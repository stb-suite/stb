#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Energy-resolved COOP/COHP (Crystal Orbital Overlap/Hamilton Population)
bonding/antibonding curves for user-selected atom pairs, from a SIESTA
full-BZ .WFSX + .HSX.

Uses sisl's own EigenstateElectron.COOP()/.COHP() (an "experimental" sisl
API per its own docstring), which operates one k-point at a time and
returns an orbital-pair-resolved, energy-resolved sparse matrix (shape
(no, no*n_s), spanning the full auxiliary supercell so periodic-image
bonds are included) -- there is no Brillouin-zone-averaging convenience in
sisl, so this tool loops every k-point itself, weighting by the .WFSX's
own k-point weights, then reduces the accumulated orbital-pair matrix down
to the user's chosen atom pairs by summing the relevant row/column blocks.

Needs a real sisl.physics.Hamiltonian (.HSX) -- unlike stb-wfdensity/
stb-sts, there is no usable approximate fallback: COOP needs the overlap
matrix (Sk()) and COHP needs the Hamiltonian matrix itself (Hk()), neither
of which exists on a bare Geometry.

Needs a full-BZ-sampled .WFSX (SIESTA: WriteWaveFunctions T + an explicit
%block WaveFuncKPoints listing the desired k-mesh -- named <label>.selected.WFSX
by SIESTA's own convention) -- NOT the band-path .WFSX stb-fatbands uses,
since COOP/COHP is a DOS-like, Brillouin-zone-integrated quantity.

Known limitations:
 - sisl's COOP/COHP API is explicitly marked "experimental" and can use
   substantial memory for many energy points / large systems -- start with
   a modest --npoints/--erange.
 - No metadata distinguishes a full-BZ-mesh WFSX from a band-path one;
   this tool cannot detect the difference and does not try to.
 - No --shift vbm/cbm (same reasoning as stb-sts: would need a second,
   expensive full pass over the whole WFSX). Use --shift fermi instead
   (see the Fermi-energy hierarchy below).
 - Not tested for non-collinear states (sisl's own COOP/COHP docstring
   caveat, passed through unchanged).

Output/report style rewritten (v1.0.0 -> v2.0.0) to match the rest of the
Analysis category (stb-wfdensity/stb-sts): a numbered [0]...[7] report,
--save-report, --save-gnuplot (previously wrote coop.dat/cohp.dat
unconditionally with no way to opt out, and no .gplot script at all --
now both are off by default and opted INTO together, with a real,
multi-pair .gplot script). --view replaces the old --no-plot -- the
matplotlib preview is now off by default and opted into, instead of on by
default and opted out of (same convention flip every rewritten Analysis
tool this session has gotten).

Also fixes the same --label + --hsx-file rejection bug already fixed in
stb-wfdensity/stb-sts (a real fdf/HSX is very often not literally named
<label>.*), and gives --shift fermi the same priority-ordered
Fermi-energy hierarchy those two tools have (--fermi > --bands-file >
--fermi-file > an auto-detected .out log, decoupled from --label) via the
shared core.siesta_bands.resolve_fermi_energy_hierarchy.

**A real, serious bug found and fixed in --bond-order**: this cross-check
used to call `H.bond_order(method="mulliken", ...)` directly on the
Hamiltonian object `H`. sisl's `bond_order()` computes
`2 * M_ij * S_ij` where `M` is meant to be a DENSITY MATRIX -- but `H` is
the HAMILTONIAN, so this silently substituted Hamiltonian matrix elements
(eV, typically tens of eV for on-site terms) for a density-matrix
population (dimensionless, typically O(0.1-2)). Verified live on a real
Sn3O4 oxide: the on-site (self) "bond order" came out as -92 (nonsensical
for what should be a population-like number), and a genuinely bonded Sn-O
pair (2.08 Ang) gave -38 -- barely different in magnitude from a
NON-bonded Sn-Sn pair 3.8 Ang apart, which gave -29. Two pairs with
wildly different bonding character giving near-identical "bond orders" is
the smoking gun that the numbers weren't measuring bonding at all.

Fixed by reading a REAL density matrix (new --dm-file, or an
auto-detected <label>.DM) and calling bond_order on THAT instead. One
more wrinkle: SIESTA's .DM file format doesn't store the overlap matrix
at all (confirmed live: DM.Sk() gives an exactly-zero overlap trace right
after reading a bare .DM) -- only .HSX/.TSHS do. So the already-loaded
Hamiltonian's own overlap column is spliced into the freshly-read density
matrix's sparse structure before calling bond_order (`DM._csr._D[:, -1] =
H._csr._D[:, -1]`, the standard sisl idiom for this exact combination;
confirmed by re-checking DM.Sk()'s trace afterward, which then matches
H's own ~182 for this fixture's 182-orbital basis). After the fix: the
same Sn-O bonds (~2.08-2.13 Ang) give a physically sensible bond order of
~0.51-0.55, and the non-bonded Sn-Sn pair gives essentially zero (-0.01).

This also gave the first live, independent cross-check of this suite's
own COOP sign convention: the energy-integrated COOP curve itself (a
completely different computation path, from eigenstate wavefunction
coefficients rather than the density matrix) agrees in sign with the
fixed bond order for both pairs (+0.133 vs +0.527 for the bonded pair,
-0.0008 vs -0.010 for the non-bonded one) -- confirming COOP > 0 really
does mean bonding here, not just by convention/assertion.

COHP sign: sisl's COHP() uses the raw H(k), i.e. Dronskowski-Bloechl's convention, in
which a bonding interaction is NEGATIVE (-ICOHP > 0). The report labels each quantity by
its own convention (bonding_character); before 2.0.1 the COOP rule was applied to COHP too,
and every Sn-O/X-C bond came out "antibonding".
"""

from stb import __version__ as VERSION
from stb.core.cli import new_report_path

import argparse
import os
from datetime import datetime

import numpy as np

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_sisl
from stb.core.siesta_bands import resolve_fermi_energy_hierarchy
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, iter_wfsx_states_by_k, read_k_weights
from stb.core.broadening import sigma_ev_from_args

REPORT_FILE = "stb_coop_report.txt"
BIB_FILE = "references.bib"


def build_atom_pairs(geometry, pair_indices, pair_species, image=None):
    """Returns a list of (label, atoms_i, atoms_j, image) tuples. --pair gives a
    single-index-vs-single-index pair; --pair-species aggregates ALL
    atom-index pairs between the two species sets into one combined
    curve.

    `image` (n1, n2, n3) or None selects which periodic images of atoms_j are
    summed: None = every image within the orbital range (atoms_j folded back to
    the unit cell), a tuple = that one image only (one bond). --image applies to
    --pair only. A --pair I I with no --image is the on-site term (image 0, as
    documented); I I with --image is I and its own periodic copy.

    Note: a same-species --pair-species (e.g. --pair-species Sn Sn) sets
    atoms_i == atoms_j, so the aggregate curve includes each atom's
    on-site (self) term alongside the genuine inter-atomic Sn-Sn terms --
    there's no way to separate "bonding between Sn atoms" from "each Sn
    atom's own on-site character" in that aggregate. Use --pair I I
    (single on-site term) or --pair I J (single inter-atomic term) for a
    curve that isn't a mix of the two."""
    pairs = []
    na = geometry.na
    if image is not None:
        image = tuple(int(n) for n in image)
        nsc = geometry.lattice.nsc
        if any(abs(n) > (m - 1) // 2 for n, m in zip(image, nsc)):
            raise ValueError(f"--image {' '.join(map(str, image))} is outside the auxiliary supercell "
                             f"of the Hamiltonian (nsc = {' '.join(map(str, nsc))}): no orbital of that "
                             "image interacts with the unit cell.")
    for i, j in pair_indices or []:
        if not (0 <= i < na) or not (0 <= j < na):
            raise ValueError(f"--pair {i} {j} out of range (0-{na - 1}).")
        if image is not None:
            pairs.append((f"{i}-{j}[{','.join(map(str, image))}]", [i], [j], image))
        elif i == j:
            pairs.append((f"{i}-{j}", [i], [j], (0, 0, 0)))
        else:
            pairs.append((f"{i}-{j}", [i], [j], None))
    for sp_a, sp_b in pair_species or []:
        idx_a = [ia for ia in range(na) if geometry.atoms[ia].symbol == sp_a]
        idx_b = [ia for ia in range(na) if geometry.atoms[ia].symbol == sp_b]
        if not idx_a:
            raise ValueError(f"--pair-species: no atoms of species '{sp_a}' found.")
        if not idx_b:
            raise ValueError(f"--pair-species: no atoms of species '{sp_b}' found.")
        pairs.append((f"{sp_a}-{sp_b}", idx_a, idx_b, None))
    return pairs


def pair_columns(geometry, atoms_j, image=None):
    """Column indices (orbitals of the auxiliary supercell, no * n_s) of atoms_j: all of
    their periodic images (image None) or the one image (n1, n2, n3). geometry.o2a gives
    SUPERCELL atom indices (0 .. na*n_s - 1), so they are folded with asc2uc before the
    comparison; before 2.0.1 they were not, and only image 0 of J was ever summed -- a bond
    to a neighbouring cell (the C-C bridge of XC8, 1.5 Ang) came out ICOOP ~ 0 while its
    bond order was 0.88."""
    if image is not None:
        isc = geometry.lattice.sc_index(image)
        return np.asarray(geometry.a2o(atoms_j, all=True)) + isc * geometry.no
    col_atoms = geometry.asc2uc(geometry.o2a(np.arange(geometry.no * geometry.n_s)))
    return np.where(np.isin(col_atoms, atoms_j))[0]


def bond_order_columns(geometry, atoms_j, image=None):
    """The same selection as pair_columns, in supercell ATOM indices (for the bond-order matrix)."""
    if image is not None:
        isc = geometry.lattice.sc_index(image)
        return np.asarray(atoms_j) + isc * geometry.na
    col_atoms = geometry.asc2uc(np.arange(geometry.na * geometry.n_s))
    return np.where(np.isin(col_atoms, atoms_j))[0]


def image_breakdown(total_oplist, geometry, E, i, j, tol=1e-4):
    """[(image, distance Ang, integrated value)] of every periodic image of atom j that
    interacts with atom i, sorted by distance -- what a folded --pair I J sums, one image
    at a time. Images whose integral is below `tol` in magnitude are left out."""
    rows = []
    orb_i = geometry.a2o([i], all=True)
    xyz_i = geometry.xyz[i]
    # Energy integral of every column summed over the rows of atom i, in one pass over the
    # energy points (trapezoid weights), then one block of columns per image.
    weights = np.zeros(len(E))
    if len(E) > 1:
        dE = np.diff(E)
        weights[:-1] += dE / 2
        weights[1:] += dE / 2
    integrated_cols = np.zeros(geometry.no * geometry.n_s)
    for w, m in zip(weights, total_oplist):
        integrated_cols += w * np.real(np.asarray(m[orb_i, :].sum(axis=0))).ravel()
    orb_j = np.asarray(geometry.a2o([j], all=True))
    for isc, off in enumerate(geometry.lattice.sc_off):
        if i == j and isc == 0:
            continue
        value = float(integrated_cols[orb_j + isc * geometry.no].sum())
        if abs(value) < tol:
            continue
        r = geometry.xyz[j] + off @ geometry.cell - xyz_i
        rows.append((tuple(int(n) for n in off), float(np.linalg.norm(r)), value))
    rows.sort(key=lambda row: row[1])
    return rows


def population_factor(H):
    """sisl's COOP/COHP is per spin channel (occupation 1 per state) and one-sided (i -> j);
    the Mulliken bond order counts both directions and, for a spin-unpolarized run, both
    spins. Returns 4 (unpolarized) or 2 (polarized, non-collinear, spin-orbit): bond order =
    factor x ICOOP for a window from the bottom of the bands to E_F (Sn3O4: 3.98; the sum
    over every orbital pair is N_electrons / 2)."""
    spin = H.spin
    return 4 if spin.is_unpolarized else 2


def bonding_character(quantity, integrated):
    """'bonding'/'antibonding'/'non-bonding' of an energy-integrated curve. COOP > 0 is
    bonding (Hughbanks-Hoffmann); sisl's COHP is built from the raw H(k), Dronskowski's
    convention, so COHP < 0 is bonding (-ICOHP > 0). The old label applied the COOP rule
    to both and called every COHP bond 'antibonding'."""
    measure = integrated if quantity == "coop" else -integrated
    return "bonding" if measure > 0 else "antibonding" if measure < 0 else "non-bonding"


def extract_pair_curve(total_oplist, geometry, atoms_i, atoms_j, image=None):
    """Sums the accumulated orbital-pair-resolved sparse matrices' entries
    whose row belongs to atoms_i and whose column belongs to atoms_j: every
    periodic image of atoms_j (image None; folded back to the unit cell with
    asc2uc, see pair_columns) or the one image given.

    A single mat[orb_i, :][:, col_idx_j] sum is the one-sided (i -> j) term of
    the interaction; the j -> i term is numerically identical (ratio
    1.0000000007 on a real 8-k-point test), so adding both (an earlier,
    incorrect version of this function did) doubled the value. The Mulliken
    bond order counts both directions and both spins: for a spin-unpolarized
    run, bond order = 4 x ICOOP integrated to E_F (checked on Sn3O4: 3.98)."""
    orb_i = geometry.a2o(atoms_i, all=True)
    col_idx_j = pair_columns(geometry, atoms_j, image)

    curve = np.empty(len(total_oplist))
    for ie, mat in enumerate(total_oplist):
        curve[ie] = np.real(mat[orb_i, :][:, col_idx_j].sum())
    return curve


def resolve_density_matrix(H, label, dm_file):
    """Returns (DM, source) or (None, None) -- a real sisl DensityMatrix
    for --bond-order, built from an explicit --dm-file or an
    auto-detected <label>.DM.

    SIESTA's .DM file format does not store the overlap matrix (verified
    live: DM.Sk() is exactly zero right after a bare read) -- only
    .HSX/.TSHS do. So the already-loaded Hamiltonian `H`'s own overlap
    column is spliced into the freshly-read DM's sparse data before
    returning it; without this, DM.bond_order() would silently use a
    zero overlap and return all-zero bond orders instead of erroring.
    """
    path = dm_file
    if path is None and label:
        candidate = f"{label}.DM"
        if os.path.isfile(candidate):
            path = candidate
    if path is None or not os.path.isfile(path):
        return None, None

    sisl = require_sisl()
    DM = sisl.get_sile(path).read_density_matrix(geometry=H.geometry)
    if DM._csr.nnz != H._csr.nnz:
        raise ValueError(
            f"'{path}' and the loaded Hamiltonian have different sparse-matrix sizes "
            f"({DM._csr.nnz} vs {H._csr.nnz} non-zeros) -- they likely aren't from the same "
            "calculation/geometry, so the overlap matrix can't be safely spliced in."
        )
    DM._csr._D[:, -1] = H._csr._D[:, -1]
    return DM, path


def main():
    parser = argparse.ArgumentParser(
        description="Energy-resolved COOP/COHP bonding/antibonding curves for selected atom "
                    "pairs, from a full-BZ SIESTA .WFSX + .HSX.",
        epilog="Example usage:\n"
               "  stb-coop --label siesta --quantity coop --pair 0 1 --erange -10 5 --sigma 200\n"
               "  stb-coop --label siesta --quantity cohp --pair-species Sn O --erange -10 5 --fwhm 300\n"
               "  stb-coop --label siesta --quantity coop --pair 0 1 --erange -10 5 --sigma 200 \\\n"
               "      --bond-order --save-gnuplot --view\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and, unless "
                             "--hsx-file overrides it, <label>.HSX (mandatory -- no approximate "
                             "fallback). Mutually exclusive with --wfsx (not --hsx-file -- see "
                             "its own help).")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Explicit .HSX path. Required if --wfsx is used instead of --label; "
                             "optional (and common) together WITH --label too, since the real "
                             "file is often not literally named <label>.HSX -- --label still "
                             "auto-detects the .WFSX in that case, only the Hamiltonian source "
                             "changes.")

    parser.add_argument("--quantity", type=str, choices=["coop", "cohp"], required=True,
                        help="'coop': Crystal Orbital Overlap Population. "
                             "'cohp': Crystal Orbital Hamilton Population.")

    parser.add_argument("--pair", type=int, nargs=2, action="append", default=None,
                        metavar=("I", "J"),
                        help="0-based atom-index pair (repeatable). Sums EVERY periodic image of "
                             "J within the orbital range (the report lists each image with its "
                             "distance); use --image for one bond. I == J is allowed (on-site "
                             "term).")
    parser.add_argument("--image", type=int, nargs=3, default=None, metavar=("N1", "N2", "N3"),
                        help="With --pair: only the periodic image (N1, N2, N3) of J, i.e. the one "
                             "bond between I in the unit cell and J in cell N1 a1 + N2 a2 + N3 a3 "
                             "(e.g. --pair 4 5 --image 0 -1 0 for a bond across the cell edge). "
                             "The report's image table gives the image of each neighbour.")
    parser.add_argument("--pair-species", type=str, nargs=2, action="append", default=None,
                        metavar=("A", "B"),
                        help="Species-pair (repeatable), e.g. --pair-species Sn O -- "
                             "aggregates ALL atom-index pairs between the two species into "
                             "one combined curve. If A == B (e.g. Sn Sn), the aggregate also "
                             "includes each atom's own on-site term, mixed in with the "
                             "inter-atomic Sn-Sn terms -- use --pair I J for a purely "
                             "inter-atomic curve instead.")

    parser.add_argument("--erange", type=float, nargs=2, required=True, metavar=("EMIN", "EMAX"),
                        help="Energy window (eV), relative to --shift's reference.")
    parser.add_argument("--npoints", type=int, default=300,
                        help="Number of energy samples (default: 300). sisl's COOP/COHP is "
                             "memory-heavy for many points -- start modest.")
    parser.add_argument("--shift", type=str, choices=["none", "fermi", "manual"], default="none",
                        help="Energy reference (default: none). No vbm/cbm option -- see the "
                             "module docstring.")
    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy (eV) for --shift fermi. Highest-priority source if given.")
    parser.add_argument("--bands-file", type=str, default=None,
                        help="Companion .bands file to read the Fermi energy from, for --shift "
                             "fermi.")
    parser.add_argument("--fermi-file", type=str, default=None,
                        help="Explicit SIESTA .out log to read the Fermi energy from, for --shift "
                             "fermi -- an alternative to --bands-file for a run with no saved "
                             ".bands file. Not assumed to be named after --label.")
    parser.add_argument("--manual-value", type=float, default=None,
                        help="Custom shift value (eV), required if --shift manual.")

    broadening = parser.add_mutually_exclusive_group(required=True)
    broadening.add_argument("--sigma", type=float, default=None,
                        help="Gaussian broadening standard deviation, in meV.")
    broadening.add_argument("--fwhm", type=float, default=None,
                        help="Gaussian broadening full width at half maximum, in meV.")
    parser.add_argument("--atol", type=float, default=1e-10,
                        help="Per-state energy-window cutoff passed straight to sisl's "
                             "COOP/COHP (default: 1e-10) -- a performance knob, not a physics "
                             "tolerance.")

    parser.add_argument("--bond-order", action="store_true",
                        help="Also compute Hamiltonian.bond_order(method='mulliken', "
                             "projection='atom') per pair from a REAL density matrix -- a "
                             "cheap, energy-INTEGRATED complementary sanity check, not a "
                             "substitute for the curve. Needs --dm-file or an auto-detected "
                             "<label>.DM (see below).")
    parser.add_argument("--dm-file", type=str, default=None,
                        help="Explicit .DM density-matrix path, for --bond-order. Falls back to "
                             "an auto-detected <label>.DM if omitted. SIESTA's .DM format doesn't "
                             "store the overlap matrix -- it's spliced in from the already-loaded "
                             "Hamiltonian, so the .DM must be from the same calculation as the "
                             ".HSX/--hsx-file.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write coop.dat/cohp.dat (and, with --save-gnuplot/"
                             "the .gplot script/references.bib) into (the report goes to the current folder) "
                             "(default: current directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a real .gplot script alongside coop.dat/cohp.dat. Off by "
                             "default -- this tool used to write the .dat file unconditionally "
                             "but never wrote a .gplot script at all.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview before finishing. Off by "
                             "default (replaces the old --no-plot, which was on by default).")

    parser.add_argument("-v", "--version", action="version", version=f"stb-coop {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if not args.pair and not args.pair_species:
        parser.error("at least one --pair or --pair-species is required.")
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is 'manual'.")

    try:
        sigma_eV = sigma_ev_from_args(args.sigma, args.fwhm)
    except ValueError as e:
        parser.error(str(e))

    if args.label:
        if args.wfsx:
            parser.error("--label cannot be combined with --wfsx.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.hsx_file:
        parser.error("--hsx-file is required when using --wfsx instead of --label.")
    # --label + --hsx-file together is valid and common -- see the module docstring (the same
    # fix already made for stb-wfdensity/stb-sts): load_parent() below already prefers an
    # explicit hsx_file over <label>.HSX on its own.

    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.WFSX', '{args.label}.selected.WFSX' and "
                 f"'{args.label}.bands.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    fermi_energy = fermi_source = None
    if args.shift == "fermi":
        fermi_energy, fermi_source = resolve_fermi_energy_hierarchy(
            args.fermi, args.bands_file, args.fermi_file, args.label)
        if fermi_energy is None:
            parser.error(
                "--shift fermi needs a Fermi energy -- none found. Pass --fermi <value>, "
                "--bands-file <path>, --fermi-file <path>, or leave a SIESTA .out log "
                "(<label>.out, or the sole *.out) in the current directory."
            )

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text(f"{args.quantity.upper()}: reading input data", 'bold'))
    print("-" * 60)

    sisl = require_sisl()

    print("[INFO] Resolving Hamiltonian (mandatory -- COOP/COHP need a real overlap/Hamiltonian) ...")
    H, source, _ = load_parent(args.label, args.hsx_file, None, mode="hamiltonian_required")
    if H is None:
        parser.error("No .HSX found -- stb-coop needs --hsx-file, or --label with a "
                     "<label>.HSX next to it. There is no approximate fallback for COOP/COHP.")
    print(f"[INFO] Using '{source}' ({H.geometry.no} orbitals, {H.geometry.na} atoms).")
    geometry = H.geometry

    try:
        if args.image is not None and not args.pair:
            raise ValueError("--image needs --pair (it selects one periodic image of J).")
        pairs = build_atom_pairs(geometry, args.pair, args.pair_species, args.image)
    except ValueError as e:
        parser.error(str(e))
    print(f"[INFO] Selected pairs: {[p[0] for p in pairs]}")

    DM = dm_source = None
    if args.bond_order:
        try:
            DM, dm_source = resolve_density_matrix(H, args.label, args.dm_file)
        except ValueError as e:
            parser.error(str(e))
        if DM is None:
            parser.error(
                "--bond-order needs a real density matrix -- pass --dm-file <path>, or leave "
                "<label>.DM next to your calculation. (SIESTA's .DM file doesn't store the "
                "overlap matrix itself; it's spliced in from the already-loaded Hamiltonian, "
                "so bond_order can never be computed from --hsx-file/<label>.HSX alone.)"
            )
        print(f"[INFO] Using '{dm_source}' for --bond-order (overlap spliced in from '{source}').")

    print(f"[INFO] Reading k-point weights from '{args.wfsx}' ...")
    k_weights = read_k_weights(args.wfsx, H)
    nk = len(k_weights)
    mesh_warning = None
    if nk <= 4 and np.allclose(k_weights, k_weights[0]):
        mesh_warning = (
            f"Only {nk} k-point(s), all equally weighted -- this doesn't look like a "
            "converged full-BZ mesh. stb-coop needs a genuine k-mesh WFSX, not a band-path one "
            "like stb-fatbands uses -- the resulting curve may not be physically meaningful.")
        print(color_text(f"[WARNING] {mesh_warning}", 'yellow'))
    k_weights = k_weights / k_weights.sum()

    E = np.linspace(args.erange[0], args.erange[1], args.npoints)
    if args.shift == "fermi":
        E_query = E + fermi_energy
    elif args.shift == "manual":
        E_query = E + args.manual_value
    else:
        E_query = E
    dist = sisl.physics.get_distribution("gaussian", smearing=sigma_eV)

    print(f"[INFO] Accumulating {args.quantity.upper()}(E) over {nk} k-points x {args.npoints} "
          "energy points (this may take a moment) ...")
    total = None
    lowest_eig = np.inf
    for k_index, block in iter_wfsx_states_by_k(args.wfsx, H):
        kw = k_weights[k_index]
        for spin, state in block.items():
            lowest_eig = min(lowest_eig, float(np.min(state.eig)))
            if args.quantity == "coop":
                result = state.COOP(E_query, distribution=dist, atol=args.atol)
            else:
                result = state.COHP(E_query, distribution=dist, atol=args.atol)
            if total is None:
                total = [kw * m for m in result]
            else:
                total = [t + kw * m for t, m in zip(total, result)]

    print("[INFO] Reducing to selected atom pairs ...")
    curves = {}
    breakdowns = {}
    for label, atoms_i, atoms_j, image in pairs:
        curves[label] = extract_pair_curve(total, geometry, atoms_i, atoms_j, image)
        if image is None and len(atoms_i) == 1 and len(atoms_j) == 1:
            breakdowns[label] = image_breakdown(total, geometry, E, atoms_i[0], atoms_j[0])
    pair_labels = list(curves.keys())

    bond_orders = {}
    if DM is not None:
        bo = DM.bond_order(method="mulliken", projection="atom")
        for label, atoms_i, atoms_j, image in pairs:
            col_idx = bond_order_columns(geometry, atoms_j, image)
            bond_orders[label] = float(sum(bo[i, :][col_idx].sum() for i in atoms_i))

    os.makedirs(args.output_dir, exist_ok=True)
    out_name = f"{args.quantity}.dat"
    out_path = os.path.join(args.output_dir, out_name)
    with open(out_path, "w") as f:
        f.write(f"# Generated by STB (stb-coop). Quantity: {args.quantity.upper()}, shift = {args.shift}\n")
        f.write(f"# Broadening: sigma = {sigma_eV * 1000:.3f} meV\n")
        header = "#" + f"{'Energy(eV)':<14}" + "".join(f"{lbl:<16}" for lbl in pair_labels)
        f.write(header + "\n")
        for ie, e in enumerate(E):
            row = f"{e:<15.6f}" + "".join(f"{curves[lbl][ie]:<16.6e}" for lbl in pair_labels)
            f.write(row + "\n")

    gplot_path = None
    if args.save_gnuplot:
        gplot_path = out_path.rsplit('.', 1)[0] + ".gplot"
        dat_base = os.path.basename(out_path)
        stem = os.path.splitext(dat_base)[0]
        ncols = len(pair_labels) + 1
        lines = [
            'set terminal pdfcairo enhanced font "Arial,14" size 8,6\n',
            f'set output "{stem}.pdf"\n',
            'set xlabel "Energy (eV)"\n',
            f'set ylabel "{args.quantity.upper()} (arb. units)"\n',
            'set key outside\n',
            'set grid\n',
            'set arrow from graph 0, first 0 to graph 1, first 0 nohead lc rgb "gray" dt 2\n',
            f'plot for [i=2:{ncols}] "{dat_base}" using 1:i with lines lw 2 title columnheader(i)\n',
        ]
        with open(gplot_path, 'w') as f:
            f.writelines(lines)

    # --- From here on: the numbered, save-able report --------------------
    report_path = new_report_path("coop")
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text(f"\n===== STB-COOP REPORT ({args.quantity.upper()}) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}" if args.label else "Label          : (explicit --wfsx/--hsx-file)", f_out)
    print_dual(f"Quantity       : {args.quantity.upper()}", f_out)
    print_dual(f"WFSX file      : {args.wfsx}", f_out)
    print_dual(f"Hamiltonian src: {source}", f_out)
    print_dual(f"Energy shift   : {args.shift}" + (f" (Fermi source: {fermi_source})" if args.shift == "fermi" else ""), f_out)
    print_dual(f"Broadening     : sigma = {sigma_eV * 1000:.3f} meV", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Orbitals (basis)", f"{geometry.no}"], None),
        (["Atoms", f"{geometry.na}"], None),
        (["k-points in .WFSX", f"{nk}"], None),
        (["Energy samples", f"{args.npoints}"], None),
    ], f_out)
    if mesh_warning:
        print_dual(color_text(f"[WARNING] {mesh_warning}", 'yellow'), f_out)

    print_section("[2] PAIR SELECTION", f_out)
    print_dual(f"Selected pairs: {pair_labels}", f_out)

    print_section(f"[3] {args.quantity.upper()} CURVE", f_out)
    g_factor = population_factor(H)
    rows = []
    for lbl in pair_labels:
        curve = curves[lbl]
        peak_idx = int(np.argmax(np.abs(curve)))
        integrated = float(np.trapezoid(curve, E))
        character = bonding_character(args.quantity, integrated)
        rows.append(([lbl, f"{curve[peak_idx]:.6e} at E={E[peak_idx]:.4f} eV",
                      f"{integrated:.6e} ({character})", f"{g_factor * integrated:+.4f}"], None))
    print_table(["Pair", "Largest |value|", "Integrated over --erange", f"x{g_factor:g} (both spins, i<->j)"],
                rows, f_out)
    print_dual(f"Normalization: the curve and {args.quantity}.dat are sisl's value per spin channel and one "
               f"direction (i -> j). The last column multiplies by {g_factor:g} = "
               f"{'2 spins x ' if g_factor == 4 else ''}2 directions (i -> j and j -> i): for COOP it is "
               "on the scale of the Mulliken bond order (equal to it when the window runs from the bottom of "
               "the bands to E_F), and it is the scale of tables that count both spins and both directions "
               "(compare with care: codes differ in what they count).", f_out)
    for lbl, rows_img in breakdowns.items():
        if len(rows_img) < 1:
            continue
        print_dual(f"Pair {lbl}: every periodic image of atom {lbl.split('-')[1]} that interacts with "
                   f"atom {lbl.split('-')[0]} (the curve above sums them all; --pair {lbl.replace('-', ' ')} "
                   "--image N1 N2 N3 gives one):", f_out)
        print_table(["Image (N1 N2 N3)", "Distance (Ang)", f"Integrated {args.quantity.upper()}"],
                    [([" ".join(f"{n:+d}" for n in img), f"{dist:.3f}", f"{val:+.6e}"], None)
                     for img, dist, val in rows_img], f_out)
    print_dual("(Integrated over --erange, NOT necessarily just the occupied states. Sign convention: "
               "COOP > 0 = bonding (Hughbanks-Hoffmann); COHP < 0 = bonding (Dronskowski-Bloechl), "
               "i.e. -ICOHP > 0 for a bond.)", f_out)

    print_section("[4] BOND ORDER (Mulliken, energy-integrated cross-check)", f_out)
    if bond_orders:
        print_table(["Pair", "Mulliken bond order"], [
            ([lbl, f"{bond_orders[lbl]:.4f}"], None) for lbl in pair_labels
        ], f_out)
        print_dual(f"Density matrix : {dm_source}", f_out)
        name = "COOP" if args.quantity == "coop" else "-COHP"
        print_dual(f"Cross-check against the integrated {name} above (same sign expected for a "
                   "converged full-BZ mesh spanning the true occupied bandwidth):", f_out)
        for lbl in pair_labels:
            integrated = float(np.trapezoid(curves[lbl], E))
            bonding_measure = integrated if args.quantity == "coop" else -integrated
            agree = (bonding_measure > 0) == (bond_orders[lbl] > 0)
            print_dual(f"  {lbl}: integrated {name} {'positive' if bonding_measure > 0 else 'negative'}, "
                       f"bond order {'positive' if bond_orders[lbl] > 0 else 'negative'} "
                       f"-> {'[OK] signs agree' if agree else '[WARNING] signs disagree'}", f_out)
        if args.quantity == "coop":
            # The sign alone missed a bond across the cell edge (ICOOP ~ 0, bond order 0.88): compare
            # the magnitude too, x g_factor ICOOP / bond order, which is 1 for a window from the bottom
            # of the bands to E_F.
            e_top = E_query[-1]
            full_window = (args.shift == "fermi" and abs(args.erange[1]) <= 0.05
                           and E_query[0] < lowest_eig - 3 * sigma_eV)
            print_dual(f"Magnitude: x{g_factor:g} ICOOP / bond order (1 for a window from the bottom of the "
                       "bands to E_F" + ("):" if full_window else
                                         f"; this window, {E_query[0]:.2f} to {e_top:.2f} eV absolute, is not "
                                         "that one -- --shift fermi --erange <below the lowest band> 0 -- so "
                                         "the ratio is shown, not judged):"), f_out)
            for lbl in pair_labels:
                integrated = float(np.trapezoid(curves[lbl], E))
                bo_val = bond_orders[lbl]
                if abs(bo_val) < 0.02:
                    print_dual(f"  {lbl}: bond order {bo_val:+.4f} -- too small for a ratio", f_out)
                    continue
                ratio = g_factor * integrated / bo_val
                off = abs(ratio - 1) > 0.1
                tag = ("[WARNING] differs from 1 by more than 10 %" if (off and full_window)
                       else "[OK]" if full_window else "")
                print_dual(f"  {lbl}: {g_factor * integrated:+.4f} / {bo_val:+.4f} = {ratio:.3f} {tag}", f_out)
    elif args.bond_order:
        print_dual("Not computed (no usable density matrix -- see the error above).", f_out)
    else:
        print_dual("Not requested (pass --bond-order, plus --dm-file or an auto-detected "
                   "<label>.DM, for this cross-check).", f_out)

    print_section("[5] OUTPUT DATA & PLOTS", f_out)
    print_dual(color_text(f"[OK] Data written to '{out_path}'.", 'green'), f_out)
    if args.save_gnuplot:
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Gnuplot script not written (off by default -- pass --save-gnuplot to write it).", f_out)

    print_section("[6] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    bib_entries.append(citations.COOP_HUGHBANKS if args.quantity == "coop" else citations.COHP_DRONSKOWSKI)
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Data           : {out_path}", f_out)
    if gplot_path:
        print_dual(f"Gnuplot script : {gplot_path}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 5))
        for lbl in pair_labels:
            plt.plot(E, curves[lbl], lw=2, label=lbl)
        plt.axhline(0, color="gray", linestyle="--", linewidth=1)
        plt.xlabel("Energy (eV)")
        plt.ylabel(f"{args.quantity.upper()} (arb. units)")
        plt.title("COOP (positive = bonding)" if args.quantity == "coop"
                  else "COHP (negative = bonding)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
