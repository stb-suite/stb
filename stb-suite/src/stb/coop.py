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
   expensive full pass over the whole WFSX). Pass --fermi directly.
 - Sign convention: COOP > 0 is conventionally bonding, < 0 antibonding.
   sisl's raw COHP sign has NOT been independently verified against a
   textbook example in this suite -- some literature plots "-COHP" so
   that positive also means bonding; treat the raw sign from this tool
   as sisl's own convention until cross-checked against a known case.
 - Not tested for non-collinear states (sisl's own COOP/COHP docstring
   caveat, passed through unchanged).
"""

VERSION = "1.0.0"

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from stb.core.cli import color_text, show_intro
from stb.core.deps import require_sisl
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, iter_wfsx_states_by_k, read_k_weights
from stb.core.broadening import sigma_ev_from_args


def build_atom_pairs(geometry, pair_indices, pair_species):
    """Returns a list of (label, atoms_i, atoms_j) triples. --pair gives a
    single-index-vs-single-index pair; --pair-species aggregates ALL
    atom-index pairs between the two species sets into one combined
    curve.

    Note: a same-species --pair-species (e.g. --pair-species Sn Sn) sets
    atoms_i == atoms_j, so the aggregate curve includes each atom's
    on-site (self) term alongside the genuine inter-atomic Sn-Sn terms --
    there's no way to separate "bonding between Sn atoms" from "each Sn
    atom's own on-site character" in that aggregate. Use --pair I I
    (single on-site term) or --pair I J (single inter-atomic term) for a
    curve that isn't a mix of the two."""
    pairs = []
    na = geometry.na
    for i, j in pair_indices or []:
        if not (0 <= i < na) or not (0 <= j < na):
            raise ValueError(f"--pair {i} {j} out of range (0-{na - 1}).")
        pairs.append((f"{i}-{j}", [i], [j]))
    for sp_a, sp_b in pair_species or []:
        idx_a = [ia for ia in range(na) if geometry.atoms[ia].symbol == sp_a]
        idx_b = [ia for ia in range(na) if geometry.atoms[ia].symbol == sp_b]
        if not idx_a:
            raise ValueError(f"--pair-species: no atoms of species '{sp_a}' found.")
        if not idx_b:
            raise ValueError(f"--pair-species: no atoms of species '{sp_b}' found.")
        pairs.append((f"{sp_a}-{sp_b}", idx_a, idx_b))
    return pairs


def extract_pair_curve(total_oplist, geometry, atoms_i, atoms_j):
    """Sums the accumulated orbital-pair-resolved sparse matrices' entries
    whose row belongs to atoms_i and whose column (folded back to the unit
    cell across the full auxiliary supercell) belongs to atoms_j.

    A single mat[orb_i, :][:, col_idx_j] sum is already the COMPLETE
    bilateral interaction between atoms_i and atoms_j -- the (no, no*n_s)
    matrix's column range already spans every periodic image of atoms_j,
    including the direct (image-0) unit-cell-to-unit-cell bond, and that
    submatrix is itself Hermitian/symmetric. Verified empirically: for
    atoms_i != atoms_j, mat[orb_i,:][:,col_idx_j].sum() and
    mat[orb_j,:][:,col_idx_i].sum() are numerically identical (ratio
    1.0000000007 on a real 8-k-point test), so adding both (an earlier,
    incorrect version of this function did) doubles the true value rather
    than adding a distinct contribution."""
    orb_i = geometry.a2o(atoms_i, all=True)
    col_atoms = geometry.o2a(np.arange(geometry.no * geometry.n_s))
    col_idx_j = np.where(np.isin(col_atoms, atoms_j))[0]

    curve = np.empty(len(total_oplist))
    for ie, mat in enumerate(total_oplist):
        curve[ie] = np.real(mat[orb_i, :][:, col_idx_j].sum())
    return curve


def main():
    parser = argparse.ArgumentParser(
        description="Energy-resolved COOP/COHP bonding/antibonding curves for selected atom "
                    "pairs, from a full-BZ SIESTA .WFSX + .HSX.",
        epilog="Example usage:\n"
               "  stb-coop --label siesta --quantity coop --pair 0 1 --erange -10 5 --sigma 200\n"
               "  stb-coop --label siesta --quantity cohp --pair-species Sn O --erange -10 5 --fwhm 300\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and <label>.HSX "
                             "(mandatory -- no approximate fallback). Mutually exclusive with "
                             "--wfsx/--hsx-file.")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Explicit .HSX path. Required if --wfsx is used instead of --label.")

    parser.add_argument("--quantity", type=str, choices=["coop", "cohp"], required=True,
                        help="'coop': Crystal Orbital Overlap Population. "
                             "'cohp': Crystal Orbital Hamilton Population.")

    parser.add_argument("--pair", type=int, nargs=2, action="append", default=None,
                        metavar=("I", "J"),
                        help="0-based atom-index pair (repeatable). I == J is allowed "
                             "(on-site term).")
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
                        help="Fermi energy (eV), required if --shift fermi.")
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
                        help="Also print Hamiltonian.bond_order(method='mulliken', "
                             "projection='atom') per pair -- a cheap, energy-INTEGRATED "
                             "complementary sanity check, not a substitute for the curve.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write coop.dat/cohp.dat into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="Skip the matplotlib preview (the data file is still written).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-coop {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if not args.pair and not args.pair_species:
        parser.error("at least one --pair or --pair-species is required.")
    if args.shift == "fermi" and args.fermi is None:
        parser.error("--fermi is required when --shift is 'fermi'.")
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is 'manual'.")

    try:
        sigma_eV = sigma_ev_from_args(args.sigma, args.fwhm)
    except ValueError as e:
        parser.error(str(e))

    if args.label:
        if args.wfsx or args.hsx_file:
            parser.error("--label cannot be combined with --wfsx/--hsx-file.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.hsx_file:
        parser.error("--hsx-file is required when using --wfsx instead of --label.")

    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.WFSX', '{args.label}.selected.WFSX' and "
                 f"'{args.label}.bands.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text(f"{args.quantity.upper()} BONDING ANALYSIS:", 'bold'))
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
        pairs = build_atom_pairs(geometry, args.pair, args.pair_species)
    except ValueError as e:
        parser.error(str(e))
    print(f"[INFO] Selected pairs: {[p[0] for p in pairs]}")

    print(f"[INFO] Reading k-point weights from '{args.wfsx}' ...")
    k_weights = read_k_weights(args.wfsx, H)
    nk = len(k_weights)
    if nk <= 4 and np.allclose(k_weights, k_weights[0]):
        print(color_text(
            f"[WARNING] Only {nk} k-point(s), all equally weighted -- this doesn't look like a "
            "converged full-BZ mesh. stb-coop needs a genuine k-mesh WFSX, not a band-path one "
            "like stb-fatbands uses -- the resulting curve may not be physically meaningful.",
            'yellow'))
    k_weights = k_weights / k_weights.sum()

    E = np.linspace(args.erange[0], args.erange[1], args.npoints)
    if args.shift == "fermi":
        E_query = E + args.fermi
    elif args.shift == "manual":
        E_query = E + args.manual_value
    else:
        E_query = E
    dist = sisl.physics.get_distribution("gaussian", smearing=sigma_eV)

    print(f"[INFO] Accumulating {args.quantity.upper()}(E) over {nk} k-points x {args.npoints} "
          "energy points (this may take a moment) ...")
    total = None
    for k_index, block in iter_wfsx_states_by_k(args.wfsx, H):
        kw = k_weights[k_index]
        for spin, state in block.items():
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
    for label, atoms_i, atoms_j in pairs:
        curves[label] = extract_pair_curve(total, geometry, atoms_i, atoms_j)

    if args.bond_order:
        print(f"[INFO] Complementary bond_order (Mulliken, energy-integrated) sanity check:")
        bo = H.bond_order(method="mulliken", projection="atom")
        # bo is a SparseAtom, shape (na, na*n_s) -- sum over ALL periodic
        # images of atom j (not just its own unit-cell copy), same
        # supercell-folding treatment as extract_pair_curve() above, via
        # Geometry.asc2uc (atom-level analogue of o2a).
        col_atoms = geometry.asc2uc(np.arange(geometry.na * geometry.n_s))
        for label, atoms_i, atoms_j in pairs:
            col_idx = np.where(np.isin(col_atoms, atoms_j))[0]
            value = sum(bo[i, :][col_idx].sum() for i in atoms_i)
            print(f"       {label}: {value:.4f}")

    os.makedirs(args.output_dir, exist_ok=True)
    out_name = f"{args.quantity}.dat"
    out_path = os.path.join(args.output_dir, out_name)
    pair_labels = list(curves.keys())
    with open(out_path, "w") as f:
        f.write(f"# Generated by STB. Quantity: {args.quantity.upper()}, shift = {args.shift}\n")
        f.write(f"# Broadening: sigma = {sigma_eV * 1000:.3f} meV\n")
        header = "#" + f"{'Energy(eV)':<14}" + "".join(f"{lbl:<16}" for lbl in pair_labels)
        f.write(header + "\n")
        for ie, e in enumerate(E):
            row = f"{e:<15.6f}" + "".join(f"{curves[lbl][ie]:<16.6e}" for lbl in pair_labels)
            f.write(row + "\n")
    print(f"[INFO] Wrote {out_path}")

    if args.plot:
        plt.figure(figsize=(7, 5))
        for lbl in pair_labels:
            plt.plot(E, curves[lbl], lw=2, label=lbl)
        plt.axhline(0, color="gray", linestyle="--", linewidth=1)
        plt.xlabel("Energy (eV)")
        plt.ylabel(f"{args.quantity.upper()} (arb. units)")
        plt.title(f"{args.quantity.upper()} (positive = bonding, by convention)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
