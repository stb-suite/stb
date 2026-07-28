#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Orbital-projected ("fatbands") band structure: combines a SIESTA .bands
file (k-path, eigenvalues) with a .WFSX wavefunction-coefficient file to
color/size each (k, band) point by its orbital, atom, or species character.

Needs sisl (like stb-bands' --eig-file path), and a .WFSX file written
alongside the same calculation that wrote .bands (`WFS.Write.For.Bands T`
in the fdf -- SIESTA/sisl's own convention names it <label>.bands.WFSX, not
<label>.WFSX). If a companion .HSX (Hamiltonian + overlap) is also
available, weights are computed with SIESTA's true non-orthogonal-basis
overlap matrix (physically correct); otherwise they fall back to an
implicit-orthogonal approximation (raw |c|^2 per orbital) with a printed
warning, since sisl's norm2() silently uses an identity overlap whenever
its `parent` is a bare Geometry rather than a Hamiltonian.

Known limitations (mirrors dos.py's own documented orbital-naming caveats,
since both tools share stb.core.orbitals):
 - l > 3 (g-orbitals and beyond) are excluded from --projection l/ml/
   species_l (counted and warned about once), but still contribute fully
   to --projection atom/species, which don't depend on l.
 - non-collinear/SOC (nspin 4/8) inputs are not specially handled beyond
   whatever sisl's own norm2()/Sk() already account for internally; spin
   labeling stays collinear-only (same gap dos.py documents for PDOS.xml).
 - There is no file metadata tying a given .WFSX to a given .bands file.
   This tool trusts file-order correspondence, validated only by a
   k-point/orbital-count guard and a --k-tol eigenvalue cross-check at a
   few sampled k-points -- always let --label auto-detect paired files from
   the same run rather than mixing --file/--wfsx from different
   calculations.

Output/report style rewritten (v1.0.0 -> v2.0.0) to match the rest of the
Analysis category (stb-bands/stb-dos): a numbered [0]...[6] report,
--save-report to persist it, --save-gnuplot to opt into writing
fatbands_<category>.dat/fatbands.gplot (previously unconditional on every
run), and --view to opt into the matplotlib preview (previously always
shown, blocking, with no way to skip it). Also adds --projection
species_l, combining chemical species AND angular momentum in one category
(e.g. 'Sn-s', 'O-p') -- --projection species alone cannot show, on its own,
which orbital character within a species dominates a given band; l/ml
alone cannot show which species that character belongs to. Marker sizes in
both the matplotlib and gnuplot output are smaller than
core/band_scatter.py's shared defaults (stb-ipr/stb-stm are unaffected --
see that module's own docstring): a fatbands plot scatters every
(k, band[, category]) point at once, several times denser than those
tools' single-series plots, and read as a cluttered blob at the shared
default size.

Spin-polarized (nspin=2) case verified live (v2.0.0 -> v2.1.0) against a
real SIESTA calculation (a spin-polarized isolated O atom, converged to its
physical 2 Bohr-magneton triplet moment, with WFS.Write.For.Bands + SaveHS):
the original weight loop merged both spin channels into one category with
no spin label at all -- silently combining two physically very different
band sets (here, spin-up/spin-down CBM differed by ~29 eV) into one
indistinguishable scatter series/data file. Every category is now split
into "<category>_up"/"<category>_down" for nspin=2 (same convention
dos.py already uses for its own spin-resolved PDOS columns) -- e.g.
--projection species_l on a spin-polarized run yields 'Sn-s_up',
'Sn-s_down', etc. nspin=1 is completely unaffected (no suffix). Also
default --projection changed from 'l' to 'species_l' -- the more
informative combined category is now what a plain run shows without
having to know to ask for it.
"""

VERSION = "2.1.0"

import argparse
import os
import sys
from datetime import datetime

import numpy as np

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.siesta_bands import read_data, _is_gamma, shift_bands, cbm_vbm
from stb.core.orbitals import ORBITAL_ORDER, get_orbital_name, get_detailed_orbital_name
from stb.core.siesta_wfsx import (
    resolve_wfsx_path, load_parent, _geometry_of, read_wfsx_states,
)
from stb.core.band_scatter import (
    write_scalar_data, write_scalar_gplot, plot_scalar_on_bands, plot_multi_series_on_bands,
)

REPORT_FILE = "stb_fatbands_report.txt"
BIB_FILE = "references.bib"

CATEGORY_CHOICES = ("l", "ml", "atom", "species", "species_l")

# Deliberately smaller than core/band_scatter.py's own defaults (10 / 200)
# -- see the module docstring above for why.
MARKER_SIZE_BASE = 2.0
MARKER_SIZE_SCALE = 40.0
GNUPLOT_POINT_SCALE = 2.5


def build_orbital_table(geometry, projection):
    """Per-global-orbital-index parallel lists: category name (per
    `projection`), atom index, species symbol. Orbitals with l > 3 get
    category=None for l/ml/species_l (excluded from those sums), but atom
    index/species are always recorded since --projection atom/species
    don't depend on l at all. Also returns (n_excluded, excluded_l) so the
    caller can report the same exclusion count both live and in the saved
    report, instead of only a one-off stderr print."""
    categories = []
    atoms = []
    species = []
    n_excluded = 0
    excluded_l = set()
    for ia, io in geometry.iter_orbitals(local=True):
        atom = geometry.atoms[ia]
        orb = atom.orbitals[io]
        atoms.append(ia)
        species.append(atom.symbol)
        if projection == "atom":
            categories.append(str(ia))
        elif projection == "species":
            categories.append(atom.symbol)
        elif projection == "l":
            name = get_orbital_name(orb.l)
            if name is None:
                n_excluded += 1
                excluded_l.add(orb.l)
            categories.append(name)
        elif projection == "species_l":
            name = get_orbital_name(orb.l)
            if name is None:
                n_excluded += 1
                excluded_l.add(orb.l)
                categories.append(None)
            else:
                categories.append(f"{atom.symbol}-{name}")
        else:  # ml
            name = get_detailed_orbital_name(orb.l, orb.m)
            if name is None:
                n_excluded += 1
                excluded_l.add(orb.l)
            categories.append(name)

    if n_excluded and projection in ("l", "ml", "species_l"):
        print(f"[WARNING] {n_excluded} orbital(s) with l={sorted(excluded_l)} "
              "(l > 3, i.e. g-orbitals or beyond, are not supported) excluded "
              f"from --projection {projection} -- still counted in --projection "
              "atom/species, which don't depend on l.", file=sys.stderr)

    return categories, atoms, species, n_excluded, sorted(excluded_l)


def _category_sort_key(projection):
    if projection in ("l", "ml"):
        def key(c):
            return (ORBITAL_ORDER.index(c) if c in ORBITAL_ORDER else 999, c)
        return key
    if projection == "atom":
        return lambda c: int(c)
    if projection == "species_l":
        def key(c):
            sp, _, orb = c.partition("-")
            return (sp, ORBITAL_ORDER.index(orb) if orb in ORBITAL_ORDER else 999)
        return key
    return lambda c: c  # species: plain alphabetical


def main():
    parser = argparse.ArgumentParser(
        description="Orbital-projected (fatbands) band structure from .bands + .WFSX.",
        epilog="Example usage:\n"
               "  stb-fatbands --label siesta --shift fermi --projection l\n"
               "  stb-fatbands --label siesta --shift fermi --projection species --category O Sn\n"
               "  stb-fatbands --label siesta --shift fermi --projection species_l\n"
               "  stb-fatbands --file siesta.bands --wfsx siesta.bands.WFSX --hsx-file siesta.HSX \\\n"
               "      --shift fermi --save-report --save-gnuplot --view",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Shorthand for --file <label>.bands, "
                             "auto-detects <label>.bands.WFSX (falling back to <label>.WFSX "
                             "with a note), and auto-detects <label>.HSX (optional, upgrades "
                             "weight accuracy) then <label>.fdf (geometry+basis fallback -- NOT "
                             "<label>.XV, which has no per-orbital basis information). "
                             "Mutually exclusive with --file/--wfsx/--hsx-file/--geometry-file.")

    parser.add_argument("--file", dest="input_file", type=str, default=None,
                        help="Explicit path to the .bands file (alternative to --label).")

    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file. Required if --file is used "
                             "instead of --label.")

    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Optional explicit .HSX path (Hamiltonian+overlap) for physically "
                             "correct overlap-weighted norms. Falls back to --geometry-file / "
                             "auto-detected .fdf with an accuracy warning if omitted.")

    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Optional explicit .fdf geometry path (must have the basis's "
                             ".ion/.ion.xml files alongside it), used as the fallback parent "
                             "when no .HSX is available. A bare .XV has no per-orbital basis "
                             "information and will fail the orbital-count check against the "
                             ".WFSX -- use --hsx-file or a .fdf instead. Requires sisl.")

    parser.add_argument("--shift", type=str, choices=["vbm", "cbm", "fermi", "manual"], required=True,
                        help="Reference energy shift, same vocabulary as stb-bands:\n"
                             "  - 'vbm'    : Valence Band Maximum\n"
                             "  - 'cbm'    : Conduction Band Minimum\n"
                             "  - 'fermi'  : Fermi level\n"
                             "  - 'manual' : Custom shift value (requires --manual-value).")

    parser.add_argument("--manual-value", type=float,
                        help="Custom energy shift value (required if --shift manual is used).")

    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification when --shift "
                             "is vbm/cbm (default: 0.01). Same meaning as stb-bands.")

    parser.add_argument("--projection", type=str, choices=CATEGORY_CHOICES, default="species_l",
                        help="Category to color/size markers by (default: species_l):\n"
                             "  l:         aggregate by angular momentum (s, p, d, f)\n"
                             "  ml:        detailed orbital (s, px, py, pz, dxy, ...)\n"
                             "  atom:      per atom index\n"
                             "  species:   per chemical species\n"
                             "  species_l: chemical species AND angular momentum combined\n"
                             "             (e.g. 'Sn-s', 'Sn-p', 'O-s', 'O-p') -- shows which\n"
                             "             orbital character within a species dominates a band,\n"
                             "             which plain 'species' or plain 'l' cannot on their own\n"
                             "Single choice only -- unlike stb-dos's --type, the plot needs one "
                             "weight per marker; run again for another projection.")

    parser.add_argument("--category", type=str, nargs="+", default=None,
                        help="Restrict output to specific category values (e.g. --projection "
                             "species --category O Sn, or --projection species_l --category "
                             "Sn-s O-p). Default: all categories found.")

    parser.add_argument("--k-tol", type=float, default=0.001,
                        help="Cross-check tolerance (eV) between WFSX-stored eigenvalues and "
                             ".bands eigenvalues at sampled k-points (default: 0.001). This is a "
                             "file-order sanity guard, not a physics tolerance like --gap-tol -- "
                             "raising it only weakens the check that .bands/.WFSX truly "
                             "correspond to the same calculation.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write fatbands_<category>.dat/fatbands.gplot (with "
                             "--save-gnuplot) and stb_fatbands_report.txt/references.bib into "
                             "(default: current directory). Created if it doesn't exist.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")

    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write fatbands_<category>.dat (one file per category, "
                             "k_position/energy(eV)/value columns) and fatbands.gplot together. "
                             "Off by default -- this tool used to write both unconditionally on "
                             "every run; that's no longer the case.")

    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of the fatbands plot before "
                             "finishing. Off by default -- this tool used to always show it, "
                             "with no way to skip.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-fatbands {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is set to 'manual'.")

    if args.label:
        if args.input_file or args.wfsx or args.hsx_file or args.geometry_file:
            parser.error("--label cannot be combined with --file/--wfsx/--hsx-file/--geometry-file.")
        args.input_file = f"{args.label}.bands"
        args.wfsx = resolve_wfsx_path(args.label)
    elif not args.input_file:
        parser.error("one of --label or --file is required.")
    elif not args.wfsx:
        parser.error("--wfsx is required when using --file instead of --label.")

    if not os.path.isfile(args.input_file):
        parser.error(f"'{args.input_file}' not found.")
    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = f"'{args.label}.bands.WFSX' and '{args.label}.WFSX'" if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}). stb-fatbands needs the wavefunction "
                     "coefficients written for the band-structure k-path (SIESTA fdf option "
                     "'WFS.Write.For.Bands T').")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    shift_desc = args.shift if args.shift != "manual" else f"manual ({args.manual_value} eV)"

    print("\n" + color_text("FATBANDS: reading input data", 'bold'))
    print("-" * 60)

    print("[INFO] Reading band structure data ...")
    fermi_energy, high_sym, dic_bands, nspin = read_data(args.input_file)
    bands_keys = list(dic_bands.keys())

    print("[INFO] Resolving geometry/Hamiltonian ...")
    parent, parent_source, is_approx = load_parent(args.label, args.hsx_file, args.geometry_file)
    if parent is None:
        parser.error("No geometry/Hamiltonian could be resolved -- pass --hsx-file, "
                      "--geometry-file, or use --label with a <label>.HSX/.fdf next to it.")
    if is_approx:
        print(color_text(
            "[WARNING] No .HSX found -- orbital weights use an implicit-orthogonal-basis "
            "approximation (|c|^2), not SIESTA's true non-orthogonal overlap-weighted norm. "
            "Pass --hsx-file <label>.HSX for physically correct weights.", 'yellow'))
        accuracy_desc = "approximate (|c|^2, no .HSX)"
    else:
        print(f"[INFO] Using '{parent_source}' for overlap-aware orbital weights.")
        accuracy_desc = "accurate (overlap-weighted norm, real Hamiltonian)"

    geometry = _geometry_of(parent)

    print(f"[INFO] Reading wavefunction coefficients from '{args.wfsx}' ...")
    sizes, states_by_k = read_wfsx_states(args.wfsx, parent)

    if sizes.nk != len(bands_keys):
        parser.error(
            f"k-point count mismatch: '{args.wfsx}' has {sizes.nk} k-points but "
            f"'{args.input_file}' has {len(bands_keys)} -- they must be from the same "
            "calculation/k-path (see stb-fatbands --help for the k-order assumption this "
            "tool relies on)."
        )
    if sizes.no_u != geometry.no:
        parser.error(
            f"Orbital count mismatch: '{args.wfsx}' has {sizes.no_u} orbitals but the "
            f"loaded geometry ('{parent_source}') has {geometry.no} -- they must be from "
            "the same calculation."
        )

    # Sanity cross-check: sample first/middle/last k, compare WFSX-stored
    # eigenvalues against .bands' own -- a hard error (not a warning) on
    # mismatch, since a silently-mis-zipped plot would look plausible but
    # be wrong.
    n = len(bands_keys)
    sample_idx = sorted(set([0, n // 2, n - 1]))
    for i in sample_idx:
        k_key = bands_keys[i]
        band_eigs = dic_bands[k_key]  # shape (nspin, nbands)
        wfsx_block = states_by_k[i]
        for s in range(nspin):
            state = wfsx_block.get(s)
            if state is None:
                continue
            wfsx_eigs = np.asarray(state.eig)
            m = min(len(wfsx_eigs), band_eigs.shape[1])
            diff = np.abs(wfsx_eigs[:m] - band_eigs[s][:m])
            if diff.size and np.max(diff) > args.k_tol:
                parser.error(
                    f"Eigenvalue mismatch between '{args.wfsx}' and '{args.input_file}' at "
                    f"k-index {i} (spin {s}): max diff {np.max(diff):.6f} eV exceeds --k-tol "
                    f"{args.k_tol} eV -- these files likely don't correspond to the same "
                    "calculation/k-path order."
                )
    print("[INFO] .bands / .WFSX correspondence check passed.")

    print(f"[INFO] Categorizing {geometry.no} orbitals by --projection {args.projection} ...")
    categories, _atoms, _species, n_excluded, excluded_l = build_orbital_table(geometry, args.projection)
    categories = np.array(categories, dtype=object)

    all_categories = sorted({c for c in categories if c is not None}, key=_category_sort_key(args.projection))
    if args.category:
        missing = set(args.category) - set(all_categories)
        if missing:
            parser.error(f"--category value(s) not found for --projection {args.projection}: "
                         f"{sorted(missing)}. Available: {all_categories}")
        selected_categories = [c for c in all_categories if c in args.category]
    else:
        selected_categories = all_categories

    result = cbm_vbm(fermi_energy, dic_bands, nspin, args.gap_tol)
    vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k = result["combined"]
    if args.shift == "vbm":
        rshift = vbm
    elif args.shift == "cbm":
        rshift = cbm
    elif args.shift == "fermi":
        rshift = fermi_energy
    else:
        rshift = args.manual_value
    shifted_bands = shift_bands(dic_bands, rshift)

    # Spin-resolved output categories: verified live on a real spin-polarized
    # calculation (an isolated O atom, SIESTA-converged to its physical 2
    # Bohr-magneton triplet moment) that merging both spin channels into one
    # category -- the original behavior -- silently combines two physically
    # very different band sets (here, spin-up/spin-down CBM differed by
    # ~29 eV) into one indistinguishable scatter series/data file, with no
    # way to tell which point belongs to which spin. nspin=2 now gets each
    # category split into "<category>_up"/"<category>_down" (same convention
    # dos.py already uses for spin-resolved PDOS columns); nspin=1 is
    # unaffected (no suffix). Any other nspin (non-collinear/SOC) falls back
    # to a raw "_s{n}" suffix, same as dos.py -- not specially handled or
    # verified beyond that, matching this tool's own documented limitation.
    spin_suffixes = {1: [""], 2: ["_up", "_down"]}.get(
        nspin, [f"_s{s + 1}" for s in range(nspin)])
    output_categories = [f"{cat}{suf}" for cat in selected_categories for suf in spin_suffixes]
    category_of_output = {f"{cat}{suf}": cat for cat in selected_categories for suf in spin_suffixes}

    print("[INFO] Computing orbital-projected weights (this may take a moment) ...")
    masks = {cat: (categories == cat) for cat in selected_categories}
    weights_by_category = {out_cat: [] for out_cat in output_categories}

    for i, k_key in enumerate(bands_keys):
        wfsx_block = states_by_k[i]
        energies = shifted_bands[k_key]
        for s in range(nspin):
            state = wfsx_block.get(s)
            if state is None:
                continue
            # norm2(projection="hadamard") is a sesquilinear (basis-resolved
            # inner product) quantity -- real by construction, but sisl
            # returns it as a complex dtype (matching the complex
            # wavefunction coefficients) with a numerically ~0 imaginary
            # part. Take .real explicitly rather than letting float()
            # silently (and noisily -- ComplexWarning) drop it.
            weights = state.norm2(projection="hadamard").real  # (nwf, norb)
            nwf = weights.shape[0]
            nbands_avail = min(nwf, energies.shape[1])
            suf = spin_suffixes[s]
            for cat in selected_categories:
                cat_weight_per_band = weights[:nbands_avail, masks[cat]].sum(axis=1)
                out_cat = f"{cat}{suf}"
                for b in range(nbands_avail):
                    weights_by_category[out_cat].append(
                        (k_key, float(energies[s][b]), float(cat_weight_per_band[b]))
                    )

    # --- From here on: the numbered, save-able report -------------------
    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("\n===== STB-FATBANDS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time       : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label           : {args.label}" if args.label else "Label           : (explicit --file/--wfsx, no --label)", f_out)
    print_dual(f"Bands file      : {args.input_file}", f_out)
    print_dual(f"WFSX file       : {args.wfsx}", f_out)
    print_dual(f"Shift mode      : {shift_desc}", f_out)
    print_dual(f"Gap tolerance   : {args.gap_tol} eV", f_out)
    print_dual(f"k-tol (cross-check): {args.k_tol} eV", f_out)
    print_dual(f"Projection      : {args.projection}", f_out)
    print_dual(f"Category filter : {', '.join(args.category) if args.category else 'all'}", f_out)
    print_dual(f"Output dir      : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot    : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    nk = len(bands_keys)
    nbands = next(iter(dic_bands.values())).shape[1]

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Fermi energy", f"{fermi_energy:.6f} eV"], None),
        (["Spin channels", f"{nspin} ({'non-polarized' if nspin == 1 else 'polarized'})"], None),
        (["Bands x k-points", f"{nbands} x {nk}"], None),
        (["Orbitals (basis)", f"{geometry.no}"], None),
        (["Geometry/Hamiltonian source", parent_source], None),
        (["Weight accuracy", accuracy_desc], None),
    ], f_out)
    if is_approx:
        print_dual(color_text(
            "[WARNING] No .HSX found -- orbital weights use an implicit-orthogonal-basis "
            "approximation (|c|^2), not SIESTA's true non-orthogonal overlap-weighted norm. "
            "Pass --hsx-file <label>.HSX for physically correct weights.", 'yellow'), f_out)
    print_dual(color_text(
        f"[OK] .bands / .WFSX k-point and orbital-count correspondence verified "
        f"(sampled eigenvalue cross-check, tol={args.k_tol} eV).", 'green'), f_out)

    print_section("[2] BAND GAP ANALYSIS (k-path)", f_out)
    print_table(["Quantity", "Value"], [
        (["VBM", f"{vbm:.6f} eV (k = {vbm_k:.6f})"], None),
        (["CBM", f"{cbm:.6f} eV (k = {cbm_k:.6f})"], None),
        (["Indirect gap", f"{indirect_gap:.6f} eV (fundamental: CBM - VBM, any k)"], None),
        (["Direct gap", f"{direct_gap:.6f} eV (same-k minimum, k = {direct_k:.6f})"], None),
        (["Gap type", gap_type], 'yellow' if gap_type == "Metallic" else None),
    ], f_out)
    if nspin == 2:
        for s, spin_result in enumerate(result["spins"]):
            spin_name = "Spin up" if s == 0 else "Spin down"
            sv, sc_, svk, sck, sig, sgt, sdg, sdk = spin_result
            print_dual(f"{spin_name}:", f_out)
            print_table(["Quantity", "Value"], [
                (["VBM", f"{sv:.6f} eV (k = {svk:.6f})"], None),
                (["CBM", f"{sc_:.6f} eV (k = {sck:.6f})"], None),
                (["Indirect gap", f"{sig:.6f} eV"], None),
                (["Direct gap", f"{sdg:.6f} eV (k = {sdk:.6f})"], None),
                (["Gap type", sgt], 'yellow' if sgt == "Metallic" else None),
            ], f_out)
        half_metallic = result["half_metallic"]
        line = f"Half-metallic : {'Yes' if half_metallic else 'No'}"
        print_dual(color_text(line, 'green') if half_metallic else line, f_out)

    print_section("[3] ORBITAL PROJECTION", f_out)
    print_dual(f"Categories found : {len(all_categories)} ({', '.join(all_categories)})", f_out)
    if args.category:
        print_dual(f"Categories used  : {', '.join(selected_categories)}", f_out)
    if n_excluded:
        print_dual(color_text(
            f"[WARNING] {n_excluded} orbital(s) with l={excluded_l} (l > 3, i.e. g-orbitals or "
            f"beyond, are not supported) excluded from --projection {args.projection}.", 'yellow'), f_out)
    if nspin == 2:
        print_dual(f"Spin-resolved    : yes -- each category split into <category>_up/"
                   f"<category>_down ({len(output_categories)} series total).", f_out)

    category_rows = []
    for cat in output_categories:
        n_orb = int(np.sum(categories == category_of_output[cat]))
        vals = np.array([w for _, _, w in weights_by_category[cat]])
        category_rows.append((
            [cat, str(n_orb), f"{vals.mean():.4f}", f"{vals.max():.4f}"], None
        ))
    print_table(["Category", "Orbitals", "Mean weight", "Max weight"], category_rows, f_out)

    print_section("[4] WRITING OUTPUT FILES", f_out)
    written_data = []
    gplot_path = None
    if args.save_gnuplot:
        data_names = [f"fatbands_{cat}" for cat in output_categories]
        for cat, name in zip(output_categories, data_names):
            path = write_scalar_data(args.output_dir, name, weights_by_category[cat])
            written_data.append(path)
        gplot_path = write_scalar_gplot(args.output_dir, "fatbands.gplot", "fatbands.pdf",
                                        high_sym, data_names, ps_scale=GNUPLOT_POINT_SCALE)
        print_dual(color_text(
            f"[OK] Data written: {len(written_data)} file(s) "
            f"({', '.join(os.path.basename(p) for p in written_data)}).", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write "
                   "fatbands_<category>.dat + fatbands.gplot).", f_out)

    print_section("[5] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status          : OK", f_out)
    print_dual(f"Best (indirect) gap: {indirect_gap:.6f} eV ({gap_type}).", f_out)
    if nspin == 2 and result["half_metallic"]:
        print_dual(color_text(
            "Half-metallic character detected (one spin channel metallic, the other has a gap).",
            'green'), f_out)
    for path in written_data:
        print_dual(f"Data file       : {path}", f_out)
    if gplot_path:
        print_dual(f"Gnuplot script  : {gplot_path}", f_out)
    print_dual(f"References      : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report          : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        if len(output_categories) == 1:
            cat = output_categories[0]
            k_pos, energy, weight = zip(*weights_by_category[cat])
            plot_scalar_on_bands(high_sym, k_pos, energy, weight, f"{cat} weight", _is_gamma,
                                  size_base=MARKER_SIZE_BASE, size_scale=MARKER_SIZE_SCALE)
        else:
            plot_multi_series_on_bands(high_sym, weights_by_category, _is_gamma,
                                        size_base=MARKER_SIZE_BASE, size_scale=MARKER_SIZE_SCALE)


if __name__ == "__main__":
    main()
