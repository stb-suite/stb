#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import sys
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.calc_directives import build_optical_block, OPTICAL_DIRECTION_VECTORS, OPTICAL_OFFDIAG_PAIRS
from stb.core.adsorption_sites import SINGLE_POINT_BLOCK
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo

REPORT_FILE = "optical_stage1.txt"
CONFIG_EXTRA_FILE = "config_extra.fdf"

# See core/calc_directives.py's own docstring for the full diagonal-vs-
# -biaxial derivation -- shared with stb-opticalAnalysis (Stage 2), which
# performs the actual eps_ij reconstruction from the biaxial folders this
# Stage 1 writes.
_AXIS_VECTORS = OPTICAL_DIRECTION_VECTORS
_OFFDIAG_REQUIRES = OPTICAL_OFFDIAG_PAIRS


def detect_structure_vacuum(fdf_structure, vacuum_gap):
    """Detects (vacuum_axes, is_2d) from `fdf_structure`'s own lattice/
    coordinates -- same pattern as gqca.py's own detect_structure_vacuum,
    duplicated here (not extracted to core/) since this is currently the
    second, not third, consumer of this exact 4-line combination and the
    two callers' surrounding context differs enough (gqca.py re-detects
    per representative-cluster structure; here it's a single detection
    on the one input structure, reused informationally in the report and
    NOT used to change which/how many direction folders get written --
    see stb-optical --help for why no direction is skipped for 2D
    inputs).
    """
    positions = np.array([pos for _, pos in fdf_structure.atoms])
    is_cartesian = fdf_structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, fdf_structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, fdf_structure.lattice, vacuum_gap)
    return vacuum_axes, any(vacuum_axes)


def write_direction_folder(out_dir, fdf_structure, calc_text, extra_fdf_text, pp_path):
    """Writes structure.fdf (unchanged -- no new geometry, this workflow
    never touches atomic positions/lattice) + calc.fdf (the user's own
    --calc template, untouched except for the %include at its very top --
    see structure_io.prepend_include) + config_extra.fdf (the forced
    single-point directives plus this one direction's own Optical.* block,
    `extra_fdf_text`, varying per folder) + copied pseudos for one
    direction folder. Same config_extra.fdf convention as every other
    WORKFLOW_TOOLS prep stage in this suite (see CLAUDE.md's
    "config_extra.fdf" convention note) -- never edits the user's own
    template text directly, and SIESTA's fdf reader is first-occurrence
    -wins for a duplicate label, so %include'ing this ahead of the
    template correctly overrides any conflicting MD.*/Optical.* setting
    the template itself might already carry. Deliberately does NOT write/
    recompute a kgrid.MonkhorstPack block -- the ground-state SCF k-grid
    is left entirely to the user's own -c/--calc template, matching
    raman_modes.py's own precedent for this exact Optical.* use case (the
    sampling density that actually matters for the optical spectrum is
    Optical.Mesh, already exposed via --optical-mesh, not the SCF k-grid).
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
        f.write(extra_fdf_text)
    symbols = sorted({sym for sym, _ in fdf_structure.atoms})
    for sym in symbols:
        copy_pseudo(pp_path, sym, out_dir)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 2: writes one SIESTA folder per requested "
        "Cartesian/biaxial polarization direction, ready for a linear-response Optical.Calculate run.", 'bold')}
Takes an ALREADY-RELAXED structure (--file) and a calc.fdf template (--calc, default: calc.fdf in the
current directory), and for each requested direction (--directions, default xx yy zz) writes
'<output-dir>/dir_<axis>/' with the SAME structure (no new geometry -- unlike this suite's other
WORKFLOW_TOOLS prep stages, there's nothing to deform/search/enumerate here). Every forced directive
-- MD.TypeOfRun CG / MD.Steps 0 (single-point, no relaxation) plus that folder's own
Optical.Mesh/Optical.Vector/Optical.Broaden/Optical.NumberOfBands block (SIESTA's own linear
-response/RPA dielectric-function machinery -- see the SIESTA manual's Optical.* directives) -- is
written into a config_extra.fdf sidecar, %include'd at the very top of calc.fdf, never edited into
your own --calc template text directly (this suite's standard config_extra.fdf convention -- see
CLAUDE.md). Doesn't run SIESTA -- run each folder yourself, then use stb-opticalAnalysis (which
aggregates every direction folder present into one combined report/spectrum, so run it once against
the whole --output-dir root, not once per folder).

[DIAGONAL vs. BIAXIAL DIRECTIONS] 'xx'/'yy'/'zz' set Optical.Vector to the pure Cartesian unit axis
(e_x/e_y/e_z) -- the diagonal dielectric-tensor component eps_ii comes directly from that one folder.
'xy'/'xz'/'yz' set Optical.Vector to the NORMALIZED BISECTOR of the two axes involved (e.g. xy ->
(e_x+e_y)/sqrt(2)) -- SIESTA's Optical.Vector always computes a single-direction response n^T.eps.n,
never a genuine off-diagonal component directly, so the bisector is the standard trick to reach one:
for a symmetric dielectric tensor, response along that bisector expands to
(eps_ii + eps_jj)/2 + eps_ij, so once eps_ii, eps_jj AND the bisector's own response are all known,
  eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2
recovers the true off-diagonal component. Requesting 'xy' without also requesting 'xx'/'yy' in the
same (or an earlier) run still writes a valid folder -- it's simply not reconstructible into eps_xy
until the 2 diagonal folders it needs also exist; the Stage 1 report flags this per direction.

Stage 2 (stb-opticalAnalysis) auto-detects each folder's direction from its own Optical.Vector value
(not the folder name) and, whenever a biaxial direction's 2 matching diagonal folders are ALSO
present in the same run, applies the reconstruction formula above automatically -- eps_ij(E) is
written to its own dedicated output alongside the raw per-direction eps_ii(E)/eps_jj(E) results, no
extra flag needed. A biaxial folder without its matching diagonal pair is still read normally by
Stage 2 (its raw eps_(bisector)(E) response is real data), just not reconstructible into eps_ij yet.

[KNOWN LIMITATION, 2D/slab inputs] No direction is skipped or specially treated for a vacuum-padded
input -- every requested direction (including out-of-plane and any biaxial one involving it) is still
written and computed. For a slab, the periodic-supercell dielectric response along ANY direction is
diluted by the vacuum region (a well-known methodological caveat for slab dielectric/polarizability
DFT calculations) -- no vacuum-thickness rescaling is applied here. Interpret out-of-plane 2D results
with this in mind.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f relaxed.fdf -p dojo\n"
               "  %(prog)s -f relaxed.fdf -p dojo --directions xx yy zz xy xz yz\n"
    )

    parser.add_argument("-f", "--file", type=str, required=True,
                         help="Already-relaxed input structure.fdf.")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                         help="calc.fdf template for the Optical calculations (default: calc.fdf).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--directions", type=str, nargs='+', default=["xx", "yy", "zz"],
                         choices=["xx", "yy", "zz", "xy", "xz", "yz"],
                         help="Direction(s) to compute (default: xx yy zz). 'xx'/'yy'/'zz' are the "
                              "diagonal Cartesian components (Optical.Vector = a unit axis, "
                              "eps_ii read directly). 'xy'/'xz'/'yz' are BIAXIAL directions "
                              "(Optical.Vector = the normalized bisector of the 2 axes, e.g. "
                              "(1,1,0)/sqrt(2) for xy) -- the off-diagonal eps_ij component is only "
                              "reconstructible from these once the matching xx/yy/zz pair also "
                              "exists; see the description above for the exact formula.")
    parser.add_argument("--optical-mesh", type=int, nargs=3, default=[10, 10, 10], metavar=("NX", "NY", "NZ"),
                         help="Optical.Mesh k-point sampling for the optical-transition integral "
                              "(default: 10 10 10).")
    parser.add_argument("--optical-broaden", type=float, default=0.2, metavar="EV",
                         help="Optical.Broaden Gaussian broadening in eV (default: 0.2).")
    parser.add_argument("--optical-nbands", type=int, default=None, metavar="N",
                         help="Optical.NumberOfBands (default: SIESTA's own default -- all bands).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang, for the report's "
                              "dimensionality note only (default: 10.0).")
    parser.add_argument("-O", "--output-dir", type=str, default="optical_study",
                         help="Root directory for the whole optical-properties workflow "
                              "(default: optical_study).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-optical {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OPTICAL PROPERTIES WORKFLOW -- STAGE 1: DIRECTION FOLDERS", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.file):
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
        sys.exit(1)
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    structure = structure_io.read_fdf(args.file)
    vacuum_axes, is_2d = detect_structure_vacuum(structure, args.vacuum_gap)

    with open(args.calc) as f:
        calc_template = f.read()

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None

    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== OPTICAL STAGE 1 REPORT (DIRECTION FOLDERS) =====', 'magenta')}", f_out)

    species = structure_io.species_list(structure)
    counts = structure_io.atom_counts(structure)
    formula = " ".join(f"{sym}{counts[sym]}" for sym in species)
    lattice = structure_io.lattice_only(structure)
    cell_volume = float(abs(np.linalg.det(lattice)))
    n_atoms = len(structure.atoms)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Structure       : {args.file}", f_out)
    print_dual(f"Formula         : {formula}  ({n_atoms} atom(s), {len(species)} species: {', '.join(species)})", f_out)
    print_dual(f"Cell volume     : {cell_volume:.4f} Ang^3", f_out)
    print_dual(f"Calc template   : {args.calc}", f_out)
    print_dual(f"Directions      : {' '.join(args.directions)}", f_out)
    print_dual(f"Pseudo source   : {args.pseudo_dir if args.pseudo_dir else '(none -- pseudos not copied)'}", f_out)
    print_dual(f"Output root     : {output_root}", f_out)

    print_section('[1] DIMENSIONALITY', f_out)
    print_dual(f"Detected : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    if is_2d:
        print_dual(color_text(
            "[KNOWN LIMITATION] Vacuum-padded (2D/slab) input -- the periodic-supercell "
            "dielectric response along ANY direction (including out-of-plane) is diluted by "
            "the vacuum region, a known methodological caveat for slab dielectric/"
            "polarizability DFT calculations. No vacuum-thickness rescaling is applied. Every "
            "requested direction is still computed (no direction is skipped for 2D inputs).",
            'yellow'), f_out)

    print_section('[2] OPTICAL PARAMETERS (applied identically to every direction)', f_out)
    print_dual(f"Optical.Mesh          : {args.optical_mesh[0]} {args.optical_mesh[1]} {args.optical_mesh[2]}"
                "  (k-point sampling for the optical-transition integral -- independent of the "
                "ground-state SCF k-grid, which comes entirely from your own --calc template)", f_out)
    print_dual(f"Optical.Broaden       : {args.optical_broaden} eV  (Gaussian broadening applied to "
                "every interband transition before summing into eps2(E))", f_out)
    print_dual(f"Optical.NumberOfBands : "
                f"{args.optical_nbands if args.optical_nbands else '(SIESTA default -- all bands)'}", f_out)
    print_dual(f"Single-point forced   : MD.TypeOfRun CG / MD.Steps 0 (no relaxation -- the input "
                "geometry is evaluated exactly as given in every direction folder)", f_out)
    print_dual(f"\nEvery direction folder gets its own {CONFIG_EXTRA_FILE} (%include'd at the top "
                "of calc.fdf, never edited into your own --calc template -- SIESTA's fdf reader is "
                "first-occurrence-wins for a duplicate label, so this correctly overrides any "
                "conflicting MD.*/Optical.* setting the template itself carries). The single-point "
                "part, shared by every direction, is:", f_out)
    for line in SINGLE_POINT_BLOCK.rstrip("\n").split("\n"):
        print_dual(f"  {line}", f_out)
    print_dual("followed by that folder's own Optical.* block (varies per direction -- see [3] "
                "below for each one's exact Optical.Vector).", f_out)

    print_section('[3] DIRECTIONS & DIELECTRIC-TENSOR MAPPING', f_out)
    print_dual(f"{'Direction':<10} {'Kind':<12} {'Optical.Vector':<24} {'Gives access to':<28} {'Needs also':<12}", f_out)
    print_dual("-" * 60, f_out)
    for axis in args.directions:
        axis_vec = _AXIS_VECTORS[axis]
        vec_str = f"{axis_vec[0]:.4f} {axis_vec[1]:.4f} {axis_vec[2]:.4f}"
        if axis in _OFFDIAG_REQUIRES:
            need_a, need_b = _OFFDIAG_REQUIRES[axis]
            have = (need_a in args.directions, need_b in args.directions)
            needs_str = f"{need_a}+{need_b}"
            if not all(have):
                missing = [n for n, ok in zip((need_a, need_b), have) if not ok]
                needs_str += color_text(f" (missing {'/'.join(missing)} in this run)", 'yellow')
            print_dual(f"{axis:<10} {'biaxial':<12} {vec_str:<24} "
                        f"{'eps_' + axis + ' (needs reconstruction)':<28} {needs_str}", f_out)
        else:
            print_dual(f"{axis:<10} {'diagonal':<12} {vec_str:<24} {'eps_' + axis + ' (direct)':<28} {'--':<12}", f_out)
    if any(axis in _OFFDIAG_REQUIRES for axis in args.directions):
        print_dual(color_text(
            "\nReconstruction formula for any biaxial direction above: "
            "eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2 -- stb-opticalAnalysis applies this "
            "automatically once all 3 folders involved (ii, jj, ij) are present and done; no extra "
            "flag needed.", 'cyan'), f_out)

    print_section('[4] FOLDERS WRITTEN', f_out)

    forced_calc_text = structure_io.prepend_include(calc_template, CONFIG_EXTRA_FILE)
    for axis in args.directions:
        axis_vec = _AXIS_VECTORS[axis]
        optical_block = build_optical_block(args.optical_mesh, args.optical_broaden, axis_vec,
                                             args.optical_nbands)
        extra_fdf_text = SINGLE_POINT_BLOCK + "\n" + optical_block
        out_dir = os.path.join(output_root, f"dir_{axis}")
        write_direction_folder(out_dir, structure, forced_calc_text, extra_fdf_text, args.pseudo_dir)
        kind = "biaxial" if axis in _OFFDIAG_REQUIRES else "diagonal"
        print_dual(f"  {color_text('[OK]', 'green')} {out_dir}  ({kind}, Optical.Vector = "
                    f"{axis_vec[0]:.4f} {axis_vec[1]:.4f} {axis_vec[2]:.4f}, {n_atoms} atoms, "
                    f"species: {', '.join(species)})", f_out)

    print_section('[5] SUMMARY & NEXT STEPS', f_out)
    n_diag = sum(1 for a in args.directions if a not in _OFFDIAG_REQUIRES)
    n_off = len(args.directions) - n_diag
    print_dual(f"{len(args.directions)} direction folder(s) written under '{output_root}' "
                f"({n_diag} diagonal, {n_off} biaxial).", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{output_root}/dir_*/' folder.", f_out)
    print_dual(f"  2. Once they're done, run: stb-opticalAnalysis --directory {output_root}", f_out)
    if n_off:
        print_dual("  3. stb-opticalAnalysis reconstructs eps_ij(E) automatically for any biaxial "
                    "direction whose matching diagonal pair is also present -- see Section [3] above.",
                    f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Direction folders ready for Stage 2 (stb-opticalAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
