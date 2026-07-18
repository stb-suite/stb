#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.calc_directives import force_single_point, build_optical_block
from stb.core.pseudopotentials import resolve_pseudo_source, link_pseudo

REPORT_FILE = "optical_stage1.txt"

_AXIS_VECTORS = {
    "x": (1.0, 0.0, 0.0),
    "y": (0.0, 1.0, 0.0),
    "z": (0.0, 0.0, 1.0),
}


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


def write_direction_folder(out_dir, fdf_structure, calc_text, pp_path):
    """Writes structure.fdf (unchanged -- no new geometry, this workflow
    never touches atomic positions/lattice) + calc.fdf (single-point +
    the Optical.* block for one Cartesian direction, already baked into
    calc_text by the caller) + linked pseudos for one direction folder.
    Deliberately does NOT write/recompute a kgrid.MonkhorstPack block --
    the ground-state SCF k-grid is left entirely to the user's own -c/
    --calc template, matching raman_modes.py's own precedent for this
    exact Optical.* use case (the sampling density that actually matters
    for the optical spectrum is Optical.Mesh, already exposed via
    --optical-mesh, not the SCF k-grid).
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = sorted({sym for sym, _ in fdf_structure.atoms})
    for sym in symbols:
        link_pseudo(pp_path, sym, out_dir)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 2: writes one SIESTA folder per requested "
        "Cartesian polarization direction, ready for a linear-response Optical.Calculate run.", 'bold')}
Takes an ALREADY-RELAXED structure (--file) and a calc.fdf template (--calc), and for each requested
direction (--directions, default x y z) writes '<output-dir>/dir_<axis>/' with the SAME structure
(no new geometry -- unlike this suite's other WORKFLOW_TOOLS prep stages, there's nothing to deform/
search/enumerate here) plus a forced single-point calc.fdf carrying the Optical.Mesh/Optical.Vector/
Optical.Broaden block for that one axis (SIESTA's own linear-response/RPA dielectric-function
machinery -- see the SIESTA manual's Optical.* directives). Doesn't run SIESTA -- run each folder
yourself, then use stb-opticalAnalysis (which aggregates every direction folder present into one
combined report/spectrum, so run it once against the whole --output-dir root, not once per folder).

[KNOWN LIMITATION, 2D/slab inputs] No direction is skipped or specially treated for a vacuum-padded
input -- every requested direction (including out-of-plane) is still written and computed. For a
slab, the periodic-supercell dielectric response along ANY direction is diluted by the vacuum region
(a well-known methodological caveat for slab dielectric/polarizability DFT calculations) -- no
vacuum-thickness rescaling is applied here. Interpret out-of-plane 2D results with this in mind.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s -f relaxed.fdf -c calc.fdf -p dojo --directions x y z\n"
    )

    parser.add_argument("-f", "--file", type=str, required=True,
                         help="Already-relaxed input structure.fdf.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template for the Optical calculations.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--directions", type=str, nargs='+', default=["x", "y", "z"],
                         choices=["x", "y", "z"],
                         help="Cartesian Optical.Vector direction(s) to compute (default: x y z). "
                              "Only the diagonal response along these axes is computed -- no off-"
                              "diagonal dielectric-tensor components, no point-group-based direction "
                              "reduction in this version.")
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
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, 'w') as f_out:
        print_dual(f"{color_text('===== OPTICAL STAGE 1 REPORT (DIRECTION FOLDERS) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Structure       : {args.file}", f_out)
        print_dual(f"Calc template   : {args.calc}", f_out)
        print_dual(f"Directions      : {' '.join(args.directions)}", f_out)
        print_dual(f"Output root     : {output_root}", f_out)

        print_dual(f"\n{color_text('[1] DIMENSIONALITY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Detected : {kspace.dimensionality_label(vacuum_axes)}", f_out)
        if is_2d:
            print_dual(color_text(
                "[KNOWN LIMITATION] Vacuum-padded (2D/slab) input -- the periodic-supercell "
                "dielectric response along ANY direction (including out-of-plane) is diluted by "
                "the vacuum region, a known methodological caveat for slab dielectric/"
                "polarizability DFT calculations. No vacuum-thickness rescaling is applied. Every "
                "requested direction is still computed (no direction is skipped for 2D inputs).",
                'yellow'), f_out)

        print_dual(f"\n{color_text('[2] OPTICAL PARAMETERS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Optical.Mesh     : {args.optical_mesh[0]} {args.optical_mesh[1]} {args.optical_mesh[2]}", f_out)
        print_dual(f"Optical.Broaden  : {args.optical_broaden} eV", f_out)
        print_dual(f"Optical.NumberOfBands : {args.optical_nbands if args.optical_nbands else '(SIESTA default)'}", f_out)

        print_dual(f"\n{color_text('[3] FOLDERS WRITTEN', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        calc_text_base = force_single_point(calc_template)
        for axis in args.directions:
            axis_vec = _AXIS_VECTORS[axis]
            optical_block = build_optical_block(args.optical_mesh, args.optical_broaden, axis_vec,
                                                 args.optical_nbands)
            calc_text = calc_text_base + "\n" + optical_block
            out_dir = os.path.join(output_root, f"dir_{axis}")
            write_direction_folder(out_dir, structure, calc_text, args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} {out_dir} (Optical.Vector = "
                        f"{axis_vec[0]:.0f} {axis_vec[1]:.0f} {axis_vec[2]:.0f})", f_out)

        print_dual(f"\n{color_text('[4] SUMMARY & NEXT STEPS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"{len(args.directions)} direction folder(s) written under '{output_root}'.", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual(f"  1. Run SIESTA in every '{output_root}/dir_*/' folder.", f_out)
        print_dual(f"  2. Once they're done, run: stb-opticalAnalysis --directory {output_root}", f_out)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Direction folders ready for Stage 2 (stb-opticalAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
