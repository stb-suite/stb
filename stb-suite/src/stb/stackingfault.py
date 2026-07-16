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
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.pseudopotentials import resolve_pseudo_source, link_pseudo
from stb.core.calc_directives import force_single_point
from stb.core.heterostructure import find_zsl_match, build_stacked_structure

REPORT_FILE = "stackingfault_setup.txt"


def read_layer(path):
    """Reads one monolayer .fdf into a pymatgen Structure, exiting with a
    clear error on any parse failure -- same "fail loud, don't propagate a
    cryptic traceback" convention as adsorb.py/neb.py's own file reading.
    """
    try:
        return structure_io.to_pymatgen(structure_io.read_fdf(path))
    except Exception as e:
        print(color_text(f"[ERROR] Could not read '{path}': {e}", 'red'))
        sys.exit(1)


def min_interlayer_distance(pmg_structure, n_layer1):
    """Minimum periodic distance between any layer-1 atom (index <
    n_layer1) and any layer-2 atom (index >= n_layer1) -- same style as
    adsorb.py's min_adsorbate_slab_distance, generalized from
    adsorbate-vs-slab to layer-vs-layer. Purely informational here (unlike
    adsorb.py's overlap WARNING): a close contact at some grid point is
    often exactly the physics a stacking-fault sweep is meant to sample
    (an eclipsed, high-energy registry), not a mistake to flag.
    """
    n_total = len(pmg_structure)
    if n_layer1 == 0 or n_layer1 >= n_total:
        return None
    dm = pmg_structure.distance_matrix
    return float(dm[:n_layer1, n_layer1:].min())


def write_grid_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path):
    """Writes structure.fdf + calc.fdf + linked pseudopotentials for one
    shift_II_JJ/ grid-point folder -- same shape as adsorb.py's
    write_reference_folder / neb.py's write_image_folder. `species_meta`
    only needs to cover layer1's declared species -- any symbol unique to
    layer2 gets a fresh id automatically via structure_io.from_pymatgen's
    ensure_species_id, same as adsorb.py passing the slab's species_meta
    alone when the adsorbate introduces a new element.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta,
                                                coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = {site.specie.symbol for site in pmg_structure}
    for sym in sorted(symbols):
        link_pseudo(pp_path, sym, out_dir)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Prepares SIESTA folders for a 2D stacking-fault (generalized "
        "stacking fault energy / gamma-surface) study: rigidly slides one layer of a bilayer "
        "across a 2D grid of lateral offsets and writes one single-point SIESTA folder per grid "
        "point.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -n 7\n"
               "  %(prog)s -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf -n 9 -g 3.4\n"
               "\n"
               "Pass the SAME file for -l1/-l2 for the canonical use case (a material sliding "
               "against itself, e.g. graphite ABA vs. ABC stacking) -- already the precedented "
               "'identical layers' case for stb-2Dstacking's own ZSL match. Different files "
               "study an interlayer/heterostructure sliding energy landscape instead.\n"
    )

    parser.add_argument("-l1", "--layer1", type=str, required=True,
                         help="Bottom monolayer .fdf (stays fixed).")
    parser.add_argument("-l2", "--layer2", type=str, required=True,
                         help="Top monolayer .fdf (rigidly slid across the grid). Pass the same "
                              "file as --layer1 for a material sliding against itself.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template (kgrid, basis, XC, %%include structure.fdf, "
                              "etc.) -- copied into every shift_II_JJ/ folder with MD.TypeOfRun/"
                              "MD.NumCGsteps forced to a single-point evaluation.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (optional): a bundled bank or a folder path.")

    parser.add_argument("-a", "--max_area", type=float, default=150.0,
                         help="Maximum ZSL commensurate-supercell area, Ang^2 (default: 150.0, "
                              "same as stb-2Dstacking).")
    parser.add_argument("-s", "--max_strain", type=float, default=0.05,
                         help="Maximum allowed ZSL match strain fraction (default: 0.05).")
    parser.add_argument("-id", "--match_id", type=int, default=0,
                         help="Which ZSL match to use, 0-based, best (lowest strain) first "
                              "(default: 0).")
    parser.add_argument("-g", "--gap", type=float, default=3.2,
                         help="Interlayer gap in Ang, FIXED across the whole grid (default: 3.2, "
                              "same as stb-2Dstacking) -- a rigid-shift protocol, standard in the "
                              "literature; per-point interlayer-distance relaxation is not done "
                              "here (some high-energy registries may end up with a non-relaxed "
                              "distance -- see the [INFO] closest-contact line in the report).")
    parser.add_argument("--vacuum", type=float, default=None,
                         help="Target vacuum space in Ang. Inherits --layer1's by default.")
    parser.add_argument("-sm", "--strain_mode", choices=["top", "bottom", "sym"], default="top",
                         help="Strain distribution mode for the ZSL match (default: top).")
    parser.add_argument("-t", "--twist", type=float, default=0.0,
                         help="Twist angle of layer 2 in degrees (default: 0.0), FIXED for the "
                              "whole grid -- a different twist is a different physical system, "
                              "not swept here.")

    parser.add_argument("-n", "--grid-n", type=int, default=7,
                         help="Grid resolution: an N x N grid of lateral shifts covering [0, 1) "
                              "in each direction (default: 7).")

    parser.add_argument("-O", "--output-dir", type=str, default=".",
                         help="Root directory (default: current directory) for every "
                              "shift_II_JJ/.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-stackingfault {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.grid_n < 2:
        parser.error("-n/--grid-n must be >= 2 (need at least 2 points per axis for a sweep).")

    print("\n" + color_text("Prepare a 2D stacking-fault study:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.layer1):
        print(color_text(f"[ERROR] Layer 1 file '{args.layer1}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.layer2):
        print(color_text(f"[ERROR] Layer 2 file '{args.layer2}' not found.", 'red'))
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

    layer1_structure = structure_io.read_fdf(args.layer1)
    species_meta = structure_io.species_dict(layer1_structure)
    layer1_pmg = read_layer(args.layer1)
    layer2_pmg = read_layer(args.layer2)

    with open(args.calc) as f:
        calc_text = f.read()
    single_point_calc_text = force_single_point(calc_text)

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, 'w') as f_out:
        print_dual(f"\n{color_text('[0] RUN METADATA', 'bold')}", f_out)
        print_dual(f"  Layer 1         : {args.layer1}", f_out)
        print_dual(f"  Layer 2         : {args.layer2}", f_out)
        print_dual(f"  Calc template   : {args.calc}", f_out)
        print_dual(f"  Pseudo dir      : {args.pseudo_dir or '(none)'}", f_out)
        print_dual(f"  Output dir      : {output_root}", f_out)
        print_dual(f"  Grid            : {args.grid_n} x {args.grid_n}", f_out)
        print_dual(f"  Gap (fixed)     : {args.gap} Ang", f_out)
        print_dual(f"  Twist (fixed)   : {args.twist} deg", f_out)
        print_dual(f"  max_area/max_strain/match_id: {args.max_area} / {args.max_strain} / {args.match_id}", f_out)

        print_dual(f"\n{color_text('[1] ZSL MATCH', 'bold')}", f_out)
        t_mat1, t_mat2, best_match_data = find_zsl_match(
            layer1_pmg, layer2_pmg, max_area=args.max_area, max_strain=args.max_strain,
            match_id=args.match_id, interactive=False, twist_angle=args.twist)
        print_dual(f"  Selected match ID {args.match_id}: area {best_match_data['area']:.2f} Ang^2, "
                    f"strain {best_match_data['strain']:.2f}%, angular strain "
                    f"{best_match_data['angle_strain']:.2f} deg.", f_out)
        n_layer1_supercell = len(layer1_pmg) * round(np.linalg.det(t_mat1))
        n_layer2_supercell = len(layer2_pmg) * round(np.linalg.det(t_mat2))
        n_total = n_layer1_supercell + n_layer2_supercell
        n_grid_points = args.grid_n * args.grid_n
        print_dual(f"  Each grid point: {n_total} atoms ({n_layer1_supercell} layer 1 + "
                    f"{n_layer2_supercell} layer 2). {n_grid_points} grid point(s) total -- "
                    f"{n_grid_points} independent single-point SIESTA runs at {n_total} atoms each.",
                    f_out)

        shifts = list(np.linspace(0.0, 1.0, args.grid_n, endpoint=False))

        print_dual(f"\n{color_text('[2] GRID FOLDERS', 'bold')}", f_out)
        report_rows = []  # (label, i, j, shift_x, shift_y, dir)
        closest_contact = None
        closest_label = None
        for i, shift_x in enumerate(shifts):
            for j, shift_y in enumerate(shifts):
                label = f"shift_{i:02d}_{j:02d}"
                hetero, n_layer1_atoms, _max_strain_val = build_stacked_structure(
                    layer1_pmg, layer2_pmg, t_mat1, t_mat2, shift_x, shift_y, args.gap,
                    target_vacuum=args.vacuum, strain_mode=args.strain_mode)

                contact = min_interlayer_distance(hetero, n_layer1_atoms)
                if contact is not None and (closest_contact is None or contact < closest_contact):
                    closest_contact = contact
                    closest_label = label

                grid_dir = os.path.join(output_root, label)
                write_grid_folder(grid_dir, hetero, single_point_calc_text, species_meta,
                                   args.pseudo_dir)
                print_dual(f"  {color_text('[OK]', 'green')} {grid_dir} "
                            f"(shift: {shift_x:.4f}, {shift_y:.4f})", f_out)
                report_rows.append((label, i, j, shift_x, shift_y, grid_dir))

        print_dual(f"\n{color_text('[3] SUMMARY', 'bold')}", f_out)
        print_dual(f"  {len(report_rows)} grid folder(s) written under '{output_root}'.", f_out)
        if closest_contact is not None:
            print_dual(f"  [INFO] Closest interlayer contact anywhere in the grid: "
                        f"{closest_contact:.3f} Ang at {closest_label} -- expected to be small at "
                        "eclipsed/high-energy registries, not necessarily a problem.", f_out)
        print_dual("  Run SIESTA (single-point: MD.NumCGsteps forced to 0) in every "
                    "shift_II_JJ/ folder, then use stb-stackingfaultAnalysis.", f_out)

        f_out.write("\n# GRID_TABLE -- parsed by stb-stackingfaultAnalysis, do not reorder the columns\n")
        f_out.write(f"# {'label':<16}{'i':<5}{'j':<5}{'shift_x':<12}{'shift_y':<12}{'dir'}\n")
        for label, i, j, shift_x, shift_y, grid_dir in report_rows:
            f_out.write(f"{label:<18}{i:<5}{j:<5}{shift_x:<12.6f}{shift_y:<12.6f}{grid_dir}\n")

    print(f"\n{color_text('Success:', 'green')} {len(report_rows)} grid folder(s) written under "
          f"'{output_root}'.")
    print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
