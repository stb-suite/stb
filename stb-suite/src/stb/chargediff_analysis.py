#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import os
import sys
import glob
import json
import argparse
from datetime import datetime

import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.rho_io import read_total, try_read_net_spin, report_quantity, find_one_rho

REPORT_FILE = "chargediff_analysis_report.txt"
BIB_FILE = "references.bib"
CHARGEDIFF_MANIFEST_FILE = "chargediff_manifest.json"


def find_geometry_file(folder):
    """Returns a (path, sisl Geometry) for --cube's geometry, read from
    whichever finished '*.XV' is in `folder` (preferred, the actual
    converged geometry) or its 'structure.fdf' (written by
    stb-chargediffPrep) otherwise. Returns (None, None) if neither exists.
    """
    xv_matches = sorted(glob.glob(os.path.join(folder, "*.XV")))
    if xv_matches:
        path = xv_matches[0]
    else:
        path = os.path.join(folder, "structure.fdf")
        if not os.path.isfile(path):
            return None, None
    return path, sisl.get_sile(path).read_geometry()


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of the Charge Density Difference workflow: given just "
        "--dir (the same folder stb-chargediffPrep wrote into), auto-discovers the combined system "
        "and every fragment from chargediff_manifest.json, reads each one's own .RHO (once SIESTA "
        "has been run in every folder), and computes Delta rho = rho_combined - "
        "sum(rho_fragment_i) -- the standard way to visualize charge accumulation/depletion "
        "regions.", 'bold')}
Produces the same 2D slice / 3D point cloud / planar-averaged 1D profile /
Gaussian cube output as stb-density, so the visualization conventions match
the rest of the suite. --spin computes the spin-density difference instead
(only available if the combined system and every fragment are
spin-polarized).""",
        epilog="Example usage:\n"
               "  %(prog)s\n"
               "  %(prog)s --dir chargediff_run --profile --axis 2\n"
               "  %(prog)s --dir chargediff_run --view\n"
               "  %(prog)s --dir chargediff_run --no-cube\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--dir", default="chargediff_run",
                         help="Directory stb-chargediffPrep wrote into -- must contain "
                              "'chargediff_manifest.json' plus the 'combined/' and fragment "
                              "subfolders it lists, each already run in SIESTA. Default: "
                              "chargediff_run (stb-chargediffPrep's own default --output-dir).")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write the data/gnuplot/cube/report files and "
                             "references.bib into (default: current directory). Created if "
                             "it doesn't exist.")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--3d", action="store_true", dest="full_3d", help="Export full 3D volume")
    mode_group.add_argument("--profile", action="store_true",
                        help="Write the planar average of Delta rho along --axis as a 1D profile "
                             "(average taken over the other two axes) instead of a 2D/3D map -- "
                             "useful for bilayers, slabs, and interfaces")
    parser.add_argument("-a", "--axis", type=int, default=2, choices=[0,1,2],
                        help="For a slice/profile: the axis normal to the cut plane (slice) or "
                             "the axis the profile varies along (profile) (0=X, 1=Y, 2=Z)")
    parser.add_argument("-p", "--pos", type=float, help="Position (Angstrom) of the 2D cut (ignored in --profile mode)")
    parser.add_argument("--spin", action="store_true",
                        help="Compute the SPIN density difference instead of the charge density "
                             "difference -- only available if the combined system AND every "
                             "fragment are spin-polarized.")
    parser.add_argument("--iso-min", type=float, default=None, metavar="RHO",
                        help="For --3d: only export/preview points with |Delta rho| >= this "
                             "threshold (e/Ang^3) -- keeps the point cloud a manageable size "
                             "and visually meaningful instead of dominated by near-zero points")
    parser.add_argument("--no-cube", dest="cube", action="store_false", default=True,
                        help="Skip writing Delta rho's full 3D grid as a Gaussian .cube file -- "
                             "written by default (needs a '*.XV' or 'structure.fdf' in the "
                             "combined folder for the geometry), for real 3D isosurface "
                             "rendering in VESTA/VMD/Avogadro.")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Fix the colorbar/palette lower bound manually (e/Ang^3), e.g. to "
                             "compare several slices on the same scale.")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Fix the colorbar/palette upper bound manually (e/Ang^3).")
    parser.add_argument("--contour", action="store_true",
                        help="Overlay contour lines on a 2D slice map/preview (ignored for --3d/--profile)")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a .gplot script next to each .dat file this run "
                             "generates. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of Delta rho (slice heatmap "
                             "/ profile line / 3D scatter, matching the mode). Off by default.")

    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-chargediffAnalysis {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Charge Density Difference (Analysis)",
            "Delta rho = rho_combined - sum(rho_fragment_i)",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CHARGE DENSITY DIFFERENCE -- STAGE 2: ANALYSIS:", 'bold'))
    print("-" * 60)

    manifest_path = os.path.join(args.dir, CHARGEDIFF_MANIFEST_FILE)
    if not os.path.isfile(manifest_path):
        print(color_text(
            f"[ERROR] '{manifest_path}' not found -- run stb-chargediffPrep with "
            f"--output-dir {args.dir!r} first (this file is where it records the combined/"
            "fragment folder names Stage 2 needs; it's never typed by hand).", 'red'))
        sys.exit(1)
    with open(manifest_path) as f:
        manifest = json.load(f)
    if "combined_folder" not in manifest or "fragment_folders" not in manifest:
        print(color_text(
            f"[ERROR] '{manifest_path}' is from an older stb-chargediffPrep (schema "
            f"{manifest.get('schema', 'unknown')!r}) that didn't yet write a 'combined/' folder "
            "or record 'combined_folder'/'fragment_folders' -- re-run stb-chargediffPrep with "
            f"--output-dir {args.dir!r} to regenerate it before using this Stage 2.", 'red'))
        sys.exit(1)
    combined_folder = manifest["combined_folder"]
    fragment_folders = manifest["fragment_folders"]

    folders = {"combined": os.path.join(args.dir, combined_folder)}
    for name in fragment_folders:
        folders[name] = os.path.join(args.dir, name)

    os.makedirs(args.output_dir, exist_ok=True)

    filenames = {}
    try:
        for name, folder in folders.items():
            filenames[name] = find_one_rho(folder)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    grids = {}
    for name, filename in filenames.items():
        print(f"[INFO] Reading: {color_text(filename, 'cyan')}")
        try:
            sile = sisl.get_sile(filename)
            grid_obj = try_read_net_spin(sile, filename, required=True) if args.spin else read_total(sile)
        except Exception as e:
            print(color_text(f"[FATAL] {filename}: {e}", 'red'))
            sys.exit(1)
        grids[name] = grid_obj

    combined_grid_obj = grids["combined"]
    lattice = combined_grid_obj.lattice.cell
    origin = combined_grid_obj.origin
    combined_shape = combined_grid_obj.grid.shape
    print(f"[INFO] Combined grid: {combined_shape[0]}x{combined_shape[1]}x{combined_shape[2]}")

    for name in fragment_folders:
        frag_shape = grids[name].grid.shape
        if frag_shape != combined_shape:
            print(color_text(
                f"[ERROR] Grid shape mismatch: '{filenames['combined']}' is "
                f"{combined_shape} but '{filenames[name]}' is {frag_shape} -- every fragment "
                "must be on the same grid as the combined system to take a difference (see "
                "stb-chargediffPrep's --mesh-cutoff).", 'red'))
            sys.exit(1)
        if not np.allclose(grids[name].lattice.cell, lattice, atol=1e-3):
            print(color_text(
                f"[ERROR] Lattice mismatch: '{filenames[name]}' does not share the combined "
                "system's cell -- Delta rho is only meaningful point-by-point on an identical "
                "cell/grid.", 'red'))
            sys.exit(1)

    diff_data = combined_grid_obj.grid.copy()
    for name in fragment_folders:
        diff_data = diff_data - grids[name].grid

    quantity_label = "Delta Spin Density" if args.spin else "Delta Charge Density"

    if args.profile:
        mode = 'profile'
    elif args.full_3d:
        mode = '3d'
    else:
        mode = 'slice'

    if args.contour and mode != 'slice':
        print(color_text(
            "[WARNING] --contour only applies to the default 2D slice mode -- ignoring.", 'yellow'))
    if args.pos is not None and mode == 'profile':
        print(color_text(
            "[WARNING] --pos has no effect in --profile mode -- the profile covers the "
            "whole axis.", 'yellow'))

    stem = f"{os.path.basename(os.path.normpath(args.dir)) or 'chargediff'}_chargediff"

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-CHARGEDIFF ANALYSIS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Prep dir       : {args.dir}", f_out)
    print_dual(f"Combined folder: {folders['combined']} ({filenames['combined']})", f_out)
    for name in fragment_folders:
        print_dual(f"Fragment folder: {folders[name]} ({filenames[name]})", f_out)
    print_dual(f"Mode           : {mode}" + (f" (axis {args.axis})" if mode != '3d' else ""), f_out)
    print_dual(f"Quantity       : {quantity_label}", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)

    files_written = {}
    files_written["charge"] = report_quantity(
        f"[1] {quantity_label.upper()}", diff_data, lattice, origin, quantity_label, True,
        stem, mode, args, args.output_dir, f_out)
    print_dual(
        "Note: the GLOBAL integral above is expected to be close to 0 (no net charge leaves "
        "the simulation cell) -- the point of this plot is the nonzero LOCAL accumulation/"
        "depletion, not the global integral.", f_out)

    cube_path = None
    if args.cube:
        print_section("[2] CUBE FILE", f_out)
        geo_file, geometry = find_geometry_file(folders["combined"])
        if geometry is None:
            print_dual(color_text(
                f"[ERROR] --cube requires a geometry file -- neither a '*.XV' nor "
                f"'structure.fdf' was found in '{folders['combined']}'.", 'red'), f_out)
        elif not np.allclose(geometry.cell, lattice, atol=1e-2):
            print_dual(color_text(
                f"[ERROR] The cell in '{geo_file}' does not match the density grid's own "
                "cell -- refusing to write a cube file with mismatched geometry.", 'red'), f_out)
        else:
            cube_path = os.path.join(args.output_dir, f"{stem}.cube")
            combined_grid_obj.grid = diff_data
            combined_grid_obj.set_geometry(geometry)
            combined_grid_obj.write(cube_path)
            print_dual(color_text(f"[OK] Cube file saved to '{cube_path}' (open in VESTA/VMD/"
                                  "Avogadro for real 3D isosurfaces).", 'green'), f_out)

    print_section("[3] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[4] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Data           : {files_written['charge']['dat']}", f_out)
    if files_written['charge']['gplot']:
        print_dual(f"Gnuplot        : {files_written['charge']['gplot']}", f_out)
    if cube_path:
        print_dual(f"Cube file      : {cube_path}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

if __name__ == "__main__":
    main()
