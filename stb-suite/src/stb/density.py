#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import os
import sys
import argparse
from datetime import datetime
import numpy as np

from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
sisl = require_sisl()

# ANSI Colors for terminal
from stb.core import citations
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.rho_io import (
    read_total, try_read_net_spin, compute_colorbar_range, report_quantity,
)

REPORT_FILE = "stb_density_report.txt"
BIB_FILE = "references.bib"


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Converts a SIESTA charge-density grid (.RHO) to 2D maps, 3D point clouds, or planar-averaged 1D profiles.", 'bold')}
Plots the total charge density by default. If the .RHO file is
spin-polarized, the spin density (up-down) is now detected and reported
automatically alongside it -- --spin no longer needs to be remembered just
to notice a magnetic calculation; pass --spin to look at ONLY the spin
density (skip the charge section) instead. --rho2 subtracts a second .RHO
file's density from the first (Delta rho = rho1 - rho2), the standard way
to visualize charge transfer/bonding (e.g. adsorbate+substrate vs. isolated
fragments) -- both files must be on the same grid. --iso-min filters low
-|density| points out of a --3d export, since an unfiltered production-size
grid can produce a multi-GB file dominated by near-zero vacuum/interstitial
points. --profile writes the planar average along the chosen axis instead
of a map (useful for slabs/interfaces). --cube additionally writes a
standard Gaussian cube file (needs a <label>.XV or <label>.fdf for the
geometry) for real isosurface rendering in VESTA/VMD/Avogadro. --vmin/
--vmax fix the colorbar range manually for the primary quantity, e.g. to
compare several slices on the same scale.""",
        epilog="Example usage:\n"
               "  %(prog)s --label siesta\n"
               "  %(prog)s --label siesta --3d --iso-min 0.01\n"
               "  %(prog)s --label siesta --spin\n"
               "  %(prog)s --label adsorbed --rho2 isolated.RHO\n"
               "  %(prog)s --label siesta --profile --axis 2\n"
               "  %(prog)s --label siesta --cube --view\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", required=True, help="System Label (looks for .RHO file)")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write the data/gnuplot/cube/report files and "
                             "references.bib into (default: current directory). Created if "
                             "it doesn't exist.")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--3d", action="store_true", dest="full_3d", help="Export full 3D volume")
    mode_group.add_argument("--profile", action="store_true",
                        help="Write the planar average of the density along --axis as a 1D "
                             "profile (average taken over the other two axes) instead of a "
                             "2D/3D map -- useful for slabs, interfaces, and superlattices")
    parser.add_argument("-a", "--axis", type=int, default=2, choices=[0,1,2],
                        help="For a slice/profile: the axis normal to the cut plane (slice) or "
                             "the axis the profile varies along (profile) (0=X, 1=Y, 2=Z)")
    parser.add_argument("-p", "--pos", type=float, help="Position (Angstrom) of the 2D cut (ignored in --profile mode)")
    parser.add_argument("--spin", action="store_true",
                        help="Process ONLY the spin density (up-down), skipping the charge "
                             "density section entirely -- only available for spin-polarized "
                             ".RHO files. Without this flag, spin density is now detected and "
                             "reported automatically alongside the charge density whenever the "
                             ".RHO file has it, no flag needed.")
    parser.add_argument("--rho2", metavar="FILE", default=None,
                        help="A second .RHO file to subtract from the first (Delta rho = "
                             "rho1 - rho2), e.g. to visualize charge transfer upon adsorption "
                             "or defect formation. Must be on the same grid (same shape) as "
                             "<label>.RHO. If both files are spin-polarized, the spin "
                             "difference is computed too.")
    parser.add_argument("--iso-min", type=float, default=None, metavar="RHO",
                        help="For --3d: only export/preview points with |density| >= this "
                             "threshold (e/Ang^3) -- keeps the point cloud a manageable size "
                             "and visually meaningful instead of dominated by near-zero points")
    parser.add_argument("--cube", action="store_true",
                        help="Also write the primary quantity's full 3D grid as a Gaussian "
                             ".cube file (needs <label>.XV or <label>.fdf for the geometry) -- "
                             "for real 3D isosurface rendering in VESTA/VMD/Avogadro, which "
                             "handle volumetric data better than a point cloud. Applies to the "
                             "charge density, or the spin density if --spin was given -- not "
                             "both at once.")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Fix the colorbar/palette lower bound manually (e/Ang^3) for the "
                             "primary quantity, e.g. to compare several slices on the same "
                             "scale. Auto-detected spin density always uses its own "
                             "zero-symmetric range regardless of this.")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Fix the colorbar/palette upper bound manually (e/Ang^3), for the "
                             "primary quantity (see --vmin).")
    parser.add_argument("--contour", action="store_true",
                        help="Overlay contour lines on a 2D slice map/preview (ignored for --3d/--profile)")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a .gplot script next to each .dat file this run "
                             "generates. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of each quantity plotted "
                             "this run (slice heatmap / profile line / 3D scatter, matching the "
                             "mode). Off by default.")

    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-density {VERSION}")

    args = parser.parse_args()
    filename = f"{args.label}.RHO"
    stem = f"{args.label}_density"

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Density Plotter",
            "Exports RHO to 2D Maps, 3D Clouds, and matplotlib previews",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("DENSITY PLOTTER:", 'bold'))
    print("-" * 60)

    if not os.path.exists(filename):
        print(f"{COLORS['red']}[ERROR] File '{filename}' not found.{COLORS['reset']}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"[INFO] Reading: {color_text(filename, 'cyan')}")
    try:
        sile = sisl.get_sile(filename)
        if args.spin:
            grid_obj = try_read_net_spin(sile, filename, required=True)
            charge_grid_obj = None
        else:
            grid_obj = read_total(sile)
            charge_grid_obj = grid_obj
    except Exception as e:
        print(f"{COLORS['red']}[FATAL] {e}{COLORS['reset']}")
        sys.exit(1)

    lattice = grid_obj.lattice.cell
    origin = grid_obj.origin
    nx, ny, nz = grid_obj.grid.shape
    print(f"[INFO] Grid: {nx}x{ny}x{nz}")

    charge_data = charge_grid_obj.grid.copy() if charge_grid_obj is not None else None
    spin_data = None
    has_spin = False
    if not args.spin:
        print("[INFO] Checking for a spin component...")
        spin_grid_obj = try_read_net_spin(sile, filename)
        if spin_grid_obj is not None and spin_grid_obj.grid.shape == charge_data.shape:
            spin_data = spin_grid_obj.grid.copy()
            has_spin = True
            print(color_text(
                "[INFO] Spin-polarized .RHO detected -- the net spin (magnetization) density "
                "will be reported automatically alongside the total charge density.", 'cyan'))
        elif spin_grid_obj is not None:
            print(color_text(
                "[WARNING] A spin component exists but its shape doesn't match the charge "
                "grid -- ignoring it.", 'yellow'))
    else:
        spin_data = grid_obj.grid.copy()
        has_spin = True

    if args.rho2:
        if not os.path.exists(args.rho2):
            print(f"{COLORS['red']}[ERROR] File '{args.rho2}' not found.{COLORS['reset']}")
            sys.exit(1)
        print(f"[INFO] Reading second grid for difference: {color_text(args.rho2, 'cyan')}")
        try:
            sile2 = sisl.get_sile(args.rho2)
        except Exception as e:
            print(f"{COLORS['red']}[FATAL] {e}{COLORS['reset']}")
            sys.exit(1)

        if charge_data is not None:
            grid2 = read_total(sile2)
            if grid2.grid.shape != charge_data.shape:
                print(color_text(
                    f"[ERROR] Grid shape mismatch: {filename} is {charge_data.shape} but "
                    f"{args.rho2} is {grid2.grid.shape} -- both files must be on the same "
                    "grid to take a difference.", 'red'))
                sys.exit(1)
            charge_data = charge_data - grid2.grid

        if has_spin:
            spin2_obj = try_read_net_spin(sile2, args.rho2)
            if spin2_obj is not None and spin2_obj.grid.shape == spin_data.shape:
                spin_data = spin_data - spin2_obj.grid
            else:
                print(color_text(
                    f"[WARNING] --rho2: '{args.rho2}' has no matching spin component -- "
                    "the spin density section below is for the FIRST file alone, not a "
                    "difference.", 'yellow'))

    quantity_prefix = "Delta " if args.rho2 else ""

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

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-DENSITY REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}", f_out)
    print_dual(f"Grid file      : {filename}", f_out)
    print_dual(f"Mode           : {mode}" + (f" (axis {args.axis})" if mode != '3d' else ""), f_out)
    print_dual(f"Spin-polarized : {'yes' if has_spin else 'no'}"
               + (" (auto-detected)" if has_spin and not args.spin else ""), f_out)
    print_dual(f"--rho2         : {args.rho2 if args.rho2 else 'no'}", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    files_written = {}

    if charge_data is not None:
        label = f"{quantity_prefix}Charge Density"
        files_written["charge"] = report_quantity(
            "[1] CHARGE DENSITY", charge_data, lattice, origin, label, bool(args.rho2),
            stem, mode, args, args.output_dir, f_out)
    else:
        print_section("[1] CHARGE DENSITY", f_out)
        print_dual("Skipped (--spin: only the spin density is processed this run).", f_out)

    if has_spin and spin_data is not None:
        label = f"{quantity_prefix}Spin Density"
        spin_stem = stem if args.spin else f"{stem}_spin"
        files_written["spin"] = report_quantity(
            "[2] SPIN DENSITY", spin_data, lattice, origin, label, True,
            spin_stem, mode, args, args.output_dir, f_out, apply_manual_range=args.spin)

    cube_path = None
    if args.cube:
        print_section("[3] CUBE FILE", f_out)
        geometry, geo_file = read_sisl_geometry_xv_or_fdf(args.label)
        if geometry is None:
            print_dual(color_text(
                f"[ERROR] --cube requires a geometry file -- neither '{args.label}.XV' nor "
                f"'{args.label}.fdf' was found.", 'red'), f_out)
        elif not np.allclose(geometry.cell, lattice, atol=1e-2):
            print_dual(color_text(
                f"[ERROR] The cell in '{geo_file}' does not match the density grid's own "
                "cell -- refusing to write a cube file with mismatched geometry.", 'red'), f_out)
        else:
            cube_data = spin_data if args.spin else charge_data
            cube_path = os.path.join(args.output_dir, f"{stem}.cube")
            grid_obj.grid = cube_data
            grid_obj.set_geometry(geometry)
            grid_obj.write(cube_path)
            print_dual(color_text(f"[OK] Cube file saved to '{cube_path}' (open in VESTA/VMD/"
                                  "Avogadro for real 3D isosurfaces).", 'green'), f_out)

    print_section("[4] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    if "charge" in files_written:
        print_dual(f"Charge data    : {files_written['charge']['dat']}", f_out)
        if files_written['charge']['gplot']:
            print_dual(f"Charge gnuplot : {files_written['charge']['gplot']}", f_out)
    if "spin" in files_written:
        print_dual(f"Spin data      : {files_written['spin']['dat']}", f_out)
        if files_written['spin']['gplot']:
            print_dual(f"Spin gnuplot   : {files_written['spin']['gplot']}", f_out)
    if cube_path:
        print_dual(f"Cube file      : {cube_path}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

if __name__ == "__main__":
    main()
