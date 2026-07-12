#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.22.0"

import os
import sys
import argparse
import numpy as np

from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
sisl = require_sisl()

# ANSI Colors for terminal
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.grid_export import (
    NON_ORTHOGONAL_WARN_DEG, check_planar_orthogonality, write_gnuplot_script,
    integrated_charge, write_profile_data_file, write_data_file,
)


def read_component(sile, component_index, source_desc):
    """Reads one grid component (0=charge/total, 1=spin) via sisl, giving a
    clean error instead of sisl's own bare IndexError when --spin is
    requested on a non-spin-polarized .RHO (which only has component 0).
    """
    try:
        return sile.read_grid(index=component_index)
    except IndexError:
        if component_index == 0:
            raise
        print(color_text(
            f"[ERROR] --spin requested but {source_desc} has no spin component -- "
            "this doesn't look like a spin-polarized calculation.", 'red'))
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Converts a SIESTA charge-density grid (.RHO) to Gnuplot 2D maps, 3D point clouds, or planar-averaged 1D profiles.", 'bold')}
By default plots the total charge density. --spin plots the spin density
(up-down) instead, if the .RHO file is spin-polarized. --rho2 subtracts a
second .RHO file's density from the first (Delta rho = rho1 - rho2), the
standard way to visualize charge transfer/bonding (e.g. adsorbate+substrate
vs. isolated fragments) -- both files must be on the same grid. --iso-min
filters low-|density| points out of a --3d export, since an unfiltered
production-size grid can produce a multi-GB file dominated by
near-zero vacuum/interstitial points. --profile writes the planar average
along the chosen axis instead of a map (useful for slabs/interfaces). --cube
additionally writes a standard Gaussian cube file (needs a <label>.XV or
<label>.fdf for the geometry) for real isosurface rendering in VESTA/VMD/
Avogadro. --vmin/--vmax fix the colorbar range manually, e.g. to compare
several slices on the same scale.""",
        epilog="Example usage:\n"
               "  %(prog)s --label siesta\n"
               "  %(prog)s --label siesta --3d --iso-min 0.01\n"
               "  %(prog)s --label siesta --spin\n"
               "  %(prog)s --label adsorbed --rho2 isolated.RHO\n"
               "  %(prog)s --label siesta --profile --axis 2\n"
               "  %(prog)s --label siesta --cube\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", required=True, help="System Label (looks for .RHO file)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output data filename (default: <label>_density.dat)")
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
                        help="Plot the spin density (up-down) instead of the total charge "
                             "density -- only available for spin-polarized .RHO files")
    parser.add_argument("--rho2", metavar="FILE", default=None,
                        help="A second .RHO file to subtract from the first (Delta rho = "
                             "rho1 - rho2), e.g. to visualize charge transfer upon adsorption "
                             "or defect formation. Must be on the same grid (same shape) as "
                             "<label>.RHO.")
    parser.add_argument("--iso-min", type=float, default=None, metavar="RHO",
                        help="For --3d: only export points with |density| >= this threshold "
                             "(e/Ang^3) -- keeps the point cloud a manageable size and "
                             "visually meaningful instead of dominated by near-zero points")
    parser.add_argument("--cube", action="store_true",
                        help="Also write the full 3D grid as a Gaussian .cube file (needs "
                             "<label>.XV or <label>.fdf for the geometry) -- for real 3D "
                             "isosurface rendering in VESTA/VMD/Avogadro, which handle "
                             "volumetric data better than a gnuplot point cloud")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Fix the colorbar/palette lower bound manually (e/Ang^3), e.g. "
                             "to compare several slices on the same scale")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Fix the colorbar/palette upper bound manually (e/Ang^3)")
    parser.add_argument("--contour", action="store_true",
                        help="Overlay contour lines on a 2D slice map (ignored for --3d/--profile)")

    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-density {VERSION}")

    args = parser.parse_args()
    filename = f"{args.label}.RHO"
    output = args.output if args.output else f"{args.label}_density.dat"

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Density Plotter",
            "Exports RHO to Gnuplot (2D Maps & 3D Clouds)",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("DENSITY PLOTTER:", 'bold'))
    print("-" * 60)

    if not os.path.exists(filename):
        print(f"{COLORS['red']}[ERROR] File '{filename}' not found.{COLORS['reset']}")
        sys.exit(1)

    component_index = 1 if args.spin else 0

    print(f"[INFO] Reading: {color_text(filename, 'cyan')}")
    try:
        sile = sisl.get_sile(filename)
        grid_obj = read_component(sile, component_index, filename)
    except Exception as e:
        print(f"{COLORS['red']}[FATAL] {e}{COLORS['reset']}")
        sys.exit(1)

    grid_data = grid_obj.grid
    lattice = grid_obj.lattice.cell
    origin = grid_obj.origin
    nx, ny, nz = grid_data.shape

    print(f"[INFO] Grid: {nx}x{ny}x{nz}")

    is_signed = args.spin
    quantity_label = "Spin Density" if args.spin else "Charge Density"

    if args.rho2:
        if not os.path.exists(args.rho2):
            print(f"{COLORS['red']}[ERROR] File '{args.rho2}' not found.{COLORS['reset']}")
            sys.exit(1)
        print(f"[INFO] Reading second grid for difference: {color_text(args.rho2, 'cyan')}")
        try:
            sile2 = sisl.get_sile(args.rho2)
            grid_obj2 = read_component(sile2, component_index, args.rho2)
        except Exception as e:
            print(f"{COLORS['red']}[FATAL] {e}{COLORS['reset']}")
            sys.exit(1)

        if grid_obj2.grid.shape != grid_data.shape:
            print(color_text(
                f"[ERROR] Grid shape mismatch: {filename} is {grid_data.shape} but "
                f"{args.rho2} is {grid_obj2.grid.shape} -- both files must be on the same "
                "grid to take a difference.", 'red'))
            sys.exit(1)

        grid_data = grid_data - grid_obj2.grid
        quantity_label = f"Delta {quantity_label}"
        is_signed = True

    total_charge = integrated_charge(grid_data, lattice)
    print(f"[INFO] Integrated {quantity_label}: {color_text(f'{total_charge:.4f} e', 'cyan')} "
          "over the cell (sanity check against the expected valence electron count)")

    if args.profile:
        mode = 'profile'
    elif args.full_3d:
        mode = '3d'
    else:
        mode = 'slice'

    if args.contour and mode != 'slice':
        print(color_text(
            "[WARNING] --contour only applies to the default 2D slice mode -- ignoring.",
            'yellow'))

    if args.pos is not None and mode == 'profile':
        print(color_text(
            "[WARNING] --pos has no effect in --profile mode -- the profile covers the "
            "whole axis.", 'yellow'))

    slice_idx = 0
    pos_val = 0.0
    plane_angle_deg = 90.0
    data_vmin = data_vmax = None

    if mode == 'slice':
        dim_size = [nx, ny, nz][args.axis]
        if args.pos is not None:
            axis_vec = lattice[args.axis]
            axis_len = np.linalg.norm(axis_vec)
            slice_idx = int(round((args.pos / axis_len) * dim_size))
            pos_val = args.pos
        else:
            slice_idx = dim_size // 2
            pos_val = (slice_idx / dim_size) * np.linalg.norm(lattice[args.axis])
            print(f"{COLORS['yellow']}[NOTICE] Using center of axis {args.axis}: {pos_val:.2f} A{COLORS['reset']}")

        if slice_idx < 0 or slice_idx >= dim_size:
            print(f"{COLORS['red']}[ERROR] Position out of bounds.{COLORS['reset']}")
            sys.exit(1)

        plane_angle_deg = check_planar_orthogonality(lattice, args.axis)
        if abs(plane_angle_deg - 90.0) > NON_ORTHOGONAL_WARN_DEG:
            print(color_text(
                f"[WARNING] The cut plane is skewed ({plane_angle_deg:.1f} deg between its "
                "two in-plane lattice vectors, not 90) -- gnuplot's pm3d map assumes a "
                "roughly rectangular grid and may render this visibly distorted even though "
                "the data file's coordinates are correct (common for hexagonal cells).",
                'yellow'))

    if mode == 'profile':
        data_vmin, data_vmax = write_profile_data_file(grid_data, lattice, args.axis, output)
    else:
        data_vmin, data_vmax = write_data_file(grid_data, lattice, origin, output, mode, slice_idx,
                        args.axis, iso_min=args.iso_min if mode == '3d' else None)

    # Colorbar range: an explicit --vmin/--vmax always wins; otherwise a signed quantity
    # (spin density, --rho2 difference) gets a range symmetric around zero so the diverging
    # palette's white midpoint actually lands on zero (gnuplot's autorange would center white
    # at (data_min+data_max)/2 instead, which is only 0 for coincidentally symmetric data); a
    # non-negative quantity is anchored at 0 for a consistent "white = no charge" baseline.
    cb_range = None
    if mode in ('slice', '3d'):
        if args.vmin is not None or args.vmax is not None:
            cb_min = args.vmin if args.vmin is not None else data_vmin
            cb_max = args.vmax if args.vmax is not None else data_vmax
        elif is_signed:
            m = max(abs(data_vmin), abs(data_vmax))
            cb_min, cb_max = -m, m
        else:
            cb_min, cb_max = 0.0, data_vmax
        # A degenerate range (e.g. an all-zero --rho2 diff of identical files) would write
        # a meaningless "set cbrange [0:0]" -- leave it to gnuplot's own autorange instead.
        if cb_max > cb_min:
            cb_range = (cb_min, cb_max)

    contour = args.contour and mode == 'slice'
    write_gnuplot_script(output, output, mode, quantity_label, is_signed, args.axis, pos_val,
                        plane_angle_deg, cb_range=cb_range, contour=contour)

    if args.cube:
        geometry, geo_file = read_sisl_geometry_xv_or_fdf(args.label)
        if geometry is None:
            print(color_text(
                f"[ERROR] --cube requires a geometry file -- neither '{args.label}.XV' nor "
                f"'{args.label}.fdf' was found.", 'red'))
            sys.exit(1)
        if not np.allclose(geometry.cell, lattice, atol=1e-2):
            print(color_text(
                f"[ERROR] The cell in '{geo_file}' does not match the density grid's own "
                "cell -- refusing to write a cube file with mismatched geometry.", 'red'))
            sys.exit(1)
        cube_filename = output.rsplit('.', 1)[0] + ".cube"
        grid_obj.grid = grid_data
        grid_obj.set_geometry(geometry)
        grid_obj.write(cube_filename)
        print(f"[INFO] Cube file saved to '{color_text(cube_filename, 'green')}' "
              "(open in VESTA/VMD/Avogadro for real 3D isosurfaces)")

    print("\n[INFO] Complete job!")
    print("\n"+"-"*60)
    print(color_text("Density mapped. I am your density (destiny).\n\n", 'bold'))

if __name__ == "__main__":
    main()
