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
import matplotlib.pyplot as plt

from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
sisl = require_sisl()

# ANSI Colors for terminal
from stb.core import citations
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.grid_export import (
    NON_ORTHOGONAL_WARN_DEG, check_planar_orthogonality, write_gnuplot_script,
    integrated_charge, write_profile_data_file, write_data_file,
)

REPORT_FILE = "stb_density_report.txt"
BIB_FILE = "references.bib"

# Above this many points, an unfiltered --3d matplotlib scatter tends to be slow to render
# and to produce an overplotted, visually meaningless cloud -- a hint to use --iso-min,
# not a hard limit (the export/gnuplot path is unaffected, only the --view preview).
MPL_3D_POINT_WARN = 200_000


def read_total(sile):
    """Reads the TOTAL charge density via sisl's string-based `index='total'`
    convention (== 'up' + 'down' for a 2-component collinear .RHO, or just
    the only component for a non-spin-polarized one -- see the module note
    below for why this ISN'T the same as read_grid(index=0))."""
    return sile.read_grid(index='total')


def try_read_net_spin(sile, source_desc, required=False):
    """Reads the NET SPIN (magnetization) density via sisl's `index='z'`
    convention (== 'up' - 'down'), or returns None if this .RHO isn't
    spin-polarized (only 1 component -- sisl raises ValueError requesting 2
    factors ([1,-1]) from a single-component grid). If `required` (--spin
    was explicitly given), a missing spin component is a fatal, clearly
    -worded error instead of sisl's own bare exception; otherwise (auto
    -detection) it's silently treated as "not spin-polarized", the same
    best-effort convention as stb-bader's own read_spin_density.

    ** Why `index='z'`/`index='total'`, not the raw integers 0/1 **: for a
    genuinely spin-polarized (nspin=2) SIESTA .RHO, sisl's read_grid(index=N)
    with a plain integer returns the RAW stored component directly --
    component 0 is the up-spin channel, component 1 is the down-spin
    channel, NOT "total charge" and "spin density" as an earlier version of
    this tool (and stb-bader's own read_spin_density) assumed. That earlier
    assumption was a real, verified bug: on a real spin-polarized O2
    calculation (test/6-utils/3-cube/o2.RHO, textbook triplet ground state,
    SIESTA's own log reports |S| = 2.0), naively reading index=0 as "charge"
    integrates to 7.0 e (the up channel alone) instead of the correct 12.0 e
    (2 O atoms x 6 valence electrons), and index=1 as "spin" integrates to
    5.0 e (the down channel alone) instead of the correct 2.0 e net moment.
    `index='total'` (sums components with weight [1, 1]) and `index='z'`
    (weight [1, -1]) are sisl's own combination for exactly this -- the
    same convention already used, correctly, by stb-cube's own SPIN_INDEX
    (`cube.py`: `{'total': 'total', 'up': 0, 'down': 1, 'diff': 'z'}`).
    """
    try:
        return sile.read_grid(index='z')
    except (IndexError, ValueError):
        if required:
            print(color_text(
                f"[ERROR] --spin requested but {source_desc} has no spin component -- "
                "this doesn't look like a spin-polarized calculation.", 'red'))
            sys.exit(1)
        return None


def compute_colorbar_range(data_vmin, data_vmax, is_signed, vmin_arg=None, vmax_arg=None):
    """Colorbar range: an explicit --vmin/--vmax always wins; otherwise a
    signed quantity (spin density, --rho2 difference) gets a range symmetric
    around zero so the diverging palette's white midpoint actually lands on
    zero (autorange would center white at (data_min+data_max)/2 instead,
    only 0 for coincidentally symmetric data); a non-negative quantity is
    anchored at 0 for a consistent "white = no charge" baseline. Returns
    None for a degenerate range (e.g. an all-zero --rho2 diff of identical
    files), leaving it to autorange instead of a meaningless [0:0].
    """
    if vmin_arg is not None or vmax_arg is not None:
        cb_min = vmin_arg if vmin_arg is not None else data_vmin
        cb_max = vmax_arg if vmax_arg is not None else data_vmax
    elif is_signed:
        m = max(abs(data_vmin), abs(data_vmax))
        cb_min, cb_max = -m, m
    else:
        cb_min, cb_max = 0.0, data_vmax
    if cb_max > cb_min:
        return (cb_min, cb_max)
    return None


def plot_matplotlib_slice(lattice, origin, grid_data, axis, slice_idx, quantity_label,
                          is_signed, cb_range, pos_val, contour):
    """2D heatmap preview of one slice, in real-space coordinates -- same
    blue-white-red (signed) / white-yellow-red (non-negative) convention
    as the gnuplot palette in write_gnuplot_script, so the two views read
    as the same plot."""
    nx, ny, nz = grid_data.shape
    if axis == 2:
        arr, u_n, v_n, ulab, vlab = grid_data[:, :, slice_idx], nx, ny, "X", "Y"
    elif axis == 1:
        arr, u_n, v_n, ulab, vlab = grid_data[:, slice_idx, :], nx, nz, "X", "Z"
    else:
        arr, u_n, v_n, ulab, vlab = grid_data[slice_idx, :, :], ny, nz, "Y", "Z"

    cmap = "RdBu_r" if is_signed else "YlOrRd"
    vmin, vmax = cb_range if cb_range else (None, None)

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(arr.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    if contour:
        ax.contour(arr.T, colors="k", linewidths=0.5, levels=10)
    fig.colorbar(im, ax=ax, label=f"{quantity_label} (e/Ang^3)")
    ax.set_xlabel(f"{ulab} (grid index, {u_n} pts)")
    ax.set_ylabel(f"{vlab} (grid index, {v_n} pts)")
    ax.set_title(f"{quantity_label} slice ({'XYZ'[axis]}={pos_val:.2f} Ang)")
    fig.tight_layout()
    plt.show()


def plot_matplotlib_profile(lattice, axis, grid_data, quantity_label, is_signed):
    """Line plot of the planar-averaged profile along `axis` -- same
    convention as stb-workfunction's own plot_matplotlib."""
    avg_axes = tuple(a for a in (0, 1, 2) if a != axis)
    profile = grid_data.mean(axis=avg_axes)
    n = profile.shape[0]
    axis_len = np.linalg.norm(lattice[axis])
    positions = (np.arange(n) / n) * axis_len

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(positions, profile, color="#cc5522", linewidth=2)
    if is_signed:
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel(f"{'XYZ'[axis]} (Ang)")
    ax.set_ylabel(f"{quantity_label} (e/Ang^3)")
    ax.set_title(f"Planar-averaged {quantity_label} along {'XYZ'[axis]}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()


def plot_matplotlib_3d(lattice, origin, grid_data, quantity_label, is_signed, cb_range, iso_min):
    """3D scatter preview of the point cloud -- same iso_min filtering as
    the exported point cloud (write_data_file), recomputed here rather than
    threaded through as a return value, matching stb-workfunction's own
    convention of the matplotlib preview independently recomputing its
    display arrays instead of sharing internals with the file writer."""
    nx, ny, nz = grid_data.shape
    ix, iy, iz = np.meshgrid(range(nx), range(ny), range(nz), indexing='ij')
    frac_coords = np.vstack([ix.flatten() / nx, iy.flatten() / ny, iz.flatten() / nz]).T
    values = grid_data.flatten()

    if iso_min is not None:
        mask = np.abs(values) >= iso_min
        frac_coords = frac_coords[mask]
        values = values[mask]

    if len(values) == 0:
        print(color_text("[WARNING] --view: no points survived --iso-min -- skipping the 3D preview.", 'yellow'))
        return
    if len(values) > MPL_3D_POINT_WARN:
        print(color_text(
            f"[NOTE] --view: {len(values)} points to render in 3D -- this may be slow/cluttered. "
            "Consider a tighter --iso-min for a clearer preview.", 'yellow'))

    real_coords = origin + np.dot(frac_coords, lattice)
    cmap = "RdBu_r" if is_signed else "YlOrRd"
    vmin, vmax = cb_range if cb_range else (None, None)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(real_coords[:, 0], real_coords[:, 1], real_coords[:, 2],
                    c=values, cmap=cmap, vmin=vmin, vmax=vmax, s=4, depthshade=True)
    fig.colorbar(sc, ax=ax, label=f"{quantity_label} (e/Ang^3)", shrink=0.7)
    ax.set_xlabel("X (Ang)")
    ax.set_ylabel("Y (Ang)")
    ax.set_zlabel("Z (Ang)")
    ax.set_title(f"3D {quantity_label} (point cloud)")
    fig.tight_layout()
    plt.show()


def report_quantity(section_label, grid_data, lattice, origin, quantity_label, is_signed,
                    stem, mode, args, output_dir, f_out, apply_manual_range=True):
    """Computes + writes everything for ONE quantity (charge or spin, or
    their --rho2 difference) -- shared by both the always-present charge
    section and the conditional (auto-detected) spin section below, so the
    slice/profile/3d dispatch logic isn't duplicated for each. Returns a
    dict of the files actually written, for the closing SUMMARY section.
    """
    print_section(section_label, f_out)

    total_value = integrated_charge(grid_data, lattice)
    print_dual(f"Integrated {quantity_label}: {total_value:.4f} e "
               "(sanity check against the expected valence electron count/magnetic moment)",
               f_out)

    nx, ny, nz = grid_data.shape
    slice_idx = 0
    pos_val = 0.0
    plane_angle_deg = 90.0

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
            print_dual(f"Using center of axis {args.axis}: {pos_val:.2f} Ang (no --pos given)", f_out)

        if slice_idx < 0 or slice_idx >= dim_size:
            print_dual(color_text("[ERROR] Position out of bounds.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)

        plane_angle_deg = check_planar_orthogonality(lattice, args.axis)
        if abs(plane_angle_deg - 90.0) > NON_ORTHOGONAL_WARN_DEG:
            print_dual(color_text(
                f"[WARNING] The cut plane is skewed ({plane_angle_deg:.1f} deg between its two "
                "in-plane lattice vectors, not 90) -- gnuplot's pm3d map / matplotlib's imshow "
                "both assume a roughly rectangular grid and may render this visibly distorted "
                "even though the data file's coordinates are correct (common for hexagonal "
                "cells).", 'yellow'), f_out)

    dat_path = os.path.join(output_dir, f"{stem}.dat")
    if mode == 'profile':
        data_vmin, data_vmax = write_profile_data_file(grid_data, lattice, args.axis, dat_path)
    else:
        data_vmin, data_vmax = write_data_file(
            grid_data, lattice, origin, dat_path, mode, slice_idx, args.axis,
            iso_min=args.iso_min if mode == '3d' else None)

    cb_range = None
    if mode in ('slice', '3d'):
        vmin_arg = args.vmin if apply_manual_range else None
        vmax_arg = args.vmax if apply_manual_range else None
        cb_range = compute_colorbar_range(data_vmin, data_vmax, is_signed, vmin_arg, vmax_arg)

    gplot_path = None
    if args.save_gnuplot:
        contour = args.contour and mode == 'slice'
        write_gnuplot_script(dat_path, dat_path, mode, quantity_label, is_signed, args.axis,
                             pos_val, plane_angle_deg, cb_range=cb_range, contour=contour)
        gplot_path = dat_path.rsplit('.', 1)[0] + ".gplot"

    print_dual(color_text(f"[OK] Data written to '{dat_path}'.", 'green'), f_out)
    if gplot_path:
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)

    if args.view:
        if mode == 'profile':
            plot_matplotlib_profile(lattice, args.axis, grid_data, quantity_label, is_signed)
        elif mode == 'slice':
            plot_matplotlib_slice(lattice, origin, grid_data, args.axis, slice_idx, quantity_label,
                                  is_signed, cb_range, pos_val, args.contour)
        else:
            plot_matplotlib_3d(lattice, origin, grid_data, quantity_label, is_signed, cb_range,
                               args.iso_min)

    return {"total": total_value, "dat": dat_path, "gplot": gplot_path}


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
