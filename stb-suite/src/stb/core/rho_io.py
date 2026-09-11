"""Shared SIESTA .RHO reading/reporting helpers, extracted from density.py
once chargediff_analysis.py became a second consumer of the exact same
"read a total/spin density grid via sisl, then slice/3D/profile/gnuplot/
matplotlib it" need (same extract-on-second-use policy as the rest of
core/). density.py itself still owns the "diff against a second file"
detection (--rho2) and its own CLI/report structure -- only the pieces
generic over an arbitrary grid_data array moved here.
"""

import os
import sys

import numpy as np

from stb.core.cli import color_text, print_dual, print_section
from stb.core.grid_export import (
    NON_ORTHOGONAL_WARN_DEG, check_planar_orthogonality, write_gnuplot_script,
    integrated_charge, write_profile_data_file, write_data_file,
    plot_matplotlib_slice, plot_matplotlib_profile,
)


def read_total(sile):
    """Reads the TOTAL charge density via sisl's string-based `index='total'`
    convention (== 'up' + 'down' for a 2-component collinear .RHO, or just
    the only component for a non-spin-polarized one -- see try_read_net_spin's
    own docstring below for why this ISN'T the same as read_grid(index=0))."""
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
    density.py (and stb-bader's own read_spin_density) assumed. That earlier
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
    signed quantity (spin density, a difference density) gets a range
    symmetric around zero so the diverging palette's white midpoint
    actually lands on zero (autorange would center white at
    (data_min+data_max)/2 instead, only 0 for coincidentally symmetric
    data); a non-negative quantity is anchored at 0 for a consistent
    "white = no charge" baseline. Returns None for a degenerate range
    (e.g. an all-zero difference of identical files), leaving it to
    autorange instead of a meaningless [0:0].
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


def plot_matplotlib_3d(lattice, origin, grid_data, quantity_label, is_signed, cb_range, iso_min):
    """3D scatter preview of the point cloud -- same iso_min filtering as
    the exported point cloud (write_data_file), recomputed here rather than
    threaded through as a return value, matching stb-workfunction's own
    convention of the matplotlib preview independently recomputing its
    display arrays instead of sharing internals with the file writer."""
    import matplotlib.pyplot as plt

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
    if len(values) > 200_000:
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
    """Computes + writes everything for ONE quantity (charge or spin, or a
    difference of either) -- shared by every section that needs the
    slice/profile/3d dispatch logic, so it isn't duplicated per caller.
    Returns a dict of the files actually written, for the closing SUMMARY
    section.
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
