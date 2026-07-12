#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""STM image simulator (Tersoff-Hamann approximation): turns a SIESTA
energy-integrated local density of states grid (<label>.LDOS, written by
`%block LocalDensityOfStates Emin Emax` in the fdf -- the bias window is a
property of that SIESTA run, not something this tool controls after the
fact) into a simulated STM image, the same way a real tip does under the
Tersoff-Hamann approximation: tunneling current is proportional to the
LDOS at the tip position.

Two modes:
  --mode height  : the classic 2D LDOS map at a fixed height above the
                   surface (constant-height STM).
  --mode current : for each (x, y), finds the height where the LDOS first
                   drops to the --iso threshold coming in from vacuum --
                   the topography a real tip follows to hold tunneling
                   current constant (constant-current STM, the usual
                   textbook STM image).

Needs sisl (like stb-density) and a geometry file (<label>.XV/.fdf) to
locate the topmost surface atom -- heights are always reported relative to
it, not to the raw cell origin.

Known limitations:
 - Only images the surface exposed in the +axis direction beyond the
   topmost atom. A slab with two exposed faces (vacuum on both sides) only
   gets its "top" one imaged this way -- flip the structure first (e.g.
   via a mirrored geometry) to image the other face.
 - Assumes the chosen vacuum axis is Cartesian-aligned (true for every
   slab this suite's own tools produce, e.g. stb-slab); a sheared cell
   along that axis gets a printed warning but is not corrected for.
 - --iso is a relative LDOS threshold in the .LDOS file's own units
   (e/Bohr^3), not a calibrated absolute tunneling current -- Tersoff-Hamann
   is a proportionality, not an absolute current model.
"""

VERSION = "1.0.0"

import argparse
import os
import sys
import numpy as np

from stb.core.cli import color_text, show_intro
from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
from stb.core import kspace

VACUUM_GAP_ANG = 10.0  # same threshold used by stb-workfunction/stb-mlrelax


def _resolve_ldos_file(label):
    primary = f"{label}.LDOS"
    if os.path.isfile(primary):
        return primary
    fallback = f"{label}.STM.LDOS"
    if os.path.isfile(fallback):
        print(f"[INFO] No '{primary}' found; using '{fallback}' instead.")
        return fallback
    return None


def resolve_axis(geometry, requested_axis, grid_cell):
    """Returns (axis, vacuum_axes). Unlike stb-workfunction's analogous
    helper (which degrades gracefully to axis=2 since a work function can
    still be computed from an assumed axis), STM imaging without a genuine,
    confirmed vacuum axis and a matching geometry is meaningless -- both
    are hard requirements here, raising ValueError instead of guessing."""
    if geometry is None:
        raise ValueError(
            "No geometry file found -- stb-stm needs <label>.XV or <label>.fdf "
            "(or --geometry-file) to locate the topmost surface atom."
        )
    if not np.allclose(geometry.cell, grid_cell, atol=1e-2):
        raise ValueError(
            "The geometry file's cell does not match the .LDOS grid's own cell -- "
            "they likely don't belong to the same calculation."
        )
    vacuum_axes = [bool(v) for v in kspace.detect_vacuum_axes(geometry.fxyz, geometry.cell, VACUUM_GAP_ANG)]

    if requested_axis is not None:
        return requested_axis, vacuum_axes

    n_vacuum = sum(vacuum_axes)
    if n_vacuum != 1:
        raise ValueError(
            f"Expected exactly 1 vacuum-padded axis to auto-detect the surface normal, "
            f"found {n_vacuum} ({vacuum_axes}). STM imaging needs a slab (bulk/wire/"
            "molecule geometries don't have a well-defined exposed surface) -- pass "
            "--axis explicitly if this really is a slab."
        )
    return vacuum_axes.index(True), vacuum_axes


def check_axis_alignment(lattice, axis):
    """Warns (doesn't block) if the chosen axis's lattice vector isn't
    dominantly aligned with a single global Cartesian direction -- this
    tool treats grid index -> position along `axis` as a plain fractional
    scaling of that vector's own length, which only means "height above the
    surface" in the usual sense if the vector doesn't leak into the other
    two Cartesian directions."""
    vec = lattice[axis]
    norm = np.linalg.norm(vec)
    dominant = np.max(np.abs(vec))
    if dominant < norm * 0.999:
        print(color_text(
            f"[WARNING] Lattice vector along --axis {axis} is not aligned with a single "
            "Cartesian direction -- 'height above the surface' along this axis may not "
            "mean what you expect for a sheared cell.", 'yellow'))


def write_stm_data(output_file, u_vals, v_vals, data2d, value_label):
    """pm3d-blocked (X Y Value) data file, same blank-line-per-row layout
    as stb-density's slice mode. NaN entries (constant-current points that
    never reached --iso) are written as the literal text 'NaN', which
    gnuplot treats as missing data rather than a real value."""
    with open(output_file, 'w') as f:
        f.write(f"# Generated by STB. Format: {value_label}\n")
        for i, u in enumerate(u_vals):
            for j, v in enumerate(v_vals):
                val = data2d[i, j]
                val_str = "NaN" if np.isnan(val) else f"{val:.6e}"
                f.write(f"{u:.6f} {v:.6f} {val_str}\n")
            f.write("\n")


def write_stm_gplot(output_file, gplot_file, title, cb_label, u_label, v_label, is_signed=False):
    pdf_name = gplot_file.rsplit('.', 1)[0] + ".pdf"
    lines = []
    lines.append('set terminal pdfcairo enhanced color font "Arial,14" size 6,5\n')
    lines.append(f'set output "{pdf_name}"\n')
    lines.append(f'set title "{title}"\n')
    lines.append(f'set xlabel "{u_label}"\n')
    lines.append(f'set ylabel "{v_label}"\n')
    lines.append(f'set cblabel "{cb_label}"\n')
    lines.append('set size ratio -1\n')
    lines.append('set tics out nomirror\n')
    lines.append('set pm3d map\n')
    lines.append('set pm3d interpolate 4,4\n')
    if is_signed:
        lines.append('set palette defined (0 "#000099", 1 "#2255cc", 2 "#ffffff", 3 "#cc5522", 4 "#990000")\n')
    else:
        lines.append('set palette defined (0 "#000000", 1 "#4b0082", 2 "#b8860b", 3 "#ffd700", 4 "#ffffff")\n')
    lines.append(f'splot "{os.path.basename(output_file)}" using 1:2:3 with pm3d\n')
    with open(gplot_file, 'w') as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Simulate an STM image (Tersoff-Hamann approximation) from a SIESTA "
                    "energy-integrated LDOS grid (<label>.LDOS).",
        epilog="Example usage:\n"
               "  stb-stm --label siesta --mode current --iso 0.001\n"
               "  stb-stm --label siesta --mode height --z 3.0\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.LDOS (falling back to "
                             "<label>.STM.LDOS with a note) and <label>.XV/<label>.fdf for the "
                             "topmost-atom reference. Mutually exclusive with --file/--geometry-file.")
    parser.add_argument("--file", dest="ldos_file", type=str, default=None,
                        help="Explicit path to the .LDOS/.STM.LDOS file (alternative to --label).")
    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Explicit .XV/.fdf geometry path (alternative/override to --label's "
                             "auto-detect). Required if --file is used instead of --label.")

    parser.add_argument("--axis", type=int, default=None, choices=[0, 1, 2],
                        help="Surface-normal axis (0=X, 1=Y, 2=Z). Auto-detected from the "
                             "geometry's vacuum padding if omitted and exactly one vacuum axis "
                             "is found.")

    parser.add_argument("--mode", type=str, choices=["height", "current"], default="current",
                        help="'current' (default): classic constant-current STM image -- the "
                             "topography a tip follows to hold the LDOS at --iso. 'height': "
                             "constant-height LDOS map at --z above the surface.")

    parser.add_argument("--z", type=float, default=None,
                        help="Height (Ang) above the topmost surface atom for --mode height. "
                             "Required for --mode height.")
    parser.add_argument("--iso", type=float, default=None,
                        help="LDOS threshold (same units as the .LDOS file, e/Bohr^3) defining "
                             "the constant-current contour for --mode current. Required for "
                             "--mode current.")
    parser.add_argument("--z-max", type=float, default=None,
                        help="Upper bound (Ang above the topmost atom) of the search/plot "
                             "window along --axis. Default: distance to the cell boundary.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write stm_<mode>.dat and stm_<mode>.gplot into "
                             "(default: current directory). Created if it doesn't exist.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-stm {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.mode == "height" and args.z is None:
        parser.error("--z is required when --mode is 'height'.")
    if args.mode == "current" and args.iso is None:
        parser.error("--iso is required when --mode is 'current'.")

    if args.label:
        if args.ldos_file or args.geometry_file:
            parser.error("--label cannot be combined with --file/--geometry-file.")
        args.ldos_file = _resolve_ldos_file(args.label)
    elif not args.ldos_file:
        parser.error("one of --label or --file is required.")
    elif not args.geometry_file:
        parser.error("--geometry-file is required when using --file instead of --label.")

    if args.ldos_file is None or not os.path.isfile(args.ldos_file):
        tried = f"'{args.label}.LDOS' and '{args.label}.STM.LDOS'" if args.label else f"'{args.ldos_file}'"
        parser.error(f"No LDOS file found ({tried}). Generate one with SIESTA's "
                     "'%block LocalDensityOfStates Emin Emax' in the fdf.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STM SIMULATOR:", 'bold'))
    print("-" * 60)

    sisl = require_sisl()

    print(f"[INFO] Reading LDOS grid from '{args.ldos_file}' ...")
    grid_obj = sisl.get_sile(args.ldos_file).read_grid()
    grid_data = grid_obj.grid
    lattice = grid_obj.lattice.cell
    print(f"[INFO] Grid: {grid_data.shape[0]}x{grid_data.shape[1]}x{grid_data.shape[2]}, "
          f"LDOS range [{grid_data.min():.6e}, {grid_data.max():.6e}]")

    print("[INFO] Resolving geometry / surface-normal axis ...")
    if args.geometry_file:
        geometry = sisl.get_sile(args.geometry_file).read_geometry()
        geo_source = args.geometry_file
    else:
        geometry, geo_source = read_sisl_geometry_xv_or_fdf(args.label)

    try:
        axis, vacuum_axes = resolve_axis(geometry, args.axis, lattice)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)
    print(f"[INFO] Using axis {axis} (from '{geo_source}') as the surface normal.")
    check_axis_alignment(lattice, axis)

    z_top = float(geometry.xyz[:, axis].max())
    axis_len = float(np.linalg.norm(lattice[axis]))
    n_axis = grid_data.shape[axis]
    z_max_default = axis_len - z_top  # distance from the topmost atom to the far cell boundary
    z_max = args.z_max if args.z_max is not None else z_max_default
    print(f"[INFO] Topmost surface atom at {z_top:.3f} Ang along axis {axis}; "
          f"search window up to {z_max:.3f} Ang above it.")

    # Move the chosen axis to the end so the rest of the code is axis-agnostic.
    data = np.moveaxis(grid_data, axis, -1)
    other_axes = [i for i in range(3) if i != axis]
    nu, nv = data.shape[0], data.shape[1]
    u_len = float(np.linalg.norm(lattice[other_axes[0]]))
    v_len = float(np.linalg.norm(lattice[other_axes[1]]))
    u_vals = (np.arange(nu) / nu) * u_len
    v_vals = (np.arange(nv) / nv) * v_len

    z_positions = (np.arange(n_axis) / n_axis) * axis_len
    heights = z_positions - z_top  # height above the topmost atom, can be negative

    os.makedirs(args.output_dir, exist_ok=True)

    axis_names = ['X', 'Y', 'Z']
    u_label = f"{axis_names[other_axes[0]]} (Angstrom)"
    v_label = f"{axis_names[other_axes[1]]} (Angstrom)"

    if args.mode == "height":
        idx = int(np.argmin(np.abs(heights - args.z)))
        actual_height = float(heights[idx])
        print(f"[INFO] Mode: constant-height. Requested {args.z:.3f} Ang, using nearest grid "
              f"height {actual_height:.3f} Ang (index {idx}/{n_axis}).")
        image = data[:, :, idx]
        print(f"[INFO] LDOS at this height: min={image.min():.6e}, max={image.max():.6e}, "
              f"mean={image.mean():.6e}")

        out_data = os.path.join(args.output_dir, "stm_height.dat")
        out_gplot = os.path.join(args.output_dir, "stm_height.gplot")
        write_stm_data(out_data, u_vals, v_vals, image, "LDOS (e/Bohr^3)")
        write_stm_gplot(out_data, out_gplot,
                        f"STM constant-height image (z = {actual_height:.2f} Ang above surface)",
                        "LDOS (e/Bohr^3)", u_label, v_label, is_signed=False)

    else:
        window_mask = (heights >= 0) & (heights <= z_max)
        if not window_mask.any():
            print(color_text(
                "[ERROR] --z-max is smaller than the grid's own resolution along --axis -- "
                "no grid points fall inside the search window.", 'red'))
            sys.exit(1)
        window_idx = np.where(window_mask)[0]
        # Order from farthest-out (largest height) down to the surface, so the
        # first crossing found is the outermost point still at/above --iso --
        # exactly where a real tip retreating from the surface would stop.
        window_idx_desc = window_idx[np.argsort(-heights[window_idx])]
        window = data[:, :, window_idx_desc]  # (nu, nv, nwindow), ordered outside-in
        mask = window >= args.iso
        found = mask.any(axis=-1)
        first_idx = np.argmax(mask, axis=-1)  # first True along the outside-in order

        height_map = np.full((nu, nv), np.nan)
        flat_heights = heights[window_idx_desc]
        height_map[found] = flat_heights[first_idx[found]]

        n_missing = int((~found).sum())
        print(f"[INFO] Mode: constant-current (iso = {args.iso:.6e}).")
        if n_missing:
            pct = 100.0 * n_missing / found.size
            print(color_text(
                f"[WARNING] {n_missing}/{found.size} points ({pct:.1f}%) never reached --iso "
                "within the search window -- consider a larger --z-max or a lower --iso. "
                "These points are written as NaN.", 'yellow'))
        valid = height_map[found]
        if valid.size:
            print(f"[INFO] Height map: min={valid.min():.3f} Ang, max={valid.max():.3f} Ang, "
                  f"mean={valid.mean():.3f} Ang (corrugation: {valid.max() - valid.min():.3f} Ang)")

        out_data = os.path.join(args.output_dir, "stm_current.dat")
        out_gplot = os.path.join(args.output_dir, "stm_current.gplot")
        write_stm_data(out_data, u_vals, v_vals, height_map, "Height (Angstrom)")
        write_stm_gplot(out_data, out_gplot,
                        f"STM constant-current image (iso = {args.iso:.3e})",
                        "Height above surface (Angstrom)", u_label, v_label, is_signed=False)

    print(f"[INFO] Data written to '{color_text(out_data, 'green')}'")
    print(f"[INFO] Gnuplot script written to '{color_text(out_gplot, 'green')}'")
    print(f"[INFO] To plot: {color_text(f'gnuplot {out_gplot}', 'cyan')}")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
