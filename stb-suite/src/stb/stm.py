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
 - The .LDOS file itself carries no record of the Emin/Emax bias window
   used to build it (that's a property of the SIESTA run, not the grid
   file) -- know your own %block LocalDensityOfStates settings if you need
   to report/reproduce the bias window.

Output/report style rewritten (v1.0.0 -> v2.0.0) to match the rest of the
Analysis category (stb-workfunction/stb-density): a numbered [0]...[6]
report, --save-report, --save-gnuplot (the .dat/.gplot files used to be
written unconditionally on every run), and --view (a matplotlib preview --
this tool previously had none at all). --iso/--z now default to 0.001
e/Bohr^3 / 3.0 Ang (a representative, literature-typical choice for each --
see the module docstring above on why --iso is relative, not absolute)
instead of being required with no default. The gnuplot script's own
`set output`/`splot` filenames are now always bare basenames, never
prefixed with --output-dir -- the user is expected to `cd` into that
directory before running gnuplot there directly, so an embedded directory
prefix pointed at the wrong place (a real, verified bug: with -o, the old
`set output "<output-dir>/stm_current.pdf"` line was wrong relative to
that convention).

A second real bug (v2.0.0 -> v2.1.0), found while trying this tool on a
real, externally-fetched CrS monolayer structure: the "topmost atom" used
to be a naive `xyz[:, axis].max()`, silently picking the WRONG bounding
atom whenever a structure's atoms straddle the periodic cell boundary with
the real vacuum gap in the *middle* of the cell instead of padded after
the atoms (this CrS cell's atoms sit at fractional z = 0, 0, 0.066, 0.934
-- the real ~87%-of-cell vacuum gap is between the 0.066/0.934 atoms, but
`.max()` picked the 0.934 one and searched into the tiny ~7% wraparound
sliver beyond it instead, collapsing the whole search window to ~1.5 Ang
instead of the genuine ~20 Ang vacuum region). Fixed with
`core/kspace.py::find_surface_reference`, which locates the atom
immediately below the LARGEST circular (wrap-aware) gap, the same gap
-finding logic `detect_vacuum_axes` already uses to decide THAT an axis is
vacuum-padded, now also used to find WHERE. Verified live: identical
numbers as before on every existing (conventional) fixture, and a real,
usable ~20 Ang search window (instead of a nonsensical, saturated ~1.5 Ang
one) on the CrS structure.
"""

VERSION = "2.1.0"

import argparse
import os
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
from stb.core import kspace

REPORT_FILE = "stb_stm_report.txt"
BIB_FILE = "references.bib"

VACUUM_GAP_ANG = 10.0  # same threshold used by stb-workfunction/stb-mlrelax

# Defaults chosen to just work out of the box on a typical slab, not as a
# universally "correct" value for every material/bias window:
#  - DEFAULT_ISO (e/Bohr^3): --iso is a RELATIVE threshold in the .LDOS
#    file's own units, not a calibrated absolute tunneling current (see the
#    module docstring). 1e-3 sits comfortably inside the vacuum-decay tail
#    for a typical slab LDOS grid (verified against this tool's own
#    graphene test fixture: max LDOS ~0.65, with only ~0.1% of grid points
#    above 0.44 -- 1e-3 is well past the near-atom, near-saturated region
#    without being so small that noise/grid resolution dominates).
#  - DEFAULT_Z (Ang): 3.0 Ang above the topmost atom is a commonly cited
#    representative constant-height STM tip-sample separation in the
#    Tersoff-Hamann literature -- close enough to the surface for real
#    corrugation contrast, far enough to sit in the decaying vacuum tail
#    rather than inside the atomic core region.
DEFAULT_ISO = 0.001
DEFAULT_Z = 3.0


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
        return (f"Lattice vector along --axis {axis} is not aligned with a single "
                "Cartesian direction -- 'height above the surface' along this axis may not "
                "mean what you expect for a sheared cell.")
    return None


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
    """Both filenames referenced INSIDE the script (the data file in
    `splot`, and the PDF in `set output`) are stripped to bare basenames --
    the user is expected to `cd` into --output-dir and run `gnuplot
    <name>.gplot` directly from there (same convention every other tool's
    .gplot script in this suite follows), so an embedded directory prefix
    would point at the wrong place relative to that cwd. Verified live:
    with -o some_dir, the old code left the directory prefix in the
    `set output` line (but not in `splot`, an existing inconsistency) --
    e.g. `set output "some_dir/stm_current.pdf"` while sitting inside
    some_dir already, which gnuplot would try (and typically fail) to
    write to some_dir/some_dir/stm_current.pdf."""
    gplot_basename = os.path.basename(gplot_file)
    pdf_name = gplot_basename.rsplit('.', 1)[0] + ".pdf"
    data_basename = os.path.basename(output_file)
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
    lines.append(f'splot "{data_basename}" using 1:2:3 with pm3d\n')
    with open(gplot_file, 'w') as f:
        f.writelines(lines)


def plot_matplotlib(u_vals, v_vals, image, u_label, v_label, cb_label, title):
    """Interactive preview mirroring the .gplot's own pm3d map: a filled
    2D color map, u/v axes in Angstrom, shared 'afmhot'-like sequential
    palette convention (this tool has no signed quantity -- both LDOS and
    height-above-surface are non-negative by construction)."""
    fig, ax = plt.subplots(figsize=(7, 6))
    mesh = ax.pcolormesh(u_vals, v_vals, image.T, shading='auto', cmap='afmhot')
    fig.colorbar(mesh, ax=ax, label=cb_label)
    ax.set_xlabel(u_label)
    ax.set_ylabel(v_label)
    ax.set_title(title)
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Simulate an STM image (Tersoff-Hamann approximation) from a SIESTA "
                    "energy-integrated LDOS grid (<label>.LDOS).",
        epilog="Example usage:\n"
               "  stb-stm --label siesta --mode current --iso 0.001\n"
               "  stb-stm --label siesta --mode height --z 3.0\n"
               "  stb-stm --label siesta --save-report --save-gnuplot --view\n",
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

    parser.add_argument("--z", type=float, default=DEFAULT_Z,
                        help=f"Height (Ang) above the topmost surface atom for --mode height "
                             f"(default: {DEFAULT_Z}, a commonly used representative constant-"
                             "height tip-sample separation).")
    parser.add_argument("--iso", type=float, default=DEFAULT_ISO,
                        help=f"LDOS threshold (same units as the .LDOS file, e/Bohr^3) defining "
                             f"the constant-current contour for --mode current (default: "
                             f"{DEFAULT_ISO} -- a relative threshold, not a calibrated absolute "
                             "current; tune it for your own system's LDOS magnitude if the "
                             "default leaves too many/few points reaching it, see [3] in the "
                             "report).")
    parser.add_argument("--z-max", type=float, default=None,
                        help="Upper bound (Ang above the topmost atom) of the search/plot "
                             "window along --axis. Default: distance to the cell boundary.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write stm_<mode>.dat/.gplot into (with --save-gnuplot) "
                             "and stb_stm_report.txt/references.bib (default: current directory). "
                             "Created if it doesn't exist.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write stm_<mode>.dat + stm_<mode>.gplot together. Off by "
                             "default -- this tool used to write both unconditionally on every run.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of the STM image before "
                             "finishing. Off by default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-stm {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

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

    print("\n" + color_text("STM SIMULATOR: reading input data", 'bold'))
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
    axis_warning = check_axis_alignment(lattice, axis)
    if axis_warning:
        print(color_text(f"[WARNING] {axis_warning}", 'yellow'))

    axis_len = float(np.linalg.norm(lattice[axis]))
    n_axis = grid_data.shape[axis]
    # frac_start is the atom immediately below the largest real vacuum gap --
    # NOT necessarily the same as a naive xyz[:, axis].max() (see
    # find_surface_reference's own docstring: a structure whose atoms
    # straddle the periodic cell boundary, with the real vacuum in the
    # middle of the cell instead of padded after the atoms, needs the gap
    # -aware reference or the search window silently collapses to a tiny,
    # wrong-direction sliver).
    frac_start, gap_size_frac = kspace.find_surface_reference(geometry.fxyz[:, axis])
    z_top = frac_start * axis_len
    z_max_default = gap_size_frac * axis_len
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

    axis_names = ['X', 'Y', 'Z']
    u_label = f"{axis_names[other_axes[0]]} (Angstrom)"
    v_label = f"{axis_names[other_axes[1]]} (Angstrom)"

    print("[INFO] Computing the STM image ...")
    image_stats = {}
    if args.mode == "height":
        idx = int(np.argmin(np.abs(heights - args.z)))
        actual_height = float(heights[idx])
        image = data[:, :, idx]
        cb_label = "LDOS (e/Bohr^3)"
        title = f"STM constant-height image (z = {actual_height:.2f} Ang above surface)"
        image_stats = {
            "requested_z": args.z, "actual_height": actual_height, "idx": idx,
            "min": float(image.min()), "max": float(image.max()), "mean": float(image.mean()),
        }
        n_missing = 0
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

        image = np.full((nu, nv), np.nan)
        flat_heights = heights[window_idx_desc]
        image[found] = flat_heights[first_idx[found]]
        cb_label = "Height above surface (Angstrom)"
        title = f"STM constant-current image (iso = {args.iso:.3e})"

        n_missing = int((~found).sum())
        valid = image[found]
        image_stats = {
            "iso": args.iso, "n_missing": n_missing, "n_total": int(found.size),
        }
        if valid.size:
            image_stats.update({
                "min": float(valid.min()), "max": float(valid.max()), "mean": float(valid.mean()),
                "corrugation": float(valid.max() - valid.min()),
            })

    os.makedirs(args.output_dir, exist_ok=True)
    out_data = os.path.join(args.output_dir, f"stm_{args.mode}.dat")
    out_gplot = os.path.join(args.output_dir, f"stm_{args.mode}.gplot")

    # --- From here on: the numbered, save-able report --------------------
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("\n===== STB-STM REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"LDOS file      : {args.ldos_file}", f_out)
    print_dual(f"Geometry source: {geo_source}", f_out)
    print_dual(f"Mode           : {args.mode}", f_out)
    if args.mode == "height":
        print_dual(f"Requested z    : {args.z} Ang", f_out)
    else:
        print_dual(f"Iso threshold  : {args.iso:.6e} e/Bohr^3", f_out)
        print_dual(f"z-max (search) : {z_max:.3f} Ang above the topmost atom", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Grid shape", f"{grid_data.shape[0]} x {grid_data.shape[1]} x {grid_data.shape[2]}"], None),
        (["LDOS range", f"[{grid_data.min():.6e}, {grid_data.max():.6e}] e/Bohr^3"], None),
        (["Surface-normal axis", f"{axis} ({axis_names[axis]})"], None),
        (["Vacuum-padded axes", str(vacuum_axes)], None),
        (["Topmost atom", f"{z_top:.3f} Ang along axis {axis}"], None),
        (["Axis length", f"{axis_len:.3f} Ang"], None),
    ], f_out)
    if axis_warning:
        print_dual(color_text(f"[WARNING] {axis_warning}", 'yellow'), f_out)
    else:
        print_dual("Requested/detected axis is Cartesian-aligned -- heights along it mean "
                   "what you'd expect.", f_out)

    print_section("[2] STM IMAGE", f_out)
    if args.mode == "height":
        print_table(["Quantity", "Value"], [
            (["Requested height", f"{image_stats['requested_z']:.3f} Ang"], None),
            (["Actual grid height used", f"{image_stats['actual_height']:.3f} Ang "
              f"(index {image_stats['idx']}/{n_axis})"], None),
            (["LDOS min", f"{image_stats['min']:.6e} e/Bohr^3"], None),
            (["LDOS max", f"{image_stats['max']:.6e} e/Bohr^3"], None),
            (["LDOS mean", f"{image_stats['mean']:.6e} e/Bohr^3"], None),
        ], f_out)
    else:
        rows = [
            (["Iso threshold", f"{image_stats['iso']:.6e} e/Bohr^3"], None),
            (["Points reaching iso", f"{image_stats['n_total'] - image_stats['n_missing']}/"
              f"{image_stats['n_total']}"], None),
        ]
        if image_stats["n_missing"]:
            pct = 100.0 * image_stats["n_missing"] / image_stats["n_total"]
            rows.append((["Points never reaching iso", f"{image_stats['n_missing']} ({pct:.1f}%)"], 'yellow'))
        if "min" in image_stats:
            rows.extend([
                (["Height min", f"{image_stats['min']:.3f} Ang"], None),
                (["Height max", f"{image_stats['max']:.3f} Ang"], None),
                (["Height mean", f"{image_stats['mean']:.3f} Ang"], None),
                (["Corrugation (max-min)", f"{image_stats['corrugation']:.3f} Ang"], 'green'),
            ])
        print_table(["Quantity", "Value"], rows, f_out)
        if image_stats["n_missing"]:
            pct = 100.0 * image_stats["n_missing"] / image_stats["n_total"]
            print_dual(color_text(
                f"[WARNING] {image_stats['n_missing']}/{image_stats['n_total']} points ({pct:.1f}%) "
                "never reached --iso within the search window -- consider a larger --z-max or a "
                "lower --iso. These points are written as NaN.", 'yellow'), f_out)

    print_section("[3] OUTPUT DATA & PLOTS", f_out)
    if args.save_gnuplot:
        write_stm_data(out_data, u_vals, v_vals, image, cb_label)
        write_stm_gplot(out_data, out_gplot, title, cb_label, u_label, v_label, is_signed=False)
        print_dual(color_text(f"[OK] Data written to '{out_data}'.", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{out_gplot}' "
                   "(run gnuplot from inside its own folder).", 'green'), f_out)
    else:
        print_dual(f"Not written (off by default -- pass --save-gnuplot to write "
                   f"stm_{args.mode}.dat/stm_{args.mode}.gplot).", f_out)

    print_section("[4] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.TERSOFF_HAMANN]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    if args.mode == "current" and "corrugation" in image_stats:
        print_dual(f"Corrugation    : {image_stats['corrugation']:.3f} Ang", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if args.save_gnuplot:
        print_dual(f"Data           : {out_data}", f_out)
        print_dual(f"Gnuplot script : {out_gplot}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        plot_matplotlib(u_vals, v_vals, image, u_label, v_label, cb_label, title)


if __name__ == "__main__":
    main()
