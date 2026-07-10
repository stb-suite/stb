#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.18.0"

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from stb.core import siesta_log
from stb.core import kspace

from stb.core.deps import require_sisl
sisl = require_sisl()

# ANSI Colors for terminal
from stb.core.cli import COLORS, color_text, show_intro

VACUUM_GAP_ANG = 10.0  # min empty span (wrapped) along an axis to treat it as vacuum-padded --
                        # matches the same constant used by stb-kgrid/stb-mlrelax/etc.

# Fraction of grid points (by lowest |gradient|) considered for vacuum-plateau detection --
# same value the original single-region heuristic used.
VACUUM_PLATEAU_FRACTION = 0.15

# A genuine vacuum region should span at least this much (real vacuum gaps are usually
# 10+ Ang, see VACUUM_GAP_ANG below) -- a "flat" run narrower than this is more likely a
# brief near-zero-gradient moment at a local extremum of an otherwise-varying potential
# (e.g. deep inside a bulk/periodic region with no vacuum at all) than a real plateau.
# Expressed as a physical width, not a fixed point count, so it scales with grid
# resolution instead of becoming meaningless on a finer or coarser mesh.
MIN_PLATEAU_WIDTH_ANG = 2.0

# Two detected plateaus whose mean potential differs by more than this are treated as
# genuinely distinct vacuum levels (e.g. an asymmetric slab with different terminations on
# each side) rather than the same region split by a small numerical artifact.
ASYMMETRIC_SLAB_TOL = 0.05


def read_grid_data(grid_file):
    """Reads the grid using sisl (.VT grids are already returned in eV by sisl)."""
    try:
        grid = sisl.get_sile(grid_file).read_grid()
        potential_ev = grid.grid
        return grid, potential_ev
    except Exception as e:
        print(f"{COLORS['red']}[ERROR] Reading grid: {e}{COLORS['reset']}")
        sys.exit(1)


def calculate_planar_avg(grid, potential_ev, axis):
    """Calculates the planar average along the specified axis."""
    avg_axes = tuple(a for a in (0, 1, 2) if a != axis)
    L = grid.lattice.length[axis]
    N = grid.shape[axis]

    v_planar = np.mean(potential_ev, axis=avg_axes)
    # endpoint=False: the grid is periodic, so point N would just be point 0 again --
    # including it as a distinct, distance-L coordinate would misrepresent the spacing.
    z_vals = np.linspace(0, L, N, endpoint=False)
    return z_vals, v_planar, L


def detect_vacuum_axis(label, requested_axis):
    """Reads <label>.XV or <label>.fdf (whichever exists) and reuses this
    suite's own core.kspace.detect_vacuum_axes() to check which lattice
    directions are vacuum-padded.

    Returns (axis, vacuum_axes) where `axis` is `requested_axis` if given,
    otherwise the auto-detected single vacuum axis (falling back to 2/z if
    detection is unavailable or ambiguous); `vacuum_axes` is the raw length
    -3 bool list from detect_vacuum_axes, or None if no geometry file was
    found or reading/detection failed -- callers use it to warn if the
    chosen axis doesn't actually look vacuum-padded, but this is never fatal
    since it's an optional cross-check, not a hard requirement.
    """
    file_xv = f"{label}.XV"
    file_fdf = f"{label}.fdf"
    geom_file = file_xv if os.path.exists(file_xv) else (file_fdf if os.path.exists(file_fdf) else None)

    vacuum_axes = None
    if geom_file is not None:
        try:
            geometry = sisl.get_sile(geom_file).read_geometry()
            vacuum_axes = kspace.detect_vacuum_axes(geometry.fxyz, geometry.cell, VACUUM_GAP_ANG)
        except Exception:
            vacuum_axes = None

    if requested_axis is not None:
        return requested_axis, vacuum_axes

    if vacuum_axes is not None and sum(vacuum_axes) == 1:
        return vacuum_axes.index(True), vacuum_axes

    return 2, vacuum_axes


def _periodic_gradient(v):
    """Centered difference using periodic wrap-around. Plain np.gradient()
    doesn't know the grid is periodic and uses one-sided differences at the
    array edges, which can spuriously mark a genuine vacuum plateau that
    straddles the cell boundary -- the single most common slab setup, where
    the vacuum gap sits at the top/bottom of the cell -- as "not flat" right
    where it wraps around.
    """
    return (np.roll(v, -1) - np.roll(v, 1)) / 2.0


def _contiguous_runs(is_flat):
    """Returns a list of index arrays, each a maximal contiguous run of True
    values in `is_flat`, treating the array as circular (the grid wraps from
    the last index back to the first).
    """
    n = len(is_flat)
    if not is_flat.any():
        return []
    if is_flat.all():
        return [np.arange(n)]

    doubled = np.concatenate([is_flat, is_flat])
    diff = np.diff(doubled.astype(int))
    run_starts = list(np.where(diff == 1)[0] + 1)
    run_ends = list(np.where(diff == -1)[0] + 1)
    if doubled[0]:
        run_starts = [0] + run_starts
    if doubled[-1]:
        run_ends = run_ends + [len(doubled)]

    runs = []
    seen_starts = set()
    for s, e in zip(run_starts, run_ends):
        if s >= n:
            continue  # a repeat of a run already captured from the first copy
        start_mod = s % n
        if start_mod in seen_starts:
            continue
        seen_starts.add(start_mod)
        length = min(e - s, n)
        runs.append((np.arange(s, s + length) % n))
    return runs


def find_vacuum_plateaus(v_planar, cell_length):
    """Finds contiguous (periodic-boundary-aware) low-gradient regions in the
    planar-averaged potential -- each one a candidate vacuum plateau.

    Unlike picking the flattest points anywhere in the cell regardless of
    where they are, requiring them to form a connected region means two
    physically different flat regions (e.g. the two vacuum gaps on either
    side of an asymmetric slab, sitting at different potential levels, or a
    coincidentally-flat spot inside the material itself) can't get silently
    averaged together into a physically meaningless number. A minimum
    physical width (MIN_PLATEAU_WIDTH_ANG) further rejects a merely brief
    near-zero-gradient moment near a local extremum of an otherwise-varying
    bulk-like potential -- a real vacuum plateau is wide, not just locally
    flat for an instant.

    Returns a list of dicts (one per distinct plateau, sorted by position
    along the axis), each with 'indices', 'mean', 'std', 'center_idx',
    'size' -- or an empty list if no genuine plateau is found at all (e.g.
    the chosen axis has no real vacuum region to speak of).
    """
    grad = np.abs(_periodic_gradient(v_planar))
    n_points = len(grad)
    n_flat = max(1, int(n_points * VACUUM_PLATEAU_FRACTION))
    threshold = np.partition(grad, n_flat - 1)[n_flat - 1]
    is_flat = grad <= threshold

    point_width = cell_length / n_points
    min_points = max(1, int(np.ceil(MIN_PLATEAU_WIDTH_ANG / point_width)))
    raw_runs = [run for run in _contiguous_runs(is_flat) if len(run) >= min_points]
    if not raw_runs:
        return []

    raw = [{'indices': run, 'mean': float(np.mean(v_planar[run])),
            'center_idx': float(np.mean(run))} for run in raw_runs]

    # Merge runs whose mean potential agrees within tolerance -- the same physical
    # plateau can get split into two runs by a small numerical artifact right at
    # its edge (see _periodic_gradient's docstring), which isn't evidence of two
    # distinct vacuum levels.
    merged = []
    used = [False] * len(raw)
    for i, p in enumerate(raw):
        if used[i]:
            continue
        group_idx = [i]
        used[i] = True
        for j in range(i + 1, len(raw)):
            if not used[j] and abs(raw[j]['mean'] - p['mean']) <= ASYMMETRIC_SLAB_TOL:
                group_idx.append(j)
                used[j] = True
        combined = np.concatenate([raw[k]['indices'] for k in group_idx])
        merged.append({
            'indices': combined,
            'mean': float(np.mean(v_planar[combined])),
            'std': float(np.std(v_planar[combined])),
            'center_idx': float(np.mean([raw[k]['center_idx'] for k in group_idx])),
            'size': len(combined),
        })

    merged.sort(key=lambda p: p['center_idx'])
    return merged


def write_gnuplot_wf(z_vals, v_planar, E_f, plateaus, label):
    """Generates a .gplot file and a data file for plotting."""
    data_filename = "workfunction_data.dat"

    try:
        with open(data_filename, 'w') as f:
            f.write("# Distance(Ang)  Potential(eV)\n")
            f.write(f"# Fermi level (eV): {E_f:.6f}\n")
            for i, p in enumerate(plateaus, start=1):
                z_c = z_vals[int(round(p['center_idx'])) % len(z_vals)]
                f.write(f"# Vacuum region {i}: E_vac={p['mean']:.6f} eV, "
                        f"WF={p['mean'] - E_f:.6f} eV, near z={z_c:.4f} Ang\n")
            for z, v in zip(z_vals, v_planar):
                f.write(f"{z:.6f}     {v:.6f}\n")
    except OSError as e:
        print(color_text(f"[ERROR] Could not write '{data_filename}': {e}", 'red'))
        return

    fileout = []
    fileout.append('# Set terminal and output\n')
    fileout.append('set terminal pdfcairo enhanced font "Arial,20" size 8,6\n')
    fileout.append(f'set output "{label}_WF.pdf"\n')
    fileout.append('\n')
    fileout.append('# Labels and Styles\n')
    fileout.append('set ylabel "Potential (eV)" font "Arial,22"\n')
    fileout.append('set xlabel "Position (Angstrom)" font "Arial,22"\n')
    fileout.append('set grid xtics ytics lt 0 lw 1 lc rgb "#bbbbbb"\n')
    fileout.append(f'set title "Work Function: {label}"\n')
    fileout.append('\n')

    fileout.append('# Fermi level\n')
    fileout.append(f'set label "Ef = {E_f:.2f} eV" at graph 0.05, graph 0.1 textcolor rgb "forest-green"\n')
    fileout.append(f'set arrow from graph 0, first {E_f} to graph 1, first {E_f} nohead dt 2 '
                    'lc rgb "forest-green" lw 2\n')

    colors = ["red", "orange", "purple", "brown"]
    for i, p in enumerate(plateaus):
        color = colors[i % len(colors)]
        y_offset = 0.15 + 0.05 * i
        wf_val = p['mean'] - E_f
        z_c = z_vals[int(round(p['center_idx'])) % len(z_vals)]
        fileout.append(f'set label "Evac{i+1} = {p["mean"]:.2f} eV" at graph 0.05, graph {y_offset:.2f} '
                        f'textcolor rgb "{color}"\n')
        fileout.append(f'set arrow from graph 0, first {p["mean"]} to graph 1, first {p["mean"]} '
                        f'nohead dt 2 lc rgb "{color}" lw 2\n')
        fileout.append(f'set arrow from {z_c}, {E_f} to {z_c}, {p["mean"]} heads lc rgb "black" lw 2\n')
        fileout.append(f'set label "WF{i+1} = {wf_val:.3f} eV" at {z_c}, {(E_f + p["mean"])/2} '
                        'center font ",14"\n')

    fileout.append('\n')
    fileout.append(f'plot "{data_filename}" using 1:2 with lines lw 3 lc rgb "navy" title "Planar Avg"\n')

    try:
        with open('workfunction.gplot', 'w') as file:
            file.writelines(fileout)
    except OSError as e:
        print(color_text(f"[ERROR] Could not write 'workfunction.gplot': {e}", 'red'))
        return

    print(f"[INFO] Gnuplot script saved to 'workfunction.gplot'")
    print(f"[INFO] Data saved to '{data_filename}'")


def plot_matplotlib(z_vals, v_planar, E_f, plateaus, label, axis):
    """Generates the Matplotlib preview."""
    try:
        plt.figure(figsize=(8, 6))
        plt.plot(z_vals, v_planar, label='Planar Avg Potential', color='navy', linewidth=2)
        plt.axhline(y=E_f, color='forestgreen', linestyle='--', label=f'$E_F$ = {E_f:.2f}')

        colors = ["firebrick", "darkorange", "purple", "saddlebrown"]
        for i, p in enumerate(plateaus):
            color = colors[i % len(colors)]
            wf_val = p['mean'] - E_f
            z_c = z_vals[int(round(p['center_idx'])) % len(z_vals)]
            plt.axhline(y=p['mean'], color=color, linestyle='--',
                        label=rf"$E_{{vac,{i+1}}}$ = {p['mean']:.2f} ($\Phi_{{{i+1}}}$={wf_val:.2f})")
            plt.annotate('', xy=(z_c, E_f), xytext=(z_c, p['mean']),
                         arrowprops=dict(arrowstyle='<->', color=color, lw=1.5))

        plt.xlabel(rf"Position along Axis {axis} ($\AA$)", fontsize=14)
        plt.ylabel("Potential (eV)", fontsize=14)
        plt.title(f"Work Function: {label}", fontsize=16)
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"{COLORS['yellow']}[WARNING] Could not create interactive plot: {e}{COLORS['reset']}")


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Calculates the work function from a SIESTA electrostatic potential grid.", 'bold')}
Reads <label>.VT (total potential, planar-averaged along --axis) plus the
Fermi energy from <label>.out (or --fermi to override), and locates the
vacuum plateau(s) in the averaged potential to compute Phi = E_vac - E_F.
If <label>.XV or <label>.fdf is present, --axis is auto-detected from the
structure's actual vacuum-padded direction (core.kspace.detect_vacuum_axes,
same mechanism as stb-kgrid/stb-mlrelax); otherwise it defaults to z (2).
An asymmetric slab (different terminations on each side) genuinely has two
different vacuum levels -- this is detected automatically and reported as
two separate work functions rather than one physically meaningless average.""",
        epilog="Example usage:\n"
               "  stb-workfunction -l graphene\n"
               "  stb-workfunction -l slab --axis 2\n"
               "  stb-workfunction -l slab --fermi -4.5 --no-plot",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", help="System Label (required).", required=True)
    parser.add_argument("-f", "--file", help="Output file (.out). Default: label.out", default=None)
    parser.add_argument("-g", "--grid", help="Potential grid file (.VT recommended). Default: label.VT", default=None)
    parser.add_argument("-z", "--axis", type=int, choices=[0, 1, 2], default=None,
                         help="Axis normal to surface (0=x, 1=y, 2=z). Default: auto-detected from "
                              "label.XV/label.fdf if available, else z (2)")

    parser.add_argument("--fermi", type=float, help="Manually force Fermi Energy (eV).", default=None)
    parser.add_argument("--no-plot", action="store_true", help="Disable automatic plotting.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    parser.add_argument("-v", "--version", action="version", version=f"stb-workfunction {VERSION}")

    args = parser.parse_args()

    # Filenames
    out_file = args.file if args.file else f"{args.label}.out"
    grid_file = args.grid if args.grid else f"{args.label}.VT"

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Work Function",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("WORK FUNCTION CALCULATOR:", 'bold'))
    print("-" * 60)

    # --- 1. Fermi Energy ---
    print(f"\n[INFO] Detecting Fermi Energy...")
    E_f = args.fermi

    if E_f is None:
        if os.path.exists(out_file):
            E_f = siesta_log.get_fermi_energy(out_file)
        else:
            print(f"{COLORS['yellow']}[WARNING] File {out_file} not found. Cannot auto-detect Fermi.{COLORS['reset']}")

    if E_f is None:
        print("\n" + color_text("FATAL ERROR: Could not determine Fermi Energy.", 'red'))
        print("Please run with manual value found in your logs:")
        print(f"  python {sys.argv[0]} -l {args.label} --fermi -3.635")
        sys.exit(1)

    print(f"[INFO] Fermi Energy (Ef): {E_f:.6f} eV")

    # --- 2. Axis ---
    axis, vacuum_axes = detect_vacuum_axis(args.label, args.axis)
    if args.axis is None:
        if vacuum_axes is not None and sum(vacuum_axes) == 1:
            print(f"[INFO] Auto-detected vacuum axis: {axis} (from {args.label}.XV/.fdf)")
        else:
            print(color_text(
                f"   [WARNING] Could not auto-detect a unique vacuum axis (no geometry file "
                f"found, or {'more than one' if vacuum_axes and sum(vacuum_axes) > 1 else 'no'} "
                f"axis looks vacuum-padded) -- defaulting to axis {axis}. Pass --axis explicitly "
                "if this is wrong.", 'yellow'))
    elif vacuum_axes is not None and not vacuum_axes[axis]:
        print(color_text(
            f"   [WARNING] Axis {axis} doesn't look vacuum-padded in {args.label}.XV/.fdf "
            f"(no >= {VACUUM_GAP_ANG} Ang gap found there) -- the work function on a genuinely "
            "periodic direction isn't physically meaningful. Double-check --axis.", 'yellow'))

    # --- 3. Grid Reading ---
    if not os.path.exists(grid_file):
        print(f"{COLORS['red']}[ERROR] Grid file '{grid_file}' not found.{COLORS['reset']}")
        sys.exit(1)

    print(f"[INFO] Reading potential from {grid_file}...")
    grid, potential_ev = read_grid_data(grid_file)

    # --- 4. Planar Average ---
    print(f"[INFO] Calculating planar average along axis {axis}...")
    z_vals, v_planar, cell_length = calculate_planar_avg(grid, potential_ev, axis)

    # --- 5. Vacuum Level(s) & Work Function(s) ---
    plateaus = find_vacuum_plateaus(v_planar, cell_length)
    if not plateaus:
        print("\n" + color_text(
            f"FATAL ERROR: No vacuum plateau found along axis {axis} -- this direction may not "
            "actually have vacuum. Try a different --axis, or check the structure.", 'red'))
        sys.exit(1)

    # --- 6. Reporting ---
    print("-" * 40)
    print(color_text(f"RESULTS for {args.label}:", 'cyan'))
    print(f"  Fermi Level    = {E_f:8.4f} eV")
    if len(plateaus) == 1:
        WF = plateaus[0]['mean'] - E_f
        print(f"  Vacuum Level   = {plateaus[0]['mean']:8.4f} eV")
        print(color_text(f"  Work Function  = {WF:8.4f} eV", 'green'))
        if plateaus[0]['std'] > 0.05:
            print(f"{COLORS['yellow']}[WARNING] Vacuum plateau is noisy (std={plateaus[0]['std']:.3f} "
                  f"eV). Results might be inaccurate.{COLORS['reset']}")
    else:
        print(color_text(
            f"  {len(plateaus)} distinct vacuum levels detected -- this looks like an "
            "asymmetric slab (different terminations on each side); reporting each "
            "separately instead of one averaged (and physically meaningless) value.",
            'yellow'))
        wf_values = []
        for i, p in enumerate(plateaus, start=1):
            z_c = z_vals[int(round(p['center_idx'])) % len(z_vals)]
            wf = p['mean'] - E_f
            wf_values.append(wf)
            print(f"  Vacuum Level {i} (near z={z_c:6.2f} Ang) = {p['mean']:8.4f} eV")
            print(color_text(f"  Work Function {i}                    = {wf:8.4f} eV", 'green'))
            if p['std'] > 0.05:
                print(f"{COLORS['yellow']}[WARNING] Vacuum region {i} is noisy (std={p['std']:.3f} "
                      f"eV).{COLORS['reset']}")
        print(f"  Average Work Function              = {np.mean(wf_values):8.4f} eV")
    print("-" * 40)

    # --- 7. Output Files & Plotting ---
    print(f"[INFO] Writing output files...")
    write_gnuplot_wf(z_vals, v_planar, E_f, plateaus, args.label)

    if not args.no_plot:
        plot_matplotlib(z_vals, v_planar, E_f, plateaus, args.label, axis)

    print("\n[INFO] Complete job!")
    print("\n"+"-"*60)
    print(color_text("Work Function calculated. Now I need a vacation.\n\n", 'bold'))

if __name__ == "__main__":
    main()
