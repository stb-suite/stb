#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.19.0"

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from stb.core import siesta_log
from stb.core import kspace

from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
sisl = require_sisl()

# ANSI Colors for terminal
from stb.core.cli import COLORS, color_text, show_intro

VACUUM_GAP_ANG = 10.0  # min empty span (wrapped) along an axis to treat it as vacuum-padded --
                        # matches the same constant used by stb-kgrid/stb-mlrelax/etc.

# Fraction of grid points (by lowest |gradient|) considered for vacuum-plateau detection --
# same value the original single-region heuristic used.
VACUUM_PLATEAU_FRACTION = 0.15

# A genuine vacuum region should span at least this much (real vacuum gaps are usually
# 10+ Ang, see VACUUM_GAP_ANG above) -- a "flat" run narrower than this is more likely a
# brief near-zero-gradient moment at a local extremum of an otherwise-varying potential
# (e.g. deep inside a bulk/periodic region with no vacuum at all) than a real plateau.
# Expressed as a physical width, not a fixed point count, so it scales with grid
# resolution instead of becoming meaningless on a finer or coarser mesh.
MIN_PLATEAU_WIDTH_ANG = 2.0

# A short non-flat run (a few points of elevated gradient) inside an otherwise-flat
# region is bridged -- i.e. treated as flat too -- if it's narrower than this AND the
# potential level on both sides agrees (see MERGE_LEVEL_TOL). Real DFT data is noisier
# than a clean synthetic step function, and a couple of stray points can otherwise
# fragment one genuine, wide vacuum plateau into pieces each below
# MIN_PLATEAU_WIDTH_ANG. The level check is what keeps this from bridging an actual
# transition between two different potential levels, even a short/sharp one.
MAX_GAP_ANG = 1.0

# Bridging tolerance for MAX_GAP_ANG (above) and the tolerance for merging two detected
# plateaus that are the same physical vacuum level split by noise, vs. two GENUINELY
# distinct levels (e.g. an asymmetric slab's different terminations on each side).
MERGE_LEVEL_TOL = 0.05

# A per-plateau standard deviation above this is flagged as "noisy" in the report --
# independent from MERGE_LEVEL_TOL above even though they default to the same value:
# one decides whether two regions are the same physical plateau, the other just flags
# how much a single already-accepted plateau's potential wobbles internally.
NOISY_PLATEAU_STD_TOL = 0.05


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


def detect_vacuum_axis(label, requested_axis, grid_cell):
    """Reads <label>.XV or <label>.fdf (whichever exists) and reuses this
    suite's own core.kspace.detect_vacuum_axes() to check which lattice
    directions are vacuum-padded.

    `grid_cell` (the .VT grid's own lattice, already read by the caller) is
    cross-checked against the geometry file's cell -- if they disagree, the
    geometry almost certainly doesn't belong to the same calculation as the
    grid actually being planar-averaged, so it's ignored (treated the same
    as "no geometry file found") rather than trusted for axis detection.

    Returns (axis, vacuum_axes, geometry) where `axis` is `requested_axis`
    if given, otherwise the auto-detected single vacuum axis (falling back
    to 2/z if detection is unavailable or ambiguous); `vacuum_axes` is the
    raw length-3 bool list from detect_vacuum_axes, or None if no geometry
    was found/usable; `geometry` is the sisl Geometry if one was
    successfully read and matched the grid, else None. None of this is ever
    fatal -- it's an optional cross-check, not a hard requirement.
    """
    vacuum_axes = None
    geometry = None
    try:
        geometry, geom_file = read_sisl_geometry_xv_or_fdf(label)
    except Exception as e:
        print(color_text(
            f"   [WARNING] Found a geometry file for '{label}' but couldn't read it ({e}) "
            "-- ignoring it for axis auto-detection.", 'yellow'))
        geometry, geom_file = None, None

    if geometry is not None:
        if not np.allclose(geometry.cell, grid_cell, atol=1e-2):
            print(color_text(
                f"   [WARNING] {geom_file}'s cell doesn't match {label}.VT's -- they likely "
                "don't belong to the same calculation. Ignoring it for axis auto-detection.",
                'yellow'))
            geometry = None
        else:
            vacuum_axes = kspace.detect_vacuum_axes(geometry.fxyz, geometry.cell, VACUUM_GAP_ANG)

    if requested_axis is not None:
        return requested_axis, vacuum_axes, geometry

    if vacuum_axes is not None and sum(vacuum_axes) == 1:
        return vacuum_axes.index(True), vacuum_axes, geometry

    return 2, vacuum_axes, geometry


def atoms_in_plateau(plateau, geometry, axis, n_points):
    """Returns the number of atoms whose position projects into `plateau`
    along `axis` -- a genuine vacuum region should contain no atoms at all;
    if one does, this "plateau" may not be real vacuum (e.g. an adsorbate
    sitting alone on one side of a slab, creating its own locally-flat
    region) rather than a numerical artifact this tool can safely ignore.
    Returns 0 if no geometry is available (nothing to check against).
    """
    if geometry is None:
        return 0
    plateau_idx = set(int(i) for i in plateau['indices'])
    atom_idx = np.round((geometry.fxyz[:, axis] % 1.0) * n_points).astype(int) % n_points
    return sum(1 for i in atom_idx if int(i) in plateau_idx)


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

    Rotates the array to start at a guaranteed False index first, so the
    circular case reduces to a single ordinary linear scan -- no separate
    wraparound bookkeeping, and no risk of the same physical run coming back
    twice (as a spurious short "prefix" fragment plus the real run) the way
    a naive "double the array and de-duplicate by start index" approach can.
    """
    n = len(is_flat)
    if not is_flat.any():
        return []
    if is_flat.all():
        return [np.arange(n)]

    false_idx = int(np.where(~is_flat)[0][0])
    rotated = np.roll(is_flat, -false_idx)

    runs = []
    idx = 0
    while idx < n:
        if rotated[idx]:
            start = idx
            while idx < n and rotated[idx]:
                idx += 1
            runs.append((np.arange(start, idx) + false_idx) % n)
        else:
            idx += 1
    return runs


def _circular_mean_index(indices, n):
    """Mean position of a set of indices on a circular array of length n.
    A plain arithmetic mean is wrong whenever the indices wrap around the
    boundary -- e.g. mean([9, 0, 1]) = 10/3 = 3.3, nowhere near the actual
    tight cluster of points right at the wrap point. Computed the standard
    way for circular data: average the unit vectors at each index's angle,
    then take the angle of the resultant.
    """
    angles = 2 * np.pi * np.asarray(indices) / n
    mean_angle = np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
    return (mean_angle / (2 * np.pi)) % 1.0 * n


def _close_small_gaps(is_flat, v_planar, max_gap_points, level_tol):
    """Bridges (marks as flat) short runs of non-flat points whose
    surrounding flat levels agree within `level_tol` -- these are noise
    blips inside one physical plateau, not a real second region. A gap is
    left alone -- even a short one -- if the levels on either side of it
    genuinely differ, since that's a real transition (e.g. a sharp step
    between two different vacuum levels on an asymmetric slab), not noise.
    """
    n = len(is_flat)
    closed = is_flat.copy()
    if max_gap_points <= 0 or not is_flat.any() or is_flat.all():
        return closed
    for gap in _contiguous_runs(~is_flat):
        if len(gap) > max_gap_points:
            continue
        before_idx = (gap[0] - 1) % n
        after_idx = (gap[-1] + 1) % n
        if not is_flat[before_idx] or not is_flat[after_idx]:
            continue
        if abs(v_planar[before_idx] - v_planar[after_idx]) <= level_tol:
            closed[gap] = True
    return closed


def find_vacuum_plateaus(v_planar, cell_length, merge_tol=MERGE_LEVEL_TOL):
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
    flat for an instant. Small gaps are bridged (see _close_small_gaps)
    BEFORE the width filter runs, so realistic numerical noise can't
    fragment one genuinely wide plateau into pieces that each individually
    fail the width check.

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
    max_gap_points = max(0, int(np.floor(MAX_GAP_ANG / point_width)))
    is_flat = _close_small_gaps(is_flat, v_planar, max_gap_points, merge_tol)

    min_points = max(1, int(np.ceil(MIN_PLATEAU_WIDTH_ANG / point_width)))
    raw_runs = [run for run in _contiguous_runs(is_flat) if len(run) >= min_points]
    if not raw_runs:
        return []

    raw = [{'indices': run, 'mean': float(np.mean(v_planar[run]))} for run in raw_runs]

    # Group runs whose mean potential agrees within tolerance, sorted by mean first so
    # each run is only ever compared against its immediate neighbor in sorted order --
    # unlike comparing every run against a single group anchor, this correctly chains
    # together a plateau split into 3+ fragments by a gradual mean drift.
    raw.sort(key=lambda r: r['mean'])
    groups = [[raw[0]]]
    for r in raw[1:]:
        if abs(r['mean'] - groups[-1][-1]['mean']) <= merge_tol:
            groups[-1].append(r)
        else:
            groups.append([r])

    plateaus = []
    for group in groups:
        combined = np.unique(np.concatenate([g['indices'] for g in group]))
        plateaus.append({
            'indices': combined,
            'mean': float(np.mean(v_planar[combined])),
            'std': float(np.std(v_planar[combined])),
            'center_idx': _circular_mean_index(combined, n_points),
            'size': len(combined),
        })

    plateaus.sort(key=lambda p: p['center_idx'])
    return plateaus


PLATEAU_COLORS = [
    {"gnuplot": "red", "mpl": "firebrick"},
    {"gnuplot": "orange", "mpl": "darkorange"},
    {"gnuplot": "purple", "mpl": "purple"},
    {"gnuplot": "brown", "mpl": "saddlebrown"},
]


def annotate_plateaus(plateaus, z_vals, E_f):
    """Attaches the derived 'z' (position along the axis) and 'wf' (work
    function) values to each plateau dict once, so main()/write_gnuplot_wf()
    /plot_matplotlib() all read the same precomputed numbers instead of each
    re-deriving them from 'center_idx'.
    """
    for p in plateaus:
        p['z'] = z_vals[int(round(p['center_idx'])) % len(z_vals)]
        p['wf'] = p['mean'] - E_f
    return plateaus


def write_gnuplot_wf(z_vals, v_planar, E_f, plateaus, label):
    """Generates a .gplot file and a data file for plotting."""
    data_filename = "workfunction_data.dat"

    try:
        with open(data_filename, 'w') as f:
            f.write("# Distance(Ang)  Potential(eV)\n")
            f.write(f"# Fermi level (eV): {E_f:.6f}\n")
            for i, p in enumerate(plateaus, start=1):
                f.write(f"# Vacuum region {i}: E_vac={p['mean']:.6f} eV, "
                        f"WF={p['wf']:.6f} eV, near z={p['z']:.4f} Ang\n")
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

    for i, p in enumerate(plateaus):
        color = PLATEAU_COLORS[i % len(PLATEAU_COLORS)]["gnuplot"]
        y_offset = 0.15 + 0.05 * i
        fileout.append(f'set label "Evac{i+1} = {p["mean"]:.2f} eV" at graph 0.05, graph {y_offset:.2f} '
                        f'textcolor rgb "{color}"\n')
        fileout.append(f'set arrow from graph 0, first {p["mean"]} to graph 1, first {p["mean"]} '
                        f'nohead dt 2 lc rgb "{color}" lw 2\n')
        fileout.append(f'set arrow from {p["z"]}, {E_f} to {p["z"]}, {p["mean"]} heads lc rgb "black" lw 2\n')
        fileout.append(f'set label "WF{i+1} = {p["wf"]:.3f} eV" at {p["z"]}, {(E_f + p["mean"])/2} '
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

        for i, p in enumerate(plateaus):
            color = PLATEAU_COLORS[i % len(PLATEAU_COLORS)]["mpl"]
            plt.axhline(y=p['mean'], color=color, linestyle='--',
                        label=rf"$E_{{vac,{i+1}}}$ = {p['mean']:.2f} ($\Phi_{{{i+1}}}$={p['wf']:.2f})")
            plt.annotate('', xy=(p['z'], E_f), xytext=(p['z'], p['mean']),
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
If <label>.XV or <label>.fdf is present (and its cell matches <label>.VT's),
--axis is auto-detected from the structure's actual vacuum-padded direction
(core.kspace.detect_vacuum_axes, same mechanism as stb-kgrid/stb-mlrelax);
otherwise it defaults to z (2). An asymmetric slab (different terminations
on each side) genuinely has two different vacuum levels -- this is detected
automatically and reported as two separate work functions rather than one
physically meaningless average.""",
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
                              "label.XV/label.fdf if available and matching label.VT's cell, else z (2)")

    parser.add_argument("--fermi", type=float, help="Manually force Fermi Energy (eV).", default=None)
    parser.add_argument("--asymmetric-tol", type=float, default=MERGE_LEVEL_TOL, metavar="EV",
                         help=f"Vacuum levels within this many eV of each other are treated as the "
                              f"same physical plateau rather than a genuinely asymmetric slab "
                              f"(default: {MERGE_LEVEL_TOL})")
    parser.add_argument("--no-plot", action="store_true", help="Disable automatic plotting.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    parser.add_argument("-v", "--version", action="version", version=f"stb-workfunction {VERSION}")

    args = parser.parse_args()

    if args.asymmetric_tol <= 0:
        parser.error("--asymmetric-tol must be positive.")

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

    # --- 2. Grid Reading ---
    if not os.path.exists(grid_file):
        print(f"{COLORS['red']}[ERROR] Grid file '{grid_file}' not found.{COLORS['reset']}")
        sys.exit(1)

    print(f"[INFO] Reading potential from {grid_file}...")
    grid, potential_ev = read_grid_data(grid_file)

    # --- 3. Axis (cross-checked against the grid actually being read, not just the
    # geometry file in isolation) ---
    axis, vacuum_axes, geometry = detect_vacuum_axis(args.label, args.axis, grid.cell)
    if args.axis is None:
        if vacuum_axes is not None and sum(vacuum_axes) == 1:
            print(f"[INFO] Auto-detected vacuum axis: {axis} (from {args.label}.XV/.fdf)")
        else:
            print(color_text(
                f"   [WARNING] Could not auto-detect a unique vacuum axis (no usable geometry "
                f"file found, or {'more than one' if vacuum_axes and sum(vacuum_axes) > 1 else 'no'} "
                f"axis looks vacuum-padded) -- defaulting to axis {axis}. Pass --axis explicitly "
                "if this is wrong.", 'yellow'))
    elif vacuum_axes is not None and not vacuum_axes[axis]:
        print(color_text(
            f"   [WARNING] Axis {axis} doesn't look vacuum-padded in {args.label}.XV/.fdf "
            f"(no >= {VACUUM_GAP_ANG} Ang gap found there) -- the work function on a genuinely "
            "periodic direction isn't physically meaningful. Double-check --axis.", 'yellow'))

    # --- 4. Planar Average ---
    print(f"[INFO] Calculating planar average along axis {axis}...")
    z_vals, v_planar, cell_length = calculate_planar_avg(grid, potential_ev, axis)

    # --- 5. Vacuum Level(s) & Work Function(s) ---
    plateaus = find_vacuum_plateaus(v_planar, cell_length, merge_tol=args.asymmetric_tol)
    if not plateaus:
        print("\n" + color_text(
            f"FATAL ERROR: No vacuum plateau found along axis {axis} -- this direction may not "
            "actually have vacuum. Try a different --axis, or check the structure.", 'red'))
        sys.exit(1)
    annotate_plateaus(plateaus, z_vals, E_f)

    for p in plateaus:
        n_atoms_in = atoms_in_plateau(p, geometry, axis, len(z_vals))
        if n_atoms_in:
            print(color_text(
                f"   [WARNING] {n_atoms_in} atom(s) project into the vacuum region near "
                f"z={p['z']:.2f} Ang -- this may not be genuine vacuum (e.g. an adsorbate) "
                "rather than a real, empty plateau.", 'yellow'))

    # --- 6. Reporting ---
    print("-" * 40)
    print(color_text(f"RESULTS for {args.label}:", 'cyan'))
    print(f"  Fermi Level    = {E_f:8.4f} eV")
    if len(plateaus) == 1:
        print(f"  Vacuum Level   = {plateaus[0]['mean']:8.4f} eV")
        print(color_text(f"  Work Function  = {plateaus[0]['wf']:8.4f} eV", 'green'))
        if plateaus[0]['std'] > NOISY_PLATEAU_STD_TOL:
            print(f"{COLORS['yellow']}[WARNING] Vacuum plateau is noisy (std={plateaus[0]['std']:.3f} "
                  f"eV). Results might be inaccurate.{COLORS['reset']}")
    else:
        print(color_text(
            f"  {len(plateaus)} distinct vacuum levels detected -- this looks like an "
            "asymmetric slab (different terminations on each side); reporting each "
            "separately instead of one averaged (and physically meaningless) value.",
            'yellow'))
        for i, p in enumerate(plateaus, start=1):
            print(f"  Vacuum Level {i} (near z={p['z']:6.2f} Ang) = {p['mean']:8.4f} eV")
            print(color_text(f"  Work Function {i}                    = {p['wf']:8.4f} eV", 'green'))
            if p['std'] > NOISY_PLATEAU_STD_TOL:
                print(f"{COLORS['yellow']}[WARNING] Vacuum region {i} is noisy (std={p['std']:.3f} "
                      f"eV). Results might be inaccurate.{COLORS['reset']}")
        print(f"  Average Work Function              = {np.mean([p['wf'] for p in plateaus]):8.4f} eV")
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
