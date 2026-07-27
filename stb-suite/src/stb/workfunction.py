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
from stb.core import siesta_log
from stb.core import kspace
from stb.core import citations

from stb.core.deps import require_sisl, read_sisl_geometry_xv_or_fdf
sisl = require_sisl()

# ANSI Colors for terminal
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table

REPORT_FILE = "stb_workfunction_report.txt"
BIB_FILE = "references.bib"

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

# A plateau whose potential varies (in a systematic, monotonic sense -- not just
# scatter) by more than this across its own width is flagged as sloped rather than
# flat. A slope like this on an asymmetric slab is the classic symptom of a missing
# dipole correction (SIESTA's SlabDipoleCorrection): without it, the artificial
# periodic-image field shows up as a ramp across the vacuum instead of a plateau.
# Detected from the potential profile itself (a physical symptom) rather than by
# parsing the .out for whether the flag was set, since that keeps working regardless
# of exactly how a given SIESTA version logs it.
SLOPE_WARN_TOL = 0.05


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


def _unwrap_plateau_positions(indices, n):
    """Returns (original_indices, unwrapped_positions) for a plateau's index
    set, in a continuous monotonically-increasing order -- needed for any
    per-point analysis (line-fitting a slope, drawing a shaded z-range) that
    a wraparound plateau's indices would otherwise break. A wraparound
    plateau's indices, sorted numerically (e.g. [0,1,...,20,280,...,299] on
    a 300-point grid), have one large gap where they jump across the array
    boundary; positions before that gap are shifted up by `n` so the whole
    sequence becomes continuous (representing the true circular order:
    280,281,...,299,300(=0),301(=1),...,320(=20)), leaving positions
    possibly >= n -- callers convert to physical z via `position * point_width`
    and wrap back into [0, cell_length) themselves if needed.
    """
    idx = np.sort(np.asarray(indices))
    if len(idx) < 2:
        return idx, idx.astype(float)
    gaps = np.diff(idx)
    max_gap_pos = int(np.argmax(gaps))
    if gaps[max_gap_pos] <= 1:
        return idx, idx.astype(float)  # already contiguous, no wraparound
    unwrapped = idx.astype(float).copy()
    unwrapped[:max_gap_pos + 1] += n
    order = np.argsort(unwrapped)
    return idx[order], unwrapped[order]


def detect_plateau_slope(plateau, v_planar, n_points, point_width):
    """Fits a line to a plateau's potential vs. position (wraparound-aware
    via _unwrap_plateau_positions) and returns (slope_ev_per_ang,
    total_variation_ev) -- the latter being the fitted line's value change
    across the plateau's own width, which is what actually matters for
    deciding whether this looks like a real systematic ramp (see
    SLOPE_WARN_TOL) rather than random noise around a flat mean. Returns
    None if the plateau has too few points to fit meaningfully.
    """
    orig_idx, pos = _unwrap_plateau_positions(plateau['indices'], n_points)
    if len(orig_idx) < 3:
        return None
    z_pos = pos * point_width
    values = v_planar[orig_idx]
    slope, _intercept = np.polyfit(z_pos, values, 1)
    variation = float(slope * (z_pos[-1] - z_pos[0]))
    return float(slope), variation


def plateau_z_ranges(plateau, n_points, point_width, cell_length):
    """Returns a list of (z_start, z_end) tuples spanning a plateau's actual
    extent along the axis -- 2 tuples for a plateau that wraps around the
    cell boundary, 1 otherwise. For shading the detected region in a plot.
    """
    _orig_idx, pos = _unwrap_plateau_positions(plateau['indices'], n_points)
    z_start = pos[0] * point_width
    z_end = pos[-1] * point_width
    if z_end <= cell_length:
        return [(z_start, z_end)]
    return [(z_start, cell_length), (0.0, z_end - cell_length)]


def slab_and_vacuum_extent(geometry, axis, cell_length, plateaus, point_width):
    """Returns (slab_thickness_ang, vacuum_size_ang) along `axis` --
    vacuum_size is just the total length of every detected plateau (no
    geometry needed for that); slab_thickness is the span between the
    extreme atom positions along the axis, or None if no geometry is
    available. A quick sanity-check number for "is my vacuum gap actually
    generous:" one glance, not a substitute for checking the real structure.

    Wraparound-aware: a real slab is very commonly centered near z=0 (i.e.
    straddling the periodic boundary, e.g. fractional z of 0.93 and 0.07 for
    the two faces of a thin layer) -- a naive max()-min() on the raw wrapped
    fractional coordinates would then read this as spanning almost the whole
    cell instead of the true, narrow thickness. Fixed the same way
    _circular_mean_index/_unwrap_plateau_positions above handle the
    identical issue for a vacuum plateau: shift into a frame centered on the
    atoms' own circular mean position before taking the span. Verified live
    on a real CrS monolayer (2D-fetched structure, atoms at fractional
    z = 0, 0, 0.0659, 0.9341): the old naive calculation reported a
    slab thickness of ~21.8 Ang (the full cell minus the true ~1.5 Ang
    buckling), the fixed one correctly reports ~3.1 Ang.
    """
    vacuum_size = sum(p['size'] for p in plateaus) * point_width
    slab_thickness = None
    if geometry is not None:
        frac = geometry.fxyz[:, axis] % 1.0
        angles = 2 * np.pi * frac
        mean_frac = (np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles))) / (2 * np.pi)) % 1.0
        shifted = ((frac - mean_frac + 0.5) % 1.0 - 0.5) * cell_length
        slab_thickness = float(shifted.max() - shifted.min())
    return slab_thickness, vacuum_size


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


def write_gnuplot_wf(z_vals, v_planar, E_f, plateaus, label, output_dir):
    """Generates a .gplot file and a data file for plotting, both inside
    `output_dir`. Returns (data_path, gplot_path), or (None, None) if either
    write failed."""
    data_filename = os.path.join(output_dir, "workfunction_data.dat")
    gplot_filename = os.path.join(output_dir, "workfunction.gplot")

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
        return None, None

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
    # Relative to the .gplot file's own location (gnuplot is run from
    # output_dir), not the full joined path -- matches dos.py's own
    # dat_base/stem convention for the same reason.
    data_basename = os.path.basename(data_filename)
    fileout.append(f'plot "{data_basename}" using 1:2 with lines lw 3 lc rgb "navy" title "Planar Avg"\n')

    try:
        with open(gplot_filename, 'w') as file:
            file.writelines(fileout)
    except OSError as e:
        print(color_text(f"[ERROR] Could not write '{gplot_filename}': {e}", 'red'))
        return data_filename, None

    return data_filename, gplot_filename


def plot_matplotlib(z_vals, v_planar, E_f, plateaus, label, axis, cell_length):
    """Generates the Matplotlib preview."""
    try:
        n_points = len(z_vals)
        point_width = cell_length / n_points

        plt.figure(figsize=(8, 6))
        plt.plot(z_vals, v_planar, label='Planar Avg Potential', color='navy', linewidth=2)
        plt.axhline(y=E_f, color='forestgreen', linestyle='--', label=f'$E_F$ = {E_f:.2f}')

        for i, p in enumerate(plateaus):
            color = PLATEAU_COLORS[i % len(PLATEAU_COLORS)]["mpl"]
            plt.axhline(y=p['mean'], color=color, linestyle='--',
                        label=rf"$E_{{vac,{i+1}}}$ = {p['mean']:.2f} ($\Phi_{{{i+1}}}$={p['wf']:.2f})")
            plt.annotate('', xy=(p['z'], E_f), xytext=(p['z'], p['mean']),
                         arrowprops=dict(arrowstyle='<->', color=color, lw=1.5))
            for z0, z1 in plateau_z_ranges(p, n_points, point_width, cell_length):
                plt.axvspan(z0, z1, color=color, alpha=0.12, lw=0)

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
physically meaningless average, along with the vacuum-level difference
(proportional to the surface dipole moment). Also reports the slab
thickness/vacuum size along the axis (when geometry is available), and
warns if a "vacuum" plateau has a systematic slope (a common symptom of a
missing dipole correction) rather than being genuinely flat.""",
        epilog="Example usage:\n"
               "  stb-workfunction -l graphene\n"
               "  stb-workfunction -l slab --axis 2\n"
               "  stb-workfunction -l slab --fermi -4.5 --save-report --save-gnuplot",
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

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                         help="Directory to write references.bib into (and the report/data+gnuplot "
                              "files, with --save-report/--save-gnuplot) (default: current "
                              "directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                         help="Also write workfunction_data.dat and workfunction.gplot (the planar"
                              "-averaged potential profile, annotated with the detected vacuum "
                              "level(s)/work function(s)). Off by default.")
    parser.add_argument("--view", action="store_true",
                         help="Show an interactive matplotlib preview of the potential profile "
                              "before exiting. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    parser.add_argument("-v", "--version", action="version", version=f"stb-workfunction {VERSION}")

    args = parser.parse_args()

    if args.asymmetric_tol <= 0:
        parser.error("--asymmetric-tol must be positive.")

    # Filenames
    out_file = args.file if args.file else f"{args.label}.out"
    grid_file = args.grid if args.grid else f"{args.label}.VT"
    fermi_source = None

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Work Function",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("WORK FUNCTION CALCULATOR:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)

    # Remove stale workfunction_data.dat/workfunction.gplot left by an older
    # version of this tool (which always wrote them unconditionally to the
    # current directory -- there was no --output-dir/--save-gnuplot concept
    # before this migration) so they're never mistaken for current output.
    for stale_name in ("workfunction_data.dat", "workfunction.gplot"):
        if os.path.exists(stale_name):
            os.remove(stale_name)

    # --- 1. Fermi Energy ---
    print(f"\n[INFO] Detecting Fermi Energy...")
    E_f = args.fermi

    if E_f is None:
        if os.path.exists(out_file):
            E_f = siesta_log.get_fermi_energy(out_file)
            fermi_source = out_file
        else:
            print(f"{COLORS['yellow']}[WARNING] File {out_file} not found. Cannot auto-detect Fermi.{COLORS['reset']}")
    else:
        fermi_source = "--fermi (manual override)"

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
    axis_note = None
    if args.axis is None:
        if vacuum_axes is not None and sum(vacuum_axes) == 1:
            axis_note = f"auto-detected from {args.label}.XV/.fdf"
        else:
            axis_note = color_text(
                f"could not auto-detect a unique vacuum axis (no usable geometry file found, or "
                f"{'more than one' if vacuum_axes and sum(vacuum_axes) > 1 else 'no'} axis looks "
                f"vacuum-padded) -- defaulting to axis {axis}. Pass --axis explicitly if this is "
                "wrong.", 'yellow')
    elif vacuum_axes is not None and not vacuum_axes[axis]:
        axis_note = color_text(
            f"axis {axis} doesn't look vacuum-padded in {args.label}.XV/.fdf (no >= "
            f"{VACUUM_GAP_ANG} Ang gap found there) -- the work function on a genuinely periodic "
            "direction isn't physically meaningful. Double-check --axis.", 'yellow')

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

    point_width = cell_length / len(z_vals)
    plateau_notes = {}
    for p in plateaus:
        notes = []
        n_atoms_in = atoms_in_plateau(p, geometry, axis, len(z_vals))
        if n_atoms_in:
            notes.append(color_text(
                f"[WARNING] {n_atoms_in} atom(s) project into this vacuum region -- may not be "
                "genuine vacuum (e.g. an adsorbate) rather than a real, empty plateau.", 'yellow'))

        slope_result = detect_plateau_slope(p, v_planar, len(z_vals), point_width)
        if slope_result is not None:
            slope, variation = slope_result
            if abs(variation) > SLOPE_WARN_TOL:
                notes.append(color_text(
                    f"[WARNING] Systematic slope (~{variation:+.3f} eV across the plateau, "
                    f"{slope*1000:+.2f} meV/Ang) instead of a flat plateau -- a common symptom of "
                    "a missing dipole correction (SIESTA's SlabDipoleCorrection) or a real applied "
                    "field. Treat this region's vacuum level as approximate.", 'yellow'))
        if p['std'] > NOISY_PLATEAU_STD_TOL:
            notes.append(color_text(
                f"[WARNING] Noisy plateau (std={p['std']:.3f} eV) -- results might be inaccurate.",
                'yellow'))
        plateau_notes[id(p)] = notes

    slab_thickness, vacuum_size = slab_and_vacuum_extent(geometry, axis, cell_length, plateaus, point_width)

    # --- 6. Data + gnuplot (opt-in) ---
    data_path = gplot_path = None
    if args.save_gnuplot:
        data_path, gplot_path = write_gnuplot_wf(z_vals, v_planar, E_f, plateaus, args.label, args.output_dir)

    # --- 7. Report ---
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-WORKFUNCTION REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}", f_out)
    print_dual(f"Fermi source   : {fermi_source}", f_out)
    print_dual(f"Grid file      : {grid_file}", f_out)
    print_dual(f"Axis           : {axis}" + (f" ({axis_note})" if axis_note else ""), f_out)
    print_dual(f"Asymmetric tol : {args.asymmetric_tol} eV", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] FERMI ENERGY", f_out)
    print_dual(f"E_F = {E_f:.6f} eV (source: {fermi_source})", f_out)

    print_section("[2] VACUUM AXIS DETECTION", f_out)
    print_dual(f"Axis used: {axis} ({'x' if axis == 0 else 'y' if axis == 1 else 'z'})", f_out)
    if axis_note:
        print_dual(axis_note, f_out)
    else:
        print_dual("Requested axis matches the geometry's own detected vacuum-padded direction.", f_out)

    print_section("[3] VACUUM PLATEAU DETECTION", f_out)
    rows = []
    for i, p in enumerate(plateaus, start=1):
        rows.append(([str(i), f"{p['z']:.4f}", f"{p['mean']:.6f}", f"{p['std']:.6f}",
                      "yes" if p['std'] > NOISY_PLATEAU_STD_TOL else "no"], None))
    print_table(["Region", "z (Å)", "Mean potential (eV)", "Std (eV)", "Noisy?"], rows, f_out)
    for i, p in enumerate(plateaus, start=1):
        for note in plateau_notes[id(p)]:
            print_dual(f"Region {i}: {note}", f_out)

    print_section("[4] WORK FUNCTION RESULTS", f_out)
    print_dual(f"Vacuum size (axis {axis}) ~ {vacuum_size:.2f} Å"
               + (f"  |  Slab thickness ~ {slab_thickness:.2f} Å" if slab_thickness is not None else ""),
               f_out)
    if len(plateaus) == 1:
        rows = [(["Vacuum level", f"{plateaus[0]['mean']:.4f} eV"], None),
                (["Work Function", f"{plateaus[0]['wf']:.4f} eV"], 'green')]
        print_table(["Quantity", "Value"], rows, f_out)
    else:
        print_dual(color_text(
            f"{len(plateaus)} distinct vacuum levels detected -- this looks like an asymmetric "
            "slab (different terminations on each side); reporting each separately instead of "
            "one averaged (and physically meaningless) value.", 'yellow'), f_out)
        rows = []
        for i, p in enumerate(plateaus, start=1):
            rows.append(([f"Vacuum level {i} (z={p['z']:.2f} Å)", f"{p['mean']:.4f} eV"], None))
            rows.append(([f"Work Function {i}", f"{p['wf']:.4f} eV"], 'green'))
        rows.append((["Average Work Function", f"{np.mean([p['wf'] for p in plateaus]):.4f} eV"], 'green'))
        print_table(["Quantity", "Value"], rows, f_out)
        if len(plateaus) == 2:
            delta_v = plateaus[1]['mean'] - plateaus[0]['mean']
            print_dual(color_text(
                f"Vacuum level difference (dV) = {delta_v:+.4f} eV (proportional to the surface "
                "dipole moment, dV ~ 4*pi*mu/A -- reported directly rather than converting to an "
                "absolute dipole moment in e*Ang, since that needs a convention-specific prefactor; "
                "see e.g. Bengtsson, PRB 59, 12301 (1999)).", 'cyan'), f_out)

    print_section("[5] OUTPUT DATA & PLOTS", f_out)
    if data_path:
        print_dual(color_text(f"[OK] Data written to '{data_path}'.", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write "
                   "workfunction_data.dat/workfunction.gplot).", f_out)

    print_section("[6] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if data_path:
        print_dual(f"Data           : {data_path}", f_out)
        print_dual(f"Gnuplot script : {gplot_path}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        plot_matplotlib(z_vals, v_planar, E_f, plateaus, args.label, axis, cell_length)

if __name__ == "__main__":
    main()
