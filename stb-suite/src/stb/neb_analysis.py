#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import time
import shutil
import argparse
from collections import Counter
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import siesta_log
from stb.core import structure_io
from stb.core.siesta_log import check_scf_and_force
from stb.core.cli import color_text, show_intro, print_dual, print_section, capture_library_noise
from stb.core.ase_view import view_structure_interactive
from stb.neb import MODE_DESCRIPTIONS, write_path_trajectory

REPORT_FILE = "neb_report.txt"
SETUP_REPORT_FILE = "neb_setup.txt"
_CYCLE_DIR_RE = re.compile(r"^cycle_(\d+)$")
_TOTAL_PATH_LENGTH_RE = re.compile(r"Total path length\s*:\s*([0-9.]+)\s*Ang")


def read_image_table(root_dir):
    """Parses the '# IMAGE_TABLE' section and the '# ML_NEB_USED: yes|no'
    marker stb-neb writes in neb_setup.txt. Returns (rows, ml_neb_used)
    where rows is a list of (label, index, reaction_coord) sorted by
    index, or (None, False) if the report (or its table section) isn't
    found -- the caller then falls back to a sorted glob of 'image_*'
    folders. Mirrors adsorb_analysis.py::read_site_table -- the image
    directory is recomputed by the caller as os.path.join(root_dir,
    label), not trusted from the table's own 'dir' column (written
    relative to stb-neb's --output-dir at prep time, which may not be the
    same as --dir here).
    """
    report_path = os.path.join(root_dir, SETUP_REPORT_FILE)
    if not os.path.isfile(report_path):
        return None, False
    rows = []
    ml_neb_used = False
    in_table = False
    with open(report_path) as f:
        for line in f:
            if line.startswith("# ML_NEB_USED:"):
                ml_neb_used = line.split(":", 1)[1].strip().lower() == "yes"
                continue
            if line.startswith("# IMAGE_TABLE"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            label, index_str, coord_str = parts[0], parts[1], parts[2]
            try:
                index = int(index_str)
                reaction_coord = float(coord_str)
            except ValueError:
                continue
            rows.append((label, index, reaction_coord))
    if not rows:
        return None, False
    rows.sort(key=lambda r: r[1])
    return rows, ml_neb_used


def read_mode(root_dir):
    """Reads '# MODE: N' from SETUP_REPORT_FILE (stb-neb v2's mode marker,
    see neb.py's own writer). Returns 1 (the original, pre-mode single-
    point behavior) if no marker is found -- either an old neb_setup.txt
    from before this workflow gained modes, or no report at all.
    """
    report_path = os.path.join(root_dir, SETUP_REPORT_FILE)
    if os.path.isfile(report_path):
        with open(report_path) as f:
            for line in f:
                if line.startswith("# MODE:"):
                    try:
                        return int(line.split(":", 1)[1].strip())
                    except ValueError:
                        break
    return 1


def find_analysis_cycle(root_dir, out_filename):
    """For modes 2/3: returns (cycle_dir, converged). (None, False) if no
    'cycle_NN' exists yet (stb-neb --mode 2/3 hasn't been run, or SIESTA
    hasn't finished cycle_00 yet).

    Once NEB_CONVERGED exists, stb-nebCycle has stopped writing new
    cycles, so the highest-numbered cycle_NN IS the converged one --
    plain `find_latest_cycle` is correct there. But while NOT YET
    CONVERGED, the highest-numbered cycle_NN is stb-nebCycle's most
    recent WRITE TARGET (the next geometry it wants SIESTA to run, built
    by reading the previous cycle's results and writing cycle_{N+1}/,
    see neb_cycle.py's own `next_dir = ... cycle_{cycle_num + 1:02d}`) --
    it has NO calc.out at all yet by construction, live-reproduced: right
    after any stb-nebCycle call, the latest folder is always one step
    AHEAD of the last cycle SIESTA actually finished. Picking it
    unconditionally (the previous behavior here) meant analyzing a cycle
    with zero results every single time between two stb-nebCycle calls,
    surfacing as a hard "No valid image results found" error instead of
    the previous cycle's real numbers. Fixed by scanning cycle_NN in
    descending order and picking the highest-numbered one that actually
    has at least one image's out_filename on disk (a cycle only needs
    ONE image with results to be worth reporting -- partial-completion
    tolerance already exists per-image below, via the SKIP/[WARNING]
    machinery).
    """
    from stb.neb_cycle import find_latest_cycle, CONVERGED_SENTINEL
    converged = os.path.isfile(os.path.join(root_dir, CONVERGED_SENTINEL))
    if not converged:
        candidates = sorted(
            (m for m in (_CYCLE_DIR_RE.match(e) for e in os.listdir(root_dir))
             if m and os.path.isdir(os.path.join(root_dir, m.group(0)))),
            key=lambda m: int(m.group(1)), reverse=True)
        for m in candidates:
            cycle_dir = os.path.join(root_dir, m.group(0))
            has_results = any(
                os.path.isfile(os.path.join(cycle_dir, d, out_filename))
                for d in os.listdir(cycle_dir)
                if os.path.isdir(os.path.join(cycle_dir, d)) and d.startswith("image_"))
            if has_results:
                return cycle_dir, converged
        # No cycle has any results at all (e.g. cycle_00 fresh out of
        # stb-neb, before SIESTA/stb-nebCycle has ever run) -- fall
        # through to plain latest-folder behavior so the existing
        # downstream "No valid image results found" error still fires
        # with its normal wording.
    _cycle_num, cycle_dir = find_latest_cycle(root_dir)
    return cycle_dir, converged


class ImageRow:
    """One analyzed image: label, band index, reaction coordinate (Ang,
    cumulative displacement from image_00 -- None if the table fallback
    was used), and its computed energy/quality diagnostics. Plain
    attribute-holder, same style as adsorb_analysis.py's SiteRow.
    """
    def __init__(self, label, index, reaction_coord, energy, scf_ok, max_force):
        self.label = label
        self.index = index
        self.reaction_coord = reaction_coord
        self.energy = energy
        self.scf_ok = scf_ok
        self.max_force = max_force


def fit_energy_spline(rows):
    """Fits a cubic spline (scipy's default not-a-knot boundary condition,
    not a "natural" spline -- not-a-knot avoids forcing zero curvature at
    the two endpoints, which a natural spline would, biasing the fitted
    curve near the edges) through the (reaction_coord, energy) points,
    sorted by reaction_coord. Returns (x, y, spline) or None if there are
    fewer than 4 points (not enough for a meaningful cubic fit) or two
    images share the same reaction coordinate (a degenerate/duplicate
    image -- see stb-neb's check_path_quality). Factored out as its own
    function so the barrier-fitting search (find_spline_barrier) and the
    plotting code (build_energy_profile_figure, write_curve_plot's
    smooth-line decision) all reuse the SAME fitted spline instead of
    each re-fitting it independently.
    """
    if len(rows) < 4:
        return None
    rows_sorted = sorted(rows, key=lambda r: r.reaction_coord)
    x = [r.reaction_coord for r in rows_sorted]
    y = [r.energy for r in rows_sorted]
    if len(set(x)) != len(x):
        return None
    from scipy.interpolate import CubicSpline
    spline = CubicSpline(x, y)
    return x, y, spline


def find_spline_barrier(energy_spline):
    """Given fit_energy_spline's (x, y, spline) tuple (or None), returns
    (ts_coord, spline_barrier, spline_dE) -- a smoother, less grid-
    dependent transition-state estimate than "whichever discrete image
    happens to be highest", the DFT-side analogue of
    ase.mep.neb.NEBTools.get_barrier(fit=True) (already used by stb-neb
    for the MACE band -- that one additionally has per-image forces for
    a Hermite-style fit; here we only have a scalar SIESTA single-point
    energy per image, so a plain energy-only cubic spline is the fair
    equivalent). Returns None if energy_spline is None.

    Finding the maximum is deliberately NOT a single scipy.optimize.
    minimize_scalar(bounds=(x[0], x[-1])) call: Brent's method (its
    'bounded' method) converges to the FIRST local maximum its bracketing
    search happens to find, not necessarily the global one over the whole
    domain -- a real risk here specifically, since these energies come
    from independent SIESTA single-points on a fixed, possibly imperfect
    geometric path (no inter-image relaxation coupling at the DFT level),
    which is noisier than a converged, coupled NEB profile and more prone
    to a spline with more than one local hump. Evaluating the spline on a
    dense grid first to find the global maximum, then polishing with a
    bounded local search over just the grid cell it falls in, costs
    nothing (spline evaluation is cheap) and is robust to that.
    """
    if energy_spline is None:
        return None
    x, y, spline = energy_spline
    from scipy.optimize import minimize_scalar
    x_dense = np.linspace(x[0], x[-1], 500)
    y_dense = spline(x_dense)
    i_max = int(np.argmax(y_dense))
    lo = x_dense[max(i_max - 1, 0)]
    hi = x_dense[min(i_max + 1, len(x_dense) - 1)]
    result = minimize_scalar(lambda xv: -spline(xv), bounds=(lo, hi), method='bounded')
    ts_coord = float(result.x)
    ts_energy = float(spline(ts_coord))
    return ts_coord, ts_energy - y[0], y[-1] - y[0]


def read_band_structures(images_root, image_rows_sorted):
    """Reads every image's structure.fdf (modes 1/2/3) into pymatgen
    Structures, in band order, via structure_io.read_fdf/to_pymatgen --
    the same primitives neb.py itself uses to write these same folders.
    Returns a list of (label, Structure_or_None) -- None for any image
    whose structure.fdf is missing/unreadable, reported by the caller as
    a [WARNING], never fatal here (consistency checks/--view-path already
    tolerate a partial band).
    """
    result = []
    for label, _index, _coord in image_rows_sorted:
        path = os.path.join(images_root, label, "structure.fdf")
        try:
            structure = structure_io.to_pymatgen(structure_io.read_fdf(path))
        except (OSError, ValueError):
            structure = None
        result.append((label, structure))
    return result


def print_consistency_checks(f_out, image_rows, band_structures, use_reaction_coord,
                              barrier_fit, setup_report_path, images_root, table_found,
                              is_refined_cycle):
    """Prints '[2] CONSISTENCY CHECKS': re-validates the band as it
    actually exists on disk right now (not just what stb-neb reported at
    prep time), each check ending with an explicit [OK]/[WARNING] line --
    never silent on success, matching neb.py's own reporting model.
    """
    print_section('[2] CONSISTENCY CHECKS', f_out)

    # (a) declared vs. found image count
    found_dirs = sorted(d for d in os.listdir(images_root)
                         if os.path.isdir(os.path.join(images_root, d)) and d.startswith("image_"))
    declared_labels = {label for label, _i, _c in image_rows}
    found_labels = set(found_dirs)
    if declared_labels == found_labels:
        print_dual(f"  [OK] Image count: {len(found_labels)} folder(s) on disk, matching "
                    f"{'the declared image table' if table_found else 'the glob fallback'}.", f_out)
    else:
        missing = declared_labels - found_labels
        extra = found_labels - declared_labels
        detail = []
        if missing:
            detail.append(f"missing: {', '.join(sorted(missing))}")
        if extra:
            detail.append(f"unexpected extra: {', '.join(sorted(extra))}")
        print_dual(color_text(
            f"  [WARNING] Declared/found image folders disagree ({'; '.join(detail)}).",
            'yellow'), f_out)

    # (b) reaction-coordinate ordering
    if use_reaction_coord:
        coords = [r[2] for r in sorted(image_rows, key=lambda r: r[1])]
        bad_pairs = [(i, i + 1) for i in range(len(coords) - 1) if coords[i + 1] <= coords[i]]
        if not bad_pairs:
            print_dual(f"  [OK] Reaction coordinate is strictly increasing along the band "
                        f"(0.0 to {coords[-1]:.4f} Ang).", f_out)
        else:
            pairs_str = ", ".join(f"image_{i:02d}/image_{j:02d}" for i, j in bad_pairs)
            print_dual(color_text(
                f"  [WARNING] Reaction coordinate is not strictly increasing at: {pairs_str} -- "
                "this is why a spline-fitted barrier may be unavailable below.", 'yellow'), f_out)
    else:
        print_dual("  [NOTE] No reaction-coordinate data available (index-only fallback) -- "
                    "ordering check skipped.", f_out)

    # (c) cross-image structural consistency
    valid_structures = [(label, s) for label, s in band_structures if s is not None]
    missing_structs = [label for label, s in band_structures if s is None]
    if missing_structs:
        print_dual(color_text(
            f"  [WARNING] Could not read structure.fdf for {len(missing_structs)} image(s): "
            f"{', '.join(missing_structs)} -- structural checks below only cover the rest.",
            'yellow'), f_out)
    if len(valid_structures) >= 2:
        ref_label, ref_structure = valid_structures[0]
        ref_counts = Counter(str(s.specie) for s in ref_structure)
        mismatches = []
        for label, structure in valid_structures[1:]:
            counts = Counter(str(s.specie) for s in structure)
            if counts != ref_counts:
                mismatches.append(label)
        if not mismatches:
            n_atoms = sum(ref_counts.values())
            formula = " ".join(f"{sym}{n}" for sym, n in sorted(ref_counts.items()))
            print_dual(f"  [OK] All {len(valid_structures)} image(s) share the same composition "
                        f"({formula}, {n_atoms} atoms) -- re-validated from the folders as they "
                        "exist now, not just at prep time.", f_out)
        else:
            print_dual(color_text(
                f"  [WARNING] Composition mismatch vs. {ref_label} in: {', '.join(mismatches)} -- "
                "these folders may have been hand-edited or partially regenerated since stb-neb "
                "wrote them.", 'yellow'), f_out)
    else:
        print_dual("  [NOTE] Fewer than 2 readable structures -- structural consistency check "
                    "skipped.", f_out)

    # (d) path-length cross-check
    if len(valid_structures) >= 2:
        total_length = 0.0
        for (_l0, s0), (_l1, s1) in zip(valid_structures[:-1], valid_structures[1:]):
            total_length += float(np.linalg.norm(s1.cart_coords - s0.cart_coords))
        declared_length = None
        if setup_report_path and os.path.isfile(setup_report_path):
            with open(setup_report_path) as f:
                m = _TOTAL_PATH_LENGTH_RE.search(f.read())
            if m:
                declared_length = float(m.group(1))
        if declared_length is None:
            print_dual("  [NOTE] No 'Total path length' line found in neb_setup.txt (older "
                        "report format, or the file is missing) -- path-length cross-check "
                        "skipped.", f_out)
        elif is_refined_cycle:
            delta = total_length - declared_length
            print_dual(f"  [OK] Path length has evolved by {delta:+.4f} Ang under real-DFT "
                        f"refinement since the original interpolated guess ({declared_length:.4f} "
                        f"Ang -> {total_length:.4f} Ang, informational).", f_out)
        else:
            delta = total_length - declared_length
            tol = 0.05
            if abs(delta) <= tol:
                print_dual(f"  [OK] Recomputed path length ({total_length:.4f} Ang) matches "
                            f"stb-neb's own record ({declared_length:.4f} Ang, from "
                            "neb_setup.txt).", f_out)
            else:
                print_dual(color_text(
                    f"  [WARNING] Recomputed path length ({total_length:.4f} Ang) differs from "
                    f"stb-neb's own record ({declared_length:.4f} Ang) by {delta:+.4f} Ang -- "
                    "the image folders may have moved since stb-neb wrote them.", 'yellow'), f_out)
    else:
        print_dual("  [NOTE] Fewer than 2 readable structures -- path-length cross-check "
                    "skipped.", f_out)

    # (e) spline-fitted TS plausibility
    if barrier_fit is not None and use_reaction_coord:
        ts_coord, _barrier, _dE = barrier_fit
        coords = sorted(r[2] for r in image_rows)
        span = coords[-1] - coords[0]
        frac = (ts_coord - coords[0]) / span if span > 0 else 0.5
        if frac < 0.05 or frac > 0.95:
            print_dual(color_text(
                f"  [WARNING] The spline-fitted transition state sits very close to an endpoint "
                f"({frac * 100:.1f}% along the path) -- the barrier may be underestimated, or "
                "the path may need extending.", 'yellow'), f_out)
        else:
            print_dual(f"  [OK] Spline-fitted transition state sits {frac * 100:.1f}% along the "
                        "path (not degenerately close to an endpoint).", f_out)
    else:
        print_dual("  [NOTE] No spline-fitted transition state available -- plausibility check "
                    "skipped (see BARRIER ANALYSIS below for why).", f_out)


def write_curve_plot(dat_path, rows, use_reaction_coord, smooth):
    """Writes <dat_path> plus a companion .gplot: energy vs. reaction
    coordinate (or plain image index, if reaction_coord wasn't available)
    -- a continuous path, so linespoints, not a scatter (unlike
    adsorb_analysis.py's per-site curve, which is a discrete set of
    candidates). Same <name>.dat + <name>.gplot convention as the rest of
    the suite. `smooth`: True adds a second 'smooth csplines' series (a
    genuinely new convention for this suite, requested to visually match
    the matplotlib PNG's own interpolated curve -- csplines, not
    acsplines, since NEB image energies are exact values, not noisy
    samples) -- caller passes True only when the same >=4-point/no-
    duplicate-x guard fit_energy_spline itself uses holds, so this never
    hits gnuplot's own failure modes for smooth on non-monotonic/
    duplicate x.

    Energies are plotted RELATIVE to the first (lowest-index) point --
    E - E(first point) -- not the raw absolute SIESTA total energy, purely
    for readability (an absolute LCAO total energy is a large negative
    number with all the interesting structure squeezed into its last few
    digits). The report's own [1] IMAGE ENERGIES table still prints the
    real, absolute per-image energy -- only this plot is shifted.
    """
    rows_sorted = sorted(rows, key=lambda r: r.index)
    e0 = rows_sorted[0].energy
    xlabel_col = "ReactionCoord(Ang)" if use_reaction_coord else "ImageIndex"
    with open(dat_path, 'w') as f:
        f.write("# NEB energy profile\n")
        f.write(f"# Energy relative to the first point ({rows_sorted[0].label}, "
                f"E = {e0:.6f} eV)\n")
        f.write(f"# 1:{xlabel_col} 2:E-E0(eV) 3:Label\n")
        for r in rows_sorted:
            x = r.reaction_coord if use_reaction_coord else r.index
            f.write(f"{x:.6f}  {r.energy - e0:.6f}  {r.label}\n")

    base = os.path.splitext(dat_path)[0]
    base_name = os.path.basename(base)
    gplot_path = f"{base}.gplot"
    dat_name = os.path.basename(dat_path)
    xlabel = "Reaction coordinate (Ang)" if use_reaction_coord else "Image index"
    if smooth:
        plot_lines = [
            f'plot "{dat_name}" using 1:2 with points pt 7 ps 1.5 lc rgb "#2255cc" '
            'title "E (images)", \\\n',
            '     "" using 1:2 smooth csplines lw 2 lc rgb "#cc5522" '
            'title "Cubic spline (interpolated)"\n',
        ]
    else:
        plot_lines = [
            '# Not enough points (or duplicate/non-monotonic reaction coordinates) for a\n',
            '# smooth interpolated curve -- plotting the raw path only.\n',
            f'plot "{dat_name}" using 1:2 with linespoints lw 2 pt 7 lc rgb "#2255cc" title "E"\n',
        ]
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-nebAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            'set title "NEB energy profile"\n',
            f'set xlabel "{xlabel}"\n',
            'set ylabel "E - E(first point) (eV)"\n',
            'set grid\n',
            'set key top right\n',
        ] + plot_lines)
    return gplot_path


def build_energy_profile_figure(rows, use_reaction_coord, energy_spline, barrier_fit, title):
    """Builds (does not save/close/show) a matplotlib Figure: raw image
    markers connected by a light line, plus -- when energy_spline is
    available -- a denser cubic-spline curve (the SAME spline already fit
    for the barrier estimate, reused here rather than refit) with the
    fitted transition-state point marked. Returns the Figure; the caller
    decides save/close/show timing (see main()'s [4] ENERGY PROFILE PLOT
    section).

    Energies are plotted RELATIVE to the first (lowest-index) point --
    E - E(first point) -- not the raw absolute SIESTA total energy, purely
    for readability (same reasoning/convention as write_curve_plot's own
    gnuplot .dat; the report's own [1] IMAGE ENERGIES table still prints
    the real, absolute per-image energy).
    """
    rows_sorted = sorted(rows, key=lambda r: r.index)
    xs = [r.reaction_coord if use_reaction_coord else r.index for r in rows_sorted]
    e0 = rows_sorted[0].energy
    ys = [r.energy - e0 for r in rows_sorted]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(xs, ys, '-', color='#b0b0b0', linewidth=1, zorder=1)
    ax.plot(xs, ys, 'o', color='#2255cc', markersize=7, label='E (images)', zorder=2)

    if energy_spline is not None:
        x, y0, spline = energy_spline
        x_dense = np.linspace(x[0], x[-1], 400)
        ax.plot(x_dense, spline(x_dense) - e0, '-', color='#cc5522', linewidth=2,
                label='Cubic spline (interpolated)', zorder=1.5)
        if barrier_fit is not None:
            ts_coord, spline_barrier, _dE = barrier_fit
            ts_energy = y0[0] + spline_barrier - e0
            ax.plot([ts_coord], [ts_energy], '*', color='#22aa44', markersize=16,
                    label=f'TS (fitted): {spline_barrier:.4f} eV', zorder=3)

    ax.set_xlabel("Reaction coordinate (Ang)" if use_reaction_coord else "Image index")
    ax.set_ylabel("E - E(first point) (eV)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig


def view_band_path(band_structures):
    """Converts every (label, Structure) into an ase.Atoms (via
    AseAtomsAdaptor.get_atoms, same conversion neb.py itself performs)
    and opens them as a multi-frame band in ASE's interactive viewer, in
    band order -- step through the reaction path with the frame slider.
    Structures that couldn't be read are silently skipped (already
    reported as a [WARNING] in [2] CONSISTENCY CHECKS).
    """
    atoms_list = [AseAtomsAdaptor.get_atoms(s) for _label, s in band_structures if s is not None]
    if not atoms_list:
        print(color_text("[FAIL] No readable structures to view.", 'red'))
        return
    view_structure_interactive(atoms_list)


def collect_cycle_barrier_history(root_dir, out_filename):
    """For modes 2/3: forward barrier (highest-energy image minus the
    first image) per COMPLETE cycle_NN found under root_dir -- 'complete'
    meaning every image_* in that cycle already has a readable energy.
    Incomplete cycles (SIESTA still queued/running, or a stale partial
    cycle) are silently skipped -- this is a bonus convergence-history
    view on top of the main analysis (which only ever looks at the
    single latest/converged cycle), not a second required data source.
    Returns a sorted list of (cycle_num, forward_barrier_eV).
    """
    history = []
    for entry in sorted(os.listdir(root_dir)):
        m = _CYCLE_DIR_RE.match(entry)
        if not m:
            continue
        cycle_dir = os.path.join(root_dir, entry)
        image_dirs = sorted(d for d in os.listdir(cycle_dir)
                             if os.path.isdir(os.path.join(cycle_dir, d)) and d.startswith("image_"))
        if not image_dirs:
            continue
        energies = []
        for label in image_dirs:
            out_path = os.path.join(cycle_dir, label, out_filename)
            e = siesta_log.get_free_energy(out_path) if os.path.isfile(out_path) else None
            if e is None:
                energies = None
                break
            energies.append(e)
        if energies is None:
            continue
        history.append((int(m.group(1)), max(energies) - energies[0]))
    return history


def write_cycle_convergence_plot(history, out_path):
    """Forward barrier vs. cycle number -- shows the real-DFT refinement
    (modes 2/3) actually settling down (or not) as stb-nebCycle iterates.
    """
    cycles = [c for c, _b in history]
    barriers = [b for _c, b in history]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(cycles, barriers, marker='o', color='tab:blue')
    ax.set_xlabel("Cycle")
    ax.set_ylabel("Forward barrier (eV)")
    ax.set_title("Real-DFT NEB refinement: barrier vs. cycle")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _print_library_warnings_section(f_out, library_warnings, section_label):
    """Prints the shared 'LIBRARY WARNINGS' section, mirroring neb.py's
    own _print_library_warnings_section -- collected noise from
    matplotlib/scipy/ASE instead of interleaving it with the report.
    """
    print_section(section_label, f_out)
    if library_warnings:
        print_dual(color_text(
            "  Messages emitted by external libraries (matplotlib/scipy/ASE) during this run -- "
            "collected here instead of interleaved with the report above; harmless in almost "
            "every case (deprecation-style warnings), but worth a glance.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(f"  {entry}", f_out)
    else:
        print_dual("  No library warnings.", f_out)


def _default_analysis_dir():
    """Smart default for --dir: prefers 'neb_run' (Stage 1's own
    self-contained run folder), but falls back to '.' when the CURRENT
    directory itself already looks like a run folder (has its own
    neb_setup.txt) -- otherwise defaulting to 'neb_run' would try to
    descend into a non-existent 'neb_run/neb_run' and fail with a
    misleading "no image folders found". This exact situation happens
    live: stb-neb's own printed cluster-submission snippet (--mode 2/3)
    does `cd "$run_root"` (i.e. into neb_run/) for its SIESTA loop, so a
    user who runs stb-nebAnalysis right after that loop finishes,
    without cd'ing back out first, is standing INSIDE neb_run/ already.
    """
    if os.path.isfile(SETUP_REPORT_FILE):
        return "."
    return "neb_run"


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes a NEB energy profile from an stb-neb band: reads "
        "each image_NN/'s SIESTA single-point energy and estimates the reaction barrier from "
        "the highest-energy image.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s                              # analyzes ./neb_run by default\n"
               "  %(prog)s --dir . --file calc.out\n"
               "  %(prog)s --apply ts_guess.fdf\n"
               "  %(prog)s --save-gnuplot --view\n"
               "  %(prog)s --view-path\n"
               "  %(prog)s --save-path-xyz\n"
    )

    parser.add_argument("--dir", type=str, default=_default_analysis_dir(),
                         help="Root directory containing every 'image_NN/' (default: 'neb_run', "
                              "the self-contained run folder stb-neb --Stage 1 always writes its "
                              "output into -- auto-detected as '.' instead when the current "
                              "directory already has its own neb_setup.txt, e.g. right after "
                              "stb-neb's own printed cluster-submission snippet left you inside "
                              "neb_run/ already).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-o", "--output", type=str, default="neb_curve.dat",
                         help="Base filename for the gnuplot data file (default: neb_curve.dat), "
                              "written under '<dir>/plot/' -- only relevant with --save-gnuplot.")
    parser.add_argument("--apply", type=str, default=None, metavar="STRUCTURE_FDF",
                         help="Copy the highest-energy image's structure.fdf (the approximate "
                              "transition-state guess) to this path.")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05, same as "
                              "stb-adsorbAnalysis/stb-cohesiveAnalysis) above which an image's "
                              "calc.out is flagged as possibly not relaxed. Advisory only, never "
                              "blocks the result.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                         help="Also save the energy-profile curve as gnuplot .dat + .gplot "
                              "scripts, under '<dir>/plot/'. Off by default.")
    parser.add_argument("--view", action="store_true",
                         help="View the energy-profile plot interactively via matplotlib now, in "
                              "addition to always saving it as a PNG. Off by default. Needs a "
                              "display.")
    parser.add_argument("--view-path", dest="view_path", action="store_true",
                         help="Open every image's structure (in band order) in ASE's interactive "
                              "multi-frame 3D viewer, to step through the reaction path -- after "
                              "everything else has finished. Needs a display. Independent of "
                              "--view (which shows the energy-vs-reaction-coordinate plot, not "
                              "structures).")
    parser.add_argument("--save-path-xyz", dest="save_path_xyz", action="store_true",
                         help="Write the currently-analyzed band (the same images used for "
                              "[1] IMAGE ENERGIES above -- the latest complete cycle for "
                              "--mode 2/3, converged or not) as a single multi-frame XYZ "
                              "trajectory, '<dir>/neb_path_current.xyz'. Off by default. "
                              "Distinct from stb-neb's own 'neb_path.xyz', which is written "
                              "once at prep time from the pre-refinement path and never "
                              "updated afterwards.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nebAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()
    run_start = time.monotonic()
    library_warnings = []

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze a NEB reaction-path study:", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.dir):
        print(color_text(f"[ERROR] '{args.dir}' not found.", 'red'))
        sys.exit(1)

    mode = read_mode(args.dir)

    images_root = args.dir
    cycle_note = None
    is_refined_cycle = False
    if mode in (2, 3):
        cycle_dir, cycle_converged = find_analysis_cycle(args.dir, args.file)
        if cycle_dir is None:
            print(color_text(
                f"[ERROR] No 'cycle_NN' folder found in '{args.dir}'. Did you run stb-neb "
                "--mode 2/3 and at least one stb-nebCycle round?", 'red'))
            sys.exit(1)
        images_root = cycle_dir
        is_refined_cycle = os.path.basename(cycle_dir) != "cycle_00"
        cycle_num = int(_CYCLE_DIR_RE.match(os.path.basename(cycle_dir)).group(1))
        next_dir = os.path.join(args.dir, f"cycle_{cycle_num + 1:02d}")
        pending_note = (
            f" A newer 'cycle_{cycle_num + 1:02d}' folder already exists, waiting for SIESTA "
            "to finish there before the next stb-nebCycle round."
            if not cycle_converged and os.path.isdir(next_dir) else "")
        cycle_note = (
            f"Analyzing '{os.path.basename(cycle_dir)}'"
            + (" -- NEB_CONVERGED." if cycle_converged
               else " -- NOT YET CONVERGED, this is the latest cycle with SIESTA results so "
                    f"far; re-run once stb-nebCycle's loop finishes.{pending_note}")
        )

    image_rows, ml_neb_used = read_image_table(args.dir)
    table_found = image_rows is not None
    if not table_found:
        dirs = sorted(d for d in os.listdir(images_root)
                       if os.path.isdir(os.path.join(images_root, d)) and d.startswith("image_"))
        image_rows = [(d, i, None) for i, d in enumerate(dirs)]

    if not image_rows:
        print(color_text(f"[ERROR] No 'image_*' folders found in '{images_root}'. Did you run "
                          "stb-neb?", 'red'))
        sys.exit(1)

    image_rows_sorted = sorted(image_rows, key=lambda r: r[1])
    band_structures = read_band_structures(images_root, image_rows_sorted)

    report_path = os.path.join(args.dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== NEB ENERGY PROFILE REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Started    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory  : {args.dir}", f_out)
    print_dual(f"Output file: {args.file}", f_out)
    print_dual(f"Mode       : {mode} ({MODE_DESCRIPTIONS[mode]})", f_out)
    if band_structures and band_structures[0][1] is not None:
        counts = Counter(str(s.specie) for s in band_structures[0][1])
        formula = " ".join(f"{sym}{n}" for sym, n in sorted(counts.items()))
        print_dual(f"Composition: {formula} ({sum(counts.values())} atoms)", f_out)
    print_dual(f"Force tolerance: {args.force_tolerance} eV/Ang (residual-force advisory "
                "threshold)", f_out)
    if cycle_note:
        print_dual(color_text(f"[NOTE] {cycle_note}", 'cyan'), f_out)
    if not table_found:
        print_dual(color_text(
            "[NOTE] No 'neb_setup.txt' image table found -- falling back to a sorted glob of "
            "'image_*' folders (index = sort position, no reaction-coordinate data available; "
            "the x-axis below uses plain image index instead).", 'yellow'), f_out)

    print_section('[1] IMAGE ENERGIES', f_out)
    rows = []
    n_skipped = 0
    scf_warn_labels = []
    force_warn_labels = []
    header = f"{'Image':<14}{'Index':<7}{'ReactionCoord':<16}{'E(eV)':<16}{'SCF':<6}{'MaxF(eV/A)':<12}"
    print_dual(header, f_out)
    print_dual("-" * len(header), f_out)
    for label, index, reaction_coord in image_rows:
        image_dir = os.path.join(images_root, label)
        out_path = os.path.join(image_dir, args.file)
        if not os.path.exists(out_path):
            n_skipped += 1
            print_dual(f"{label:<14}{color_text('SKIP', 'yellow')} (missing {args.file})", f_out)
            continue
        energy = siesta_log.get_free_energy(out_path)
        if energy is None:
            n_skipped += 1
            print_dual(f"{label:<14}{color_text('SKIP', 'yellow')} (could not parse energy)", f_out)
            continue
        scf_ok, max_force = check_scf_and_force(out_path)
        if not scf_ok:
            scf_warn_labels.append(label)
        if max_force is not None and max_force > args.force_tolerance:
            force_warn_labels.append(label)

        rows.append(ImageRow(label, index, reaction_coord, energy, scf_ok, max_force))
        coord_str = f"{reaction_coord:<16.4f}" if reaction_coord is not None else f"{'--':<16}"
        scf_str = color_text("WARN", 'yellow') if not scf_ok else "OK"
        force_str = f"{max_force:.4f}" if max_force is not None else "--"
        print_dual(f"{label:<14}{index:<7}{coord_str}{energy:<16.6f}{scf_str:<6}{force_str:<12}", f_out)
    print_dual("-" * len(header), f_out)
    if scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(scf_warn_labels)} image(s) never confirmed SCF convergence -- "
            f"their energy may be unreliable: {', '.join(scf_warn_labels)}.", 'yellow'), f_out)
    if force_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(force_warn_labels)} image(s) have residual force above "
            f"--force-tolerance ({args.force_tolerance} eV/Ang), possibly not single-point-"
            f"converged: {', '.join(force_warn_labels)}.", 'yellow'), f_out)

    if not rows:
        print_dual(color_text("\n[ERROR] No valid image results found.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    use_reaction_coord = all(r.reaction_coord is not None for r in rows)
    with capture_library_noise(library_warnings, "scipy cubic spline fit"):
        energy_spline = fit_energy_spline(rows) if use_reaction_coord else None
        barrier_fit = find_spline_barrier(energy_spline)

    setup_report_path = os.path.join(args.dir, SETUP_REPORT_FILE)
    print_consistency_checks(f_out, image_rows, band_structures, use_reaction_coord,
                              barrier_fit, setup_report_path, images_root, table_found,
                              is_refined_cycle)

    print_section('[3] BARRIER ANALYSIS', f_out)
    max_index = max(idx for _, idx, _ in image_rows)
    initial_row = next((r for r in rows if r.index == 0), None)
    final_row = next((r for r in rows if r.index == max_index), None)
    if initial_row is None or final_row is None:
        print_dual(color_text(
            "[ERROR] Could not read the initial and/or final endpoint image's energy -- "
            "cannot compute a barrier.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    max_row = max(rows, key=lambda r: r.energy)
    forward_barrier = max_row.energy - initial_row.energy
    backward_barrier = max_row.energy - final_row.energy
    delta_e = final_row.energy - initial_row.energy

    print_dual(f"Highest-energy image (approx. TS) : {max_row.label}  "
                f"(E = {max_row.energy:.6f} eV)", f_out)
    print_dual(f"Forward barrier  (TS - initial)   : {forward_barrier:.6f} eV", f_out)
    print_dual(f"Backward barrier (TS - final)     : {backward_barrier:.6f} eV", f_out)
    rxn_verdict = "exothermic (favorable)" if delta_e < 0 else "endothermic (unfavorable)"
    print_dual(f"Reaction energy  (final - initial) : {delta_e:.6f} eV, {rxn_verdict}", f_out)

    if barrier_fit:
        ts_coord, spline_barrier, spline_dE = barrier_fit
        print_dual(f"Spline-fitted barrier (smoothed)  : {spline_barrier:.6f} eV at reaction "
                    f"coordinate {ts_coord:.4f} Ang (cubic spline through the {len(rows)} "
                    "energies)", f_out)
    elif use_reaction_coord and len(rows) < 4:
        print_dual("[NOTE] Not enough images (need >= 4) for a spline-fitted barrier estimate.", f_out)

    if mode in (2, 3):
        caveat = ("This cycle has gone through real-DFT NEB refinement (stb-nebCycle) -- a "
                   "genuine coupled-spring optimization using SIESTA forces, not just "
                   "independent single points on a fixed guess -- but see the convergence "
                   "status noted in [0] before trusting it as fully converged.")
    elif ml_neb_used:
        caveat = ("Stage 1 already converged a real climbing-image band on MACE-MP-0 first, "
                   "so this is a stronger (ML-relaxed-path) estimate -- but still not an "
                   "independent DFT confirmation of a converged saddle point.")
    else:
        caveat = ("Approximate, interpolated-path estimate: no real reaction coordinate was "
                   "optimized between these DFT single points (no inter-image spring "
                   "coupling at the DFT level) -- the true saddle point may lie off this "
                   "fixed path. Re-run stb-neb with --mode 1/2/3 for a physically relaxed band "
                   "shape before the DFT single points.")
    print_dual(color_text(f"[NOTE] {caveat}", 'yellow'), f_out)

    if mode in (2, 3):
        print_section('[3b] BARRIER VS. CYCLE (real-DFT refinement history)', f_out)
        cycle_history = collect_cycle_barrier_history(args.dir, args.file)
        if len(cycle_history) >= 2:
            convergence_plot_path = os.path.join(args.dir, "neb_cycle_convergence.png")
            with capture_library_noise(library_warnings, "matplotlib cycle-convergence plot"):
                write_cycle_convergence_plot(cycle_history, convergence_plot_path)
            print_dual(f"{color_text('[Saved]', 'cyan')} {convergence_plot_path} "
                        f"({len(cycle_history)} complete cycle(s))", f_out)
            for cyc, barrier in cycle_history:
                print_dual(f"  cycle_{cyc:02d}: forward barrier = {barrier:.6f} eV", f_out)
        else:
            print_dual("  Not enough complete cycles yet for a convergence curve (need >= 2 "
                        "cycles with every image's calc.out readable).", f_out)

    print_section('[4] ENERGY PROFILE PLOT', f_out)
    png_path = os.path.join(args.dir, "neb_energy_profile.png")
    with capture_library_noise(library_warnings, "matplotlib energy-profile plot"):
        fig = build_energy_profile_figure(
            rows, use_reaction_coord, energy_spline, barrier_fit, title="NEB energy profile")
        fig.savefig(png_path, dpi=150)
    print_dual(f"{color_text('[Saved]', 'cyan')} {png_path}", f_out)
    if args.save_gnuplot:
        plot_dir = os.path.join(args.dir, "plot")
        os.makedirs(plot_dir, exist_ok=True)
        dat_path = os.path.join(plot_dir, os.path.basename(args.output))
        gplot_path = write_curve_plot(dat_path, rows, use_reaction_coord,
                                       smooth=energy_spline is not None)
        print_dual(f"{color_text('[Saved]', 'cyan')} {dat_path}, {gplot_path} "
                    f"(cd {plot_dir} && gnuplot {os.path.basename(gplot_path)})", f_out)
    else:
        print_dual("Gnuplot data+script : not written (off by default -- pass --save-gnuplot to "
                    "write it under a 'plot/' subfolder).", f_out)
    if not args.view:
        plt.close(fig)

    elapsed = time.monotonic() - run_start
    print_section('[5] SUMMARY', f_out)
    print_dual(f"Images analyzed : {len(rows)} (skipped: {n_skipped})", f_out)
    print_dual(f"Total elapsed   : {elapsed:.1f}s", f_out)
    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report -> {report_path}", f_out)

    if args.apply:
        print_section('[6] APPLY', f_out)
        src = os.path.join(images_root, max_row.label, "structure.fdf")
        try:
            shutil.copy(src, args.apply)
        except OSError as e:
            print_dual(color_text(f"[ERROR] Could not copy '{src}' to '{args.apply}': {e}", 'red'), f_out)
        else:
            print_dual(f"{color_text('[Applied]', 'green')} {max_row.label} -> {args.apply}", f_out)

    if args.save_path_xyz:
        print_section('[6b] PATH EXPORT', f_out)
        valid_band = [(label, s) for label, s in band_structures if s is not None]
        if len(valid_band) < 2:
            print_dual(color_text(
                "[WARNING] Fewer than 2 structures available (missing structure.fdf files) -- "
                "not enough to write a path trajectory.", 'yellow'), f_out)
        else:
            if mode in (2, 3):
                status = ("NEB_CONVERGED" if cycle_converged
                           else f"NOT YET CONVERGED, latest complete cycle "
                                f"({os.path.basename(images_root)})")
            elif mode == 1:
                status = "single-point on Stage 1's MACE-MP-0 path (no DFT-refinement cycles)"
            else:
                status = "unrefined interpolated path"
            xyz_path = os.path.join(args.dir, "neb_path_current.xyz")
            write_path_trajectory([s for _l, s in valid_band], xyz_path)
            print_dual(f"{color_text('[Saved]', 'cyan')} {xyz_path} ({len(valid_band)} frame(s), "
                        f"{status}).", f_out)

    _print_library_warnings_section(f_out, library_warnings, '[7] LIBRARY WARNINGS')

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} mode {mode} analysis complete "
          f"({len(rows)} image(s), forward barrier {forward_barrier:.4f} eV).")
    if report_path:
        print(f"Full report: {report_path}")

    if args.view:
        plt.show()
    if args.view_path:
        print(f"\n{color_text('--view-path:', 'cyan')} opening {len(band_structures)} frame(s) "
              "in ASE's interactive viewer (use the frame slider to step through the reaction "
              "path):")
        for i, (label, _s) in enumerate(band_structures):
            print(f"  {i} = {label}")
        view_band_path(band_structures)


if __name__ == "__main__":
    main()
