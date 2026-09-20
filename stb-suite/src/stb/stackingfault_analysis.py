#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import json
import shutil
import argparse
from datetime import datetime
import numpy as np
from ase.io import write as ase_write
from stb.core import siesta_log, structure_io
from stb.core.siesta_log import check_scf_and_force
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.grid_export import write_gnuplot_script, check_planar_orthogonality
from stb.core.ase_view import view_structure_interactive

# Same self-contained run-folder convention as stb-neb's neb_run/ -- Stage 1
# (stb-stackingfault) always writes into <output-dir>/sf_run/, so this is
# also RUN_SUBDIR's default value for --dir (see _default_analysis_dir()).
RUN_SUBDIR = "sf_run"
MANIFEST_FILE = "sf_manifest.json"
REPORT_FILE = "stackingfault_report.txt"
_SHIFT_DIR_RE = re.compile(r'^shift_(\d+)_(\d+)$')


def read_grid_manifest(root_dir):
    """Reads Stage 1's always-on sf_manifest.json -- the primary,
    machine-readable source of the grid layout (label, i, j, shift_x,
    shift_y per point, plus grid_nx/grid_ny/scan/scan_points/mode/
    n_layer1_atoms), independent of whether --save-report was used for
    the narrative report. Returns (rows, grid_nx, grid_ny, scan,
    scan_points, mode, n_layer1_atoms), or None if the manifest is
    missing/unreadable, in which case the caller falls back to a sorted
    glob of 'shift_II_JJ' folders (see fallback_grid_rows()). Mirrors
    adsorb_analysis.py::read_site_table / neb_analysis.py::read_image_table
    -- the grid directory is recomputed by the caller as
    os.path.join(root_dir, label), not trusted from the manifest.
    """
    manifest_path = os.path.join(root_dir, MANIFEST_FILE)
    if not os.path.isfile(manifest_path):
        return None
    try:
        with open(manifest_path) as f:
            manifest = json.load(f)
        rows = [(r["label"], r["i"], r["j"], r["shift_x"], r["shift_y"])
                for r in manifest["rows"]]
        if not rows:
            return None
        # "mode"/"n_layer1_atoms" are absent in manifests written before
        # stb-stackingfault's --mode {1,2,3} existed -- treat those as mode
        # 3 (the previous, only behavior: fixed gap, plain single-point).
        # "scan"/"scan_points" are absent in manifests written before
        # --scan {surface,x,y,xy} existed -- treat those as 'surface' (the
        # previous, only behavior).
        mode = manifest.get("mode", 3)
        n_layer1_atoms = manifest.get("n_layer1_atoms")
        scan = manifest.get("scan", "surface")
        scan_points = manifest.get("scan_points")
        return rows, manifest["grid_nx"], manifest["grid_ny"], scan, scan_points, mode, n_layer1_atoms
    except (OSError, ValueError, KeyError):
        return None


def fallback_grid_rows(root_dir):
    """Recovers the grid layout from a sorted glob of 'positions/shift_II_JJ/'
    folders when sf_manifest.json is missing -- unlike a generic "no
    metadata at all" fallback, the folder-naming convention here IS the
    (i, j) grid index, and stb-stackingfault always builds shifts via
    linspace(0, 1, grid_nx/grid_ny, endpoint=False) independently per
    axis, so shift_x/shift_y are exactly recoverable as i/grid_nx,
    j/grid_ny once grid_nx/grid_ny are inferred from the max index
    actually present along each axis (asymmetric grids recoverable too,
    same as a square one) -- no data is lost, unlike stb-neb's index-only
    fallback (which can't recover reaction_coord at all). Returns (rows,
    grid_nx, grid_ny), or None if no 'shift_II_JJ' folder is found under
    '<root_dir>/positions/'.

    Does NOT recover --scan (x/y/xy/surface) -- the caller always treats a
    manifest-less run as 'surface'. A genuine 1D scan missing its manifest
    would still reconstruct correct shift_x/shift_y here (its folder names
    already encode a degenerate grid_ny=1/grid_nx=1/diagonal shape that
    this same math handles fine), just attempted as a 'surface' plot
    instead of the right 1D line -- acceptable given the manifest is
    always-on in Stage 1 and only a manually-deleted manifest hits this.
    """
    positions_dir = os.path.join(root_dir, "positions")
    if not os.path.isdir(positions_dir):
        return None
    dirs = [d for d in os.listdir(positions_dir) if os.path.isdir(os.path.join(positions_dir, d))]
    parsed = []
    for d in dirs:
        m = _SHIFT_DIR_RE.match(d)
        if m:
            parsed.append((d, int(m.group(1)), int(m.group(2))))
    if not parsed:
        return None
    grid_nx = max(i for _, i, j in parsed) + 1
    grid_ny = max(j for _, i, j in parsed) + 1
    parsed.sort(key=lambda r: (r[1], r[2]))
    rows = [(d, i, j, i / grid_nx, j / grid_ny) for d, i, j in parsed]
    return rows, grid_nx, grid_ny


class GridRow:
    """One analyzed grid point: label, (i, j) grid index, fractional
    shift, and its computed energy/quality diagnostics. Plain
    attribute-holder, same style as adsorb_analysis.py's SiteRow /
    neb_analysis.py's ImageRow. `final_gap` is None unless mode 1
    (SIESTA-relaxed z) and a relaxed structure was actually found.
    """
    def __init__(self, label, i, j, shift_x, shift_y, energy, scf_ok, max_force, final_gap=None):
        self.label = label
        self.i = i
        self.j = j
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.energy = energy
        self.scf_ok = scf_ok
        self.max_force = max_force
        self.final_gap = final_gap


def effective_gap(pmg_structure, n_layer1):
    """Cartesian z-gap (Ang) between layer 1's topmost atom and layer 2's
    bottommost atom -- same definition/formula as stackingfault.py's own
    effective_gap (duplicated here rather than cross-imported, same
    "duplicated per self-contained workflow stage" convention documented
    for read_site_theory_flags elsewhere in the suite). Used to report the
    REAL gap a --mode 1 (SIESTA) restricted relaxation actually reached.
    """
    n_total = len(pmg_structure)
    if n_layer1 is None or n_layer1 <= 0 or n_layer1 >= n_total:
        return None
    z_max_l1 = max(pmg_structure[i].coords[2] for i in range(n_layer1))
    z_min_l2 = min(pmg_structure[i].coords[2] for i in range(n_layer1, n_total))
    return float(z_min_l2 - z_max_l1)


def cartesian_shift(lattice, shift_x, shift_y):
    """Converts a fractional (shift_x, shift_y) into the actual Cartesian
    in-plane displacement (Ang), using the structure's own a/b lattice
    vectors -- needed for physically meaningful axis units in the
    gamma-surface plot (a plain fractional 0-1 axis would be mislabeled by
    write_gnuplot_script's hardcoded "(Angstrom)" slice-mode labels).
    """
    a_vec = np.array(lattice[0][:2])
    b_vec = np.array(lattice[1][:2])
    return shift_x * a_vec + shift_y * b_vec


def cell_area(lattice):
    """In-plane (a x b) cross-product magnitude, Ang^2 -- the SAME a_vec/
    b_vec (lattice[0][:2]/lattice[1][:2]) vectors cartesian_shift() uses.
    Needed to convert a relative stacking-fault energy from eV into the
    physically standard, supercell-size-independent meV/Ang^2 unit.
    """
    a_vec = np.array(lattice[0][:2])
    b_vec = np.array(lattice[1][:2])
    return float(abs(a_vec[0] * b_vec[1] - a_vec[1] * b_vec[0]))


def to_mev(e_ev):
    """eV -> meV. The base unit-of-account for every relative energy
    quantity in this report (gamma-surface dE, corrugation) -- never
    applied to an absolute per-point SCF energy, which stays in eV. Area
    normalization (mev_per_area) is a distinct, final step, applied only
    when a quantity is actually written/printed as an areal energy
    density -- never baked into this conversion itself.
    """
    return e_ev * 1000.0


def mev_per_area(e_mev, area_ang2):
    """meV -> meV/Ang^2. Applied ONLY at the point a gamma-surface/
    corrugation value is actually written to a file or printed -- see
    to_mev's own docstring.
    """
    return e_mev / area_ang2


def write_surface_data(dat_path, rows, grid_nx, grid_ny, lattice, area, e_min):
    """Writes the gamma-surface data file in gnuplot PM3D block format
    (blank line after each row of constant i) with columns
    dx(Ang) dy(Ang) 0.0(unused) dE(meV/Ang^2) -- matches
    core/grid_export.py::write_data_file's 'slice'-mode convention exactly
    (4 columns, "1:2:4" used by write_gnuplot_script's axis_idx=2 branch),
    so that function can be reused verbatim for the heatmap script instead
    of writing a new one. grid_nx/grid_ny may differ (an asymmetric grid).

    `area`/`e_min` are computed once in main() and passed in, so the
    plotted dE=0 point is numerically identical to what [2]'s report
    calls "equilibrium".
    """
    by_ij = {(r.i, r.j): r for r in rows}
    with open(dat_path, 'w') as f:
        f.write("# Stacking-fault (gamma-surface) energy landscape\n")
        f.write("# 1:dx(Ang) 2:dy(Ang) 3:(unused) 4:dE(meV/Ang^2)\n")
        for i in range(grid_nx):
            for j in range(grid_ny):
                row = by_ij.get((i, j))
                if row is None:
                    continue
                dx, dy = cartesian_shift(lattice, row.shift_x, row.shift_y)
                f.write(f"{dx:.6f} {dy:.6f} 0.0 {mev_per_area(to_mev(row.energy - e_min), area):.6f}\n")
            f.write("\n")


def write_scan_line_data(dat_path, rows, lattice, area, e_min):
    """Writes a --scan x/y/xy line-plot data file: 2 plain columns,
    distance traveled along the scan path (Ang, from the origin) and
    dE (meV/Ang^2) -- the 1D analog of write_surface_data's PM3D block,
    read by write_gnuplot_script's mode='profile' template ('plot ...
    using 1:2'). Distance is the Cartesian norm of cartesian_shift(shift_x,
    shift_y) -- works uniformly for all 3 scan shapes (shift_y is always 0
    for 'x', shift_x is always 0 for 'y', so the norm collapses to the
    single swept component in both cases; for 'xy' it's the true diagonal
    distance). Rows are sorted by this distance (not just manifest order)
    so the line always plots left-to-right even if the caller's row order
    doesn't already match it. See write_surface_data for area/e_min.
    """
    points = []
    for r in rows:
        dx, dy = cartesian_shift(lattice, r.shift_x, r.shift_y)
        distance = float(np.hypot(dx, dy))
        points.append((distance, mev_per_area(to_mev(r.energy - e_min), area)))
    points.sort(key=lambda p: p[0])
    with open(dat_path, 'w') as f:
        f.write("# Stacking-fault energy along a 1D scan path\n")
        f.write("# 1:distance(Ang) 2:dE(meV/Ang^2)\n")
        for distance, e_rel in points:
            f.write(f"{distance:.6f} {e_rel:.6f}\n")


def view_surface_plot(rows, grid_nx, grid_ny, area, e_min):
    """Interactive matplotlib heatmap preview of the gamma-surface -- the
    on-screen counterpart to write_surface_data's saved gnuplot .dat/
    .gplot pair, same meV/Ang^2 unit. Fractional shift_x/shift_y axes
    (like stackingfault.py's own Stage 1 write_ml_preview_plot), NOT the
    saved file's Cartesian dx/dy -- this is a quick on-screen check, never
    written to disk (see CLAUDE.md: WORKFLOW_TOOLS writes gnuplot pairs,
    matplotlib here is preview-only). Degrades gracefully (a NaN gap) on a
    partial grid, unlike gnuplot's pm3d map -- still usable even when
    --save-gnuplot's plot is skipped for an incomplete grid. matplotlib is
    imported lazily, only when --view was actually passed.
    """
    import matplotlib.pyplot as plt
    dE = np.full((grid_nx, grid_ny), np.nan)
    for r in rows:
        dE[r.i, r.j] = mev_per_area(to_mev(r.energy - e_min), area)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(dE.T, origin='lower', extent=(0, 1, 0, 1), cmap='YlOrRd', aspect='equal')
    ax.set_xlabel("shift_x (fractional)")
    ax.set_ylabel("shift_y (fractional)")
    ax.set_title("Stacking-fault energy landscape")
    fig.colorbar(im, ax=ax, label="dE (meV/Ang^2)")
    fig.tight_layout()
    plt.show()


def view_scan_line_plot(rows, lattice, area, e_min, scan):
    """Interactive matplotlib line preview for --scan x/y/xy -- the
    on-screen counterpart to write_scan_line_data, same Cartesian-distance
    x-axis (Ang) as the saved file (Stage 2 always uses Cartesian Ang,
    unlike Stage 1's own fractional-axis preview). Never written to disk,
    matplotlib imported lazily only when --view was actually passed.
    """
    import matplotlib.pyplot as plt
    points = []
    for r in rows:
        dx, dy = cartesian_shift(lattice, r.shift_x, r.shift_y)
        points.append((float(np.hypot(dx, dy)), mev_per_area(to_mev(r.energy - e_min), area)))
    points.sort(key=lambda p: p[0])
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    scan_desc = {"x": "shift_x", "y": "shift_y", "xy": "the shift_x=shift_y diagonal"}
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, ys, marker='o', color='#cc5522')
    ax.set_xlabel("Distance along scan path (Angstrom)")
    ax.set_ylabel("dE (meV/Ang^2)")
    ax.set_title(f"Stacking-fault energy along {scan_desc[scan]}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()


def build_animation_frames(dir_path, rows, lattice, scan, area, e_min):
    """Builds one ase.Atoms per successfully-analyzed grid point
    (structure_io.read_relaxed_or_input + to_pymatgen +
    AseAtomsAdaptor.get_atoms -- same conversion neb.py's own
    write_path_trajectory/view_band_path use), for the always-on
    'stackingfault_animation.xyz' export and the optional --view-animation
    preview -- built once, reused for both.

    Frame order matches the corresponding saved plot's own iteration
    order: row-major (i then j) for scan == 'surface' (same as
    write_surface_data), sorted by Cartesian distance along the path
    otherwise (same as write_scan_line_data) -- so the animation always
    steps through the grid in the same visual order as the saved/
    previewed plot.

    Each frame's .info carries: label, shift_x, shift_y, energy_eV (raw,
    diagnostic, same convention as [1]'s own E(eV) column), dE_meV_per_A2
    (relative to the grid's e_min, same value as [1]'s dE(meV/A2) column).

    A point whose structure can't be re-read is skipped (best-effort,
    never blocks the numeric analysis -- same policy as final_gap) and
    its label collected in the returned `unreadable` list. Returns
    (frames: list[ase.Atoms], unreadable: list[str]).
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    if scan == "surface":
        ordered = sorted(rows, key=lambda r: (r.i, r.j))
    else:
        ordered = sorted(rows, key=lambda r: float(np.hypot(
            *cartesian_shift(lattice, r.shift_x, r.shift_y))))
    frames, unreadable = [], []
    for r in ordered:
        grid_dir = os.path.join(dir_path, "positions", r.label)
        try:
            structure, _used_relaxed = structure_io.read_relaxed_or_input(grid_dir)
        except (OSError, ValueError):
            unreadable.append(r.label)
            continue
        atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure))
        atoms.info.update({
            "label": r.label, "shift_x": r.shift_x, "shift_y": r.shift_y,
            "energy_eV": r.energy, "dE_meV_per_A2": mev_per_area(to_mev(r.energy - e_min), area),
        })
        frames.append(atoms)
    return frames, unreadable


def _default_analysis_dir():
    """Smart default for --dir: prefers 'sf_run' (Stage 1's own
    self-contained run folder), but falls back to '.' when the CURRENT
    directory itself already looks like a run folder (has its own
    sf_manifest.json) -- otherwise defaulting to 'sf_run' would try to
    descend into a non-existent 'sf_run/sf_run' and fail with a misleading
    "no shift_II_JJ folders found". Same convention as stb-nebAnalysis's
    own _default_analysis_dir().
    """
    if os.path.isfile(MANIFEST_FILE):
        return "."
    return RUN_SUBDIR


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes a stacking-fault (gamma-surface) energy landscape "
        "from an stb-stackingfault grid: reads each positions/shift_II_JJ/'s SIESTA energy and "
        "reports the equilibrium stacking, the highest-energy registry, and the corrugation "
        "energy.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --dir . --file calc.out\n"
               "  %(prog)s --dir . --apply equilibrium.fdf\n"
               "  %(prog)s --save-gnuplot --view\n"
               "  %(prog)s --view-animation\n"
    )

    parser.add_argument("--dir", type=str, default=_default_analysis_dir(),
                         help="Root directory containing 'positions/shift_II_JJ/' (default: "
                              "'sf_run', the self-contained run folder stb-stackingfault Stage 1 "
                              "always writes its output into -- auto-detected as '.' instead when "
                              "the current directory already has its own sf_manifest.json).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-o", "--output", type=str, default="stackingfault_surface.dat",
                         help="Base filename for the gnuplot data file (default: "
                              "stackingfault_surface.dat), written under '<dir>/plot/' -- only "
                              "relevant with --save-gnuplot.")
    parser.add_argument("--apply", type=str, default=None, metavar="STRUCTURE_FDF",
                         help="Copy the equilibrium (lowest-energy) grid point's structure.fdf to "
                              "this path.")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05, same as "
                              "stb-adsorbAnalysis/stb-nebAnalysis) above which a grid point's "
                              "calc.out is flagged as possibly not single-point-converged. "
                              "Advisory only, never blocks the result.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                         help="Also save the gamma-surface (or line-scan) data as gnuplot .dat "
                              "+ .gplot scripts, under '<dir>/plot/'. Off by default.")
    parser.add_argument("--view", action="store_true",
                         help="View the gamma-surface heatmap (or 1D line, for --scan x/y/xy) "
                              "interactively via matplotlib now. Off by default. Needs a "
                              "display. Independent of --save-gnuplot -- matplotlib here is "
                              "only an on-screen preview, never written to disk.")
    parser.add_argument("--view-animation", dest="view_animation", action="store_true",
                         help="Open every analyzed grid point's structure (same order as the "
                              "always-saved animation XYZ) in ASE's interactive multi-frame 3D "
                              "viewer. Off by default. Needs a display. Independent of --view "
                              "(which shows the energy plot, not structures).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-stackingfaultAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze a 2D stacking-fault study:", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.dir):
        print(color_text(f"[ERROR] '{args.dir}' not found.", 'red'))
        sys.exit(1)

    manifest_result = read_grid_manifest(args.dir)
    manifest_found = manifest_result is not None
    if manifest_found:
        grid_rows, grid_nx, grid_ny, scan, scan_points, mode, n_layer1_atoms = manifest_result
    else:
        fallback_result = fallback_grid_rows(args.dir)
        if fallback_result is None:
            print(color_text(f"[ERROR] No 'shift_II_JJ' folders found under '{args.dir}/positions/'. "
                              "Did you run stb-stackingfault?", 'red'))
            sys.exit(1)
        grid_rows, grid_nx, grid_ny = fallback_result
        scan, scan_points = "surface", None
        mode, n_layer1_atoms = 3, None

    report_path = os.path.join(args.dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== STACKING-FAULT (GAMMA-SURFACE) REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory  : {args.dir}", f_out)
    print_dual(f"Output file: {args.file}", f_out)
    if scan == "surface":
        print_dual(f"Grid       : {grid_nx} x {grid_ny} (surface)"
                    + (" (asymmetric)" if grid_nx != grid_ny else ""), f_out)
    else:
        scan_desc = {"x": "shift_x only", "y": "shift_y only", "xy": "diagonal shift_x=shift_y"}
        print_dual(f"Grid       : {scan_points} points, 1D ({scan_desc[scan]})", f_out)
    mode_labels = {
        1: "SIESTA relaxed z for real (restricted CG, x/y frozen)",
        2: "MACE-MP-0 relaxed z, then SIESTA single-point",
        3: "fixed gap, plain SIESTA single-point",
    }
    print_dual(f"Mode       : {mode} ({mode_labels.get(mode, 'unknown')})", f_out)
    if not manifest_found:
        print_dual(color_text(
            f"[NOTE] No '{MANIFEST_FILE}' found -- recovered the grid layout from "
            "'shift_II_JJ' folder names instead (shift_x/shift_y reconstructed as "
            "i/grid_nx, j/grid_ny).", 'yellow'), f_out)

    print_section('[1] GRID ENERGIES', f_out)
    rows = []
    entries = []  # ("skip", message) | ("row", GridRow) -- printed once
                  # area/e_min are known (need the full scan first)
    n_skipped = 0
    scf_warn_labels = []
    force_warn_labels = []
    lattice = None
    show_final_gap = mode == 1 and n_layer1_atoms is not None
    for label, i, j, shift_x, shift_y in grid_rows:
        grid_dir = os.path.join(args.dir, "positions", label)
        out_path = os.path.join(grid_dir, args.file)
        if not os.path.exists(out_path):
            n_skipped += 1
            entries.append(("skip", f"{label:<16}{color_text('SKIP', 'yellow')} (missing {args.file})"))
            continue
        energy = siesta_log.get_free_energy(out_path)
        if energy is None:
            n_skipped += 1
            entries.append(("skip", f"{label:<16}{color_text('SKIP', 'yellow')} (could not parse energy)"))
            continue
        if lattice is None:
            struct_path = os.path.join(grid_dir, "structure.fdf")
            if os.path.isfile(struct_path):
                lattice = structure_io.read_fdf(struct_path).lattice
        scf_ok, max_force = check_scf_and_force(out_path)
        if not scf_ok:
            scf_warn_labels.append(label)
        if max_force is not None and max_force > args.force_tolerance:
            force_warn_labels.append(label)

        final_gap = None
        if show_final_gap:
            try:
                relaxed_structure, _used_relaxed = structure_io.read_relaxed_or_input(grid_dir)
                final_gap = effective_gap(structure_io.to_pymatgen(relaxed_structure), n_layer1_atoms)
            except (OSError, ValueError):
                pass

        row = GridRow(label, i, j, shift_x, shift_y, energy, scf_ok, max_force, final_gap)
        rows.append(row)
        entries.append(("row", row))

    # The dE(meV/A2) column needs the whole grid's own e_min/area -- a
    # whole-grid quantity only known once every point has been scanned --
    # so [1]'s table is now printed here, after the loop, instead of
    # interleaved with it. e_min/area are computed once and reused
    # everywhere downstream ([1]'s column, [2]'s report,
    # write_surface_data/write_scan_line_data, both --view preview
    # functions, and build_animation_frames) for one consistent number.
    area = cell_area(lattice) if lattice is not None else None
    if rows and area is not None:
        min_row = min(rows, key=lambda r: r.energy)
        max_row = max(rows, key=lambda r: r.energy)
        e_min = min_row.energy
        e_max = max_row.energy
        corrugation = e_max - e_min
        corrugation_mev = mev_per_area(to_mev(corrugation), area)
    else:
        min_row = max_row = None
        e_min = e_max = corrugation = corrugation_mev = None

    header = f"{'Point':<16}{'Shift(x,y)':<20}{'E(eV)':<16}{'dE(meV/A2)':<14}{'SCF':<6}{'MaxF(eV/A)':<12}"
    if show_final_gap:
        header += f"{'GapFinal(A)':<12}"
    print_dual(header, f_out)
    print_dual("-" * len(header), f_out)
    for entry in entries:
        if entry[0] == "skip":
            print_dual(entry[1], f_out)
            continue
        _, row = entry
        dE_str = (f"{mev_per_area(to_mev(row.energy - e_min), area):<14.3f}" if e_min is not None
                   else f"{'--':<14}")
        shift_str = f"({row.shift_x:.3f},{row.shift_y:.3f})"
        scf_str = color_text("WARN", 'yellow') if not row.scf_ok else "OK"
        force_str = f"{row.max_force:.4f}" if row.max_force is not None else "--"
        row_line = (f"{row.label:<16}{shift_str:<20}{row.energy:<16.6f}{dE_str}"
                    f"{scf_str:<6}{force_str:<12}")
        if show_final_gap:
            row_line += f"{row.final_gap:<12.4f}" if row.final_gap is not None else f"{'--':<12}"
        print_dual(row_line, f_out)
    print_dual("-" * len(header), f_out)
    if scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(scf_warn_labels)} grid point(s) never confirmed SCF convergence "
            f"-- their energy may be unreliable: {', '.join(scf_warn_labels)}.", 'yellow'), f_out)
    if force_warn_labels:
        force_context = ("possibly not converged (mode 1's restricted relaxation may not have "
                          "reached --relax-z-steps' target)" if mode == 1 else
                          "the rigid/pre-relaxed geometry at this point may not be physically "
                          "reasonable (not itself a convergence failure)")
        print_dual(color_text(
            f"[WARNING] {len(force_warn_labels)} grid point(s) have residual force above "
            f"--force-tolerance ({args.force_tolerance} eV/Ang) -- {force_context}: "
            f"{', '.join(force_warn_labels)}.", 'yellow'), f_out)

    if not rows:
        print_dual(color_text("\n[ERROR] No valid grid results found.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    if lattice is None:
        print_dual(color_text(
            "\n[ERROR] Could not read a lattice from any grid point's structure.fdf -- "
            "cannot compute Cartesian shift distances for the gamma-surface plot.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_section('[2] STACKING FAULT ANALYSIS', f_out)
    print_dual(f"[INFO] Interface area (in-plane): {area:.4f} Ang^2 (used to convert relative "
                "energies to meV/Ang^2)", f_out)
    print_dual(f"Equilibrium stacking (min)    : {min_row.label}  (shift "
                f"{min_row.shift_x:.4f}, {min_row.shift_y:.4f})  "
                f"(E = {e_min:.6f} eV, dE = 0.000 meV/Ang^2)", f_out)
    print_dual(f"Highest-energy registry (max) : {max_row.label}  (shift "
                f"{max_row.shift_x:.4f}, {max_row.shift_y:.4f})  "
                f"(E = {e_max:.6f} eV, dE = {corrugation_mev:.3f} meV/Ang^2)", f_out)
    print_dual(f"Corrugation (stacking-fault) energy : {corrugation_mev:.3f} meV/Ang^2 "
                f"(max - min over the sampled grid; {corrugation:.6f} eV over a {area:.4f} "
                "Ang^2 cell)", f_out)
    if mode == 3:
        print_dual(color_text(
            "[NOTE] The interlayer gap was fixed across the whole grid (rigid-shift protocol, "
            "--mode 3) -- high-energy registries may have a non-relaxed interlayer distance; "
            "the true corrugation with a fully relaxed distance at each point may be somewhat "
            "lower (see --mode 1/2 in stb-stackingfault).", 'yellow'), f_out)
    elif mode == 2:
        print_dual(
            "[NOTE] The interlayer gap was relaxed per point with MACE-MP-0 before the SIESTA "
            "single-point (--mode 2) -- a cheaper but different level of theory than the "
            "reported DFT energy; see --mode 1 for a real-SIESTA-relaxed gap.", f_out)
    final_gaps = [r.final_gap for r in rows if r.final_gap is not None]
    if mode == 1 and final_gaps:
        print_dual(f"[INFO] SIESTA-relaxed gap ranged {min(final_gaps):.3f} - "
                    f"{max(final_gaps):.3f} Ang across the grid.", f_out)
    if n_skipped:
        print_dual(color_text(
            f"[NOTE] {n_skipped} grid point(s) were skipped -- the reported min/max above are "
            "only over the points actually computed; the true equilibrium stacking or highest-"
            "energy registry may be among the missing ones.", 'yellow'), f_out)

    print_section('[3] SUMMARY', f_out)
    print_dual(f"Grid points analyzed : {len(rows)} (skipped: {n_skipped})", f_out)

    # A gnuplot pm3d map needs EVERY point of the rectangular grid present
    # to render at all -- verified live: even a single missing/NaN point
    # in an otherwise-complete N x N block doesn't just leave a local gap,
    # it makes pm3d render the WHOLE plot blank (no error, no warning from
    # gnuplot itself -- a silent, easy-to-miss failure). A 1D --scan x/y/xy
    # line plot has the same all-or-nothing requirement for a different
    # reason: a gap in the middle of `plot ... with lines` just draws a
    # straight connector across the missing point, silently hiding it
    # rather than rendering blank -- still wrong, so the same "only plot
    # when complete" gate applies to both. Otherwise the numeric analysis
    # above still stands (it tolerates a partial grid/scan fine), but the
    # plot is skipped with a clear reason instead of silently writing a
    # broken/misleading one.
    n_expected = grid_nx * grid_ny if scan == "surface" else scan_points
    grid_complete = len(rows) == n_expected
    if args.save_gnuplot:
        if grid_complete:
            plot_dir = os.path.join(args.dir, "plot")
            os.makedirs(plot_dir, exist_ok=True)
            dat_path = os.path.join(plot_dir, os.path.basename(args.output))
            if scan == "surface":
                write_surface_data(dat_path, rows, grid_nx, grid_ny, lattice, area, e_min)
                plane_angle = check_planar_orthogonality(lattice, axis_idx=2)
                cb_range = (0.0, corrugation_mev) if corrugation_mev > 0 else None
                write_gnuplot_script(dat_path, dat_path, mode='slice',
                                      quantity_label="Stacking Fault Energy", is_signed=False,
                                      axis_idx=2, plane_angle_deg=plane_angle, cb_range=cb_range,
                                      contour=True, generator_name="stb-stackingfaultAnalysis",
                                      units="meV/A^2")
                plot_label = "Surface data"
            else:
                scan_desc = {"x": "shift_x", "y": "shift_y", "xy": "the shift_x=shift_y diagonal"}
                write_scan_line_data(dat_path, rows, lattice, area, e_min)
                write_gnuplot_script(dat_path, dat_path, mode='profile',
                                      quantity_label="Stacking Fault Energy", is_signed=False,
                                      generator_name="stb-stackingfaultAnalysis", units="meV/A^2",
                                      title=f"Stacking-Fault Energy along {scan_desc[scan]}",
                                      xlabel="Distance along scan path (Angstrom)")
                plot_label = "Line data"
            gplot_path = dat_path.rsplit('.', 1)[0] + ".gplot"
            print_dual(f"{color_text('[Saved]', 'cyan')} {plot_label} -> {dat_path}, {gplot_path} "
                        f"(cd {plot_dir} && gnuplot {os.path.basename(gplot_path)})", f_out)
        elif scan == "surface":
            print_dual(color_text(
                f"[WARNING] Skipping the gamma-surface plot: {n_skipped} grid point(s) are missing "
                f"out of {n_expected} -- a gnuplot pm3d map needs every point of the grid to "
                "render correctly (a single missing point can make the ENTIRE map render blank, "
                "not just a local gap). Complete SIESTA in the remaining folders and re-run to get "
                "the map.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                f"[WARNING] Skipping the line plot: {n_skipped} grid point(s) are missing out of "
                f"{n_expected} -- a gap in the scan would just draw a straight connector across the "
                "missing point instead of showing it's missing. Complete SIESTA in the remaining "
                "folders and re-run to get the plot.", 'yellow'), f_out)
    else:
        print_dual("Gnuplot data+script : not written (off by default -- pass --save-gnuplot to "
                    "write it, under '<dir>/plot/').", f_out)

    anim_path = os.path.join(args.dir, "stackingfault_animation.xyz")
    anim_frames, anim_unreadable = build_animation_frames(args.dir, rows, lattice, scan, area, e_min)
    if anim_frames:
        # format= is REQUIRED -- ASE's bare-".xyz"-extension auto-detect
        # defaults to plain (non-extended) XYZ, which drops per-frame
        # Lattice/.info -- same gotcha ani2traj.py's own
        # ase_write(..., format=ase_format) call works around.
        ase_write(anim_path, anim_frames, format="extxyz")
        print_dual(f"{color_text('[Saved]', 'cyan')} Animation     -> {anim_path} "
                    f"({len(anim_frames)} frame(s)).", f_out)
        if anim_unreadable:
            print_dual(color_text(
                f"[NOTE] {len(anim_unreadable)} grid point(s) omitted from the animation "
                f"(structure unreadable): {', '.join(anim_unreadable)}.", 'yellow'), f_out)
    else:
        print_dual(color_text("[WARNING] No readable structures -- animation not written.",
                               'yellow'), f_out)

    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report        -> {report_path}", f_out)

    if args.apply:
        print_section('[4] APPLY', f_out)
        src = os.path.join(args.dir, "positions", min_row.label, "structure.fdf")
        try:
            shutil.copy(src, args.apply)
        except OSError as e:
            print_dual(color_text(f"[ERROR] Could not copy '{src}' to '{args.apply}': {e}", 'red'), f_out)
        else:
            print_dual(f"{color_text('[Applied]', 'green')} {min_row.label} -> {args.apply}", f_out)

    if f_out:
        f_out.close()

    if args.view:
        if scan == "surface":
            view_surface_plot(rows, grid_nx, grid_ny, area, e_min)
        else:
            view_scan_line_plot(rows, lattice, area, e_min, scan)
    if args.view_animation:
        if anim_frames:
            print(f"\n{color_text('--view-animation:', 'cyan')} opening {len(anim_frames)} "
                  "frame(s) in ASE's interactive viewer (use the frame slider):")
            view_structure_interactive(anim_frames)
        else:
            print(color_text("[FAIL] No readable structures to view.", 'red'))


if __name__ == "__main__":
    main()
