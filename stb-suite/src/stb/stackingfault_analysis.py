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
import json
import shutil
import argparse
from datetime import datetime
import numpy as np
from stb.core import siesta_log, structure_io
from stb.core.siesta_log import check_scf_and_force
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.grid_export import write_gnuplot_script, check_planar_orthogonality

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
    `energy_corrected` is None unless stb-stackingfaultBsse's own
    'bsse/<label>/' folders were found AND both ghost-fragment energies
    were readable for this point -- see read_bsse_energy.
    """
    def __init__(self, label, i, j, shift_x, shift_y, energy, scf_ok, max_force, final_gap=None,
                 energy_corrected=None):
        self.label = label
        self.i = i
        self.j = j
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.energy = energy
        self.scf_ok = scf_ok
        self.max_force = max_force
        self.final_gap = final_gap
        self.energy_corrected = energy_corrected

    @property
    def energy_used(self):
        """The BSSE-corrected energy when available, else the raw one --
        what equilibrium/corrugation should actually be computed from."""
        return self.energy_corrected if self.energy_corrected is not None else self.energy


def read_bsse_energy(bsse_dir, file_name):
    """Reads the two ghost-fragment reference energies stb-stackingfaultBsse
    writes under 'bsse/<label>/' -- a tree parallel to the grid's own
    'shift_II_JJ/' folders, not nested inside them -- '<bsse_dir>/
    bsse_layer1/<file_name>' and '<bsse_dir>/bsse_layer2/<file_name>' --
    and returns (e_bsse_layer1, e_bsse_layer2), or (None, None) if either
    is missing/unreadable. Best-effort per point (never blocks the
    uncorrected result for that point or any other) -- same contract as
    adsorb_analysis.py's own read_bsse_energy, adapted from
    bsse_slab/bsse_adsorbate to bsse_layer1/bsse_layer2. Callers are
    expected to have already confirmed '<bsse_dir>/bsse_layer1' exists
    before calling this, same convention as the adsorption sibling.
    """
    layer1_path = os.path.join(bsse_dir, "bsse_layer1", file_name)
    layer2_path = os.path.join(bsse_dir, "bsse_layer2", file_name)
    e_layer1 = siesta_log.get_free_energy(layer1_path)
    e_layer2 = siesta_log.get_free_energy(layer2_path)
    if e_layer1 is None or e_layer2 is None:
        return None, None
    return e_layer1, e_layer2


def effective_gap(pmg_structure, n_layer1):
    """Cartesian z-gap (Ang) between layer 1's topmost atom and layer 2's
    bottommost atom -- same definition/formula as stackingfault.py's own
    effective_gap (duplicated here rather than cross-imported: Stage 1 and
    Stage 2 are meant to stay self-contained, same "duplicated per
    self-contained workflow stage" convention documented for
    read_site_theory_flags elsewhere in the suite). Used to report the
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


def write_surface_data(dat_path, rows, grid_nx, grid_ny, lattice):
    """Writes the gamma-surface data file in gnuplot PM3D block format
    (blank line after each row of constant i) with columns
    dx(Ang) dy(Ang) 0.0(unused) E_rel(eV) -- matches
    core/grid_export.py::write_data_file's 'slice'-mode convention exactly
    (4 columns, "1:2:4" used by write_gnuplot_script's axis_idx=2 branch),
    so that function can be reused verbatim for the heatmap script instead
    of writing a new one. grid_nx/grid_ny may differ (an asymmetric grid).
    """
    by_ij = {(r.i, r.j): r for r in rows}
    e_min = min(r.energy for r in rows)
    with open(dat_path, 'w') as f:
        f.write("# Stacking-fault (gamma-surface) energy landscape\n")
        f.write("# 1:dx(Ang) 2:dy(Ang) 3:(unused) 4:E_rel(eV)\n")
        for i in range(grid_nx):
            for j in range(grid_ny):
                row = by_ij.get((i, j))
                if row is None:
                    continue
                dx, dy = cartesian_shift(lattice, row.shift_x, row.shift_y)
                f.write(f"{dx:.6f} {dy:.6f} 0.0 {row.energy - e_min:.6f}\n")
            f.write("\n")


def write_scan_line_data(dat_path, rows, lattice):
    """Writes a --scan x/y/xy line-plot data file: 2 plain columns,
    distance traveled along the scan path (Ang, from the origin) and
    E_rel (eV) -- the 1D analog of write_surface_data's PM3D block, read
    by write_gnuplot_script's mode='profile' template ('plot ... using
    1:2'). Distance is the Cartesian norm of cartesian_shift(shift_x,
    shift_y) -- works uniformly for all 3 scan shapes (shift_y is always 0
    for 'x', shift_x is always 0 for 'y', so the norm collapses to the
    single swept component in both cases; for 'xy' it's the true diagonal
    distance). Rows are sorted by this distance (not just manifest order)
    so the line always plots left-to-right even if the caller's row order
    doesn't already match it.
    """
    e_min = min(r.energy for r in rows)
    points = []
    for r in rows:
        dx, dy = cartesian_shift(lattice, r.shift_x, r.shift_y)
        distance = float(np.hypot(dx, dy))
        points.append((distance, r.energy - e_min))
    points.sort(key=lambda p: p[0])
    with open(dat_path, 'w') as f:
        f.write("# Stacking-fault energy along a 1D scan path\n")
        f.write("# 1:distance(Ang) 2:E_rel(eV)\n")
        for distance, e_rel in points:
            f.write(f"{distance:.6f} {e_rel:.6f}\n")


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
    )

    parser.add_argument("--dir", type=str, default=_default_analysis_dir(),
                         help="Root directory containing 'positions/shift_II_JJ/' (default: "
                              "'sf_run', the self-contained run folder stb-stackingfault Stage 1 "
                              "always writes its output into -- auto-detected as '.' instead when "
                              "the current directory already has its own sf_manifest.json).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-o", "--output", type=str, default="stackingfault_surface.dat",
                         help="Output data file name (default: stackingfault_surface.dat).")
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
    bsse_root = os.path.join(args.dir, "bsse")
    show_bsse_column = os.path.isdir(bsse_root)
    print_dual(f"BSSE       : {'stb-stackingfaultBsse folder found -- will use where available' if show_bsse_column else 'not run (see stb-stackingfaultBsse, Stage 2)'}",
                f_out)
    if not manifest_found:
        print_dual(color_text(
            f"[NOTE] No '{MANIFEST_FILE}' found -- recovered the grid layout from "
            "'shift_II_JJ' folder names instead (shift_x/shift_y reconstructed as "
            "i/grid_nx, j/grid_ny).", 'yellow'), f_out)

    print_section('[1] GRID ENERGIES', f_out)
    rows = []
    n_skipped = 0
    scf_warn_labels = []
    force_warn_labels = []
    lattice = None
    show_final_gap = mode == 1 and n_layer1_atoms is not None
    bsse_pending_labels = []  # bsse/<label>/ exists but energies aren't readable yet
    header = f"{'Point':<16}{'Shift(x,y)':<20}{'E(eV)':<16}{'SCF':<6}{'MaxF(eV/A)':<12}"
    if show_final_gap:
        header += f"{'GapFinal(A)':<12}"
    if show_bsse_column:
        header += f"{'BSSE':<9}"
    print_dual(header, f_out)
    print_dual("-" * len(header), f_out)
    for label, i, j, shift_x, shift_y in grid_rows:
        grid_dir = os.path.join(args.dir, "positions", label)
        out_path = os.path.join(grid_dir, args.file)
        if not os.path.exists(out_path):
            n_skipped += 1
            print_dual(f"{label:<16}{color_text('SKIP', 'yellow')} (missing {args.file})", f_out)
            continue
        energy = siesta_log.get_free_energy(out_path)
        if energy is None:
            n_skipped += 1
            print_dual(f"{label:<16}{color_text('SKIP', 'yellow')} (could not parse energy)", f_out)
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

        energy_corrected = None
        bsse_str = "--"
        if show_bsse_column:
            point_bsse_dir = os.path.join(bsse_root, label)
            if os.path.isdir(os.path.join(point_bsse_dir, "bsse_layer1")):
                e_layer1, e_layer2 = read_bsse_energy(point_bsse_dir, args.file)
                if e_layer1 is not None and e_layer2 is not None:
                    energy_corrected = energy - e_layer1 - e_layer2
                    bsse_str = "yes"
                else:
                    bsse_pending_labels.append(label)
                    bsse_str = "pending"

        rows.append(GridRow(label, i, j, shift_x, shift_y, energy, scf_ok, max_force, final_gap,
                             energy_corrected))
        shift_str = f"({shift_x:.3f},{shift_y:.3f})"
        scf_str = color_text("WARN", 'yellow') if not scf_ok else "OK"
        force_str = f"{max_force:.4f}" if max_force is not None else "--"
        row_line = f"{label:<16}{shift_str:<20}{energy:<16.6f}{scf_str:<6}{force_str:<12}"
        if show_final_gap:
            row_line += f"{final_gap:<12.4f}" if final_gap is not None else f"{'--':<12}"
        if show_bsse_column:
            row_line += f"{bsse_str:<9}"
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
    if bsse_pending_labels:
        print_dual(color_text(
            f"[NOTE] {len(bsse_pending_labels)} grid point(s) have BSSE folders (stb-"
            f"stackingfaultBsse) but SIESTA hasn't finished the ghost-fragment calculations yet "
            f"-- using uncorrected energies for those: {', '.join(bsse_pending_labels)}.",
            'yellow'), f_out)

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
    bsse_count = sum(1 for r in rows if r.energy_corrected is not None)
    use_bsse = bsse_count == len(rows) and bsse_count > 0
    min_row = min(rows, key=lambda r: r.energy_used) if use_bsse else min(rows, key=lambda r: r.energy)
    max_row = max(rows, key=lambda r: r.energy_used) if use_bsse else max(rows, key=lambda r: r.energy)
    e_min = min_row.energy_used if use_bsse else min_row.energy
    e_max = max_row.energy_used if use_bsse else max_row.energy
    corrugation = e_max - e_min

    bsse_tag = " (BSSE-corrected)" if use_bsse else ""
    print_dual(f"Equilibrium stacking (min)    : {min_row.label}  (shift "
                f"{min_row.shift_x:.4f}, {min_row.shift_y:.4f})  "
                f"(E = {e_min:.6f} eV{bsse_tag})", f_out)
    print_dual(f"Highest-energy registry (max) : {max_row.label}  (shift "
                f"{max_row.shift_x:.4f}, {max_row.shift_y:.4f})  "
                f"(E = {e_max:.6f} eV{bsse_tag})", f_out)
    print_dual(f"Corrugation (stacking-fault) energy : {corrugation:.6f} eV{bsse_tag} "
                "(max - min over the sampled grid)", f_out)
    if show_bsse_column:
        if use_bsse:
            print_dual(color_text(
                "[INFO] BSSE (counterpoise) correction applied at every grid point analyzed -- "
                "the values above already account for it.", 'green'), f_out)
        elif bsse_count > 0:
            print_dual(color_text(
                f"[NOTE] BSSE correction is only available for {bsse_count} of {len(rows)} "
                "analyzed grid point(s) -- using UNCORRECTED energies above for consistency "
                "(a mix of corrected/uncorrected points would not be a meaningful comparison). "
                "Finish stb-stackingfaultBsse's ghost-fragment SIESTA runs for the remaining "
                "points and re-run for the corrected result.", 'yellow'), f_out)
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
    if len(rows) == n_expected:
        if scan == "surface":
            write_surface_data(args.output, rows, grid_nx, grid_ny, lattice)
            plane_angle = check_planar_orthogonality(lattice, axis_idx=2)
            e_max_rel = corrugation
            cb_range = (0.0, e_max_rel) if e_max_rel > 0 else None
            write_gnuplot_script(args.output, args.output, mode='slice',
                                  quantity_label="Stacking Fault Energy", is_signed=False,
                                  axis_idx=2, plane_angle_deg=plane_angle, cb_range=cb_range,
                                  contour=True, generator_name="stb-stackingfaultAnalysis", units="eV")
            plot_label = "Surface data"
        else:
            scan_desc = {"x": "shift_x", "y": "shift_y", "xy": "the shift_x=shift_y diagonal"}
            write_scan_line_data(args.output, rows, lattice)
            write_gnuplot_script(args.output, args.output, mode='profile',
                                  quantity_label="Stacking Fault Energy", is_signed=False,
                                  generator_name="stb-stackingfaultAnalysis", units="eV",
                                  title=f"Stacking-Fault Energy along {scan_desc[scan]}",
                                  xlabel="Distance along scan path (Angstrom)")
            plot_label = "Line data"
        gplot_path = args.output.rsplit('.', 1)[0] + ".gplot"
        print_dual(f"{color_text('[Saved]', 'cyan')} {plot_label} -> {args.output}, {gplot_path} "
                    f"(cd {os.path.dirname(args.output) or '.'} && gnuplot {os.path.basename(gplot_path)})",
                    f_out)
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
    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report       -> {report_path}", f_out)

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


if __name__ == "__main__":
    main()
