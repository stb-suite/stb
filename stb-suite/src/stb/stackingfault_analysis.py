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
import shutil
import argparse
from datetime import datetime
import numpy as np
from stb.core import siesta_log, structure_io
from stb.core.siesta_log import check_scf_and_force
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.grid_export import write_gnuplot_script, check_planar_orthogonality

REPORT_FILE = "stackingfault_report.txt"
SETUP_REPORT_FILE = "stackingfault_setup.txt"
_SHIFT_DIR_RE = re.compile(r'^shift_(\d+)_(\d+)$')


def read_grid_table(root_dir):
    """Parses the '# GRID_TABLE' section stb-stackingfault writes in
    stackingfault_setup.txt into a list of (label, i, j, shift_x,
    shift_y), or None if the report (or its table section) isn't found --
    the caller then falls back to a sorted glob of 'shift_II_JJ' folders
    (see fallback_grid_rows()). Mirrors adsorb_analysis.py::
    read_site_table / neb_analysis.py::read_image_table -- the grid
    directory is recomputed by the caller as os.path.join(root_dir,
    label), not trusted from the table's own 'dir' column.
    """
    report_path = os.path.join(root_dir, SETUP_REPORT_FILE)
    if not os.path.isfile(report_path):
        return None
    rows = []
    in_table = False
    with open(report_path) as f:
        for line in f:
            if line.startswith("# GRID_TABLE"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            label, i_str, j_str, x_str, y_str = parts[0], parts[1], parts[2], parts[3], parts[4]
            try:
                i, j = int(i_str), int(j_str)
                shift_x, shift_y = float(x_str), float(y_str)
            except ValueError:
                continue
            rows.append((label, i, j, shift_x, shift_y))
    return rows or None


def fallback_grid_rows(root_dir):
    """Recovers the grid layout from a sorted glob of 'shift_II_JJ'
    folders when stackingfault_setup.txt is missing -- unlike a generic
    "no metadata at all" fallback, the folder-naming convention here IS
    the (i, j) grid index, and stb-stackingfault always builds shifts via
    linspace(0, 1, grid_n, endpoint=False), so shift_x/shift_y are exactly
    recoverable as i/grid_n, j/grid_n once grid_n is inferred from the max
    index actually present -- no data is lost, unlike stb-neb's index-only
    fallback (which can't recover reaction_coord at all).
    """
    dirs = [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]
    parsed = []
    for d in dirs:
        m = _SHIFT_DIR_RE.match(d)
        if m:
            parsed.append((d, int(m.group(1)), int(m.group(2))))
    if not parsed:
        return None
    grid_n = max(max(i, j) for _, i, j in parsed) + 1
    parsed.sort(key=lambda r: (r[1], r[2]))
    return [(d, i, j, i / grid_n, j / grid_n) for d, i, j in parsed]


class GridRow:
    """One analyzed grid point: label, (i, j) grid index, fractional
    shift, and its computed energy/quality diagnostics. Plain
    attribute-holder, same style as adsorb_analysis.py's SiteRow /
    neb_analysis.py's ImageRow.
    """
    def __init__(self, label, i, j, shift_x, shift_y, energy, scf_ok, max_force):
        self.label = label
        self.i = i
        self.j = j
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.energy = energy
        self.scf_ok = scf_ok
        self.max_force = max_force


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


def write_surface_data(dat_path, rows, grid_n, lattice):
    """Writes the gamma-surface data file in gnuplot PM3D block format
    (blank line after each row of constant i) with columns
    dx(Ang) dy(Ang) 0.0(unused) E_rel(eV) -- matches
    core/grid_export.py::write_data_file's 'slice'-mode convention exactly
    (4 columns, "1:2:4" used by write_gnuplot_script's axis_idx=2 branch),
    so that function can be reused verbatim for the heatmap script instead
    of writing a new one.
    """
    by_ij = {(r.i, r.j): r for r in rows}
    e_min = min(r.energy for r in rows)
    with open(dat_path, 'w') as f:
        f.write("# Stacking-fault (gamma-surface) energy landscape\n")
        f.write("# 1:dx(Ang) 2:dy(Ang) 3:(unused) 4:E_rel(eV)\n")
        for i in range(grid_n):
            for j in range(grid_n):
                row = by_ij.get((i, j))
                if row is None:
                    continue
                dx, dy = cartesian_shift(lattice, row.shift_x, row.shift_y)
                f.write(f"{dx:.6f} {dy:.6f} 0.0 {row.energy - e_min:.6f}\n")
            f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes a stacking-fault (gamma-surface) energy landscape "
        "from an stb-stackingfault grid: reads each shift_II_JJ/'s SIESTA single-point energy and "
        "reports the equilibrium stacking, the highest-energy registry, and the corrugation "
        "energy.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --dir . --file calc.out\n"
               "  %(prog)s --dir . --apply equilibrium.fdf\n"
    )

    parser.add_argument("--dir", type=str, default=".",
                         help="Root directory containing every 'shift_II_JJ/' (default: current directory).")
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
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
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

    grid_rows = read_grid_table(args.dir)
    table_found = grid_rows is not None
    if not table_found:
        grid_rows = fallback_grid_rows(args.dir)

    if not grid_rows:
        print(color_text(f"[ERROR] No 'shift_II_JJ' folders found in '{args.dir}'. Did you run "
                          "stb-stackingfault?", 'red'))
        sys.exit(1)

    grid_n = max(max(i, j) for _, i, j, _, _ in grid_rows) + 1

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== STACKING-FAULT (GAMMA-SURFACE) REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory  : {args.dir}", f_out)
    print_dual(f"Output file: {args.file}", f_out)
    print_dual(f"Grid       : {grid_n} x {grid_n}", f_out)
    if not table_found:
        print_dual(color_text(
            "[NOTE] No 'stackingfault_setup.txt' grid table found -- recovered the grid "
            "layout from 'shift_II_JJ' folder names instead (shift_x/shift_y reconstructed "
            "as i/grid_n, j/grid_n).", 'yellow'), f_out)

    print_section('[1] GRID ENERGIES', f_out)
    rows = []
    n_skipped = 0
    scf_warn_labels = []
    force_warn_labels = []
    lattice = None
    header = f"{'Point':<16}{'Shift(x,y)':<20}{'E(eV)':<16}{'SCF':<6}{'MaxF(eV/A)':<12}"
    print_dual(header, f_out)
    print_dual("-" * len(header), f_out)
    for label, i, j, shift_x, shift_y in grid_rows:
        grid_dir = os.path.join(args.dir, label)
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

        rows.append(GridRow(label, i, j, shift_x, shift_y, energy, scf_ok, max_force))
        shift_str = f"({shift_x:.3f},{shift_y:.3f})"
        scf_str = color_text("WARN", 'yellow') if not scf_ok else "OK"
        force_str = f"{max_force:.4f}" if max_force is not None else "--"
        print_dual(f"{label:<16}{shift_str:<20}{energy:<16.6f}{scf_str:<6}{force_str:<12}", f_out)
    print_dual("-" * len(header), f_out)
    if scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(scf_warn_labels)} grid point(s) never confirmed SCF convergence "
            f"-- their energy may be unreliable: {', '.join(scf_warn_labels)}.", 'yellow'), f_out)
    if force_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(force_warn_labels)} grid point(s) have residual force above "
            f"--force-tolerance ({args.force_tolerance} eV/Ang), possibly not single-point-"
            f"converged: {', '.join(force_warn_labels)}.", 'yellow'), f_out)

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
    min_row = min(rows, key=lambda r: r.energy)
    max_row = max(rows, key=lambda r: r.energy)
    corrugation = max_row.energy - min_row.energy

    print_dual(f"Equilibrium stacking (min)    : {min_row.label}  (shift "
                f"{min_row.shift_x:.4f}, {min_row.shift_y:.4f})  "
                f"(E = {min_row.energy:.6f} eV)", f_out)
    print_dual(f"Highest-energy registry (max) : {max_row.label}  (shift "
                f"{max_row.shift_x:.4f}, {max_row.shift_y:.4f})  "
                f"(E = {max_row.energy:.6f} eV)", f_out)
    print_dual(f"Corrugation (stacking-fault) energy : {corrugation:.6f} eV "
                "(max - min over the sampled grid)", f_out)
    print_dual(color_text(
        "[NOTE] The interlayer gap was fixed across the whole grid (rigid-shift protocol) -- "
        "high-energy registries may have a non-relaxed interlayer distance; the true "
        "corrugation with a fully relaxed distance at each point may be somewhat lower.",
        'yellow'), f_out)
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
    # gnuplot itself -- a silent, easy-to-miss failure). So the surface
    # plot is only attempted when every point was actually computed;
    # otherwise the numeric analysis above still stands (it tolerates a
    # partial grid fine), but the map is skipped with a clear reason
    # instead of silently writing a broken/blank one.
    if len(rows) == grid_n * grid_n:
        write_surface_data(args.output, rows, grid_n, lattice)
        plane_angle = check_planar_orthogonality(lattice, axis_idx=2)
        e_max_rel = corrugation
        cb_range = (0.0, e_max_rel) if e_max_rel > 0 else None
        write_gnuplot_script(args.output, args.output, mode='slice',
                              quantity_label="Stacking Fault Energy", is_signed=False,
                              axis_idx=2, plane_angle_deg=plane_angle, cb_range=cb_range,
                              contour=True, generator_name="stb-stackingfaultAnalysis", units="eV")
        gplot_path = args.output.rsplit('.', 1)[0] + ".gplot"
        print_dual(f"{color_text('[Saved]', 'cyan')} Surface data -> {args.output}, {gplot_path} "
                    f"(cd {os.path.dirname(args.output) or '.'} && gnuplot {os.path.basename(gplot_path)})",
                    f_out)
    else:
        print_dual(color_text(
            f"[WARNING] Skipping the gamma-surface plot: {n_skipped} grid point(s) are missing "
            f"out of {grid_n * grid_n} -- a gnuplot pm3d map needs every point of the grid to "
            "render correctly (a single missing point can make the ENTIRE map render blank, "
            "not just a local gap). Complete SIESTA in the remaining folders and re-run to get "
            "the map.", 'yellow'), f_out)
    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report       -> {report_path}", f_out)

    if args.apply:
        print_section('[4] APPLY', f_out)
        src = os.path.join(args.dir, min_row.label, "structure.fdf")
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
