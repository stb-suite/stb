#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import shutil
import argparse
from datetime import datetime
import numpy as np
from stb.core import siesta_log
from stb.core.siesta_log import check_scf_and_force
from stb.core.cli import color_text, show_intro, print_dual

REPORT_FILE = "neb_report.txt"
SETUP_REPORT_FILE = "neb_setup.txt"


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


def fit_spline_barrier(rows):
    """Fits a cubic spline (scipy's default not-a-knot boundary condition,
    not a "natural" spline -- not-a-knot avoids forcing zero curvature at
    the two endpoints, which a natural spline would, biasing the fitted
    curve near the edges) through the (reaction_coord, energy) points
    (sorted by reaction_coord) and returns
    (ts_coord, spline_barrier, spline_dE) -- a smoother, less
    grid-dependent transition-state estimate than "whichever discrete
    image happens to be highest", the DFT-side analogue of
    ase.mep.neb.NEBTools.get_barrier(fit=True) (already used by stb-neb
    --ml-neb for the MACE band -- that one additionally has per-image
    forces for a Hermite-style fit; here we only have a scalar SIESTA
    single-point energy per image, so a plain energy-only cubic spline is
    the fair equivalent). Returns None if there are fewer than 4 points
    (not enough for a meaningful cubic fit) or if two images share the
    same reaction coordinate (a degenerate/duplicate image -- see
    stb-neb's check_path_quality).

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
    if len(rows) < 4:
        return None
    rows_sorted = sorted(rows, key=lambda r: r.reaction_coord)
    x = [r.reaction_coord for r in rows_sorted]
    y = [r.energy for r in rows_sorted]
    if len(set(x)) != len(x):
        return None
    from scipy.interpolate import CubicSpline
    from scipy.optimize import minimize_scalar
    spline = CubicSpline(x, y)

    x_dense = np.linspace(x[0], x[-1], 500)
    y_dense = spline(x_dense)
    i_max = int(np.argmax(y_dense))
    lo = x_dense[max(i_max - 1, 0)]
    hi = x_dense[min(i_max + 1, len(x_dense) - 1)]
    result = minimize_scalar(lambda xv: -spline(xv), bounds=(lo, hi), method='bounded')
    ts_coord = float(result.x)
    ts_energy = float(spline(ts_coord))
    return ts_coord, ts_energy - y[0], y[-1] - y[0]


def write_curve_plot(dat_path, rows, use_reaction_coord):
    """Writes <dat_path> plus a companion .gplot: energy vs. reaction
    coordinate (or plain image index, if reaction_coord wasn't available)
    -- a continuous path, so linespoints, not a scatter (unlike
    adsorb_analysis.py's per-site curve, which is a discrete set of
    candidates). Same <name>.dat + <name>.gplot convention as the rest of
    the suite.
    """
    rows_sorted = sorted(rows, key=lambda r: r.index)
    xlabel_col = "ReactionCoord(Ang)" if use_reaction_coord else "ImageIndex"
    with open(dat_path, 'w') as f:
        f.write("# NEB energy profile\n")
        f.write(f"# 1:{xlabel_col} 2:E(eV) 3:Label\n")
        for r in rows_sorted:
            x = r.reaction_coord if use_reaction_coord else r.index
            f.write(f"{x:.6f}  {r.energy:.6f}  {r.label}\n")

    base = os.path.splitext(dat_path)[0]
    base_name = os.path.basename(base)
    gplot_path = f"{base}.gplot"
    dat_name = os.path.basename(dat_path)
    xlabel = "Reaction coordinate (Ang)" if use_reaction_coord else "Image index"
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-nebAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            'set title "NEB energy profile"\n',
            f'set xlabel "{xlabel}"\n',
            'set ylabel "E (eV)"\n',
            'set grid\n',
            'set key top right\n',
            f'plot "{dat_name}" using 1:2 with linespoints lw 2 pt 7 lc rgb "#2255cc" title "E"\n',
        ])
    return gplot_path


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes a NEB energy profile from an stb-neb band: reads "
        "each image_NN/'s SIESTA single-point energy and estimates the reaction barrier from "
        "the highest-energy image.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --dir . --file calc.out\n"
               "  %(prog)s --dir . --apply ts_guess.fdf\n"
    )

    parser.add_argument("--dir", type=str, default=".",
                         help="Root directory containing every 'image_NN/' (default: current directory).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-o", "--output", type=str, default="neb_curve.dat",
                         help="Output data file name (default: neb_curve.dat).")
    parser.add_argument("--apply", type=str, default=None, metavar="STRUCTURE_FDF",
                         help="Copy the highest-energy image's structure.fdf (the approximate "
                              "transition-state guess) to this path.")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05, same as "
                              "stb-adsorbAnalysis/stb-cohesiveAnalysis) above which an image's "
                              "calc.out is flagged as possibly not relaxed. Advisory only, never "
                              "blocks the result.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nebAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

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

    image_rows, ml_neb_used = read_image_table(args.dir)
    table_found = image_rows is not None
    if not table_found:
        dirs = sorted(d for d in os.listdir(args.dir)
                       if os.path.isdir(os.path.join(args.dir, d)) and d.startswith("image_"))
        image_rows = [(d, i, None) for i, d in enumerate(dirs)]

    if not image_rows:
        print(color_text(f"[ERROR] No 'image_*' folders found in '{args.dir}'. Did you run stb-neb?", 'red'))
        sys.exit(1)

    with open(REPORT_FILE, 'w') as f_out:
        print_dual(f"{color_text('===== NEB ENERGY PROFILE REPORT =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Directory  : {args.dir}", f_out)
        print_dual(f"Output file: {args.file}", f_out)
        if not table_found:
            print_dual(color_text(
                "[NOTE] No 'neb_setup.txt' image table found -- falling back to a sorted glob of "
                "'image_*' folders (index = sort position, no reaction-coordinate data available; "
                "the x-axis below uses plain image index instead).", 'yellow'), f_out)

        print_dual(f"\n{color_text('[1] IMAGE ENERGIES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        rows = []
        n_skipped = 0
        scf_warn_labels = []
        force_warn_labels = []
        header = f"{'Image':<14}{'Index':<7}{'ReactionCoord':<16}{'E(eV)':<16}{'SCF':<6}{'MaxF(eV/A)':<12}"
        print_dual(header, f_out)
        print_dual("-" * len(header), f_out)
        for label, index, reaction_coord in image_rows:
            image_dir = os.path.join(args.dir, label)
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
            sys.exit(1)

        print_dual(f"\n{color_text('[2] BARRIER ANALYSIS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        max_index = max(idx for _, idx, _ in image_rows)
        initial_row = next((r for r in rows if r.index == 0), None)
        final_row = next((r for r in rows if r.index == max_index), None)
        if initial_row is None or final_row is None:
            print_dual(color_text(
                "[ERROR] Could not read the initial and/or final endpoint image's energy -- "
                "cannot compute a barrier.", 'red'), f_out)
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

        use_reaction_coord = all(r.reaction_coord is not None for r in rows)
        spline_fit = fit_spline_barrier(rows) if use_reaction_coord else None
        if spline_fit:
            ts_coord, spline_barrier, spline_dE = spline_fit
            print_dual(f"Spline-fitted barrier (smoothed)  : {spline_barrier:.6f} eV at reaction "
                        f"coordinate {ts_coord:.4f} Ang (cubic spline through the {len(rows)} "
                        "energies)", f_out)
        elif use_reaction_coord and len(rows) < 4:
            print_dual("[NOTE] Not enough images (need >= 4) for a spline-fitted barrier estimate.", f_out)

        if ml_neb_used:
            caveat = ("Stage 1 already converged a real climbing-image band on MACE-MP-0 first, "
                       "so this is a stronger (ML-relaxed-path) estimate -- but still not an "
                       "independent DFT confirmation of a converged saddle point.")
        else:
            caveat = ("Approximate, interpolated-path estimate: no real reaction coordinate was "
                       "optimized between these DFT single points (no inter-image spring "
                       "coupling at the DFT level) -- the true saddle point may lie off this "
                       "fixed path. Re-run stb-neb with --ml-neb for a physically relaxed band "
                       "shape before the DFT single points.")
        print_dual(color_text(f"[NOTE] {caveat}", 'yellow'), f_out)

        print_dual(f"\n{color_text('[3] SUMMARY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Images analyzed : {len(rows)} (skipped: {n_skipped})", f_out)

        gplot_path = write_curve_plot(args.output, rows, use_reaction_coord)
        print_dual(f"{color_text('[Saved]', 'cyan')} Curve data -> {args.output}, {gplot_path} "
                    f"(cd {os.path.dirname(args.output) or '.'} && gnuplot {os.path.basename(gplot_path)})",
                    f_out)
        print_dual(f"{color_text('[Saved]', 'cyan')} Report     -> {REPORT_FILE}", f_out)

        if args.apply:
            print_dual(f"\n{color_text('[4] APPLY', 'magenta')}", f_out)
            print_dual("-" * 60, f_out)
            src = os.path.join(args.dir, max_row.label, "structure.fdf")
            try:
                shutil.copy(src, args.apply)
            except OSError as e:
                print_dual(color_text(f"[ERROR] Could not copy '{src}' to '{args.apply}': {e}", 'red'), f_out)
            else:
                print_dual(f"{color_text('[Applied]', 'green')} {max_row.label} -> {args.apply}", f_out)


if __name__ == "__main__":
    main()
