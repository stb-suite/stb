#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

# 2.0.1: numbered-section report (matching the other Stage 2 analysis tools);
# --input-dir now defaults to 'xrd_search' (stb-xrdsearch's own default output
# folder); --raw-experimental peak-list auto-broadening fix (same bug already
# fixed in stb-xrd's own --compare-to); gnuplot .dat+.gplot output (ranking bar
# chart + best-candidate simulated/experimental overlay) via --save-gnuplot,
# written under its own 'plot/' subfolder; matplotlib preview (same two plots,
# one figure) via --view; space group added to the ranking table; final
# dedicated LIBRARY WARNINGS section.
VERSION = "2.0.1"

import os
import sys
import glob
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
from stb.core import structure_io
from stb.core import xrd as xrd_core
from stb.core.symmetry import space_group_label
from stb.core.cli import (color_text, show_intro, run_with_spinner, print_dual, print_section,
                          print_table, capture_library_noise)
from stb.core.deps import require_pyxtal

REPORT_FILE = "xrd_rank_report.txt"


def find_candidates(input_dir):
    """Finds one candidate structure per subfolder of input_dir, preferring a
    post-DFT relaxed '*.STRUCT_OUT' file (if the user has already run SIESTA
    there) over the original 'structure.fdf' written by stb-xrdsearch. Also
    matches flat .fdf files directly in input_dir, for candidate sets built
    by hand rather than by stb-xrdsearch.
    Returns a list of (name, path, fmt, source) tuples, fmt in {"fdf",
    "struct_out"} and source in {"raw", "relaxed"} (for reporting only).
    """
    candidates = []
    seen_names = set()
    for entry in sorted(os.listdir(input_dir)):
        folder = os.path.join(input_dir, entry)
        if not os.path.isdir(folder):
            continue
        struct_out = sorted(glob.glob(os.path.join(folder, "*.STRUCT_OUT")))
        if struct_out:
            candidates.append((entry, struct_out[0], "struct_out", "relaxed"))
            seen_names.add(entry)
            continue
        fdf_path = os.path.join(folder, "structure.fdf")
        if os.path.isfile(fdf_path):
            candidates.append((entry, fdf_path, "fdf", "raw"))
            seen_names.add(entry)
    for path in sorted(glob.glob(os.path.join(input_dir, "*.fdf"))):
        name = os.path.splitext(os.path.basename(path))[0]
        if name in seen_names:
            print(color_text(
                f"  Warning: skipping '{path}' -- its name collides with the '{name}' "
                "subfolder candidate above, and there's no way to tell them apart.", 'yellow'))
            continue
        candidates.append((name, path, "fdf", "raw"))
        seen_names.add(name)
    return candidates


def write_rank_gnuplot(dat_path, gplot_path, shown, wavelength_label):
    """Writes the ranking bar chart as a .dat (rank, similarity, name) +
    .gplot pair -- the gnuplot-workflow analog of the [3] RANKING table,
    matching this WORKFLOW_TOOLS stage's own sibling analysis tools
    (strain_analysis.py/elastic_analysis.py/convergence_analysis.py), which
    all pair a plain-text .dat with a .gplot script rather than the newer
    matplotlib convention the ML Simulations tools use.
    """
    with open(dat_path, "w") as f:
        f.write("# ============================================================\n")
        f.write("# STB-XRDRANK -- Candidate similarity ranking\n")
        f.write("# ============================================================\n")
        f.write(f"# Wavelength : {wavelength_label}\n")
        f.write(f"# Candidates : {len(shown)}\n")
        f.write("# ------------------------------------------------------------\n")
        f.write("# Columns: rank  similarity  name\n")
        for rank, (name, similarity, formula, source, sg) in enumerate(shown, start=1):
            f.write(f"{rank:<5d} {similarity:<12.6f} {name}\n")

    dat_basename = os.path.basename(dat_path)
    lines = [
        'set terminal pdfcairo enhanced font "Arial,12" size 10,6\n',
        'set output "xrd_rank.pdf"\n',
        'set title "XRD candidate ranking"\n',
        'set ylabel "Similarity" font "Arial,14"\n',
        'set yrange [0:1]\n',
        'set style data histograms\n',
        'set style fill solid 0.7 border -1\n',
        'set xtics rotate by -45 right\n',
        'set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
        'unset key\n',
        '\n',
        f'plot "{dat_basename}" using 2:xtic(3) with boxes lc rgb "navy"\n',
    ]
    with open(gplot_path, "w") as f:
        f.writelines(lines)


def write_overlay_gnuplot(sim_dat_path, exp_dat_path, gplot_path, sim_xy, exp_xy,
                           best_name, similarity, two_theta_range):
    """Writes the best candidate's simulated pattern vs. the experimental
    one as two small (2theta, normalized intensity) .dat files (different,
    generally unequal grids -- pyxtal's own broadened-profile resolution
    for one, the experimental file's own sampling for the other -- so two
    separate data files plotted on the same axes, not one shared-grid file)
    plus one .gplot overlay script.
    """
    lo, hi = two_theta_range

    def write_xy(path, x, y, label):
        with open(path, "w") as f:
            f.write(f"# {label} pattern, normalized intensity (best match: {best_name})\n")
            f.write("# 2theta(deg)  intensity(normalized)\n")
            for xv, yv in zip(x, y):
                f.write(f"{xv:12.6f} {yv:12.6f}\n")

    write_xy(sim_dat_path, sim_xy[0], sim_xy[1], "Simulated")
    write_xy(exp_dat_path, exp_xy[0], exp_xy[1], "Experimental")

    sim_b, exp_b = os.path.basename(sim_dat_path), os.path.basename(exp_dat_path)
    lines = [
        'set terminal pdfcairo enhanced font "Arial,14" size 10,5\n',
        'set output "xrd_rank_overlay.pdf"\n',
        f'set title "Best match: {best_name} (similarity {similarity:.4f})"\n',
        'set xlabel "2{/Symbol Q} (deg)" font "Arial,16"\n',
        'set ylabel "Normalized intensity" font "Arial,16"\n',
        f'set xrange [{lo}:{hi}]\n',
        'set grid xtics ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
        '\n',
        f'plot "{sim_b}" using 1:2 with lines lw 2 lc rgb "navy" title "Simulated", \\\n'
        f'     "{exp_b}" using 1:2 with lines lw 2 lc rgb "orange" title "Experimental"\n',
    ]
    with open(gplot_path, "w") as f:
        f.writelines(lines)


def plot_matplotlib(shown, sim_xy, exp_xy, best_name, similarity, two_theta_range):
    """Interactive two-panel preview: the ranking as a horizontal bar chart
    (best match highlighted, on top) and the best candidate's simulated vs.
    experimental pattern overlay below it -- one figure instead of two
    separate windows, so both pieces of information (who wins, and does the
    winner's pattern actually look right) are visible at a glance. Blocking
    plt.show(), same convention already used by stb-xrd/stb-density/
    stb-workfunction's own --view.
    """
    lo, hi = two_theta_range
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))

    names = [r[0] for r in shown][::-1]
    scores = [r[1] for r in shown][::-1]
    colors = ['C2' if i == len(shown) - 1 else 'C0' for i in range(len(shown))]
    ax1.barh(names, scores, color=colors)
    ax1.set_xlabel("Similarity")
    ax1.set_xlim(0, 1)
    ax1.set_title(f"Candidate ranking ({len(shown)} shown, best match highlighted)")
    ax1.grid(True, axis='x', alpha=0.3)

    ax2.plot(sim_xy[0], sim_xy[1], label="Simulated", color="C0")
    ax2.plot(exp_xy[0], exp_xy[1], label="Experimental", color="C1", alpha=0.7)
    ax2.set_xlabel("2-theta (deg)")
    ax2.set_ylabel("Normalized intensity")
    ax2.set_xlim(lo, hi)
    ax2.set_title(f"Best match: {best_name} (similarity {similarity:.4f})")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    plt.show()


def main():
    library_warnings = []

    # require_pyxtal() prints its own [CRITICAL ERROR]+install hint and exits if missing --
    # must run OUTSIDE capture_library_noise, or that message would be captured into a
    # buffer that's never flushed (sys.exit() skips capture_library_noise's post-yield
    # flush step), leaving the user with no output at all instead of the install hint.
    require_pyxtal()
    with capture_library_noise(library_warnings, "pyxtal import"):
        from pyxtal.XRD import Similarity

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Ranks candidate structures by similarity to an experimental XRD pattern (Stage 2 - Analysis).", 'bold')}
Part of the Structure Solution (XRD) workflow: reads every candidate structure
in --input-dir (default 'xrd_search', stb-xrdsearch's own default output
folder -- one candidate per subfolder, picking up a relaxed '*.STRUCT_OUT'
over the raw 'structure.fdf' if you've already run SIESTA there), simulates
each one's powder XRD pattern, and ranks them by similarity to an
experimental pattern -- a fast pre-screen for which candidate space group/
arrangement is most likely correct. Run it once on the raw candidates to
prioritize which to relax with real DFT, then again on the relaxed results
to confirm. Each comparison takes 10s of seconds (pyxtal's Similarity() has
no faster mode), so this can take a while for a large --input-dir.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --experimental measured.dat\n"
               "  %(prog)s --input-dir candidates --experimental measured.dat \\\n"
               "      --wavelength MoKa --top 5\n"
               "  %(prog)s --experimental measured.dat --save-gnuplot --view --save-report\n"
    )

    parser.add_argument("--input-dir", type=str, default="xrd_search",
                        help="Folder of candidate structures, one per subfolder (e.g. written "
                             "by stb-xrdsearch as '<subfolder>/structure.fdf') -- defaults to "
                             "'xrd_search', stb-xrdsearch's own default --output-dir, so Stage 2 "
                             "runs against Stage 1's output with no extra flag in the common "
                             "case. If a subfolder also has a '*.STRUCT_OUT' file (from running "
                             "SIESTA there), the relaxed structure is used instead of the raw "
                             "one. Flat .fdf files directly in --input-dir are also accepted.")
    parser.add_argument("--experimental", type=str, required=True, metavar="PATH",
                        help="Experimental XRD pattern: a plain text file with two columns "
                             "(2theta intensity), whitespace- or comma-separated, '#' comments "
                             "and blank lines allowed. If the file looks like a sparse list of "
                             "discrete peaks (e.g. a literature peak table, or stb-xrd's own "
                             "stick-pattern output) rather than a continuous scan, it is "
                             "automatically Gaussian-broadened before comparing -- see "
                             "--raw-experimental.")
    parser.add_argument("--raw-experimental", action="store_true",
                        help="Never auto-broaden a peak-list-looking --experimental file -- "
                             "compare it exactly as read. Off by default (auto-detection "
                             "normally handles this correctly; matches stb-xrd's own "
                             "--compare-to flag of the same name).")
    parser.add_argument("--wavelength", type=str, default="CuKa",
                        help="X-ray source: a known name (CuKa, CuKa1, MoKa, CrKa, FeKa, CoKa, "
                             "AgKa, ...) or a wavelength in Ang as a plain number "
                             "(default: CuKa, 1.54184 Ang).")
    parser.add_argument("--two-theta-range", type=float, nargs=2, default=[0.0, 90.0],
                        metavar=("MIN", "MAX"),
                        help="2-theta range in degrees to compute (default: 0 90).")
    parser.add_argument("--top", type=int, default=None,
                        help="Only show/plot the N best-matching candidates (default: show "
                             "all). The full ranking is always written to --output regardless.")
    parser.add_argument("-o", "--output", type=str, default="xrd_rank.txt",
                        help="Output ranking table file name (default: xrd_rank.txt). Always "
                             "lists every candidate, regardless of --top. Gnuplot/report files "
                             "(--save-gnuplot/--save-report) are written alongside it, in the "
                             "same directory.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write the ranking as a gnuplot bar chart (xrd_rank.dat + "
                             "xrd_rank.gplot) and the best candidate's simulated-vs-experimental "
                             "pattern overlay (xrd_rank_top_sim.dat/xrd_rank_top_exp.dat + "
                             "xrd_rank_overlay.gplot), all under a 'plot/' subfolder alongside "
                             "--output. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview (ranking bar chart + best "
                             "candidate's simulated-vs-experimental overlay, one figure) before "
                             "finishing. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full numbered run report to {REPORT_FILE} "
                             "(written alongside --output). Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-xrdrank {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    lo, hi = args.two_theta_range
    if lo >= hi:
        parser.error("--two-theta-range MIN must be less than MAX.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if not os.path.isdir(args.input_dir):
        print(color_text(
            f"Error: '{args.input_dir}' is not a directory -- run stb-xrdsearch first "
            "(default output folder: 'xrd_search'), or point --input-dir at your own "
            "candidate folder.", 'red'))
        sys.exit(1)

    output_dir = os.path.dirname(args.output) or "."
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-XRDRANK STAGE 2 REPORT (STRUCTURE SOLUTION RANKING) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time              : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input directory        : {args.input_dir}", f_out)
    print_dual(f"Experimental pattern   : {args.experimental}", f_out)

    try:
        wavelength, wavelength_label = xrd_core.resolve_wavelength(args.wavelength)
    except ValueError as e:
        fail(str(e))
    print_dual(f"Wavelength             : {wavelength_label} ({wavelength:.5f} Ang)", f_out)
    print_dual(f"2-theta range          : {lo} - {hi} deg", f_out)
    print_dual(f"Top shown/plotted      : {'all' if args.top is None else args.top}", f_out)
    print_dual(f"Save gnuplot           : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib)      : {'yes' if args.view else 'no'}", f_out)
    print_dual(f"Report                 : {report_path if report_path else '(not saved)'}", f_out)

    print_section("[1] CANDIDATES FOUND", f_out)
    candidate_files = find_candidates(args.input_dir)
    if not candidate_files:
        fail(f"No candidate structures found in '{args.input_dir}' (expected "
             "'<subfolder>/structure.fdf', '<subfolder>/*.STRUCT_OUT', or flat .fdf files).")
    n_relaxed = sum(1 for *_, source in candidate_files if source == "relaxed")
    print_dual(f"Found {len(candidate_files)} candidate(s) "
               f"({n_relaxed} relaxed, {len(candidate_files) - n_relaxed} raw).", f_out)
    print_table(["Name", "Source", "File"], [
        ([name, source, path], 'green' if source == 'relaxed' else None)
        for name, path, fmt, source in candidate_files
    ], f_out)

    print_section("[2] EXPERIMENTAL PATTERN", f_out)
    try:
        experimental = xrd_core.read_experimental_pattern(args.experimental)
    except (FileNotFoundError, ValueError) as e:
        fail(str(e))
    is_peak_list = xrd_core.looks_like_peak_list(experimental[0])
    broadened = False
    if is_peak_list and not args.raw_experimental:
        with capture_library_noise(library_warnings, "pyxtal broaden_peak_list"):
            experimental = xrd_core.broaden_peak_list(experimental[0], experimental[1], (lo, hi))
        broadened = True
    print_dual(f"Points                 : {experimental.shape[1]}", f_out)
    print_dual(f"2-theta span (file)    : {experimental[0].min():.3f} - "
               f"{experimental[0].max():.3f} deg", f_out)
    if broadened:
        print_dual(color_text(
            "Peak-list detected     : yes -- auto-broadened (Gaussian, FWHM=0.1 deg) before "
            "comparing, so it's on the same footing as each candidate's own broadened "
            "simulated profile (pass --raw-experimental to disable).", 'cyan'), f_out)
    elif is_peak_list:
        print_dual(color_text(
            "Peak-list detected     : yes -- but --raw-experimental was given, so it is "
            "compared unbroadened (low similarity scores may just reflect that mismatch, "
            "not a real structural difference).", 'yellow'), f_out)
    else:
        print_dual("Peak-list detected     : no -- compared as a continuous scan, unbroadened.", f_out)

    print_section("[3] RANKING", f_out)
    results = []
    for name, path, fmt, source in candidate_files:
        # capture_library_noise wraps only the third-party calls (structure read, pattern
        # computation, Similarity) -- the print_dual status line below stays outside it, so
        # this tool's own report output is never redirected/delayed. Called once per
        # candidate, so the same underlying warning can be appended once per candidate that
        # triggers it; the final [6] LIBRARY WARNINGS section dedups by exact text before
        # printing, so this stays clean regardless of how many candidates were scanned.
        try:
            with capture_library_noise(library_warnings, "pyxtal XRD ranking"):
                pmg = structure_io.read_siesta_structure(path, fmt)
                pattern = xrd_core.compute_pattern(pmg, wavelength=wavelength,
                                                    two_theta_range=(lo, hi))
                sim_profile = pattern.get_profile()
                similarity = run_with_spinner(
                    lambda sp=sim_profile: Similarity(sp, experimental).value,
                    label=f"Comparing {name}")
                # space_group_label (pymatgen -> spglib) is called INSIDE the same capture
                # block -- found live (against this exact real Si data) that spglib emits its
                # own 'Set OLD_ERROR_HANDLING...' DeprecationWarning here, which leaked
                # straight to the terminal when this call sat outside any capture_library_
                # noise wrap (Python's own once-per-location warning dedup made it invisible
                # in some runs and a raw leak in others, depending on whether pyxtal's own
                # internal spglib use had already "used up" that one-time default warning
                # earlier in the same process -- wrapping it removes that dependency on
                # ordering entirely). A generic Exception here (unrelated to the comparison
                # itself) is caught locally so a symmetry-detection hiccup never discards an
                # otherwise-valid similarity score.
                try:
                    sg = space_group_label(pmg)
                except Exception:
                    sg = "N/A"
        except (FileNotFoundError, ValueError) as e:
            print_dual(color_text(f"  Warning: skipping '{name}' ({e}).", 'yellow'), f_out)
            continue
        results.append((name, similarity, pmg.composition.reduced_formula, source, sg))
        print_dual(f"  {color_text('->', 'green')} {name}: similarity {similarity:.4f} "
                   f"({pmg.composition.reduced_formula}, {sg}, {source})", f_out)

    if not results:
        fail("None of the candidates could be compared -- see warnings above.")

    results.sort(key=lambda r: r[1], reverse=True)
    shown = results if args.top is None else results[:args.top]

    print_table(["Rank", "Name", "Similarity", "Formula", "Space group", "Source"], [
        ([str(rank), name, f"{similarity:.4f}", formula, sg, source],
         'green' if rank == 1 else None)
        for rank, (name, similarity, formula, source, sg) in enumerate(shown, start=1)
    ], f_out)
    if len(results) >= 2:
        gap = results[0][1] - results[1][1]
        print_dual(f"Score gap (#1 vs #2)   : {gap:.4f}" +
                   ("  [close call -- consider relaxing both before deciding]"
                    if gap < 0.02 else ""), f_out)

    best_name, best_similarity, best_formula, best_source, best_sg = results[0]

    print_section("[4] OUTPUT DATA & PLOTS", f_out)
    sim_xy = exp_xy = None
    if args.save_gnuplot or args.view:
        # Recomputed for the best candidate only -- compute_pattern() itself is fast (the
        # slow step, Similarity(), already ran once per candidate above and isn't repeated).
        best_path, best_fmt = next((path, fmt) for name, path, fmt, _ in candidate_files
                                    if name == best_name)
        with capture_library_noise(library_warnings, "pyxtal overlay recompute"):
            best_pmg = structure_io.read_siesta_structure(best_path, best_fmt)
            best_pattern = xrd_core.compute_pattern(best_pmg, wavelength=wavelength,
                                                      two_theta_range=(lo, hi))
            sim_x, sim_y = best_pattern.get_profile()
        sim_y = sim_y / sim_y.max() if sim_y.max() > 0 else sim_y
        exp_x, exp_y = experimental
        exp_y = exp_y / exp_y.max() if exp_y.max() > 0 else exp_y
        sim_xy, exp_xy = (sim_x, sim_y), (exp_x, exp_y)

    if args.save_gnuplot:
        # All plot data/scripts (and the PDFs gnuplot itself produces from them, via each
        # .gplot's own relative "set output") live together under their own 'plot/'
        # subfolder, kept separate from the plain-text ranking/report files right beside it.
        plot_dir = os.path.join(output_dir, "plot")
        os.makedirs(plot_dir, exist_ok=True)

        rank_dat = os.path.join(plot_dir, "xrd_rank.dat")
        rank_gplot = os.path.join(plot_dir, "xrd_rank.gplot")
        write_rank_gnuplot(rank_dat, rank_gplot, shown, wavelength_label)
        print_dual(color_text(f"[OK] Ranking bar chart written to '{rank_dat}' / '{rank_gplot}'.",
                               'green'), f_out)

        sim_dat = os.path.join(plot_dir, "xrd_rank_top_sim.dat")
        exp_dat = os.path.join(plot_dir, "xrd_rank_top_exp.dat")
        overlay_gplot = os.path.join(plot_dir, "xrd_rank_overlay.gplot")
        write_overlay_gnuplot(sim_dat, exp_dat, overlay_gplot, sim_xy, exp_xy,
                               best_name, best_similarity, (lo, hi))
        print_dual(color_text(
            f"[OK] Best-match overlay written to '{sim_dat}' / '{exp_dat}' / '{overlay_gplot}'.",
            'green'), f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write the ranking "
                   "bar chart and best-match overlay as .dat/.gplot pairs, under a 'plot/' "
                   "subfolder).", f_out)

    print_section("[5] SUMMARY & NEXT STEPS", f_out)
    with open(args.output, "w") as f:
        name_width = max(28, max((len(r[0]) for r in results), default=0) + 2)
        header = f"{'Rank':<5} {'File':<{name_width}}{'Similarity':<12}{'Formula':<15}{'Space group':<16}{'Source':<10}"
        f.write(header + "\n")
        for rank, (name, similarity, formula, source, sg) in enumerate(results, start=1):
            f.write(f"{rank:<5} {name:<{name_width}}{similarity:<12.4f}{formula:<15}{sg:<16}{source:<10}\n")
    print_dual(color_text(f"[OK] Full ranking ({len(results)} candidate(s)) written to "
                          f"'{args.output}'.", 'green'), f_out)
    print_dual(f"Best match             : {best_name} (similarity {best_similarity:.4f}, "
               f"{best_formula}, {best_sg})", f_out)
    print_dual(color_text(
        "\nNote: a similarity ranking from simulated vs. experimental XRD, not a full "
        "Rietveld refinement -- use it to prioritize which candidate(s) to investigate "
        "further (e.g. relax with stb-mlrelax, then real DFT), not as definitive proof.",
        'yellow'), f_out)
    print_dual(f"\n{color_text('Next:', 'cyan')} relax the top candidate(s) with real SIESTA "
              "(or stb-mlrelax first for a quick pre-screen), then rerun stb-xrdrank on the "
              "relaxed .STRUCT_OUT results to confirm.", f_out)

    print_section("[6] LIBRARY WARNINGS", f_out)
    library_warnings = list(dict.fromkeys(library_warnings))  # exact-text dedup, order preserved
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (pyxtal/scipy/numpy/pymatgen/spglib) "
            "during this run -- collected here instead of interleaved with the report above; "
            "harmless in almost every case (import-time notices, deprecation-style warnings), "
            "but worth a glance.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a blocking matplotlib
    # window never delays or hides it (same convention as stb-xrd's own --view).
    if args.view:
        plot_matplotlib(shown, sim_xy, exp_xy, best_name, best_similarity, (lo, hi))


if __name__ == "__main__":
    main()
