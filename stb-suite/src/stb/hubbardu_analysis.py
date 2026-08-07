#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

# 2.0.0: numbered-section report (matching stages 1/2's convention); --save-report;
# gnuplot .dat+.gplot output moved under its own 'plot/' subfolder and made opt-in
# via --save-gnuplot (previously written unconditionally, directly in the cwd);
# new --view matplotlib preview of the same response curve; a per-run READING RUNS
# table instead of one print() line per run; a final LIBRARY WARNINGS section
# (numpy's own RankWarning on a poorly-conditioned polyfit is a real, relevant
# signal for this tool's own core fit -- now captured there instead of leaking).
VERSION = "2.0.0"

import os
import re
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from stb.core.cli import (COLORS, color_text, show_intro, print_dual, print_section,
                          print_table, capture_library_noise)
from stb.core.dftu_data import REFERENCE_U, ldau_proj_block, load_manifest

REPORT_FILE = "hubbardu_stage3.txt"

_OCC_RE = re.compile(r'Occupations:\s*([0-9eE+.\-\s]+)')


def parse_occupation(out_path):
    """Last 'Occupations:' line's TOTAL shell occupation (summed over spin)
    from a SIESTA .out file produced with LDAU.PotentialShift T.

    The line has either 2 or 3 numbers depending on nspin (verified against
    SIESTA's own source, Src/dftu.F, format
    write(6,'(a,/,a,3f12.6)') 'hubbard_term: Total projector shell',
    'Occupations: ', (oc(ispin),ispin=1,nspin), sum(oc)
    -- i.e. it always prints nspin per-spin values followed by sum(oc)):

    - nspin=2 (spin-polarized): 3 numbers (up, down, sum) -- the last one
      IS the genuine total, no adjustment needed.
    - nspin=1 (non-polarized): only 2 numbers, and BOTH are oc(1) --
      sum(oc) degenerates to oc(1) itself when there's only one spin
      channel to sum. Critically, occu() (and so oc(1)) is computed
      upstream as Dij*Sik*Sjk/(3-nspin) -- i.e. divided by 2 for nspin=1,
      because SIESTA's non-polarized Dscf already stores the FULL
      (both-spin) density in its single channel. So oc(1) here is only
      the PER-SPIN occupation, and the true total electron count in the
      shell is 2x that -- verified directly: for a non-polarized Mn run,
      the printed value was ~2.47 (implausibly low for a d-shell), and
      halving worked backwards from the /(3-nspin) factor confirms 2x
      that (~4.94) is the physically-expected total. Not doubling this
      silently halved every non-polarized run's occupation and slope,
      doubling the computed U.

    Returns None on any parse error.
    """
    last = None
    try:
        with open(out_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                match = _OCC_RE.search(line)
                if match:
                    nums = match.group(1).split()
                    if nums:
                        try:
                            total = float(nums[-1])
                            if len(nums) == 2:
                                total *= 2.0
                            last = total
                        except ValueError:
                            pass
    except Exception:
        return None
    return last


def check_scf_converged(out_path):
    """True if the LAST convergence marker in the file is 'SCF Convergence
    by ...', False if it's 'SCF_NOT_CONV', None if neither marker is found.
    Uses whichever marker occurs LAST in the file, not just "is either
    present anywhere" -- a file that's the concatenation of a failed attempt
    followed by a successful restart (e.g. SIESTA re-run with `>>` instead of
    overwriting calc.out) has both markers, and only the last one reflects
    the run's actual final state. A not-converged run can still have a
    usable occupation value (SIESTA prints it before aborting), but the
    caller should be able to tell converged points from not-quite-converged
    ones instead of only finding out indirectly via a poor R^2.
    """
    result = None
    try:
        with open(out_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "SCF_NOT_CONV" in line:
                    result = False
                elif "SCF Convergence by" in line:
                    result = True
    except Exception:
        return None
    return result


def linear_fit_r2(alphas, occs):
    """np.polyfit degree-1 fit, returning (slope, intercept, r_squared)."""
    alphas = np.array(alphas)
    occs = np.array(occs)
    slope, intercept = np.polyfit(alphas, occs, 1)
    predicted = slope * alphas + intercept
    ss_res = np.sum((occs - predicted) ** 2)
    ss_tot = np.sum((occs - np.mean(occs)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(intercept), float(r2)


# How far the fitted alpha=0 intercept may drift from the actually-measured
# alpha=0 occupation before it's flagged -- a decent overall R^2 can still
# hide a systematically-off zero point (e.g. the fit tilting to average out
# an intercept mismatch rather than a genuinely noisy slope).
INTERCEPT_DEVIATION_TOL = 0.01


def write_response_plot(data_path, gplot_path, scf_points, frozen_points,
                        chi, chi_intercept, chi0, chi0_intercept, species):
    """Writes a gnuplot data+script pair plotting occupation vs. perturbation
    strength for both branches, with the fitted lines overlaid -- the same
    data-file + .gplot convention used throughout the suite (stb-kgrid,
    stb-density, stb-workfunction, stb-xrd). `data_path`/`gplot_path` are
    full paths (typically under a 'plot/' subfolder); the script itself
    references the data file by BASENAME only (matching stb-xrdrank's own
    plot/ convention), so it works as-is when run from inside that folder
    regardless of where the caller's own cwd was.
    """
    with open(data_path, 'w') as f:
        f.write("# index 0: self-consistent (screened) response\n")
        f.write("# alpha(eV) occupation\n")
        for alpha, occ in scf_points:
            f.write(f"{alpha:.6f}  {occ:.6f}\n")
        f.write("\n\n")
        f.write("# index 1: frozen-density (bare) response\n")
        f.write("# alpha(eV) occupation\n")
        for alpha, occ in frozen_points:
            f.write(f"{alpha:.6f}  {occ:.6f}\n")

    data_basename = os.path.basename(data_path)
    gplot_basename = os.path.basename(gplot_path)
    lines = []
    lines.append('# --- STB Plot Configuration ---\n')
    lines.append('# Generated by stb-hubbarduAnalysis\n')
    lines.append('set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n')
    lines.append(f'set output "{gplot_basename.rsplit(".", 1)[0]}.pdf"\n\n')
    lines.append(f'set title "Hubbard U linear response ({species})"\n')
    lines.append('set xlabel "Perturbation alpha (eV)"\n')
    lines.append('set ylabel "Shell occupation"\n')
    lines.append('set grid\n')
    lines.append('set key top left\n')
    lines.append(f'f_scf(x) = {chi:.8f}*x + {chi_intercept:.8f}\n')
    lines.append(f'f_frozen(x) = {chi0:.8f}*x + {chi0_intercept:.8f}\n')
    lines.append(
        f'plot "{data_basename}" index 0 using 1:2 with points pt 7 ps 1.5 lc rgb "#cc5522" '
        'title "Self-consistent (screened)", \\\n'
        '     f_scf(x) with lines lc rgb "#cc5522" dt 2 notitle, \\\n'
        f'     "{data_basename}" index 1 using 1:2 with points pt 7 ps 1.5 lc rgb "#2255cc" '
        'title "Frozen-density (bare)", \\\n'
        '     f_frozen(x) with lines lc rgb "#2255cc" dt 2 notitle\n'
    )
    with open(gplot_path, 'w') as f:
        f.writelines(lines)


def plot_matplotlib(scf_points, frozen_points, chi, chi_intercept, chi0, chi0_intercept,
                     species, U):
    """Interactive preview of the same response curve write_response_plot()
    writes as a gnuplot pair: both branches' points plus their fitted lines,
    on one figure. Blocking plt.show(), same convention as stb-xrd/
    stb-density/stb-workfunction/stb-xrdrank's own --view.
    """
    fig, ax = plt.subplots(figsize=(9, 6))

    scf_x = [a for a, _ in scf_points]
    scf_y = [o for _, o in scf_points]
    frozen_x = [a for a, _ in frozen_points]
    frozen_y = [o for _, o in frozen_points]

    x_span = [min(scf_x + frozen_x), max(scf_x + frozen_x)]
    ax.scatter(scf_x, scf_y, color="#cc5522", zorder=3, label="Self-consistent (screened)")
    ax.plot(x_span, [chi * x + chi_intercept for x in x_span],
            color="#cc5522", linestyle="--", alpha=0.7)
    ax.scatter(frozen_x, frozen_y, color="#2255cc", zorder=3, label="Frozen-density (bare)")
    ax.plot(x_span, [chi0 * x + chi0_intercept for x in x_span],
            color="#2255cc", linestyle="--", alpha=0.7)

    ax.set_xlabel("Perturbation alpha (eV)")
    ax.set_ylabel("Shell occupation")
    ax.set_title(f"Hubbard U linear response ({species}) -- U = {U:.4f} eV")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()


def main():
    library_warnings = []

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 3 of 3: analyzes a stb-hubbardu linear-response sweep and computes the Hubbard U.", 'bold')}
Fits the self-consistent (screened) and frozen-density (bare) occupation
responses to the perturbation strength, then applies the standard
Cococcioni & de Gironcoli formula: U = 1/chi0 - 1/chi (Phys. Rev. B 71,
035105, 2005). Reports a convergence summary, both fits' R^2 and intercepts
(flagging a fit whose zero point drifted from what was actually measured
there), and warns if chi/chi0 have opposite signs or U comes out negative --
all signs of a bad fit rather than a real result. Writes a ready-to-use
%%block LDAU.proj snippet with the computed U, plus (opt-in) a gnuplot
data+script pair and/or a live matplotlib preview to inspect the fit
visually instead of trusting R^2 alone.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --dir hubbardu_runs --file calc.out\n"
               "  %(prog)s --dir hubbardu_runs --save-gnuplot --view --save-report\n"
    )

    parser.add_argument("--dir", type=str, default="hubbardu_runs",
                        help="Directory containing the run folders and run_manifest.json "
                             "written by stb-hubbardu/stb-hubbarduAlphas (stages 1/2) -- "
                             "defaults to 'hubbardu_runs', their own default --output-dir/--dir, "
                             "so stage 3 runs against their output with no extra flag in the "
                             "common case.")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--r2-tolerance", type=float, default=0.98,
                        help="Warn if either fit's R^2 falls below this (default: 0.98) -- a poor "
                             "linear fit means the perturbation range is probably too wide.")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output .fdf snippet filename (default: <species>_LDAU.fdf).")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write the response curve as a gnuplot data+script pair "
                             "(<species>_hubbardu_response.dat/.gplot), under a 'plot/' "
                             "subfolder inside --dir. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of the response curve "
                             "(both branches + fitted lines) before finishing. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbarduAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    # os.path.isdir never expands '~' itself (that's shell-only, and argv
    # passed through subprocess.run() never goes through a shell), so a path
    # like '~/hubbardu_runs' would otherwise be checked completely literally.
    args.dir = os.path.expanduser(args.dir)
    if args.output:
        args.output = os.path.expanduser(args.output)

    if not os.path.isdir(args.dir):
        sys.exit(color_text(f"Error: Directory '{args.dir}' not found.", 'red'))

    manifest_path = os.path.join(args.dir, "run_manifest.json")
    if not os.path.isfile(manifest_path):
        sys.exit(color_text(
            f"Error: 'run_manifest.json' not found in '{args.dir}' -- this directory "
            "wasn't generated by stb-hubbardu.", 'red'))

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = os.path.join(args.dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-HUBBARDUANALYSIS STAGE 3 REPORT (HUBBARD U COMPUTATION) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Directory              : {args.dir}", f_out)
    print_dual(f"SIESTA output filename : {args.file}", f_out)
    print_dual(f"R^2 tolerance          : {args.r2_tolerance}", f_out)
    print_dual(f"Save gnuplot           : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib)      : {'yes' if args.view else 'no'}", f_out)
    print_dual(f"Report                 : {report_path if report_path else '(not saved)'}", f_out)

    try:
        manifest = load_manifest(manifest_path)
    except (ValueError, json.JSONDecodeError) as e:
        fail(str(e))
    species, n, l, j = manifest["species"], manifest["n"], manifest["l"], manifest["j"]
    print_dual(f"Species / shell        : {species} (n={n}, l={l}), J={j:.3f} eV", f_out)

    print_section("[1] READING RUNS", f_out)
    scf_points, frozen_points = [], []
    frozen_iterations_seen = set()
    n_skipped = 0
    n_converged = n_not_converged = n_unknown = 0
    n_frozen_unknown = 0
    n_manifest_runs = len(manifest["runs"])
    run_rows = []
    for folder, info in manifest["runs"].items():
        if "kind" not in info or "alpha" not in info:
            run_rows.append(([folder, "?", "?", "--", "SKIP (malformed manifest entry)"], 'yellow'))
            n_skipped += 1
            continue

        out_path = os.path.join(args.dir, folder, args.file)
        occ = parse_occupation(out_path)
        if occ is None:
            run_rows.append(([folder, info["kind"], f"{info['alpha']:.4f}", "--",
                              f"SKIP (missing/unparseable {args.file})"], 'yellow'))
            n_skipped += 1
            continue

        kind, alpha = info["kind"], info["alpha"]
        if kind in ("scf", "reference"):
            scf_points.append((alpha, occ))
        if kind == "frozen":
            frozen_points.append((alpha, occ))
            if "frozen_iterations" in info:
                frozen_iterations_seen.add(info["frozen_iterations"])

        converged = check_scf_converged(out_path)
        row_color = None
        if kind == "frozen":
            # SCF_NOT_CONV is the EXPECTED outcome here, not a problem: the
            # frozen branch deliberately caps at MaxSCFIterations (1 by
            # default) to measure the bare/unrelaxed response, before the
            # density has a chance to reach real self-consistency.
            if converged is None:
                n_frozen_unknown += 1
                conv_str, row_color = "status unknown", 'yellow'
            elif converged is False:
                conv_str = "frozen eval. complete (stopped at MaxSCFIterations, as designed)"
            else:
                conv_str = "frozen eval. converged before hitting MaxSCFIterations"
        else:
            if converged is True:
                n_converged += 1
                conv_str = "converged"
            elif converged is False:
                n_not_converged += 1
                conv_str, row_color = "NOT CONVERGED, used anyway", 'yellow'
            else:
                n_unknown += 1
                conv_str, row_color = "status unknown", 'yellow'
        run_rows.append(([folder, kind, f"{alpha:.4f}", f"{occ:.6f}", conv_str], row_color))

    print_table(["Folder", "Kind", "Alpha (eV)", "Occupation", "SCF status"], run_rows, f_out)

    n_read = n_manifest_runs - n_skipped
    n_scf_read = len(scf_points)
    print_dual(f"\nConvergence summary    : {n_read}/{n_manifest_runs} run(s) read"
              + (f" ({n_skipped} skipped)" if n_skipped else "") + f"; of the {n_scf_read} "
              f"self-consistent (scf/reference) run(s), {n_converged}/{n_scf_read} fully converged"
              + (f", {n_not_converged} NOT converged, {n_unknown} unknown"
                 if n_not_converged or n_unknown else "")
              + (f". Frozen-density runs: {len(frozen_points) - n_frozen_unknown}/{len(frozen_points)} "
                 "completed as expected" + (f", {n_frozen_unknown} unknown" if n_frozen_unknown else "")
                 if frozen_points else ""), f_out)

    if len(frozen_iterations_seen) > 1:
        print_dual(color_text(
            f"[WARNING] The frozen-density runs used different --frozen-iterations values across "
            f"this sweep ({sorted(frozen_iterations_seen)}) -- the frozen branch is only physically "
            "meaningful when every point uses the exact same recipe. Regenerate them with a "
            "consistent --frozen-iterations in stb-hubbarduAlphas.", 'yellow'), f_out)

    if len(scf_points) < 2 or len(frozen_points) < 2:
        fail("Not enough valid runs to fit a response (need at least 2 scf points, including "
             "the reference, and 2 frozen points).")

    for name, points in (("self-consistent", scf_points), ("frozen-density", frozen_points)):
        if len(points) == 2:
            print_dual(color_text(
                f"[WARNING] The {name} branch has only 2 points -- a line always fits 2 points "
                "exactly, so its R^2 will read 1.0 regardless of whether the response is actually "
                "linear. Add more --alphas in stb-hubbarduAlphas for a meaningful check.", 'yellow'), f_out)

    print_section("[2] LINEAR RESPONSE FIT", f_out)
    scf_points.sort()
    frozen_points.sort()
    with capture_library_noise(library_warnings, "numpy polyfit"):
        chi, chi_intercept, chi_r2 = linear_fit_r2(*zip(*scf_points))
        chi0, chi0_intercept, chi0_r2 = linear_fit_r2(*zip(*frozen_points))

    print_table(["Branch", "chi (1/eV)", "Intercept", "R^2"], [
        (["Self-consistent (screened)", f"{chi:.6f}", f"{chi_intercept:.6f}", f"{chi_r2:.4f}"],
         'yellow' if chi_r2 < args.r2_tolerance else None),
        (["Frozen-density (bare)", f"{chi0:.6f}", f"{chi0_intercept:.6f}", f"{chi0_r2:.4f}"],
         'yellow' if chi0_r2 < args.r2_tolerance else None),
    ], f_out)

    for name, r2 in (("self-consistent", chi_r2), ("frozen-density", chi0_r2)):
        if r2 < args.r2_tolerance:
            print_dual(color_text(
                f"[WARNING] The {name} response fit has R^2={r2:.4f}, below the "
                f"{args.r2_tolerance} tolerance -- the perturbation range may be too wide "
                "for the linear regime. Consider narrower --alphas in stb-hubbarduAlphas.", 'yellow'), f_out)

    # A decent overall R^2 can still hide a fit whose zero point drifted from
    # what was actually measured there (the fit tilting to average it out
    # rather than a genuinely noisy slope) -- compare the fitted intercept
    # against the actual alpha=0 occupation directly.
    scf_zero = next((occ for alpha, occ in scf_points if alpha == 0.0), None)
    frozen_zero = next((occ for alpha, occ in frozen_points if alpha == 0.0), None)
    for name, zero_occ, intercept in (("self-consistent", scf_zero, chi_intercept),
                                       ("frozen-density", frozen_zero, chi0_intercept)):
        if zero_occ is not None and abs(zero_occ - intercept) > INTERCEPT_DEVIATION_TOL:
            print_dual(color_text(
                f"[WARNING] The {name} fit's intercept ({intercept:.6f}) differs from the "
                f"actually-measured alpha=0 occupation ({zero_occ:.6f}) by more than "
                f"{INTERCEPT_DEVIATION_TOL} -- the response may not be as linear as the R^2 "
                "alone suggests.", 'yellow'), f_out)

    if chi == 0 or chi0 == 0:
        fail("A zero response slope was fitted -- cannot compute U (division by zero). Check "
             "that the occupation actually responds to the perturbation.")

    print_section("[3] COMPUTED U", f_out)
    if (chi > 0) != (chi0 > 0):
        print_dual(color_text(
            "[WARNING] chi and chi0 have OPPOSITE signs -- physically, both describe the same "
            "direction of response (just screened vs. bare), so this usually means a bad run, "
            "a mislabeled point, or a fit dominated by noise. Treat the computed U with caution.",
            'yellow'), f_out)

    U = 1.0 / chi0 - 1.0 / chi
    if U < 0:
        print_dual(color_text(
            "[WARNING] Computed U is NEGATIVE -- the Hubbard U is a repulsive on-site term and "
            "should always be positive. This strongly suggests a bad fit or bad data rather "
            "than a real result.", 'yellow'), f_out)
    print_dual(f"Formula                : U = 1/chi0 - 1/chi (Cococcioni & de Gironcoli, "
              "PRB 71, 035105, 2005)", f_out)
    print_dual(color_text(f"Computed U             : {U:.4f} eV", 'green'), f_out)

    if species in REFERENCE_U:
        print_dual(color_text(
            f"[INFO] Literature reference (GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 "
            f"(2006)) for {species}: {REFERENCE_U[species]:.2f} eV -- a sanity-check comparison "
            "only, not a validation (functional/pseudopotential/basis-dependent).", 'cyan'), f_out)

    print_section("[4] OUTPUT FILES", f_out)
    output = args.output or f"{species}_LDAU.fdf"
    block = ldau_proj_block([{"species": species, "n": n, "l": l, "u": U, "j": j}])
    with open(output, 'w') as f:
        f.write(block)
    print_dual(color_text(f"[OK] Ready-to-use DFT+U block written to '{output}'.", 'green'), f_out)

    if args.save_gnuplot:
        plot_dir = os.path.join(args.dir, "plot")
        os.makedirs(plot_dir, exist_ok=True)
        data_path = os.path.join(plot_dir, f"{species}_hubbardu_response.dat")
        gplot_path = os.path.join(plot_dir, f"{species}_hubbardu_response.gplot")
        write_response_plot(data_path, gplot_path, scf_points, frozen_points,
                            chi, chi_intercept, chi0, chi0_intercept, species)
        print_dual(color_text(
            f"[OK] Response plot written to '{data_path}' / '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Gnuplot data+script : not written (off by default -- pass --save-gnuplot "
                   "to write it under a 'plot/' subfolder).", f_out)

    print_section("[5] SUMMARY & NEXT STEPS", f_out)
    print_dual(color_text(f"Computed U for {species}: {U:.4f} eV", 'green'), f_out)
    print_dual(f"DFT+U block            : {output}", f_out)
    if report_path:
        print_dual(f"Report                 : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Review [2]/[3] above for any [WARNING] before trusting this U value.", f_out)
    print_dual(f"  2. Copy the %block LDAU.proj from '{output}' into your production calc.fdf.", f_out)

    print_section("[6] LIBRARY WARNINGS", f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (numpy) during this run -- collected here "
            "instead of interleaved with the report above; harmless in almost every case, but a "
            "RankWarning specifically means the polyfit itself was poorly conditioned and the "
            "fitted chi/chi0 above may not be trustworthy.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a blocking matplotlib
    # window never delays or hides it (same convention as stb-xrd/stb-xrdrank's own --view).
    if args.view:
        plot_matplotlib(scf_points, frozen_points, chi, chi_intercept, chi0, chi0_intercept,
                        species, U)


if __name__ == "__main__":
    main()
