#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.3.0"

import os
import re
import sys
import json
import argparse
import numpy as np
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.dftu_data import REFERENCE_U, ldau_proj_block, load_manifest

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


def write_response_plot(data_filename, gplot_filename, scf_points, frozen_points,
                        chi, chi_intercept, chi0, chi0_intercept, species):
    """Writes a gnuplot data+script pair plotting occupation vs. perturbation
    strength for both branches, with the fitted lines overlaid -- the same
    data-file + .gplot convention used throughout the suite (stb-kgrid,
    stb-density, stb-workfunction, stb-xrd), so the fit can be inspected
    visually instead of trusting R^2 alone.
    """
    with open(data_filename, 'w') as f:
        f.write("# index 0: self-consistent (screened) response\n")
        f.write("# alpha(eV) occupation\n")
        for alpha, occ in scf_points:
            f.write(f"{alpha:.6f}  {occ:.6f}\n")
        f.write("\n\n")
        f.write("# index 1: frozen-density (bare) response\n")
        f.write("# alpha(eV) occupation\n")
        for alpha, occ in frozen_points:
            f.write(f"{alpha:.6f}  {occ:.6f}\n")

    lines = []
    lines.append('# --- STB Plot Configuration ---\n')
    lines.append('# Generated by stb-hubbarduAnalysis\n')
    lines.append('set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n')
    lines.append(f'set output "{gplot_filename.rsplit(".", 1)[0]}.pdf"\n\n')
    lines.append(f'set title "Hubbard U linear response ({species})"\n')
    lines.append('set xlabel "Perturbation alpha (eV)"\n')
    lines.append('set ylabel "Shell occupation"\n')
    lines.append('set grid\n')
    lines.append('set key top left\n')
    lines.append(f'f_scf(x) = {chi:.8f}*x + {chi_intercept:.8f}\n')
    lines.append(f'f_frozen(x) = {chi0:.8f}*x + {chi0_intercept:.8f}\n')
    lines.append(
        f'plot "{data_filename}" index 0 using 1:2 with points pt 7 ps 1.5 lc rgb "#cc5522" '
        'title "Self-consistent (screened)", \\\n'
        '     f_scf(x) with lines lc rgb "#cc5522" dt 2 notitle, \\\n'
        f'     "{data_filename}" index 1 using 1:2 with points pt 7 ps 1.5 lc rgb "#2255cc" '
        'title "Frozen-density (bare)", \\\n'
        '     f_frozen(x) with lines lc rgb "#2255cc" dt 2 notitle\n'
    )
    with open(gplot_filename, 'w') as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Analyzes a stb-hubbardu linear-response sweep and computes the Hubbard U.", 'bold')}
Fits the self-consistent (screened) and frozen-density (bare) occupation
responses to the perturbation strength, then applies the standard
Cococcioni & de Gironcoli formula: U = 1/chi0 - 1/chi (Phys. Rev. B 71,
035105, 2005). Reports a convergence summary, both fits' R^2 and intercepts
(flagging a fit whose zero point drifted from what was actually measured
there), and warns if chi/chi0 have opposite signs or U comes out negative --
all signs of a bad fit rather than a real result. Writes a ready-to-use
%%block LDAU.proj snippet with the computed U, plus a gnuplot data+script
pair to inspect the fit visually.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --dir hubbardu_runs --file calc.out\n"
    )

    parser.add_argument("--dir", type=str, default="hubbardu_runs",
                        help="Directory containing the run folders and run_manifest.json "
                             "written by stb-hubbardu (default: hubbardu_runs).")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--r2-tolerance", type=float, default=0.98,
                        help="Warn if either fit's R^2 falls below this (default: 0.98) -- a poor "
                             "linear fit means the perturbation range is probably too wide.")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output .fdf snippet filename (default: <species>_LDAU.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbarduAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze a Hubbard U linear-response sweep:", 'bold'))
    print("-" * 60)

    # os.path.isdir never expands '~' itself (that's shell-only, and argv
    # passed through subprocess.run() never goes through a shell), so a path
    # like '~/hubbardu_runs' would otherwise be checked completely literally.
    args.dir = os.path.expanduser(args.dir)
    if args.output:
        args.output = os.path.expanduser(args.output)

    if not os.path.isdir(args.dir):
        print(color_text(f"Error: Directory '{args.dir}' not found.", 'red'))
        sys.exit(1)

    manifest_path = os.path.join(args.dir, "run_manifest.json")
    if not os.path.isfile(manifest_path):
        print(color_text(
            f"Error: 'run_manifest.json' not found in '{args.dir}' -- this directory "
            "wasn't generated by stb-hubbardu.", 'red'))
        sys.exit(1)

    try:
        manifest = load_manifest(manifest_path)
    except (ValueError, json.JSONDecodeError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)
    species, n, l, j = manifest["species"], manifest["n"], manifest["l"], manifest["j"]

    scf_points, frozen_points = [], []
    frozen_iterations_seen = set()
    n_read = n_skipped = 0
    n_converged = n_not_converged = n_unknown = 0
    n_frozen_unknown = 0
    n_manifest_runs = len(manifest["runs"])
    print(f"\n{color_text('READING RUNS:', 'bold')}")
    for folder, info in manifest["runs"].items():
        if "kind" not in info or "alpha" not in info:
            print(f"   -> {folder:<24} : {color_text('SKIP', 'yellow')} (malformed manifest entry, "
                  "missing 'kind'/'alpha')")
            n_skipped += 1
            continue

        out_path = os.path.join(args.dir, folder, args.file)
        occ = parse_occupation(out_path)
        if occ is None:
            print(f"   -> {folder:<24} : {color_text('SKIP', 'yellow')} (missing or unparseable {args.file})")
            n_skipped += 1
            continue
        n_read += 1

        kind, alpha = info["kind"], info["alpha"]
        if kind in ("scf", "reference"):
            scf_points.append((alpha, occ))
        if kind == "frozen":
            frozen_points.append((alpha, occ))
            if "frozen_iterations" in info:
                frozen_iterations_seen.add(info["frozen_iterations"])

        converged = check_scf_converged(out_path)
        if kind == "frozen":
            # SCF_NOT_CONV is the EXPECTED outcome here, not a problem: the
            # frozen branch deliberately caps at MaxSCFIterations (1 by
            # default) to measure the bare/unrelaxed response, before the
            # density has a chance to reach real self-consistency. Treating
            # it the same as an scf_alpha_*/reference run not converging
            # (a genuine problem there) made every frozen run show up as
            # "NOT CONVERGED, used anyway" despite having completed exactly
            # as designed.
            if converged is None:
                n_frozen_unknown += 1
                conv_str, conv_color = "SCF: status unknown", 'yellow'
            elif converged is False:
                conv_str = "SCF: frozen evaluation complete (stopped at MaxSCFIterations, as designed)"
                conv_color = 'green'
            else:
                conv_str = "SCF: frozen evaluation converged before hitting MaxSCFIterations"
                conv_color = 'green'
        else:
            if converged is True:
                n_converged += 1
                conv_str, conv_color = "SCF: converged", 'green'
            elif converged is False:
                n_not_converged += 1
                conv_str, conv_color = "SCF: NOT CONVERGED, used anyway", 'yellow'
            else:
                n_unknown += 1
                conv_str, conv_color = "SCF: status unknown", 'yellow'
        print(f"   -> {folder:<24} : {color_text('OK', 'green')} (occupation={occ:.6f}, "
              f"{color_text(conv_str, conv_color)})")

    n_scf_read = len(scf_points)
    print(f"\n{color_text('Convergence summary:', 'cyan')} {n_read}/{n_manifest_runs} run(s) read"
          + (f" ({n_skipped} skipped)" if n_skipped else "") + f"; of the {n_scf_read} "
          f"self-consistent (scf/reference) run(s), {n_converged}/{n_scf_read} fully converged"
          + (f", {n_not_converged} NOT converged, {n_unknown} unknown"
             if n_not_converged or n_unknown else "")
          + (f". Frozen-density runs: {len(frozen_points) - n_frozen_unknown}/{len(frozen_points)} "
             "completed as expected" + (f", {n_frozen_unknown} unknown" if n_frozen_unknown else "")
             if frozen_points else ""))

    if len(frozen_iterations_seen) > 1:
        print(color_text(
            f"[WARNING] The frozen-density runs used different --frozen-iterations values across "
            f"this sweep ({sorted(frozen_iterations_seen)}) -- the frozen branch is only physically "
            "meaningful when every point uses the exact same recipe. Regenerate them with a "
            "consistent --frozen-iterations in stb-hubbarduAlphas.", 'yellow'))

    if len(scf_points) < 2 or len(frozen_points) < 2:
        print(color_text(
            "\nError: Not enough valid runs to fit a response (need at least 2 scf points, "
            "including the reference, and 2 frozen points).", 'red'))
        sys.exit(1)

    for name, points in (("self-consistent", scf_points), ("frozen-density", frozen_points)):
        if len(points) == 2:
            print(color_text(
                f"[WARNING] The {name} branch has only 2 points -- a line always fits 2 points "
                "exactly, so its R^2 will read 1.0 regardless of whether the response is actually "
                "linear. Add more --alphas in stb-hubbarduAlphas for a meaningful check.", 'yellow'))

    scf_points.sort()
    frozen_points.sort()
    chi, chi_intercept, chi_r2 = linear_fit_r2(*zip(*scf_points))
    chi0, chi0_intercept, chi0_r2 = linear_fit_r2(*zip(*frozen_points))

    print(f"\n{color_text('Self-consistent (screened) response chi:', 'cyan')} "
          f"{chi:.6f} 1/eV  (intercept={chi_intercept:.6f}, R^2={chi_r2:.4f})")
    print(f"{color_text('Frozen-density (bare) response chi0:', 'cyan')} "
          f"{chi0:.6f} 1/eV  (intercept={chi0_intercept:.6f}, R^2={chi0_r2:.4f})")

    for name, r2 in (("self-consistent", chi_r2), ("frozen-density", chi0_r2)):
        if r2 < args.r2_tolerance:
            print(color_text(
                f"[WARNING] The {name} response fit has R^2={r2:.4f}, below the "
                f"{args.r2_tolerance} tolerance -- the perturbation range may be too wide "
                "for the linear regime. Consider narrower --alphas in stb-hubbarduAlphas.", 'yellow'))

    # A decent overall R^2 can still hide a fit whose zero point drifted from
    # what was actually measured there (the fit tilting to average it out
    # rather than a genuinely noisy slope) -- compare the fitted intercept
    # against the actual alpha=0 occupation directly.
    scf_zero = next((occ for alpha, occ in scf_points if alpha == 0.0), None)
    frozen_zero = next((occ for alpha, occ in frozen_points if alpha == 0.0), None)
    for name, zero_occ, intercept in (("self-consistent", scf_zero, chi_intercept),
                                       ("frozen-density", frozen_zero, chi0_intercept)):
        if zero_occ is not None and abs(zero_occ - intercept) > INTERCEPT_DEVIATION_TOL:
            print(color_text(
                f"[WARNING] The {name} fit's intercept ({intercept:.6f}) differs from the "
                f"actually-measured alpha=0 occupation ({zero_occ:.6f}) by more than "
                f"{INTERCEPT_DEVIATION_TOL} -- the response may not be as linear as the R^2 "
                "alone suggests.", 'yellow'))

    if chi == 0 or chi0 == 0:
        print(color_text(
            "\nError: A zero response slope was fitted -- cannot compute U (division by "
            "zero). Check that the occupation actually responds to the perturbation.", 'red'))
        sys.exit(1)

    if (chi > 0) != (chi0 > 0):
        print(color_text(
            "[WARNING] chi and chi0 have OPPOSITE signs -- physically, both describe the same "
            "direction of response (just screened vs. bare), so this usually means a bad run, "
            "a mislabeled point, or a fit dominated by noise. Treat the computed U with caution.",
            'yellow'))

    U = 1.0 / chi0 - 1.0 / chi
    if U < 0:
        print(color_text(
            "[WARNING] Computed U is NEGATIVE -- the Hubbard U is a repulsive on-site term and "
            "should always be positive. This strongly suggests a bad fit or bad data rather "
            "than a real result.", 'yellow'))
    print(f"\n{color_text('Computed U (Cococcioni & de Gironcoli, PRB 71, 035105, 2005):', 'bold')} "
          f"{color_text(f'{U:.4f} eV', 'green')}")

    if species in REFERENCE_U:
        print(color_text(
            f"[INFO] Literature reference (GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 "
            f"(2006)) for {species}: {REFERENCE_U[species]:.2f} eV -- a sanity-check comparison "
            "only, not a validation (functional/pseudopotential/basis-dependent).", 'cyan'))

    output = args.output or f"{species}_LDAU.fdf"
    block = ldau_proj_block([{"species": species, "n": n, "l": l, "u": U, "j": j}])
    with open(output, 'w') as f:
        f.write(block)
    print(f"\n{color_text('[Saved]', 'cyan')} Ready-to-use DFT+U block -> {output}")

    data_filename = f"{species}_hubbardu_response.dat"
    gplot_filename = f"{species}_hubbardu_response.gplot"
    write_response_plot(data_filename, gplot_filename, scf_points, frozen_points,
                        chi, chi_intercept, chi0, chi0_intercept, species)
    print(f"{color_text('[Saved]', 'cyan')} Response plot -> {gplot_filename} "
          f"(gnuplot {gplot_filename})")


if __name__ == "__main__":
    main()
