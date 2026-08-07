#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

# 1.3.0: numbered-section report (matching stage 1/stb-hubbardu and the other
# Stage 1 prep tools' convention); --save-report. --dir's default
# ("hubbardu_runs") already matched stage 1's own default --output-dir, so
# the common case (stage 2 right after stage 1, same folder) needed no
# extra flag before this change either -- unchanged, just confirmed/reported
# explicitly now in [0] RUN METADATA.
VERSION = "1.3.0"

import os
import sys
import json
import shutil
import argparse
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.dftu_data import load_manifest
from stb.hubbardu import (
    ldau_perturbation_block, write_run_folder, check_no_existing_ldau, copy_pseudopotentials,
    TEMPLATE_STRUCTURE_NAME, TEMPLATE_CALC_NAME, PSEUDO_SNAPSHOT_DIRNAME,
)

REPORT_FILE = "hubbardu_stage2.txt"

# Default perturbation strengths (eV) for the linear-response fit -- symmetric
# around zero so the fit isn't skewed by an accidental asymmetry in the
# response. Narrowed from +-0.15 to +-0.10: on a real MnO run, including the
# +-0.15 points pulled the fitted U from ~4.4 eV (using only up to +-0.10) to
# ~9 eV -- a sign the response is no longer linear that far out. +-0.10 was
# the widest range that stayed consistent there. Still just a default --
# widen via --alphas for a system with a weak response (signal near numerical
# noise at +-0.10), and check stage 3's R^2/intercept warnings either way.
DEFAULT_ALPHAS = [-0.10, -0.05, 0.05, 0.10]

# With DM.UseSaveDM T seeding a converged reference .DM, a single SCF
# iteration is enough: SIESTA still prints the DFT+U occupation diagnostic
# and exits normally -- verified against real serial SIESTA 5.4.2 output
# (frozen_alpha_0.0000 with MaxSCFIterations 1 finished with 0_NORMAL_EXIT
# and the hubbard_term occupations printed). Exposed as --frozen-iterations
# (not just this hardcoded default) in case a different SIESTA version/
# parallel run needs more than one iteration to trigger the same print.
DEFAULT_FROZEN_ITERATIONS = 1


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 3: generates the perturbed run folders for a Hubbard U linear-response calculation.", 'bold')}
Reads the 'reference/' folder and template snapshot written by stb-hubbardu
(stage 1) -- 'reference/' must already be converged, since its .DM is copied
into every generated folder here. Writes one 'scf_alpha_<v>' folder per
perturbation strength (full self-consistent response) and one
'frozen_alpha_<v>' folder per strength INCLUDING v=0 (single-iteration,
frozen-density response) -- the v=0 frozen point is evaluated by the exact
same recipe as the others rather than reusing reference's own fully-converged
value, so every point on the frozen-response line is mutually consistent.
Also copies in any pseudopotentials stb-hubbardu was given via --pseudo-dir.
Doesn't run SIESTA itself -- run each folder's calculation yourself, then use
stb-hubbarduAnalysis.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --dir hubbardu_runs\n"
               "  %(prog)s --dir hubbardu_runs --alphas -0.2 -0.1 0.1 0.2\n"
               "  %(prog)s --dir hubbardu_runs --save-report\n"
    )

    parser.add_argument("--dir", type=str, default="hubbardu_runs",
                        help="Directory containing stage 1's 'reference/' folder and manifest -- "
                             "defaults to 'hubbardu_runs', stb-hubbardu's own default "
                             "--output-dir, so stage 2 runs against stage 1's output with no "
                             "extra flag in the common case.")
    parser.add_argument("--alphas", type=float, nargs='+', default=DEFAULT_ALPHAS,
                        help=f"Perturbation strengths in eV (default: {DEFAULT_ALPHAS}).")
    parser.add_argument("--frozen-iterations", type=int, default=DEFAULT_FROZEN_ITERATIONS,
                        help=f"MaxSCFIterations for the frozen-density runs (default: "
                             f"{DEFAULT_FROZEN_ITERATIONS}; raise this only if "
                             "stb-hubbarduAnalysis finds nothing to parse in those folders).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbarduAlphas {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.frozen_iterations < 1:
        sys.exit(color_text("Error: --frozen-iterations must be >= 1.", 'red'))

    if len(set(args.alphas)) != len(args.alphas):
        sys.exit(color_text("Error: --alphas contains duplicate values.", 'red'))
    # Folder names round to 4 decimals (f"{alpha:.4f}"), so two alphas that
    # are exact-distinct floats can still collide on disk/in the manifest --
    # e.g. 0.00001 and the implicit alpha=0.0 frozen point both become
    # "frozen_alpha_0.0000", silently overwriting one of them.
    rounded = [f"{a:.4f}" for a in args.alphas]
    if len(set(rounded)) != len(rounded):
        sys.exit(color_text(
            "Error: two or more --alphas round to the same folder name at 4 decimal places -- "
            "use more distinct values.", 'red'))
    if "0.0000" in rounded:
        sys.exit(color_text(
            "Error: an alpha rounds to 0.0000 eV -- already covered automatically (the frozen "
            "branch gets its own alpha=0 evaluation, the scf branch reuses 'reference').", 'red'))

    # os.path.isdir never expands '~' itself (that's shell-only, and argv
    # passed through subprocess.run() never goes through a shell), so a path
    # like '~/hubbardu_runs' would otherwise be checked completely literally.
    args.dir = os.path.expanduser(args.dir)
    if not os.path.isdir(args.dir):
        sys.exit(color_text(f"Error: Directory '{args.dir}' not found.", 'red'))

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
        "===== STB-HUBBARDUALPHAS STAGE 2 REPORT (HUBBARD U PERTURBATION PREP) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Directory              : {args.dir}", f_out)
    print_dual(f"Perturbations (eV)     : {args.alphas}", f_out)
    print_dual(f"Frozen iterations      : {args.frozen_iterations}", f_out)
    print_dual(f"Report                 : {report_path if report_path else '(not saved)'}", f_out)

    print_section("[1] VALIDATION & REFERENCE CHECK", f_out)
    manifest_path = os.path.join(args.dir, "run_manifest.json")
    if not os.path.isfile(manifest_path):
        fail(f"'run_manifest.json' not found in '{args.dir}' -- run stb-hubbardu (stage 1) first.")
    try:
        manifest = load_manifest(manifest_path)
    except (ValueError, json.JSONDecodeError) as e:
        fail(str(e))

    if "reference" not in manifest["runs"]:
        fail(f"'{manifest_path}' has no 'reference' entry.")

    species, n, l, j, label = manifest["species"], manifest["n"], manifest["l"], manifest["j"], manifest["label"]
    # Falls back to `species` for manifests written before stage 1 supported
    # --atom-index (no aliasing then, so the species itself was always the
    # one actually carrying the LDAU.proj perturbation).
    perturbed_species = manifest.get("perturbed_species", species)
    print_dual(f"Species / shell        : {species} (n={n}, l={l}), J={j:.3f} eV", f_out)
    if perturbed_species != species:
        print_dual(f"Perturbed atom         : isolated via alias species '{perturbed_species}' "
                   "(see stage 1's --atom-index)", f_out)

    ref_dm = os.path.join(args.dir, "reference", f"{label}.DM")
    if not os.path.isfile(ref_dm):
        fail(f"'{ref_dm}' not found -- run SIESTA in '{args.dir}/reference/' first and make "
             "sure it converged (it must write a .DM file).")
    print_dual(f"Reference DM           : {ref_dm} (found, will be copied into every folder)", f_out)

    structure_path = os.path.join(args.dir, TEMPLATE_STRUCTURE_NAME)
    calc_path = os.path.join(args.dir, TEMPLATE_CALC_NAME)
    if not os.path.isfile(structure_path) or not os.path.isfile(calc_path):
        fail(f"Template snapshot(s) not found in '{args.dir}' -- this directory wasn't "
             "generated by stb-hubbardu (stage 1).")
    with open(calc_path, 'r') as f:
        calc_template = f.read()
    # Re-checked here, not just in stage 1: a user might hand-edit this saved
    # snapshot (e.g. to fix an SCF convergence problem in reference/) and
    # accidentally reintroduce a DFT+U block, which would otherwise silently
    # double up once stage 2 prepends its own. Its own plain print+exit (not
    # routed through f_out) is a known, accepted minor gap shared with stage
    # 1 -- see hubbardu.py's own check_no_existing_ldau docstring.
    check_no_existing_ldau(calc_template, calc_path)
    print_dual("Template snapshot      : OK (no pre-existing DFT+U block)", f_out)

    print_section("[2] PSEUDOPOTENTIALS", f_out)
    pseudo_snapshot = os.path.join(args.dir, PSEUDO_SNAPSHOT_DIRNAME)
    have_pseudos = os.path.isdir(pseudo_snapshot)
    if have_pseudos:
        print_dual(f"Source                 : {pseudo_snapshot} (saved by stage 1's --pseudo-dir)", f_out)
        print_dual(color_text(
            "[OK] Will be copied into every scf_alpha_*/frozen_alpha_* folder.", 'green'), f_out)
    else:
        print_dual(color_text(
            "[NOTE] No pseudopotentials saved by stage 1 (--pseudo-dir wasn't given) -- copy "
            "them into every scf_alpha_*/frozen_alpha_* folder yourself before running SIESTA.",
            'yellow'), f_out)

    print_section("[3] PERTURBATION FOLDERS", f_out)
    for alpha in args.alphas:
        scf_name = f"scf_alpha_{alpha:.4f}"
        scf_extra = (ldau_perturbation_block(perturbed_species, n, l, alpha, j, potential_shift=True)
                     + "DM.UseSaveDM T\n")
        scf_folder, _ = write_run_folder(args.dir, scf_name, calc_template, structure_path, scf_extra)
        shutil.copy(ref_dm, os.path.join(scf_folder, f"{label}.DM"))
        if have_pseudos:
            copy_pseudopotentials(pseudo_snapshot, scf_folder)
        manifest["runs"][scf_name] = {"kind": "scf", "alpha": alpha}
        print_dual(color_text(
            f"[OK] {scf_folder} (self-consistent response, seeded with reference DM)", 'green'), f_out)

    for alpha in [0.0] + list(args.alphas):
        frozen_name = f"frozen_alpha_{alpha:.4f}"
        frozen_extra = (ldau_perturbation_block(perturbed_species, n, l, alpha, j, potential_shift=True)
                        + f"MaxSCFIterations {args.frozen_iterations}\n" + "DM.UseSaveDM T\n")
        frozen_folder, _ = write_run_folder(args.dir, frozen_name, calc_template, structure_path, frozen_extra)
        shutil.copy(ref_dm, os.path.join(frozen_folder, f"{label}.DM"))
        if have_pseudos:
            copy_pseudopotentials(pseudo_snapshot, frozen_folder)
        manifest["runs"][frozen_name] = {"kind": "frozen", "alpha": alpha,
                                          "frozen_iterations": args.frozen_iterations}
        print_dual(color_text(
            f"[OK] {frozen_folder} (frozen-density response, reference DM copied in)", 'green'), f_out)

    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print_section("[4] SUMMARY & NEXT STEPS", f_out)
    n_runs = 2 * len(args.alphas) + 1
    print_dual(color_text(f"Success: Generated {n_runs} run(s) in '{args.dir}'", 'green'), f_out)
    if report_path:
        print_dual(f"Report                 : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual("  1. Run SIESTA in every 'scf_alpha_*' and 'frozen_alpha_*' folder.", f_out)
    print_dual(f"  2. Run stb-hubbarduAnalysis --dir {args.dir}", f_out)

    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
