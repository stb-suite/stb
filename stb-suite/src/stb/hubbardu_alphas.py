#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.2.0"

import os
import sys
import json
import shutil
import argparse
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.dftu_data import load_manifest
from stb.hubbardu import (
    ldau_perturbation_block, write_run_folder, check_no_existing_ldau, copy_pseudopotentials,
    TEMPLATE_STRUCTURE_NAME, TEMPLATE_CALC_NAME, PSEUDO_SNAPSHOT_DIRNAME,
)

# Default perturbation strengths (eV) for the linear-response fit -- small
# enough to stay in the linear regime, symmetric around zero so the fit isn't
# skewed by an accidental asymmetry in the response.
DEFAULT_ALPHAS = [-0.15, -0.10, -0.05, 0.05, 0.10, 0.15]

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
    )

    parser.add_argument("--dir", type=str, default="hubbardu_runs",
                        help="Directory containing stage 1's 'reference/' folder and manifest "
                             "(default: hubbardu_runs).")
    parser.add_argument("--alphas", type=float, nargs='+', default=DEFAULT_ALPHAS,
                        help=f"Perturbation strengths in eV (default: {DEFAULT_ALPHAS}).")
    parser.add_argument("--frozen-iterations", type=int, default=DEFAULT_FROZEN_ITERATIONS,
                        help=f"MaxSCFIterations for the frozen-density runs (default: "
                             f"{DEFAULT_FROZEN_ITERATIONS}; raise this only if "
                             "stb-hubbarduAnalysis finds nothing to parse in those folders).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbarduAlphas {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Generate the Hubbard U perturbation sweep (stage 2/3):", 'bold'))
    print("-" * 60)

    # os.path.isdir never expands '~' itself (that's shell-only, and argv
    # passed through subprocess.run() never goes through a shell), so a path
    # like '~/hubbardu_runs' would otherwise be checked completely literally.
    args.dir = os.path.expanduser(args.dir)

    if args.frozen_iterations < 1:
        print(color_text("Error: --frozen-iterations must be >= 1.", 'red'))
        sys.exit(1)

    if len(set(args.alphas)) != len(args.alphas):
        print(color_text("Error: --alphas contains duplicate values.", 'red'))
        sys.exit(1)
    # Folder names round to 4 decimals (f"{alpha:.4f}"), so two alphas that
    # are exact-distinct floats can still collide on disk/in the manifest --
    # e.g. 0.00001 and the implicit alpha=0.0 frozen point both become
    # "frozen_alpha_0.0000", silently overwriting one of them.
    rounded = [f"{a:.4f}" for a in args.alphas]
    if len(set(rounded)) != len(rounded):
        print(color_text(
            "Error: two or more --alphas round to the same folder name at 4 decimal places -- "
            "use more distinct values.", 'red'))
        sys.exit(1)
    if "0.0000" in rounded:
        print(color_text(
            "Error: an alpha rounds to 0.0000 eV -- already covered automatically (the frozen "
            "branch gets its own alpha=0 evaluation, the scf branch reuses 'reference').", 'red'))
        sys.exit(1)

    if not os.path.isdir(args.dir):
        print(color_text(f"Error: Directory '{args.dir}' not found.", 'red'))
        sys.exit(1)

    manifest_path = os.path.join(args.dir, "run_manifest.json")
    if not os.path.isfile(manifest_path):
        print(color_text(
            f"Error: 'run_manifest.json' not found in '{args.dir}' -- run stb-hubbardu "
            "(stage 1) first.", 'red'))
        sys.exit(1)
    try:
        manifest = load_manifest(manifest_path)
    except (ValueError, json.JSONDecodeError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    if "reference" not in manifest["runs"]:
        print(color_text(f"Error: '{manifest_path}' has no 'reference' entry.", 'red'))
        sys.exit(1)

    species, n, l, j, label = manifest["species"], manifest["n"], manifest["l"], manifest["j"], manifest["label"]

    ref_dm = os.path.join(args.dir, "reference", f"{label}.DM")
    if not os.path.isfile(ref_dm):
        print(color_text(
            f"Error: '{ref_dm}' not found -- run SIESTA in '{args.dir}/reference/' first "
            "and make sure it converged (it must write a .DM file).", 'red'))
        sys.exit(1)

    structure_path = os.path.join(args.dir, TEMPLATE_STRUCTURE_NAME)
    calc_path = os.path.join(args.dir, TEMPLATE_CALC_NAME)
    if not os.path.isfile(structure_path) or not os.path.isfile(calc_path):
        print(color_text(
            f"Error: Template snapshot(s) not found in '{args.dir}' -- this directory "
            "wasn't generated by stb-hubbardu (stage 1).", 'red'))
        sys.exit(1)
    with open(calc_path, 'r') as f:
        calc_template = f.read()
    # Re-checked here, not just in stage 1: a user might hand-edit this saved
    # snapshot (e.g. to fix an SCF convergence problem in reference/) and
    # accidentally reintroduce a DFT+U block, which would otherwise silently
    # double up once stage 2 prepends its own.
    check_no_existing_ldau(calc_template, calc_path)

    print(f"  {color_text('Species / shell:', 'cyan')} {species} (n={n}, l={l}), J={j:.3f} eV")
    print(f"  {color_text('Perturbations (eV):', 'cyan')} {args.alphas}")

    pseudo_snapshot = os.path.join(args.dir, PSEUDO_SNAPSHOT_DIRNAME)
    have_pseudos = os.path.isdir(pseudo_snapshot)
    if not have_pseudos:
        print(color_text(
            "  [NOTE] No pseudopotentials saved by stage 1 (--pseudo-dir wasn't given) -- copy "
            "them into every scf_alpha_*/frozen_alpha_* folder yourself before running SIESTA.",
            'yellow'))

    for alpha in args.alphas:
        scf_name = f"scf_alpha_{alpha:.4f}"
        scf_extra = (ldau_perturbation_block(species, n, l, alpha, j, potential_shift=True)
                     + "DM.UseSaveDM T\n")
        scf_folder = write_run_folder(args.dir, scf_name, calc_template, structure_path, scf_extra)
        shutil.copy(ref_dm, os.path.join(scf_folder, f"{label}.DM"))
        if have_pseudos:
            copy_pseudopotentials(pseudo_snapshot, scf_folder)
        manifest["runs"][scf_name] = {"kind": "scf", "alpha": alpha}
        print(f"  {color_text('[OK]', 'green')} {scf_folder} (self-consistent response, seeded with reference DM)")

    for alpha in [0.0] + list(args.alphas):
        frozen_name = f"frozen_alpha_{alpha:.4f}"
        frozen_extra = (ldau_perturbation_block(species, n, l, alpha, j, potential_shift=True)
                        + f"MaxSCFIterations {args.frozen_iterations}\n" + "DM.UseSaveDM T\n")
        frozen_folder = write_run_folder(args.dir, frozen_name, calc_template, structure_path, frozen_extra)
        shutil.copy(ref_dm, os.path.join(frozen_folder, f"{label}.DM"))
        if have_pseudos:
            copy_pseudopotentials(pseudo_snapshot, frozen_folder)
        manifest["runs"][frozen_name] = {"kind": "frozen", "alpha": alpha,
                                          "frozen_iterations": args.frozen_iterations}
        print(f"  {color_text('[OK]', 'green')} {frozen_folder} (frozen-density response, reference DM copied in)")

    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    n_runs = 2 * len(args.alphas) + 1
    print(f"\n{color_text('Success:', 'green')} Generated {n_runs} run(s) in "
          f"'{color_text(args.dir, 'bold')}'")
    print(color_text("\nNext steps:", 'yellow'))
    print(f"  1. Run SIESTA in every 'scf_alpha_*' and 'frozen_alpha_*' folder.")
    print(f"  2. Run stb-hubbarduAnalysis --dir {args.dir}")


if __name__ == "__main__":
    main()
