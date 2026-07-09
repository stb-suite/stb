#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import os
import sys
import glob
import argparse
from stb.core import structure_io
from stb.core import xrd as xrd_core
from stb.core.cli import color_text, show_intro, run_with_spinner
from stb.core.deps import require_pyxtal


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


def main():
    require_pyxtal()
    from pyxtal.XRD import Similarity

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Ranks candidate structures by similarity to an experimental XRD pattern (Stage 2 - Analysis).", 'bold')}
Part of the Structure Solution (XRD) workflow: reads every candidate structure
in --input-dir (e.g. from stb-xrdsearch, one per subfolder -- picking up a
relaxed '*.STRUCT_OUT' over the raw 'structure.fdf' if you've already run
SIESTA there), simulates each one's powder XRD pattern, and ranks them by
similarity to an experimental pattern -- a fast pre-screen for which
candidate space group/arrangement is most likely correct. Run it once on
the raw candidates to prioritize which to relax with real DFT, then again
on the relaxed results to confirm. Each comparison takes 10s of seconds
(pyxtal's Similarity() has no faster mode), so this can take a while for a
large --input-dir.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --input-dir xrd_search --experimental measured.dat\n"
               "  %(prog)s --input-dir xrd_search --experimental measured.dat \\\n"
               "      --wavelength MoKa --top 5\n"
    )

    parser.add_argument("--input-dir", type=str, required=True,
                        help="Folder of candidate structures, one per subfolder (e.g. written "
                             "by stb-xrdsearch as '<subfolder>/structure.fdf'). If a subfolder "
                             "also has a '*.STRUCT_OUT' file (from running SIESTA there), the "
                             "relaxed structure is used instead of the raw one. Flat .fdf files "
                             "directly in --input-dir are also accepted.")
    parser.add_argument("--experimental", type=str, required=True, metavar="PATH",
                        help="Experimental XRD pattern: a plain text file with two columns "
                             "(2theta intensity), whitespace- or comma-separated, '#' comments "
                             "and blank lines allowed.")
    parser.add_argument("--wavelength", type=str, default="CuKa",
                        help="X-ray source: a known name (CuKa, CuKa1, MoKa, CrKa, FeKa, CoKa, "
                             "AgKa, ...) or a wavelength in Ang as a plain number "
                             "(default: CuKa, 1.54184 Ang).")
    parser.add_argument("--two-theta-range", type=float, nargs=2, default=[0.0, 90.0],
                        metavar=("MIN", "MAX"),
                        help="2-theta range in degrees to compute (default: 0 90).")
    parser.add_argument("--top", type=int, default=None,
                        help="Only show the N best-matching candidates (default: show all).")
    parser.add_argument("-o", "--output", type=str, default="xrd_rank.txt",
                        help="Output ranking report file name (default: xrd_rank.txt). Always "
                             "lists every candidate, regardless of --top.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-xrdrank {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Rank candidate structures against an experimental XRD pattern:", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.input_dir):
        print(color_text(f"Error: '{args.input_dir}' is not a directory.", 'red'))
        sys.exit(1)

    candidate_files = find_candidates(args.input_dir)
    if not candidate_files:
        print(color_text(
            f"Error: no candidate structures found in '{args.input_dir}' (expected "
            "'<subfolder>/structure.fdf', '<subfolder>/*.STRUCT_OUT', or flat .fdf files).",
            'red'))
        sys.exit(1)

    try:
        wavelength, wavelength_label = xrd_core.resolve_wavelength(args.wavelength)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    lo, hi = args.two_theta_range
    if lo >= hi:
        print(color_text("Error: --two-theta-range MIN must be less than MAX.", 'red'))
        sys.exit(1)

    try:
        experimental = xrd_core.read_experimental_pattern(args.experimental)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    print(f"  {color_text('Candidates found:', 'cyan')} {len(candidate_files)} in '{args.input_dir}'")
    print(f"  {color_text('Experimental pattern:', 'cyan')} {args.experimental}")
    print(f"  {color_text('Wavelength:', 'cyan')} {wavelength_label} ({wavelength:.5f} Ang)")
    print()

    results = []
    for name, path, fmt, source in candidate_files:
        try:
            pmg = structure_io.read_siesta_structure(path, fmt)
            pattern = xrd_core.compute_pattern(pmg, wavelength=wavelength, two_theta_range=(lo, hi))
        except (FileNotFoundError, ValueError) as e:
            print(color_text(f"  Warning: skipping '{name}' ({e}).", 'yellow'))
            continue

        sim_profile = pattern.get_profile()
        similarity = run_with_spinner(
            lambda sp=sim_profile: Similarity(sp, experimental).value, label=f"Comparing {name}")
        results.append((name, similarity, pmg.composition.reduced_formula, source))
        print(f"  {color_text('->', 'green')} {name}: similarity {similarity:.4f} "
              f"({pmg.composition.reduced_formula}, {source})")

    if not results:
        print(color_text(
            "\nError: none of the candidates could be compared -- see warnings above.", 'red'))
        sys.exit(1)

    results.sort(key=lambda r: r[1], reverse=True)
    shown = results if args.top is None else results[:args.top]

    # Sized to the longest candidate name -- --input-dir also accepts hand-built candidate
    # sets with arbitrary folder names, so a fixed width can glue the name into the next
    # column (the same bug crystalcast.py's own ranked table had to fix for the same reason).
    name_width = max(28, max((len(name) for name, *_ in results), default=0) + 2)

    def format_row(rank, name, similarity, formula, source):
        return f"{rank:<5} {name:<{name_width}}{similarity:<12.4f}{formula:<15}{source:<10}"

    header = f"{'Rank':<5} {'File':<{name_width}}{'Similarity':<12}{'Formula':<15}{'Source':<10}"

    print(f"\n{color_text(f'Ranking ({len(shown)} of {len(results)} shown, best match first):', 'bold')}")
    print(f"  {header}")
    for rank, (name, similarity, formula, source) in enumerate(shown, start=1):
        print(f"  {format_row(rank, name, similarity, formula, source)}")

    with open(args.output, "w") as f:
        f.write(header + "\n")
        for rank, (name, similarity, formula, source) in enumerate(results, start=1):
            f.write(format_row(rank, name, similarity, formula, source) + "\n")

    print(f"\n{color_text('Success:', 'green')} ranking ({len(results)} candidate(s)) written to "
          f"'{color_text(args.output, 'bold')}'")
    print(color_text(
        "\n  Note: a similarity ranking from simulated vs. experimental XRD, not a full "
        "Rietveld refinement -- use it to prioritize which candidate(s) to investigate "
        "further (e.g. relax with stb-mlrelax, then real DFT), not as definitive proof.",
        'yellow'))


if __name__ == "__main__":
    main()
