#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import os
import re
import sys
import shutil
import argparse
import subprocess
from stb.core.cli import color_text, show_intro
from stb.core.deps import require_pyxtal, require_mace


def parse_groups(spec):
    """Parses '225,227,229' or '225 227 229' into [225, 227, 229]."""
    try:
        return [int(g) for g in spec.replace(",", " ").split()]
    except ValueError:
        raise ValueError(f"--groups must be a list of space group numbers, got '{spec}'.")


def strip_ansi(text):
    """Removes ANSI color codes -- stb-crystalcast's own stdout is already
    colored, and nesting it inside this tool's own color_text() calls
    without stripping first makes the two color-reset codes collide and the
    terminal output look broken.
    """
    return re.sub(r'\x1b\[[0-9;]*m', '', text)


def parse_ml_ranking(stdout):
    """Parses stb-crystalcast --ml-rank's own ranked-table stdout (see its
    'ML-ranked structures ...' block) into an ordered list of file
    basenames, most stable (lowest energy) first. Returns None if no ranking
    table is found in the output (e.g. every attempt in this batch failed
    before ranking, so --ml-rank never got to run).
    """
    lines = strip_ansi(stdout).splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if "ML-ranked structures" in line)
    except StopIteration:
        return None

    # The File column echoes crystalcast's own -o path verbatim, which can itself contain
    # spaces (e.g. a caller-supplied --output-dir with a space in it) -- greedily capture
    # everything up to the last two numeric columns instead of assuming \S+.
    row_re = re.compile(r'^\s*\d+\s+(.*\S)\s+[-\d.]+\s+[-\d.]+\s*$')
    ranked = []
    for line in lines[start + 2:]:  # +2 skips the "ML-ranked structures" line and its column header
        match = row_re.match(line)
        if not match:
            break
        ranked.append(match.group(1))
    return ranked or None


def main():
    require_pyxtal()  # fail fast here rather than after spawning stb-crystalcast subprocesses

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Casts candidate structures across a set of space groups (Stage 1 - Prep).", 'bold')}
Part of the Structure Solution (XRD) workflow: given a composition and a
list of space groups to try, casts --count-per-group random candidates in
each (by driving stb-crystalcast once per group), and puts each candidate
in its own folder as '<folder>/structure.fdf' -- ready for you to add your
own calc.fdf and run a real SIESTA relaxation in each one. Follow up with
stb-xrdrank to score every candidate against an experimental XRD pattern
(either now, on these raw candidates, or later, on your relaxed
.STRUCT_OUT results).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --species Ni O --num-ions 4 8 --groups 225,227,229\n"
               "  %(prog)s --species Na Cl --num-ions 4 4 --groups 225,221,62 \\\n"
               "      --count-per-group 3 -o candidates\n"
               "  %(prog)s --species Na Cl --num-ions 4 4 --groups 225,221,62 \\\n"
               "      --count-per-group 5 --ml-rank --keep-top-per-group 1\n"
    )

    parser.add_argument("--species", nargs="+", required=True,
                        help="Element symbols in the structure, e.g. --species Ni O.")
    parser.add_argument("--num-ions", nargs="+", type=int, required=True,
                        help="Number of atoms of each --species, same order and count, "
                             "e.g. --num-ions 4 8.")
    parser.add_argument("--groups", type=str, required=True,
                        help="Comma- or space-separated space group numbers to try, e.g. "
                             "'225,227,229' or '225 227 229'.")
    parser.add_argument("--count-per-group", type=int, default=1,
                        help="Random candidates to generate per space group (default: 1).")
    parser.add_argument("--molecular", action="store_true",
                        help="Passed through to stb-crystalcast: treat --species as whole "
                             "rigid molecules instead of bare elements.")
    parser.add_argument("--volume-factor", type=float, default=1.1,
                        help="Passed through to stb-crystalcast (default: 1.1).")
    parser.add_argument("--max-attempts", type=int, default=10,
                        help="Passed through to stb-crystalcast (default: 10).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Base random seed; each space group's batch uses seed+group so "
                             "groups don't share identical draws (default: not seeded).")
    parser.add_argument("--ml-rank", action="store_true",
                        help="Passed through to stb-crystalcast: quickly relax each candidate "
                             "with MACE-MP-0 and rank by energy (needs the optional 'ml' extra: "
                             "pip install stb_suite[ml]) -- a fast, local pre-screen, not a "
                             "substitute for the real SIESTA relaxation. Folders are then named "
                             "'<group>_rank<N>' (most stable first) instead of by attempt order.")
    parser.add_argument("--keep-top-per-group", type=int, default=None, metavar="N",
                        help="Only keep the N most stable candidates per group (requires "
                             "--ml-rank; default: keep all --count-per-group candidates).")
    parser.add_argument("-o", "--output-dir", type=str, default="xrd_search",
                        help="Folder to write candidate structures into (default: xrd_search); "
                             "each candidate gets its own '<output-dir>/<name>/structure.fdf'. "
                             "Created if it doesn't exist; refuses to run if it already has any "
                             "content, to avoid mixing batches.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-xrdsearch {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if len(args.species) != len(args.num_ions):
        parser.error(f"--species has {len(args.species)} entries but --num-ions has "
                     f"{len(args.num_ions)} -- they must match one-to-one.")
    if args.count_per_group < 1:
        parser.error("--count-per-group must be at least 1.")
    if args.keep_top_per_group is not None:
        if not args.ml_rank:
            parser.error("--keep-top-per-group requires --ml-rank (there's no ranking to keep "
                         "the top of otherwise).")
        if args.keep_top_per_group < 1:
            parser.error("--keep-top-per-group must be at least 1.")

    try:
        groups = parse_groups(args.groups)
    except ValueError as e:
        parser.error(str(e))
    if not groups:
        parser.error("--groups must list at least one space group number.")

    deduped_groups = list(dict.fromkeys(groups))
    if len(deduped_groups) != len(groups):
        print(color_text(
            "Warning: --groups had duplicate entries -- each space group is only cast once "
            "(duplicates would otherwise silently overwrite each other's candidate folder).",
            'yellow'))
        groups = deduped_groups

    if args.ml_rank:
        require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Cast candidate structures across space groups:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    if os.listdir(args.output_dir):
        print(color_text(
            f"Error: '{args.output_dir}' already has content -- move it aside or choose a "
            "different -o to avoid mixing batches.", 'red'))
        sys.exit(1)

    composition = " ".join(f"{s}{n}" for s, n in zip(args.species, args.num_ions))
    print(f"  {color_text('Composition:', 'cyan')} {composition}")
    print(f"  {color_text('Space groups:', 'cyan')} {', '.join(str(g) for g in groups)}")
    print(f"  {color_text('Candidates per group:', 'cyan')} {args.count_per_group}")
    if args.ml_rank:
        keep_note = f", keeping top {args.keep_top_per_group}" if args.keep_top_per_group else ""
        print(f"  {color_text('ML pre-screen:', 'cyan')} enabled{keep_note}")
    print()

    written = []
    failed = []
    for group in groups:
        tmp_stem = os.path.join(args.output_dir, f"_tmp_group_{group}")
        tmp_out = f"{tmp_stem}.fdf"
        cmd = [
            "stb-crystalcast", "--group", str(group), "--species", *args.species,
            "--num-ions", *[str(n) for n in args.num_ions],
            "--count", str(args.count_per_group),
            "--volume-factor", str(args.volume_factor),
            "--max-attempts", str(args.max_attempts),
            "-o", tmp_out, "--no-intro",
        ]
        if args.molecular:
            cmd.append("--molecular")
        if args.ml_rank:
            cmd.append("--ml-rank")
        if args.seed is not None:
            cmd.extend(["--seed", str(args.seed + group)])

        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            print(color_text(
                "\nError: 'stb-crystalcast' was not found on PATH -- make sure the stb_suite "
                "package is installed (pip install -e .) so its console scripts are available.",
                'red'))
            sys.exit(1)

        if args.count_per_group == 1:
            raw_candidates = [(1, tmp_out)] if os.path.exists(tmp_out) else []
        else:
            raw_candidates = [(i, f"{tmp_stem}_{i}.fdf") for i in range(1, args.count_per_group + 1)
                              if os.path.exists(f"{tmp_stem}_{i}.fdf")]

        if not raw_candidates:
            tail = (strip_ansi(result.stdout).strip().splitlines() or ["(no output)"])[-1]
            print(color_text(f"  -> group {group}: FAILED ({tail})", 'yellow'))
            failed.append(group)
            continue

        # With --ml-rank, name folders by rank (most stable first) instead of attempt
        # order. Any raw candidate that doesn't make it into the final set -- it failed
        # to relax, or --keep-top-per-group pruned it -- is deleted here rather than left
        # behind as an unexplained flat file directly under --output-dir.
        if args.ml_rank:
            ranking = parse_ml_ranking(result.stdout)
            if ranking:
                by_basename = {os.path.basename(path): (i, path) for i, path in raw_candidates}
                ordered = [by_basename[os.path.basename(name)] for name in ranking
                          if os.path.basename(name) in by_basename]
                if args.keep_top_per_group is not None:
                    ordered = ordered[:args.keep_top_per_group]
                folder_names = [f"group_{group}_rank{rank}" for rank in range(1, len(ordered) + 1)]
                candidates = [path for _, path in ordered]
                kept_paths = set(candidates)
                for _, path in raw_candidates:
                    if path not in kept_paths:
                        os.remove(path)
            else:
                if args.keep_top_per_group is not None:
                    print(color_text(
                        f"  Warning: group {group}: could not parse an ML ranking from "
                        "stb-crystalcast's output -- --keep-top-per-group not applied, "
                        f"keeping all {len(raw_candidates)} candidate(s).", 'yellow'))
                folder_names = [f"group_{group}" if args.count_per_group == 1 else f"group_{group}_{i}"
                                for i, _ in raw_candidates]
                candidates = [path for _, path in raw_candidates]
        else:
            folder_names = [f"group_{group}" if args.count_per_group == 1 else f"group_{group}_{i}"
                            for i, _ in raw_candidates]
            candidates = [path for _, path in raw_candidates]

        if not candidates:
            print(color_text(
                f"  -> group {group}: FAILED (every candidate was excluded by ML ranking/"
                "relaxation -- none left to keep)", 'yellow'))
            failed.append(group)
            continue

        for folder_name, tmp_path in zip(folder_names, candidates):
            folder_path = os.path.join(args.output_dir, folder_name)
            os.makedirs(folder_path, exist_ok=True)
            shutil.move(tmp_path, os.path.join(folder_path, "structure.fdf"))
            written.append(folder_path)

        note = "" if len(candidates) == args.count_per_group else \
            f" ({args.count_per_group - len(candidates)} of {args.count_per_group} skipped/dropped)"
        print(f"  {color_text('->', 'green')} group {group}: {len(candidates)} folder(s){note} "
              f"({', '.join(folder_names)})")

    if not written:
        print(color_text("\nError: no candidates were generated -- see failures above.", 'red'))
        sys.exit(1)

    if failed:
        print(color_text(
            f"\nPartial success: {len(written)} candidate folder(s) written to '{args.output_dir}' "
            f"({len(failed)} of {len(groups)} group(s) failed, see above).", 'yellow'))
        print(f"  {color_text('Next:', 'cyan')} run SIESTA in each folder, then "
              f"stb-xrdrank --input-dir {args.output_dir} --experimental <file>")
        sys.exit(1)

    print(f"\n{color_text('Success:', 'green')} {len(written)} candidate folder(s) written to "
          f"'{color_text(args.output_dir, 'bold')}'.")
    print(f"  {color_text('Next:', 'cyan')} add calc.fdf and run SIESTA in each folder, then "
          f"stb-xrdrank --input-dir {args.output_dir} --experimental <file>")


if __name__ == "__main__":
    main()
