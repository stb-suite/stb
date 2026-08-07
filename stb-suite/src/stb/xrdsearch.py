#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

# 1.14.0: richer report -- Date/time in [0], Volume/Max displacement columns in the
# [3] MACE PRE-RELAXATION table, and a new final [7] LIBRARY WARNINGS section that
# collects every third-party warning (MACE/torch/pymatgen/spglib/stb-crystalcast)
# emitted during the run instead of letting it leak into the terminal interleaved
# with the report.
VERSION = "1.14.0"

import os
import re
import io
import sys
import shutil
import argparse
import subprocess
import contextlib
from datetime import datetime
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                          capture_library_noise)
from stb.core.deps import require_pyxtal, require_mace
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos, copy_pseudo
from stb.core import structure_io, kspace, mace_relax
# First reuse of inputfile.py's generate_calculation() as a library call (previously only
# ever invoked via the stb-inputfile CLI). It writes 'calc.fdf' (and a citations .bib file)
# relative to os.getcwd() and prints its own full report -- called below with os.chdir()
# into each candidate folder, and its stdout captured/re-summarized rather than left to
# flood this tool's own report.
from stb.inputfile import generate_calculation

REPORT_FILE = "xrdsearch_stage1.txt"
MACE_VACUUM_GAP_ANG = 10.0  # matches stb-kgrid/stb-mlrelax's own default


def parse_groups(spec):
    """Parses '225,227,229' or '225 227 229' into [225, 227, 229]."""
    try:
        return [int(g) for g in spec.replace(",", " ").split()]
    except ValueError:
        raise ValueError(f"--groups must be a list of space group numbers, got '{spec}'.")


def strip_ansi(text):
    """Removes ANSI color codes -- stb-crystalcast's own stdout (and
    inputfile.generate_calculation's own captured report) is already
    colored, and nesting it inside this tool's own color_text() calls
    without stripping first makes the two color-reset codes collide and the
    terminal output look broken.
    """
    return re.sub(r'\x1b\[[0-9;]*m', '', text)


def extract_warning_lines(text):
    """Pulls out genuine Python warnings.warn() lines (stdlib's default
    format: '<file>:<line>: <Category>Warning: <message>') from a
    subprocess's captured stderr -- distinct from stb-crystalcast's own
    '[WARNING]'-tagged domain messages (already surfaced elsewhere, e.g. in
    [2] CANDIDATE CASTING), which don't match this format. Used to route
    real third-party library noise (pymatgen/pyxtal/spglib deprecation
    -style warnings) from the stb-crystalcast subprocess into the final
    LIBRARY WARNINGS section instead of silently discarding it.
    """
    return [line.strip() for line in text.splitlines() if re.search(r'\w+Warning:', line)]


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
in its own folder as '<folder>/structure.fdf'. Optionally pre-relaxes every
final candidate with a MACE potential (--mace-relax, positions + cell) --
a fast cleanup of pyxtal's raw random placement, not a substitute for a
real SIESTA relaxation. By default, also generates the standard relaxation
calc.fdf (same as stb-inputfile -t relax) in every folder -- use --calc to
supply your own template instead -- and links pseudopotentials if
--pseudo-dir is given, so every folder is ready for a real SIESTA
relaxation as-is. Follow up with stb-xrdrank to score every candidate
against an experimental XRD pattern (either now, on these raw candidates,
or later, on your relaxed .STRUCT_OUT results).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --species Ni O --num-ions 4 8 --groups 225,227,229\n"
               "  %(prog)s --species Na Cl --num-ions 4 4 --groups 225,221,62 \\\n"
               "      --count-per-group 3 -o candidates --pseudo-dir dojo\n"
               "  %(prog)s --species Na Cl --num-ions 4 4 --groups 225,221,62 \\\n"
               "      --count-per-group 5 --ml-rank --keep-top-per-group 1\n"
               "  %(prog)s --species Ni O --num-ions 4 8 --groups 225,227 \\\n"
               "      --calc my_calc_template.fdf --pseudo-dir /path/to/pps --view\n"
               "  %(prog)s --species Ni O --num-ions 4 8 --groups 225,227 \\\n"
               "      --mace-relax --mace-model medium --pseudo-dir dojo\n"
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
    parser.add_argument("--mace-relax", action="store_true",
                        help="Pre-relax every FINAL candidate structure (positions + cell) "
                             "with a MACE potential before generating calc.fdf -- a fast "
                             "heuristic cleanup of pyxtal's raw random placement, not a "
                             "substitute for the real SIESTA relaxation. The relaxed geometry "
                             "overwrites each candidate's structure.fdf. Distinct from "
                             "--ml-rank, which only relaxes (positions-only, transiently) to "
                             "RANK/select candidates during casting -- the geometry it uses "
                             "for that is discarded, not kept. Needs the optional 'ml' extra "
                             "(pip install stb_suite[ml]).")
    parser.add_argument("--mace-model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size for --mace-relax (default: small). Ignored "
                             "if --mace-custom-model is given. Only valid together with "
                             "--mace-relax.")
    parser.add_argument("--mace-custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --mace-relax, "
                             "instead of a MACE-MP-0 foundation size. Only valid together "
                             "with --mace-relax.")
    parser.add_argument("--mace-device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the MACE model on (default: cpu) -- used for both "
                             "--mace-relax (here) and --ml-rank (passed through to the "
                             "stb-crystalcast subprocess).")
    parser.add_argument("--mace-fmax", type=float, default=0.05,
                        help="Force convergence threshold for --mace-relax, eV/Ang "
                             "(default: 0.05).")
    parser.add_argument("--mace-max-steps", type=int, default=200,
                        help="Maximum optimizer steps for --mace-relax (default: 200).")
    parser.add_argument("-c", "--calc", type=str, default=None, metavar="PATH",
                        help="Custom calc.fdf template, copied (as 'calc.fdf') into every "
                             "candidate folder -- should %%include 'structure.fdf' by name. "
                             "If omitted (default), the standard relaxation calc.fdf -- the "
                             "same default stb-inputfile -t relax produces -- is auto-generated "
                             "for every candidate instead (see --d3/--spin-polarized).")
    parser.add_argument("--d3", action="store_true",
                        help="Auto-generate mode only (--calc not given): enable the DFT-D3 "
                             "(Grimme) van der Waals correction in the generated calc.fdf. "
                             "Rejected together with --calc.")
    parser.add_argument("-s", "--spin-polarized", dest="spin_polarized", action="store_true",
                        help="Auto-generate mode only (--calc not given): enable spin "
                             "polarization in the generated calc.fdf. Rejected together with "
                             "--calc.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="", dest="pseudo_dir",
                        metavar="PATH_OR_BANK",
                        help=f"Pseudopotentials source (optional): a bundled bank "
                             f"({', '.join(BANKS)}) or a folder path. If given, the required "
                             "pseudopotential for every species actually present in the cast "
                             "candidates is copied into every candidate folder.")
    parser.add_argument("--view", action="store_true",
                        help="After generating, open one interactive 3D view (via ASE) paging "
                             "through every candidate structure. Needs a display. Off by "
                             "default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by "
                             "default.")
    parser.add_argument("-o", "--output-dir", type=str, default="xrd_search",
                        help="Folder to write candidate structures into (default: xrd_search); "
                             "each candidate gets its own '<output-dir>/<name>/'. Created if it "
                             "doesn't exist; refuses to run if it already has any content, to "
                             "avoid mixing batches.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-xrdsearch {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()
    library_warnings = []

    if args.calc and (args.d3 or args.spin_polarized):
        parser.error("--d3/--spin-polarized only apply when auto-generating calc.fdf -- they "
                     "have no effect together with --calc (a custom template is copied "
                     "verbatim); drop --calc or drop --d3/--spin-polarized.")
    if not args.mace_relax and (args.mace_custom_model or args.mace_model != "small"):
        parser.error("--mace-model/--mace-custom-model are only valid together with "
                     "--mace-relax.")

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

    if args.calc and not os.path.isfile(args.calc):
        sys.exit(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
    if args.mace_custom_model and not os.path.isfile(args.mace_custom_model):
        sys.exit(color_text(
            f"[ERROR] --mace-custom-model file not found: {args.mace_custom_model}", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))

    if args.ml_rank or args.mace_relax:
        with capture_library_noise(library_warnings, "MACE import"):
            require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    if os.listdir(args.output_dir):
        print(color_text(
            f"Error: '{args.output_dir}' already has content -- move it aside or choose a "
            "different -o to avoid mixing batches.", 'red'))
        sys.exit(1)

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-XRDSEARCH STAGE 1 REPORT (STRUCTURE SOLUTION PREP) =====", 'magenta'), f_out)

    composition = " ".join(f"{s}{n}" for s, n in zip(args.species, args.num_ions))

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time              : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Composition            : {composition}", f_out)
    print_dual(f"Space groups           : {', '.join(str(g) for g in groups)}", f_out)
    print_dual(f"Candidates per group   : {args.count_per_group}", f_out)
    if args.ml_rank:
        keep_note = f", keeping top {args.keep_top_per_group}" if args.keep_top_per_group else ""
        print_dual(f"ML pre-screen          : enabled{keep_note}", f_out)
    else:
        print_dual("ML pre-screen          : disabled", f_out)
    print_dual(f"MACE pre-relaxation    : {'enabled' if args.mace_relax else 'disabled'}", f_out)
    print_dual(f"Output directory       : {args.output_dir}", f_out)
    print_dual(f"Report                 : {report_path if report_path else '(not saved)'}", f_out)

    print_section("[1] CALC.FDF SETUP", f_out)
    if args.calc:
        print_dual(f"Mode                   : custom template ({args.calc})", f_out)
        with open(args.calc, "r") as fh:
            calc_text_check = fh.read()
        if "%include" not in calc_text_check or "structure.fdf" not in calc_text_check:
            print_dual(color_text(
                "[NOTE] The custom calc.fdf template does not appear to '%include "
                "structure.fdf' -- make sure it references the per-candidate structure file "
                "by that name (every candidate folder's structure is written as "
                "'structure.fdf').", 'yellow'), f_out)
    else:
        print_dual("Mode                   : auto-generate (standard relaxation calc.fdf, "
                   "same as stb-inputfile -t relax)", f_out)
        print_dual(f"DFT-D3 correction      : {'ENABLED' if args.d3 else 'disabled'}", f_out)
        print_dual(
            f"Spin polarization      : {'polarized' if args.spin_polarized else 'non-polarized'}",
            f_out)

    print_section("[2] CANDIDATE CASTING", f_out)
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
            cmd.extend(["--device", args.mace_device])
        if args.seed is not None:
            cmd.extend(["--seed", str(args.seed + group)])

        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            fail("'stb-crystalcast' was not found on PATH -- make sure the stb_suite package "
                 "is installed (pip install -e .) so its console scripts are available.")

        for line in extract_warning_lines(strip_ansi(result.stderr)):
            entry = f"[stb-crystalcast group {group}] {line}"
            if entry not in library_warnings:
                library_warnings.append(entry)

        if args.count_per_group == 1:
            raw_candidates = [(1, tmp_out)] if os.path.exists(tmp_out) else []
        else:
            raw_candidates = [(i, f"{tmp_stem}_{i}.fdf") for i in range(1, args.count_per_group + 1)
                              if os.path.exists(f"{tmp_stem}_{i}.fdf")]

        if not raw_candidates:
            tail = (strip_ansi(result.stdout).strip().splitlines() or ["(no output)"])[-1]
            print_dual(f"  {color_text('->', 'yellow')} group {group}: FAILED ({tail})", f_out)
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
                    print_dual(color_text(
                        f"  Warning: group {group}: could not parse an ML ranking from "
                        "stb-crystalcast's output -- --keep-top-per-group not applied, "
                        f"keeping all {len(raw_candidates)} candidate(s).", 'yellow'), f_out)
                folder_names = [f"group_{group}" if args.count_per_group == 1 else f"group_{group}_{i}"
                                for i, _ in raw_candidates]
                candidates = [path for _, path in raw_candidates]
        else:
            folder_names = [f"group_{group}" if args.count_per_group == 1 else f"group_{group}_{i}"
                            for i, _ in raw_candidates]
            candidates = [path for _, path in raw_candidates]

        if not candidates:
            print_dual(color_text(
                f"  -> group {group}: FAILED (every candidate was excluded by ML ranking/"
                "relaxation -- none left to keep)", 'yellow'), f_out)
            failed.append(group)
            continue

        for folder_name, tmp_path in zip(folder_names, candidates):
            folder_path = os.path.join(args.output_dir, folder_name)
            os.makedirs(folder_path, exist_ok=True)
            shutil.move(tmp_path, os.path.join(folder_path, "structure.fdf"))
            written.append(folder_path)

        note = "" if len(candidates) == args.count_per_group else \
            f" ({args.count_per_group - len(candidates)} of {args.count_per_group} skipped/dropped)"
        print_dual(f"  {color_text('->', 'green')} group {group}: {len(candidates)} folder(s){note} "
              f"({', '.join(folder_names)})", f_out)

    if not written:
        fail("No candidates were generated -- see [2] above.")

    print_section("[3] MACE PRE-RELAXATION", f_out)
    relax_failures = []
    if not args.mace_relax:
        print_dual("Not requested (pass --mace-relax to pre-relax every final candidate "
                   "structure -- positions + cell -- with a MACE potential before "
                   "generating calc.fdf).", f_out)
    else:
        mace_model_arg = args.mace_custom_model if args.mace_custom_model else args.mace_model
        mace_desc = (f"custom ({args.mace_custom_model})" if args.mace_custom_model
                     else f"MACE-MP-0 ({args.mace_model})")
        print_dual(f"Model                  : {mace_desc}", f_out)
        print_dual(f"Device                 : {args.mace_device}", f_out)
        print_dual(f"fmax / max steps       : {args.mace_fmax} eV/Ang / {args.mace_max_steps}", f_out)
        with capture_library_noise(library_warnings, "MACE calculator setup"):
            try:
                mace_calc = mace_relax.get_calculator(model=mace_model_arg, device=args.mace_device)
            except ValueError as e:
                fail(str(e))
        for line in mace_relax.describe_model(mace_model_arg, mace_calc):
            print_dual(line, f_out)

        relax_rows = []
        with capture_library_noise(library_warnings, "MACE relax"):
            for folder_path in written:
                folder_name = os.path.basename(folder_path)
                struct_path = os.path.join(folder_path, "structure.fdf")
                try:
                    orig_structure = structure_io.read_fdf(struct_path)
                    pmg = structure_io.to_pymatgen(orig_structure)
                    frac_coords = [site.frac_coords for site in pmg]
                    vacuum_axes = kspace.detect_vacuum_axes(
                        frac_coords, pmg.lattice.matrix, MACE_VACUUM_GAP_ANG)
                    cell_mask = None if all(vacuum_axes) else mace_relax.build_cell_mask(vacuum_axes)

                    atoms = AseAtomsAdaptor.get_atoms(pmg)
                    atoms_before = atoms.copy()
                    atoms.calc = mace_calc
                    e0 = atoms.get_potential_energy()
                    converged, steps_used = mace_relax.relax(
                        atoms, mace_calc, cell_mask=cell_mask,
                        fmax=args.mace_fmax, max_steps=args.mace_max_steps)
                    e1 = atoms.get_potential_energy()

                    vol0, vol1 = atoms_before.get_volume(), atoms.get_volume()
                    vol_change_pct = 100 * (vol1 - vol0) / vol0 if vol0 else 0.0
                    max_disp = np.linalg.norm(
                        atoms.get_positions() - atoms_before.get_positions(), axis=1).max()

                    new_pmg = AseAtomsAdaptor.get_structure(atoms)
                    new_structure = structure_io.from_pymatgen(
                        new_pmg, species_meta=orig_structure.species_meta)
                    structure_io.write_fdf(new_structure, struct_path, header_comment=[
                        f"Structure cast by stb-crystalcast, MACE pre-relaxed by stb-xrdsearch "
                        f"using {mace_desc}.",
                        f"Optimizer: FIRE, {steps_used} step(s), "
                        f"{'converged' if converged else 'NOT converged (hit step cap)'}.",
                        f"Energy: {e1:.6f} eV (delta {e1 - e0:+.6f} eV from the cast structure).",
                        f"Volume: {vol0:.4f} -> {vol1:.4f} Ang^3 ({vol_change_pct:+.2f}%). "
                        f"Max displacement: {max_disp:.4f} Ang.",
                    ])
                    status = "converged" if converged else "NOT converged"
                    relax_rows.append((
                        [folder_name, f"{e0:.4f}", f"{e1:.4f}", f"{vol_change_pct:+.2f}%",
                         f"{max_disp:.4f}", f"{steps_used} ({status})"],
                        None if converged else 'yellow'))
                except Exception as e:
                    relax_failures.append(folder_path)
                    relax_rows.append((
                        [folder_name, "--", "--", "--", "--", f"FAILED ({e})"], 'red'))

        print_table(["Folder", "E before (eV)", "E after (eV)", "Vol change",
                     "Max disp (Ang)", "Steps"], relax_rows, f_out)
        if relax_failures:
            print_dual(color_text(
                f"[WARNING] MACE pre-relaxation failed for {len(relax_failures)} candidate(s) "
                "-- kept unrelaxed (original cast geometry).", 'yellow'), f_out)

    # ---- species resolution (post-casting): --molecular treats --species as molecule
    # names resolved into real elements only during generation, so args.species is not a
    # safe proxy for "which elements need a pseudopotential". Composition is identical
    # across every candidate in one run, so reading the first written candidate's actual
    # structure.fdf once is enough. Done after MACE pre-relaxation (species identity is
    # unaffected by relaxation, but reading after keeps this simply "read whatever the
    # final structure.fdf looks like").
    try:
        first_struct = structure_io.read_fdf(os.path.join(written[0], "structure.fdf"))
        real_species = sorted(structure_io.species_list(first_struct))
    except Exception as e:
        print_dual(color_text(
            f"[WARNING] Could not read species from '{written[0]}/structure.fdf' ({e}) -- "
            "falling back to --species as a best-effort approximation (may be wrong with "
            "--molecular).", 'yellow'), f_out)
        real_species = sorted(set(args.species))

    print_section("[4] PSEUDOPOTENTIALS", f_out)
    missing = []
    if args.pseudo_dir:
        found, missing = get_required_pseudos(real_species, args.pseudo_dir)
        print_dual(f"Source                 : {args.pseudo_dir}", f_out)
        print_table(["Species", "Status"], [
            ([sp, "MISSING" if sp in missing else "found"], 'yellow' if sp in missing else None)
            for sp in real_species
        ], f_out)
        if missing:
            print_dual(color_text(
                f"[WARNING] Missing pseudopotential(s) for: {', '.join(missing)} -- these "
                "will need to be added manually to every generated folder.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[OK] All required pseudopotentials found -- will be copied into every "
                "generated folder.", 'green'), f_out)
    else:
        print_dual("Not given (pass -p/--pseudo-dir -- a bundled bank or a folder path -- to "
                   "copy the required pseudopotential for every species into every generated "
                   "folder). Pseudopotentials will need to be added manually.", f_out)

    print_section("[5] CALC.FDF & PSEUDOPOTENTIAL GENERATION", f_out)
    folder_rows = []
    calc_failures = []
    notes_by_folder = {}
    # Outer capture_library_noise wraps the whole loop to catch genuine third-party
    # warnings.warn() calls (e.g. pymatgen/spglib symmetry analysis inside
    # generate_calculation) across every folder, deduped into one entry for [7] below.
    # It does NOT interfere with the per-folder redirect_stdout right below -- that
    # inner redirect wins for stdout while this call captures warnings independently.
    with capture_library_noise(library_warnings, "calc.fdf generation (pymatgen/spglib)"):
        for folder_path in written:
            folder_name = os.path.basename(folder_path)

            if args.calc:
                shutil.copy2(args.calc, os.path.join(folder_path, "calc.fdf"))
                calc_ok, calc_note = True, "custom template"
            else:
                buf = io.StringIO()
                cwd = os.getcwd()
                try:
                    os.chdir(folder_path)
                    with contextlib.redirect_stdout(buf):
                        calc_ok = generate_calculation(
                            "structure.fdf", "relax", args.d3, args.spin_polarized, "", None)
                finally:
                    os.chdir(cwd)
                captured = strip_ansi(buf.getvalue())
                if calc_ok:
                    calc_note = "auto-generated (relax)"
                    warn_lines = [l.strip() for l in captured.splitlines()
                                  if "[WARNING]" in l or "[NOTE]" in l]
                    if warn_lines:
                        calc_note += f" ({len(warn_lines)} note/warning(s))"
                        notes_by_folder[folder_name] = warn_lines
                else:
                    calc_note = "FAILED (see notes)"
                    calc_failures.append(folder_path)
                    notes_by_folder[folder_name] = [l.strip() for l in captured.splitlines() if l.strip()]

            # Pseudopotentials are always copied here (not via generate_calculation's own
            # copy_pseudopotentials, called above with pp_path="" to skip it -- that one
            # copies straight into os.getcwd() with no per-species-missing report) so both
            # calc.fdf modes get the same handling, matching every other multi-folder prep
            # tool (strain.py/elastic_inputs.py/convergence.py).
            for sym in real_species:
                copy_pseudo(args.pseudo_dir, sym, folder_path)
            if args.pseudo_dir:
                pseudo_note = "copied" if not missing else f"missing: {', '.join(missing)}"
            else:
                pseudo_note = "not requested"

            row_color = 'yellow' if (not calc_ok or (args.pseudo_dir and missing)) else None
            folder_rows.append(([folder_name, calc_note, pseudo_note], row_color))

    print_table(["Folder", "Calc.fdf", "Pseudopotentials"], folder_rows, f_out)
    if notes_by_folder:
        print_dual("", f_out)
        print_dual(color_text("Notes from calc.fdf generation:", 'yellow'), f_out)
        for name, lines in notes_by_folder.items():
            for line in lines:
                print_dual(f"  [{name}] {line}", f_out)

    print_section("[6] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"Candidate folders written : {len(written)}", f_out)
    if failed:
        print_dual(color_text(
            f"[WARNING] {len(failed)} of {len(groups)} space group(s) failed to cast: "
            f"{', '.join(str(g) for g in failed)} (see [2] above).", 'yellow'), f_out)
    if relax_failures:
        print_dual(color_text(
            f"[WARNING] MACE pre-relaxation failed for {len(relax_failures)} candidate "
            "folder(s) (see [3] above).", 'yellow'), f_out)
    if calc_failures:
        print_dual(color_text(
            f"[WARNING] calc.fdf generation failed for {len(calc_failures)} candidate "
            "folder(s) (see notes above).", 'yellow'), f_out)
    if report_path:
        print_dual(f"Report                     : {report_path}", f_out)
    print_dual(f"\n{color_text('Next:', 'cyan')} run SIESTA in each folder, then "
              f"stb-xrdrank --input-dir {args.output_dir} --experimental <file>", f_out)

    print_section("[7] LIBRARY WARNINGS", f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (MACE/torch/pymatgen/spglib/"
            "stb-crystalcast) during this run -- collected here instead of interleaved "
            "with the report above; harmless in almost every case (import-time notices, "
            "deprecation-style warnings), but worth a glance.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every report section above has already printed, so a
    # blocking GUI window never delays or hides them (same convention as
    # inputfile.py/crystalcast.py's own --view).
    if args.view:
        from stb.core.ase_view import view_structure_interactive
        atoms_list = []
        for folder_path in written:
            try:
                struct = structure_io.read_fdf(os.path.join(folder_path, "structure.fdf"))
                atoms_list.append(AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(struct)))
            except Exception as e:
                print(color_text(f"[WARNING] Could not load '{folder_path}' for viewing: {e}",
                                  'yellow'))
        if atoms_list:
            view_structure_interactive(atoms_list)

    if failed or relax_failures or calc_failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
