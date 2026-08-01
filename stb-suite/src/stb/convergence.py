#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 1 of 2: writes one ready-to-run SIESTA folder (structure + calc.fdf
+ config_extra.fdf + pseudopotentials) per swept value of one or more
convergence-test parameters (Mesh.CutOff, PAO.EnergyShift, k-grid density),
one subfolder per parameter under -o/--output-dir. Doesn't run SIESTA itself
-- run each folder yourself, then use stb-convergenceAnalysis.

Every generated folder fully relaxes the structure (both atomic positions
AND the cell, MD.VariableCell true) at each swept value, rather than running
a fixed-geometry single-point SCF. This is deliberate, and specific to why a
Mesh.CutOff/PAO.EnergyShift/k-grid convergence test matters for SIESTA in
particular: SIESTA evaluates the Hartree/XC potential on a discrete
real-space grid (set by Mesh.CutOff), so forces and stresses carry a small,
grid-position-dependent numerical noise (the "eggbox effect") that plane
-wave codes don't have. At a coarse Mesh.CutOff (or an incomplete PAO basis,
or too sparse a k-mesh) this noise is large enough to bias where a variable
-cell relaxation actually lands -- the SCF total energy at some ARBITRARY
fixed geometry can look fine even while the relaxed lattice parameters/bond
lengths are still drifting with the parameter. Relaxing the full structure
at each swept value is what actually answers "is this setting good enough
to trust the geometry SIESTA gives me", which is normally the point of
running the calculation in the first place -- not just "is the energy at
one arbitrary geometry converged".

'%include config_extra.fdf' is PREPENDED at the very top of the generated
calc.fdf (SIESTA's fdf reader is first-occurrence-wins for duplicate
labels), forcing CG relaxation (MD.NumCGsteps/MD.Steps = --relax-steps,
MD.VariableCell true) regardless of what --calc's own template says. The
swept parameter's own value is ALSO written into config_extra.fdf (not
regex-substituted into --calc's own body) for the exact same reason -- it's
forced active via the same first-occurrence-wins mechanism, so it doesn't
need to already exist in --calc at all.

Output layout: <output-dir>/<parameter>/convergence_<parameter>_<value>/ --
one subfolder per swept parameter (so 'all' cleanly separates its 3
independent sweeps), same nested-plus-repeated-name convention
stb-elasticInputs already uses for its own <output-dir>/<direction>/
strain_<direction>_<pct>/ folders.

Per-parameter range flags (--meshcutoff-range/--energyshift-range/
--kgrid-range) let each selected parameter get its own custom sweep range
in the SAME invocation -- deliberately not a single shared --min/--max/
--step (an earlier version of this tool had that, and it could only ever
apply to exactly one parameter at a time, forcing multiple separate
invocations -- and multiple separate reports -- whenever more than one
parameter needed a non-default range). Any parameter without its own
-range flag uses its suggested default (see [2] of the report); every
selected parameter, customized or not, is always shown together in the
same single report.

KNOWN, DELIBERATE FOLLOW-UP: stb-convergenceAnalysis (Stage 2) does not yet
understand this per-parameter subfolder layout -- it still expects the old
flat <output-dir>/convergence_<parameter>_<value>/ layout (see its own
_FOLDER_RE + a non-recursive os.listdir(args.dir)). This is Stage 1 only;
Stage 2 gets its own follow-up update separately (same sequencing as this
session's phonons 4.4.1-then-4.4.2 work) -- flagged explicitly in [6] below
rather than left silently broken.
"""

VERSION = "2.2.0"

import os
import re
import sys
import shutil
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos, link_pseudo
from stb.core.structure_io import read_md_state, prepend_include

REPORT_FILE = "convergence_stage1.txt"
EXTRA_FDF_FILE = "config_extra.fdf"

# NOTE: substitute_numeric_tag/substitute_kgrid_tag/build_values/
# TAG_PATTERNS are imported directly by convergence_analysis.py
# (`from stb.convergence import substitute_numeric_tag, substitute_kgrid_tag`)
# -- kept for its own --apply feature (writing a converged value back into a
# real calc.fdf). Stage 1 itself no longer uses substitute_numeric_tag/
# substitute_kgrid_tag (see build_parameter_directive below) -- keep their
# names/signatures/behavior stable regardless.
TAG_PATTERNS = {
    "meshcutoff": re.compile(r'(Mesh\.CutOff\s+)([\d.]+)(\s*Ry)', re.IGNORECASE),
    "energyshift": re.compile(r'(PAO\.EnergyShift\s+)([\d.]+)(\s*Ry)', re.IGNORECASE),
}

# Suggested sweep ranges per parameter -- standard SIESTA convergence-test
# scales, not auto-derived from the structure itself. 'tighter' records
# which direction is MORE converged: meshcutoff is the opposite direction
# from energyshift/kgrid (a real, physical fact -- see
# convergence_analysis.py's own INCREASING_IMPROVES_ACCURACY, which encodes
# the same thing for its own plateau scan). Heavier elements/higher-accuracy
# work may need a wider meshcutoff range than the default; energyshift's
# textbook range is often quoted as a handful of specific values (0.001,
# 0.005, 0.01, 0.02, 0.05 Ry) rather than evenly spaced -- this tool's sweep
# is always linear (min/max/step), so the default below is an evenly spaced
# approximation of that same range, not the exact textbook set.
PARAMETER_DEFAULTS = {
    "meshcutoff":  {"min": 100.0, "max": 400.0, "step": 50.0, "unit": "Ry", "tighter": "larger"},
    "energyshift": {"min": 0.001, "max": 0.05, "step": 0.01, "unit": "Ry", "tighter": "smaller"},
    "kgrid":       {"min": 0.05, "max": 0.30, "step": 0.05, "unit": "1/Ang", "tighter": "smaller"},
}
PARAMETER_ORDER = ["meshcutoff", "energyshift", "kgrid"]


def substitute_numeric_tag(calc_text, parameter, value):
    """Replaces the numeric value of a 'Tag   NUMBER  Ry' line, keeping the
    tag name and unit untouched. Raises ValueError if the tag isn't found --
    a silently-unmodified calc.fdf would make every run identical.
    """
    pattern = TAG_PATTERNS[parameter]
    new_text, count = pattern.subn(rf'\g<1>{value:.4f}\g<3>', calc_text)
    if count == 0:
        raise ValueError(
            f"Could not find a '{parameter}' tag to substitute in the calc.fdf template."
        )
    return new_text


def substitute_kgrid_tag(calc_text, structure, density, vacuum_gap):
    """Computes a dimension-aware Monkhorst-Pack grid at the given density
    (same logic as stb-kgrid) and substitutes it into 'kgrid.MonkhorstPack'.
    """
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
    divisions = kspace.compute_monkhorts(
        structure.lattice[0], structure.lattice[1], structure.lattice[2], density, vacuum_axes
    )
    kgrid_str = f"kgrid.MonkhorstPack   [{divisions[0]}  {divisions[1]}  {divisions[2]}]"
    new_text, count = re.subn(r'kgrid\.MonkhorstPack\s+\[.*?\]', kgrid_str, calc_text, flags=re.IGNORECASE)
    if count == 0:
        raise ValueError(
            "Could not find a 'kgrid.MonkhorstPack' tag to substitute in the calc.fdf template."
        )
    return new_text, divisions


def build_values(min_v, max_v, step):
    n_steps = int(round((max_v - min_v) / step)) + 1
    return list(np.linspace(min_v, min_v + (n_steps - 1) * step, n_steps))


def build_parameter_directive(parameter, value, structure, vacuum_gap):
    """The swept parameter's own fdf directive line for THIS value, to be
    forced active via config_extra.fdf (not regex-substituted into --calc's
    body) -- so it doesn't need to already exist in --calc at all. Returns
    (directive_line, detail_string_for_the_report_table).
    """
    if parameter == "kgrid":
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
        divisions = kspace.compute_monkhorts(
            structure.lattice[0], structure.lattice[1], structure.lattice[2], value, vacuum_axes)
        line = f"kgrid.MonkhorstPack   [{divisions[0]}  {divisions[1]}  {divisions[2]}]"
        detail = f"grid {divisions[0]} {divisions[1]} {divisions[2]}"
    elif parameter == "meshcutoff":
        line = f"Mesh.CutOff           {value:.4f}  Ry"
        detail = f"{value:.4f} Ry"
    elif parameter == "energyshift":
        line = f"PAO.EnergyShift     {value:.4f} Ry"
        detail = f"{value:.4f} Ry"
    else:
        raise ValueError(f"Unknown parameter '{parameter}'.")
    return line, detail


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 of 2: writes one ready-to-run SIESTA folder per swept value of one "
                    "or more convergence-test parameters, one subfolder per parameter, fully "
                    "relaxing the structure (positions + cell) at each. Doesn't run SIESTA -- "
                    "run each folder yourself, then use stb-convergenceAnalysis.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo --parameter all\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo --parameter all "
               "--meshcutoff-range 150 350 50\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo --parameter kgrid "
               "--relax-steps 200 --save-report\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure.fdf, copied as 'structure.fdf' into every "
                             "generated folder (default: structure.fdf).")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                        help="calc.fdf template (default: calc.fdf) -- SCF/basis/MD settings, normally %%including "
                             "the structure file above by name. Copied into every generated "
                             "folder with one line added: '%%include config_extra.fdf', "
                             "prepended at the very top (see [3] of the report) -- this forces "
                             "full relaxation and the swept parameter's value, regardless of "
                             "what this template itself contains.")
    parser.add_argument("-p", "--pseudo-dir", default="",
                        help="Pseudopotentials source (optional): a bundled bank "
                             f"({', '.join(BANKS)}) or a folder path. If given, the required "
                             "pseudopotential for every species in the structure is linked into "
                             "every generated folder, so it's immediately ready for SIESTA.")
    parser.add_argument("--parameter", nargs='+', required=True,
                        help="Which tag(s) to sweep: meshcutoff, energyshift, kgrid, or 'all' "
                             "for all three (one independent sweep per parameter, each in its "
                             "own subfolder).")
    parser.add_argument("--meshcutoff-range", type=float, nargs=3, default=None,
                        metavar=("MIN", "MAX", "STEP"),
                        help="Custom Mesh.CutOff sweep range (Ry). Omit to use the suggested "
                             "default (see [2] of the report). Only meaningful if 'meshcutoff' "
                             "is in --parameter.")
    parser.add_argument("--energyshift-range", type=float, nargs=3, default=None,
                        metavar=("MIN", "MAX", "STEP"),
                        help="Custom PAO.EnergyShift sweep range (Ry) (see --meshcutoff-range). "
                             "Only meaningful if 'energyshift' is in --parameter.")
    parser.add_argument("--kgrid-range", type=float, nargs=3, default=None,
                        metavar=("MIN", "MAX", "STEP"),
                        help="Custom k-grid density sweep range (1/Ang) (see "
                             "--meshcutoff-range). Only meaningful if 'kgrid' is in --parameter.")
    parser.add_argument("--relax-steps", type=int, default=100,
                        help="Number of relaxation steps (MD.NumCGsteps/MD.Steps) forced in "
                             "every generated folder (default: 100). Both atomic positions and "
                             "the cell relax (MD.VariableCell true).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="--parameter kgrid only: vacuum-axis detection threshold in Angstrom, "
                             "same meaning as stb-kgrid's --vacuum-gap (default: 10.0).")
    parser.add_argument("-o", "--output-dir", type=str, default="convergence_runs",
                        help="Top-level directory (default: convergence_runs). Run folders go "
                             "into its own '<output-dir>/<parameter>/convergence_<parameter>_"
                             "<value>/' subfolder, one tree per swept parameter.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-convergence {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    # Expand 'all' / validate parameter names, de-duplicated, in a stable
    # (PARAMETER_ORDER) order regardless of how the user listed them.
    requested = []
    for p in args.parameter:
        if p == "all":
            for k in PARAMETER_ORDER:
                if k not in requested:
                    requested.append(k)
        elif p in PARAMETER_ORDER:
            if p not in requested:
                requested.append(p)
        else:
            parser.error(f"--parameter: unrecognized value '{p}' (choose from "
                         f"{', '.join(PARAMETER_ORDER)}, or 'all').")
    args.parameter = requested

    if args.relax_steps <= 0:
        parser.error("--relax-steps must be > 0.")

    range_flags = {
        "meshcutoff": args.meshcutoff_range,
        "energyshift": args.energyshift_range,
        "kgrid": args.kgrid_range,
    }
    for p, custom in range_flags.items():
        if custom is not None and p not in args.parameter:
            parser.error(f"--{p}-range was given but '{p}' is not in --parameter.")

    ranges = {}
    for p in args.parameter:
        custom = range_flags[p]
        if custom is not None:
            mn, mx, st = custom
            if st <= 0:
                parser.error(f"--{p}-range: step must be > 0.")
            if mx <= mn:
                parser.error(f"--{p}-range: max must be greater than min.")
            ranges[p] = {"min": mn, "max": mx, "step": st, "source": f"custom (--{p}-range)"}
        else:
            d = PARAMETER_DEFAULTS[p]
            ranges[p] = {"min": d["min"], "max": d["max"], "step": d["step"],
                         "source": "suggested default"}

    if not os.path.isfile(args.structure):
        sys.exit(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
    if not os.path.isfile(args.calc):
        sys.exit(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-CONVERGENCE STAGE 1 REPORT (CONVERGENCE TEST PREP) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc template     : {args.calc}", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir or '(not given)'}", f_out)
    print_dual(f"Parameter(s)      : {', '.join(args.parameter)}", f_out)
    print_dual(f"Relax steps       : {args.relax_steps}", f_out)
    print_dual(f"Output dir        : {args.output_dir}", f_out)
    print_dual(f"Save report       : {'yes' if args.save_report else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] INPUT STRUCTURE", f_out)
    try:
        structure = structure_io.read_fdf(args.structure)
    except (ValueError, FileNotFoundError) as e:
        fail(str(e))

    counts = structure_io.atom_counts(structure)
    species = sorted(sp for sp, n in counts.items() if n > 0)
    composition = ", ".join(f"{sp}{n}" for sp, n in counts.items())
    print_table(["Quantity", "Value"], [
        (["Composition", composition], None),
        (["Total atoms", str(sum(counts.values()))], None),
        (["Coordinate format", structure.coord_format], None),
    ], f_out)

    print_section("[2] PARAMETER SELECTION", f_out)
    rows = []
    for p in args.parameter:
        r = ranges[p]
        d = PARAMETER_DEFAULTS[p]
        rows.append(([p, r["source"],
                      f"{r['min']:g} to {r['max']:g}, step {r['step']:g} ({d['unit']})",
                      f"{d['tighter']} = more converged"], None))
    print_table(["Parameter", "Range source", "Sweep range", "Convergence direction"], rows, f_out)

    print_section("[3] RELAXATION & PARAMETER ENFORCEMENT", f_out)
    with open(args.calc) as f:
        original_calc_text = f.read()
    structure_basename = os.path.basename(args.structure)
    calc_basename = os.path.basename(args.calc)
    if "%include" not in original_calc_text or structure_basename not in original_calc_text:
        print_dual(color_text(
            f"[NOTE] '{calc_basename}' may not %include '{structure_basename}' -- add it if "
            "needed so SIESTA picks up the structure in each generated folder.", 'yellow'), f_out)

    before = read_md_state(original_calc_text)
    steps_label = f"{before['steps_key']}={before['steps']}" if before['steps_key'] else "(absent)"
    print_dual(f"Calc template (before forcing): MD.TypeOfRun={before['typeofrun'] or '(absent)'}  "
               f"Steps: {steps_label}  MD.VariableCell={before['variablecell'] or '(absent)'}", f_out)
    print_dual(f"Forced via '{EXTRA_FDF_FILE}' (prepended -- always wins) in every folder -- the "
               "full structure (positions + cell) relaxes at each swept value:", f_out)
    print_table(["Directive", "Forced value"], [
        (["MD.TypeOfRun", "CG"], None),
        (["MD.NumCGsteps / MD.Steps", str(args.relax_steps)], None),
        (["MD.VariableCell", "true (full relaxation: positions + cell)"], None),
        (["Swept parameter", f"{', '.join(args.parameter)} -- one value per folder, see [5]"], None),
    ], f_out)
    forced_calc_text = prepend_include(original_calc_text, EXTRA_FDF_FILE)

    print_section("[4] PSEUDOPOTENTIALS", f_out)
    if args.pseudo_dir:
        found, missing = get_required_pseudos(species, args.pseudo_dir)
        print_dual(f"Source            : {args.pseudo_dir}", f_out)
        print_table(["Species", "Status"], [
            ([sp, "MISSING" if sp in missing else "found"], 'yellow' if sp in missing else None)
            for sp in species
        ], f_out)
        if missing:
            print_dual(color_text(
                f"[WARNING] Missing pseudopotential(s) for: {', '.join(missing)} -- these will "
                "need to be added manually to every generated folder.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[OK] All required pseudopotentials found -- will be linked into every "
                "generated folder.", 'green'), f_out)
    else:
        print_dual("Not given (pass -p/--pseudo-dir -- a bundled bank or a folder path -- to "
                   "link the required pseudopotential for every species into every generated "
                   "folder). Pseudopotentials will need to be added manually.", f_out)

    print_section("[5] GENERATED CONVERGENCE FOLDERS", f_out)
    print_dual(f"Under '{args.output_dir}/<parameter>/':", f_out)
    folder_rows = []
    per_parameter_count = {}
    for p in args.parameter:
        r = ranges[p]
        values = build_values(r["min"], r["max"], r["step"])
        for value in values:
            folder_name = f"convergence_{p}_{value:.4f}"
            folder = os.path.join(args.output_dir, p, folder_name)
            os.makedirs(folder, exist_ok=True)

            try:
                directive_line, detail = build_parameter_directive(p, value, structure, args.vacuum_gap)
            except ValueError as e:
                fail(str(e))

            extra_fdf_text = (
                "# Auto-generated by stb-convergence -- forces this folder's own swept parameter\n"
                "# value plus full relaxation (positions + cell), regardless of --calc's own settings.\n"
                f"{directive_line}\n"
                "\n"
                "MD.TypeOfRun       CG\n"
                f"MD.NumCGsteps      {args.relax_steps}\n"
                f"MD.Steps           {args.relax_steps}\n"
                "MD.VariableCell    true\n"
            )

            shutil.copy(args.structure, os.path.join(folder, structure_basename))
            with open(os.path.join(folder, calc_basename), "w") as fh:
                fh.write(forced_calc_text)
            with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
                fh.write(extra_fdf_text)

            for sym in species:
                link_pseudo(args.pseudo_dir, sym, folder)

            folder_rows.append(([folder_name, p, detail,
                                f"{structure_basename}, {calc_basename}, {EXTRA_FDF_FILE}"], None))
            per_parameter_count[p] = per_parameter_count.get(p, 0) + 1

    print_table(["Folder", "Parameter", "Value", "Files"], folder_rows, f_out)

    print_section("[6] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"{len(folder_rows)} folder(s) written under '{args.output_dir}' "
               f"({len(per_parameter_count)} parameter(s)):", f_out)
    for p, n in per_parameter_count.items():
        print_dual(f"  {p}: {n} folder(s)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    print_dual(color_text(
        "\n[NOTE] stb-convergenceAnalysis (Stage 2) doesn't yet support this per-parameter "
        f"layout -- point --dir at one parameter's own subfolder (e.g. '{args.output_dir}/"
        f"{args.parameter[0]}'), not '{args.output_dir}' itself. A Stage 2 update is planned "
        "separately.", 'yellow'), f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{args.output_dir}/<parameter>/convergence_*/' folder.", f_out)
    print_dual("  2. Once they're done: stb-convergenceAnalysis --dir <one parameter's subfolder>", f_out)

    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
