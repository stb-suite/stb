#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Reads and tabulates SIESTA's own NATIVE Hirshfeld/Voronoi atomic-charge
sections from a finished .out log -- no new charge-partitioning physics,
just formatting what SIESTA already computed and printed when
'Charge.Hirshfeld end'/'Charge.Voronoi end' were set in calc.fdf. Simple
Hirshfeld does not need reimplementing anywhere in this suite: it has been
validated (see new_functions/PROPOSAL_hirshfeld-I_implementation.md) to
agree with SIESTA's own native output to <0.02 e, so this tool exists
purely to spare the user from grepping the .out by hand -- the same
"read+tabulate what's already there" role stb-status plays for run-type/
SCF/force/energy. For a genuinely new (iterative Hirshfeld-I) partitioning
method, see stb-hirshfeldPrep/Ions/Analysis instead.
"""

VERSION = "1.0.0"

import os
import sys
import glob
import statistics
import argparse
from datetime import datetime

from stb.core import siesta_log, structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table

REPORT_FILE = "stb_nativecharges_report.txt"


def _by_species(rows):
    """{species_symbol: [charge, ...]} from a get_hirshfeld_charges/
    get_voronoi_charges-shaped row list -- same aggregation shape as
    stb-bader's own by_species, reused here for the per-species summary
    table below."""
    grouped = {}
    for row in rows:
        grouped.setdefault(row['sym'], []).append(row['charge'])
    return grouped


def _print_population_table(title, rows, f_out):
    print_section(title, f_out)
    headers = ["Idx", "Elem", "Charge(e-)", "Population(e-)", "Sz(µB)"]
    table_rows = []
    for r in rows:
        sz_str = f"{r['sz']:+.4f}" if r['sz'] is not None else "N/A"
        table_rows.append(([str(r['id']), r['sym'], f"{r['charge']:+.4f}",
                             f"{r['population']:.4f}", sz_str], None))
    print_table(headers, table_rows, f_out)

    grouped = _by_species(rows)
    species_rows = []
    for sym in sorted(grouped):
        values = grouped[sym]
        mean = statistics.mean(values)
        std = statistics.pstdev(values) if len(values) > 1 else 0.0
        species_rows.append(([sym, str(len(values)), f"{mean:+.4f}", f"{std:.4f}"], None))
    print_dual("", f_out)
    print_table(["Elem", "N", "Mean(e-)", "Std(e-)"], species_rows, f_out)


def inspect_folder(path, explicit_label, f_out=None):
    """Reads one folder's .out and prints/returns its native Hirshfeld/
    Voronoi tables. Returns {'path', 'label', 'hirshfeld': list[dict] |
    None, 'voronoi': list[dict] | None} -- both fields degrade to None
    (never raise) exactly like every core.siesta_log parser, so a --path
    glob batch run can always finish."""
    label = structure_io.find_fdf_system_label(path, explicit_label)
    out_file = siesta_log.find_out_file(path, label)

    result = {"path": path, "label": label, "hirshfeld": None, "voronoi": None}

    if label:
        print_dual(f"SystemLabel : {label}", f_out)
    else:
        print_dual(color_text(
            "[WARNING] Could not determine SystemLabel (pass --label, or make sure "
            "exactly one .fdf with a SystemLabel line is present).", 'yellow'), f_out)

    if not out_file:
        print_dual(color_text(
            "[WARNING] No .out log found -- can't report native charges.", 'yellow'), f_out)
        return result

    print_dual(f".out file   : {out_file}", f_out)
    hirshfeld = siesta_log.get_hirshfeld_charges(out_file)
    voronoi = siesta_log.get_voronoi_charges(out_file)
    result.update(hirshfeld=hirshfeld, voronoi=voronoi)

    if hirshfeld:
        _print_population_table("HIRSHFELD POPULATIONS", hirshfeld, f_out)
    if voronoi:
        _print_population_table("VORONOI POPULATIONS", voronoi, f_out)
    if not hirshfeld and not voronoi:
        print_dual(color_text(
            "[WARNING] Neither a native 'Hirshfeld Atomic Populations:' nor 'Voronoi Atomic "
            "Populations:' block was found in this .out -- add 'Charge.Hirshfeld end' and/or "
            "'Charge.Voronoi end' to calc.fdf and re-run SIESTA to get this tool's output.",
            'yellow'), f_out)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Reads and tabulates SIESTA's own native Hirshfeld/Voronoi atomic-charge "
                     "sections from a finished .out log (Charge.Hirshfeld end / Charge.Voronoi "
                     "end). No new charge-partitioning method -- just formats what SIESTA "
                     "already wrote. --path can be a glob pattern to scan several folders at "
                     "once (e.g. 'strain_*')."
    )
    parser.add_argument(
        "--path", default=".",
        help="Directory to inspect, or a glob pattern (e.g. 'strain_*') matching several "
             "directories (default: current directory)."
    )
    parser.add_argument(
        "-l", "--label", default=None,
        help="SystemLabel override, applied to every matched directory. Auto-detected per "
             "directory from its .fdf file(s) if not given."
    )
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nativecharges {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Native Charges",
            "Read and tabulate SIESTA's own native Hirshfeld/Voronoi output",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    is_glob = any(c in args.path for c in '*?[')
    matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p)) if is_glob else [args.path]
    if not matches:
        parser.error(f"No directory matched '{args.path}'.")
    if not is_glob and not os.path.isdir(args.path):
        parser.error(f"'{args.path}' is not a directory.")
    batch_mode = len(matches) > 1

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-NATIVECHARGES REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if batch_mode:
        print_dual(f"Target pattern : '{args.path}' ({len(matches)} director(y/ies) matched)", f_out)
    else:
        print_dual(f"Target directory: {matches[0]}", f_out)
    if report_path:
        print_dual(f"Report file    : {report_path}", f_out)

    results = []
    for d in matches:
        print_section(f"[1] {d}" if batch_mode else "[1] NATIVE CHARGES", f_out)
        results.append(inspect_folder(d, args.label, f_out))

    print_section("[2] SUMMARY", f_out)
    n_hirshfeld = sum(1 for r in results if r['hirshfeld'])
    n_voronoi = sum(1 for r in results if r['voronoi'])
    print_dual(f"Folders inspected           : {len(results)}", f_out)
    print_dual(f"Folders with Hirshfeld data  : {n_hirshfeld}", f_out)
    print_dual(f"Folders with Voronoi data    : {n_voronoi}", f_out)
    if report_path:
        print_dual(f"Report                       : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Native charges report complete.\n", 'bold'))


if __name__ == "__main__":
    main()
