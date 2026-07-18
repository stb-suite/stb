#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import glob
import argparse
from datetime import datetime

from stb.core import siesta_log, structure_io
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_status_report.txt"

# Extensions worth checking for presence, grouped the way this suite's own
# analysis tools consume them -- not exhaustive, just "would you be able to
# run stb-bands/-dos/-cube/-ani2traj on this folder right now".
KEY_FILES = {
    "geometry": [".XV", ".STRUCT_OUT"],
    "electronic": [".bands", ".DOS", "PDOS.xml"],
    "density/potential": [".RHO", ".VT", ".VH"],
    "AIMD trajectory": [".ANI"],
}


def inspect_folder(path, explicit_label, f_out=None):
    """Prints one folder's status block (via print_dual) and returns a dict
    summary: {'path', 'label', 'category', 'steps', 'scf_converged',
    'max_force', 'free_energy', 'final_temp'}. Every field degrades to None
    (or 'unknown'/0 where a dict-table column needs a stable type) instead
    of raising, so a batch run over many folders can always finish and
    report what it found per folder rather than aborting on the first one
    missing a .out/.fdf.
    """
    label = structure_io.find_fdf_system_label(path, explicit_label)
    out_file = siesta_log.find_out_file(path, label)

    result = {
        "path": path, "label": label, "category": "unknown", "steps": 0,
        "scf_converged": None, "max_force": None, "free_energy": None, "final_temp": None,
    }

    if label:
        print_dual(f"SystemLabel       : {label}", f_out)
    else:
        print_dual(color_text(
            "[WARNING] Could not determine SystemLabel (pass --label, or make sure "
            "exactly one .fdf with a SystemLabel line is present).", 'yellow'), f_out)

    if not out_file:
        print_dual(color_text(
            "[WARNING] No .out log found -- can't report run type/convergence/energy.",
            'yellow'), f_out)
    else:
        dynamics = siesta_log.get_dynamics_type(out_file)
        category = siesta_log.categorize_dynamics(dynamics)
        steps = siesta_log.count_geometry_steps(out_file)
        scf_ok, _ = siesta_log.get_scf_convergence(out_file)
        max_force = siesta_log.get_max_force(out_file)
        free_energy = siesta_log.get_free_energy(out_file)
        result.update(category=category, steps=steps, scf_converged=scf_ok,
                       max_force=max_force, free_energy=free_energy)

        print_dual(f"Run type          : {category}" + (f" ({dynamics})" if dynamics else ""), f_out)
        if category in ("relaxation", "aimd"):
            print_dual(f"Steps             : {steps}", f_out)
        if scf_ok is not None:
            print_dual(f"SCF converged     : {'yes' if scf_ok else 'no'}", f_out)
        if max_force is not None:
            print_dual(f"Max force (eV/Ang): {max_force:.6f}", f_out)
        if free_energy is not None:
            print_dual(f"Free energy (eV)  : {free_energy:.6f}", f_out)

        if category == "aimd":
            md_steps = siesta_log.get_md_trajectory(out_file)
            if md_steps and md_steps[-1]["Temp_ion"] is not None:
                result["final_temp"] = md_steps[-1]["Temp_ion"]
                print_dual(f"Final temperature : {result['final_temp']:.3f} K", f_out)

    print_dual("Files present:", f_out)
    for group_name, exts in KEY_FILES.items():
        present = [ext for ext in exts if glob.glob(os.path.join(path, f"*{ext}"))]
        marker = color_text("yes", 'green') if present else color_text("no", 'yellow')
        detail = f" ({', '.join(present)})" if present else ""
        print_dual(f"    {group_name:<20}: {marker}{detail}", f_out)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Quick status summary of a SIESTA calc folder: run type "
                     "(single-point/relaxation/AIMD), SCF convergence, max force, final "
                     "energy, and which key output files (.XV, .bands, .DOS, .RHO, .ANI...) "
                     "are present. --path can be a glob pattern to scan several folders at "
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
    parser.add_argument("-v", "--version", action="version", version=f"stb-status {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Status Checker",
            "What's in this calc folder, and did it finish OK?",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STATUS:", 'bold'))
    print("-" * 60)

    is_glob = any(c in args.path for c in '*?[')
    matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p)) if is_glob else [args.path]
    if not matches:
        parser.error(f"No directory matched '{args.path}'.")
    if not is_glob and not os.path.isdir(args.path):
        parser.error(f"'{args.path}' is not a directory.")
    batch_mode = len(matches) > 1

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-STATUS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if batch_mode:
        print_dual(f"Target pattern    : '{args.path}' ({len(matches)} director(y/ies) matched)", f_out)
    else:
        print_dual(f"Target directory  : {matches[0]}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] STATUS", f_out)
    results = []
    for d in matches:
        if batch_mode:
            print_dual(f"\n--- {d} ---", f_out)
        results.append(inspect_folder(d, args.label, f_out))

    print_section("[2] SUMMARY & FILES", f_out)
    if batch_mode:
        print_dual(f"{'Folder':<24}{'Type':<14}{'Steps':<8}{'SCF':<6}{'Energy (eV)'}", f_out)
        for r in results:
            scf_str = "yes" if r["scf_converged"] else ("no" if r["scf_converged"] is False else "?")
            energy_str = f"{r['free_energy']:.4f}" if r["free_energy"] is not None else "N/A"
            folder_name = os.path.basename(os.path.normpath(r["path"])) or r["path"]
            print_dual(f"{folder_name:<24}{r['category']:<14}{r['steps']:<8}{scf_str:<6}{energy_str}", f_out)
    print_dual(f"Folders inspected : {len(results)}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Status check complete.\n", 'bold'))


if __name__ == "__main__":
    main()
