#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import re
import sys
import argparse
from stb.core import structure_io, siesta_log
from stb.core.cli import COLORS, color_text, show_intro

_FOLDER_RE = re.compile(r'^convergence_([a-zA-Z]+)_([\d.]+)$')


def parse_convergence_folder_name(folder_name):
    """Parses 'convergence_meshcutoff_250.0000' into ('meshcutoff', 250.0),
    matching the naming convention stb-convergence writes.
    """
    match = _FOLDER_RE.match(folder_name)
    if not match:
        return None, None
    return match.group(1), float(match.group(2))


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Analyzes a stb-convergence sweep and reports where the total energy converges.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --dir convergence_runs --file calc.out --tolerance 0.001\n"
    )

    parser.add_argument("--dir", type=str, default="convergence_runs",
                        help="Directory containing the 'convergence_*' run folders (default: convergence_runs).")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--tolerance", type=float, default=0.001,
                        help="Convergence tolerance in eV/atom (default: 0.001).")
    parser.add_argument("-o", "--output", type=str, default="convergence_curve.dat",
                        help="Output curve data file name (default: convergence_curve.dat).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-convergenceAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze a convergence-test sweep:", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.dir):
        print(color_text(f"Error: Directory '{args.dir}' not found.", 'red'))
        sys.exit(1)

    folders = sorted(
        f for f in os.listdir(args.dir)
        if os.path.isdir(os.path.join(args.dir, f)) and f.startswith("convergence_")
    )
    if not folders:
        print(color_text(f"Error: No 'convergence_*' folders found in '{args.dir}'.", 'red'))
        sys.exit(1)

    parameter_name = None
    data = []  # (value, energy, n_atoms)

    print(f"\n{color_text('READING FOLDERS:', 'bold')}")
    for folder in folders:
        param, value = parse_convergence_folder_name(folder)
        if param is None:
            continue

        folder_path = os.path.join(args.dir, folder)
        struct_path = os.path.join(folder_path, "structure.fdf")
        out_path = os.path.join(folder_path, args.file)

        if not os.path.exists(struct_path) or not os.path.exists(out_path):
            print(f"   -> {folder:<30} : {color_text('SKIP', 'yellow')} (missing structure.fdf or {args.file})")
            continue

        energy = siesta_log.get_free_energy(out_path)
        if energy is None:
            print(f"   -> {folder:<30} : {color_text('SKIP', 'yellow')} (could not parse energy)")
            continue

        structure = structure_io.read_fdf(struct_path)
        n_atoms = len(structure.atoms)

        parameter_name = parameter_name or param
        data.append((value, energy, n_atoms))
        print(f"   -> {folder:<30} : {color_text('OK', 'green')}")

    if not data:
        print(color_text("\nError: No valid data found in convergence folders.", 'red'))
        sys.exit(1)

    data.sort(key=lambda d: d[0])

    print(f"\n  {color_text('Parameter:', 'cyan')} {parameter_name}")
    print(f"  {color_text('Runs found:', 'cyan')} {len(data)}")

    header = f"{'Value':<12} {'Energy(eV)':<16} {'E/atom(eV)':<14} {'|Delta E/atom|':<16}"
    print("\n" + header)
    print("-" * len(header))

    rows = []
    converged_value = None
    prev_e_per_atom = None
    for value, energy, n_atoms in data:
        e_per_atom = energy / n_atoms
        delta = None if prev_e_per_atom is None else abs(e_per_atom - prev_e_per_atom)
        rows.append((value, energy, e_per_atom, delta))

        delta_str = f"{delta:.6f}" if delta is not None else "--"
        print(f"{value:<12.4f} {energy:<16.6f} {e_per_atom:<14.6f} {delta_str:<16}")

        if converged_value is None and delta is not None and delta < args.tolerance:
            converged_value = value

        prev_e_per_atom = e_per_atom

    print("-" * len(header))

    if converged_value is not None:
        print(f"\n{color_text('Converged at:', 'green')} {parameter_name} = {converged_value:.4f} "
              f"(|Delta E/atom| < {args.tolerance} eV)")
    else:
        print(f"\n{color_text('[WARNING]', 'yellow')} No point met the {args.tolerance} eV/atom tolerance; "
              "consider extending --max.")

    with open(args.output, 'w') as f:
        f.write(f"# Convergence curve | Parameter: {parameter_name}\n")
        f.write("# 1:Value 2:Energy(eV) 3:Energy_per_atom(eV) 4:Abs_Delta_E_per_atom(eV)\n")
        for value, energy, e_per_atom, delta in rows:
            f.write(f"{value:.6f}  {energy:.6f}  {e_per_atom:.6f}  {(delta or 0.0):.6f}\n")

    report_file = "convergence_report.txt"
    with open(report_file, 'w') as f:
        f.write("========================================\n")
        f.write("      CONVERGENCE TEST REPORT      \n")
        f.write("========================================\n")
        f.write(f"Parameter: {parameter_name}\n")
        f.write(f"Tolerance: {args.tolerance} eV/atom\n")
        if converged_value is not None:
            f.write(f"Converged at: {parameter_name} = {converged_value:.4f}\n")
        else:
            f.write("Not converged within the tested range.\n")

    print(f"\n{color_text('[Saved]', 'cyan')} Curve data -> {args.output}")
    print(f"{color_text('[Saved]', 'cyan')} Report    -> {report_file}")


if __name__ == "__main__":
    main()
