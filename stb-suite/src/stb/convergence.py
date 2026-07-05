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
import shutil
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro

TAG_PATTERNS = {
    "meshcutoff": re.compile(r'(Mesh\.CutOff\s+)([\d.]+)(\s*Ry)', re.IGNORECASE),
    "energyshift": re.compile(r'(PAO\.EnergyShift\s+)([\d.]+)(\s*Ry)', re.IGNORECASE),
}


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


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Generates a sweep of calc.fdf variants for a convergence test.", 'bold')}
Sweeps Mesh.CutOff, PAO.EnergyShift, or k-grid density across a range,
one run folder per value.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -s structure.fdf -c calc.fdf --parameter meshcutoff --min 150 --max 350 --step 50\n"
               "  %(prog)s -s structure.fdf -c calc.fdf --parameter kgrid --min 0.10 --max 0.30 --step 0.05\n"
    )

    parser.add_argument("-s", "--structure", type=str, required=True,
                        help="Input structure.fdf, copied as 'structure.fdf' into every run folder.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                        help="Existing calc.fdf template (with %%include structure.fdf and everything "
                             "else already configured); exactly one tag is substituted per folder.")
    parser.add_argument("--parameter", choices=["meshcutoff", "energyshift", "kgrid"], required=True,
                        help="Which tag to sweep.")
    parser.add_argument("--min", type=float, required=True, help="Sweep range minimum.")
    parser.add_argument("--max", type=float, required=True, help="Sweep range maximum.")
    parser.add_argument("--step", type=float, required=True, help="Sweep step.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="--parameter kgrid only: vacuum-axis detection threshold in Angstrom, "
                             "same meaning as stb-kgrid's --vacuum-gap (default: 10.0).")
    parser.add_argument("-o", "--output-dir", type=str, default="convergence_runs",
                        help="Output directory for the run folders (default: convergence_runs).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-convergence {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.step <= 0:
        parser.error("--step must be > 0.")
    if args.max <= args.min:
        parser.error("--max must be greater than --min.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Generate a convergence-test sweep:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.structure):
        print(color_text(f"Error: File '{args.structure}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"Error: File '{args.calc}' not found.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.structure)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    with open(args.calc, 'r') as f:
        calc_template = f.read()

    values = build_values(args.min, args.max, args.step)
    print(f"  {color_text('Parameter:', 'cyan')} {args.parameter}")
    print(f"  {color_text('Sweep values:', 'cyan')} {[round(float(v), 4) for v in values]}")

    os.makedirs(args.output_dir, exist_ok=True)

    for value in values:
        folder = os.path.join(args.output_dir, f"convergence_{args.parameter}_{value:.4f}")
        os.makedirs(folder, exist_ok=True)

        try:
            if args.parameter == "kgrid":
                new_calc, divisions = substitute_kgrid_tag(calc_template, structure, value, args.vacuum_gap)
                detail = f"grid {divisions[0]} {divisions[1]} {divisions[2]}"
            else:
                new_calc = substitute_numeric_tag(calc_template, args.parameter, value)
                detail = f"{value:.4f} Ry"
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

        shutil.copy(args.structure, os.path.join(folder, "structure.fdf"))
        with open(os.path.join(folder, "calc.fdf"), 'w') as f:
            f.write(new_calc)

        print(f"  {color_text('[OK]', 'green')} {folder} ({detail})")

    print(f"\n{color_text('Success:', 'green')} Generated {len(values)} run(s) in "
          f"'{color_text(args.output_dir, 'bold')}'")
    print("Run the SIESTA calculation inside each 'convergence_...' folder, "
          "then use stb-convergenceAnalysis.")


if __name__ == "__main__":
    main()
