#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.10.3"

import sys
import os
import argparse
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.symmetry import reduce_to_unitcell


def generate_txt_report(args, sga, input_formula, input_atoms, output_formula, output_atoms):
    def row(label, value):
        return f"{label:<20}: {value}"

    lines = []
    lines.append("========================================")
    lines.append("      UNIT CELL REPORT")
    lines.append("========================================")
    lines.append(row("Input file", args.filename))
    lines.append(row("Input formula", input_formula))
    lines.append(row("Input atoms", input_atoms) + "\n")

    lines.append("== Symmetry ==")
    lines.append(row("Space group", f"{sga.get_space_group_symbol()} (No. {sga.get_space_group_number()})"))
    lines.append(row("Point group", sga.get_point_group_symbol()))
    lines.append(row("Crystal system", sga.get_crystal_system()))
    lines.append(row("symprec", args.symprec))
    lines.append(row("angle-tolerance", args.angle_tolerance) + "\n")

    lines.append("== Result ==")
    lines.append(row("Mode", args.mode))
    lines.append(row("Output formula", output_formula))
    lines.append(row("Output atoms", output_atoms))
    lines.append(row("Reduction factor", f"{input_atoms / output_atoms:.3g}x"))
    lines.append(row("Output structure", args.output))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Finds the primitive/conventional unit cell of a structure, or refines it.", 'bold')}
Uses pymatgen's symmetry analyzer (spglib) to detect the crystal's true
unit cell, which may be much smaller than the literal input cell (e.g. a
4-atom conventional FCC cell reduces to a 1-atom primitive cell).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf\n"
               "  %(prog)s -f structure.fdf --mode conventional -o conventional.fdf\n"
               "  %(prog)s -f structure.fdf --mode refined -o clean.fdf\n"
               "\n"
               "Notes:\n"
               "  - This operates on the literal 3D periodic cell as given. For "
               "slabs/2D materials/nanotubes with a vacuum direction, spot-check that "
               "the vacuum axis wasn't reoriented before using the output.\n"
               "  - The output's atom order and coordinate origin are NOT guaranteed to "
               "match the input, in any mode -- spglib rebuilds the cell from the detected "
               "symmetry operations from scratch, free to pick any symmetry-equivalent "
               "origin/ordering. For highly symmetric structures (e.g. a monoatomic FCC "
               "lattice) even tiny input noise can flip which equivalent origin it picks. "
               "This is expected: the crystal is the same, just relabeled.\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf).")
    parser.add_argument("--mode", choices=["primitive", "conventional", "refined"], default="primitive",
                        help="primitive: smallest possible cell (default).\n"
                             "conventional: standardized, usually larger, cell.\n"
                             "refined: conventional cell with atomic positions snapped to the\n"
                             "detected symmetry -- cleans up numerical noise from relaxations or\n"
                             "hand-built/CIF structures without changing which atoms are present.")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision (default: 1e-3, matches stb-symmetry).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry angle tolerance in degrees (default: 5.0, pymatgen's own default).")
    parser.add_argument("-o", "--output", type=str, default="unitcell.fdf",
                        help="Output .fdf file name (default: unitcell.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-unitcell {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Find the primitive/conventional unit cell:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    pmg = structure_io.to_pymatgen(structure)
    input_formula = pmg.composition.reduced_formula
    input_atoms = len(pmg)
    print(f"  {color_text('Input formula:', 'cyan')} {input_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {input_atoms}")

    try:
        new_pmg, sga = reduce_to_unitcell(pmg, args.mode, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
    except ValueError as e:
        print(color_text(f"Error: symmetry detection failed: {e}", 'red'))
        sys.exit(1)

    print(f"  {color_text('Space group:', 'cyan')} {sga.get_space_group_symbol()} (No. {sga.get_space_group_number()})")
    print(f"  {color_text('Point group:', 'cyan')} {sga.get_point_group_symbol()}")
    print(f"  {color_text('Crystal system:', 'cyan')} {sga.get_crystal_system()}")
    print(f"  {color_text('Mode:', 'cyan')} {args.mode}")

    output_formula = new_pmg.composition.reduced_formula
    output_atoms = len(new_pmg)

    print(f"\n  {color_text('Output formula:', 'cyan')} {output_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {output_atoms}")

    if args.mode != "refined" and output_atoms == input_atoms:
        print(color_text(
            f"  Note: input is already the {args.mode} cell at this symprec (no reduction).", 'yellow'))

    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output)

    report_file = "unitcell_report.txt"
    report_text = generate_txt_report(args, sga, input_formula, input_atoms, output_formula, output_atoms)
    with open(report_file, "w") as f:
        f.write(report_text)

    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")
    print(f"{color_text('[Saved]', 'cyan')} Report    -> {report_file}")


if __name__ == "__main__":
    main()
