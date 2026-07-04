#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import sys
import os
import argparse
import numpy as np
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro


def parse_matrix(values: list[float]) -> np.ndarray:
    """Builds a 3x3 supercell matrix from either 3 numbers (diagonal) or 9
    numbers (full matrix, row-major) -- same convention as phonopy's DIM tag.
    """
    if len(values) == 3:
        return np.diag(values).astype(float)
    if len(values) == 9:
        return np.array(values, dtype=float).reshape(3, 3)
    raise ValueError(f"Expected 3 (diagonal) or 9 (full matrix) numbers, got {len(values)}.")


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Builds a supercell from a SIESTA FDF structure file.", 'bold')}
Give -d 3 numbers for a diagonal supercell, or 9 for a full row-major matrix.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf -d 2 2 2\n"
               "  %(prog)s -f structure.fdf -d 2 0 0 0 2 0 0 0 1 -o big_cell.fdf\n"
    )

    parser.add_argument(
        "-f", "--file", dest="filename", type=str, required=True,
        help="Path to the input structure file (.fdf)."
    )
    parser.add_argument(
        "-d", "--dim", type=float, nargs='+', required=True,
        help="Supercell dimensions: 3 numbers for a diagonal supercell (e.g. 2 2 2),\n"
             "or 9 numbers for a full row-major 3x3 matrix (e.g. 2 0 0 0 2 0 0 0 1)."
    )
    parser.add_argument(
        "-o", "--output", type=str, default="supercell.fdf",
        help="Output .fdf file name (default: supercell.fdf)."
    )
    parser.add_argument("-v", "--version", action="version", version=f"stb-supercell {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Build a supercell from a structure file:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    try:
        matrix = parse_matrix(args.dim)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    rounded = np.round(matrix)
    if not np.allclose(matrix, rounded, atol=1e-6):
        print(color_text(f"Error: supercell matrix must have integer entries, got:\n{matrix}", 'red'))
        sys.exit(1)
    matrix = rounded.astype(int)

    det = int(round(np.linalg.det(matrix)))
    if det == 0:
        print(color_text("Error: supercell matrix is singular (determinant is zero).", 'red'))
        sys.exit(1)
    if det < 0:
        print(color_text(
            f"Warning: supercell matrix has a negative determinant ({det}) -- "
            "the resulting cell is mirrored (still geometrically valid).", 'yellow'
        ))

    try:
        structure = structure_io.read_fdf(args.filename)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    pmg_structure = structure_io.to_pymatgen(structure)
    print(f"  {color_text('Input formula:', 'cyan')} {pmg_structure.composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {len(pmg_structure)}")
    print(f"  {color_text('Supercell matrix:', 'cyan')}")
    for row in matrix:
        print(f"    {row[0]:3d} {row[1]:3d} {row[2]:3d}")
    print(f"  {color_text('Determinant (cell multiplication factor):', 'cyan')} {abs(det)}")

    supercell = pmg_structure.copy()
    supercell.make_supercell(matrix)

    print("\n" + color_text("--- Supercell Built ---", 'bold'))
    print(f"  {color_text('Output formula:', 'cyan')} {supercell.composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {len(supercell)}")
    print(f"  {color_text('Output volume:', 'cyan')} {supercell.volume:.4f} Ang^3")

    new_structure = structure_io.from_pymatgen(supercell, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Supercell written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
