#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.13.0"

import sys
import os
import argparse
from pymatgen.core.periodic_table import Element
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.passivation import passivate_dangling_bonds


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Caps dangling bonds on a cut surface with a passivating atom.", 'bold')}
Detects undercoordinated atoms (e.g. from stb-slab) via local coordination
geometry and places a passivant (default H) along each missing-bond
direction. Only single-missing-bond sites are auto-passivated -- sites
missing 2+ bonds are geometrically underdetermined from local coordination
alone and are reported instead of guessed.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f slab.fdf\n"
               "  %(prog)s -f slab.fdf --passivant F --bond-length 1.4 -o slab_F.fdf\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf), typically a cut slab.")
    parser.add_argument("--passivant", type=str, default="H",
                        help="Element to cap dangling bonds with (default: H).")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="Neighbor-search radius in Angstrom, used to determine each atom's "
                             "coordination number. Default: auto-detected as 1.25x the shortest "
                             "pairwise distance in the structure.")
    parser.add_argument("--bond-length", type=float, default=None,
                        help="Passivant bond length in Angstrom. Default: auto, per surface-"
                             "species/passivant pair, as the sum of their pymatgen atomic radii "
                             "(an approximation -- e.g. Si+H = 1.35 Ang vs. the experimental "
                             "Si-H ~1.48 Ang).")
    parser.add_argument("-o", "--output", type=str, default="passivated.fdf",
                        help="Output .fdf file name (default: passivated.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-passivate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    try:
        Element(args.passivant)
    except ValueError as e:
        parser.error(str(e))

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Passivate dangling bonds on a surface:", 'bold'))
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
    print(f"  {color_text('Input formula:', 'cyan')} {pmg.composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {len(pmg)}")
    print(f"  {color_text('Passivant:', 'cyan')} {args.passivant}")

    new_pmg, report = passivate_dangling_bonds(
        pmg, passivant=args.passivant, cutoff=args.cutoff, bond_length=args.bond_length)

    n_passivated = len(report["passivated"])
    n_unresolved = len(report["unresolved"])
    print(f"\n  {color_text('Dangling bonds found:', 'cyan')} {n_passivated + n_unresolved}")
    print(f"  {color_text('Auto-passivated:', 'cyan')} {n_passivated}")

    if n_unresolved:
        print(color_text(
            f"\n  Warning: {n_unresolved} atom(s) are missing 2+ bonds -- the missing directions "
            "are geometrically underdetermined from local coordination alone, so they were left "
            "unpassivated. Review these manually:", 'yellow'))
        for idx, symbol, pos, deficit in report["unresolved"]:
            print(f"    #{idx + 1:<4} {symbol:<3} deficit={deficit}  at {pos}")

    species_meta = structure_io.ensure_species_id(dict(structure.species_meta), args.passivant)

    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=species_meta)
    structure_io.write_fdf(new_structure, args.output)

    print(f"\n  {color_text('Output formula:', 'cyan')} {new_pmg.composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {len(new_pmg)}")
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
