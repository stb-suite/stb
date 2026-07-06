#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.14.0"

import sys
import argparse
import difflib
from ase.build import molecule
from ase.collections import g2
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro


def print_available_names():
    names = sorted(g2.names)
    print(color_text(f"Available molecules ({len(names)}):", 'bold'))
    columns = 6
    width = max(len(n) for n in names) + 2
    for i in range(0, len(names), columns):
        row = names[i:i + columns]
        print("  " + "".join(f"{n:<{width}}" for n in row))


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Builds a reference molecule from ASE's G2 database.", 'bold')}
Places the molecule in a vacuum box, ready as an isolated-molecule
reference (e.g. for stb-cohesive) or a quick starting structure.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --name H2O\n"
               "  %(prog)s --name CH3OH --vacuum 12 -o methanol.fdf\n"
               "  %(prog)s --list\n"
    )

    parser.add_argument("--name", type=str, default=None,
                        help="Exact (case-sensitive) name from ASE's G2 database, "
                             "e.g. 'H2O', 'CO2', 'C6H6'. See --list for all names. "
                             "Required unless --list is given.")
    parser.add_argument("--list", dest="list_names", action="store_true",
                        help="Print all available molecule names and exit.")
    parser.add_argument("--vacuum", type=float, default=10.0,
                        help="Vacuum padding in Angstrom, added in every direction "
                             "around the molecule (default: 10.0).")
    parser.add_argument("-o", "--output", type=str, default="molecule.fdf",
                        help="Output .fdf file name (default: molecule.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-molecule {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.list_names:
        print_available_names()
        sys.exit(0)

    if args.name is None:
        parser.error("--name is required (unless --list is given).")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Build a reference molecule:", 'bold'))
    print("-" * 60)
    print(f"  {color_text('Name:', 'cyan')} {args.name}")

    try:
        atoms = molecule(args.name, vacuum=args.vacuum)
    except KeyError:
        case_insensitive_hit = next(
            (n for n in g2.names if n.lower() == args.name.lower()), None)
        msg = f"Error: '{args.name}' is not a known molecule name (names are case-sensitive)."
        if case_insensitive_hit:
            msg += f" Did you mean '{case_insensitive_hit}'?"
        else:
            suggestions = difflib.get_close_matches(args.name, g2.names, n=3)
            if suggestions:
                msg += f" Did you mean: {', '.join(suggestions)}?"
        msg += " Run with --list to see all available names."
        print(color_text(msg, 'red'))
        sys.exit(1)

    pmg = AseAtomsAdaptor.get_structure(atoms)
    print(f"  {color_text('Formula:', 'cyan')} {pmg.composition.reduced_formula}")
    print(f"  {color_text('Atoms:', 'cyan')} {len(pmg)}")
    print(f"  {color_text('Cell:', 'cyan')} {pmg.lattice.abc[0]:.2f} x {pmg.lattice.abc[1]:.2f} x "
          f"{pmg.lattice.abc[2]:.2f} Ang (vacuum={args.vacuum})")

    species_meta = {}
    for symbol in dict.fromkeys(site.specie.symbol for site in pmg):
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(pmg, species_meta=species_meta, coord_format="cartesian")
    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
