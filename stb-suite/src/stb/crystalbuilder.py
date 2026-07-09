#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import sys
import argparse
from pymatgen.core import Structure, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro


def parse_sites(raw_sites):
    """Parses [['Ni','0','0','0'], ['O','0.25','0.25','0.25']] (argparse's
    --site nargs=4, all strings) into [('Ni', [0.0, 0.0, 0.0]), ...].
    """
    sites = []
    for symbol, x, y, z in raw_sites:
        try:
            coords = [float(x), float(y), float(z)]
        except ValueError:
            raise ValueError(f"--site '{symbol} {x} {y} {z}': x/y/z must be numbers.")
        sites.append((symbol, coords))
    return sites


def resolve_spacegroup(spec):
    """Returns spec as an int if it looks like one, else as-is (a symbol string)."""
    try:
        return int(spec)
    except ValueError:
        return spec


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Builds a structure from a space group and Wyckoff positions.", 'bold')}
Give only the symmetrically-distinct sites (e.g. one per Wyckoff letter) --
pymatgen expands the rest via the space group's own symmetry operations.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0\n"
               "  %(prog)s --spacegroup 225 --a 3.52 --site Ni 0 0 0 -o fcc_ni.fdf\n"
               "  %(prog)s --spacegroup Fd-3m --a 8.396 \\\n"
               "      --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 \\\n"
               "      --site O 0.2548 0.2548 0.2548\n"
    )

    parser.add_argument("--spacegroup", type=str, required=True,
                        help="Space group as an international symbol (e.g. 'Fm-3m') or "
                             "number (e.g. '225').")
    parser.add_argument("--a", type=float, required=True, help="Lattice constant a (Ang).")
    parser.add_argument("--b", type=float, default=None, help="Lattice constant b (Ang). Default: same as --a.")
    parser.add_argument("--c", type=float, default=None, help="Lattice constant c (Ang). Default: same as --a.")
    parser.add_argument("--alpha", type=float, default=90.0, help="Lattice angle alpha, degrees (default: 90).")
    parser.add_argument("--beta", type=float, default=90.0, help="Lattice angle beta, degrees (default: 90).")
    parser.add_argument("--gamma", type=float, default=90.0, help="Lattice angle gamma, degrees (default: 90).")

    parser.add_argument("--site", dest="sites", action="append", nargs=4, required=True,
                        metavar=("SYMBOL", "X", "Y", "Z"),
                        help="One symmetrically-distinct Wyckoff site: element symbol and "
                             "fractional x y z. Repeat --site for each distinct site.")

    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision for the post-build verification step "
                             "(default: 1e-3, matches stb-symmetry/stb-unitcell).")
    parser.add_argument("-o", "--output", type=str, default="crystal.fdf",
                        help="Output .fdf file name (default: crystal.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-crystalbuilder {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Build a structure from a space group:", 'bold'))
    print("-" * 60)

    try:
        sites = parse_sites(args.sites)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    for symbol, _ in sites:
        try:
            Element(symbol)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

    b = args.b if args.b is not None else args.a
    c = args.c if args.c is not None else args.a
    lattice = Lattice.from_parameters(args.a, b, c, args.alpha, args.beta, args.gamma)

    spacegroup = resolve_spacegroup(args.spacegroup)
    print(f"  {color_text('Requested space group:', 'cyan')} {args.spacegroup}")
    print(f"  {color_text('Lattice:', 'cyan')} a={args.a} b={b} c={c} "
          f"alpha={args.alpha} beta={args.beta} gamma={args.gamma}")
    for symbol, coords in sites:
        print(f"  {color_text('Site:', 'cyan')} {symbol} at {coords}")

    species = [symbol for symbol, _ in sites]
    coords = [xyz for _, xyz in sites]

    try:
        structure = Structure.from_spacegroup(spacegroup, lattice, species, coords)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    min_dist = structure_io.min_pairwise_distance(structure)
    if min_dist is not None and min_dist < 0.5:
        print(color_text(
            f"  Warning: some atoms are unusually close ({min_dist:.3f} Ang) -- "
            "check for overlapping --site coordinates.", 'yellow'))

    sga = SpacegroupAnalyzer(structure, symprec=args.symprec)
    detected_number = sga.get_space_group_number()
    detected_symbol = sga.get_space_group_symbol()
    try:
        requested_number = SpaceGroup(str(spacegroup)).int_number if isinstance(spacegroup, str) \
            else SpaceGroup.from_int_number(spacegroup).int_number
    except ValueError:
        requested_number = None

    print(f"\n  {color_text('Output formula:', 'cyan')} {structure.composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {len(structure)}")
    print(f"  {color_text('Detected space group:', 'cyan')} {detected_symbol} (No. {detected_number})")

    if requested_number is not None and requested_number != detected_number:
        print(color_text(
            f"  Note: requested space group (No. {requested_number}) differs from the "
            f"detected one (No. {detected_number}) -- the given sites likely sit on a "
            "higher-symmetry special position than requested. Still valid, just be aware.",
            'yellow'))

    species_meta = {}
    for symbol in species:
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(structure, species_meta=species_meta)
    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
