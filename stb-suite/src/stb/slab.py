#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.10.0"

import sys
import os
import argparse
from pymatgen.core.periodic_table import Element
from pymatgen.core.surface import SlabGenerator
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.passivation import passivate_dangling_bonds


def slab_metrics(slab):
    """(thickness, vacuum) in Angstrom. reorient_lattice=True (SlabGenerator's
    own default, kept here) puts the vacuum-normal on the c-axis, so this
    matches what kspace.detect_vacuum_axes will later find on the written file.
    """
    cart_z = slab.cart_coords[:, 2]
    thickness = float(cart_z.max() - cart_z.min())
    vacuum = float(slab.lattice.c - thickness)
    return thickness, vacuum


def sort_key(slab, symprec):
    """Non-polar and symmetric terminations first -- pymatgen's own get_slabs()
    order carries no physical meaning, so index 0 needs a meaningful default.
    """
    return (slab.is_polar(), not slab.is_symmetric(symprec=symprec))


def print_termination_table(slabs, symprec):
    header = (f"{'ID':<4} | {'Formula':<12} | {'Atoms':<6} | {'Thickness(A)':<13} | "
              f"{'Vacuum(A)':<10} | {'Polar':<6} | {'Symmetric':<9}")
    print(color_text(header, 'blue'))
    print("-" * len(header))
    for i, slab in enumerate(slabs):
        thickness, vacuum = slab_metrics(slab)
        polar = slab.is_polar()
        symmetric = slab.is_symmetric(symprec=symprec)
        row = (f"{i:<4} | {slab.composition.reduced_formula:<12} | {len(slab):<6} | "
               f"{thickness:<13.2f} | {vacuum:<10.2f} | "
               f"{'Yes' if polar else 'No':<6} | {'Yes' if symmetric else 'No':<9}")
        if not polar and symmetric:
            print(color_text(row, 'green'))
        else:
            print(row)


def prompt_termination(slabs, hkl, symprec):
    print("\n" + color_text(f"--- Terminations for hkl={hkl} ({len(slabs)} found) ---", 'bold'))
    print_termination_table(slabs, symprec)
    while True:
        user_input = input(color_text(
            f"\nSelect a termination ID (0 to {len(slabs) - 1}) or 'q' to quit: ", 'yellow'))
        if user_input.lower() == 'q':
            sys.exit(0)
        try:
            selected = int(user_input)
            if 0 <= selected < len(slabs):
                return selected
        except ValueError:
            pass


def print_slab_summary(idx, slab, hkl, symprec):
    thickness, vacuum = slab_metrics(slab)
    polar = slab.is_polar()
    symmetric = slab.is_symmetric(symprec=symprec)

    print("\n" + color_text("=" * 60, 'blue'))
    print(color_text(f"  Slab Summary: [Termination {idx}]", 'bold'))
    print(color_text("=" * 60, 'blue'))
    print(f"{'Formula':<22}: {color_text(slab.composition.reduced_formula, 'yellow')}")
    print(f"{'Total Atoms':<22}: {color_text(str(len(slab)), 'yellow')}")
    print(f"{'Miller Index (hkl)':<22}: {color_text(str(hkl), 'yellow')}")
    print(f"{'Slab Thickness':<22}: {color_text(f'{thickness:.2f} Ang', 'yellow')}")
    print(f"{'Vacuum Thickness':<22}: {color_text(f'{vacuum:.2f} Ang', 'yellow')}")
    print(f"{'Polar':<22}: {color_text('Yes' if polar else 'No', 'cyan')}")
    print(f"{'Symmetric':<22}: {color_text('Yes' if symmetric else 'No', 'cyan')}")
    print(color_text("-" * 60, 'blue'))


def write_slab(slab, species_meta, out_path):
    new_structure = structure_io.from_pymatgen(slab, species_meta=species_meta)
    structure_io.write_fdf(new_structure, out_path)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Cuts a Miller-index slab from a bulk SIESTA FDF structure.", 'bold')}
Adds vacuum along the surface normal; the vacuum axis lands on c.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f bulk.fdf --hkl 1 0 0\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 -i\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 --all -o surface.fdf\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 --passivate\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input bulk structure file (.fdf).")
    parser.add_argument("--hkl", type=int, nargs=3, required=True, metavar=("H", "K", "L"),
                        help="Miller index of the surface to cut, e.g. --hkl 1 0 0")
    parser.add_argument("--min-slab-size", type=float, default=10.0,
                        help="Minimum slab thickness in Angstrom (default: 10.0).")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0,
                        help="Minimum vacuum thickness in Angstrom (default: 15.0).")
    parser.add_argument("--lll-reduce", action="store_true",
                        help="LLL-reduce the slab lattice (pymatgen SlabGenerator option).")
    parser.add_argument("--center-slab", action="store_true",
                        help="Center the slab in the middle of the vacuum (pymatgen SlabGenerator option).")
    parser.add_argument("--primitive", action="store_true",
                        help="Reduce the bulk to its primitive cell before cutting the slab. "
                             "Off by default so the input cell is used as given.")
    parser.add_argument("--symmetrize", action="store_true",
                        help="Ask pymatgen to symmetrize polar/asymmetric terminations.")
    parser.add_argument("--symprec", type=float, default=0.1,
                        help="Symmetry tolerance used for the polar/symmetric diagnostics (default: 0.1).")

    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("-i", "--interactive", action="store_true",
                           help="Show all terminations found and prompt for one to keep.")
    selection.add_argument("--termination", type=int, default=None,
                           help="Directly select a termination index (0-based). "
                                "Default: index 0 after sorting non-polar/symmetric terminations first.")
    parser.add_argument("--all", action="store_true",
                        help="Write every termination found, instead of just one.")

    parser.add_argument("--passivate", action="store_true",
                        help="Cap dangling bonds on each written termination with a passivating "
                             "atom (same logic as stb-passivate). Only single-missing-bond sites "
                             "are auto-passivated; sites missing 2+ bonds are reported instead of "
                             "guessed.")
    parser.add_argument("--passivant", type=str, default="H",
                        help="Element to cap dangling bonds with, only with --passivate (default: H).")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="Neighbor-search radius in Angstrom for --passivate. "
                             "Default: auto-detected (see stb-passivate --help).")
    parser.add_argument("--bond-length", type=float, default=None,
                        help="Passivant bond length in Angstrom for --passivate. "
                             "Default: auto per species pair (see stb-passivate --help).")

    parser.add_argument("-o", "--output", type=str, default="slab.fdf",
                        help="Output .fdf file name (default: slab.fdf). When more than one "
                             "termination is written, it is suffixed '_term0', '_term1', ...")
    parser.add_argument("-v", "--version", action="version", version=f"stb-slab {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if not args.passivate:
        if args.passivant != "H":
            parser.error("--passivant is only valid with --passivate.")
        if args.cutoff is not None:
            parser.error("--cutoff is only valid with --passivate.")
        if args.bond_length is not None:
            parser.error("--bond-length is only valid with --passivate.")
    else:
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

    print("\n" + color_text("Build a Miller-index slab from a bulk structure:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    hkl = tuple(args.hkl)
    if hkl == (0, 0, 0):
        print(color_text("Error: Miller index (0, 0, 0) is not valid.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    pmg_structure = structure_io.to_pymatgen(structure)
    print(f"  {color_text('Input formula:', 'cyan')} {pmg_structure.composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {len(pmg_structure)}")
    print(f"  {color_text('Miller index (hkl):', 'cyan')} {hkl}")

    try:
        gen = SlabGenerator(
            pmg_structure, hkl, args.min_slab_size, args.min_vacuum_size,
            lll_reduce=args.lll_reduce, center_slab=args.center_slab, primitive=args.primitive,
        )
        slabs = gen.get_slabs(symmetrize=args.symmetrize)
    except Exception as e:
        print(color_text(f"Error: Could not generate slabs for hkl={hkl}: {e}", 'red'))
        sys.exit(1)

    if not slabs:
        print(color_text(f"Error: No slabs could be generated for Miller index {hkl}.", 'red'))
        sys.exit(1)

    slabs.sort(key=lambda s: sort_key(s, args.symprec))
    print(f"\n  {color_text('Terminations found:', 'cyan')} {len(slabs)}")

    if args.interactive:
        chosen_indices = [prompt_termination(slabs, hkl, args.symprec)]
    elif args.termination is not None:
        if not (0 <= args.termination < len(slabs)):
            print(color_text(
                f"Error: --termination {args.termination} is out of range "
                f"(0 to {len(slabs) - 1}).", 'red'))
            sys.exit(1)
        chosen_indices = [args.termination]
    elif args.all:
        chosen_indices = list(range(len(slabs)))
    else:
        chosen_indices = [0]

    base, ext = os.path.splitext(args.output)
    if not ext:
        ext = ".fdf"
    multi = len(chosen_indices) > 1

    for idx in chosen_indices:
        slab = slabs[idx]
        print_slab_summary(idx, slab, hkl, args.symprec)
        species_meta = structure.species_meta

        if args.passivate:
            slab, report = passivate_dangling_bonds(
                slab, passivant=args.passivant, cutoff=args.cutoff, bond_length=args.bond_length)
            n_passivated = len(report["passivated"])
            n_unresolved = len(report["unresolved"])
            print(f"  {color_text('Dangling bonds found:', 'cyan')} {n_passivated + n_unresolved}")
            print(f"  {color_text('Auto-passivated:', 'cyan')} {n_passivated} with {args.passivant}")
            if n_unresolved:
                print(color_text(
                    f"  Warning: {n_unresolved} atom(s) missing 2+ bonds -- left unpassivated "
                    "(geometrically underdetermined from local coordination alone). Review:", 'yellow'))
                for atom_idx, symbol, pos, deficit in report["unresolved"]:
                    print(f"    #{atom_idx + 1:<4} {symbol:<3} deficit={deficit}  at {pos}")
            species_meta = structure_io.ensure_species_id(dict(species_meta), args.passivant)

        out_path = f"{base}_term{idx}{ext}" if multi else args.output
        write_slab(slab, species_meta, out_path)
        print(f"  {color_text('Success:', 'green')} Slab written to '{color_text(out_path, 'bold')}'")

    print("\n" + color_text("Complete job!", 'bold'))


if __name__ == "__main__":
    main()
