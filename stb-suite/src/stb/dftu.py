#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import sys
import argparse
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.dftu_data import SHELL_NAMES, DEFAULT_SHELL, REFERENCE_U, ldau_proj_block


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Generates a ready-to-use %block LDAU.proj snippet for one or more species.", 'bold')}
You supply U (and optionally J) explicitly -- this tool never guesses or
estimates U for you. --list-reference prints a small table of literature GGA+U
values (citation included) as a convenience starting point; --suggest prints
that value for one element with a clear disclaimer. Neither is ever used
automatically -- --u must always be given to produce a block. For a
first-principles U computed from your own system via the Cococcioni & de
Gironcoli linear-response method, see the Workflow menu's "Hubbard U (Linear
Response)" entry (stb-hubbardu / stb-hubbarduAnalysis) instead.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --species Mn --u 3.9\n"
               "  %(prog)s --species Fe Co --u 5.3 3.32 --shell 3d 3d\n"
               "  %(prog)s --list-reference\n"
               "  %(prog)s --suggest Ni\n"
    )

    parser.add_argument("--species", type=str, nargs='+',
                         help="Species label(s) to add a DFT+U correction to.")
    parser.add_argument("--u", type=float, nargs='+',
                         help="Hubbard U (eV), one per --species, in the same order.")
    parser.add_argument("--j", type=float, nargs='+', default=None,
                         help="Exchange J (eV), one per --species (default: 0.0 for each).")
    parser.add_argument("--shell", type=str, nargs='+', default=None, choices=sorted(SHELL_NAMES),
                         help="Correlated shell per species (3d/4d/5d/4f/5f). Default: the standard "
                              "shell for each species (transition metals -> (n)d, lanthanides -> 4f, "
                              "actinides -> 5f).")
    parser.add_argument("-o", "--output", type=str, default=None,
                         help="Also save the block to this file (in addition to printing it).")
    parser.add_argument("--list-reference", action="store_true",
                         help="Print the literature GGA+U reference table (with citation) and exit.")
    parser.add_argument("--suggest", type=str, metavar="ELEMENT",
                         help="Print the literature reference U for one element (with a disclaimer) "
                              "and exit -- does not generate a block.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-dftu {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("DFT+U / Hubbard Block Generator:", 'bold'))
    print("-" * 60)

    if args.list_reference:
        print(color_text(
            "Reference GGA+U values (eV), oxides, Wang/Maxisch/Ceder, "
            "Phys. Rev. B 73, 195107 (2006), as tabulated by the Materials Project:", 'cyan'))
        for el, u in REFERENCE_U.items():
            print(f"  {el:<3} {u:.2f} eV")
        print(color_text(
            "\nThese are a starting-point sanity check, not a validation -- the actual U "
            "depends on your functional, pseudopotential, and basis set.", 'yellow'))
        return

    if args.suggest:
        el = args.suggest
        if el not in REFERENCE_U:
            print(color_text(f"Error: No reference value tabulated for '{el}'.", 'red'))
            sys.exit(1)
        print(color_text(
            f"[INFO] Literature reference (GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 "
            f"(2006)) for {el}: {REFERENCE_U[el]:.2f} eV -- a starting-point sanity check only, "
            "not a validation (functional/pseudopotential/basis-dependent). Use --u to set your "
            "own value explicitly.", 'cyan'))
        return

    if not args.species or not args.u:
        parser.error("--species and --u are required (unless using --list-reference or --suggest).")

    if len(args.u) != len(args.species):
        print(color_text(
            f"Error: --u has {len(args.u)} value(s) but --species has {len(args.species)} -- "
            "must match one-to-one.", 'red'))
        sys.exit(1)

    j_values = args.j if args.j is not None else [0.0] * len(args.species)
    if len(j_values) != len(args.species):
        print(color_text(
            f"Error: --j has {len(j_values)} value(s) but --species has {len(args.species)} -- "
            "must match one-to-one.", 'red'))
        sys.exit(1)

    shell_values = args.shell if args.shell is not None else [None] * len(args.species)
    if len(shell_values) != len(args.species):
        print(color_text(
            f"Error: --shell has {len(shell_values)} value(s) but --species has "
            f"{len(args.species)} -- must match one-to-one.", 'red'))
        sys.exit(1)

    entries = []
    for species, u, j, shell in zip(args.species, args.u, j_values, shell_values):
        shell_name = shell or DEFAULT_SHELL.get(species)
        if shell_name is None:
            print(color_text(
                f"Error: No default correlated shell known for '{species}' -- "
                f"pass --shell explicitly ({', '.join(sorted(SHELL_NAMES))}).", 'red'))
            sys.exit(1)
        n, l = SHELL_NAMES[shell_name]
        entries.append({"species": species, "n": n, "l": l, "u": u, "j": j})
        print(f"  {color_text(species, 'cyan')}: shell={shell_name} (n={n}, l={l}), "
              f"U={u:.3f} eV, J={j:.3f} eV")
        if species in REFERENCE_U:
            print(color_text(
                f"    [INFO] Literature reference for {species}: {REFERENCE_U[species]:.2f} eV "
                "(GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 (2006)) -- sanity check only.",
                'yellow'))

    block = ldau_proj_block(entries)
    print("\n" + block)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(block)
        print(f"{color_text('[Saved]', 'cyan')} {args.output}")

    print(color_text("\nCopy this block into your .fdf's DFT+U section.", 'green'))


if __name__ == "__main__":
    main()
