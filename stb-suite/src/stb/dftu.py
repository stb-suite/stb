#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import os
import sys
import argparse
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.dftu_data import SHELL_NAMES, DEFAULT_SHELL, REFERENCE_U, ldau_proj_block
from stb.core.structure_io import read_fdf, species_list
from stb.core import citations

REPORT_FILE = "stb_dftu_report.txt"
BIB_FILE = "references.bib"

# The literature GGA+U reference table this tool tabulates (see
# core/dftu_data.py::REFERENCE_U) -- cited whenever it's actually shown or
# used (--list-reference, --suggest, --use-reference, or the informational
# sanity-check note on an explicitly-given --u that happens to match a
# tabulated species). Kept local since this is currently the only consumer.
_BIB_WANG_MAXISCH_CEDER = ("WangMaxischCeder2006", """@article{WangMaxischCeder2006,
  author  = {Wang, Lei and Maxisch, Thomas and Ceder, Gerbrand},
  title   = {Oxidation energies of transition metal oxides within the {GGA+U} framework},
  journal = {Physical Review B},
  year    = {2006},
  volume  = {73},
  number  = {19},
  pages   = {195107},
  doi     = {10.1103/PhysRevB.73.195107}
}""")


def _fail(message, f_out, report_path):
    """Prints a red [ERROR] line, closes the report file if one is open,
    and exits with status 1 -- the single error exit point every check
    below funnels through, so none of them has to repeat the close/exit
    boilerplate."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Generates a ready-to-use %block LDAU.proj snippet for one or more species.", 'bold')}
You supply U (and optionally J) explicitly -- this tool never guesses or
estimates U for you. --list-reference prints a small table of literature GGA+U
values (citation included) as a convenience starting point; --suggest prints
that value for one element with a clear disclaimer. Neither is ever used
automatically -- --u must always be given to produce a block, unless you pass
--use-reference to explicitly opt in to auto-filling U from that same table.
--fdf reads the species list off an existing structure file instead of typing
--species by hand -- combine it with --use-reference to build a block straight
from a structure: only metal species with a tabulated value go into the block,
with a loud warning; non-metals are skipped silently and untabulated metals
are skipped with a warning. For a
first-principles U computed from your own system via the Cococcioni & de
Gironcoli linear-response method, see the Workflow menu's "Hubbard U (Linear
Response)" entry (stb-hubbardu / stb-hubbarduAnalysis) instead.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --species Mn --u 3.9\n"
               "  %(prog)s --species Fe Co --u 5.3 3.32 --shell 3d 3d\n"
               "  %(prog)s --fdf calc.fdf --use-reference\n"
               "  %(prog)s --fdf calc.fdf --species Mn --u 3.9\n"
               "  %(prog)s --list-reference\n"
               "  %(prog)s --suggest Ni\n"
    )

    parser.add_argument("--fdf", type=str, metavar="PATH",
                         help="Read the species list off this .fdf structure file instead of "
                              "typing --species by hand. If --species is also given, it filters "
                              "down to that subset (each one must be present in the file).")
    parser.add_argument("--species", type=str, nargs='+',
                         help="Species label(s) to add a DFT+U correction to. Required unless --fdf "
                              "is given.")
    parser.add_argument("--u", type=float, nargs='+',
                         help="Hubbard U (eV), one per species, in the same order. Required unless "
                              "--use-reference is given (and every species has a tabulated value).")
    parser.add_argument("--use-reference", action="store_true",
                         help="For any species without an explicit --u, auto-fill U from the "
                              "tabulated Materials Project reference table (see --list-reference). "
                              "A non-metal species (no default correlated shell, e.g. O, Si) is "
                              "silently skipped -- DFT+U doesn't apply to it. A metal species with "
                              "no tabulated value is skipped too, but with a warning. The block is "
                              "built from whatever is left; errors out only if nothing qualifies.")
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
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-dftu {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()
    if args.output:
        args.output = os.path.expanduser(args.output)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-DFTU REPORT =====", 'magenta'), f_out)

    # --- --list-reference: print the table and exit ---
    if args.list_reference:
        print_section("[0] REFERENCE TABLE", f_out)
        print_dual(
            "Reference GGA+U values (eV), oxides, Wang/Maxisch/Ceder, "
            "Phys. Rev. B 73, 195107 (2006), as tabulated by the Materials Project:", f_out)
        for el, u in REFERENCE_U.items():
            print_dual(f"  {el:<3} {u:.2f} eV", f_out)
        print_dual(color_text(
            "These are a starting-point sanity check, not a validation -- the actual U "
            "depends on your functional, pseudopotential, and basis set.", 'yellow'), f_out)

        print_section("[1] REFERENCES", f_out)
        citations.write_bib_file(BIB_FILE, [_BIB_WANG_MAXISCH_CEDER])
        print_dual(color_text(
            f"[OK] Citation for this reference table written to '{BIB_FILE}'.", 'green'), f_out)

        if report_path:
            print_dual(f"Report : {report_path}", f_out)
        if f_out:
            f_out.close()
        return

    # --- --suggest: print one element's reference value and exit ---
    if args.suggest:
        print_section("[0] SUGGESTED VALUE", f_out)
        el = args.suggest
        if el not in REFERENCE_U:
            _fail(f"No reference value tabulated for '{el}'.", f_out, report_path)
        print_dual(color_text(
            f"[INFO] Literature reference (GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 "
            f"(2006)) for {el}: {REFERENCE_U[el]:.2f} eV -- a starting-point sanity check only, "
            "not a validation (functional/pseudopotential/basis-dependent). Use --u to set your "
            "own value explicitly.", 'cyan'), f_out)

        print_section("[1] REFERENCES", f_out)
        citations.write_bib_file(BIB_FILE, [_BIB_WANG_MAXISCH_CEDER])
        print_dual(color_text(
            f"[OK] Citation for this reference value written to '{BIB_FILE}'.", 'green'), f_out)

        if report_path:
            print_dual(f"Report : {report_path}", f_out)
        if f_out:
            f_out.close()
        return

    # --- Main path: generate a %block LDAU.proj ---
    used_reference_table = False

    print_section("[0] RUN METADATA", f_out)
    if args.fdf:
        print_dual(f"Structure file : {args.fdf}", f_out)
    if args.species:
        print_dual(f"Species requested : {', '.join(args.species)}", f_out)
    print_dual(f"U source : {'--use-reference (auto-fill)' if args.use_reference else '--u (explicit)'}", f_out)

    print_section("[1] SPECIES RESOLUTION", f_out)
    if args.fdf:
        fdf_path = os.path.expanduser(args.fdf)
        try:
            fdf_structure = read_fdf(fdf_path)
        except (FileNotFoundError, ValueError) as e:
            _fail(f"Reading '{fdf_path}': {e}", f_out, report_path)
        fdf_species = species_list(fdf_structure)
        print_dual(f"Species found in '{fdf_path}' : {', '.join(fdf_species)}", f_out)
        if args.species:
            unknown = [s for s in args.species if s not in fdf_species]
            if unknown:
                _fail(
                    f"Species {unknown} passed via --species not found in '{fdf_path}' "
                    f"(found: {fdf_species}).", f_out, report_path)
            resolved_species = args.species
        else:
            resolved_species = fdf_species
    elif args.species:
        resolved_species = args.species
    else:
        _fail("--species or --fdf is required (unless using --list-reference or --suggest).", f_out, report_path)

    if args.u is not None:
        if len(args.u) != len(resolved_species):
            _fail(
                f"--u has {len(args.u)} value(s) but {len(resolved_species)} species are "
                "being processed -- must match one-to-one.", f_out, report_path)
        u_values = list(args.u)
    elif args.use_reference:
        used_reference_table = True
        untabulated_metals = [s for s in resolved_species if s in DEFAULT_SHELL and s not in REFERENCE_U]
        if untabulated_metals:
            print_dual(color_text(
                f"[WARNING] No tabulated reference U for metal species {untabulated_metals} -- "
                "skipped (pass --u explicitly for these if you need a DFT+U correction on them).",
                'yellow'), f_out)
        tabulated_species = [s for s in resolved_species if s in REFERENCE_U]
        if not tabulated_species:
            _fail(
                "None of the requested species have a tabulated reference U value -- "
                "nothing to generate. Pass --u explicitly, or restrict to species with a "
                f"tabulated value ({', '.join(sorted(REFERENCE_U))}).", f_out, report_path)
        resolved_species = tabulated_species
        u_values = [REFERENCE_U[s] for s in resolved_species]
        print_dual(color_text(
            f"[WARNING] Using tabulated Materials Project reference U values (GGA+U, oxides, "
            f"Wang/Maxisch/Ceder PRB 73, 195107 (2006)) for {resolved_species} -- a starting-point "
            "sanity check, not a validation. Verify against your own functional/pseudopotential/"
            "basis before production runs.", 'yellow'), f_out)
    else:
        _fail("--u is required unless --use-reference is given (with every species having "
              "a tabulated reference value).", f_out, report_path)

    j_values = args.j if args.j is not None else [0.0] * len(resolved_species)
    if len(j_values) != len(resolved_species):
        _fail(
            f"--j has {len(j_values)} value(s) but {len(resolved_species)} species are "
            "being processed -- must match one-to-one.", f_out, report_path)

    shell_values = args.shell if args.shell is not None else [None] * len(resolved_species)
    if len(shell_values) != len(resolved_species):
        _fail(
            f"--shell has {len(shell_values)} value(s) but {len(resolved_species)} species "
            "are being processed -- must match one-to-one.", f_out, report_path)

    print_section("[2] BLOCK GENERATION", f_out)
    entries = []
    for species, u, j, shell in zip(resolved_species, u_values, j_values, shell_values):
        shell_name = shell or DEFAULT_SHELL.get(species)
        if shell_name is None:
            _fail(
                f"No default correlated shell known for '{species}' -- "
                f"pass --shell explicitly ({', '.join(sorted(SHELL_NAMES))}).", f_out, report_path)
        n, l = SHELL_NAMES[shell_name]
        entries.append({"species": species, "n": n, "l": l, "u": u, "j": j})
        print_dual(f"  {species}: shell={shell_name} (n={n}, l={l}), U={u:.3f} eV, J={j:.3f} eV", f_out)
        if species in REFERENCE_U:
            used_reference_table = True
            print_dual(color_text(
                f"    [INFO] Literature reference for {species}: {REFERENCE_U[species]:.2f} eV "
                "(GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 (2006)) -- sanity check only.",
                'yellow'), f_out)

    block = ldau_proj_block(entries)
    print_dual("\n" + block, f_out)

    if args.output:
        print_section("[3] WRITING OUTPUT FILE", f_out)
        with open(args.output, 'w') as f:
            f.write(block)
        print_dual(color_text(f"[OK] Block saved to '{args.output}'.", 'green'), f_out)

    print_section("[4] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if used_reference_table:
        bib_entries.append(_BIB_WANG_MAXISCH_CEDER)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Species        : {', '.join(resolved_species)}", f_out)
    print_dual(f"Output file    : {args.output if args.output else '(not saved -- pass -o/--output)'}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    print_dual(color_text("\nCopy the block above into your .fdf's DFT+U section.", 'green'), f_out)

    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
