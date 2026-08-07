#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.structure_io import read_fdf, species_list
from stb.core.pseudopotentials import (
    BANKS, resolve_pseudo_source, detect_pseudo_bank, list_bank_elements,
    get_required_pseudos, copy_pseudo,
)
from stb.core import citations

REPORT_FILE = "stb_pseudo_report.txt"
BIB_FILE = "references.bib"


def _fail(message, f_out, report_path):
    """Prints a red [ERROR] line, closes the report file if one is open,
    and exits with status 1 -- same single error-exit funnel as
    dftu.py::_fail."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def resolve_species(structure_file, species_arg, f_out, report_path):
    """Returns the sorted, deduplicated list of element symbols required,
    either read off a .fdf structure or taken directly from --species."""
    if structure_file:
        path = os.path.expanduser(structure_file)
        try:
            structure = read_fdf(path)
        except (FileNotFoundError, ValueError) as e:
            _fail(f"Reading '{path}': {e}", f_out, report_path)
        found = species_list(structure)
        print_dual(f"Structure file : {path}", f_out)
        print_dual(f"Species found in structure : {', '.join(found)}", f_out)
        return sorted(set(found))
    print_dual(f"Species requested : {', '.join(species_arg)}", f_out)
    return sorted(set(species_arg))


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Resolves the pseudopotentials a structure needs and copies them into a run folder.", 'bold')}
Give a structure file (-f/--file, .fdf) or an explicit element list (--species),
plus a pseudopotential source (-p/--pp-path -- a bundled bank name or a
directory of your own .psf/.psml files). Reports which elements were found
and which are missing, and copies the resolved files into -o/--output
(default: the current directory) unless --dry-run is given. --fallback-dir
supplies a second source consulted only for elements missing from the
primary one. --list-elements browses a bundled bank without needing a
structure at all.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf -p dojo -o .\n"
               "  %(prog)s --species Fe O -p dojo --fallback-dir virtual_vault\n"
               "  %(prog)s -f structure.fdf -p dojo --dry-run\n"
               "  %(prog)s --list-elements dojo\n"
    )

    parser.add_argument("-f", "--file", type=str, dest="structure_file", metavar="PATH",
                         help="Path to a SIESTA structure file (.fdf) to read the required "
                              "species from. Required unless --species or --list-elements "
                              "is given.")
    parser.add_argument("--species", type=str, nargs='+',
                         help="Element symbol(s) to resolve directly, instead of reading a "
                              "structure file. Mutually exclusive with -f/--file.")
    parser.add_argument("-p", "--pp-path", type=str, dest="pp_path", metavar="BANK_OR_PATH",
                         help=f"Pseudopotential source: a bundled bank name ({', '.join(sorted(BANKS))}) "
                              "or a directory of your own .psf/.psml files. Required unless "
                              "--list-elements is given.")
    parser.add_argument("--fallback-dir", type=str, dest="fallback_dir", metavar="BANK_OR_PATH", default=None,
                         help="A second pseudopotential source, consulted only for elements "
                              "missing from -p/--pp-path -- fills gaps automatically instead of "
                              "just reporting them as missing.")
    parser.add_argument("--list-elements", type=str, dest="list_elements", metavar="BANK",
                         choices=sorted(BANKS),
                         help="Print the elements available in a bundled bank and exit -- no "
                              "structure or -p/--pp-path needed.")
    parser.add_argument("-o", "--output", type=str, default=".",
                         help="Directory to copy the resolved pseudopotentials into. Default: "
                              "the current directory.")
    parser.add_argument("--dry-run", action="store_true",
                         help="Only run the availability check and report -- do not copy any "
                              "files.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-pseudo {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                         help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-PSEUDO REPORT =====", 'magenta'), f_out)

    # --- --list-elements: browse a bank and exit, no structure needed ---
    if args.list_elements:
        print_section("[0] BANK ELEMENTS", f_out)
        bank = BANKS[args.list_elements]
        print_dual(f"Bank : {args.list_elements} -- {bank['description']}", f_out)
        elements = list_bank_elements(args.list_elements)
        print_dual(f"{len(elements)} element(s) available:", f_out)
        for i in range(0, len(elements), 10):
            print_dual("  " + " ".join(f"{e:<3}" for e in elements[i:i + 10]), f_out)
        if report_path:
            print_dual(f"Report : {report_path}", f_out)
        if f_out:
            f_out.close()
        return

    if not args.pp_path:
        _fail("-p/--pp-path is required (unless using --list-elements).", f_out, report_path)
    if not args.structure_file and not args.species:
        _fail("Either -f/--file or --species is required (unless using --list-elements).",
              f_out, report_path)
    if args.structure_file and args.species:
        _fail("Give either -f/--file or --species, not both.", f_out, report_path)

    ok = True
    try:
        print_section("[0] RUN METADATA", f_out)
        species = resolve_species(args.structure_file, args.species, f_out, report_path)
        print_dual(f"Pseudopotential source (primary)  : {args.pp_path}", f_out)
        if args.fallback_dir:
            print_dual(f"Pseudopotential source (fallback) : {args.fallback_dir}", f_out)
        print_dual(f"Output directory : {os.path.abspath(os.path.expanduser(args.output))}", f_out)
        print_dual(f"Dry run          : {'yes' if args.dry_run else 'no'}", f_out)

        print_section("[1] REQUIRED SPECIES", f_out)
        print_dual(f"{len(species)} unique element(s): {', '.join(species)}", f_out)

        print_section("[2] PSEUDOPOTENTIAL SOURCE(S)", f_out)
        try:
            pp_path_resolved = resolve_pseudo_source(args.pp_path)
        except ValueError as e:
            _fail(str(e), f_out, report_path)
        print_dual(f"Primary   : {args.pp_path} -> {pp_path_resolved}", f_out)
        primary_bank = detect_pseudo_bank(pp_path_resolved)

        fallback_resolved = None
        fallback_bank = None
        if args.fallback_dir:
            try:
                fallback_resolved = resolve_pseudo_source(args.fallback_dir)
            except ValueError as e:
                _fail(str(e), f_out, report_path)
            print_dual(f"Fallback  : {args.fallback_dir} -> {fallback_resolved}", f_out)
            fallback_bank = detect_pseudo_bank(fallback_resolved)

        print_section("[3] AVAILABILITY CHECK", f_out)
        resolved = {}  # element -> ("primary"/"fallback", file_path)
        found_primary, missing = get_required_pseudos(species, pp_path_resolved)
        for path in found_primary:
            el = os.path.splitext(os.path.basename(path))[0]
            resolved[el] = ("primary", path)

        if fallback_resolved and missing:
            found_fallback, missing = get_required_pseudos(missing, fallback_resolved)
            for path in found_fallback:
                el = os.path.splitext(os.path.basename(path))[0]
                resolved[el] = ("fallback", path)

        for el in species:
            if el in resolved:
                source, path = resolved[el]
                ext = os.path.splitext(path)[1].lstrip(".")
                print_dual(f"  {el:<3} FOUND    ({source}, .{ext})  {path}", f_out)
            else:
                print_dual(color_text(f"  {el:<3} MISSING", 'yellow'), f_out)

        print_dual(f"Resolved : {len(resolved)}/{len(species)}", f_out)
        if missing:
            print_dual(color_text(
                f"[WARNING] {len(missing)} element(s) not found in any given source: "
                f"{', '.join(missing)}.", 'yellow'), f_out)
            ok = False

        print_section("[4] ACTION TAKEN", f_out)
        if args.dry_run:
            print_dual("[INFO] --dry-run: no files copied.", f_out)
        else:
            output_dir = os.path.abspath(os.path.expanduser(args.output))
            os.makedirs(output_dir, exist_ok=True)
            for el, (source, path) in resolved.items():
                src_dir = pp_path_resolved if source == "primary" else fallback_resolved
                copy_pseudo(src_dir, el, output_dir)
                print_dual(f"[OK] Copied {el} pseudopotential to '{output_dir}'.", f_out)
            if not resolved:
                print_dual("[INFO] Nothing to copy -- no elements resolved.", f_out)

        print_section("[5] REFERENCES", f_out)
        bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
        if primary_bank:
            bib_entries.append(citations.PSEUDO_BANK_CITATIONS[primary_bank])
        if fallback_bank and fallback_bank != primary_bank:
            bib_entries.append(citations.PSEUDO_BANK_CITATIONS[fallback_bank])
        citations.write_bib_file(BIB_FILE, bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)

        print_section("[6] SUMMARY & NEXT STEPS", f_out)
        print_dual(f"Status     : {'OK' if ok else 'INCOMPLETE (missing elements)'}", f_out)
        print_dual(f"Resolved   : {len(resolved)}/{len(species)}", f_out)
        if not args.dry_run:
            print_dual(f"Copied to  : {os.path.abspath(os.path.expanduser(args.output))}", f_out)
        print_dual(f"References : {BIB_FILE}", f_out)
        if report_path:
            print_dual(f"Report     : {report_path}", f_out)
        print_dual(color_text(
            "\nPoint stb-inputfile's -p/--pp-path (or your calc.fdf's pseudopotential setup) "
            "at the output directory above.", 'green'), f_out)

    except Exception as e:
        print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
        ok = False

    if f_out:
        f_out.close()

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
