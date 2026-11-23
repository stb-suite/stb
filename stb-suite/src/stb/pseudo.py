#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by the STB-SUITE team               #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION
from stb.core.cli import new_report_path

import os
import sys
import argparse
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.structure_io import read_fdf, species_list
from stb.core.pseudopotentials import (
    BANK_HELP, bank_names, bank_info, bank_manifest, verify_bank, resolve_pseudo_source,
    detect_pseudo_bank, list_bank_elements, get_required_pseudos, copy_pseudo, pseudo_valence,
    read_pseudo_header, warn_pseudo_level,
)
from stb.core import citations

REPORT_FILE = "stb_pseudo_report.txt"


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
structure at all, and --list-banks lists every bundled table.

PseudoDojo tables: -p dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent],
words in any order, a missing word takes its default (sr, pbe, standard):
  dojo                 scalar-relativistic PBE standard (the default table)
  dojo-pbesol          PBEsol;  dojo-pw (or dojo-lda): LDA (PW92)
  dojo-stringent       more semicore states (groups 13-16: the s,p shell below the d;
                       Hf/Ta: 4f) -- more electrons, more expensive
  dojo-fr[-pbesol]     full-relativistic, required for 'Spin spin-orbit'
  dojo-sr3plus         lanthanides Ce-Lu with 4f in the core (Ln3+ compounds only)
The table's functional must match the calculation's XC.Functional/XC.Authors.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf -p dojo -o .\n"
               "  %(prog)s --species Fe O -p dojo --fallback-dir virtual_vault\n"
               "  %(prog)s -f structure.fdf -p dojo --dry-run\n"
               "  %(prog)s --list-elements dojo\n"
               "  %(prog)s -f structure.fdf -p dojo-fr-pbesol --dry-run\n"
               "  %(prog)s -f structure.fdf -p dojo-sr3plus --fallback-dir dojo\n"
               "  %(prog)s --list-banks\n"
    )

    parser.add_argument("-f", "--file", type=str, dest="structure_file", metavar="PATH",
                         help="Path to a SIESTA structure file (.fdf) to read the required "
                              "species from. Required unless --species or --list-elements "
                              "is given.")
    parser.add_argument("--species", type=str, nargs='+',
                         help="Element symbol(s) to resolve directly, instead of reading a "
                              "structure file. Mutually exclusive with -f/--file.")
    parser.add_argument("-p", "--pp-path", type=str, dest="pp_path", metavar="BANK_OR_PATH",
                         help=f"Pseudopotential source: a bundled bank name ({BANK_HELP}) "
                              "or a directory of your own .psf/.psml files. Required unless "
                              "--list-elements is given.")
    parser.add_argument("--fallback-dir", type=str, dest="fallback_dir", metavar="BANK_OR_PATH", default=None,
                         help="A second pseudopotential source, consulted only for elements "
                              "missing from -p/--pp-path -- fills gaps automatically instead of "
                              "just reporting them as missing.")
    parser.add_argument("--list-elements", type=str, dest="list_elements", metavar="BANK",
                         help="Print the elements available in a bundled bank (any name -p "
                              "accepts, e.g. dojo-fr-pbesol) and exit -- no structure or "
                              "-p/--pp-path needed.")
    parser.add_argument("--list-banks", action="store_true",
                         help="Print every bundled bank (each PseudoDojo table with its relativity, "
                              "functional, accuracy and element count) and exit.")
    parser.add_argument("-o", "--output", type=str, default=".",
                         help="Directory to copy the resolved pseudopotentials into. Default: "
                              "the current directory.")
    parser.add_argument("--dry-run", action="store_true",
                         help="Only run the availability check and report -- do not copy any "
                              "files.")
    parser.add_argument("--save-report", action="store_true",
                         help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action="version", version=f"stb-pseudo {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                         help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by the STB-SUITE team"
        ])

    report_path = new_report_path("pseudo")
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-PSEUDO REPORT =====", 'magenta'), f_out)

    # --- --list-banks: every bundled table and exit ---
    if args.list_banks:
        print_section("[0] BUNDLED BANKS", f_out)
        print_dual(f"  {'name (-p)':<27} {'table':<29} {'format':<6} {'elements':>8}", f_out)
        bank_problems = {}
        for name in bank_names():
            info = bank_info(name)
            n = len(list_bank_elements(name))
            fmt = "PSF" if name == "virtual_vault" else "PSML"
            problems = verify_bank(name)
            if problems:
                bank_problems[name] = problems
            print_dual(f"  {name:<27} {info['table'] or '-':<29} {fmt:<6} {n:>8}", f_out)
        print_dual("", f_out)
        print_dual("Origin (manifest.json of each folder; every file checked against its md5):", f_out)
        for name in bank_names():
            manifest = bank_manifest(name)
            if manifest is None:
                continue
            src = manifest.get("source", {})
            state = "OK" if name not in bank_problems else f"{len(bank_problems[name])} problem(s)"
            print_dual(f"  {name:<27} {src.get('archive', '?')}, taken {src.get('taken', '?')}"
                       f"{', verified ' + src['verified'] if src.get('verified') else ''}: {state}", f_out)
        for name, problems in bank_problems.items():
            for problem in problems:
                print_dual(color_text(f"[WARNING] Bank {name}: {problem}", 'yellow'), f_out)
        print_dual("", f_out)
        print_dual("Short names: words in any order, a missing word takes its default (sr, pbe, standard):", f_out)
        print_dual("  dojo = dojo-sr-pbe-standard, dojo-fr = dojo-fr-pbe-standard, dojo-lda = dojo-sr-pw-standard,", f_out)
        print_dual("  dojo-soc-pbesol = dojo-fr-pbesol-standard, dojo-sr3plus = dojo-sr3plus-pbe-standard.", f_out)
        print_dual("Notes:", f_out)
        print_dual("  - The functional must match XC.Functional/XC.Authors of the calculation.", f_out)
        print_dual("  - stringent: groups 13-16 (Ga-Se, In-Te, Tl-Po) add the s,p shell below the d (Ge/Sn z 14 -> 22);", f_out)
        print_dual("    Hf/Ta add the 4f shell (z 12 -> 26, 13 -> 27). More electrons and orbitals: more expensive.", f_out)
        print_dual("  - fr (full-relativistic): required for 'Spin spin-orbit'. Without SOC, the fr-pbe files run (do not", f_out)
        print_dual("    mix absolute energies of fr and sr files); of the fr-pbesol files only Ge, Sb and Te do -- the", f_out)
        print_dual("    others have no scalar-relativistic part. No full-relativistic LDA table exists.", f_out)
        print_dual("  - sr3plus: lanthanides Ce-Lu with 4f frozen in the core, for Ln3+ compounds only (no 4f physics);", f_out)
        print_dual("    use --fallback-dir dojo for the other elements.", f_out)
        print_dual(f"Reference: {bank_info('dojo')['citation']}", f_out)
        if f_out:
            f_out.close()
        return

    # --- --list-elements: browse a bank and exit, no structure needed ---
    if args.list_elements:
        print_section("[0] BANK ELEMENTS", f_out)
        try:
            bank = bank_info(args.list_elements)
        except ValueError as e:
            _fail(str(e), f_out, report_path)
        print_dual(f"Bank : {bank['name']} -- {bank['description']}", f_out)
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
        if primary_bank:
            print_dual(f"            {bank_info(primary_bank)['description']}", f_out)

        fallback_resolved = None
        fallback_bank = None
        if args.fallback_dir:
            try:
                fallback_resolved = resolve_pseudo_source(args.fallback_dir)
            except ValueError as e:
                _fail(str(e), f_out, report_path)
            print_dual(f"Fallback  : {args.fallback_dir} -> {fallback_resolved}", f_out)
            fallback_bank = detect_pseudo_bank(fallback_resolved)
            if fallback_bank:
                print_dual(f"            {bank_info(fallback_bank)['description']}", f_out)

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

        level_warnings = []
        headers = {}
        for el in species:
            if el in resolved:
                source, path = resolved[el]
                ext = os.path.splitext(path)[1].lstrip(".")
                header = read_pseudo_header(path) or {}
                headers[el] = header
                zval, _note = pseudo_valence(el, os.path.dirname(path))
                zval_txt = f"z={zval:g}" if zval is not None and ext == "psml" else ""
                level = f"{header.get('xc') or '?'}, {header.get('relativity') or '?'}"
                print_dual(f"  {el:<3} FOUND    ({source}, .{ext}, {level}{', ' + zval_txt if zval_txt else ''})  {path}",
                           f_out)
            else:
                print_dual(color_text(f"  {el:<3} MISSING", 'yellow'), f_out)
        xcs = sorted({h.get('xc') for h in headers.values() if h.get('xc')})
        if len(xcs) > 1:
            level_warnings.append(f"the resolved pseudopotentials mix functionals ({', '.join(xcs)}): one calculation "
                                  "needs a single functional.")
        ljs = {bool(h.get('lj')) for h in headers.values() if h}
        if len(ljs) > 1:
            level_warnings.append("the resolved pseudopotentials mix full-relativistic and scalar-relativistic files: "
                                  "spin-orbit needs every file full-relativistic.")
        for msg in level_warnings:
            print_dual(color_text(f"[WARNING] {msg}", 'yellow'), f_out)

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
            # The output folder may already hold the calculation's .fdf: compare its functional
            # and spin-orbit setting with the files just copied (warning only).
            level_warnings += warn_pseudo_level([output_dir], recursive=False)

        print_section("[5] REFERENCES", f_out)
        bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
        for bank in (primary_bank, fallback_bank):
            entry = citations.pseudo_bank_citation(bank)
            if entry and entry not in bib_entries:
                bib_entries.append(entry)
        citations.cite(bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{citations.run_references_path()}' "
            f" (with the libraries the run loads; written when it ends).", 'green'), f_out)

        print_section("[6] SUMMARY & NEXT STEPS", f_out)
        print_dual(f"Status     : {'OK' if ok else 'INCOMPLETE (missing elements)'}", f_out)
        print_dual(f"Resolved   : {len(resolved)}/{len(species)}", f_out)
        if not args.dry_run:
            print_dual(f"Copied to  : {os.path.abspath(os.path.expanduser(args.output))}", f_out)
        print_dual(f"References : {citations.run_references_path()}", f_out)
        if level_warnings:
            print_dual(color_text(f"Warnings   : {len(level_warnings)} (level of theory, see above)", 'yellow'), f_out)
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
