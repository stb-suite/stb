#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import argparse
import sys
import numpy as np
from stb.core import structure_io, kspace, citations
from stb.core.cli import color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_kgrid_report.txt"
BIB_FILE = "references.bib"


def build_bib_entries():
    """stb-kgrid always computes a genuine Monkhorst-Pack grid -- unlike
    stb-inputfile (which skips it for AIMD), there's no mode here where the
    k-grid isn't actually the point of the run."""
    return [citations.SIESTA, citations.SIESTA_RECENT, citations.MONKHORST_PACK]


def main():
    parser = argparse.ArgumentParser(
        description="Compute a Monkhorst-Pack k-point grid from a SIESTA "
                    "structure file (.fdf) and a target k-point density.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf -d 0.2\n"
               "  %(prog)s -f slab.fdf -d 0.2 --vacuum-gap 15\n"
    )
    parser.add_argument(
        "-f", "--file", type=str, required=True, dest="structure_file",
        help="Path to the SIESTA structure file (.fdf)."
    )
    parser.add_argument(
        "-d", "--density", type=float, required=True,
        help="Target k-point density (in 1/Ang). Example: 0.2"
    )
    parser.add_argument(
        "--vacuum-gap", type=float, default=10.0,
        help="Minimum empty span (in Angstrom) between atoms along an axis, wrapped "
             "periodically, to treat that axis as vacuum-padded (non-periodic) and force "
             "a single k-point there regardless of density. Default: 10.0"
    )
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kgrid {VERSION}")
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

    print_dual(color_text("===== STB-KGRID REPORT =====", 'magenta'), f_out)

    ok = True
    try:
        print_section("[0] RUN METADATA", f_out)
        print_dual(f"Structure file : {args.structure_file}", f_out)
        print_dual(f"Target density : {args.density} 1/Ang", f_out)
        print_dual(f"Vacuum-gap threshold : {args.vacuum_gap} Ang", f_out)

        structure = structure_io.read_fdf(args.structure_file)
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'

        print_section("[1] K-POINT DENSITY GUIDE", f_out)
        kspace.print_density_recommendation(f_out)

        print_section("[2] K-GRID CALCULATION", f_out)
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
        print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes)}", f_out)

        divisions = kspace.compute_monkhorts(
            structure.lattice[0], structure.lattice[1], structure.lattice[2],
            args.density, vacuum_axes
        )
        print_dual(
            f"Suggested Monkhorst-Pack grid : {divisions[0]} {divisions[1]} {divisions[2]}",
            f_out
        )

        vacuum_count = sum(vacuum_axes)
        if vacuum_count > 0:
            axis_word = "axis" if vacuum_count == 1 else "axes"
            print_dual(
                f"The '1' division(s) above correspond to the {vacuum_count} vacuum-padded "
                f"{axis_word} detected -- no real periodicity there, so a single k-point "
                "is enough.",
                f_out
            )

        print_section("[3] REFERENCES", f_out)
        bib_entries = build_bib_entries()
        citations.write_bib_file(BIB_FILE, bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)

        print_section("[4] SUMMARY & FILES", f_out)
        print_dual("Status         : OK", f_out)
        print_dual(f"Suggested grid : {divisions[0]} {divisions[1]} {divisions[2]}", f_out)
        print_dual(f"References     : {BIB_FILE}", f_out)

    except FileNotFoundError:
        print_dual(
            color_text(f"[ERROR] Structure file '{args.structure_file}' not found.", 'red'),
            f_out
        )
        ok = False
    except Exception as e:
        print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
        ok = False

    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
