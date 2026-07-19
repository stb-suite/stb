#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import argparse
import glob
from datetime import datetime
import numpy as np
from phonopy.interface.siesta import read_siesta, get_physical_units
from stb.core.cli import color_text, show_intro
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos
from stb.core import kspace
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders

REPORT_FILE = "raman_stage1.txt"


def print_dual(text, file_handle=None):
    """Prints to stdout with color, writes to file without color. Same
    duplicated-per-tool helper as phonons_create.py."""
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 3: generates the phonon-displacement folders for a "
        "Raman spectrum workflow.", 'bold')}
Workflow: (1) this tool generates 'raman_study/phonon_disp/disp-*/' -- run
SIESTA in every folder yourself; (2) once they're done, run stb-ramanModes
to build FORCE_SETS from them, identify the Gamma-point vibrational modes,
and generate the +/-delta Optical-calculation displacement folders needed
for the Raman tensor; (3) run SIESTA in those, then stb-ramanAnalysis to
get the Raman-active frequencies/activities and spectrum. This whole
workflow computes its own phonons -- it does NOT need item 4's
stb-phononsCreate/stb-phononsPos run first (though the underlying Phonopy
machinery is shared, via core/phonon_workflow.py).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p ~/pseudos -dim 2 2 2 -d 0.02\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure file containing the unit cell (default: structure.fdf)")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                        help="Input calculation parameters file (default: calc.fdf)")
    parser.add_argument("-p", "--pseudo-dir", type=str, default=".",
                        help=f"Pseudopotentials source: a bundled bank ({', '.join(BANKS)}) or a "
                             "folder path (default: current directory).")
    parser.add_argument("-dim", type=int, nargs=3, default=[2, 2, 2],
                        help="Supercell dimensions for the phonon force-constant calculation "
                             "(default: 2 2 2).")
    parser.add_argument("-d", "--distance", type=float, default=0.02,
                        help="Phonon displacement distance in Angstroms (default: 0.02).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum gap (Ang) along an axis to consider it vacuum-padded, "
                             "for the supercell-dimension advisory (default: 10.0)")
    parser.add_argument("-O", "--output-dir", type=str, default="raman_study",
                        help="Root directory for the whole Raman workflow (default: raman_study) -- "
                             "phonon displacements are written under <output-dir>/phonon_disp/.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-raman {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    try:
        args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("RAMAN SPECTRUM WORKFLOW -- STAGE 1: PHONON DISPLACEMENTS", 'bold'))
    print("-" * 60)

    print("\n[INFO] Validating input files ...")
    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calculation file '{args.calc}' not found.", 'red'))
        sys.exit(1)

    print(f"[INFO] Reading unit cell from '{args.structure}' ...")
    try:
        unitcell = read_siesta(args.structure)
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read {args.structure}. Make sure it's properly "
                          f"formatted.\nDetails: {e}", 'red'))
        sys.exit(1)

    output_root = args.output_dir
    phonon_disp_dir = os.path.join(output_root, "phonon_disp")
    existing_disps = sorted(glob.glob(os.path.join(phonon_disp_dir, "disp-*")))
    if existing_disps:
        print(color_text(
            f"\n[CRITICAL ERROR] '{phonon_disp_dir}' already contains {len(existing_disps)} "
            "displacement folder(s) from a previous run.", 'red'))
        print(color_text(
            "Regenerating on top of them can leave stale disp-* folders mixed in with the new "
            f"ones, silently corrupting FORCE_SETS during Stage 2. Remove or move aside "
            f"'{output_root}' and rerun.", 'yellow'))
        sys.exit(1)
    os.makedirs(phonon_disp_dir, exist_ok=True)

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== RAMAN STAGE 1 REPORT (PHONON DISPLACEMENTS) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Structure file    : {args.structure}", f_out)
        print_dual(f"Calc file         : {args.calc}", f_out)
        print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
        print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
        print_dual(f"Pseudopotentials  : {args.pseudo_dir}", f_out)
        print_dual(f"Output root       : {output_root}", f_out)

        bohr_to_angstrom = get_physical_units().Bohr
        lattice_ang = np.array(unitcell.cell) * bohr_to_angstrom
        frac_coords = np.array(unitcell.scaled_positions)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
        print_dual(f"Dimensionality    : {kspace.dimensionality_label(vacuum_axes)}", f_out)

        vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                                  if is_vac and args.dim['abc'.index(axis)] > 1]
        if vacuum_dims_requested:
            print_dual(color_text(
                f"[WARNING] -dim requests more than 1 repetition along vacuum-padded "
                f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
                "Replicating a supercell across vacuum only multiplies the SIESTA cost without "
                "adding real periodicity -- consider -dim 1 on that axis.", 'yellow'), f_out)

        print_dual(f"\n{color_text('[1] PSEUDOPOTENTIALS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        symbols = unitcell.symbols
        unique_elements = list(set(symbols))
        print_dual(f"Elements in unit cell : {', '.join(unique_elements)}", f_out)
        print_dual(f"Searching in          : {args.pseudo_dir}", f_out)
        pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)

        if missing:
            print_dual(color_text(
                f"[CRITICAL ERROR] Missing pseudopotentials for the following elements: "
                f"{', '.join(missing)}", 'red'), f_out)
            print_dual(color_text(
                f"Action required: add the necessary '{missing[0]}.psf' or '{missing[0]}.psml' "
                f"files into '{args.pseudo_dir}' and rerun.", 'yellow'), f_out)
            sys.exit(1)

        print_dual(f"Found all required    : "
                    f"{', '.join(os.path.basename(p) for p in pseudos_to_copy)}", f_out)

        print_dual(f"\n{color_text('[2] PHONON DISPLACEMENTS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print(f"[INFO] Generating supercell {args.dim} with {args.distance} Ang displacements ...")
        supercell_matrix = [
            [args.dim[0], 0, 0],
            [0, args.dim[1], 0],
            [0, 0, args.dim[2]]
        ]
        # phonopy's SIESTA interface keeps the structure internally in bohr --
        # distance conversion is handled inside build_phonon_displacements()
        # (see core/phonon_workflow.py for the ~1.89x pitfall this avoids).
        phonon, supercells = build_phonon_displacements(unitcell, supercell_matrix, args.distance)

        n_used = len(phonon.dataset['first_atoms'])
        n_naive = phonon.dataset['natom'] * 6
        print_dual(f"Displacements needed  : {color_text(str(n_used), 'green')} "
                    f"(vs. {n_naive} without symmetry reduction)", f_out)

        print_dual(f"Building {len(supercells)} displacement folders in '{phonon_disp_dir}' ...", f_out)
        folders, yaml_path = write_displacement_folders(
            phonon_disp_dir, phonon, supercells, args.structure, args.calc, pseudos_to_copy)
        print_dual(f"Saved Phonopy metadata to '{yaml_path}'", f_out)

        print_dual(f"\n{color_text('[3] SUMMARY & NEXT STEPS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Displacement folders : {len(folders)} (disp-001 .. disp-{len(folders):03d})", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(f"Files                : {yaml_path}, {phonon_disp_dir}/disp-*/", f_out)
        print_dual(color_text(
            f"\n[NOTE] '{os.path.basename(args.calc)}' was copied as-is into every disp-* folder. "
            f"Its k-grid was tuned for the {args.dim[0]}x{args.dim[1]}x{args.dim[2]} times smaller "
            "unit cell -- review it for the generated supercell (roughly kgrid_unitcell / dim per "
            "direction gives the same sampling density at much lower cost).", 'yellow'), f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual(f"  1. Run SIESTA in every '{phonon_disp_dir}/disp-*/' folder.", f_out)
        print_dual(f"  2. Once they're done, run: stb-ramanModes --directory {output_root}", f_out)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Phonon displacement folders ready for Stage 2 (stb-ramanModes).\n", 'bold'))


if __name__ == "__main__":
    main()
