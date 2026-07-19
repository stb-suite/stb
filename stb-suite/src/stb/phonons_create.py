#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import os
import re
import sys
import argparse
from time import sleep
from datetime import datetime
import glob
import numpy as np
from ase import Atoms
from phonopy.interface.siesta import read_siesta, write_siesta, get_physical_units
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos
from stb.core import kspace, mace_relax
from stb.core.deps import require_mace
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders

REPORT_FILE = "phonon_prep_properties.txt"


def print_dual(text, file_handle=None):
    """Prints to stdout with color, writes to file without color. Same
    helper as phonons_pos.py -- duplicated per tool, not factored into
    core/ (presentational, not computational logic)."""
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


def print_symmetry_table(phonon, f_out=None):
    """Reports the symmetry reduction Phonopy already applied (via spglib)
    when generating the displacement dataset -- how many finite-difference
    displacements were actually needed vs. the naive count with no symmetry
    reduction at all (3 Cartesian directions x 2 signs per supercell atom).
    Same table-style presentation as elastic_inputs.py/strain.py's symmetry
    tables, but this is Phonopy's own symmetry analysis, not core/symmetry.py
    (that module solves a different problem -- strain-tensor equivalence,
    not atomic-displacement-pattern reduction).
    """
    sym = phonon.symmetry
    space_group = sym.get_international_table() or "unknown"
    point_group = sym.pointgroup_symbol
    n_ops = len(sym.symmetry_operations['rotations'])
    n_used = len(phonon.dataset['first_atoms'])
    n_naive = phonon.dataset['natom'] * 6

    print_dual("\n" + color_text("--- Symmetry Reduction ---", 'bold'), f_out)
    print_dual(f"Detected symmetry : space group {space_group}, point group {point_group} "
               f"-- {n_ops} operation(s)", f_out)
    print_dual(f"Displacements needed : {color_text(str(n_used), 'green')} "
               f"(vs. {n_naive} without symmetry reduction)", f_out)
    if n_naive > 0:
        print_dual(f"Reduction : {100 * (1 - n_used / n_naive):.1f}% fewer SIESTA runs", f_out)
    print_dual("-" * 60, f_out)


def main():
    parser = argparse.ArgumentParser(
        description="Automate SIESTA phonon displacement folders with Phonopy.",
        epilog="Example usage:\n"
               "  stb_phonon --structure structure.fdf --calc calc.fdf --dim 2 2 2\n"
               "  stb_phonon --pseudo-dir ~/pseudos --distance 0.015",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-dim", type=int, nargs=3, default=[2, 2, 2],
                        help="Supercell dimensions (default: 2 2 2)")

    parser.add_argument("-d", "--distance", type=float, default=0.01,
                        help="Displacement distance in Angstroms (default: 0.01)")

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure file containing the unit cell (default: structure.fdf)")

    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                        help="Input calculation parameters file (default: calc.fdf)")

    parser.add_argument("-p", "--pseudo-dir", type=str, default=".",
                        help=f"Pseudopotentials source: a bundled bank ({', '.join(BANKS)}) or a "
                             "folder path (default: current directory).")

    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum gap (Ang) along an axis to consider it vacuum-padded, "
                             "for the supercell-dimension advisory (default: 10.0)")

    parser.add_argument("--ml-prerelax", action="store_true",
                        help="Optional pre-flight check (opt-in): compute the residual force "
                             "on the input structure with the MACE-MP-0 foundation potential, "
                             "and if it's large, run a quick ML relax for reference. An "
                             "unrelaxed reference structure is a common cause of spurious "
                             "imaginary phonon modes. Diagnostic only -- this run still "
                             "generates displacements from the original --structure file "
                             "unchanged; the ML-relaxed structure is written out separately "
                             "for you to review/use in a follow-up run.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size for --ml-prerelax (default: small).")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu",
                        help="Device for --ml-prerelax (default: cpu).")
    parser.add_argument("--ml-fmax", type=float, default=0.05,
                        help="Force threshold (eV/Ang) for --ml-prerelax: above this, a "
                             "warning is shown and a relax is offered (default: 0.05).")

    parser.add_argument("-v", "--version", action="version", version=f"stb-phononsCreate {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    try:
        args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    # Exibe a introdução se a flag não for acionada
    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("PHONON CALCULATION (VIA PHONOPY):", 'bold'))
    print("-"*60)

    # 1. Validação de arquivos de entrada
    print("\n[INFO] Validating input files ...")
    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found in the current directory.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calculation file '{args.calc}' not found in the current directory.", 'red'))
        sys.exit(1)

    # 2. Leitura da estrutura
    print(f"[INFO] Reading unit cell from '{args.structure}' ...")
    try:
        unitcell = read_siesta(args.structure)
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read {args.structure}. Make sure it's properly formatted.\nDetails: {e}", 'red'))
        sys.exit(1)

    # Diretório de saída + guarda contra re-execução em cima de dados antigos
    # -- movido pra cá (antes só existia no passo 5) pra abrir o relatório
    # persistente cedo e capturar tudo que segue (dimensionalidade, checagem
    # ML, simetria), não só a criação de pastas.
    output_root = "phonon_runs"
    existing_disps = sorted(glob.glob(os.path.join(output_root, "disp-*")))
    if existing_disps:
        print(color_text(
            f"\n[CRITICAL ERROR] '{output_root}' already contains {len(existing_disps)} "
            "displacement folder(s) from a previous run.", 'red'))
        print(color_text(
            "Regenerating on top of them can leave stale disp-* folders (from a "
            "different --dim/--distance/structure) mixed in with the new ones, "
            "silently corrupting FORCE_SETS during post-processing. Remove or move "
            f"aside '{output_root}' and rerun.", 'yellow'))
        sys.exit(1)
    os.makedirs(output_root, exist_ok=True)

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== PHONON PREP REPORT =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Structure file    : {args.structure}", f_out)
        print_dual(f"Calc file         : {args.calc}", f_out)
        print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
        print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
        print_dual(f"Pseudopotentials  : {args.pseudo_dir}", f_out)

        # phonopy's SIESTA interface keeps cell/positions internally in bohr (see
        # the distance-unit note near generate_displacements below) -- convert to
        # Angstrom before handing off to the suite's shared, Angstrom-based
        # vacuum-axis detector.
        bohr_to_angstrom = get_physical_units().Bohr
        lattice_ang = np.array(unitcell.cell) * bohr_to_angstrom
        frac_coords = np.array(unitcell.scaled_positions)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
        print_dual(f"Dimensionality    : {kspace.dimensionality_label(vacuum_axes)}", f_out)

        vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                                  if is_vac and args.dim['abc'.index(axis)] > 1]
        if vacuum_dims_requested:
            print_dual(color_text(
                f"[WARNING] --dim requests more than 1 repetition along vacuum-padded "
                f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
                "Replicating a supercell across vacuum only multiplies the SIESTA cost "
                "without adding real periodicity -- consider --dim 1 on that axis.", 'yellow'), f_out)

        if args.ml_prerelax:
            require_mace()
            print_dual(f"\n{color_text('[0b] ML PRE-FLIGHT CHECK', 'magenta')}", f_out)
            print_dual("-" * 60, f_out)
            print_dual(f"Running MACE-MP-0 ({args.ml_model}) pre-flight check on "
                       f"'{args.structure}' ...", f_out)
            atoms = Atoms(numbers=unitcell.numbers,
                           positions=np.array(unitcell.positions) * bohr_to_angstrom,
                           cell=lattice_ang, pbc=True)
            calc = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
            atoms.calc = calc
            f0 = np.abs(atoms.get_forces()).max()
            print_dual(f"Max residual force on input structure: {f0:.4f} eV/Ang "
                       f"(threshold: {args.ml_fmax} eV/Ang)", f_out)

            if f0 <= args.ml_fmax:
                print_dual(color_text(
                    "Looks relaxed -- a good sign for the finite-difference phonon calculation "
                    "below (which assumes ~zero net force at the reference geometry).", 'green'), f_out)
            else:
                print_dual(color_text(
                    "[WARNING] Residual force is above threshold -- the reference structure may "
                    "not be at a real energy minimum, a common cause of spurious imaginary "
                    "phonon modes. Running a quick MACE relax (positions only) for reference "
                    "...", 'yellow'), f_out)
                converged, steps_used = mace_relax.relax(
                    atoms, calc, cell_mask=None, fmax=args.ml_fmax, max_steps=200)
                f1 = np.abs(atoms.get_forces()).max()
                print_dual(f"After ML relax: max|F| = {f1:.4f} eV/Ang "
                           f"({'converged' if converged else 'hit step cap, not fully converged'}, "
                           f"{steps_used} steps)", f_out)

                relaxed = unitcell.copy()
                relaxed.positions = atoms.get_positions() / bohr_to_angstrom
                relaxed_path = f"{os.path.splitext(args.structure)[0]}_mlrelaxed.fdf"
                write_siesta(relaxed_path, relaxed)
                max_disp = np.linalg.norm(
                    atoms.get_positions() - np.array(unitcell.positions) * bohr_to_angstrom,
                    axis=1).max()
                print_dual(color_text(
                    f"Wrote ML-relaxed structure to '{relaxed_path}' (max atomic displacement "
                    f"from '{args.structure}': {max_disp:.4f} Ang) -- for reference only. This "
                    f"run continues using '{args.structure}' unchanged; rerun with "
                    f"-s {relaxed_path} if you want to use it instead.", 'yellow'), f_out)

        # 3. Extração e validação de Pseudopotenciais
        print_dual(f"\n{color_text('[1] PSEUDOPOTENTIALS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        symbols = unitcell.symbols
        unique_elements = list(set(symbols))
        print_dual(f"Elements in unit cell : {', '.join(unique_elements)}", f_out)
        print_dual(f"Searching in          : {args.pseudo_dir}", f_out)
        pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)

        if missing:
            print_dual(color_text(f"[CRITICAL ERROR] Missing pseudopotentials for the following elements: {', '.join(missing)}", 'red'), f_out)
            print_dual(color_text(f"Action required: Please add the necessary '{missing[0]}.psf' or '{missing[0]}.psml' files into the '{args.pseudo_dir}' directory and rerun the script.", 'yellow'), f_out)
            sys.exit(1)

        print_dual(f"Found all required    : {', '.join([os.path.basename(p) for p in pseudos_to_copy])}", f_out)

        # 4. Inicialização do Phonopy
        print_dual(f"\n{color_text('[2] SYMMETRY REDUCTION', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print(f"[INFO] Generating supercell {args.dim} with {args.distance} Å displacements ...")
        supercell_matrix = [
            [args.dim[0], 0, 0],
            [0, args.dim[1], 0],
            [0, 0, args.dim[2]]
        ]

        # phonopy's SIESTA interface keeps the structure internally in bohr (not
        # Angstrom, unlike every other calculator interface) -- distance has to be
        # converted to that same unit, or the real Cartesian displacement ends up
        # ~1.89x smaller than requested (0.01 Ang asked -> 0.00529 Ang actually
        # applied), verified numerically against this tool's own output. Handled
        # inside build_phonon_displacements() (bohr_to_angstrom computed above,
        # right after reading the structure, is the same physical constant).
        phonon, supercells = build_phonon_displacements(unitcell, supercell_matrix, args.distance)

        print_symmetry_table(phonon, f_out)

        # 5. Criação dos diretórios e cópia
        print_dual(f"\n{color_text('[3] DISPLACEMENT FOLDERS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Building {len(supercells)} displacement folders in '{output_root}' ...", f_out)

        _folders, yaml_path = write_displacement_folders(
            output_root, phonon, supercells, args.structure, args.calc, pseudos_to_copy)

        print_dual(f"Saved Phonopy metadata to '{yaml_path}'", f_out)

        print_dual(f"\n{color_text('[4] SUMMARY & FILES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Displacement folders : {len(supercells)} (disp-001 .. disp-{len(supercells):03d})", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(f"Files                : {yaml_path}, {output_root}/disp-*/", f_out)
        print_dual(color_text(
            f"\n[NOTE] '{os.path.basename(args.calc)}' was copied as-is into every disp-* "
            f"folder. Its k-grid was tuned for the {args.dim[0]}x{args.dim[1]}x{args.dim[2]} "
            "times smaller unit cell -- review it for the generated supercell (roughly "
            "kgrid_unitcell / dim per direction gives the same sampling density at much "
            "lower cost).", 'yellow'), f_out)

    print("\n[INFO] Complete job!")
    print("\n"+"-"*60)
    print(color_text("Phonon folders ready! Let the atoms shake, rattle and roll.\n\n", 'bold'))

if __name__ == "__main__":
    main()
