#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.10.0"

import os
import sys
import shutil
import argparse
from time import sleep
import glob
import numpy as np
from ase import Atoms
from phonopy import Phonopy
from phonopy.interface.siesta import read_siesta, write_siesta, get_physical_units
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source
from stb.core import kspace, mace_relax
from stb.core.deps import require_mace

def get_required_pseudos(symbols: list, pseudo_dir: str):
    """
    Verifica se os pseudopotenciais existem no diretório fornecido.
    Retorna a lista de caminhos encontrados e a lista de elementos ausentes.
    """
    unique_elements = set(symbols)
    found_pseudos = []
    missing_elements = []

    for element in unique_elements:
        psf_path = os.path.join(pseudo_dir, f"{element}.psf")
        psml_path = os.path.join(pseudo_dir, f"{element}.psml")

        if os.path.exists(psf_path):
            found_pseudos.append(psf_path)
        elif os.path.exists(psml_path):
            found_pseudos.append(psml_path)
        else:
            missing_elements.append(element)

    return found_pseudos, missing_elements

def print_symmetry_table(phonon):
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

    print("\n" + color_text("--- Symmetry Reduction ---", 'bold'))
    print(f"Detected symmetry : space group {space_group}, point group {point_group} "
          f"-- {n_ops} operation(s)")
    print(f"Displacements needed : {color_text(str(n_used), 'green')} "
          f"(vs. {n_naive} without symmetry reduction)")
    if n_naive > 0:
        print(f"Reduction : {100 * (1 - n_used / n_naive):.1f}% fewer SIESTA runs")
    print("-" * 60)

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

    # phonopy's SIESTA interface keeps cell/positions internally in bohr (see
    # the distance-unit note near generate_displacements below) -- convert to
    # Angstrom before handing off to the suite's shared, Angstrom-based
    # vacuum-axis detector.
    bohr_to_angstrom = get_physical_units().Bohr
    lattice_ang = np.array(unitcell.cell) * bohr_to_angstrom
    frac_coords = np.array(unitcell.scaled_positions)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
    print(f"[INFO] Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}")

    vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                              if is_vac and args.dim['abc'.index(axis)] > 1]
    if vacuum_dims_requested:
        print(color_text(
            f"[WARNING] --dim requests more than 1 repetition along vacuum-padded "
            f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
            "Replicating a supercell across vacuum only multiplies the SIESTA cost "
            "without adding real periodicity -- consider --dim 1 on that axis.", 'yellow'))

    if args.ml_prerelax:
        require_mace()
        print(f"\n[INFO] Running MACE-MP-0 ({args.ml_model}) pre-flight check on "
              f"'{args.structure}' ...")
        atoms = Atoms(numbers=unitcell.numbers,
                       positions=np.array(unitcell.positions) * bohr_to_angstrom,
                       cell=lattice_ang, pbc=True)
        calc = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
        atoms.calc = calc
        f0 = np.abs(atoms.get_forces()).max()
        print(f"[INFO] Max residual force on input structure: {f0:.4f} eV/Ang "
              f"(threshold: {args.ml_fmax} eV/Ang)")

        if f0 <= args.ml_fmax:
            print(color_text(
                "Looks relaxed -- a good sign for the finite-difference phonon calculation "
                "below (which assumes ~zero net force at the reference geometry).", 'green'))
        else:
            print(color_text(
                "[WARNING] Residual force is above threshold -- the reference structure may "
                "not be at a real energy minimum, a common cause of spurious imaginary "
                "phonon modes. Running a quick MACE relax (positions only) for reference "
                "...", 'yellow'))
            converged, steps_used = mace_relax.relax(
                atoms, calc, cell_mask=None, fmax=args.ml_fmax, max_steps=200)
            f1 = np.abs(atoms.get_forces()).max()
            print(f"[INFO] After ML relax: max|F| = {f1:.4f} eV/Ang "
                  f"({'converged' if converged else 'hit step cap, not fully converged'}, "
                  f"{steps_used} steps)")

            relaxed = unitcell.copy()
            relaxed.positions = atoms.get_positions() / bohr_to_angstrom
            relaxed_path = f"{os.path.splitext(args.structure)[0]}_mlrelaxed.fdf"
            write_siesta(relaxed_path, relaxed)
            max_disp = np.linalg.norm(
                atoms.get_positions() - np.array(unitcell.positions) * bohr_to_angstrom,
                axis=1).max()
            print(color_text(
                f"Wrote ML-relaxed structure to '{relaxed_path}' (max atomic displacement "
                f"from '{args.structure}': {max_disp:.4f} Ang) -- for reference only. This "
                f"run continues using '{args.structure}' unchanged; rerun with "
                f"-s {relaxed_path} if you want to use it instead.", 'yellow'))

    # 3. Extração e validação de Pseudopotenciais
    symbols = unitcell.symbols
    unique_elements = list(set(symbols))
    print(f"[INFO] Elements found in unit cell: {', '.join(unique_elements)}")
    
    print(f"[INFO] Searching for pseudopotentials in '{args.pseudo_dir}' ...")
    pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)

    if missing:
        print(color_text(f"\n[CRITICAL ERROR] Missing pseudopotentials for the following elements: {', '.join(missing)}", 'red'))
        print(color_text(f"Action required: Please add the necessary '{missing[0]}.psf' or '{missing[0]}.psml' files into the '{args.pseudo_dir}' directory and rerun the script.", 'yellow'))
        sys.exit(1)
        
    print(f"[INFO] Found all required pseudopotentials: {', '.join([os.path.basename(p) for p in pseudos_to_copy])}")

    # 4. Inicialização do Phonopy
    print(f"[INFO] Generating supercell {args.dim} with {args.distance} Å displacements ...")
    supercell_matrix = [
        [args.dim[0], 0, 0],
        [0, args.dim[1], 0],
        [0, 0, args.dim[2]]
    ]

    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix, calculator="siesta")
    # phonopy's SIESTA interface keeps the structure internally in bohr (not
    # Angstrom, unlike every other calculator interface) -- distance has to be
    # converted to that same unit, or the real Cartesian displacement ends up
    # ~1.89x smaller than requested (0.01 Ang asked -> 0.00529 Ang actually
    # applied), verified numerically against this tool's own output.
    # (bohr_to_angstrom computed above, right after reading the structure.)
    phonon.generate_displacements(distance=args.distance / bohr_to_angstrom)
    supercells = phonon.supercells_with_displacements

    print_symmetry_table(phonon)

    # 5. Criação dos diretórios e cópia
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

    print(f"[INFO] Building {len(supercells)} displacement folders in '{output_root}' ...")

    for i, scell in enumerate(supercells):
        if scell is None:
            continue
            
        folder_name = os.path.join(output_root, f"disp-{i+1:03d}")
        os.makedirs(folder_name, exist_ok=True)
        
        # A. Escreve a supercélula deslocada
        disp_struct_path = os.path.join(folder_name, os.path.basename(args.structure))
        write_siesta(disp_struct_path, scell)

        # B. Copia o calc.fdf
        shutil.copy(args.calc, os.path.join(folder_name, os.path.basename(args.calc)))
        
        # C. Copia apenas os pseudopotenciais exigidos
        for pseudo_path in pseudos_to_copy:
            pseudo_filename = os.path.basename(pseudo_path)
            shutil.copy(pseudo_path, os.path.join(folder_name, pseudo_filename))

    # 6. Salvar metadados do Phonopy
    yaml_path = os.path.join(output_root, "phonopy_disp.yaml")
    phonon.save(yaml_path)
    print(f"[INFO] Saved Phonopy metadata to '{yaml_path}'")
    
    print("\n[INFO] Complete job!")
    print(color_text(
        f"\n[NOTE] '{os.path.basename(args.calc)}' was copied as-is into every disp-* "
        f"folder. Its k-grid was tuned for the {args.dim[0]}x{args.dim[1]}x{args.dim[2]} "
        "times smaller unit cell -- review it for the generated supercell (roughly "
        "kgrid_unitcell / dim per direction gives the same sampling density at much "
        "lower cost).", 'yellow'))
    print("\n"+"-"*60)
    print(color_text("Phonon folders ready! Let the atoms shake, rattle and roll.\n\n", 'bold'))

if __name__ == "__main__":
    main()
