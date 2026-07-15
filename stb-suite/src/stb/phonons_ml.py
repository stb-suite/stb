#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.1.0"

import os
import re
import sys
import argparse
from datetime import datetime
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io, kspace, mace_relax, mace_phonons
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.deps import require_mace

REPORT_FILE = "phonon_prep_properties.txt"


def print_dual(text, file_handle=None):
    """Prints to stdout with color, writes to file without color. Same
    helper as phonons_create.py/phonons_pos.py/phonons_qha.py -- duplicated
    per tool, not factored into core/ (presentational, not computational
    logic)."""
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


def print_symmetry_table(phonon, f_out=None):
    """Reports the symmetry reduction Phonopy applied (via spglib) when
    generating the displacement dataset. Same content/style as
    phonons_create.py's print_symmetry_table -- duplicated rather than
    imported across sibling CLI modules, same precedent as print_dual (a
    presentational helper, not core computational logic, kept bespoke per
    tool in this workflow).
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
        print_dual(f"Reduction : {100 * (1 - n_used / n_naive):.1f}% fewer MACE evaluations", f_out)
    print_dual("-" * 60, f_out)


def main():
    parser = argparse.ArgumentParser(
        description="Fast phonon force-constant calculation using the MACE-MP-0 foundation "
                     "potential -- no SIESTA calculation needed. A heuristic preview "
                     "(dispersion, DOS, thermal properties, all via the existing "
                     "stb-phononsPos) to sanity-check a structure or explore trends before "
                     "committing to a real SIESTA phonon calculation -- same philosophy as "
                     "stb-mlrelax: NOT a DFT replacement.",
        epilog="Example usage:\n"
               "  stb-phononsML -s structure.fdf -dim 2 2 2\n"
               "  stb-phononsML -s structure.fdf --model medium --device cuda\n"
               "\n"
               "Then analyze exactly like a SIESTA-based run:\n"
               "  stb-phononsPos -dir phonon_ml_runs --bands --dos\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure file containing the unit cell (default: structure.fdf)")
    parser.add_argument("-dim", type=int, nargs=3, default=[2, 2, 2],
                        help="Supercell dimensions (default: 2 2 2)")
    parser.add_argument("-d", "--distance", type=float, default=0.01,
                        help="Displacement distance in Angstroms (default: 0.01). Applied "
                             "directly -- unlike stb-phononsCreate's SIESTA path, this tool "
                             "builds the structure natively in Angstrom (ASE/MACE's own "
                             "convention), so there's no bohr-vs-Angstrom internal-unit "
                             "conversion needed here.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum gap (Ang) along an axis to consider it vacuum-padded, "
                             "for the supercell-dimension advisory (default: 10.0)")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size: speed/accuracy tradeoff (default: small).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the MACE pre-relax (positions only) before generating "
                             "displacements. On by default: a phonon calculation assumes "
                             "~zero net force at the reference geometry, and skipping this "
                             "on an unrelaxed structure is a common cause of spurious "
                             "imaginary modes.")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the pre-relax, eV/Ang (default: 0.05).")
    parser.add_argument("-o", "--output-dir", type=str, default="phonon_ml_runs",
                        help="Output directory for phonopy_disp.yaml (default: phonon_ml_runs) "
                             "-- point stb-phononsPos -dir here afterwards.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-phononsML {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("PHONON CALCULATION (VIA MACE-MP-0 -- NO SIESTA NEEDED):", 'bold'))
    print(color_text("Heuristic preview, not a DFT replacement -- same caveat as stb-mlrelax.", 'yellow'))
    print("-" * 60)

    # 1. Validação e leitura da estrutura
    print("\n[INFO] Validating input files ...")
    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found in the current directory.", 'red'))
        sys.exit(1)

    print(f"[INFO] Reading unit cell from '{args.structure}' ...")
    try:
        structure = structure_io.read_fdf(args.structure)
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read {args.structure}. Make sure it's properly formatted.\nDetails: {e}", 'red'))
        sys.exit(1)

    pmg = structure_io.to_pymatgen(structure)
    atoms = AseAtomsAdaptor.get_atoms(pmg)

    # Diretório de saída + guarda contra re-execução em cima de dados antigos
    # -- feito cedo pra abrir o relatório persistente antes de qualquer
    # trabalho pesado (mesmo padrão de phonons_create.py).
    output_root = args.output_dir
    existing_yaml = os.path.join(output_root, "phonopy_disp.yaml")
    if os.path.exists(existing_yaml):
        print(color_text(
            f"\n[CRITICAL ERROR] '{existing_yaml}' already exists from a previous run.", 'red'))
        print(color_text(
            f"Remove or move aside '{output_root}' and rerun, to avoid silently mixing "
            "data from a different --dim/--distance/structure.", 'yellow'))
        sys.exit(1)
    os.makedirs(output_root, exist_ok=True)

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== PHONON PREP REPORT (MACE-MP-0) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Structure file    : {args.structure}", f_out)
        print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
        print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
        print_dual(f"MACE-MP-0 model   : {args.model} ({args.device})", f_out)
        print_dual(f"Elements in unit cell : {', '.join(sorted(set(atoms.get_chemical_symbols())))}", f_out)

        # Dimensionalidade / vácuo -- já nativamente em Å (sem conversão bohr,
        # ao contrário do caminho SIESTA deste workflow).
        lattice_ang = np.array(atoms.get_cell())
        frac_coords = atoms.get_scaled_positions()
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
        print_dual(f"Dimensionality    : {kspace.dimensionality_label(vacuum_axes)}", f_out)

        vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                                  if is_vac and args.dim['abc'.index(axis)] > 1]
        if vacuum_dims_requested:
            print_dual(color_text(
                f"[WARNING] --dim requests more than 1 repetition along vacuum-padded "
                f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
                "Replicating a supercell across vacuum only multiplies the MACE cost "
                "without adding real periodicity -- consider --dim 1 on that axis.", 'yellow'), f_out)

        # 4. Calculadora MACE + relaxação opcional (só posições -- a célula já
        # define a supercélula abaixo; relaxá-la mudaria a periodicidade que os
        # deslocamentos vão assumir)
        print_dual(f"\n{color_text('[0b] ML PRE-RELAX', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print(f"[INFO] Loading MACE-MP-0 ({args.model}) on {args.device} ...")
        calc = mace_relax.get_calculator(model=args.model, device=args.device)
        atoms.calc = calc

        f0 = np.abs(atoms.get_forces()).max()
        print_dual(f"Residual force on input structure: {f0:.4f} eV/Ang", f_out)
        if not args.relax:
            print_dual(color_text(
                "[WARNING] --no-relax: generating displacements from the structure as given. "
                "If it isn't already at a MACE-relaxed minimum, expect spurious imaginary "
                "modes downstream.", 'yellow'), f_out)
        else:
            print(f"[INFO] Relaxing positions (fmax={args.fmax}) before generating displacements ...")

        # 5. Deslocamentos (core.mace_phonons -- reaproveitado também pelo
        # stb-phononsQHA)
        print(f"\n[INFO] Generating supercell {args.dim} with {args.distance} Å displacements ...")
        phonon, _, relax_info = mace_phonons.generate_ml_displacements(
            atoms, args.dim, args.distance, calc, relax=args.relax, fmax=args.fmax)

        if relax_info is not None:
            _, f1, converged, steps_used = relax_info
            print_dual(f"After relax: max|F| = {f1:.4f} eV/Ang "
                       f"({'converged' if converged else 'hit step cap, not fully converged'}, "
                       f"{steps_used} steps)", f_out)

        print_dual(f"\n{color_text('[1] SYMMETRY REDUCTION', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_symmetry_table(phonon, f_out)

        # 6. Forças via MACE para cada supercélula deslocada + constantes de força
        print_dual(f"\n{color_text('[2] FORCE CONSTANTS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        n_disp = len(phonon.supercells_with_displacements)
        print_dual(f"Computing MACE forces for {n_disp} displaced supercells ...", f_out)
        report_every = max(1, n_disp // 10)

        def _progress(i, n):
            if i % report_every == 0 or i == n:
                print(f"  ... {i}/{n}")

        mace_phonons.compute_force_constants(phonon, calc, progress_callback=_progress)

        # 7. Salvar phonopy_disp.yaml já com as constantes de força embutidas --
        # stb-phononsPos detecta isso e pula a extração de força via SIESTA.
        yaml_path = os.path.join(output_root, "phonopy_disp.yaml")
        phonon.save(yaml_path, settings={'force_constants': True})
        print_dual(f"Saved force constants to '{yaml_path}'", f_out)

        print_dual(f"\n{color_text('[3] SUMMARY & FILES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Displacements computed : {n_disp}", f_out)
        print_dual(f"Report                 : {report_path}", f_out)
        print_dual(f"Files                  : {yaml_path}", f_out)

    print("\n" + "-" * 60)
    print(color_text("Phonon force constants ready (MACE-MP-0, no SIESTA)!", 'bold'))
    print(f"Next: {color_text(f'stb-phononsPos -dir {output_root}', 'cyan')} "
          "(add --bands/--dos/--pdos/--thermal-displacements/--freeze-unstable-mode "
          "exactly like a SIESTA-based run).\n")


if __name__ == "__main__":
    main()
