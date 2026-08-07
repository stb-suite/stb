#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.13.0"  # --symprec (default 0.01, pymatgen's own default) now threaded through
                    # to Phonopy itself -- a real bug: Phonopy's own raw default (1e-5) is far
                    # tighter than any DFT relaxation's real numerical noise floor and was
                    # silently misdetecting the true space group (verified live on a real
                    # relaxed AlP structure: 1e-5 -> wrong R3m, 0.01 -> correct F-43m). Also:
                    # library prints/warnings (MACE/torch/phonopy) captured and moved to a new
                    # final [LIBRARY WARNINGS] section instead of interleaving with the report;
                    # [0b] ML PRE-FLIGHT CHECK gets an explicit Verdict line; symmetry section
                    # rewritten as a print_table with symprec shown and the result highlighted

import os
import sys
import argparse
from time import sleep
from datetime import datetime
import glob
import numpy as np
from ase import Atoms
from phonopy.interface.siesta import read_siesta, write_siesta, get_physical_units
from stb.core.cli import (COLORS, color_text, show_intro, print_dual, print_section,
                          print_table, capture_library_noise)
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos
from stb.core import kspace, mace_relax
from stb.core.deps import require_mace
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders
from stb.core.structure_io import read_md_state, prepend_include

REPORT_FILE = "phonon_prep_properties.txt"
EXTRA_FDF_FILE = "config_extra.fdf"


def print_symmetry_table(phonon, symprec, f_out=None):
    """Reports the symmetry reduction Phonopy already applied (via spglib)
    when generating the displacement dataset -- how many finite-difference
    displacements were actually needed vs. the naive count with no symmetry
    reduction at all (3 Cartesian directions x 2 signs per supercell atom).
    This is Phonopy's own symmetry analysis (its own `symprec`, shown
    explicitly here), not core/symmetry.py's (that module solves a
    different problem -- strain-tensor equivalence, not atomic
    -displacement-pattern reduction) -- the two can legitimately disagree
    if their tolerances differ.
    """
    sym = phonon.symmetry
    space_group = sym.get_international_table() or "unknown"
    point_group = sym.pointgroup_symbol
    n_ops = len(sym.symmetry_operations['rotations'])
    n_used = len(phonon.dataset['first_atoms'])
    n_naive = phonon.dataset['natom'] * 6
    reduction_pct = 100 * (1 - n_used / n_naive) if n_naive > 0 else 0.0

    print_dual(color_text(f"\n>>> Detected space group: {space_group}  "
                          f"(symprec={symprec:g} Ang) <<<", 'bold'), f_out)
    print_table(["Quantity", "Value"], [
        (["Symmetry precision (symprec)", f"{symprec:g} Ang"], None),
        (["Space group", space_group], None),
        (["Point group", point_group], None),
        (["Symmetry operations", str(n_ops)], None),
        (["Displacements needed", f"{n_used} (of {n_naive} without symmetry reduction)"], 'green'),
        (["Reduction", f"{reduction_pct:.1f}% fewer SIESTA runs"], 'green'),
    ], f_out)


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

    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry-detection tolerance (Ang) Phonopy uses internally to "
                             "reduce how many displacements are actually needed (default: 0.01, "
                             "pymatgen's own default -- matches the rest of the suite, and "
                             "deliberately NOT Phonopy's own raw default of 1e-5, which is far "
                             "too tight for a real DFT-relaxed structure and can misdetect the "
                             "true space group: verified live on a real relaxed AlP zincblende "
                             "structure, 1e-5 reports the wrong R3m while 0.01 correctly reports "
                             "F-43m). Loosen further (e.g. 0.02-0.05) for a structure relaxed "
                             "with a looser force tolerance.")

    parser.add_argument("--kgrid", type=int, nargs=3, default=None,
                        help="Explicit Monkhorst-Pack k-grid 'X Y Z' for the SUPERCELL's own "
                             "single-point SCF (e.g. --kgrid 2 2 2). Overrides the automatic "
                             "--kgrid-density suggestion. Written into config_extra.fdf (see "
                             "[1]) as the first directive, so it takes precedence over --calc's "
                             "own k-grid, which was tuned for the smaller unit cell.")
    parser.add_argument("--kgrid-density", type=float, default=0.2,
                        help="K-point density (1/Ang) used to auto-suggest the supercell's own "
                             "k-grid when --kgrid isn't given (default: 0.2, same convention as "
                             "stb-kgrid/stb-strain).")

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

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-phononsCreate {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()
    library_warnings = []

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

    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== PHONON PREP REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc file         : {args.calc}", f_out)
    print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
    print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
    print_dual(f"Symmetry precision: {args.symprec:g} Ang (symprec, see [3])", f_out)
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
        with capture_library_noise(library_warnings, "MACE import"):
            require_mace()
        print_section('[0b] ML PRE-FLIGHT CHECK', f_out)
        print_dual(f"Model             : MACE-MP-0 ({args.ml_model}, device={args.ml_device})", f_out)
        print_dual(f"Force threshold   : {args.ml_fmax} eV/Ang", f_out)
        atoms = Atoms(numbers=unitcell.numbers,
                       positions=np.array(unitcell.positions) * bohr_to_angstrom,
                       cell=lattice_ang, pbc=True)
        with capture_library_noise(library_warnings, "MACE calculator setup"):
            calc = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
        atoms.calc = calc
        f0 = np.abs(atoms.get_forces()).max()
        print_dual(f"Max residual force on input structure: {f0:.4f} eV/Ang", f_out)

        if f0 <= args.ml_fmax:
            print_dual(color_text(
                f"\nVerdict: OK -- structure looks relaxed ({f0:.4f} <= {args.ml_fmax} eV/Ang), "
                "a good sign for the finite-difference phonon calculation below (which assumes "
                "~zero net force at the reference geometry).", 'green'), f_out)
        else:
            print_dual(color_text(
                f"\n[WARNING] Residual force ({f0:.4f} eV/Ang) is above the {args.ml_fmax} "
                "eV/Ang threshold -- the reference structure may not be at a real energy "
                "minimum, a common cause of spurious imaginary phonon modes. Running a quick "
                "MACE relax (positions only) for reference ...", 'yellow'), f_out)
            with capture_library_noise(library_warnings, "MACE relax"):
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
            print_dual(f"Wrote ML-relaxed structure to '{relaxed_path}' (max atomic "
                       f"displacement from '{args.structure}': {max_disp:.4f} Ang) -- for "
                       "reference only.", f_out)
            print_dual(color_text(
                f"\nVerdict: RELAXATION RECOMMENDED -- this run continues using "
                f"'{args.structure}' unchanged; rerun with -s {relaxed_path} first if you "
                "want the ML pre-relaxed geometry instead.", 'yellow'), f_out)

    # --- [1] SINGLE-POINT SCF ENFORCEMENT ---
    print_section('[1] SINGLE-POINT SCF ENFORCEMENT', f_out)
    with open(args.calc) as f:
        original_calc_text = f.read()
    structure_basename = os.path.basename(args.structure)
    calc_basename = os.path.basename(args.calc)
    if "%include" not in original_calc_text or structure_basename not in original_calc_text:
        print_dual(color_text(
            f"[NOTE] '{calc_basename}' may not %include '{structure_basename}' -- add it if "
            "needed so SIESTA picks up each displaced supercell's geometry.", 'yellow'), f_out)

    before = read_md_state(original_calc_text)
    steps_label = f"{before['steps_key']}={before['steps']}" if before['steps_key'] else "(absent)"
    print_dual(f"Calc template (before forcing): MD.TypeOfRun={before['typeofrun'] or '(absent)'}  "
               f"Steps: {steps_label}  MD.VariableCell={before['variablecell'] or '(absent)'}", f_out)
    print_dual(f"Forced to single-point SCF via '%include {EXTRA_FDF_FILE}'.", f_out)

    # Supercell k-grid: --calc's own k-grid was tuned for the (smaller)
    # unit cell, so it over-samples the supercell if left as-is -- compute
    # (or take explicitly via --kgrid) a grid sized for the supercell
    # itself, at --kgrid-density (default 0.2, same convention as
    # stb-kgrid/stb-strain), and put it FIRST in config_extra.fdf (before
    # the MD block below) so it's the first directive SIESTA reads.
    supercell_lattice = np.diag(args.dim) @ lattice_ang
    if args.kgrid is not None:
        kgrid = tuple(args.kgrid)
        kgrid_source = "explicit --kgrid"
    else:
        try:
            kgrid = tuple(kspace.compute_monkhorts(
                supercell_lattice[0], supercell_lattice[1], supercell_lattice[2],
                args.kgrid_density, vacuum_axes))
        except ValueError as e:
            print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        kgrid_source = f"auto-suggested, density={args.kgrid_density:g} 1/Ang"
    print_dual(f"\nSupercell k-grid  : {kgrid[0]} {kgrid[1]} {kgrid[2]} ({kgrid_source})", f_out)

    extra_fdf_text = (
        "# Auto-generated by stb-phononsCreate.\n"
        f"# Supercell k-grid ({kgrid_source}) -- takes precedence over --calc's own k-grid,\n"
        "# which was tuned for the smaller unit cell.\n"
        f"kgrid.MonkhorstPack   [{kgrid[0]}  {kgrid[1]}  {kgrid[2]}]\n"
        "\n"
        "# Forces a pure single-point SCF (no ionic or cell relaxation) at each displaced\n"
        "# supercell's exact geometry, regardless of --calc's own settings.\n"
        "MD.TypeOfRun       CG\n"
        "MD.Steps           0\n"
        "MD.NumCGsteps      0\n"
        "MD.VariableCell    false\n"
    )
    print_dual(f"\n{EXTRA_FDF_FILE} (written into every generated folder):", f_out)
    for line in extra_fdf_text.rstrip("\n").split("\n"):
        print_dual(f"  {line}", f_out)
    forced_calc_text = prepend_include(original_calc_text, EXTRA_FDF_FILE)

    # 3. Extração e validação de Pseudopotenciais
    print_section('[2] PSEUDOPOTENTIALS', f_out)
    symbols = unitcell.symbols
    unique_elements = list(set(symbols))
    print_dual(f"Elements in unit cell : {', '.join(unique_elements)}", f_out)
    print_dual(f"Searching in          : {args.pseudo_dir}", f_out)
    pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)

    if missing:
        print_dual(color_text(f"[CRITICAL ERROR] Missing pseudopotentials for the following elements: {', '.join(missing)}", 'red'), f_out)
        print_dual(color_text(f"Action required: Please add the necessary '{missing[0]}.psf' or '{missing[0]}.psml' files into the '{args.pseudo_dir}' directory and rerun the script.", 'yellow'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(f"Found all required    : {', '.join([os.path.basename(p) for p in pseudos_to_copy])}", f_out)

    # 4. Inicialização do Phonopy
    print_section('[3] SYMMETRY REDUCTION', f_out)
    print_dual(f"Generating supercell {args.dim} with {args.distance} Ang displacements ...", f_out)
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
    with capture_library_noise(library_warnings, "Phonopy"):
        phonon, supercells = build_phonon_displacements(
            unitcell, supercell_matrix, args.distance, symprec=args.symprec)

    print_symmetry_table(phonon, args.symprec, f_out)

    # 5. Criação dos diretórios e cópia
    print_section('[4] DISPLACEMENT FOLDERS', f_out)
    print_dual(f"Building {len(supercells)} displacement folders in '{output_root}' ...", f_out)

    _folders, yaml_path = write_displacement_folders(
        output_root, phonon, supercells, args.structure, args.calc, pseudos_to_copy)

    # write_displacement_folders copies --calc verbatim (it's shared with
    # her_refs.py/ir.py/oer_refs.py/raman.py, which don't force single-point)
    # -- overwrite each folder's copy with the forced version, and add the
    # config_extra.fdf sidecar, right after (see [1] above).
    for folder in _folders:
        with open(os.path.join(folder, calc_basename), "w") as fh:
            fh.write(forced_calc_text)
        with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
            fh.write(extra_fdf_text)

    print_dual(f"Saved Phonopy metadata to '{yaml_path}'", f_out)

    print_section('[5] SUMMARY & FILES', f_out)
    print_dual(f"Displacement folders : {len(supercells)} (disp-001 .. disp-{len(supercells):03d})", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(f"Files                : {yaml_path}, {output_root}/disp-*/", f_out)
    print_dual(color_text(
        f"\n[NOTE] '{calc_basename}' was forced to single-point SCF with its own supercell "
        f"k-grid ({kgrid[0]} {kgrid[1]} {kgrid[2]}, see [1]) in every disp-* folder.",
        'yellow'), f_out)

    print_section('[6] LIBRARY WARNINGS', f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (MACE/torch/phonopy) during this run -- "
            "collected here instead of interleaved with the report above; harmless in almost "
            "every case (import-time notices, deprecation-style warnings), but worth a glance.",
            'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n"+"-"*60)

if __name__ == "__main__":
    main()
