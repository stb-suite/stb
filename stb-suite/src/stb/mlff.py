#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import shutil
import argparse
from datetime import datetime

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io
from stb.core.pseudopotentials import get_required_pseudos
from stb.core.cli import color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_mlff_report.txt"


def build_rattled_configs(atoms, coord_format, n_configs, stdevs, seed):
    """Returns a list of `n_configs` position arrays (same coordinate system
    as `coord_format`, matching the reference calc.fdf's own
    AtomicCoordinatesFormat), each a Gaussian-displaced ("rattled") copy of
    `atoms`'s equilibrium positions. Cycles through `stdevs` (one or more
    displacement amplitudes) round-robin across the `n_configs` -- using
    several amplitudes instead of just one is standard MLIP training-set
    practice, so the fitted potential sees more than one noise scale.
    """
    configs = []
    for i in range(n_configs):
        stdev = stdevs[i % len(stdevs)]
        a = atoms.copy()
        a.rattle(stdev=stdev, seed=seed + i)
        if coord_format == "fractional":
            configs.append(a.get_scaled_positions(wrap=False))
        else:
            configs.append(a.get_positions())
    return configs


def sample_aimd_configs(label, coord_format, lattice_constant, n_samples, f_out=None):
    """Returns a list of (positions, raw_lattice_or_None) pairs, one per
    sampled AIMD frame -- `n_samples` frames spaced evenly across
    `<label>.ANI`. raw_lattice is None when the frame's cell matches the
    reference structure closely enough not to need rewriting (same
    LatticeConstant-normalized convention rewrite_fdf_lattice expects,
    i.e. already divided by `lattice_constant`), otherwise the frame's own
    cell (for a variable-cell/NPT run, so each sampled snapshot keeps its
    real cell instead of silently reusing the reference one).
    """
    from stb.core.deps import require_sisl
    sisl = require_sisl()
    from stb.core.md_traj import read_frame_lattices

    ani_file = f"{label}.ANI"
    if not os.path.exists(ani_file):
        raise FileNotFoundError(f"'{ani_file}' not found.")

    sile = sisl.get_sile(ani_file)
    all_frames = sile.read_geometry[:]()
    if not all_frames:
        raise ValueError(f"'{ani_file}' contains no frames.")

    all_cells, _ = read_frame_lattices(label, len(all_frames), f_out)
    if all_cells is None:
        raise ValueError(
            f"Could not find a readable '{label}.out', '{label}.XV' or '{label}.fdf' "
            f"-- a lattice is required ('{ani_file}' doesn't carry one of its own).")

    n_samples = min(n_samples, len(all_frames))
    indices = np.linspace(0, len(all_frames) - 1, n_samples).round().astype(int)
    indices = sorted(set(indices.tolist()))

    configs = []
    for idx in indices:
        cell = np.array(all_cells[idx])
        cart = all_frames[idx].xyz
        if coord_format == "fractional":
            pos = cart @ np.linalg.inv(cell)
        else:
            pos = cart
        raw_lattice = cell / lattice_constant
        configs.append((pos, raw_lattice))
    return configs


def write_config_folder(output_root, index, ref_fdf, ref_basename, positions,
                         raw_lattice, pseudos_to_copy):
    """Writes mlff_config_NNN/<ref_basename> (positions rewritten via
    rewrite_fdf_positions, lattice ALSO rewritten via rewrite_fdf_lattice
    when raw_lattice is given -- e.g. a variable-cell AIMD snapshot) plus
    the copied pseudopotentials, ready for the user to run SIESTA in.
    Returns the folder path.
    """
    folder = os.path.join(output_root, f"mlff_config_{index:03d}")
    os.makedirs(folder, exist_ok=True)
    out_fdf = os.path.join(folder, ref_basename)

    structure_io.rewrite_fdf_positions(ref_fdf, positions, out_fdf)
    if raw_lattice is not None:
        structure_io.rewrite_fdf_lattice(out_fdf, raw_lattice, out_fdf)

    for pseudo_path in pseudos_to_copy:
        shutil.copy(pseudo_path, os.path.join(folder, os.path.basename(pseudo_path)))

    return folder


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 (Prep) of the Custom ML Force Field workflow: generates "
                     "training configurations (rattled displacements and/or AIMD-sampled "
                     "snapshots) from a reference SIESTA calculation, ready for the user to "
                     "run through SIESTA. Stage 2 (stb-mlffAnalysis) aggregates the finished "
                     "calculations and fine-tunes a MACE-MP foundation model on them.",
        epilog="Example usage:\n"
               "  stb-mlff --file calc.fdf --n-configs 10 --stdev 0.05 0.1\n"
               "  stb-mlff --file calc.fdf --n-configs 6 --from-aimd aimd --n-samples 6",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True,
                        help="Reference SIESTA .fdf (basis/mesh cutoff/pseudos/SCF settings "
                             "are reused verbatim in every generated configuration).")
    parser.add_argument("--n-configs", type=int, default=10,
                        help="Number of rattled configurations to generate (default: 10). "
                             "Ignored if --from-aimd is given without also wanting rattled "
                             "configs (see --n-configs 0 to disable rattling entirely).")
    parser.add_argument("--stdev", type=float, nargs="+", default=[0.05],
                        help="One or more rattling displacement standard deviations, in "
                             "Angstrom (default: 0.05). Cycled round-robin across --n-configs.")
    parser.add_argument("--pseudo-dir", default=".",
                        help="Directory to search for <element>.psf/.psml pseudopotentials "
                             "(default: current directory).")
    parser.add_argument("--from-aimd", default=None, metavar="LABEL",
                        help="Also sample starting geometries from an existing AIMD "
                             "trajectory ('<LABEL>.ANI' + '.out'), in addition to (or instead "
                             "of, with --n-configs 0) rattled configurations -- more "
                             "physically realistic diversity than rattling alone.")
    parser.add_argument("--n-samples", type=int, default=5,
                        help="Number of evenly-spaced AIMD frames to sample when --from-aimd "
                             "is given (default: 5).")
    parser.add_argument("-o", "--output-dir", default=".",
                        help="Directory under which mlff_config_NNN/ folders are created "
                             "(default: current directory).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for rattling (default: 42).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlff {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Prep)",
            "Generate SIESTA training configurations for MACE fine-tuning",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CUSTOM ML FORCE FIELD -- PREP:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLFF REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Reference file    : {args.file}", f_out)
    print_dual(f"Rattled configs   : {args.n_configs} (stdev {args.stdev})", f_out)
    if args.from_aimd:
        print_dual(f"AIMD samples      : {args.n_samples} (from '{args.from_aimd}.ANI')", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING REFERENCE STRUCTURE", f_out)

    if not os.path.isfile(args.file):
        fail(f"Reference file '{args.file}' not found.")

    try:
        ref_structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")

    unique_elements = sorted(ref_structure.species)
    print_dual(f"[OK] Read {len(ref_structure.atoms)} atom(s), species: {', '.join(unique_elements)}", f_out)

    pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)
    if missing:
        fail(f"Missing pseudopotential(s) for: {', '.join(missing)} (searched '{args.pseudo_dir}').")
    print_dual(f"[OK] Pseudopotentials found: {', '.join(os.path.basename(p) for p in pseudos_to_copy)}", f_out)

    pmg_structure = structure_io.to_pymatgen(ref_structure)
    ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    ref_basename = os.path.basename(args.file)

    os.makedirs(args.output_dir, exist_ok=True)

    print_section("[2] GENERATING CONFIGURATIONS", f_out)
    folders = []
    index = 0

    if args.n_configs > 0:
        rattled = build_rattled_configs(
            ase_atoms, ref_structure.coord_format, args.n_configs, args.stdev, args.seed)
        for positions in rattled:
            index += 1
            folder = write_config_folder(
                args.output_dir, index, args.file, ref_basename, positions, None, pseudos_to_copy)
            folders.append(folder)
        print_dual(f"[OK] Generated {len(rattled)} rattled configuration(s).", f_out)

    if args.from_aimd:
        try:
            aimd_configs = sample_aimd_configs(
                args.from_aimd, ref_structure.coord_format, ref_structure.lattice_constant,
                args.n_samples, f_out)
        except Exception as e:
            fail(f"Could not sample from AIMD trajectory '{args.from_aimd}': {e}")
        for positions, raw_lattice in aimd_configs:
            index += 1
            folder = write_config_folder(
                args.output_dir, index, args.file, ref_basename, positions, raw_lattice, pseudos_to_copy)
            folders.append(folder)
        print_dual(f"[OK] Generated {len(aimd_configs)} AIMD-sampled configuration(s).", f_out)

    if not folders:
        fail("No configurations were generated (both --n-configs 0 and no --from-aimd).")

    print_section("[3] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Configurations    : {len(folders)}", f_out)
    for folder in folders:
        print_dual(f"  - {folder}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Run SIESTA in each mlff_config_NNN/ folder, then use stb-mlffAnalysis.\n", 'bold'))


if __name__ == "__main__":
    main()
