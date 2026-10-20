#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import sys
import shutil
import argparse
from datetime import datetime

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io
from stb.core.adsorption_sites import SINGLE_POINT_BLOCK, CONFIG_EXTRA_FILE
from stb.core.pseudopotentials import get_required_pseudos
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core import mlff_workflow

REPORT_FILE = "stb_mlff_report.txt"

# Every training folder is ONE fixed sample point: SIESTA must evaluate it, not
# move it. A relaxation template (MD.TypeOfRun CG + MD.Steps > 0) would drift the
# geometry towards equilibrium, and MD.VariableCell T would relax a strained cell
# back to the reference volume -- either way the labels end up describing a
# different structure than the one sampled. Forced through config_extra.fdf
# (first-occurrence-wins, see structure_io.prepend_include), same mechanism as
# the other workflows' single-point folders.
_MLFF_EXTRA_BLOCK = (
    SINGLE_POINT_BLOCK
    + "MD.NumCGsteps         0\n"
    + "MD.VariableCell       false\n"
)


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


def build_strained_configs(atoms, coord_format, raw_lattice, strains_pct, stdevs, seed):
    """Returns one (positions, raw_lattice) pair per entry of `strains_pct`
    (volumetric strain, %): the cell scaled isotropically by
    (1 + strain/100)**(1/3) -- same convention as stb-eosInputs/stb-mleos --
    and the atoms carried with it (same fractional coordinates), then rattled
    with the `stdevs` amplitudes cycled round-robin.

    Rattled rather than left ideal: at ideal positions of a high-symmetry
    crystal every force vanishes by symmetry, so such a folder would teach
    the model energy and stress only. Rattled, it also adds forces at a
    volume the position-only configurations never reach. Without at least a
    few of these, the training set has a single volume and the fitted model
    has nothing to learn the volume response (E(V), bulk modulus, pressure)
    from. `raw_lattice` is the reference cell in the file's own
    LatticeConstant basis (structure_io.raw_lattice_vectors).
    """
    configs = []
    scaled0 = atoms.get_scaled_positions(wrap=False)
    for i, strain in enumerate(strains_pct):
        factor = (1.0 + strain / 100.0) ** (1.0 / 3.0)
        a = atoms.copy()
        a.set_cell(atoms.get_cell() * factor, scale_atoms=False)
        a.set_scaled_positions(scaled0)
        a.rattle(stdev=stdevs[i % len(stdevs)], seed=seed + 10000 + i)
        if coord_format == "fractional":
            positions = a.get_scaled_positions(wrap=False)
        else:
            positions = a.get_positions()
        configs.append((positions, np.asarray(raw_lattice) * factor))
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


TEST_SEED_OFFSET = 777000
TEST_DIR = "mlff_test"


def build_test_configs(atoms, coord_format, raw_lattice, n_test, stdevs, strains_pct, seed):
    """`n_test` held-out configurations drawn like the training ones but with an
    independent seed: the kinds -- one rattle per --stdev, one strained cell per
    --volume-strains entry -- are cycled round-robin. Returns a list of
    (positions, raw_lattice or None, kind, origin). The point is a test set
    from the same distribution the model is meant for, never trained on, so
    stb-mlffValidate can tell how a round really does and whether the next one
    helped."""
    kinds = [("rattle", s) for s in stdevs] + [("strain", p) for p in strains_pct]
    out = []
    for i in range(n_test):
        kind, value = kinds[i % len(kinds)]
        stdev = stdevs[i % len(stdevs)]
        if kind == "rattle":
            pos = build_rattled_configs(atoms, coord_format, 1, [value], seed + TEST_SEED_OFFSET + i)[0]
            out.append((pos, None, "rattle", f"rattle {value:g} Ang"))
        else:
            pos, lat = build_strained_configs(atoms, coord_format, raw_lattice, [value], [stdev],
                                              seed + TEST_SEED_OFFSET + i)[0]
            out.append((pos, lat, "strain", f"volume {value:+g} % (rattle {stdev:g} Ang)"))
    return out


def write_config_folder(output_root, index, ref_fdf, ref_basename, positions,
                         raw_lattice, pseudos_to_copy, extra_block="", prefix="mlff_config"):
    """Writes mlff_config_NNN/<ref_basename> (positions rewritten via
    rewrite_fdf_positions, lattice ALSO rewritten via rewrite_fdf_lattice
    when raw_lattice is given -- e.g. a variable-cell AIMD snapshot) plus
    the copied pseudopotentials, ready for the user to run SIESTA in, plus a
    config_extra.fdf forcing a single-point evaluation (%include'd at the top
    of the .fdf; see _MLFF_EXTRA_BLOCK). `extra_block` (e.g. a supercell's own
    k-grid, stb-mlffActiveLearning --supercell) is written first in that file,
    so it wins over the template. `prefix` names the folder (<prefix>_NNN;
    the held-out test set uses "test"). Returns the folder path.
    """
    folder = os.path.join(output_root, f"{prefix}_{index:03d}")
    os.makedirs(folder, exist_ok=True)
    out_fdf = os.path.join(folder, ref_basename)

    structure_io.rewrite_fdf_positions(ref_fdf, positions, out_fdf)
    if raw_lattice is not None:
        structure_io.rewrite_fdf_lattice(out_fdf, raw_lattice, out_fdf)

    with open(os.path.join(folder, CONFIG_EXTRA_FILE), "w") as f:
        f.write(extra_block + _MLFF_EXTRA_BLOCK)
    with open(out_fdf) as f:
        fdf_text = f.read()
    if f"%include {CONFIG_EXTRA_FILE}" not in fdf_text:
        with open(out_fdf, "w") as f:
            f.write(structure_io.prepend_include(fdf_text, CONFIG_EXTRA_FILE))

    for pseudo_path in pseudos_to_copy:
        shutil.copy(pseudo_path, os.path.join(folder, os.path.basename(pseudo_path)))

    return folder


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 (Prep) of the Custom ML Force Field workflow: generates "
                     "training configurations (rattled displacements and/or AIMD-sampled "
                     "snapshots) from a reference SIESTA calculation, ready for the user to "
                     "run through SIESTA, plus an optional held-out test set (--n-test). "
                     "Creates mlff_workflow.json, through which the next stages -- Train "
                     "(stb-mlffAnalysis), Validate (stb-mlffValidate), Active learning "
                     "(stb-mlffActiveLearning) -- find what the earlier ones did.",
        epilog="Example usage:\n"
               "  stb-mlff --file calc.fdf --n-configs 10 --stdev 0.05 0.1\n"
               "  stb-mlff --file calc.fdf --n-configs 30 --stdev 0.05 0.1 --volume-strains -5 -3 -1 1 3 5\n"
               "  stb-mlff --file calc.fdf --n-configs 6 --from-aimd aimd --n-samples 6\n"
               "  stb-mlff --file calc.fdf --n-configs 30 --stdev 0.05 0.1 --volume-strains -3 3 --n-test 8",
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
    parser.add_argument("--volume-strains", type=float, nargs="+", default=[], metavar="PCT",
                        help="Also write one configuration per volumetric strain, in %% (e.g. "
                             "-5 -3 -1 1 3 5): the cell scaled isotropically by "
                             "(1 + PCT/100)^(1/3), atoms rattled with --stdev. Without them "
                             "every configuration has the reference volume and the fitted "
                             "model cannot learn the volume response (E(V), bulk modulus, "
                             "pressure). Default: none. Bulk only -- in a slab this also "
                             "scales the vacuum.")
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
    parser.add_argument("--n-test", type=int, default=0,
                        help="Also write N held-out test configurations to "
                             f"{TEST_DIR}/test_NNN/ (default: 0): drawn like the training ones "
                             "(each --stdev, each --volume-strains entry, cycled) with an "
                             "independent seed, never trained on. Stage 3 (stb-mlffValidate) "
                             "scores every round against them -- needed to know when to stop.")
    parser.add_argument("-o", "--output-dir", default=".",
                        help="Directory under which mlff_config_NNN/ folders are created "
                             "(default: current directory).")
    parser.add_argument("--workflow", default=None, metavar="JSON",
                        help=f"Workflow state file shared by the 4 stages (default: "
                             f"<output-dir>/{mlff_workflow.WORKFLOW_FILE}). Created here; an "
                             "existing one is copied to .bak first.")
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
    if args.volume_strains:
        print_dual(f"Volume strains    : {', '.join(f'{x:+g}' for x in args.volume_strains)} %", f_out)
    if args.from_aimd:
        print_dual(f"AIMD samples      : {args.n_samples} (from '{args.from_aimd}.ANI')", f_out)
    print_dual(f"Test configs      : {args.n_test}" + (f" (in {TEST_DIR}/)" if args.n_test else ""), f_out)
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

    if args.n_test < 0:
        fail(f"--n-test must be >= 0 (got {args.n_test}).")
    workflow_path = args.workflow or os.path.join(args.output_dir, mlff_workflow.WORKFLOW_FILE)
    wf_base = mlff_workflow.base_dir(workflow_path)

    print_section("[2] GENERATING CONFIGURATIONS", f_out)
    folders = []
    entries = []
    index = 0
    n_ref = len(ref_structure.atoms)

    def _entry(folder, kind, origin, which="train"):
        entries.append({"folder": folder, "set": which, "kind": kind, "origin": origin,
                        "stage": "prep", "round": 0, "n_atoms": n_ref})

    if args.n_configs > 0:
        rattled = build_rattled_configs(
            ase_atoms, ref_structure.coord_format, args.n_configs, args.stdev, args.seed)
        for i, positions in enumerate(rattled):
            index += 1
            folder = write_config_folder(
                args.output_dir, index, args.file, ref_basename, positions, None, pseudos_to_copy)
            folders.append(folder)
            _entry(folder, "rattle", f"rattle {args.stdev[i % len(args.stdev)]:g} Ang")
        print_dual(f"[OK] Generated {len(rattled)} rattled configuration(s).", f_out)

    if args.volume_strains:
        bad = [x for x in args.volume_strains if x <= -100.0]
        if bad:
            fail(f"--volume-strains must be > -100 % (got {bad}).")
        strained = build_strained_configs(
            ase_atoms, ref_structure.coord_format,
            structure_io.raw_lattice_vectors(ref_structure),
            args.volume_strains, args.stdev, args.seed)
        for i, (positions, raw_lattice) in enumerate(strained):
            index += 1
            folder = write_config_folder(
                args.output_dir, index, args.file, ref_basename, positions, raw_lattice, pseudos_to_copy)
            folders.append(folder)
            _entry(folder, "strain", f"volume {args.volume_strains[i]:+g} % "
                                     f"(rattle {args.stdev[i % len(args.stdev)]:g} Ang)")
        print_dual(f"[OK] Generated {len(strained)} strained configuration(s) "
                   f"(volume {', '.join(f'{x:+g}' for x in args.volume_strains)} %).", f_out)

    if args.from_aimd:
        try:
            aimd_configs = sample_aimd_configs(
                args.from_aimd, ref_structure.coord_format, ref_structure.lattice_constant,
                args.n_samples, f_out)
        except Exception as e:
            fail(f"Could not sample from AIMD trajectory '{args.from_aimd}': {e}")
        for i, (positions, raw_lattice) in enumerate(aimd_configs):
            index += 1
            folder = write_config_folder(
                args.output_dir, index, args.file, ref_basename, positions, raw_lattice, pseudos_to_copy)
            folders.append(folder)
            _entry(folder, "aimd", f"AIMD frame {i + 1}/{len(aimd_configs)} of '{args.from_aimd}'")
        print_dual(f"[OK] Generated {len(aimd_configs)} AIMD-sampled configuration(s).", f_out)

    if not folders:
        fail("No configurations were generated (--n-configs 0, no --volume-strains and no --from-aimd).")

    test_folders = []
    if args.n_test:
        tests = build_test_configs(
            ase_atoms, ref_structure.coord_format, structure_io.raw_lattice_vectors(ref_structure),
            args.n_test, args.stdev, args.volume_strains, args.seed)
        test_root = os.path.join(args.output_dir, TEST_DIR)
        os.makedirs(test_root, exist_ok=True)
        for i, (positions, raw_lattice, kind, origin) in enumerate(tests, start=1):
            folder = write_config_folder(test_root, i, args.file, ref_basename, positions,
                                         raw_lattice, pseudos_to_copy, prefix="test")
            test_folders.append(folder)
            _entry(folder, kind, origin, which="test")
        print_dual(f"[OK] Generated {len(test_folders)} held-out test configuration(s) in "
                   f"{test_root}/ (never trained on).", f_out)

    print_dual(f"[OK] Each folder forces a single-point run via {CONFIG_EXTRA_FILE} "
               "(no relaxation, no cell change).", f_out)

    bak = mlff_workflow.backup(workflow_path)
    if bak:
        print_dual(color_text(f"[WARNING] {workflow_path} already existed -- copied to {bak} "
                              "and started fresh.", 'yellow'), f_out)
    state = mlff_workflow.new_state(args.file, n_ref, unique_elements, wf_base)
    # A folder is stored relative to the workflow file; a bundled bank name
    # (dojo, virtual_vault) as given.
    state["reference"]["pseudo_dir"] = (mlff_workflow.rel(args.pseudo_dir, wf_base)
                                        if os.path.isdir(args.pseudo_dir) else args.pseudo_dir)
    mlff_workflow.add_configs(state, entries, wf_base)
    mlff_workflow.record(state, "prep", "stb-mlff", sys.argv[1:])
    mlff_workflow.save(state, workflow_path)

    print_section("[3] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Configurations    : {len(folders)} training" +
               (f" + {len(test_folders)} test" if test_folders else ""), f_out)
    for folder in folders + test_folders:
        print_dual(f"  - {folder}", f_out)
    print_dual(f"Workflow file     : {workflow_path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Run SIESTA in each mlff_config_NNN/ folder" +
                     (f" and in {TEST_DIR}/test_NNN/" if test_folders else "") +
                     ", then train with stb-mlffAnalysis (4.17.2).\n", 'bold'))


if __name__ == "__main__":
    main()
