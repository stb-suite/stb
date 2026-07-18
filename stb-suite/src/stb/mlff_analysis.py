#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import glob
import shutil
import subprocess
import argparse
from datetime import datetime

import numpy as np
from ase import Atoms
from ase.io import write as ase_write
from ase.calculators.calculator import CalculationFailed

from stb.core import structure_io, siesta_log
from stb.core.deps import require_sisl
from stb.core.cli import color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_mlffAnalysis_report.txt"


def read_fa_forces(fa_path):
    """Reads ALL atomic forces (eV/Ang) from a SIESTA .FA file -- format:
    first line = atom count, then one row per atom '<1-based index> Fx Fy
    Fz'. Generalizes her_analysis.py's read_fa_force() (which only reads
    ONE specific atom, for HER's local-ZPE Hessian) to the whole array;
    first consumer of the "all atoms at once" shape, so this stays local
    here rather than moving into core/ yet (extract-on-second-use policy).
    Returns an (natoms, 3) array, or None on any read/parse failure.
    """
    try:
        with open(fa_path) as f:
            lines = f.readlines()
    except OSError:
        return None
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError):
        return None
    if len(lines) < n + 1:
        return None
    try:
        forces = np.array([[float(x) for x in lines[i + 1].split()[1:4]] for i in range(n)])
    except (IndexError, ValueError):
        return None
    return forces


def collect_training_configs(config_dirs, f_out=None):
    """Scans each config_dirs entry for a finished SIESTA calculation
    (SystemLabel auto-detected from its .fdf) and builds one ase.Atoms per
    folder, with REF_energy (info scalar), REF_forces (per-atom array), and
    REF_stress (info, 3x3 in eV/Ang^3 -- best-effort, see the caller's own
    note on this) when available. Folders missing any of .out/.FA/.XV are
    skipped with a warning, never aborting the whole batch.

    Geometry (positions + cell) is read from the folder's own <label>.XV --
    NOT the original reference structure -- since SIESTA relaxes/updates it
    even in what's nominally a single-point run's internal SCF geometry
    handling, and because an AIMD-sampled config may carry a different cell
    than the reference altogether.

    Returns (atoms_list, n_skipped, stress_available).
    """
    sisl = require_sisl()
    atoms_list = []
    n_skipped = 0
    stress_available = True

    for d in sorted(config_dirs):
        label = structure_io.find_fdf_system_label(d)
        if not label:
            print_dual(color_text(f"[WARNING] {d}: no SystemLabel found -- skipping.", 'yellow'), f_out)
            n_skipped += 1
            continue

        out_file = siesta_log.find_out_file(d, label)
        fa_file = os.path.join(d, f"{label}.FA")
        xv_file = os.path.join(d, f"{label}.XV")

        if not (out_file and os.path.isfile(fa_file) and os.path.isfile(xv_file)):
            print_dual(color_text(
                f"[WARNING] {d}: missing .out/.FA/.XV (calculation not finished?) -- skipping.",
                'yellow'), f_out)
            n_skipped += 1
            continue

        energy = siesta_log.get_free_energy(out_file)
        forces = read_fa_forces(fa_file)
        if energy is None or forces is None:
            print_dual(color_text(
                f"[WARNING] {d}: could not parse energy/forces -- skipping.", 'yellow'), f_out)
            n_skipped += 1
            continue

        try:
            geom = sisl.get_sile(xv_file).read_geometry()
        except Exception as e:
            print_dual(color_text(f"[WARNING] {d}: could not read '{xv_file}': {e} -- skipping.", 'yellow'), f_out)
            n_skipped += 1
            continue

        if len(geom.xyz) != len(forces):
            print_dual(color_text(
                f"[WARNING] {d}: atom count mismatch between .XV ({len(geom.xyz)}) and "
                f".FA ({len(forces)}) -- skipping.", 'yellow'), f_out)
            n_skipped += 1
            continue

        atoms = Atoms(symbols=[a.symbol for a in geom.atoms], positions=geom.xyz,
                      cell=np.array(geom.cell), pbc=True)
        atoms.info['REF_energy'] = energy
        atoms.new_array('REF_forces', forces)

        stress = siesta_log.get_stress_tensor(out_file)
        if stress is not None:
            atoms.info['REF_stress'] = stress
        else:
            stress_available = False

        atoms_list.append(atoms)

    return atoms_list, n_skipped, stress_available and bool(atoms_list)


def find_mace_model(work_dir, name):
    """Best-effort glob for the fine-tuned model file(s) mace_run_train
    wrote under `work_dir` (its exact default subfolder layout isn't a
    stable contract across mace-torch versions, so this globs broadly by
    the `--name` prefix instead of hard-coding one path) and returns the
    most recently modified '<name>*.model' match, or None if none found.
    """
    matches = glob.glob(os.path.join(work_dir, "**", f"{name}*.model"), recursive=True)
    matches += glob.glob(os.path.join(work_dir, f"{name}*.model"))
    if not matches:
        return None
    return max(set(matches), key=os.path.getmtime)


def main():
    parser = argparse.ArgumentParser(
        description="Stage 2 (Analysis) of the Custom ML Force Field workflow: aggregates "
                     "the finished SIESTA calculations stb-mlff generated (energy + forces, "
                     "best-effort stress), builds a MACE-format training set, and fine-tunes "
                     "a MACE-MP foundation model on it via 'mace_run_train' -- specializing "
                     "the generic foundation potential to the user's own material using real "
                     "SIESTA data.",
        epilog="Example usage:\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --epochs 100 --device cpu\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --foundation-model medium --device cuda",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--path", default="mlff_config_*",
                        help="Glob pattern matching the stb-mlff-generated config folders "
                             "(default: 'mlff_config_*').")
    parser.add_argument("--foundation-model", default="small", choices=["small", "medium", "large"],
                        help="MACE-MP foundation checkpoint to fine-tune from (default: small "
                             "-- fastest, appropriate for the typically small dataset size of "
                             "a SIESTA-based fine-tuning run).")
    parser.add_argument("--epochs", type=int, default=100, help="Max training epochs (default: 100).")
    parser.add_argument("--batch-size", type=int, default=5, help="Training batch size (default: 5).")
    parser.add_argument("--lr", type=float, default=0.01, help="Learning rate (default: 0.01).")
    parser.add_argument("--valid-fraction", type=float, default=0.1,
                        help="Fraction of configurations held out for validation (default: 0.1).")
    parser.add_argument("--device", default=None, choices=["cpu", "cuda"],
                        help="Training device (default: auto-detect, cuda if available).")
    parser.add_argument("--name", default="mlff_model", help="Experiment/model name (default: mlff_model).")
    parser.add_argument("--seed", type=int, default=123, help="Training random seed (default: 123).")
    parser.add_argument("-o", "--work-dir", default=".",
                        help="Working directory for the training run and its outputs (default: current directory).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlffAnalysis {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Analysis)",
            "Fine-tune MACE-MP on real SIESTA data",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CUSTOM ML FORCE FIELD -- ANALYSIS:", 'bold'))
    print("-" * 60)

    if shutil.which("mace_run_train") is None:
        print(color_text("\n[CRITICAL ERROR] 'mace_run_train' not found on PATH.", 'red'))
        print("This tool needs the optional 'ml' extra. Install it with:")
        print("  pip install stb_suite[ml]")
        sys.exit(1)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLFFANALYSIS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Config pattern    : {args.path}", f_out)
    print_dual(f"Foundation model  : mace-mp-0 ({args.foundation_model})", f_out)
    print_dual(f"Epochs / batch    : {args.epochs} / {args.batch_size}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] AGGREGATING TRAINING DATA", f_out)

    config_dirs = sorted(d for d in glob.glob(args.path) if os.path.isdir(d))
    if not config_dirs:
        fail(f"No directories matched '{args.path}'.")
    print_dual(f"Config folders    : {len(config_dirs)} matched", f_out)

    atoms_list, n_skipped, stress_ok = collect_training_configs(config_dirs, f_out)
    if len(atoms_list) < 4:
        fail(f"Only {len(atoms_list)} usable configuration(s) found -- need at least a few "
             "for a train/validation split. Check that SIESTA finished in each folder.")

    print_dual(f"[OK] {len(atoms_list)} usable configuration(s) ({n_skipped} skipped).", f_out)
    if stress_ok:
        print_dual("[OK] Stress available for every configuration -- will be used as an "
                   "additional (best-effort) training target.", f_out)
    else:
        print_dual(color_text(
            "[INFO] Stress not available for every configuration -- training on "
            "energy+forces only.", 'yellow'), f_out)

    os.makedirs(args.work_dir, exist_ok=True)
    dataset_path = os.path.join(args.work_dir, "training_set.xyz")
    ase_write(dataset_path, atoms_list, format='extxyz')
    print_dual(f"[OK] Wrote training set to {dataset_path}", f_out)

    print_section("[2] FINE-TUNING (mace_run_train)", f_out)

    device = args.device or ("cuda" if _cuda_available() else "cpu")
    print_dual(f"Device            : {device}", f_out)

    cmd = [
        "mace_run_train",
        "--name", args.name,
        "--work_dir", args.work_dir,
        "--train_file", dataset_path,
        "--valid_fraction", str(args.valid_fraction),
        "--energy_key", "REF_energy",
        "--forces_key", "REF_forces",
        "--foundation_model", args.foundation_model,
        # Plain fine-tuning (--multiheads_finetuning explicitly False --
        # mace_run_train defaults it to True), with "average"
        # (per-element atomic reference energies least-squares-fit from OUR
        # OWN training set). Deliberate, verified live: SIESTA's absolute
        # total energy is on a different scale than the foundation model's
        # (different code, different pseudopotentials/core treatment), so
        # the only E0s mode that produces a sane, converging loss on SIESTA
        # data is "average" -- but mace_run_train's --multiheads_finetuning
        # (the officially recommended anti-catastrophic-forgetting mode)
        # structurally REQUIRES --E0s foundation and asserts against
        # "average" outright. Tried "foundation" first: RMSE_E_per_atom sat
        # at ~424000 meV/atom (the SIESTA/foundation absolute-scale offset)
        # and never improved. Trade-off accepted: no pretraining-replay
        # protection against forgetting the foundation's broad knowledge,
        # in exchange for training targets that are actually on the right
        # scale -- appropriate here since the goal is a potential
        # specialized to the user's own material, not a broadly-capable one.
        "--E0s", "average",
        "--multiheads_finetuning", "False",
        "--max_num_epochs", str(args.epochs),
        "--batch_size", str(args.batch_size),
        "--lr", str(args.lr),
        "--seed", str(args.seed),
        "--device", device,
    ]
    if stress_ok:
        cmd.extend(["--stress_key", "REF_stress", "--compute_stress", "True"])

    print_dual(f"Command           : {' '.join(cmd)}", f_out)
    log_path = os.path.join(args.work_dir, "mace_run_train.log")
    try:
        with open(log_path, "w") as log_f:
            result = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT)
    except FileNotFoundError:
        fail("'mace_run_train' could not be executed (not found on PATH at run time).")

    if result.returncode != 0:
        print_dual(color_text(
            f"[FAIL] mace_run_train exited with code {result.returncode} -- see {log_path} "
            "for the full training log.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"[OK] mace_run_train finished -- full log at {log_path}", f_out)

    print_section("[3] SUMMARY & FILES", f_out)
    model_path = find_mace_model(args.work_dir, args.name)
    if model_path:
        print_dual(f"[OK] Fine-tuned model: {model_path}", f_out)
        smoke_test_model(model_path, atoms_list[0], device, f_out)
    else:
        print_dual(color_text(
            "[WARNING] Training finished but no '<name>*.model' file was found under "
            f"'{args.work_dir}' -- check the log for the model's actual output path.",
            'yellow'), f_out)

    print_dual("Status            : OK", f_out)
    print_dual(f"Training set      : {dataset_path}", f_out)
    print_dual(f"Training log      : {log_path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Custom ML force field ready -- fast, and fit to your own data.\n", 'bold'))


def _cuda_available():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


def smoke_test_model(model_path, sample_atoms, device, f_out=None):
    """Loads the fine-tuned model as an ASE calculator and evaluates one
    reference configuration's energy -- a minimal "does the model actually
    load and run" check, not an accuracy validation (that's what
    mace_run_train's own validation-set RMSE, in the training log, is for).
    """
    try:
        from mace.calculators import MACECalculator
        calc = MACECalculator(model_paths=[model_path], device=device)
        atoms = sample_atoms.copy()
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        print_dual(f"[OK] Model loads and runs -- sample energy: {energy:.6f} eV", f_out)
    except (ImportError, CalculationFailed, Exception) as e:
        print_dual(color_text(f"[WARNING] Could not smoke-test the fine-tuned model: {e}", 'yellow'), f_out)


if __name__ == "__main__":
    main()
