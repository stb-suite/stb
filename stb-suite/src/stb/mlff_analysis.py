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
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import write as ase_write
from ase.calculators.calculator import CalculationFailed

from stb.core import structure_io, siesta_log
from stb.core.deps import require_sisl
from stb.core.cli import color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_mlffAnalysis_report.txt"


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
        forces = siesta_log.read_fa_forces(fa_file)
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


def split_train_valid(atoms_list, valid_fraction, seed):
    """Deterministic (seeded) train/validation split, done here in Python
    rather than left to mace_run_train's own --valid_fraction: that way we
    know exactly which configurations ended up in validation, and can
    independently evaluate the fine-tuned model (and, for comparison, the
    raw foundation model) on that same held-out set afterwards -- neither
    is possible if mace_run_train does the split internally. At least 1
    configuration is kept on each side regardless of how small
    `valid_fraction` or the dataset is (the >= 4 configs guard upstream
    already ensures at least 3 remain for training).
    """
    n = len(atoms_list)
    n_valid = min(max(1, round(n * valid_fraction)), n - 1)
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    valid_idx = set(order[:n_valid].tolist())
    train = [a for i, a in enumerate(atoms_list) if i not in valid_idx]
    valid = [a for i, a in enumerate(atoms_list) if i in valid_idx]
    return train, valid


def _rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def evaluate_on_configs(calc, atoms_list):
    """Evaluates `calc` on every atoms_list entry's OWN geometry (never
    re-relaxed) and returns (pred_e_per_atom, ref_e_per_atom, pred_forces,
    ref_forces) as flat arrays -- pred/ref pairs in the same order, ready
    for a parity plot or an RMSE. Configurations `calc` fails to evaluate
    (e.g. an element combination outside a small custom model's training
    domain) are silently skipped, not fatal -- this is a diagnostic
    report, not the main training path.
    """
    pred_e, ref_e, pred_f, ref_f = [], [], [], []
    for a in atoms_list:
        atoms = a.copy()
        atoms.calc = calc
        try:
            e_pred = atoms.get_potential_energy()
            f_pred = atoms.get_forces()
        except Exception:
            continue
        pred_e.append(e_pred / len(atoms))
        ref_e.append(a.info['REF_energy'] / len(atoms))
        pred_f.append(f_pred.flatten())
        ref_f.append(a.arrays['REF_forces'].flatten())
    pred_f = np.concatenate(pred_f) if pred_f else np.array([])
    ref_f = np.concatenate(ref_f) if ref_f else np.array([])
    return np.array(pred_e), np.array(ref_e), pred_f, ref_f


def plot_parity(series, out_path):
    """`series` is a list of (label, pred_e, ref_e, pred_f, ref_f) tuples
    (one per model evaluated, e.g. fine-tuned and/or the raw foundation
    model) -- overlays all of them on the same 2-panel energy/force parity
    figure (each with the y=x reference line) rather than one plot per
    model, so the improvement from fine-tuning is directly visible in one
    image.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for label, pred_e, ref_e, pred_f, ref_f in series:
        if len(ref_e):
            axes[0].scatter(ref_e, pred_e, s=18, alpha=0.7, label=label)
        if len(ref_f):
            axes[1].scatter(ref_f, pred_f, s=6, alpha=0.4, label=label)

    for ax, xlabel, ylabel, title in [
        (axes[0], "DFT energy (eV/atom)", "Predicted energy (eV/atom)", "Energy parity"),
        (axes[1], "DFT force (eV/Ang)", "Predicted force (eV/Ang)", "Force parity"),
    ]:
        lims = [*ax.get_xlim(), *ax.get_ylim()]
        lo, hi = min(lims), max(lims)
        ax.plot([lo, hi], [lo, hi], color='gray', linestyle='--', linewidth=0.8, zorder=0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


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
               "  stb-mlffAnalysis --path 'mlff_config_*' --foundation-model medium --device cuda\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --export-lammps",
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
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="Skip evaluating the raw (non-fine-tuned) foundation model on the "
                             "validation set for comparison. On by default (comparison runs "
                             "unless this is passed).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw parity-plot data (.dat) behind the plot. Off by default.")
    parser.add_argument("--export-lammps", action="store_true",
                        help="Also compile the fine-tuned model into a TorchScript file "
                             "loadable by LAMMPS's 'pair_style mace' (via mace_create_lammps_model, "
                             "part of mace-torch) -- '<model>-lammps.pt'. Off by default.")
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

    # Split done HERE (not via mace_run_train's own --valid_fraction) so we
    # know exactly which configurations ended up in validation, and can
    # independently evaluate the fine-tuned model (and the raw foundation
    # model, for comparison) on that same held-out set afterwards -- see
    # split_train_valid()'s docstring.
    train_atoms, valid_atoms = split_train_valid(atoms_list, args.valid_fraction, args.seed)
    train_path = os.path.join(args.work_dir, "train_set.xyz")
    valid_path = os.path.join(args.work_dir, "valid_set.xyz")
    ase_write(train_path, train_atoms, format='extxyz')
    ase_write(valid_path, valid_atoms, format='extxyz')
    print_dual(f"[OK] Wrote {len(train_atoms)} training / {len(valid_atoms)} validation "
               f"configuration(s) (train_set.xyz / valid_set.xyz)", f_out)

    print_section("[2] FINE-TUNING (mace_run_train)", f_out)

    device = args.device or ("cuda" if _cuda_available() else "cpu")
    print_dual(f"Device            : {device}", f_out)

    cmd = [
        "mace_run_train",
        "--name", args.name,
        "--work_dir", args.work_dir,
        "--train_file", train_path,
        "--valid_file", valid_path,
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

    print_section("[3] MODEL DISCOVERY", f_out)
    model_path = find_mace_model(args.work_dir, args.name)
    lammps_model_path = None
    if model_path:
        print_dual(f"[OK] Fine-tuned model: {model_path}", f_out)
        smoke_test_model(model_path, atoms_list[0], device, f_out)
        if args.export_lammps:
            lammps_model_path = export_lammps_model(model_path, f_out)
    else:
        print_dual(color_text(
            "[WARNING] Training finished but no '<name>*.model' file was found under "
            f"'{args.work_dir}' -- check the log for the model's actual output path.",
            'yellow'), f_out)
        if args.export_lammps:
            print_dual(color_text(
                "[WARNING] --export-lammps skipped: no model file to convert.", 'yellow'), f_out)

    parity_plot_path = None
    if model_path:
        print_section("[4] VALIDATION & PARITY", f_out)
        from mace.calculators import MACECalculator
        ft_calc = MACECalculator(model_paths=[model_path], device=device)
        pred_e, ref_e, pred_f, ref_f = evaluate_on_configs(ft_calc, valid_atoms)
        series = [("Fine-tuned", pred_e, ref_e, pred_f, ref_f)]
        if len(ref_e):
            print_dual(f"Fine-tuned RMSE   : E = {_rmse(pred_e, ref_e) * 1000:.2f} meV/atom, "
                       f"F = {_rmse(pred_f, ref_f) * 1000:.2f} meV/Ang "
                       f"({len(valid_atoms)} validation config(s))", f_out)
        else:
            print_dual(color_text(
                "[WARNING] Fine-tuned model could not be evaluated on the validation set.",
                'yellow'), f_out)

        if not args.skip_foundation_comparison:
            from stb.core import mace_relax
            found_calc = mace_relax.get_calculator(model=args.foundation_model, device=device)
            fpred_e, fref_e, fpred_f, fref_f = evaluate_on_configs(found_calc, valid_atoms)
            series.append((f"Foundation ({args.foundation_model})", fpred_e, fref_e, fpred_f, fref_f))
            if len(fref_e):
                print_dual(f"Foundation RMSE   : E = {_rmse(fpred_e, fref_e) * 1000:.2f} meV/atom, "
                           f"F = {_rmse(fpred_f, fref_f) * 1000:.2f} meV/Ang "
                           f"(same {len(valid_atoms)} validation config(s), no fine-tuning)", f_out)

        if len(ref_e) or len(ref_f):
            parity_plot_path = os.path.join(args.work_dir, f"{args.name}_parity.png")
            plot_parity(series, parity_plot_path)
            print_dual(f"[OK] Saved parity plot: {parity_plot_path}", f_out)
            if args.save_data:
                for label, pe, re_, pf, rf in series:
                    tag = label.split()[0].lower()
                    if len(re_):
                        e_path = os.path.join(args.work_dir, f"{args.name}_parity_{tag}_energy.dat")
                        np.savetxt(e_path, np.column_stack([re_, pe]), header="ref_E pred_E (eV/atom)")
                        print_dual(f"[OK] Saved {e_path}", f_out)
                    if len(rf):
                        f_path = os.path.join(args.work_dir, f"{args.name}_parity_{tag}_forces.dat")
                        np.savetxt(f_path, np.column_stack([rf, pf]), header="ref_F pred_F (eV/Ang)")
                        print_dual(f"[OK] Saved {f_path}", f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Training set      : {train_path} / {valid_path}", f_out)
    print_dual(f"Training log      : {log_path}", f_out)
    if parity_plot_path:
        print_dual(f"Parity plot       : {parity_plot_path}", f_out)
    if lammps_model_path:
        print_dual(f"LAMMPS model      : {lammps_model_path}", f_out)
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


def export_lammps_model(model_path, f_out=None):
    """Compiles `model_path` into a TorchScript file LAMMPS's 'pair_style
    mace' can load directly (via mace_create_lammps_model, part of
    mace-torch -- writes '<model_path>-lammps.pt', libtorch format, single
    -head models skip mace_create_lammps_model's own interactive head
    -selection prompt automatically). Requires a LAMMPS build with the MACE
    pair style (https://github.com/ACEsuit/lammps) to actually USE the
    resulting file -- this only produces it. Returns the output path, or
    None on failure (never aborts the whole run: the fine-tuned .model
    itself is already usable via ASE/Python regardless of this step).
    """
    if shutil.which("mace_create_lammps_model") is None:
        print_dual(color_text(
            "[WARNING] --export-lammps skipped: 'mace_create_lammps_model' not found on PATH.",
            'yellow'), f_out)
        return None

    out_path = model_path + "-lammps.pt"
    cmd = ["mace_create_lammps_model", model_path, "--dtype", "float64", "--format", "libtorch"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.isfile(out_path):
        print_dual(color_text(
            f"[WARNING] Could not export a LAMMPS model: {result.stderr.strip().splitlines()[-1] if result.stderr else 'unknown error'}",
            'yellow'), f_out)
        return None

    print_dual(f"[OK] Exported LAMMPS model: {out_path}", f_out)
    return out_path


if __name__ == "__main__":
    main()
