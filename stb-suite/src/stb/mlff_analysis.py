#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

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
from stb.core import mlff_workflow

REPORT_FILE = "stb_mlffAnalysis_report.txt"

# Defaults of the fine-tuning recipe, all three set on a bulk-Si benchmark
# (8-atom cell; 30 rattled + 6 volume-strained SIESTA configs + one
# active-learning round; scored on a held-out SIESTA test set of rattled,
# large-displacement, 800 K MD, E(V) and Gamma-phonon configs; 3 training
# seeds). With the earlier recipe (E0s "average", lr 0.01, stress weight 100)
# the bulk modulus came out 99-111 GPa against SIESTA's 94.4 and the force
# error at large displacements 470-790 meV/Ang; with these defaults, 94.4-96.1
# GPa and 118-146 meV/Ang, lattice parameter within 0.05 %.
#
# lr 0.001: at 0.01 the pretrained model is wrecked in the first 1-2 epochs
# (stress RMSE 1.4 -> 150-300 meV/Ang^3) and the fit never gets back below
# its own starting stress error. Stress weight 1e4: stress is in eV/Ang^3
# (~3e-3 for 0.5 GPa), so even a weight of 100 leaves its term ~2 % of the
# loss next to forces_weight=100 and the stress error stays where it started;
# 1e4 brings it down and fixes the volume response, forces unchanged.
DEFAULT_LR = 0.001
DEFAULT_STRESS_WEIGHT = 10000.0
_EV_A3_TO_GPA = 160.21766208
DEFAULT_COMMITTEE = 3


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


def foundation_residual_e0s(train_atoms, found_calc):
    """Per-element reference energies (eV) that put the foundation model's
    OWN interaction energy on the SIESTA scale from the first epoch:
    least-squares fit of E_SIESTA - E_interaction = sum_Z n_Z * E0_Z over
    the training set, where E_interaction = E_foundation - sum_Z n_Z *
    E0_Z(foundation).

    Why not mace_run_train's "average": it fits E_SIESTA = sum_Z n_Z * E0_Z,
    i.e. assumes a zero interaction energy, so the starting model is off by
    the whole cohesive-scale gap (2.6 eV/atom for Si) and the network has to
    absorb it through its weights, in the first epochs, at the cost of what
    it learned in pretraining. With these E0s the starting energy error is
    the foundation model's own error on energy differences (~24 meV/atom).
    Not "foundation" either: those are the foundation's E0s, on its own
    scale, ~1e5 meV/atom away from SIESTA's absolute energies.

    Returns (e0s, rank_ok, start_rmse, avg_e0s): `e0s` {Z: E0} (None if an
    element is outside the foundation model), `rank_ok` False when the
    training compositions don't determine one E0 per element (every config
    the same stoichiometry, for example -- only that composition's sum is
    then meaningful, same limitation as "average"), `start_rmse` the
    resulting starting energy RMSE (eV/atom), `avg_e0s` what "average"
    would have given, for the report.
    """
    model = found_calc.models[0]
    z_table = [int(z) for z in model.atomic_numbers]
    e0_found = dict(zip(z_table, model.atomic_energies_fn.atomic_energies.detach().cpu().numpy().ravel()))
    elements = sorted({int(z) for a in train_atoms for z in a.numbers})
    if any(z not in e0_found for z in elements):
        return None, False, None, None

    counts = np.array([[np.count_nonzero(a.numbers == z) for z in elements] for a in train_atoms], float)
    e_ref = np.array([a.info['REF_energy'] for a in train_atoms])
    e_int = []
    for a in train_atoms:
        b = a.copy()
        b.calc = found_calc
        e_int.append(b.get_potential_energy() - sum(e0_found[int(z)] for z in a.numbers))
    resid = e_ref - np.array(e_int)

    sol, _, rank, _ = np.linalg.lstsq(counts, resid, rcond=None)
    avg, _, _, _ = np.linalg.lstsq(counts, e_ref, rcond=None)
    n_atoms = counts.sum(axis=1)
    start_rmse = float(np.sqrt(np.mean(((counts @ sol - resid) / n_atoms) ** 2)))
    return ({z: float(v) for z, v in zip(elements, sol)}, rank == len(elements), start_rmse,
            {z: float(v) for z, v in zip(elements, avg)})


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
    ref_forces, pred_stress, ref_stress) as flat arrays -- pred/ref pairs in
    the same order, ready for a parity plot or an RMSE (the stress arrays,
    eV/Ang^3, 9 components per config, only cover configs that carry a
    REF_stress). Configurations `calc` fails to evaluate
    (e.g. an element combination outside a small custom model's training
    domain) are silently skipped, not fatal -- this is a diagnostic
    report, not the main training path.
    """
    pred_e, ref_e, pred_f, ref_f, pred_s, ref_s = [], [], [], [], [], []
    for a in atoms_list:
        atoms = a.copy()
        atoms.calc = calc
        try:
            e_pred = atoms.get_potential_energy()
            f_pred = atoms.get_forces()
            s_pred = atoms.get_stress(voigt=False) if 'REF_stress' in a.info else None
        except Exception:
            continue
        pred_e.append(e_pred / len(atoms))
        ref_e.append(a.info['REF_energy'] / len(atoms))
        pred_f.append(f_pred.flatten())
        ref_f.append(a.arrays['REF_forces'].flatten())
        if s_pred is not None:
            pred_s.append(np.asarray(s_pred).flatten())
            ref_s.append(np.asarray(a.info['REF_stress']).flatten())

    def _flat(chunks):
        return np.concatenate(chunks) if chunks else np.array([])
    return (np.array(pred_e), np.array(ref_e), _flat(pred_f), _flat(ref_f),
            _flat(pred_s), _flat(ref_s))


def energy_metrics(pred_e, ref_e):
    """Returns (rmse, offset, rmse_without_offset), eV/atom. `offset` is the
    mean pred - ref per atom; the offset-free RMSE is what says whether a
    model gets the energy DIFFERENCES between configurations right -- the
    only meaningful energy comparison against the raw foundation model,
    whose absolute energies sit on a different scale than SIESTA's
    (different code, pseudopotentials, reference energies: an offset of
    ~100 eV/atom is normal and says nothing about its accuracy).
    """
    diff = np.asarray(pred_e) - np.asarray(ref_e)
    offset = float(np.mean(diff))
    return _rmse(pred_e, ref_e), offset, float(np.sqrt(np.mean((diff - offset) ** 2)))


def report_model_metrics(label, pred_e, ref_e, pred_f, ref_f, pred_s, ref_s, n_configs,
                         f_out=None, show_absolute=True):
    """Prints one model's validation metrics block."""
    rmse_e, offset, rmse_e_rel = energy_metrics(pred_e, ref_e)
    print_dual(f"{label}:", f_out)
    if show_absolute:
        print_dual(f"  E RMSE (absolute)        : {rmse_e * 1000:9.2f} meV/atom", f_out)
    print_dual(f"  E RMSE (offset removed)  : {rmse_e_rel * 1000:9.2f} meV/atom "
               f"(mean offset {offset * 1000:+.1f} meV/atom)", f_out)
    print_dual(f"  F RMSE                   : {_rmse(pred_f, ref_f) * 1000:9.2f} meV/Ang", f_out)
    if len(ref_s):
        print_dual(f"  Stress RMSE              : {_rmse(pred_s, ref_s) * _EV_A3_TO_GPA:9.3f} GPa", f_out)
    print_dual(f"  ({n_configs} validation config(s))", f_out)


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
            # Centered on each side's own mean: a parity of ABSOLUTE energies
            # puts the foundation model ~100 eV/atom away from the diagonal
            # (scale offset, see energy_metrics) and squeezes every point
            # into one dot. Centered, the panel shows whether the energy
            # differences between configurations are right, for any model.
            axes[0].scatter(1000 * (ref_e - np.mean(ref_e)), 1000 * (pred_e - np.mean(pred_e)),
                            s=18, alpha=0.7, label=label)
        if len(ref_f):
            axes[1].scatter(ref_f, pred_f, s=6, alpha=0.4, label=label)

    for ax, xlabel, ylabel, title in [
        (axes[0], "DFT energy - mean (meV/atom)", "Predicted energy - mean (meV/atom)",
         "Energy parity (offset removed)"),
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
        description="Stage 2 (Train) of the Custom ML Force Field workflow: aggregates "
                     "the finished SIESTA calculations (energy + forces, best-effort stress), "
                     "builds a MACE-format training set, and fine-tunes a MACE-MP foundation "
                     "model on it via 'mace_run_train' -- by default a committee of 3 models "
                     "(different seeds), recorded as a new round in mlff_workflow.json for "
                     "Validate (stb-mlffValidate) and Active learning (stb-mlffActiveLearning).",
        epilog="Example usage:\n"
               "  stb-mlffAnalysis                      # configs from mlff_workflow.json, 3-model committee\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --epochs 100 --device cpu --committee 1\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --foundation-model medium --device cuda\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --export-lammps",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--path", default=None,
                        help="Glob pattern matching the config folders to train on (default: the "
                             "training configs listed in --workflow, else 'mlff_config_*').")
    parser.add_argument("--workflow", default=mlff_workflow.WORKFLOW_FILE, metavar="JSON",
                        help="Workflow state file shared by the 4 stages (default: "
                             f"{mlff_workflow.WORKFLOW_FILE}). Read for the config list, updated "
                             "with this round (committee models and metrics); created if missing.")
    parser.add_argument("--committee", type=int, default=DEFAULT_COMMITTEE,
                        help=f"Number of models to train (default: {DEFAULT_COMMITTEE}), seeds --seed, "
                             "--seed+1, ...: each gets its own train/validation split and "
                             "training seed. Their force disagreement is what stb-mlffActiveLearning "
                             "ranks candidates by, and stb-mlffValidate scores the committee mean. "
                             "1 trains a single model (that score then falls back to the largest "
                             "predicted force).")
    parser.add_argument("--foundation-model", default="small", choices=["small", "medium", "large"],
                        help="MACE-MP foundation checkpoint to fine-tune from (default: small "
                             "-- fastest, appropriate for the typically small dataset size of "
                             "a SIESTA-based fine-tuning run).")
    parser.add_argument("--epochs", type=int, default=100, help="Max training epochs (default: 100).")
    parser.add_argument("--batch-size", type=int, default=5, help="Training batch size (default: 5).")
    parser.add_argument("--lr", type=float, default=DEFAULT_LR,
                        help=f"Learning rate (default: {DEFAULT_LR:g}; 0.01 wrecks the pretrained "
                             "model in the first epochs).")
    parser.add_argument("--e0s", default="foundation-residual", choices=["foundation-residual", "average"],
                        help="Per-element reference energies (default: foundation-residual):\n"
                             "  foundation-residual - fit so that the foundation model's own interaction\n"
                             "                        energy is on the SIESTA scale from the start\n"
                             "  average             - mace_run_train's own fit, which assumes zero interaction\n"
                             "                        energy (starts eV/atom off; the earlier behavior)")
    parser.add_argument("--stress-weight", type=float, default=DEFAULT_STRESS_WEIGHT,
                        help="Weight of the stress term in the training loss, used when every "
                             f"configuration has a SIESTA stress (default: {DEFAULT_STRESS_WEIGHT:g}; "
                             "forces weigh 100, energy 1). Stress is what teaches the model the "
                             "volume response -- pair it with stb-mlff --volume-strains.")
    parser.add_argument("--no-stress-training", action="store_true",
                        help="Train on energy and forces only, even if stress is available.")
    parser.add_argument("--valid-fraction", type=float, default=0.1,
                        help="Fraction of configurations held out for validation (default: 0.1).")
    parser.add_argument("--device", default=None, choices=["cpu", "cuda"],
                        help="Training device (default: auto-detect, cuda if available).")
    parser.add_argument("--name", default="mlff_model", help="Experiment/model name (default: mlff_model).")
    parser.add_argument("--seed", type=int, default=123,
                        help="Training random seed of the first committee member (default: 123).")
    parser.add_argument("-o", "--work-dir", default=None,
                        help="Working directory for this round (default: round_NN/ next to "
                             "--workflow when that file exists, else the current directory). "
                             "With --committee > 1 each member goes to <work-dir>/m_s<seed>/.")
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

    if args.committee < 1:
        fail(f"--committee must be >= 1 (got {args.committee}).")

    # ---- where the configs come from, where this round goes
    try:
        state = mlff_workflow.load(args.workflow)
    except ValueError as e:
        fail(str(e))
    wf_base = mlff_workflow.base_dir(args.workflow)
    if args.path:
        config_dirs = sorted(d for d in glob.glob(args.path) if os.path.isdir(d))
        source = f"--path '{args.path}'"
    elif state and mlff_workflow.folders(state, wf_base, "train"):
        config_dirs = [d for d in mlff_workflow.folders(state, wf_base, "train") if os.path.isdir(d)]
        source = f"{args.workflow} (training configs)"
    else:
        config_dirs = sorted(d for d in glob.glob("mlff_config_*") if os.path.isdir(d))
        source = "'mlff_config_*'"
    round_no = mlff_workflow.next_round_number(state) if state else 1
    work_dir = args.work_dir or (os.path.relpath(os.path.join(wf_base, f"round_{round_no:02d}"))
                                 if state else ".")
    seeds = [args.seed + i for i in range(args.committee)]

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Configs from      : {source}", f_out)
    print_dual(f"Workflow file     : {args.workflow} ({'round ' + str(round_no) if state else 'new'})", f_out)
    print_dual(f"Round directory   : {work_dir}", f_out)
    print_dual(f"Committee         : {args.committee} model(s), seed(s) {', '.join(map(str, seeds))}", f_out)
    print_dual(f"Foundation model  : mace-mp-0 ({args.foundation_model})", f_out)
    print_dual(f"Epochs / batch    : {args.epochs} / {args.batch_size}", f_out)
    print_dual(f"Learning rate     : {args.lr:g}", f_out)
    print_dual(f"E0s               : {args.e0s}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] AGGREGATING TRAINING DATA", f_out)
    if not config_dirs:
        fail(f"No config folders found ({source}).")
    print_dual(f"Config folders    : {len(config_dirs)} matched", f_out)

    atoms_list, n_skipped, stress_ok = collect_training_configs(config_dirs, f_out)
    if len(atoms_list) < 4:
        fail(f"Only {len(atoms_list)} usable configuration(s) found -- need at least a few "
             "for a train/validation split. Check that SIESTA finished in each folder.")

    print_dual(f"[OK] {len(atoms_list)} usable configuration(s) ({n_skipped} skipped).", f_out)
    train_stress = stress_ok and not args.no_stress_training
    if train_stress:
        print_dual(f"[OK] Stress available for every configuration -- trained on, with "
                   f"weight {args.stress_weight:g} (--loss stress).", f_out)
    elif stress_ok:
        print_dual(color_text(
            "[INFO] --no-stress-training: training on energy+forces only (stress still "
            "reported in the validation metrics).", 'yellow'), f_out)
    else:
        print_dual(color_text(
            "[INFO] Stress not available for every configuration -- training on "
            "energy+forces only.", 'yellow'), f_out)

    device = args.device or ("cuda" if _cuda_available() else "cpu")
    print_dual(f"Device            : {device}", f_out)
    os.makedirs(work_dir, exist_ok=True)

    members = []
    for k, seed in enumerate(seeds):
        member_dir = work_dir if args.committee == 1 else os.path.join(work_dir, f"m_s{seed}")
        tag = "" if args.committee == 1 else f" -- MEMBER {k + 1}/{args.committee} (seed {seed})"
        members.append(train_member(args, atoms_list, seed, member_dir, device, train_stress,
                                    tag, compare_foundation=(k == 0 and not args.skip_foundation_comparison),
                                    f_out=f_out, fail=fail))

    if args.committee > 1:
        print_section("[5] COMMITTEE SUMMARY", f_out)
        print_dual(f"{'seed':>6}  {'train/valid':>11}  {'E RMSE':>10}  {'F RMSE':>10}  {'S RMSE':>9}", f_out)
        print_dual(f"{'':>6}  {'':>11}  {'meV/atom':>10}  {'meV/Ang':>10}  {'GPa':>9}", f_out)
        for m in members:
            v = m["valid"]
            s_txt = f"{v['S_rmse_GPa']:9.3f}" if v.get("S_rmse_GPa") is not None else f"{'-':>9}"
            print_dual(f"{m['seed']:>6}  {m['n_train']:>5}/{m['n_valid']:<5}  "
                       f"{v.get('E_rmse_meV', float('nan')):10.2f}  {v.get('F_rmse_meV', float('nan')):10.2f}  {s_txt}", f_out)
        print_dual("Each member is scored on its own validation split; the held-out test set "
                   "(stb-mlffValidate) is the common yardstick.", f_out)

    # ---- workflow state
    if state is None:
        state = mlff_workflow.new_state(base_dir=wf_base)
        mlff_workflow.add_configs(state, [{"folder": d, "set": "train", "kind": "unknown",
                                           "origin": "unknown", "stage": "unknown", "round": 0}
                                          for d in config_dirs], wf_base)
    used = [d for d in config_dirs]
    mlff_workflow.add_round(state, {
        "round": round_no,
        "dir": mlff_workflow.rel(work_dir, wf_base),
        "created": datetime.now().isoformat(timespec="seconds"),
        "configs": [mlff_workflow.rel(d, wf_base) for d in used],
        "n_usable": len(atoms_list),
        "committee": [{"seed": m["seed"], "dir": mlff_workflow.rel(m["dir"], wf_base),
                       "model": mlff_workflow.rel(m["model"], wf_base) if m["model"] else None,
                       "n_train": m["n_train"], "n_valid": m["n_valid"], "valid": m["valid"]}
                      for m in members],
        "settings": {"foundation_model": args.foundation_model, "epochs": args.epochs,
                     "batch_size": args.batch_size, "lr": args.lr, "e0s": args.e0s,
                     "stress_trained": bool(train_stress), "stress_weight": args.stress_weight,
                     "valid_fraction": args.valid_fraction},
        "validation": {}, "active_learning": [],
    })
    mlff_workflow.record(state, "train", "stb-mlffAnalysis", sys.argv[1:])
    mlff_workflow.save(state, args.workflow)

    print_section("[6] SUMMARY & FILES", f_out)
    n_models = sum(1 for m in members if m["model"])
    print_dual("Status            : OK" if n_models == len(members) else
               color_text(f"Status            : {n_models}/{len(members)} model(s) found", 'yellow'), f_out)
    print_dual(f"Round             : {round_no} ({work_dir})", f_out)
    for m in members:
        print_dual(f"Model (seed {m['seed']:>4}) : {m['model'] or '(not found)'}", f_out)
        if m["parity"]:
            print_dual(f"  parity plot     : {m['parity']}", f_out)
        if m["lammps"]:
            print_dual(f"  LAMMPS model    : {m['lammps']}", f_out)
    print_dual(f"Workflow file     : {args.workflow} (round {round_no} recorded)", f_out)
    print_dual("Next              : stb-mlffValidate (4.17.3) to score it on the test set, "
               "stb-mlffActiveLearning (4.17.4) for new configs.", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Custom ML force field ready -- fast, and fit to your own data.\n", 'bold'))


def train_member(args, atoms_list, seed, member_dir, device, train_stress, tag,
                 compare_foundation, f_out, fail):
    """Trains one committee member: its own seeded train/validation split,
    E0s, mace_run_train, model discovery, validation metrics and parity plot.
    Returns {'seed', 'dir', 'model', 'n_train', 'n_valid', 'valid', 'parity',
    'lammps'}."""
    os.makedirs(member_dir, exist_ok=True)
    # Split done HERE (not via mace_run_train's own --valid_fraction) so we
    # know exactly which configurations ended up in validation, and can
    # independently evaluate the fine-tuned model (and the raw foundation
    # model, for comparison) on that same held-out set afterwards -- see
    # split_train_valid()'s docstring.
    train_atoms, valid_atoms = split_train_valid(atoms_list, args.valid_fraction, seed)
    train_path = os.path.join(member_dir, "train_set.xyz")
    valid_path = os.path.join(member_dir, "valid_set.xyz")
    ase_write(train_path, train_atoms, format='extxyz')
    ase_write(valid_path, valid_atoms, format='extxyz')

    print_section(f"[2] FINE-TUNING (mace_run_train){tag}", f_out)
    print_dual(f"[OK] Wrote {len(train_atoms)} training / {len(valid_atoms)} validation "
               f"configuration(s) ({train_path} / {valid_path})", f_out)

    e0s_arg = "average"
    if args.e0s == "foundation-residual":
        from stb.core import mace_relax
        from ase.data import chemical_symbols
        found_calc = mace_relax.get_calculator(model=args.foundation_model, device=device)
        e0s, rank_ok, start_rmse, avg_e0s = foundation_residual_e0s(train_atoms, found_calc)
        if e0s is None:
            print_dual(color_text(
                "[WARNING] An element is outside the foundation model -- falling back to "
                "--e0s average.", 'yellow'), f_out)
        else:
            e0s_arg = "{" + ", ".join(f"{z}: {v!r}" for z, v in e0s.items()) + "}"
            for z, v in e0s.items():
                print_dual(f"E0({chemical_symbols[z]:<2})            : {v:.6f} eV "
                           f"(\"average\" would give {avg_e0s[z]:.6f} eV)", f_out)
            print_dual(f"Starting E RMSE   : {start_rmse * 1000:.2f} meV/atom (foundation model "
                       "with these E0s, before any training)", f_out)
            if not rank_ok:
                print_dual(color_text(
                    "[WARNING] The training compositions don't separate the elements' E0s "
                    "(e.g. every config has the same stoichiometry): only that composition's "
                    "total is determined, so energies of OTHER compositions are not reliable.",
                    'yellow'), f_out)

    cmd = [
        "mace_run_train",
        "--name", args.name,
        "--work_dir", member_dir,
        "--train_file", train_path,
        "--valid_file", valid_path,
        "--energy_key", "REF_energy",
        "--forces_key", "REF_forces",
        "--foundation_model", args.foundation_model,
        # Plain fine-tuning (--multiheads_finetuning explicitly False --
        # mace_run_train defaults it to True), with E0s on the SIESTA scale:
        # SIESTA's absolute total energy is on a different scale than the
        # foundation model's (different code, pseudopotentials/core
        # treatment), so "--E0s foundation" sits ~1e5 meV/atom off and never
        # converges -- and --multiheads_finetuning (the anti-forgetting mode)
        # structurally REQUIRES "--E0s foundation". Trade-off accepted: no
        # pretraining replay, targets on the right scale. The E0s are
        # foundation_residual_e0s() by default (explicit {Z: E0} dict), or
        # mace_run_train's "average" on request -- see that function for why
        # "average" is a poor start.
        "--E0s", e0s_arg,
        "--multiheads_finetuning", "False",
        "--max_num_epochs", str(args.epochs),
        "--batch_size", str(args.batch_size),
        "--lr", str(args.lr),
        "--seed", str(seed),
        "--device", device,
    ]
    if train_stress:
        # --compute_stress alone only makes the model OUTPUT a stress; it is
        # the stress loss that actually fits it (the default "weighted" loss
        # is energy+forces only, whatever keys are passed).
        cmd.extend(["--stress_key", "REF_stress", "--compute_stress", "True",
                    "--loss", "stress", "--stress_weight", str(args.stress_weight)])

    print_dual(f"Command           : {' '.join(cmd)}", f_out)
    log_path = os.path.join(member_dir, "mace_run_train.log")
    try:
        with open(log_path, "w") as log_f:
            result = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT)
    except FileNotFoundError:
        fail("'mace_run_train' could not be executed (not found on PATH at run time).")

    if result.returncode != 0:
        fail(f"mace_run_train exited with code {result.returncode} -- see {log_path} "
             "for the full training log.")
    print_dual(f"[OK] mace_run_train finished -- full log at {log_path}", f_out)

    print_section(f"[3] MODEL DISCOVERY{tag}", f_out)
    model_path = find_mace_model(member_dir, args.name)
    lammps_model_path = None
    if model_path:
        print_dual(f"[OK] Fine-tuned model: {model_path}", f_out)
        smoke_test_model(model_path, atoms_list[0], device, f_out)
        if args.export_lammps:
            lammps_model_path = export_lammps_model(model_path, f_out)
    else:
        print_dual(color_text(
            "[WARNING] Training finished but no '<name>*.model' file was found under "
            f"'{member_dir}' -- check the log for the model's actual output path.",
            'yellow'), f_out)
        if args.export_lammps:
            print_dual(color_text(
                "[WARNING] --export-lammps skipped: no model file to convert.", 'yellow'), f_out)

    parity_plot_path = None
    valid = {}
    if model_path:
        print_section(f"[4] VALIDATION & PARITY{tag}", f_out)
        from mace.calculators import MACECalculator
        ft_calc = MACECalculator(model_paths=[model_path], device=device)
        pred_e, ref_e, pred_f, ref_f, pred_s, ref_s = evaluate_on_configs(ft_calc, valid_atoms)
        series = [("Fine-tuned", pred_e, ref_e, pred_f, ref_f)]
        if len(ref_e):
            report_model_metrics("Fine-tuned", pred_e, ref_e, pred_f, ref_f, pred_s, ref_s,
                                 len(ref_e), f_out)
            rmse_e, _, rmse_e_rel = energy_metrics(pred_e, ref_e)
            valid = {"E_rmse_meV": 1000 * rmse_e, "E_rmse_rel_meV": 1000 * rmse_e_rel,
                     "F_rmse_meV": 1000 * _rmse(pred_f, ref_f),
                     "S_rmse_GPa": _rmse(pred_s, ref_s) * _EV_A3_TO_GPA if len(ref_s) else None,
                     "n": int(len(ref_e))}
        else:
            print_dual(color_text(
                "[WARNING] Fine-tuned model could not be evaluated on the validation set.",
                'yellow'), f_out)

        if compare_foundation:
            from stb.core import mace_relax
            found_calc = mace_relax.get_calculator(model=args.foundation_model, device=device)
            fpred_e, fref_e, fpred_f, fref_f, fpred_s, fref_s = evaluate_on_configs(found_calc, valid_atoms)
            series.append((f"Foundation ({args.foundation_model})", fpred_e, fref_e, fpred_f, fref_f))
            if len(fref_e):
                # No absolute E RMSE here: it is the SIESTA/MACE scale offset,
                # not an error (see energy_metrics).
                report_model_metrics(f"Foundation ({args.foundation_model}, no fine-tuning)",
                                     fpred_e, fref_e, fpred_f, fref_f, fpred_s, fref_s,
                                     len(fref_e), f_out, show_absolute=False)

        if len(ref_e) or len(ref_f):
            parity_plot_path = os.path.join(member_dir, f"{args.name}_parity.png")
            plot_parity(series, parity_plot_path)
            print_dual(f"[OK] Saved parity plot: {parity_plot_path}", f_out)
            if args.save_data:
                for label, pe, re_, pf, rf in series:
                    tag_ = label.split()[0].lower()
                    if len(re_):
                        e_path = os.path.join(member_dir, f"{args.name}_parity_{tag_}_energy.dat")
                        np.savetxt(e_path, np.column_stack([re_, pe]),
                                   header="ref_E pred_E (eV/atom, absolute -- subtract each "
                                          "column's mean to compare models)")
                        print_dual(f"[OK] Saved {e_path}", f_out)
                    if len(rf):
                        f_path = os.path.join(member_dir, f"{args.name}_parity_{tag_}_forces.dat")
                        np.savetxt(f_path, np.column_stack([rf, pf]), header="ref_F pred_F (eV/Ang)")
                        print_dual(f"[OK] Saved {f_path}", f_out)

    return {"seed": seed, "dir": member_dir, "model": model_path, "n_train": len(train_atoms),
            "n_valid": len(valid_atoms), "valid": valid, "parity": parity_plot_path,
            "lammps": lammps_model_path, "log": log_path}

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
