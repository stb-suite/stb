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
import re
import json
import time
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import write as ase_write
from ase.calculators.calculator import CalculationFailed

from stb.core import structure_io, siesta_log
from stb.core.deps import require_sisl
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                          run_with_spinner, capture_library_noise)
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


def _scf_converged(out_file):
    """True/False when the SIESTA log says whether the last SCF converged, None if it doesn't say."""
    try:
        with open(out_file, errors="replace") as f:
            text = f.read()
    except OSError:
        return None
    if "SCF_NOT_CONV" in text or "did not converge" in text:
        return False
    if "SCF cycle converged" in text or "converged after" in text:
        return True
    return None


def collect_training_configs(config_dirs, f_out=None, issues=None):
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

    Returns (atoms_list, n_skipped, stress_available). Each Atoms carries
    info['folder'] and info['scf_converged']. With `issues` (a list), the
    per-folder problems are appended to it instead of printed.
    """
    sisl = require_sisl()

    def _issue(msg):
        if issues is not None:
            issues.append(msg)
        else:
            print_dual(color_text(f"[WARNING] {msg}", 'yellow'), f_out)
    atoms_list = []
    n_skipped = 0
    stress_available = True

    for d in sorted(config_dirs):
        label = structure_io.find_fdf_system_label(d)
        if not label:
            _issue(f"{d}: no SystemLabel found -- skipping.")
            n_skipped += 1
            continue

        out_file = siesta_log.find_out_file(d, label)
        fa_file = os.path.join(d, f"{label}.FA")
        xv_file = os.path.join(d, f"{label}.XV")

        if not (out_file and os.path.isfile(fa_file) and os.path.isfile(xv_file)):
            _issue(f"{d}: missing .out/.FA/.XV (calculation not finished?) -- skipping.")
            n_skipped += 1
            continue

        energy = siesta_log.get_free_energy(out_file)
        forces = siesta_log.read_fa_forces(fa_file)
        if energy is None or forces is None:
            _issue(f"{d}: could not parse energy/forces -- skipping.")
            n_skipped += 1
            continue

        try:
            geom = sisl.get_sile(xv_file).read_geometry()
        except Exception as e:
            _issue(f"{d}: could not read '{xv_file}': {e} -- skipping.")
            n_skipped += 1
            continue

        if len(geom.xyz) != len(forces):
            _issue(f"{d}: atom count mismatch between .XV ({len(geom.xyz)}) and "
                f".FA ({len(forces)}) -- skipping.")
            n_skipped += 1
            continue

        atoms = Atoms(symbols=[a.symbol for a in geom.atoms], positions=geom.xyz,
                      cell=np.array(geom.cell), pbc=True)
        atoms.info['REF_energy'] = energy
        atoms.info['folder'] = d
        atoms.info['scf_converged'] = _scf_converged(out_file)
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


# --------------------------------------------------------------------------
# Evaluation of the fitted models on the dataset

def predict_all(calc, atoms_list):
    """One single point per config with `calc`: list of
    {'e': eV/atom, 'f': (n,3) eV/Ang, 's': 3x3 eV/Ang^3 or None}, or None
    where the model failed on that config."""
    out = []
    for a in atoms_list:
        b = a.copy()
        b.calc = calc
        try:
            e = b.get_potential_energy() / len(b)
            f = b.get_forces()
            st = np.asarray(b.get_stress(voigt=False)) if 'REF_stress' in a.info else None
        except Exception:
            out.append(None)
            continue
        out.append({"e": float(e), "f": np.asarray(f), "s": st})
    return out


def metrics(preds, atoms_list, idx):
    """Error metrics of `preds` on the configs `idx` (None when empty):
    E RMSE absolute / offset-removed / offset (meV/atom), F RMSE / MAE / max
    component error (meV/Ang), stress RMSE (GPa), n."""
    idx = [i for i in idx if preds[i] is not None]
    if not idx:
        return None
    pe = np.array([preds[i]["e"] for i in idx])
    re_ = np.array([atoms_list[i].info["REF_energy"] / len(atoms_list[i]) for i in idx])
    df = np.concatenate([(preds[i]["f"] - atoms_list[i].arrays["REF_forces"]).ravel() for i in idx])
    ds = [(preds[i]["s"] - np.asarray(atoms_list[i].info["REF_stress"])).ravel()
          for i in idx if preds[i]["s"] is not None and "REF_stress" in atoms_list[i].info]
    d = pe - re_
    off = float(np.mean(d))
    return {"n": len(idx), "E_rmse_meV": 1000 * float(np.sqrt(np.mean(d ** 2))),
            "E_rmse_rel_meV": 1000 * float(np.sqrt(np.mean((d - off) ** 2))) if len(idx) > 1 else None,
            "E_offset_meV": 1000 * off,
            "F_rmse_meV": 1000 * float(np.sqrt(np.mean(df ** 2))),
            "F_mae_meV": 1000 * float(np.mean(np.abs(df))),
            "F_maxerr_meV": 1000 * float(np.max(np.abs(df))),
            "S_rmse_GPa": float(np.sqrt(np.mean(np.concatenate(ds) ** 2))) * _EV_A3_TO_GPA if ds else None}


def _fmt(v, f="{:.1f}"):
    return "-" if v is None else f.format(v)


def read_mace_history(member_dir, seed):
    """Validation metrics per epoch from mace_run_train's results/*_train.txt
    (JSON lines), mean training loss per epoch, and the epoch whose
    checkpoint was kept. Returns {'epochs', 'val_f', 'val_e', 'val_s',
    'train_loss', 'initial', 'best_epoch'} (lists in meV units) or None."""
    files = glob.glob(os.path.join(member_dir, "results", f"*run-{seed}_train.txt"))
    if not files:
        return None
    ev, opt, initial = {}, {}, None
    with open(files[0]) as f:
        for line in f:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get("mode") == "eval":
                if r.get("epoch") is None:
                    initial = r
                else:
                    ev[int(r["epoch"])] = r
            elif r.get("mode") == "opt" and r.get("epoch") is not None:
                opt.setdefault(int(r["epoch"]), []).append(r["loss"])
    best = None
    log = os.path.join(member_dir, "mace_run_train.log")
    if os.path.isfile(log):
        with open(log, errors="replace") as f:
            m = re.findall(r"Loaded Stage one model from epoch (\d+)", f.read())
        if m:
            best = int(m[-1])
    epochs = sorted(ev)
    g = lambda r, k: 1000 * r[k] if r.get(k) is not None else None
    return {"epochs": epochs, "val_f": [g(ev[e], "rmse_f") for e in epochs],
            "val_e": [g(ev[e], "rmse_e_per_atom") for e in epochs],
            "val_s": [g(ev[e], "rmse_stress") for e in epochs],
            "train_loss": [float(np.mean(opt[e])) if e in opt else None for e in epochs],
            "initial": {"F": g(initial, "rmse_f"), "E": g(initial, "rmse_e_per_atom")} if initial else None,
            "best_epoch": best}


def dataset_table(atoms_list, kinds, f_out):
    """Prints the per-kind summary of the training data and returns the
    per-config quantities used later (E/atom, max|F|, pressure, volume)."""
    e_atom = np.array([a.info["REF_energy"] / len(a) for a in atoms_list])
    fmax = np.array([np.linalg.norm(a.arrays["REF_forces"], axis=1).max() for a in atoms_list])
    frms = np.array([np.sqrt(np.mean(np.sum(a.arrays["REF_forces"] ** 2, axis=1))) for a in atoms_list])
    press = np.array([-np.trace(np.asarray(a.info["REF_stress"])) / 3 * _EV_A3_TO_GPA
                      if "REF_stress" in a.info else np.nan for a in atoms_list])
    vol = np.array([a.get_volume() / len(a) for a in atoms_list])
    e0 = e_atom.min()
    rows = []
    for kind in sorted(set(kinds)):
        ii = [i for i, k in enumerate(kinds) if k == kind]
        rng = lambda x, f: f"{f.format(np.nanmin(x[ii]))} .. {f.format(np.nanmax(x[ii]))}" if np.isfinite(x[ii]).any() else "-"
        rows.append(([kind, str(len(ii)), rng(1000 * (e_atom - e0), "{:.1f}"), rng(frms, "{:.2f}"),
                      rng(fmax, "{:.2f}"), rng(press, "{:+.2f}"), rng(vol, "{:.3f}")], None))
    print_table(["Kind", "N", "E - E_min (meV/atom)", "F_rms (eV/Ang)", "F_max (eV/Ang)",
                 "P (GPa)", "V/atom (Ang^3)"], rows, f_out)
    return {"e_atom": e_atom, "fmax": fmax, "press": press, "vol": vol}


def plot_round_summary(out_path, histories, seeds, atoms_list, oos, found, kinds, show):
    """One figure for the round: validation learning curves (F, E) of every
    member, out-of-sample force/energy parity and error distribution against
    the raw foundation model, and the force RMSE per configuration kind."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    ax = axes[0, 0]
    for seed, h in zip(seeds, histories):
        if h and h["epochs"]:
            line, = ax.plot(h["epochs"], h["val_f"], "o-", ms=3, label=f"seed {seed}")
            if h["best_epoch"] in h["epochs"]:
                k = h["epochs"].index(h["best_epoch"])
                ax.plot([h["best_epoch"]], [h["val_f"][k]], "*", ms=14, color=line.get_color())
    ax.set_yscale("log")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Validation F RMSE (meV/Ang)")
    ax.set_title("Learning curves (star = kept checkpoint)"); ax.legend(fontsize=8)
    ax = axes[0, 1]
    for seed, h in zip(seeds, histories):
        if h and h["epochs"]:
            ax.plot(h["epochs"], h["val_e"], "o-", ms=3, label=f"seed {seed}")
    ax.set_yscale("log")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Validation E RMSE (meV/atom)")
    ax.set_title("Energy (absolute, incl. offset)"); ax.legend(fontsize=8)

    idx = [i for i in range(len(atoms_list)) if oos[i] is not None]
    ref_f = np.concatenate([atoms_list[i].arrays["REF_forces"].ravel() for i in idx]) if idx else np.array([])
    oos_f = np.concatenate([oos[i]["f"].ravel() for i in idx]) if idx else np.array([])
    fidx = [i for i in range(len(atoms_list)) if found and found[i] is not None]
    ax = axes[0, 2]
    if fidx:
        ax.scatter(np.concatenate([atoms_list[i].arrays["REF_forces"].ravel() for i in fidx]),
                   np.concatenate([found[i]["f"].ravel() for i in fidx]), s=5, alpha=0.35,
                   color="tab:orange", label="foundation")
    if idx:
        ax.scatter(ref_f, oos_f, s=6, alpha=0.6, color="tab:blue", label="fine-tuned (out-of-sample)")
    lims = [*ax.get_xlim(), *ax.get_ylim()]
    ax.plot([min(lims), max(lims)], [min(lims), max(lims)], "--", color="gray", lw=0.8)
    ax.set_xlabel("SIESTA force (eV/Ang)"); ax.set_ylabel("Predicted force (eV/Ang)")
    ax.set_title("Force parity"); ax.legend(fontsize=8)

    ax = axes[1, 0]
    re_ = lambda ii: np.array([atoms_list[i].info["REF_energy"] / len(atoms_list[i]) for i in ii])
    # Both on the same configs (the held-out ones when there are 2+), each
    # centred on its own mean there, so the two clouds are comparable.
    eidx = idx if len(idx) > 1 else fidx
    efidx = [i for i in eidx if i in set(fidx)]
    if len(efidx) > 1:
        r, pr = re_(efidx), np.array([found[i]["e"] for i in efidx])
        ax.scatter(1000 * (r - r.mean()), 1000 * (pr - pr.mean()), s=20, alpha=0.6, color="tab:orange", label="foundation")
    if len(idx) > 1:
        r, pr = re_(idx), np.array([oos[i]["e"] for i in idx])
        ax.scatter(1000 * (r - r.mean()), 1000 * (pr - pr.mean()), s=24, alpha=0.8, color="tab:blue",
                   label="fine-tuned (out-of-sample)")
    lims = [*ax.get_xlim(), *ax.get_ylim()]
    ax.plot([min(lims), max(lims)], [min(lims), max(lims)], "--", color="gray", lw=0.8)
    ax.set_xlabel("SIESTA E - mean (meV/atom)"); ax.set_ylabel("Predicted E - mean (meV/atom)")
    ax.set_title(f"Energy parity, offset removed ({len(eidx)} configs)")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)

    ax = axes[1, 1]
    errs = []
    if fidx:
        errs.append(("foundation", 1000 * np.abs(np.concatenate(
            [(found[i]["f"] - atoms_list[i].arrays["REF_forces"]).ravel() for i in fidx])), "tab:orange", 0.5))
    if idx:
        errs.append(("fine-tuned (out-of-sample)", 1000 * np.abs(oos_f - ref_f), "tab:blue", 0.7))
    if errs:
        top = max(e.max() for _, e, _, _ in errs)
        bins = np.logspace(np.log10(max(1e-1, min(max(e.min(), 1e-1) for _, e, _, _ in errs))), np.log10(top * 1.05), 35)
        for lab, e, col, al in errs:
            ax.hist(np.maximum(e, bins[0]), bins=bins, alpha=al, color=col, label=lab)
        ax.set_xscale("log")
    ax.set_xlabel("|F error| per component (meV/Ang, log)"); ax.set_ylabel("Count")
    ax.set_title("Force error distribution"); ax.legend(fontsize=8)

    ax = axes[1, 2]
    names = sorted(set(kinds))
    x = np.arange(len(names))
    for off, (lab, preds, col) in enumerate((("foundation", found, "tab:orange"),
                                             ("fine-tuned (out-of-sample)", oos, "tab:blue"))):
        if not preds:
            continue
        vals = []
        for kname in names:
            m = metrics(preds, atoms_list, [i for i, k in enumerate(kinds) if k == kname])
            vals.append(m["F_rmse_meV"] if m else np.nan)
        ax.bar(x + (off - 0.5) * 0.38, vals, 0.38, color=col, label=lab)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=20, fontsize=8)
    ax.set_ylabel("F RMSE (meV/Ang)"); ax.set_title("Force error per configuration kind"); ax.legend(fontsize=8)
    fig.suptitle("stb-mlffAnalysis -- training round summary", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    if show:
        plt.show()
    plt.close(fig)


def plot_parity(series, out_path):
    """Per-member 2-panel energy (offset removed) / force parity on its own
    validation split, fine-tuned vs the raw foundation model."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for label, pred_e, ref_e, pred_f, ref_f in series:
        if len(ref_e) > 1:
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
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    if len(series[0][2]) < 2:
        axes[0].text(0.5, 0.5, "fewer than 2 validation configs", ha="center", transform=axes[0].transAxes)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def split_indices(n, valid_fraction, seed):
    """Same seeded split as split_train_valid, as index lists (train, valid)."""
    n_valid = min(max(1, round(n * valid_fraction)), n - 1)
    order = np.random.default_rng(seed).permutation(n)
    valid = sorted(order[:n_valid].tolist())
    return [i for i in range(n) if i not in set(valid)], valid


def cv_folds(n, k, seed):
    """K-fold cross-validation: (train, valid) index lists for each fold, the
    configs shuffled once with `seed` and dealt into k near-equal folds, so
    every config is held out exactly once."""
    order = np.random.default_rng(seed).permutation(n)
    folds = [sorted(f.tolist()) for f in np.array_split(order, k)]
    return [([i for i in range(n) if i not in set(f)], f) for f in folds]


# ---- level-of-theory consistency ------------------------------------------

_THEORY_KEYS = [("XC.functional", ("xcfunctional",)), ("XC.authors", ("xcauthors",)),
                ("PAO.BasisSize", ("paobasissize",)), ("PAO.EnergyShift", ("paoenergyshift",)),
                ("Mesh.Cutoff", ("meshcutoff",)), ("ElectronicTemperature", ("electronictemperature",)),
                ("Spin", ("spin", "spinpolarized")), ("DFTD3", ("dftd3",)),
                ("SCF tolerance", ("scfdmtolerance", "dmtolerance"))]


def theory_signature(folder, atoms):
    """What SIESTA was told in `folder`, in its own reading order
    (config_extra.fdf first -- it is %included on top --, then the main .fdf,
    then the rest): the settings that change the energy surface, the k-point
    sampling as k x a per axis (Ang, so a supercell with a divided grid still
    matches its parent cell), a hash of any %block PAO.Basis and of each
    pseudopotential file. Values are normalized strings; absent = 'default'."""
    import hashlib
    from stb.core.fdf_directives import parse_directives, parse_kgrid
    files = sorted(glob.glob(os.path.join(folder, "*.fdf")))
    extra = [f for f in files if os.path.basename(f) == "config_extra.fdf"]
    main = [f for f in files if f not in extra and structure_io.get_system_label(f)]
    rest = [f for f in files if f not in extra and f not in main]
    texts = []
    for f in extra + main + rest:
        try:
            with open(f, errors="replace") as fh:
                texts.append(fh.read())
        except OSError:
            pass
    d = parse_directives(*texts)
    sig = {}
    for label, keys in _THEORY_KEYS:
        v = next((d[k] for k in keys if k in d), "default")
        v = " ".join(v.split())
        sig[label] = re.sub(r"(\d)[dD]([+-]?\d)", r"\1e\2", v)    # Fortran 1.0d-5 == 1.0e-5
    kg = parse_kgrid(*texts)
    if kg:
        lengths = atoms.cell.lengths()
        sig["k x a (Ang)"] = " ".join(f"{k * l:.0f}" if k > 1 else "G" for k, l in zip(kg, lengths))
        sig["_kgrid"] = (tuple(int(k) for k in kg), tuple(float(l) for l in lengths))
    else:
        sig["k x a (Ang)"] = f"cutoff {d['kgridcutoff']}" if "kgridcutoff" in d else "Gamma"
    blob = "\n".join(texts)
    m = re.search(r"%block\s+PAO\.Basis\b(.*?)%endblock", blob, re.I | re.S)
    sig["PAO.Basis block"] = hashlib.md5(m.group(1).encode()).hexdigest()[:8] if m else "none"
    for pp in sorted(glob.glob(os.path.join(folder, "*.psf")) + glob.glob(os.path.join(folder, "*.psml"))):
        with open(pp, "rb") as fh:
            sig[f"pseudo {os.path.basename(pp)}"] = hashlib.md5(fh.read()).hexdigest()[:8]
    return sig


def _kxa_in_reference_cell(sigs):
    """Rewrites each signature's 'k x a (Ang)' as the k-grid per reference cell
    times that cell's length: per axis, the folder's cell is round(l / l_min)
    reference cells (l_min = the shortest length of that axis in the set), so
    k x multiple x l_min. A strained cell with the reference's own grid
    (Stage 1 --volume-strains, Stage 4 strain sampling: the same grid on
    purpose) then matches, and so does a supercell with a divided grid; only a
    different k-point count per reference cell shows up. Pops '_kgrid'."""
    grids = [s_.pop("_kgrid", None) for s_ in sigs]
    known = [g for g in grids if g]
    if not known:
        return
    l_min = [min(g[1][i] for g in known) for i in range(3)]
    for s_, g in zip(sigs, grids):
        if not g:
            continue
        ks, ls = g
        mult = [max(1, int(round(l / lm))) if lm > 0 else 1 for l, lm in zip(ls, l_min)]
        s_["k x a (Ang)"] = " ".join(f"{k * m * lm:.0f}" if k * m > 1 else "G"
                                     for k, m, lm in zip(ks, mult, l_min))


def dispersion_check(theory_common, theory_odd, state, f_out, warn):
    """Prints whether the labels include DFT-D3 (the folders' own DFTD3, as
    SIESTA reads it) and checks it against the workflow file's record of
    Stage 1 (--d3). Returns True/False (D3 in most folders or not)."""
    on = mlff_workflow.fdf_true(theory_common.get("DFTD3", "default"))
    rec = mlff_workflow.dftd3_info(state)
    mixed = "DFTD3" in theory_odd
    print_dual("Dispersion        : " + ("DFT-D3 in the labels -- the model learns PBE+D3; never add a D3 "
                                         "correction on top of it" if on else "no DFT-D3 in the labels")
               + (f" (workflow: {rec['source']})" if rec.get("source") else "")
               + (" -- MIXED across folders, see [10]" if mixed else ""), f_out)
    if rec and bool(rec.get("on")) != on:
        warn(f"Dispersion: the workflow file says the labels {'include' if rec.get('on') else 'do not include'} "
             f"DFT-D3 (Stage 1: {rec.get('source')}), but the folders' own settings say "
             f"{'they do' if on else 'they do not'} -- a folder was edited or written outside the workflow; "
             "every label must come from the same energy surface.")
    return on


def check_theory(atoms_list, f_out, warn):
    """Compares every folder's theory_signature with the most common one;
    prints the common settings and warns about every folder that differs."""
    sigs = [theory_signature(a.info["folder"], a) for a in atoms_list]
    _kxa_in_reference_cell(sigs)
    keys = sorted({k for s_ in sigs for k in s_}, key=lambda k: (k.startswith("pseudo"), k))
    common, odd = {}, {}
    for k in keys:
        vals = [s_.get(k, "missing") for s_ in sigs]
        low = [v.lower() for v in vals]
        top = max(set(low), key=low.count)
        common[k] = vals[low.index(top)]
        for a, v, lv in zip(atoms_list, vals, low):
            if lv != top:
                odd.setdefault(k, []).append((os.path.basename(a.info["folder"]), v))
    main_keys = [k for k in keys if not k.startswith("pseudo")]
    shown = [f"{k} {common[k]}" for k in main_keys if common[k] != "default"]
    print_dual("Level of theory   : " + ", ".join(shown[:4]), f_out)
    for k0 in range(4, len(shown), 4):
        print_dual("                    " + ", ".join(shown[k0:k0 + 4]), f_out)
    unset = [k for k in main_keys if common[k] == "default"]
    if unset:
        print_dual(f"                    (SIESTA default for: {', '.join(unset)})", f_out)
    pp = [f"{k[7:]} {common[k]}" for k in keys if k.startswith("pseudo")]
    if pp:
        print_dual(f"Pseudopotentials  : {', '.join(pp)} (md5, first 8)", f_out)
    if odd:
        for k, lst in odd.items():
            warn(f"Data: level of theory differs in {len(lst)} folder(s) for {k} (most use '{common[k]}'): "
                 + ", ".join(f"{f} '{v}'" for f, v in lst[:6]) + (" ..." if len(lst) > 6 else "")
                 + " -- energies from different settings are not on one surface; rerun or remove them.")
    else:
        print_dual(color_text(f"[OK] Same level of theory in all {len(atoms_list)} folders "
                              "(settings, k-sampling, basis, pseudopotential files).", 'green'), f_out)
    return common, odd


# ---- model checks: MD stability and speed ---------------------------------

def stability_cell(atoms, target=9.0, max_atoms=256):
    """Supercell for the stability MD: periodic axes repeated until each is at
    least `target` Ang long (vacuum axes left alone), capped at `max_atoms`."""
    from stb.core.kspace import detect_vacuum_axes
    vac = detect_vacuum_axes(atoms.get_scaled_positions(), np.asarray(atoms.cell), 5.0)
    reps = [1 if v else max(1, int(np.ceil(target / l))) for v, l in zip(vac, atoms.cell.lengths())]
    while len(atoms) * int(np.prod(reps)) > max_atoms and max(reps) > 1:
        reps[int(np.argmax(reps))] -= 1
    return atoms.repeat(reps), reps


def _dmin(atoms, cutoff):
    from ase.neighborlist import neighbor_list
    d = neighbor_list("d", atoms, cutoff)
    return float(d.min()) if len(d) else cutoff


def stability_md(model_path, atoms0, temperature, steps, device, seed):
    """A short MD with the model: NVT (Langevin, 1 fs) for `steps`, then NVE
    for steps//2 to measure the total-energy drift. Stops early if the run
    blows up (NaN, atoms closer than 60 % of the starting shortest distance,
    T > 5x the target, |F| > 50 eV/Ang). Returns a dict of diagnostics,
    including the wall time per step (the speed benchmark)."""
    from mace.calculators import MACECalculator
    from stb.mlmd import build_dynamics
    from ase import units
    import logging
    atoms = atoms0.copy()
    level = logging.getLogger().level
    logging.getLogger().setLevel(logging.ERROR)     # "converting models to float32" notice
    try:
        atoms.calc = MACECalculator(model_paths=[model_path], device=device, default_dtype="float32")
    finally:
        logging.getLogger().setLevel(level)
    d0 = _dmin(atoms, 6.0)
    out = {"n_atoms": len(atoms), "status": "OK", "reason": "", "T_mean": None, "dmin_ratio": None,
           "fmax": None, "drift_meV_atom_ps": None, "ms_per_step": None}
    atoms.get_forces()                               # warm-up (model load, graph build)
    dyn = build_dynamics(atoms, "nvt", timestep_fs=1.0, temperature_k=temperature, friction_fs=0.01,
                         pressure_gpa=0.0, seed=seed)
    temps, dmins, fmaxs = [], [], []
    t0 = time.monotonic()
    n_done = 0

    def _check():
        f = atoms.get_forces()
        if not (np.all(np.isfinite(f)) and np.all(np.isfinite(atoms.positions))):
            return "NaN in forces/positions"
        fm = float(np.linalg.norm(f, axis=1).max())
        dm = _dmin(atoms, 6.0)
        tk = atoms.get_temperature()
        temps.append(tk); dmins.append(dm); fmaxs.append(fm)
        if dm < 0.6 * d0:
            return f"atoms at {dm:.2f} Ang (start {d0:.2f})"
        if tk > 5 * max(temperature, 50):
            return f"T = {tk:.0f} K"
        if fm > 50:
            return f"|F| = {fm:.0f} eV/Ang"
        return None
    for step in range(1, steps + 1):
        try:
            dyn.run(1)
        except Exception as e:
            out["status"], out["reason"] = "UNSTABLE", f"MD step failed: {e}"
            break
        n_done += 1
        if step % 10 == 0:
            bad = _check()
            if bad:
                out["status"], out["reason"] = "UNSTABLE", f"step {step}: {bad}"
                break
    nvt_time = time.monotonic() - t0
    out["ms_per_step"] = 1000 * nvt_time / max(1, n_done)
    if out["status"] == "OK":
        nve = build_dynamics(atoms, "nve", timestep_fs=1.0, temperature_k=temperature, friction_fs=0.01,
                             pressure_gpa=0.0, seed=seed + 1)
        etot, tt = [], []
        for step in range(1, steps // 2 + 1):
            try:
                nve.run(1)
            except Exception as e:
                out["status"], out["reason"] = "UNSTABLE", f"NVE step failed: {e}"
                break
            if step % 5 == 0:
                bad = _check()
                if bad:
                    out["status"], out["reason"] = "UNSTABLE", f"NVE step {step}: {bad}"
                    break
                etot.append(atoms.get_potential_energy() + atoms.get_kinetic_energy())
                tt.append(step * 1e-3)   # ps
        if len(etot) > 3:
            slope = np.polyfit(tt, np.array(etot) / len(atoms), 1)[0]
            out["drift_meV_atom_ps"] = 1000 * float(slope)
    if temps:
        out["T_mean"] = float(np.mean(temps))
        out["dmin_ratio"] = float(min(dmins) / d0)
        out["fmax"] = float(max(fmaxs))
    return out


def _versions():
    out = {}
    for mod in ("mace", "torch", "ase"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = "?"
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Stage 2 (Train) of the Custom ML Force Field workflow: aggregates "
                     "the finished SIESTA calculations (energy + forces, best-effort stress), "
                     "builds a MACE-format training set, and fine-tunes a MACE-MP foundation "
                     "model on it via 'mace_run_train' -- by default a committee of 3 models "
                     "(different seeds), recorded as a new round in mlff_workflow.json for "
                     "Validate (stb-mlffValidate) and Active learning (stb-mlffActiveLearning). "
                     "Reports the dataset, each member's training history and errors (train, "
                     "validation, out-of-sample per kind, vs the raw foundation model), the "
                     "committee agreement, and every warning in one final section.",
        epilog="Example usage:\n"
               "  stb-mlffAnalysis                      # configs from mlff_workflow.json, 3-model committee\n"
               "  stb-mlffAnalysis --show               # also display the summary figure\n"
               "  stb-mlffAnalysis --path 'mlff_config_*' --epochs 100 --device cpu --committee 1\n"
               "  stb-mlffAnalysis --foundation-model small --device cuda\n"
               "  stb-mlffAnalysis --export-lammps",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--path", default=None,
                        help="Glob pattern matching the config folders to train on (default: the "
                             "training configs listed in --workflow, else every mlff_config_* here or in "
                             f"{mlff_workflow.RUN_DIR}/, in configs_round_*/ subfolders or directly).")
    parser.add_argument("--workflow", default=None, metavar="JSON",
                        help="Workflow state file shared by the 4 stages (default: "
                             f"{mlff_workflow.RUN_DIR}/{mlff_workflow.WORKFLOW_FILE}, or "
                             f"./{mlff_workflow.WORKFLOW_FILE} when run inside that folder). Read for "
                             "the config list, updated with this round (committee models and "
                             "metrics); created next to the config folders if missing.")
    parser.add_argument("--skip-unlabelled", action="store_true",
                        help="Train even when the configs meant for this round (<run>/configs_round_NN/, "
                             "e.g. Stage 4's picks) have no finished SIESTA run yet -- they are left out and "
                             "train a later round. Without it the stage stops and says where to run SIESTA.")
    parser.add_argument("--committee", type=int, default=DEFAULT_COMMITTEE,
                        help=f"Number of models to train (default: {DEFAULT_COMMITTEE}), seeds --seed, "
                             "--seed+1, ...: each gets its own train/validation split and "
                             "training seed. Their force disagreement is what stb-mlffActiveLearning "
                             "ranks candidates by, and stb-mlffValidate scores the committee mean. "
                             "1 trains a single model (that score then falls back to the largest "
                             "predicted force).")
    parser.add_argument("--cv", type=int, default=0, metavar="K",
                        help="K-fold cross-validation instead of a random --valid-fraction per member: "
                             "the configs are dealt into K folds and member k validates on fold k, so "
                             "EVERY config gets one out-of-sample prediction (the committee is then the "
                             "K fold models; overrides --committee). Best with few configs. Default: off.")
    parser.add_argument("--foundation-model", default="medium", choices=["small", "medium", "large"],
                        help="MACE-MP foundation checkpoint (default: medium; small is ~2x faster). It "
                             "is the starting point unless --init-from says otherwise, and always the "
                             "baseline of the comparison.")
    parser.add_argument("--init-from", default="foundation", metavar="foundation|previous|PATH",
                        help="Starting weights: 'foundation' (default, the MACE-MP checkpoint), "
                             "'previous' (member k continues from member k of the latest round in "
                             "--workflow: faster active-learning rounds) or a .model file.")
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
                             "data for comparison. On by default (comparison runs unless this is passed).")
    parser.add_argument("--no-stability", dest="stability", action="store_false",
                        help="Skip the MD stability test (on by default: a short NVT then NVE run with "
                             "each model, in a supercell of the reference, checking it does not blow up "
                             "and how much the NVE energy drifts; it also times the model).")
    parser.add_argument("--stability-temperature", type=float, default=300.0, metavar="K",
                        help="Temperature of the stability MD (default: 300 K).")
    parser.add_argument("--stability-steps", type=int, default=400,
                        help="NVT steps of the stability MD, 1 fs each, followed by half as many NVE "
                             "steps (default: 400).")
    parser.add_argument("--view-worst", type=int, default=0, metavar="N",
                        help="Write the N worst-predicted configs to <round>/worst_configs.xyz (per-atom "
                             "force error as 'F_error' and as the initial charge, so ASE's viewer can "
                             "color by it: View > Colors > By initial charge) and open them. Default: 0.")
    parser.add_argument("--show", action="store_true",
                        help="Also display the round summary figure on screen (it is always saved).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write <round>/mlff_training_errors.dat: per config, its kind, SIESTA "
                             "max force, out-of-sample and foundation F RMSE and committee spread. Off by default.")
    parser.add_argument("--export-lammps", action="store_true",
                        help="Also compile the fine-tuned model into a TorchScript file "
                             "loadable by LAMMPS's 'pair_style mace' (via mace_create_lammps_model, "
                             "part of mace-torch) -- '<model>-lammps.pt'. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE} (in the round directory). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlffAnalysis {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Train)",
            "Fine-tune MACE-MP on real SIESTA data",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CUSTOM ML FORCE FIELD -- TRAIN:", 'bold'))
    print("-" * 60)

    if shutil.which("mace_run_train") is None:
        print(color_text("\n[CRITICAL ERROR] 'mace_run_train' not found on PATH.", 'red'))
        print("This tool needs the optional 'ml' extra. Install it with:")
        print("  pip install stb_suite[ml]")
        sys.exit(1)

    warns = []           # every warning, printed together in the last section
    library_noise = []   # third-party prints/warnings, summarized there too

    def warn(msg):
        warns.append(msg)

    f_out = None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if warns:
            print_section("WARNINGS SO FAR", f_out)
            for k, w in enumerate(warns, 1):
                print_dual(color_text(f"  {k}. {w}", 'yellow'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    if args.committee < 1:
        fail(f"--committee must be >= 1 (got {args.committee}).")

    # ---- where the configs come from, where this round goes
    explicit_workflow = args.workflow is not None
    args.workflow = mlff_workflow.default_path(args.workflow)
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
        # no workflow file: every mlff_config_* here or in the run folder, per-round
        # subfolders (configs_round_NN/) or the older flat layout
        source = "'[configs_round_*/]mlff_config_*'"
        config_dirs = mlff_workflow.find_config_dirs(".")
        if not config_dirs:
            source = f"'{mlff_workflow.RUN_DIR}/[configs_round_*/]mlff_config_*'"
            config_dirs = mlff_workflow.find_config_dirs(mlff_workflow.RUN_DIR)
    if state is None and not explicit_workflow and config_dirs:
        # No workflow file yet: start one in the run folder (above configs_round_NN/).
        parent = os.path.commonpath([os.path.abspath(os.path.dirname(os.path.abspath(d))) for d in config_dirs])
        if os.path.basename(parent).startswith(mlff_workflow.CONFIGS_PREFIX):
            parent = os.path.dirname(parent)
        args.workflow = os.path.relpath(os.path.join(parent, mlff_workflow.WORKFLOW_FILE))
        wf_base = mlff_workflow.base_dir(args.workflow)
    round_no = mlff_workflow.next_round_number(state) if state else 1
    work_dir = args.work_dir or (os.path.relpath(os.path.join(wf_base, f"round_{round_no:02d}"))
                                 if state else ".")
    if args.cv:
        if args.cv < 2:
            fail(f"--cv must be >= 2 (got {args.cv}).")
        args.committee = args.cv
    seeds = [args.seed + i for i in range(args.committee)]
    # Starting weights per member
    init_models = [None] * len(seeds)
    if args.init_from == "previous":
        prev = mlff_workflow.get_round(state) if state else None
        prev_models = []
        for m in (prev or {}).get("committee", []):
            if m.get("model"):
                path = mlff_workflow.resolve(m["model"], wf_base)
                plain = path.replace("_compiled.model", ".model")
                prev_models.append(plain if os.path.isfile(plain) else path)
        prev_models = [p_ for p_ in prev_models if os.path.isfile(p_)]
        if not prev_models:
            fail("--init-from previous: no model of an earlier round found in the workflow file.")
        init_models = [prev_models[k % len(prev_models)] for k in range(len(seeds))]
    elif args.init_from != "foundation":
        if not os.path.isfile(args.init_from):
            fail(f"--init-from: '{args.init_from}' is neither foundation/previous nor a model file.")
        init_models = [args.init_from] * len(seeds)
    os.makedirs(work_dir, exist_ok=True)
    report_path = os.path.join(work_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    device = args.device or ("cuda" if _cuda_available() else "cpu")
    t_start = time.monotonic()

    print_dual(color_text("===== STB-MLFFANALYSIS REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    vers = _versions()
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Software          : stb {VERSION}, mace-torch {vers['mace']}, torch {vers['torch']}, ase {vers['ase']}", f_out)
    dev_txt = device
    if device == "cuda":
        try:
            from stb.core.mace_relax import gpu_available
            dev_txt += f" ({gpu_available()[1]})"
        except Exception:
            pass
    print_dual(f"Device            : {dev_txt}", f_out)
    print_dual(f"Configs from      : {source}", f_out)
    print_dual(f"Workflow file     : {args.workflow} ({'round ' + str(round_no) if state else 'new'})", f_out)
    print_dual(f"Round directory   : {work_dir}", f_out)
    print_dual(f"Committee         : {args.committee} model(s), seed(s) {', '.join(map(str, seeds))}", f_out)
    print_dual(f"Foundation model  : mace-mp-0 ({args.foundation_model})", f_out)
    print_dual("Starting weights  : " + ("foundation model" if init_models[0] is None else
                                       ", ".join(sorted(set(init_models)))), f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    # ---- [1] data
    print_section("[1] TRAINING DATA", f_out)
    if not config_dirs:
        fail(f"No config folders found ({source}).")
    issues = []
    atoms_list, n_skipped, stress_ok = collect_training_configs(config_dirs, f_out, issues=issues)
    # unfinished folders are reported once per folder they live in (configs_round_NN/), with the loop
    pending = mlff_workflow.pending_by_dir(config_dirs)
    for msg in issues:
        if "missing .out/.FA/.XV" not in msg:
            warn(f"Data: {msg}")
    for parent, n, k in pending:
        if k:
            todo = [d for d in config_dirs if os.path.dirname(os.path.normpath(d)) == parent
                    and not mlff_workflow.siesta_finished(d)]
            warn(f"SIESTA: {k} of {n} folder(s) in {os.path.relpath(parent)}/ have no finished run -- run SIESTA "
                 f"there ({', '.join(os.path.basename(d) for d in todo[:5])}{' ...' if len(todo) > 5 else ''}); "
                 "this round trains without them.")
    own = os.path.normpath(mlff_workflow.configs_dir(wf_base, round_no)) if state else None
    own_pending = [(p, n, k) for p, n, k in pending if own and os.path.normpath(os.path.abspath(p)) ==
                   os.path.normpath(os.path.abspath(own)) and k == n]
    if own_pending and not args.skip_unlabelled:
        p, n, _ = own_pending[0]
        fail(f"None of the {n} configs meant for round {round_no} ({os.path.relpath(p)}/) has a finished SIESTA run -- "
             f"run SIESTA there first: {mlff_workflow.siesta_loop(p, [d for d in config_dirs if os.path.dirname(os.path.normpath(d)) == p])}"
             f"  (or pass --skip-unlabelled to train round {round_no} without them; they would then train a later "
             "round).")
    print_dual(f"Folders           : {len(config_dirs)} listed, {len(atoms_list)} usable, {n_skipped} skipped"
               + (" (see WARNINGS)" if n_skipped else ""), f_out)
    if len(pending) > 1 or any(k for _, _, k in pending):
        print_dual("By folder         : " + ", ".join(
            f"{os.path.basename(p) or p} {n - k}/{n} labelled" for p, n, k in pending), f_out)
    if len(atoms_list) < 4:
        todo = [f"{os.path.relpath(p)}/ ({k} folder(s))" for p, n, k in pending if k]
        fail(f"Only {len(atoms_list)} usable configuration(s) found -- need at least a few "
             "for a train/validation split. " + (f"Run SIESTA in {', '.join(todo)} first." if todo else
                                                  "Check that SIESTA finished in each folder."))
    kinds = []
    for a in atoms_list:
        entry = mlff_workflow.config_entry(state, a.info["folder"], wf_base) if state else None
        kinds.append((entry or {}).get("kind", "unknown") or "unknown")
    formulas = sorted({a.get_chemical_formula(mode="hill") for a in atoms_list})
    sizes = sorted({len(a) for a in atoms_list})
    species = sorted({s_ for a in atoms_list for s_ in a.get_chemical_symbols()})
    print_dual(f"Composition       : {', '.join(formulas)} ({', '.join(map(str, sizes))} atoms/config; "
               f"species {', '.join(species)})", f_out)
    n_atoms_tot = sum(len(a) for a in atoms_list)
    print_dual(f"Labels            : {len(atoms_list)} energies, {3 * n_atoms_tot} force components"
               + (f", {6 * len(atoms_list)} stress components" if stress_ok else ", stress missing in some configs"), f_out)
    ds = dataset_table(atoms_list, kinds, f_out)
    print_dual(f"Energy spread     : {1000 * (ds['e_atom'].max() - ds['e_atom'].min()):.1f} meV/atom "
               f"(E/atom {ds['e_atom'].min():.4f} .. {ds['e_atom'].max():.4f} eV)", f_out)
    print_dual(f"Force magnitude   : mean {np.mean([np.linalg.norm(a.arrays['REF_forces'], axis=1).mean() for a in atoms_list]):.3f}, "
               f"max {ds['fmax'].max():.3f} eV/Ang", f_out)
    if np.isfinite(ds["press"]).any():
        print_dual(f"Pressure range    : {np.nanmin(ds['press']):+.2f} .. {np.nanmax(ds['press']):+.2f} GPa "
                   "(from SIESTA stress)", f_out)
    print_dual(f"Volume/atom range : {ds['vol'].min():.3f} .. {ds['vol'].max():.3f} Ang^3", f_out)
    theory_common, theory_odd = check_theory(atoms_list, f_out, warn)
    labels_d3 = dispersion_check(theory_common, theory_odd, state, f_out, warn)
    not_conv = [os.path.basename(a.info["folder"]) for a in atoms_list if a.info.get("scf_converged") is False]
    if not_conv:
        warn(f"Data: SCF did not converge in {len(not_conv)} config(s) ({', '.join(not_conv)}): their energies "
             "and forces are not converged labels -- rerun or remove them.")
    else:
        print_dual(color_text("[OK] SCF converged in every config the log reports on.", 'green'), f_out)
    big = [os.path.basename(a.info["folder"]) for a, fm in zip(atoms_list, ds["fmax"])
           if fm > mlff_workflow.UNPHYSICAL_FMAX]
    if big:
        warn(f"Data: {len(big)} config(s) have SIESTA forces above {mlff_workflow.UNPHYSICAL_FMAX:g} eV/Ang "
             f"({', '.join(big)}): very far from anything MD visits; they can dominate the force loss -- "
             "remove them from the workflow file's training set, and use a smaller --stdev "
             "(Stage 1) or --max-force (Stage 4) next time.")
    if ds["vol"].max() - ds["vol"].min() < 1e-6 * ds["vol"].max():
        warn("Data: every config has the same volume -- the model has no information on E(V), the bulk "
             "modulus or pressure (stb-mlff --volume-strains).")
    if len(atoms_list) < 20 and not args.cv:
        warn(f"Data: only {len(atoms_list)} configs -- each member validates on "
             f"{max(1, round(len(atoms_list) * args.valid_fraction))} of them, so its validation numbers are "
             "noisy; --cv 5 gives every config an out-of-sample prediction.")

    train_stress = stress_ok and not args.no_stress_training
    # ---- [2] recipe
    print_section("[2] TRAINING SETUP", f_out)
    print_dual(f"Epochs / batch    : {args.epochs} / {args.batch_size}", f_out)
    print_dual(f"Learning rate     : {args.lr:g}", f_out)
    print_dual(f"Loss              : " + (f"energy 1 + forces 100 + stress {args.stress_weight:g} (--loss stress)"
                                          if train_stress else "energy 1 + forces 100 (no stress)"), f_out)
    print_dual(f"E0s               : {args.e0s}", f_out)
    if args.cv:
        if args.cv > len(atoms_list):
            fail(f"--cv {args.cv} needs at least {args.cv} configs (have {len(atoms_list)}).")
        folds = cv_folds(len(atoms_list), args.cv, args.seed)
        print_dual(f"Validation        : {args.cv}-fold cross-validation (fold sizes "
                   f"{', '.join(str(len(v)) for _, v in folds)}; every config held out once)", f_out)
    else:
        folds = [split_indices(len(atoms_list), args.valid_fraction, sd) for sd in seeds]
        print_dual(f"Validation split  : {args.valid_fraction:g} per member "
                   f"({max(1, round(len(atoms_list) * args.valid_fraction))} config(s), different for each seed)", f_out)
    if init_models[0] is not None:
        print_dual("Starting weights  : " + "; ".join(f"seed {sd} <- {m_}" for sd, m_ in zip(seeds, init_models)), f_out)
    if stress_ok and args.no_stress_training:
        warn("Setup: --no-stress-training although every config has a stress: the volume response is not fitted.")
    elif not stress_ok:
        warn("Setup: stress missing in some configs -- trained on energy + forces only; the volume "
             "response (E(V), pressure) is not fitted.")

    # ---- [3] training
    print_section("[3] TRAINING", f_out)
    members = []
    for k, seed in enumerate(seeds):
        member_dir = work_dir if args.committee == 1 else os.path.join(work_dir, f"m_s{seed}")
        members.append(train_member(args, atoms_list, seed, member_dir, device, train_stress,
                                    f"{k + 1}/{args.committee}", warn, library_noise, f_out, fail,
                                    split=folds[k], init_model=init_models[k]))
    rows = []
    for m in members:
        h = m["history"] or {}
        ini = (h.get("initial") or {})
        fin = h["val_f"][h["epochs"].index(h["best_epoch"])] if h.get("best_epoch") in h.get("epochs", []) else None
        rows.append(([str(m["seed"]), f"{m['n_train']}/{m['n_valid']}", _fmt(m["start_rmse"], "{:.1f}"),
                      _fmt(ini.get("F"), "{:.0f}"), _fmt(fin, "{:.0f}"),
                      f"{h.get('best_epoch', '-')}/{len(h.get('epochs', [])) - 1 if h.get('epochs') else '-'}",
                      f"{m['seconds']:.0f}", "ok" if m["model"] else "NO MODEL"], None if m["model"] else 'red'))
    print_table(["Seed", "Train/valid", "Start E (meV/at)", "Val F start", "Val F kept", "Kept epoch",
                 "Time (s)", "Model"], rows, f_out)
    print_dual("Start E = foundation model with the fitted E0s on the training split, before any epoch; "
               "Val F = mace_run_train's validation force RMSE (meV/Ang) before training and at the kept checkpoint.", f_out)
    for m in members:
        h = m["history"]
        if h and h.get("best_epoch") is not None and h["epochs"] and h["best_epoch"] >= h["epochs"][-1] - 1 \
                and args.epochs > 5:
            warn(f"Training: seed {m['seed']} kept epoch {h['best_epoch']} of 0-{h['epochs'][-1]} (at the end): "
                 "the validation error was still falling -- more --epochs would likely help.")
    usable = [m for m in members if m["model"]]
    if not usable:
        fail("No member produced a model file -- see the mace_run_train logs.")

    # ---- [4] errors
    print_section("[4] ERRORS ON THE DATA (independent evaluation, every model on every config)", f_out)
    from mace.calculators import MACECalculator
    all_idx = list(range(len(atoms_list)))
    for m in usable:
        with capture_library_noise(library_noise, f"load {m['model']}"):
            calc = MACECalculator(model_paths=[m["model"]], device=device)
        m["preds"] = run_with_spinner(predict_all, calc, atoms_list, label=f"Evaluating seed {m['seed']}")
        m["train_m"] = metrics(m["preds"], atoms_list, m["train_idx"])
        m["valid_m"] = metrics(m["preds"], atoms_list, m["valid_idx"])
    found = None
    if not args.skip_foundation_comparison:
        from stb.core import mace_relax
        with capture_library_noise(library_noise, "foundation model"):
            fcalc = mace_relax.get_calculator(model=args.foundation_model, device=device)
        found = run_with_spinner(predict_all, fcalc, atoms_list, label="Evaluating the foundation model")
    found_all = metrics(found, atoms_list, all_idx) if found else None

    rows = []
    for m in usable:
        for split, mm in (("train", m["train_m"]), ("valid", m["valid_m"])):
            if mm:
                rows.append(([f"seed {m['seed']}", split, str(mm["n"]), _fmt(mm["E_rmse_meV"]),
                              _fmt(mm["E_rmse_rel_meV"]), _fmt(mm["F_rmse_meV"]), _fmt(mm["F_mae_meV"]),
                              _fmt(mm["F_maxerr_meV"], "{:.0f}"), _fmt(mm["S_rmse_GPa"], "{:.3f}")], None))
    if found_all:
        rows.append(([f"foundation ({args.foundation_model})", "all", str(found_all["n"]), "(scale offset)",
                      _fmt(found_all["E_rmse_rel_meV"]), _fmt(found_all["F_rmse_meV"]), _fmt(found_all["F_mae_meV"]),
                      _fmt(found_all["F_maxerr_meV"], "{:.0f}"), _fmt(found_all["S_rmse_GPa"], "{:.3f}")], 'cyan'))
    print_table(["Model", "Set", "N", "E abs", "E rel", "F RMSE", "F MAE", "F max", "S RMSE"], rows, f_out)
    print_dual("E in meV/atom (abs = absolute; rel = mean offset removed), F in meV/Ang (max = largest "
               "component error), S in GPa. The foundation's absolute E is on another energy scale.", f_out)

    # Out-of-sample pool: each config predicted by the members that did NOT
    # train on it. Energies first lose each member's own mean offset on its
    # training split, so pooling members with different offsets does not
    # show up as a (spurious) energy-difference error.
    for m in usable:
        tr = [i for i in m["train_idx"] if m["preds"][i] is not None]
        m["e_offset"] = float(np.mean([m["preds"][i]["e"] - atoms_list[i].info["REF_energy"] / len(atoms_list[i])
                                        for i in tr])) if tr else 0.0
    oos = []
    for i in all_idx:
        ps = [(m["preds"][i], m["e_offset"]) for m in usable if i in m["valid_idx"] and m["preds"][i] is not None]
        if ps:
            ps_raw = [p for p, _ in ps]
            oos.append({"e": float(np.mean([p["e"] - off for p, off in ps])),
                        "f": np.mean([p["f"] for p in ps_raw], axis=0),
                        "s": (np.mean([p["s"] for p in ps_raw], axis=0) if all(p["s"] is not None for p in ps_raw) else None)})
        else:
            oos.append(None)
    oos_all = metrics(oos, atoms_list, all_idx)
    n_oos = sum(o is not None for o in oos)
    if oos_all:
        print_dual(f"\nOut-of-sample pool: {n_oos}/{len(atoms_list)} configs were held out by at least one "
                   "member -- each scored with the members that did not train on it:", f_out)
        rows = []
        for kname in sorted(set(kinds)):
            ii = [i for i, k in enumerate(kinds) if k == kname]
            mo = metrics(oos, atoms_list, ii)
            mf = metrics(found, atoms_list, [i for i in ii if oos[i] is not None]) if found else None
            gain = (100 * (1 - mo["F_rmse_meV"] / mf["F_rmse_meV"]) if (mo and mf and mf["F_rmse_meV"]) else None)
            rows.append(([kname, str(mo["n"] if mo else 0), _fmt(mo and mo["F_rmse_meV"]), _fmt(mo and mo["E_rmse_meV"]),
                          _fmt(mf and mf["F_rmse_meV"]), _fmt(gain, "{:+.0f} %")], None))
        mf_all = metrics(found, atoms_list, [i for i in all_idx if oos[i] is not None]) if found else None
        gain_all = (100 * (1 - oos_all["F_rmse_meV"] / mf_all["F_rmse_meV"]) if mf_all and mf_all["F_rmse_meV"] else None)
        rows.append((["ALL", str(oos_all["n"]), _fmt(oos_all["F_rmse_meV"]), _fmt(oos_all["E_rmse_meV"]),
                      _fmt(mf_all and mf_all["F_rmse_meV"]), _fmt(gain_all, "{:+.0f} %")], 'bold'))
        print_table(["Kind", "N", "F RMSE fine-tuned", "E RMSE fine-tuned*", "F RMSE foundation", "F gain"], rows, f_out)
        print_dual("* each member's mean energy offset on its own training split removed first (meV/atom).", f_out)
    else:
        gain_all = None
        print_dual("Out-of-sample pool: empty (no member held out any config).", f_out)

    # Interpretation lines
    print_dual("", f_out)
    for m in usable:
        tm, vm = m["train_m"], m["valid_m"]
        if tm and vm and tm["F_rmse_meV"] > 0:
            ratio = vm["F_rmse_meV"] / tm["F_rmse_meV"]
            m["overfit_ratio"] = ratio
            if ratio > 2.0:
                warn(f"Errors: seed {m['seed']} validation F RMSE is {ratio:.1f}x its training F RMSE -- "
                     "overfitting (few configs, or too many epochs); more data helps most.")
    ratios = [m.get("overfit_ratio") for m in usable if m.get("overfit_ratio")]
    if ratios:
        print_dual(f"Valid/train F RMSE: {', '.join(f'{r:.2f}' for r in ratios)} (1 = no overfitting; > 2 flagged)", f_out)
    if gain_all is not None:
        print_dual(f"Fine-tuning gain  : out-of-sample force error {gain_all:+.0f} % vs the raw foundation model", f_out)
        if gain_all < 10:
            warn(f"Errors: fine-tuning improved the out-of-sample force error by only {gain_all:.0f} % over the "
                 "raw foundation model -- check the data (level of theory, SCF) and the recipe.")
    if found_all and oos_all and found_all.get("E_rmse_rel_meV") and oos_all.get("E_rmse_rel_meV"):
        print_dual(f"Energy differences: {oos_all['E_rmse_rel_meV']:.1f} meV/atom out-of-sample vs "
                   f"{found_all['E_rmse_rel_meV']:.1f} for the foundation model (offset removed)", f_out)

    # Worst-predicted configs: out-of-sample error where there is one, else
    # the committee mean on its own training data (marked in-sample).
    def _cfg_err(i):
        if oos[i] is not None:
            return float(np.sqrt(np.mean((oos[i]["f"] - atoms_list[i].arrays["REF_forces"]) ** 2))), "out-of-sample", oos[i]["f"]
        ps = [m["preds"][i]["f"] for m in usable if m["preds"][i] is not None]
        if not ps:
            return None, "-", None
        fm = np.mean(ps, axis=0)
        return float(np.sqrt(np.mean((fm - atoms_list[i].arrays["REF_forces"]) ** 2))), "in-sample", fm
    cerr = [_cfg_err(i) for i in all_idx]
    ranked = sorted([i for i in all_idx if cerr[i][0] is not None], key=lambda i: -cerr[i][0])
    print_dual("\nWorst-predicted configs:", f_out)
    rows = []
    for i in ranked[:5]:
        per_atom = np.linalg.norm(cerr[i][2] - atoms_list[i].arrays["REF_forces"], axis=1)
        j = int(np.argmax(per_atom))
        rows.append(([os.path.basename(atoms_list[i].info["folder"]), kinds[i], f"{1000 * cerr[i][0]:.0f}", cerr[i][1],
                      f"{ds['fmax'][i]:.2f}", f"{atoms_list[i].get_chemical_symbols()[j]}{j + 1} ({1000 * per_atom[j]:.0f})"], None))
    print_table(["Config", "Kind", "F RMSE (meV/Ang)", "Prediction", "SIESTA F_max (eV/Ang)", "Worst atom (|dF|)"], rows, f_out)
    worst_path = None
    if args.view_worst > 0:
        frames = []
        for i in ranked[:args.view_worst]:
            b = atoms_list[i].copy()
            per_atom = np.linalg.norm(cerr[i][2] - atoms_list[i].arrays["REF_forces"], axis=1)
            b.new_array("F_error", per_atom)
            b.set_initial_charges(per_atom)
            b.info = {"folder": os.path.basename(atoms_list[i].info["folder"]), "kind": kinds[i],
                      "F_rmse_meV": round(1000 * cerr[i][0], 1), "prediction": cerr[i][1]}
            frames.append(b)
        worst_path = os.path.join(work_dir, "worst_configs.xyz")
        ase_write(worst_path, frames, format="extxyz")

    # ---- [5] committee
    comm = None
    if len(usable) > 1:
        print_section("[5] COMMITTEE", f_out)
        spreads, mean_sp = [], []
        for i in all_idx:
            fs = [m["preds"][i]["f"] for m in usable if m["preds"][i] is not None]
            if len(fs) < 2:
                spreads.append(np.nan); mean_sp.append(np.nan)
                continue
            fs = np.array(fs)
            sp = np.sqrt(np.mean(np.sum((fs - fs.mean(0)) ** 2, axis=2), axis=0))
            spreads.append(float(sp.max())); mean_sp.append(float(sp.mean()))
        spreads, mean_sp = np.array(spreads), np.array(mean_sp)
        es = np.array([np.std([m["preds"][i]["e"] for m in usable if m["preds"][i] is not None]) for i in all_idx])
        print_dual(f"Force spread      : mean {1000 * np.nanmean(mean_sp):.1f}, max {1000 * np.nanmax(spreads):.1f} meV/Ang "
                   "(per-atom RMS deviation from the committee mean)", f_out)
        print_dual(f"Energy spread     : mean {1000 * np.nanmean(es):.2f}, max {1000 * np.nanmax(es):.2f} meV/atom", f_out)
        corr = None
        err = np.array([np.sqrt(np.mean((oos[i]["f"] - atoms_list[i].arrays["REF_forces"]) ** 2)) if oos[i] is not None
                        else np.nan for i in all_idx])
        ok = np.isfinite(err) & np.isfinite(spreads)
        if ok.sum() >= 5 and np.std(spreads[ok]) > 0 and np.std(err[ok]) > 0:
            corr = float(np.corrcoef(spreads[ok], err[ok])[0, 1])
            print_dual(f"Spread vs error   : correlation {corr:+.2f} over {ok.sum()} held-out configs "
                       f"(error / spread ~ {np.nanmedian(err[ok] / np.maximum(spreads[ok], 1e-9)):.1f}x)", f_out)
            if corr < 0.3:
                warn(f"Committee: its force spread barely tracks the real error here (corr {corr:+.2f}) -- seeds "
                     "fine-tuned from one foundation model are nearly the same model; Stage 4 then ranks "
                     "candidates almost like max|F| does.")
        else:
            print_dual(f"Spread vs error   : only {int(ok.sum())} held-out config(s), too few to correlate "
                       "(stb-mlffValidate does it on the test set).", f_out)
        rows = []
        order = np.argsort(-np.nan_to_num(spreads, nan=-1))[:min(5, len(all_idx))]
        for i in order:
            rows.append(([os.path.basename(atoms_list[i].info["folder"]), kinds[i], f"{1000 * spreads[i]:.0f}",
                          "-" if not np.isfinite(err[i]) else f"{1000 * err[i]:.0f}"], None))
        print_dual("Most disputed configs (where the committee disagrees most):", f_out)
        print_table(["Config", "Kind", "Max spread (meV/Ang)", "Held-out F RMSE"], rows, f_out)
        comm = {"F_spread_mean_meV": 1000 * float(np.nanmean(mean_sp)), "F_spread_max_meV": 1000 * float(np.nanmax(spreads)),
                "spread_error_corr": corr}

    # ---- [6] round comparison
    print_section("[6] COMPARISON WITH EARLIER ROUNDS", f_out)
    prev_rounds = [r for r in (state or {}).get("rounds", []) if r["round"] < round_no]
    member_valid_f = [m["valid_m"]["F_rmse_meV"] for m in usable if m.get("valid_m")]
    cur = {"round": round_no, "n": len(atoms_list),
           "oos": oos_all["F_rmse_meV"] if oos_all else None,
           "valid": float(np.mean(member_valid_f)) if member_valid_f else None,
           "spread": comm["F_spread_mean_meV"] if comm else None, "test": None}
    table = []
    for r in prev_rounds:
        te = r.get("training_errors") or {}
        vals = [c.get("valid", {}).get("F_rmse_meV") for c in r.get("committee", [])]
        vals = [v for v in vals if v is not None]
        table.append({"round": r["round"], "n": r.get("n_usable", len(r.get("configs", []))),
                      "oos": (te.get("out_of_sample") or {}).get("F_rmse_meV"),
                      "valid": float(np.mean(vals)) if vals else None,
                      "spread": (te.get("committee") or {}).get("F_spread_mean_meV"),
                      "test": ((r.get("validation") or {}).get("overall") or {}).get("F_rmse_meV")})
    table.append(cur)
    if len(table) > 1:
        rows = []
        for k, t in enumerate(table):
            prev_t = table[k - 1] if k else None
            dv = (f"{100 * (t['valid'] / prev_t['valid'] - 1):+.0f} %" if prev_t and t["valid"] and prev_t["valid"] else "")
            rows.append(([str(t["round"]) + (" (this)" if t is cur else ""), str(t["n"]), _fmt(t["valid"]), dv,
                          _fmt(t["oos"]), _fmt(t["spread"]), _fmt(t["test"])], 'bold' if t is cur else None))
        print_table(["Round", "Configs", "Valid F RMSE", "vs prev", "Out-of-sample F", "Mean spread", "Test F RMSE"],
                    rows, f_out)
        print_dual("F in meV/Ang. Valid = mean of the members' own validation splits (different configs each "
                   "round); the test-set column (stb-mlffValidate) is the like-for-like comparison.", f_out)
    else:
        print_dual("First round -- nothing to compare yet.", f_out)

    # ---- [7] model checks: MD stability + speed
    stab = []
    speed = None
    if args.stability:
        print_section("[7] MODEL CHECKS: MD STABILITY & SPEED", f_out)
        ref_atoms = None
        ref_file = ((state or {}).get("reference") or {}).get("file")
        if ref_file:
            try:
                from pymatgen.io.ase import AseAtomsAdaptor
                ref_atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(
                    structure_io.read_fdf(mlff_workflow.resolve(ref_file, wf_base))))
            except Exception:
                ref_atoms = None
        src_txt = "reference structure"
        if ref_atoms is None:
            k0 = int(np.argmin(ds["e_atom"]))
            ref_atoms = atoms_list[k0].copy()
            src_txt = f"lowest-energy config ({os.path.basename(atoms_list[k0].info['folder'])})"
        cell, reps = stability_cell(ref_atoms)
        print_dual(f"MD cell           : {src_txt} x {reps[0]}x{reps[1]}x{reps[2]} = {len(cell)} atoms; "
                   f"NVT {args.stability_steps} fs at {args.stability_temperature:g} K, then NVE "
                   f"{args.stability_steps // 2} fs (1 fs steps)", f_out)
        for m in usable:
            plain = m["model"]
            with capture_library_noise(library_noise, f"stability MD seed {m['seed']}"):
                r = run_with_spinner(stability_md, plain, cell, args.stability_temperature, args.stability_steps,
                                     device, m["seed"], label=f"MD stability, seed {m['seed']}")
            r["seed"] = m["seed"]
            stab.append(r)
        rows = []
        for r in stab:
            rows.append(([str(r["seed"]), r["status"], _fmt(r["T_mean"], "{:.0f}"), _fmt(r["dmin_ratio"], "{:.2f}"),
                          _fmt(r["fmax"], "{:.2f}"), _fmt(r["drift_meV_atom_ps"], "{:+.2f}"),
                          _fmt(r["ms_per_step"], "{:.1f}"), r["reason"] or "-"], None if r["status"] == "OK" else 'red'))
        print_table(["Seed", "Status", "T mean (K)", "d_min/start", "F_max (eV/Ang)", "NVE drift (meV/atom/ps)",
                     "ms/step", "Note"], rows, f_out)
        print_dual("d_min/start = shortest distance seen / at the start (< 0.6 = blow-up); a good NVE drift is "
                   "|drift| < ~1 meV/atom/ps at 1 fs.", f_out)
        for r in stab:
            if r["status"] != "OK":
                warn(f"Model check: seed {r['seed']} is UNSTABLE in MD ({r['reason']}) -- it has not learned "
                     "the repulsive part well enough: add configs of this kind (Stage 4 MD sampling) or use "
                     "another member.")
            elif r["drift_meV_atom_ps"] is not None and abs(r["drift_meV_atom_ps"]) > 1.0:
                warn(f"Model check: seed {r['seed']} NVE energy drifts {r['drift_meV_atom_ps']:+.2f} meV/atom/ps at "
                     "1 fs -- the potential surface is rough; use a smaller time step in production MD or retrain.")
        ok_speed = [r["ms_per_step"] for r in stab if r["ms_per_step"]]
        if ok_speed:
            ms = float(np.median(ok_speed))
            speed = {"n_atoms": len(cell), "ms_per_step": ms, "us_per_atom_step": 1000 * ms / len(cell),
                     "ns_per_day_1fs": 86400.0 / (ms / 1000) * 1e-6, "device": device}
            print_dual(f"Speed             : {ms:.1f} ms/step for {len(cell)} atoms = {speed['us_per_atom_step']:.0f} "
                       f"us/atom/step on {device} -> ~{speed['ns_per_day_1fs']:.2f} ns/day at 1 fs "
                       f"(float32; ~{speed['ns_per_day_1fs'] * len(cell) / 1000:.2f} ns/day for 1000 atoms, linear scaling)", f_out)

    # ---- plots
    summary_png = os.path.join(work_dir, "mlff_training.png")
    show = args.show
    if show:
        import matplotlib
        if matplotlib.get_backend().lower() in ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template"):
            warn("Plots: --show needs a display (no interactive matplotlib backend here) -- the figure is only saved.")
            show = False
    try:
        plot_round_summary(summary_png, [m["history"] for m in usable], [m["seed"] for m in usable],
                           atoms_list, oos, found, kinds, show)
    except Exception as e:
        warn(f"Plots: summary figure not written ({e}).")
        summary_png = None
    for m in usable:
        try:
            vi = m["valid_idx"]
            ser = [("Fine-tuned", np.array([m["preds"][i]["e"] for i in vi if m["preds"][i]]),
                    np.array([atoms_list[i].info["REF_energy"] / len(atoms_list[i]) for i in vi if m["preds"][i]]),
                    np.concatenate([m["preds"][i]["f"].ravel() for i in vi if m["preds"][i]]),
                    np.concatenate([atoms_list[i].arrays["REF_forces"].ravel() for i in vi if m["preds"][i]]))]
            if found:
                ser.append((f"Foundation ({args.foundation_model})",
                            np.array([found[i]["e"] for i in vi if found[i]]),
                            np.array([atoms_list[i].info["REF_energy"] / len(atoms_list[i]) for i in vi if found[i]]),
                            np.concatenate([found[i]["f"].ravel() for i in vi if found[i]]),
                            np.concatenate([atoms_list[i].arrays["REF_forces"].ravel() for i in vi if found[i]])))
            m["parity"] = os.path.join(m["dir"], f"{args.name}_parity.png")
            plot_parity(ser, m["parity"])
        except Exception as e:
            warn(f"Plots: parity plot of seed {m['seed']} not written ({e}).")
    if args.save_data:
        dat = os.path.join(work_dir, "mlff_training_errors.dat")
        with open(dat, "w") as f:
            f.write("# Per-config errors, stb-mlffAnalysis\n# folder kind n_atoms E_ref(eV/atom) "
                    "Fmax_ref(eV/Ang) F_rmse_oos(meV/Ang) F_rmse_found(meV/Ang) spread_max(meV/Ang)\n")
            for i in all_idx:
                a = atoms_list[i]
                fo = (np.sqrt(np.mean((oos[i]["f"] - a.arrays["REF_forces"]) ** 2)) * 1000) if oos[i] else np.nan
                ff = (np.sqrt(np.mean((found[i]["f"] - a.arrays["REF_forces"]) ** 2)) * 1000) if (found and found[i]) else np.nan
                sp = (1000 * spreads[i]) if comm else np.nan
                f.write(f"{os.path.basename(a.info['folder'])} {kinds[i]} {len(a)} {a.info['REF_energy'] / len(a):.6f} "
                        f"{ds['fmax'][i]:.4f} {fo:.2f} {ff:.2f} {sp:.2f}\n")

    # ---- workflow state
    if state is None:
        state = mlff_workflow.new_state(base_dir=wf_base)
        mlff_workflow.add_configs(state, [{"folder": d, "set": "train", "kind": "unknown",
                                           "origin": "unknown", "stage": "unknown", "round": 0}
                                          for d in config_dirs], wf_base)
    mlff_workflow.add_round(state, {
        "round": round_no,
        "dir": mlff_workflow.rel(work_dir, wf_base),
        "created": datetime.now().isoformat(timespec="seconds"),
        "configs": [mlff_workflow.rel(d, wf_base) for d in config_dirs],
        "n_usable": len(atoms_list),
        "committee": [{"seed": m["seed"], "dir": mlff_workflow.rel(m["dir"], wf_base),
                       "model": mlff_workflow.rel(m["model"], wf_base) if m["model"] else None,
                       "n_train": m["n_train"], "n_valid": m["n_valid"],
                       "valid": _valid_entry(m.get("valid_m")), "train": _valid_entry(m.get("train_m")),
                       "best_epoch": (m["history"] or {}).get("best_epoch"), "seconds": round(m["seconds"], 1)}
                      for m in members],
        "settings": {"foundation_model": args.foundation_model, "epochs": args.epochs,
                     "batch_size": args.batch_size, "lr": args.lr, "e0s": args.e0s,
                     "stress_trained": bool(train_stress), "stress_weight": args.stress_weight,
                     "valid_fraction": args.valid_fraction, "dftd3": bool(labels_d3)},
        "training_errors": {"out_of_sample": _valid_entry(oos_all), "foundation": _valid_entry(found_all),
                            "gain_pct": gain_all, "committee": comm,
                            "cv_folds": args.cv or None},
        "model_checks": {"stability": [{k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}
                                       for r in stab], "speed": speed},
        "init_from": args.init_from,
        "validation": {}, "active_learning": [],
    })
    mlff_workflow.record(state, "train", "stb-mlffAnalysis", sys.argv[1:])
    mlff_workflow.save(state, args.workflow)

    # ---- [6] files
    print_section("[8] SUMMARY & FILES", f_out)
    n_models = len(usable)
    print_dual(f"Status            : {'OK' if n_models == len(members) else f'{n_models}/{len(members)} models'}"
               + (f" ({len(warns)} warning(s), see [10])" if warns else ""), f_out)
    print_dual(f"Round             : {round_no} ({work_dir}), {time.monotonic() - t_start:.0f} s in total", f_out)
    for m in members:
        print_dual(f"Model seed {m['seed']:<6} : {m['model'] or '(not found)'}", f_out)
        extra = [x for x in (m.get("parity"), m.get("lammps"), m["log"]) if x]
        print_dual(f"                    {', '.join(extra)}", f_out)
    if summary_png:
        print_dual(f"Summary figure    : {summary_png}", f_out)
    if args.save_data:
        print_dual(f"Per-config errors : {os.path.join(work_dir, 'mlff_training_errors.dat')}", f_out)
    if worst_path:
        print_dual(f"Worst configs     : {worst_path} (per-atom F_error; color by initial charge in ASE)", f_out)
    print_dual(f"Workflow file     : {args.workflow} (round {round_no} recorded with every number above)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    print_section("[9] NEXT STEPS", f_out)
    print_dual("1. Score the round on the held-out test set: stb-mlffValidate (4.17.3).", f_out)
    print_dual("2. Look for configs the committee is unsure of: stb-mlffActiveLearning (4.17.4).", f_out)
    print_dual(f"3. Use a model now: --custom-model {usable[0]['model']} in stb-mlmd, stb-mlphonons, "
               "stb-mlelastic, ...", f_out)
    if labels_d3:
        print_dual("   It is PBE+D3 already (the labels include DFT-D3): use it without any extra D3 "
                   "correction (MACE's dispersion=True / torch-dftd would count it twice).", f_out)
    extras = [t for t, used in (("--show (display the summary figure)", args.show),
                                ("--save-data (per-config errors)", args.save_data),
                                ("--view-worst N (worst configs in ASE's viewer)", args.view_worst),
                                ("--cv K (every config out-of-sample)", args.cv),
                                ("--init-from previous (continue from the last round)", args.init_from != "foundation"))
              if not used]
    if extras:
        print_dual(f"Also available    : {'; '.join(extras)}.", f_out)

    print_section("[10] WARNINGS", f_out)
    if warns:
        for k, w in enumerate(warns, 1):
            print_dual(color_text(f"{k:>2}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)
    if library_noise:
        print_dual(f"Library messages  : {len(library_noise)} informational message(s) from mace/torch captured "
                   "(model loading, dtype); the full training output of each member is in its "
                   "mace_run_train.log.", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Custom ML force field ready -- fast, and fit to your own data.\n", 'bold'))
    if worst_path:
        from stb.core.ase_view import view_structure_interactive
        from ase.io import read as ase_read
        print(color_text(f"[INFO] Opening the {args.view_worst} worst config(s) in ASE's viewer "
                         "(View > Colors > By initial charge = per-atom force error, eV/Ang).", 'cyan'))
        view_structure_interactive(ase_read(worst_path, ":"))


def _valid_entry(m):
    """Metrics dict as stored in the workflow file (rounded); None stays None."""
    if not m:
        return {}
    return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()}


def train_member(args, atoms_list, seed, member_dir, device, train_stress, tag, warn,
                 library_noise, f_out, fail, split=None, init_model=None):
    """Trains one committee member: its own seeded train/validation split,
    E0s, mace_run_train (with an elapsed-time spinner), model discovery and
    smoke test. Prints one compact block; problems go to `warn`. Returns a
    dict with the split indices, model path, history and timing."""
    os.makedirs(member_dir, exist_ok=True)
    # Split done HERE (not via mace_run_train's own --valid_fraction) so we
    # know exactly which configurations ended up in validation, and can
    # independently evaluate the fine-tuned model (and the raw foundation
    # model, for comparison) on that same held-out set afterwards.
    train_idx, valid_idx = split or split_indices(len(atoms_list), args.valid_fraction, seed)
    start_from = init_model or args.foundation_model
    train_atoms = [atoms_list[i] for i in train_idx]
    valid_atoms = [atoms_list[i] for i in valid_idx]
    train_path = os.path.join(member_dir, "train_set.xyz")
    valid_path = os.path.join(member_dir, "valid_set.xyz")
    ase_write(train_path, train_atoms, format='extxyz')
    ase_write(valid_path, valid_atoms, format='extxyz')
    print_dual(color_text(f"Member {tag} -- seed {seed}", 'bold') +
               f": {len(train_atoms)} train / {len(valid_atoms)} valid "
               f"(valid: {', '.join(os.path.basename(a.info['folder']) for a in valid_atoms)})", f_out)

    e0s_arg = "average"
    start_rmse = None
    if args.e0s == "foundation-residual":
        from stb.core import mace_relax
        from ase.data import chemical_symbols
        with capture_library_noise(library_noise, "starting model (E0s)"):
            found_calc = mace_relax.get_calculator(model=start_from, device=device)
        e0s, rank_ok, start_rmse, avg_e0s = foundation_residual_e0s(train_atoms, found_calc)
        if e0s is None:
            warn(f"Training: seed {seed}: an element is outside the foundation model -- E0s fell back to 'average'.")
        else:
            e0s_arg = "{" + ", ".join(f"{z}: {v!r}" for z, v in e0s.items()) + "}"
            print_dual("  E0s             : " + ", ".join(
                f"{chemical_symbols[z]} {v:.4f} eV ('average' {avg_e0s[z]:.4f})" for z, v in e0s.items())
                + f"; start E RMSE {start_rmse * 1000:.1f} meV/atom", f_out)
            if not rank_ok:
                warn("Training: the training compositions don't separate the elements' E0s (every config "
                     "has the same stoichiometry): only that composition's total is determined, so "
                     "energies of OTHER compositions are not reliable.")

    cmd = [
        "mace_run_train",
        "--name", args.name,
        "--work_dir", member_dir,
        "--train_file", train_path,
        "--valid_file", valid_path,
        "--energy_key", "REF_energy",
        "--forces_key", "REF_forces",
        "--foundation_model", start_from,
        # Plain fine-tuning (--multiheads_finetuning explicitly False --
        # mace_run_train defaults it to True), with E0s on the SIESTA scale:
        # SIESTA's absolute total energy is on a different scale than the
        # foundation model's, so "--E0s foundation" sits ~1e5 meV/atom off,
        # and --multiheads_finetuning structurally REQUIRES it. See
        # foundation_residual_e0s() for why "average" is a poor start.
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
        # the stress loss that actually fits it.
        cmd.extend(["--stress_key", "REF_stress", "--compute_stress", "True",
                    "--loss", "stress", "--stress_weight", str(args.stress_weight)])

    log_path = os.path.join(member_dir, "mace_run_train.log")
    with open(os.path.join(member_dir, "mace_run_train.cmd"), "w") as f:
        f.write(" ".join(cmd) + "\n")
    t0 = time.monotonic()
    try:
        with open(log_path, "w") as log_f:
            result = run_with_spinner(subprocess.run, cmd, stdout=log_f, stderr=subprocess.STDOUT,
                                      label=f"Training seed {seed} ({args.epochs} epochs max)")
    except FileNotFoundError:
        fail("'mace_run_train' could not be executed (not found on PATH at run time).")
    seconds = time.monotonic() - t0
    if result.returncode != 0:
        fail(f"mace_run_train exited with code {result.returncode} (seed {seed}) -- see {log_path}.")

    model_path = find_mace_model(member_dir, args.name)
    history = read_mace_history(member_dir, seed)
    lammps_model_path = None
    if model_path:
        if not smoke_test_model(model_path, atoms_list[0], device, library_noise):
            warn(f"Training: seed {seed}: the model file was written but could not be loaded/evaluated.")
        if args.export_lammps:
            lammps_model_path = export_lammps_model(model_path, warn)
    else:
        warn(f"Training: seed {seed} finished but no '{args.name}*.model' was found under '{member_dir}' "
             f"-- see {log_path}.")
    h = history or {}
    kept = h.get("best_epoch")
    print_dual(f"  done in {seconds:.0f} s, {len(h.get('epochs', []))} epochs logged, kept epoch {kept if kept is not None else '?'}"
               f" -> {model_path or 'NO MODEL'}", f_out)
    return {"seed": seed, "dir": member_dir, "model": model_path, "n_train": len(train_atoms),
            "n_valid": len(valid_atoms), "train_idx": train_idx, "valid_idx": valid_idx,
            "history": history, "start_rmse": (1000 * start_rmse if start_rmse is not None else None),
            "seconds": seconds, "lammps": lammps_model_path, "log": log_path, "parity": None}


def _cuda_available():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


def smoke_test_model(model_path, sample_atoms, device, library_noise=None):
    """Loads the fine-tuned model and evaluates one configuration: a minimal
    "does the model actually load and run" check. Returns True/False."""
    try:
        from mace.calculators import MACECalculator
        noise = library_noise if library_noise is not None else []
        with capture_library_noise(noise, f"load {model_path}"):
            calc = MACECalculator(model_paths=[model_path], device=device)
        atoms = sample_atoms.copy()
        atoms.calc = calc
        return np.isfinite(atoms.get_potential_energy())
    except (ImportError, CalculationFailed, Exception):
        return False


def export_lammps_model(model_path, warn):
    """Compiles `model_path` into a TorchScript file LAMMPS's 'pair_style
    mace' can load directly (via mace_create_lammps_model, part of
    mace-torch -- writes '<model_path>-lammps.pt'). Requires a LAMMPS build
    with the MACE pair style to USE the file. Returns the output path, or
    None on failure (never aborts the run)."""
    if shutil.which("mace_create_lammps_model") is None:
        warn("Export: --export-lammps skipped: 'mace_create_lammps_model' not found on PATH.")
        return None
    out_path = model_path + "-lammps.pt"
    cmd = ["mace_create_lammps_model", model_path, "--dtype", "float64", "--format", "libtorch"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.isfile(out_path):
        warn("Export: could not export a LAMMPS model: "
             f"{result.stderr.strip().splitlines()[-1] if result.stderr else 'unknown error'}")
        return None
    return out_path


if __name__ == "__main__":
    main()
