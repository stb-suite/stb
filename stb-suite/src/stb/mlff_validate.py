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
import time
import shutil
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from stb.core import mlff_workflow, model_card
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                          print_progress_line, finish_progress_line, capture_library_noise,
                          run_with_spinner)
from stb.mlff_analysis import (collect_training_configs, _EV_A3_TO_GPA, predict_all, metrics,
                               dataset_table, check_theory, dispersion_check, _versions, _fmt, refresh_card_validation, rebuild_round_cards)
from stb.mlff import write_config_folder

REPORT_FILE = "stb_mlffValidate_report.txt"
DATA_FILE = "mlff_validation_errors.dat"
WORST_FILE = "mlff_test_worst.xyz"
DEFAULT_PLATEAU_PCT = 5.0
# Above this ratio of true force error to committee spread the spread is
# reported as uncalibrated: still usable to RANK candidates, not as an error bar.
_CALIBRATION_WARN_RATIO = 3.0
# Test error this many times the Stage-2 out-of-sample error: the test set
# probes something the training/validation configs do not.
_GENERALIZATION_WARN_RATIO = 2.0
# A member this many times worse than the committee's median member.
_OUTLIER_MEMBER_RATIO = 1.5
# Correlations need a minimum sample to mean anything.
_MIN_CONFIGS_CORR = 3
_MIN_ATOMS_CORR = 10
_IMAGINARY_THZ = -0.05
PHONON_SYMPREC = 0.01   # as stb-phononsCreate / stb-phononsPos
# Below this force spread (meV/Ang, e.g. a config at equilibrium) error/spread is noise/noise.
_MIN_SPREAD_MEV = 1.0
EOS_FORM = "birchmurnaghan"
MD_CHECK_DIR = "md_check"
# --compare: a config counts as improved/worse beyond this change of its F RMSE.
_COMPARE_BAND_PCT = 10.0
# MD check: error at the end of the trajectory this many times the start.
_MD_GROWTH_WARN = 2.0
_N_BOOTSTRAP = 500


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _f_rmse(pred, atoms):
    return 1000 * float(np.sqrt(np.mean((pred["f"] - atoms.arrays["REF_forces"]) ** 2)))


def bootstrap_f_rmse(preds, atoms_list, n=_N_BOOTSTRAP, seed=0):
    """Standard deviation of the test F RMSE (meV/Ang) when the test configs
    are resampled with replacement: how much of a difference between two
    rounds a test set this small can resolve. None with fewer than 2 configs."""
    idx = [i for i, p in enumerate(preds) if p is not None]
    if len(idx) < 2:
        return None
    ss = np.array([np.sum((preds[i]["f"] - atoms_list[i].arrays["REF_forces"]) ** 2) for i in idx])
    nc = np.array([atoms_list[i].arrays["REF_forces"].size for i in idx])
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(idx), size=(n, len(idx)))
    return float(1000 * np.std(np.sqrt(ss[draws].sum(axis=1) / nc[draws].sum(axis=1))))


def score_committee(calcs, test_atoms, kinds):
    """Every member and the committee mean on every test config. Returns a dict
    with the per-member and mean predictions (internal, '_' keys), metrics
    overall / per kind / per member, per-config numbers and the committee
    calibration (spread vs true error)."""
    member_preds = [predict_all(c, test_atoms) for c in calcs]
    n_cfg = len(test_atoms)
    mean_preds, spread_atom = [], []
    for i in range(n_cfg):
        ok = [mp[i] for mp in member_preds if mp[i] is not None]
        if not ok:
            mean_preds.append(None)
            spread_atom.append(None)
            continue
        fs = np.array([p["f"] for p in ok])
        f_mean = fs.mean(axis=0)
        s_list = [p["s"] for p in ok if p["s"] is not None]
        mean_preds.append({"e": float(np.mean([p["e"] for p in ok])), "f": f_mean,
                           "s": np.mean(s_list, axis=0) if len(s_list) == len(ok) and s_list else None})
        spread_atom.append(np.sqrt(np.mean(np.sum((fs - f_mean) ** 2, axis=2), axis=0)) if len(ok) > 1 else None)

    all_idx = list(range(n_cfg))
    overall = metrics(mean_preds, test_atoms, all_idx)
    if overall is None:
        return None
    f_ref_all = np.concatenate([a.arrays["REF_forces"].ravel() for a in test_atoms])
    overall["F_rel_pct"] = 100 * overall["F_rmse_meV"] / max(1000 * float(np.sqrt(np.mean(f_ref_all ** 2))), 1e-9)
    overall["F_rmse_sigma_meV"] = bootstrap_f_rmse(mean_preds, test_atoms)

    groups = {}
    for g in sorted(set(kinds)):
        groups[g] = metrics(mean_preds, test_atoms, [i for i, k in enumerate(kinds) if k == g])

    members = [metrics(mp, test_atoms, all_idx) for mp in member_preds]

    per_config = []
    off = overall["E_offset_meV"]
    for i, (a, p) in enumerate(zip(test_atoms, mean_preds)):
        if p is None:
            per_config.append(None)
            continue
        err_atom = np.linalg.norm(p["f"] - a.arrays["REF_forces"], axis=1)
        sp = spread_atom[i]
        per_config.append({
            "F_rmse_meV": _f_rmse(p, a),
            "dE_meV": 1000 * (p["e"] - a.info["REF_energy"] / len(a)) - off,
            "err_atom": err_atom, "worst_atom": int(np.argmax(err_atom)),
            "err_max_meV": 1000 * float(err_atom.max()),
            "spread_atom": sp,
            "spread_max_meV": 1000 * float(sp.max()) if sp is not None else None})

    calibration = None
    if len(calcs) > 1:
        rows = [r for r in per_config if r is not None and r["spread_max_meV"] is not None]
        spread = np.array([r["spread_max_meV"] for r in rows])
        err = np.array([r["err_max_meV"] for r in rows])
        sa = np.concatenate([r["spread_atom"] for r in rows])
        ea = np.concatenate([r["err_atom"] for r in rows])
        calibration = {
            "n_configs": len(rows), "n_atoms": int(len(sa)),
            "config_corr": _corr(spread, err) if len(rows) >= _MIN_CONFIGS_CORR else None,
            "atom_corr": _corr(sa, ea) if len(sa) >= _MIN_ATOMS_CORR else None,
            "error_over_spread": (float(np.median(err[spread >= _MIN_SPREAD_MEV] / spread[spread >= _MIN_SPREAD_MEV]))
                                  if np.any(spread >= _MIN_SPREAD_MEV) else None),
            "spread_mean_meV": float(spread.mean()), "spread_max_meV": float(spread.max())}

    return {"overall": overall, "groups": groups, "members": members, "calibration": calibration,
            "per_config": per_config, "_mean_preds": mean_preds, "_member_preds": member_preds}


def foundation_calc(model_size, device):
    from stb.core import mace_relax
    return mace_relax.get_calculator(model=model_size, device=device)


# --------------------------------------------------------------------------
# Per element, equation of state, learning curve, round comparison
# --------------------------------------------------------------------------

def element_errors(atoms_list, preds):
    """Force errors per chemical element over every config: {element:
    {n_atoms, F_rmse_meV, F_mae_meV, F_maxatom_meV (largest |dF| of one atom),
    F_ref_rms_meV (RMS SIESTA force component, the scale)}}."""
    acc = {}
    for a, p in zip(atoms_list, preds):
        if p is None:
            continue
        d = p["f"] - a.arrays["REF_forces"]
        for el, di, fr in zip(a.get_chemical_symbols(), d, a.arrays["REF_forces"]):
            e = acc.setdefault(el, ([], []))
            e[0].append(di)
            e[1].append(fr)
    out = {}
    for el, (d, fr) in sorted(acc.items()):
        d, fr = np.array(d), np.array(fr)
        out[el] = {"n_atoms": len(d), "F_rmse_meV": 1000 * float(np.sqrt(np.mean(d ** 2))),
                   "F_mae_meV": 1000 * float(np.mean(np.abs(d))),
                   "F_maxatom_meV": 1000 * float(np.linalg.norm(d, axis=1).max()),
                   "F_ref_rms_meV": 1000 * float(np.sqrt(np.mean(fr ** 2)))}
    return out


def eos_kinds(atoms_list, kinds, min_points=5, min_spread=0.01):
    """Kinds whose configs span a volume range wide enough for an E(V) fit:
    {kind: [indices]} (at least `min_points` distinct volumes -- Birch-Murnaghan
    has 4 parameters, so 4 points are fitted exactly and say nothing --, a
    relative spread of at least `min_spread`)."""
    out = {}
    for g in sorted(set(kinds)):
        ii = [i for i, k in enumerate(kinds) if k == g]
        v = np.array([atoms_list[i].get_volume() / len(atoms_list[i]) for i in ii])
        if (len(np.unique(np.round(v, 4))) >= min_points
                and (v.max() - v.min()) / v.mean() >= min_spread):
            out[g] = ii
    return out


def eos_fit_set(atoms_list, idx, energies):
    """Birch-Murnaghan fit of E/atom (eV) vs V/atom (Ang^3) on configs `idx`;
    None when a prediction is missing or the fit fails."""
    from stb.core.eos_fit import fit_eos
    if energies is None or any(e is None for e in energies):
        return None
    v = np.array([atoms_list[i].get_volume() / len(atoms_list[i]) for i in idx])
    try:
        fit = fit_eos(v, np.asarray(energies, float), EOS_FORM)
    except Exception:
        return None
    return {"v0": float(fit["v0"]), "b0_gpa": float(fit["b0_gpa"]),
            "bprime": None if fit["bprime"] is None else float(fit["bprime"]),
            "r2": float(fit["r_squared"]), "outside": bool(fit["v0_outside_range"]),
            "_curve": (np.asarray(fit["curve_v"]), np.asarray(fit["curve_e"]) - float(fit["e0"])),
            "_points": (v, np.asarray(energies, float) - float(fit["e0"]))}


def _power_law(lc):
    """'a N^-b' with the sign written once (b < 0: the error grows with N)."""
    b = lc["b"]
    return f"{lc['a']:.3g} N^{'-' if b >= 0 else '+'}{abs(b):.2f}"


def far_training_configs(train_atoms, test_atoms):
    """Training configs whose energy per atom lies more than the test set's own
    energy spread above its highest energy (compared only between configs of
    the same composition, since E per atom depends on it): [(folder name,
    meV/atom above the highest test energy)], largest first. Too-hard sampling
    shows up here -- 1-2 eV/atom picks for a test set within 0.4 eV/atom."""
    def comp(a):
        sym = a.get_chemical_symbols()
        return tuple(sorted((el, round(sym.count(el) / len(sym), 4)) for el in set(sym)))
    test_e = {}
    for a in test_atoms:
        test_e.setdefault(comp(a), []).append(a.info["REF_energy"] / len(a))
    out = []
    for a in train_atoms:
        es = test_e.get(comp(a))
        if not es or len(es) < 2:
            continue
        top, spread = max(es), max(es) - min(es)
        de = a.info["REF_energy"] / len(a) - top
        if spread > 0 and de > spread:
            out.append((os.path.basename(a.info["folder"]), 1000 * de))
    return sorted(out, key=lambda x: -x[1])


def learning_curve(points, target=None):
    """Power law F = a * N^-b through the test F RMSE of the rounds against
    their number of training configs N (the best round per N). Returns
    {'points', 'a', 'b', 'r2', 'n_target'} -- the fit keys only with 3+
    distinct N; n_target (configs needed for `target`) only for b > 0."""
    by_n = {}
    for n, f in points:
        if n and f:
            by_n[n] = min(f, by_n.get(n, np.inf))
    out = {"points": sorted(by_n.items())}
    if len(by_n) < 3:
        return out
    ns = np.array(sorted(by_n), float)
    fs = np.array([by_n[n] for n in sorted(by_n)])
    slope, inter = np.polyfit(np.log(ns), np.log(fs), 1)
    resid = np.log(fs) - (inter + slope * np.log(ns))
    ss = np.sum((np.log(fs) - np.log(fs).mean()) ** 2)
    out.update({"a": float(np.exp(inter)), "b": float(-slope),
                "r2": float(1 - np.sum(resid ** 2) / ss) if ss > 0 else 1.0, "n_target": None})
    if target and out["b"] > 0.02:
        out["n_target"] = float((out["a"] / target) ** (1 / out["b"]))
    return out


# --------------------------------------------------------------------------
# MD check: does the error grow along a trajectory run with the model?
# --------------------------------------------------------------------------

def md_check_reference(state, wf_base, dims):
    """Everything needed to write SIESTA folders of MD frames, from the
    workflow file's reference (as stb-mlffActiveLearning does): the ase.Atoms
    to start from (the supercell for dims != 1 1 1), the template .fdf and its
    file name, the calc template, the config_extra block, pseudopotentials to
    copy, the coordinate format, and a temp dir to remove. Raises ValueError."""
    import tempfile
    from pymatgen.io.ase import AseAtomsAdaptor
    from stb.core import structure_io
    from stb.core.pseudopotentials import BANKS, get_required_pseudos, resolve_pseudo_source
    from stb.mlff import resolve_reference
    from stb.mlff_active_learning import supercell_kgrid_block
    rec = (state or {}).get("reference", {})
    if not rec.get("file"):
        raise ValueError("the workflow file records no reference .fdf")
    if rec.get("calc"):
        ref_file, calc, err = resolve_reference(mlff_workflow.resolve(rec["file"], wf_base),
                                                mlff_workflow.resolve(rec["calc"], wf_base), None)
    else:
        ref_file, calc, err = resolve_reference(None, None, mlff_workflow.resolve(rec["file"], wf_base))
    if err:
        raise ValueError(err)
    ref = structure_io.read_fdf(ref_file)
    pseudo_dir = rec.get("pseudo_dir") or "."
    if pseudo_dir not in BANKS:
        pseudo_dir = mlff_workflow.resolve(pseudo_dir, wf_base)
    pseudos, missing = get_required_pseudos(sorted(ref.species), resolve_pseudo_source(pseudo_dir))
    if missing:
        raise ValueError(f"missing pseudopotential(s) for {', '.join(missing)} in '{pseudo_dir}'")
    atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(ref))
    out = {"ref_file": ref_file, "template": ref_file, "basename": os.path.basename(ref_file),
           "calc": calc, "extra": "", "pseudos": pseudos, "fractional": ref.coord_format == "fractional",
           "tmp": None, "atoms": atoms}
    if list(dims) != [1, 1, 1]:
        out["tmp"] = tempfile.mkdtemp(prefix="stb_mlffval_")
        out["template"] = os.path.join(out["tmp"], out["basename"])
        n_super = structure_io.replicate_fdf_cell(ref_file, list(dims), out["template"])
        text = ""
        for path in ([calc] if calc else []) + [ref_file]:
            with open(path) as f:
                text += f.read() + "\n"
        out["extra"], _ = supercell_kgrid_block(text, list(dims))
        if calc:
            out["extra"] += f"NumberOfAtoms {n_super}\n\n"
        out["atoms"] = atoms.repeat(list(dims))   # replicate_fdf_cell uses the same atom order
    out["extra"] += mlff_workflow.dftd3_block(state)    # the same dispersion as every folder
    return out


def run_md_check(atoms0, model_path, device, temperature, steps, n_frames, seed):
    """NVT (Langevin, 1 fs) with one model; `n_frames` frames evenly spaced
    over `steps`, each wrapped into the cell. Returns ([(step, Atoms, T_inst)],
    reason the run stopped early or None)."""
    from mace.calculators import MACECalculator
    from stb.mlmd import build_dynamics
    calc = MACECalculator(model_paths=[model_path], device=device, default_dtype="float64")
    atoms = atoms0.copy()
    atoms.calc = calc
    dyn = build_dynamics(atoms, "nvt", timestep_fs=1.0, temperature_k=temperature, friction_fs=0.01,
                         pressure_gpa=0.0, seed=seed)
    marks = {max(1, int(round(steps * k / n_frames))) for k in range(1, n_frames + 1)}
    frames, stopped = [], None
    for step in range(1, steps + 1):
        try:
            dyn.run(1)
        except Exception as e:
            stopped = f"step {step}: {e}"
            break
        if not np.all(np.isfinite(atoms.get_positions())):
            stopped = f"step {step}: non-finite positions"
            break
        if step in marks:
            fr = atoms.copy()
            fr.calc = None
            fr.wrap()
            frames.append((step, fr, float(atoms.get_temperature())))
    return frames, stopped


def compare_rounds(res_a, res_b, atoms_list, kinds):
    """Config-by-config F RMSE of two scored rounds: rows (index, F_a, F_b,
    change %) and the counts improved / worse / unchanged (10 % band)."""
    rows = []
    for i, (ra, rb) in enumerate(zip(res_a["per_config"], res_b["per_config"])):
        if ra is None or rb is None:
            continue
        fa, fb = ra["F_rmse_meV"], rb["F_rmse_meV"]
        ch = 100 * (fb - fa) / fa if fa >= _MIN_SPREAD_MEV else None
        rows.append((i, fa, fb, ch))
    better = sum(1 for r in rows if r[3] is not None and r[3] < -_COMPARE_BAND_PCT)
    worse = sum(1 for r in rows if r[3] is not None and r[3] > _COMPARE_BAND_PCT)
    return rows, better, worse, len(rows) - better - worse


# --------------------------------------------------------------------------
# Phonons against a finished stb-phononsCreate run
# --------------------------------------------------------------------------

def load_siesta_phonons(phonon_dir, f_out, info=None):
    """(SIESTA phonopy object with force constants, internal_to_angstrom),
    from a finished stb-phononsCreate folder -- the same loader stb-phononsPos
    uses. Restores the working directory."""
    from stb.core.phonon_workflow import load_phonon_with_force_constants, detect_system_label
    label = detect_system_label(phonon_dir) or "siesta"
    # symprec as stb-phononsCreate/stb-phononsPos (0.01): phonopy's own 1e-5 misses the real
    # symmetry of a DFT-relaxed cell (graphene came out Cmmm), and the few displacements
    # generated under the full symmetry are then expanded wrongly (LO/TO split at Gamma).
    phonon, to_ang, original = load_phonon_with_force_constants(phonon_dir, label, False, f_out,
                                                                symprec=PHONON_SYMPREC, info=info)
    os.chdir(original)
    phonon.symmetrize_force_constants(show_drift=False)
    return phonon, to_ang


def committee_phonons(calcs, ref_phonon, to_ang):
    """Phonopy object in Angstrom/eV with the committee-mean forces on the SAME
    displacements and supercell as `ref_phonon`, so frequencies can be compared
    q-point by q-point."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from stb.core.phonon_workflow import phonopy_atoms_to_ase
    from ase import Atoms
    uc = phonopy_atoms_to_ase(ref_phonon.unitcell, to_ang)
    unitcell = PhonopyAtoms(symbols=uc.get_chemical_symbols(), cell=uc.cell[:],
                            scaled_positions=uc.get_scaled_positions())
    ml = Phonopy(unitcell, supercell_matrix=ref_phonon.supercell_matrix,
                 primitive_matrix=ref_phonon.primitive_matrix, symprec=PHONON_SYMPREC)
    ml.dataset = {"natom": ref_phonon.dataset["natom"],
                  "first_atoms": [{"number": d["number"],
                                   "displacement": (np.array(d["displacement"]) * to_ang).tolist()}
                                  for d in ref_phonon.dataset["first_atoms"]]}
    def mean_forces(cell):
        atoms = Atoms(symbols=cell.symbols, cell=cell.cell, scaled_positions=cell.scaled_positions, pbc=True)
        fs = []
        for calc in calcs:
            a = atoms.copy()
            a.calc = calc
            fs.append(a.get_forces())
        return np.mean(fs, axis=0)

    # the model's own residual forces at the SIESTA geometry (its equilibrium differs a bit):
    # subtracted like the SIESTA ones, or they add up over the supercell at Gamma
    f0 = mean_forces(ml.supercell)
    ml.forces = np.array([mean_forces(sc) - f0 for sc in ml.supercells_with_displacements])
    ml.produce_force_constants()
    ml.symmetrize_force_constants(show_drift=False)
    return ml


def phonon_errors(ref_phonon, ml_phonon, mesh):
    """Frequency errors on the q-mesh (bands sorted per q-point). The '_'
    arrays are kept for the figure, not stored."""
    ref_phonon.run_mesh(mesh)
    q = ref_phonon.mesh.qpoints
    f_ref = np.sort(ref_phonon.mesh.frequencies, axis=1)
    ml_phonon.run_qpoints(q)
    f_ml = np.sort(ml_phonon.get_qpoints_dict()["frequencies"], axis=1)
    err = f_ml - f_ref
    n_ac = 3
    return {"mesh": list(mesh), "n_q": int(len(q)),
            "rmse_THz": float(np.sqrt(np.mean(err ** 2))),
            "acoustic_rmse_THz": float(np.sqrt(np.mean(err[:, :n_ac] ** 2))),
            "optical_rmse_THz": (float(np.sqrt(np.mean(err[:, n_ac:] ** 2)))
                                 if err.shape[1] > n_ac else None),
            "mean_shift_THz": float(np.mean(err)),
            "max_abs_THz": float(np.abs(err).max()),
            "ref_max_THz": float(f_ref.max()), "model_max_THz": float(f_ml.max()),
            "ref_min_THz": float(f_ref.min()), "model_min_THz": float(f_ml.min()),
            "model_imaginary_q": int(np.sum(f_ml.min(axis=1) < _IMAGINARY_THZ)),
            "ref_imaginary_q": int(np.sum(f_ref.min(axis=1) < _IMAGINARY_THZ)),
            "_f_ref": f_ref, "_f_ml": f_ml}


# --------------------------------------------------------------------------

def verdict(history, target_force, plateau_pct, sigma=None):
    """(code, message) from the test F RMSE of the scored rounds, oldest
    first: stop / regressed / plateau / continue. `sigma` is the bootstrap
    uncertainty of the latest round's F RMSE (or None)."""
    latest_k, latest = history[-1]
    best_k, best = min(history, key=lambda kv: kv[1])
    if target_force is not None and latest <= target_force:
        return "stop", (f"Target met: test F RMSE {latest:.1f} <= {target_force:g} meV/Ang -- "
                        "stop here and use the model.")
    if len(history) < 2:
        return "continue", ("First scored round -- run an active-learning round (4.17.4), "
                            "label it, train (4.17.2) and validate again to see the trend.")
    if best_k != latest_k and latest > best * (1 + plateau_pct / 100):
        target_txt = (f" (it meets the {target_force:g} meV/Ang target)"
                      if target_force is not None and best <= target_force else "")
        return "regressed", (f"Round {latest_k} (F {latest:.1f} meV/Ang) is "
                             f"{100 * (latest - best) / best:.0f} % WORSE than round {best_k} "
                             f"({best:.1f}){target_txt}: the newest training did not help. Use round "
                             f"{best_k}'s models, and check round {latest_k}'s settings in [2] "
                             "(epochs, starting weights, training set) and its configs far above the test "
                             "set in energy ([1]).")
    prev = history[-2][1]
    gain = 100 * (prev - latest) / prev
    noise = (f" -- within the test set's own noise (+/- {sigma:.1f} meV/Ang); more test configs "
             "would tell") if sigma is not None and abs(prev - latest) < sigma else ""
    if gain < plateau_pct:
        return "plateau", (f"Round {latest_k} changed the test F RMSE by {gain:+.1f} % "
                           f"(< {plateau_pct:g} %){noise}: more configs of the SAME kind are unlikely "
                           "to help. Change the sampling (larger --stdev, higher --md-temperature, "
                           "--supercell) or stop.")
    return "continue", (f"Round {latest_k} lowered the test F RMSE by {gain:.1f} %{noise} -- "
                        "another active-learning round is likely to help.")


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

def plot_validation(scored, detail_k, found, test_atoms, kinds, target, out_path, show, extras=None):
    """Six panels: test F RMSE per round (overall, per kind, members, foundation,
    target); force and energy parity of the detailed round vs the foundation
    model; true error vs committee spread per config; F RMSE per kind; and the
    phonon frequency parity (with --phonon-ref) or the force error distribution.
    A third row holds whichever of `extras` is available: the learning curve,
    the E(V) fits, the MD check (error vs time) and the --compare scatter."""
    res = dict(scored)[detail_k]
    lab_k = "fine-tuned" if detail_k == "models" else f"round {detail_k}"
    kcol = {g: f"C{n % 10}" for n, g in enumerate(sorted(set(kinds)))}   # same color per kind in every panel
    extras = {k: v for k, v in (extras or {}).items() if v}
    nrows = 2 + (len(extras) + 2) // 3
    fig, axes = plt.subplots(nrows, 3, figsize=(16, 4.5 * nrows), squeeze=False)

    ax = axes[0, 0]
    keys = [k for k, _ in scored]
    xs = np.arange(len(keys))
    for x, (_, v) in zip(xs, scored):
        mf = [m["F_rmse_meV"] for m in v["members"] if m]
        ax.scatter([x] * len(mf), mf, s=14, color="0.65", zorder=2,
                   label="members" if x == 0 else None)
    yo = [v["overall"]["F_rmse_meV"] for _, v in scored]
    sg = [v["overall"].get("F_rmse_sigma_meV") or 0 for _, v in scored]
    ax.errorbar(xs, yo, yerr=sg, fmt="o-", color="k", lw=2, capsize=3, zorder=3,
                label="committee mean (+/- bootstrap)")
    groups = sorted({g for _, v in scored for g in v["groups"]})
    if len(groups) > 1:
        for g in groups:
            ys = [(v["groups"].get(g) or {}).get("F_rmse_meV", np.nan) for _, v in scored]
            if np.nanmax(ys) < _MIN_SPREAD_MEV:   # equilibrium configs (EOS): no forces to fit
                continue
            ax.plot(xs, ys, "o--", ms=4, color=kcol.get(g), label=g)
    if found:
        ax.axhline(found["overall"]["F_rmse_meV"], color="tab:orange", ls=":", label="foundation")
    if target:
        ax.axhline(target, color="tab:green", ls="-.", label="target")
    ax.set_xticks(xs)
    ax.set_xticklabels([str(k) for k in keys])
    ax.set_xlabel("Round")
    ax.set_ylabel("Test F RMSE (meV/Ang)")
    ax.set_yscale("log")
    ax.set_title("Force error on the test set")
    ax.legend(fontsize=7)

    mp = res["_mean_preds"]
    idx = [i for i, p in enumerate(mp) if p is not None]
    fp = found["_preds"] if found else None
    fidx = [i for i in range(len(test_atoms)) if fp and fp[i] is not None]
    ax = axes[0, 1]
    if fidx:
        ax.scatter(np.concatenate([test_atoms[i].arrays["REF_forces"].ravel() for i in fidx]),
                   np.concatenate([fp[i]["f"].ravel() for i in fidx]), s=5, alpha=0.35,
                   color="tab:orange", label="foundation")
    ax.scatter(np.concatenate([test_atoms[i].arrays["REF_forces"].ravel() for i in idx]),
               np.concatenate([mp[i]["f"].ravel() for i in idx]), s=7, alpha=0.6,
               color="tab:blue", label=lab_k)
    lims = [*ax.get_xlim(), *ax.get_ylim()]
    ax.plot([min(lims), max(lims)], [min(lims), max(lims)], "--", color="gray", lw=0.8)
    ax.set_xlabel("SIESTA force (eV/Ang)")
    ax.set_ylabel("Predicted force (eV/Ang)")
    ax.set_title(f"Force parity, test set ({lab_k})")
    ax.legend(fontsize=8)

    ax = axes[0, 2]
    ref_e = lambda ii: np.array([test_atoms[i].info["REF_energy"] / len(test_atoms[i]) for i in ii])
    for ii, preds, col, lab in ((fidx, fp, "tab:orange", "foundation"),
                                (idx, mp, "tab:blue", lab_k)):
        if len(ii) > 1:
            r, p = ref_e(ii), np.array([preds[i]["e"] for i in ii])
            ax.scatter(1000 * (r - r.mean()), 1000 * (p - p.mean()), s=24, alpha=0.7, color=col, label=lab)
    lims = [*ax.get_xlim(), *ax.get_ylim()]
    ax.plot([min(lims), max(lims)], [min(lims), max(lims)], "--", color="gray", lw=0.8)
    ax.set_xlabel("SIESTA E - mean (meV/atom)")
    ax.set_ylabel("Predicted E - mean (meV/atom)")
    ax.set_title("Energy parity, offset removed")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)

    ax = axes[1, 0]
    rows = [(i, r) for i, r in enumerate(res["per_config"]) if r is not None]
    cal = res["calibration"]
    if cal:   # equilibrium configs (zero force, zero spread) say nothing about calibration
        rows = [(i, r) for i, r in rows if (r["spread_max_meV"] or 0) >= _MIN_SPREAD_MEV] or rows
    if cal:
        sp = np.maximum([r["spread_max_meV"] for _, r in rows], 1e-2)
        er = np.maximum([r["err_max_meV"] for _, r in rows], 1e-2)
        for g in sorted(set(kinds)):
            sel = [k for k, (i, _) in enumerate(rows) if kinds[i] == g]
            if sel:
                ax.scatter(sp[sel], er[sel], s=30, alpha=0.8, color=kcol[g], label=g)
        lo, hi = max(min(sp.min(), er.min()) * 0.5, 1e-2), max(sp.max(), er.max()) * 2
        ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=0.8, label="error = spread")
        k = cal["error_over_spread"]
        if k:
            ax.plot([lo, hi], [k * lo, k * hi], ":", color="tab:red", lw=1, label=f"error = {k:.1f} x spread")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Committee spread, max over atoms (meV/Ang)")
        ax.set_ylabel("True error, max over atoms (meV/Ang)")
        ax.set_title("Calibration: is the spread an error bar?")
        ax.legend(fontsize=7)
    else:
        fmax = np.array([np.linalg.norm(test_atoms[i].arrays["REF_forces"], axis=1).max() for i, _ in rows])
        ax.scatter(fmax, [r["F_rmse_meV"] for _, r in rows], s=30)
        ax.set_xlabel("SIESTA max |F| (eV/Ang)")
        ax.set_ylabel("Test F RMSE per config (meV/Ang)")
        ax.set_title("Error vs force magnitude (single model)")

    ax = axes[1, 1]
    names = sorted(set(kinds))
    x = np.arange(len(names))
    for off, (lab, gr, col) in enumerate((("foundation", (found or {}).get("groups"), "tab:orange"),
                                          (lab_k, res["groups"], "tab:blue"))):
        if not gr:
            continue
        ax.bar(x + (off - 0.5) * 0.38, [(gr.get(n) or {}).get("F_rmse_meV", np.nan) for n in names],
               0.38, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, fontsize=8)
    ax.set_ylabel("Test F RMSE (meV/Ang)")
    ax.set_title("Force error per configuration kind")
    ax.legend(fontsize=8)

    ax = axes[1, 2]
    ph = res.get("phonons")
    if ph:
        ax.scatter(ph["_f_ref"].ravel(), ph["_f_ml"].ravel(), s=5, alpha=0.4)
        lims = [*ax.get_xlim(), *ax.get_ylim()]
        ax.plot([min(lims), max(lims)], [min(lims), max(lims)], "--", color="gray", lw=0.8)
        ax.set_xlabel("SIESTA frequency (THz)")
        ax.set_ylabel("Committee-mean frequency (THz)")
        ax.set_title(f"Phonons on a {'x'.join(map(str, ph['mesh']))} q-mesh (RMSE {ph['rmse_THz']:.2f} THz)")
    else:
        errs = []
        if fidx:
            errs.append(("foundation", 1000 * np.abs(np.concatenate(
                [(fp[i]["f"] - test_atoms[i].arrays["REF_forces"]).ravel() for i in fidx])), "tab:orange", 0.5))
        errs.append((lab_k, 1000 * np.abs(np.concatenate(
            [(mp[i]["f"] - test_atoms[i].arrays["REF_forces"]).ravel() for i in idx])), "tab:blue", 0.7))
        top = max(e.max() for _, e, _, _ in errs)
        lo = max(1e-1, min(max(e.min(), 1e-1) for _, e, _, _ in errs))
        bins = np.logspace(np.log10(lo), np.log10(max(top, lo * 10) * 1.05), 35)
        for lab, e, col, al in errs:
            ax.hist(np.maximum(e, bins[0]), bins=bins, alpha=al, color=col, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel("|F error| per component (meV/Ang, log)")
        ax.set_ylabel("Count")
        ax.set_title("Force error distribution, test set")
        ax.legend(fontsize=8)

    slots = [axes[2 + j // 3, j % 3] for j in range(3 * (nrows - 2))]
    for ax, (name, data) in zip(slots, extras.items()):
        if name == "learning":
            pts = np.array(data["points"], float)
            ax.loglog(pts[:, 0], pts[:, 1], "o", color="k", label="rounds (best per N)")
            if "b" in data:
                n_hi = max(pts[:, 0].max(), data.get("n_target") or 0) * 1.3
                nn = np.geomspace(pts[:, 0].min() * 0.8, n_hi, 50)
                ax.loglog(nn, data["a"] * nn ** -data["b"], "--", color="tab:blue",
                          label=f"F = {_power_law(data)} (R2 {data['r2']:.2f})")
            if target:
                ax.axhline(target, color="tab:green", ls="-.", label="target")
                if data.get("n_target"):
                    ax.axvline(data["n_target"], color="tab:green", ls=":", lw=0.8)
            ax.set_xlabel("Training configs N")
            ax.set_ylabel("Test F RMSE (meV/Ang)")
            ax.set_title("Learning curve")
            ax.legend(fontsize=7)
        elif name == "eos":
            kind, fits = data
            for lab, fit, col, mk in (("SIESTA", fits.get("siesta"), "k", "o"),
                                      ("foundation", fits.get("foundation"), "tab:orange", "s"),
                                      (lab_k, fits.get("detail"), "tab:blue", "^")):
                if fit:
                    ax.plot(fit["_curve"][0], 1000 * fit["_curve"][1], "-", color=col, lw=1.2,
                            label=f"{lab}: B0 {fit['b0_gpa']:.1f} GPa")
                    ax.plot(fit["_points"][0], 1000 * fit["_points"][1], mk, color=col, ms=5)
            ax.set_xlabel("V/atom (Ang^3)")
            ax.set_ylabel("E - E0 (meV/atom)")
            ax.set_title(f"Equation of state, test configs '{kind}'")
            ax.legend(fontsize=7)
        elif name == "md":
            for run, pts in data.items():
                pts = sorted(pts)
                ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=4, label=run)
            if res["overall"]:
                ax.axhline(res["overall"]["F_rmse_meV"], color="gray", ls=":", label=f"test set ({lab_k})")
            ax.set_xlabel("MD time (fs)")
            ax.set_ylabel(f"F RMSE vs SIESTA, {lab_k} (meV/Ang)")
            ax.set_title("MD check: error along the trajectory")
            ax.legend(fontsize=7)
        elif name == "compare":
            ka, kb, rows = data
            fa = np.maximum([r[1] for r in rows], 1e-2)
            fb = np.maximum([r[2] for r in rows], 1e-2)
            for g in sorted(set(kinds)):
                sel = [j for j, r in enumerate(rows) if kinds[r[0]] == g]
                if sel:
                    ax.scatter(fa[sel], fb[sel], s=28, alpha=0.8, color=kcol[g], label=g)
            lo, hi = min(fa.min(), fb.min()) * 0.7, max(fa.max(), fb.max()) * 1.4
            ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=0.8)
            ax.fill_between([lo, hi], [lo * 0.9, hi * 0.9], [lo * 1.1, hi * 1.1], color="0.9", zorder=0)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel(f"F RMSE per config, round {ka} (meV/Ang)")
            ax.set_ylabel(f"F RMSE per config, round {kb} (meV/Ang)")
            ax.set_title(f"Round {kb} vs round {ka} (below the line = better)")
            ax.legend(fontsize=7)
    for ax in slots[len(extras):]:
        ax.axis("off")

    from matplotlib.ticker import NullFormatter
    for ax in axes.ravel():   # minor-tick labels (2x10^1, 3x10^1, ...) overlap on short log axes
        for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if scale == "log":
                axis.set_minor_formatter(NullFormatter())
    fig.suptitle("stb-mlffValidate -- held-out SIESTA test set", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 1 - 0.035 * 2 / nrows])
    fig.savefig(out_path, dpi=140)
    if show:
        plt.show()
    plt.close(fig)


def can_show():
    """True when matplotlib has an interactive backend (a display) to show a figure on."""
    import matplotlib
    return matplotlib.get_backend().lower() not in ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template")


# --------------------------------------------------------------------------

def _kind_from_name(folder):
    """Kind of a test folder that has no workflow entry, from its name:
    'large_01' -> 'large', 'eos_+02' -> 'eos', 'test_003' -> 'test'."""
    import re
    m = re.match(r"^(.*?)[_-][+-]?[A-Za-z]?\d+(\.\d+)?$", os.path.basename(os.path.normpath(folder)))
    return m.group(1) if m and m.group(1) else "test"


def _model_labels(paths):
    """Short distinct labels for model files (parent folder when the file names repeat)."""
    base = [os.path.basename(p) for p in paths]
    if len(set(base)) == len(base):
        return base
    return [os.path.basename(os.path.dirname(os.path.abspath(p))) or b for p, b in zip(paths, base)]


def _siesta_loop(run_dir, example_folder):
    """The shell loop that runs SIESTA in every frame folder of `run_dir`."""
    from stb.core import structure_io
    fdfs = sorted(f for f in os.listdir(example_folder) if f.endswith(".fdf") and f != "config_extra.fdf")
    main_fdf = next((f for f in fdfs if structure_io.get_system_label(os.path.join(example_folder, f))), None)
    main_fdf = next((f for f in fdfs if f.startswith("calc")), main_fdf) or (fdfs[0] if fdfs else "calc.fdf")
    label = structure_io.find_fdf_system_label(example_folder) or "siesta"
    return (f'for d in {os.path.relpath(run_dir)}/frame_*; do (cd "$d" && siesta < {main_fdf} > {label}.out); done')


def _n_train(entry):
    """Training configs a round actually used."""
    entry = entry or {}
    return entry.get("n_usable") or len(entry.get("configs", [])) or None


def _round_label(k):
    return f"round {k}" if k != "models" else "--model"


BEST_DIR = "best_model"   # next to the workflow file; rewritten by every Stage 3 run


def export_best_model(state, rounds, best_k, wf_base, workflow_path, verdict_code, n_test, warn):
    """Copies the best model of this run, with its card, to <run folder>/best_model/ -- replacing
    whatever the previous run left there -- and records it as state['best'].

    Best = the round (among `rounds`, the (number, scores) pairs scored now) with the lowest
    committee-mean test force RMSE, and within it the single member with the lowest test force
    RMSE; a member that failed Stage 2's MD stability check only wins if every member did. The
    card in the folder is the member's own, plus a 'selection' section saying why it was chosen.
    Returns the info dict (also state['best']), or None when nothing could be exported."""
    v = dict(rounds)[best_k]
    entry = mlff_workflow.get_round(state, best_k)
    by_seed = {m["seed"]: m for m in entry.get("committee", []) if m.get("model")}
    stab = {r.get("seed"): r.get("status") for r in (entry.get("model_checks") or {}).get("stability", [])}
    cands = [(stab.get(sd) not in (None, "OK"), met["F_rmse_meV"], sd)
             for sd, met in zip(v["seeds"], v["members"]) if met and sd in by_seed]
    if not cands:
        warn(f"Best model: round {best_k} has no scored member with a model file -- {BEST_DIR}/ not written.")
        return None
    _, member_f, seed = min(cands)
    member = by_seed[seed]
    src_model = mlff_workflow.resolve(member["model"], wf_base)
    card_rel = member.get("card")
    card_src = mlff_workflow.resolve(card_rel, wf_base) if card_rel else None
    if not card_src or not os.path.isfile(card_src):
        try:
            rebuild_round_cards(state, entry, wf_base, workflow_path)
            card_src = mlff_workflow.resolve(member["card"], wf_base)
        except Exception as e:
            warn(f"Best model: round {best_k} seed {seed} has no model card and it could not be written "
                 f"({e}) -- {BEST_DIR}/ not written; run stb-mlffAnalysis --make-cards {best_k}.")
            return None
    try:
        card = model_card.load_card(card_src)
    except model_card.ModelCardError as e:
        warn(f"Best model: {e} -- {BEST_DIR}/ not written.")
        return None
    src_dir = os.path.dirname(src_model)
    dest = os.path.join(wf_base, BEST_DIR)
    os.makedirs(dest, exist_ok=True)
    for old in glob.glob(os.path.join(dest, "*.model")) + [os.path.join(dest, model_card.CARD_FILE)]:
        if os.path.isfile(old):
            os.remove(old)
    for f in card["model"]["files"]:
        src = os.path.join(src_dir, f["file"])
        if not os.path.isfile(src):
            continue
        shutil.copy2(src, os.path.join(dest, f["file"]))
        if model_card.sha256_file(os.path.join(dest, f["file"])) != f["sha256"]:
            warn(f"Best model: {f['file']} of round {best_k} seed {seed} is not the file its card describes "
                 "(retrained or replaced since the card was written?) -- run stb-mlffAnalysis --make-cards "
                 f"{best_k}.")
    others = sorted(r["round"] for r in state["rounds"] if r["round"] not in {k for k, _ in rounds})
    if others:
        warn(f"Best model: {BEST_DIR}/ was chosen among the rounds scored in this run only; round(s) "
             f"{', '.join(map(str, others))} were not scored -- run stb-mlffValidate without --rounds to compare all.")
    when = datetime.now().isoformat(timespec="seconds")
    card["selection"] = {
        "selected_by": "stb-mlffValidate", "time": when,
        "criterion": "round with the lowest committee-mean test-set force RMSE among the rounds scored in that "
                     "run; within it the member with the lowest test-set force RMSE (a member that failed the "
                     "MD stability check only if all did)",
        "round": best_k, "seed": seed, "source_model": mlff_workflow.rel(src_model, wf_base),
        "source_card": mlff_workflow.rel(card_src, wf_base), "n_test": n_test,
        "member_test": next(m for sd, m in zip(v["seeds"], v["members"]) if sd == seed),
        "committee_mean_test": v["overall"],
        "rounds_compared": [{"round": k, "F_rmse_meV": r["overall"]["F_rmse_meV"]} for k, r in rounds],
        "latest_round": rounds[-1][0], "latest_round_verdict": verdict_code,
        "note": "paths in 'committee' are relative to the workflow file, not to this folder"}
    model_card.write_card(card, os.path.join(dest, model_card.CARD_FILE))
    info = {"time": when, "round": best_k, "seed": seed, "folder": BEST_DIR,
            "model": f"{BEST_DIR}/{card['model']['file']}", "card": f"{BEST_DIR}/{model_card.CARD_FILE}",
            "source_model": mlff_workflow.rel(src_model, wf_base),
            "test_F_rmse_meV": round(float(member_f), 3),
            "committee_test_F_rmse_meV": round(float(v["overall"]["F_rmse_meV"]), 3),
            "verdict": verdict_code}
    state["best"] = info
    return info


def main():
    parser = argparse.ArgumentParser(
        description="Stage 3 (Validate) of the Custom ML Force Field workflow: scores every "
                     "training round recorded in mlff_workflow.json (its committee mean and "
                     "each member) against the held-out SIESTA test set from stb-mlff "
                     "--n-test -- energy, force and stress errors, overall, per kind and per "
                     "config, with a bootstrap error bar --, compares the test set with the "
                     "training data (level of theory, range covered), checks whether the "
                     "committee's disagreement tracks the real error, optionally compares "
                     "phonons against a finished stb-phononsCreate run, and says whether "
                     "another active-learning round is worth it. Every warning goes to one "
                     "final section. Reads finished calculations only -- except --md-check, which writes SIESTA "
                     "folders of MD frames to score on the next run.",
        epilog="Example usage:\n"
               "  stb-mlffValidate                         # every round in mlff_workflow.json\n"
               "  stb-mlffValidate --show --save-data      # display the figure, write per-config errors\n"
               "  stb-mlffValidate --target-force 30 --phonon-ref phonon_runs\n"
               "  stb-mlffValidate --view-worst 3          # worst test configs in ASE's viewer\n"
               "  stb-mlffValidate --model a.model b.model --test-path 'mlff_test/test_*'",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--workflow", default=None, metavar="JSON",
                        help=f"Workflow state file (default: {mlff_workflow.RUN_DIR}/"
                             f"{mlff_workflow.WORKFLOW_FILE}, or ./{mlff_workflow.WORKFLOW_FILE} "
                             "inside that folder). Rounds "
                             "and test set are read from it; each round's metrics are stored back.")
    parser.add_argument("--rounds", type=int, nargs="+", default=None, metavar="K",
                        help="Rounds to score (default: all rounds in --workflow).")
    parser.add_argument("--detail-round", type=int, default=None, metavar="K",
                        help="Round shown per kind, per config and per member, and in the figure "
                             "(default: the latest scored round).")
    parser.add_argument("--model", nargs="+", default=None, metavar="PATH",
                        help="Score these model files as one committee instead of the workflow's "
                             "rounds (nothing is stored back).")
    parser.add_argument("--test-path", default=None, metavar="GLOB",
                        help="Test config folders (default: the test set in --workflow, else "
                             "'mlff_test/test_*'). Quote the pattern.")
    parser.add_argument("--phonon-ref", default=None, metavar="DIR",
                        help="A finished stb-phononsCreate folder (phonopy_disp.yaml + disp-*/ with "
                             "SIESTA forces): the committee-mean model is run on the SAME "
                             "displacements and its frequencies compared q-point by q-point.")
    parser.add_argument("--phonon-mesh", type=int, nargs=3, default=[8, 8, 8], metavar=("NA", "NB", "NC"),
                        help="q-point mesh for --phonon-ref (default: 8 8 8).")
    parser.add_argument("--target-force", type=float, default=None, metavar="MEV_ANG",
                        help="Test F RMSE (meV/Ang) that is good enough for the intended use: the "
                             "verdict says stop once the latest round reaches it.")
    parser.add_argument("--plateau", type=float, default=DEFAULT_PLATEAU_PCT, metavar="PCT",
                        help="A round that lowers the test F RMSE by less than this, in %%, is "
                             f"reported as a plateau (default: {DEFAULT_PLATEAU_PCT:g}).")
    parser.add_argument("--foundation-model", default=None, choices=["small", "medium", "large"],
                        help="Foundation model for the baseline row (default: the one the latest "
                             "round was trained from, else medium).")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="Do not score the raw foundation model on the test set.")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"],
                        help="Device to run the models on (default: cpu).")
    parser.add_argument("--png", default="mlff_validation.png", metavar="FILE",
                        help="Summary figure (default: mlff_validation.png, next to --workflow).")
    parser.add_argument("--show", action="store_true",
                        help="Also display the summary figure on screen (it is always saved).")
    parser.add_argument("--save-data", action="store_true",
                        help=f"Also write {DATA_FILE} (next to --workflow): per test config of the detailed "
                             "round, its kind, SIESTA max force, model and foundation F RMSE, energy error, "
                             "committee spread, novelty. Off by default.")
    parser.add_argument("--view-worst", type=int, default=0, metavar="N",
                        help=f"Write the N worst-predicted test configs to {WORST_FILE} (per-atom force "
                             "error as 'F_error' and as the initial charge, so ASE's viewer can color by "
                             "it: View > Colors > By initial charge) and open them. Default: 0.")
    parser.add_argument("--compare", type=int, nargs=2, default=None, metavar=("A", "B"),
                        help="Compare two rounds config by config (which test configs got better or worse "
                             "from round A to round B, per kind too). Both are scored even if --rounds "
                             "leaves them out.")
    parser.add_argument("--md-check", action="store_true",
                        help="Run an NVT MD with the best round's model and write evenly spaced frames as "
                             f"SIESTA single-point folders in <run folder>/{MD_CHECK_DIR}/ (recorded in the "
                             "workflow file as their own set, not the test set). Run SIESTA in them; every "
                             "later stb-mlffValidate then reports the error along the trajectory -- the test "
                             "closest to real use. The one option that writes folders.")
    parser.add_argument("--md-check-temperature", type=float, default=300.0, metavar="K",
                        help="MD check temperature (default: 300 K).")
    parser.add_argument("--md-check-steps", type=int, default=2000,
                        help="MD check length in 1 fs steps (default: 2000).")
    parser.add_argument("--md-check-frames", type=int, default=8,
                        help="Frames written for SIESTA, evenly spaced (default: 8).")
    parser.add_argument("--md-check-supercell", type=int, nargs=3, default=[1, 1, 1], metavar=("NA", "NB", "NC"),
                        help="Run the MD check in this supercell of the reference (default: 1 1 1); the "
                             "k-grid is divided accordingly, as in stb-mlffActiveLearning --supercell.")
    parser.add_argument("--md-check-to-train", action="store_true",
                        help="Move the finished MD-check frames into the training set (they are already "
                             "SIESTA-labeled): the next stb-mlffAnalysis trains on them.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE} (next to --workflow). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlffValidate {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Validate)",
            "Score each training round against a held-out SIESTA test set",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])
    print("\n" + color_text("CUSTOM ML FORCE FIELD -- VALIDATE:", 'bold'))
    print("-" * 60)
    t_start = time.monotonic()

    warns = []           # every warning, printed together in the last section
    library_noise = []   # third-party prints/warnings, counted there too

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

    args.workflow = mlff_workflow.default_path(args.workflow)
    try:
        state = mlff_workflow.load(args.workflow)
    except ValueError as e:
        fail(str(e))
    wf_base = mlff_workflow.base_dir(args.workflow)
    out_base = wf_base if os.path.isdir(wf_base) else os.getcwd()
    near_wf = lambda name: name if os.path.isabs(name) else os.path.relpath(os.path.join(out_base, name))

    report_path = near_wf(REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    # ---- what to score
    committees = []   # (round number or "models", [model paths], round entry or None)
    if args.model:
        missing = [m for m in args.model if not os.path.isfile(m)]
        if missing:
            fail(f"Model file(s) not found: {', '.join(missing)}")
        committees.append(("models", args.model, None))
    else:
        if not state or not state["rounds"]:
            fail(f"No training rounds in '{args.workflow}' -- run stb-mlffAnalysis (4.17.2) first, "
                 "or pass --model.")
        numbers = sorted(set((args.rounds or [r["round"] for r in state["rounds"]]) + list(args.compare or [])))
        for k in numbers:
            entry = mlff_workflow.get_round(state, k)
            if entry is None:
                fail(f"Round {k} is not in '{args.workflow}' (it has round(s) "
                     f"{', '.join(str(r['round']) for r in state['rounds'])}).")
            listed = mlff_workflow.committee_models(entry, wf_base)
            models = [m for m in listed if os.path.isfile(m)]
            if not models:
                warn(f"Setup: round {k} has no model file on disk -- skipped.")
                continue
            if len(models) < len(listed):
                warn(f"Setup: round {k}: {len(listed) - len(models)} of its {len(listed)} model file(s) "
                     "are missing -- scored with the rest.")
            committees.append((k, models, entry))
        if not committees:
            fail("None of the requested rounds has a model file.")
    keys = [k for k, _, _ in committees]
    if args.compare:
        if args.model:
            fail("--compare needs the workflow's rounds, not --model.")
        missing = [k for k in args.compare if k not in keys]
        if missing:
            fail(f"--compare: round(s) {', '.join(map(str, missing))} could not be scored.")
        if args.compare[0] == args.compare[1]:
            fail("--compare needs two different rounds.")
    if (args.md_check or args.md_check_to_train) and not state:
        fail("--md-check needs the workflow file (reference structure, pseudopotentials, rounds).")
    detail_k = args.detail_round if args.detail_round is not None else keys[-1]
    if detail_k not in keys:
        fail(f"--detail-round {detail_k} is not among the scored rounds ({', '.join(map(str, keys))}).")

    if args.test_path:
        test_dirs = sorted(d for d in glob.glob(args.test_path) if os.path.isdir(d))
        test_source = f"--test-path '{args.test_path}'"
    elif state and mlff_workflow.folders(state, wf_base, "test"):
        test_dirs = [d for d in mlff_workflow.folders(state, wf_base, "test") if os.path.isdir(d)]
        test_source = f"{args.workflow} (test set)"
        on_disk = sorted(d for d in glob.glob(os.path.join(wf_base, "mlff_test", "test_*")) if os.path.isdir(d))
        listed = {os.path.abspath(d) for d in test_dirs}
        stale = [d for d in on_disk if os.path.abspath(d) not in listed]
        if stale:
            warn(f"Test set: {len(stale)} folder(s) in mlff_test/ are not in the workflow file and were not "
                 f"scored ({', '.join(os.path.basename(d) for d in stale[:6])}{' ...' if len(stale) > 6 else ''}) "
                 "-- left over from an earlier Stage 1 run? Delete them, or include them with "
                 f"--test-path '{os.path.join(os.path.relpath(wf_base), 'mlff_test', 'test_*')}'.")
    else:
        test_dirs = sorted(d for d in glob.glob(os.path.join(wf_base, "mlff_test", "test_*")) if os.path.isdir(d))
        test_source = "'mlff_test/test_*'"

    size = args.foundation_model
    if size is None:
        last = committees[-1][2]
        size = (last or {}).get("settings", {}).get("foundation_model", "medium")

    print_dual(color_text("===== STB-MLFFVALIDATE REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    vers = _versions()
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Versions          : stb {VERSION}, mace {vers['mace']}, torch {vers['torch']}, ase {vers['ase']}", f_out)
    print_dual(f"Device            : {args.device}", f_out)
    print_dual(f"Workflow file     : {args.workflow}" + ("" if state else " (not found)"), f_out)
    print_dual(f"Test set from     : {test_source}", f_out)
    print_dual("Scoring           : " + ", ".join(
        f"round {k} ({len(m)} model(s))" if k != "models" else f"{len(m)} model(s) from --model"
        for k, m, _ in committees), f_out)
    print_dual(f"Detailed round    : {_round_label(detail_k)} (per kind, per config, members, figure)", f_out)
    print_dual("Baseline          : " + ("none (--skip-foundation-comparison)" if args.skip_foundation_comparison
                                         else f"MACE-MP {size}, no fine-tuning"), f_out)
    print_dual(f"Verdict criteria  : target F RMSE {_fmt(args.target_force, '{:g} meV/Ang')}, "
               f"plateau below {args.plateau:g} % gain per round", f_out)
    if args.phonon_ref:
        print_dual(f"Phonon reference  : {args.phonon_ref} (mesh {' '.join(map(str, args.phonon_mesh))})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    # ---- [1] test set
    print_section("[1] TEST SET", f_out)
    if not test_dirs:
        fail("No test configs found -- generate them with stb-mlff --n-test N (and run SIESTA "
             "in them), or pass --test-path.")
    issues = []
    test_atoms, kinds, entries = [], [], []
    for d in test_dirs:   # one at a time, so each config keeps its folder (and kind)
        found_, _, _ = collect_training_configs([d], f_out, issues)
        for a in found_:
            entry = mlff_workflow.config_entry(state, d, wf_base) if state else None
            test_atoms.append(a)
            kinds.append(entry.get("kind", "test") if entry else _kind_from_name(d))
            entries.append(entry or {})
    for msg in issues:
        if "missing .out/.FA/.XV" not in msg:
            warn(f"Test set: {msg}")
    for parent, n, k in mlff_workflow.pending_by_dir(test_dirs):
        if k:
            warn(f"SIESTA: {k} of {n} test folder(s) in {os.path.relpath(parent)}/ have no finished run -- run SIESTA "
                 "there; they are not scored until then.")
    if not test_atoms:
        fail("No finished test calculation -- run SIESTA in the test folders first.")
    n_atoms = sum(len(a) for a in test_atoms)
    print_dual(f"Configs           : {len(test_atoms)} finished of {len(test_dirs)} folder(s), {n_atoms} atoms, "
               + ", ".join(f"{g} {kinds.count(g)}" for g in sorted(set(kinds))), f_out)
    elements = sorted({s for a in test_atoms for s in a.get_chemical_symbols()})
    print_dual(f"Elements          : {', '.join(elements)}; "
               f"{min(len(a) for a in test_atoms)}-{max(len(a) for a in test_atoms)} atoms per config", f_out)
    has_stress = all("REF_stress" in a.info for a in test_atoms)
    print_dual(f"Labels            : energy, forces" + (", stress" if has_stress else
               f" (stress in {sum('REF_stress' in a.info for a in test_atoms)}/{len(test_atoms)})"), f_out)
    ds_test = dataset_table(test_atoms, kinds, f_out)
    not_conv = [os.path.basename(a.info["folder"]) for a in test_atoms if a.info.get("scf_converged") is False]
    if not_conv:
        warn(f"Test set: SCF did not converge in {len(not_conv)} config(s) ({', '.join(not_conv)}): their "
             "reference energies/forces are unreliable -- rerun them or remove them from the test set.")
    if len(test_atoms) < 5:
        warn(f"Test set: only {len(test_atoms)} config(s) -- differences between rounds of a few % are within "
             "the noise (see the +/- column in [3]); generate more with stb-mlff --n-test.")
    if not has_stress:
        warn("Test set: stress missing in some configs -- the stress RMSE covers only the others.")

    # training data of the detailed round: overlap, level of theory, range
    train_atoms = []
    detail_entry = dict((k, e) for k, _, e in committees)[detail_k]
    if state:
        trained = {c for r in state["rounds"] for c in r.get("configs", [])}
        overlap = [os.path.basename(a.info["folder"]) for a in test_atoms
                   if mlff_workflow.rel(a.info["folder"], wf_base) in trained]
        if overlap:
            warn(f"Test set: {len(overlap)} test folder(s) were also trained on ({', '.join(overlap[:6])}) -- "
                 "their scores are optimistic.")
        if detail_entry:
            tdirs = [os.path.join(wf_base, c) for c in detail_entry.get("configs", [])]
            train_atoms, _, _ = collect_training_configs([d for d in tdirs if os.path.isdir(d)], f_out, [])
    if train_atoms:
        print_dual(f"Compared with     : the {len(train_atoms)} training config(s) of round {detail_k}", f_out)
        th_common, th_odd = check_theory(train_atoms + test_atoms, f_out,
                                         lambda m: warn(m.replace("Data: ", "Level of theory (train + test): ", 1)))
        dispersion_check(th_common, th_odd, state, f_out, warn)
        tr_vol = np.array([a.get_volume() / len(a) for a in train_atoms])
        tr_fmax = np.array([np.linalg.norm(a.arrays["REF_forces"], axis=1).max() for a in train_atoms])
        tr_e = np.array([a.info["REF_energy"] / len(a) for a in train_atoms])
        print_dual(f"Training range    : V/atom {tr_vol.min():.3f}-{tr_vol.max():.3f} Ang^3, max|F| up to "
                   f"{tr_fmax.max():.2f} eV/Ang, E spread {1000 * (tr_e.max() - tr_e.min()):.1f} meV/atom", f_out)
        outside = []
        for i, a in enumerate(test_atoms):
            why = []
            v = ds_test["vol"][i]
            if v < tr_vol.min() * 0.995 or v > tr_vol.max() * 1.005:
                why.append("V")
            if ds_test["fmax"][i] > tr_fmax.max() * 1.1:
                why.append("F")
            if ds_test["e_atom"][i] > tr_e.max() + 0.1 * max(tr_e.max() - tr_e.min(), 1e-3):
                why.append("E")
            outside.append("".join(why))
        n_out = sum(1 for o in outside if o)
        txt = (f"{len(test_atoms) - n_out} of {len(test_atoms)} test configs inside it"
               + (f"; outside: " + ", ".join(f"{os.path.basename(a.info['folder'])} ({o})"
                                            for a, o in zip(test_atoms, outside) if o)[:200] if n_out else ""))
        print_dual(f"Coverage          : {txt}  (V volume, F max force, E energy)", f_out)
        if n_out:
            warn(f"Test set: {n_out} test config(s) lie outside the training range (marked in [4]): their errors "
                 "measure extrapolation, not interpolation. Expected for a harder test set; if those "
                 "configurations matter, add similar ones to training (4.17.4 or stb-mlff).")
        far = far_training_configs(train_atoms, test_atoms)
        if far:
            print_dual(f"Far from the test : {len(far)} training config(s) lie more than one test-set energy spread "
                       f"above the highest test energy: " + ", ".join(
                           f"{n} (+{de:.0f})" for n, de in far[:6]) + (" ..." if len(far) > 6 else "")
                       + " meV/atom", f_out)
            warn(f"Training set: {len(far)} config(s) of round {detail_k} lie far above every test config in energy "
                 f"(up to +{max(de for _, de in far):.0f} meV/atom above the highest test energy; listed in [1]) -- "
                 "the model spends capacity on geometries that the tests (and phonons or MD near equilibrium) never "
                 "visit, which can make it WORSE near equilibrium. Usually an active-learning round sampled too "
                 "hard (4.17.4: a smaller --stdev or a lower --md-temperature, or --max-energy); drop them, or add "
                 "test configs of that kind if such geometries matter.")
    else:
        outside = [""] * len(test_atoms)
    novelty = [((e.get("geometry") or {}).get("novelty")) for e in entries]

    # MD-check frames written by an earlier --md-check (their own set, not the test set)
    md_atoms, md_info, md_pending = [], [], []
    if state:
        for d in mlff_workflow.folders(state, wf_base, MD_CHECK_DIR):
            entry = mlff_workflow.config_entry(state, d, wf_base) or {}
            got, _, _ = collect_training_configs([d], f_out, []) if os.path.isdir(d) else ([], 0, False)
            if got:
                md_atoms.append(got[0])
                md_info.append(entry.get("md") or {})
            else:
                md_pending.append(d)
    if md_atoms or md_pending:
        runs = sorted({i.get("run", "?") for i in md_info} |
                      {(mlff_workflow.config_entry(state, d, wf_base) or {}).get("md", {}).get("run", "?")
                       for d in md_pending})
        print_dual(f"MD-check frames   : {len(md_atoms)} finished, {len(md_pending)} awaiting SIESTA "
                   f"({', '.join(runs)}) -- scored separately in [8]", f_out)

    try:
        with capture_library_noise(library_noise, "import mace"):
            from mace.calculators import MACECalculator
    except ImportError:
        fail("mace-torch not found. This tool needs the optional 'ml' extra: pip install stb_suite[ml]")

    ref_phonon = to_ang = None
    ph_info = {}
    if args.phonon_ref:
        if not os.path.isfile(os.path.join(args.phonon_ref, "phonopy_disp.yaml")):
            fail(f"'{args.phonon_ref}' has no phonopy_disp.yaml (expected a stb-phononsCreate folder).")
        try:
            with capture_library_noise(library_noise, "phonopy"):   # its own lines are reported in [6]
                ref_phonon, to_ang = load_siesta_phonons(args.phonon_ref, None, ph_info)
        except SystemExit:
            raise
        except Exception as e:
            fail(f"Could not load the SIESTA phonons from '{args.phonon_ref}': {e}")

    # ---- scoring (quiet: everything is printed below, table by table)
    scored = []
    calcs_by_round = {}
    for n, (k, models, entry) in enumerate(committees, 1):
        print_progress_line(f"Scoring {_round_label(k)} ({n}/{len(committees)}, {len(models)} model(s))...",
                            n, len(committees))
        t0 = time.monotonic()
        with capture_library_noise(library_noise, f"load {_round_label(k)}"):
            calcs = [MACECalculator(model_paths=[m], device=args.device, default_dtype="float64")
                     for m in models]
        with capture_library_noise(library_noise, f"score {_round_label(k)}"):
            res = score_committee(calcs, test_atoms, kinds)
        if res is None:
            warn(f"Scores: {_round_label(k)}: every model failed on every test config -- skipped.")
            continue
        res["elements"] = element_errors(test_atoms, res["_mean_preds"])
        if md_atoms:
            with capture_library_noise(library_noise, f"MD check {_round_label(k)}"):
                res["md"] = score_committee(calcs, md_atoms, [i.get("run", "?") for i in md_info])
        n_fail = sum(1 for mp in res["_member_preds"] for p in mp if p is None)
        if n_fail:
            warn(f"Scores: {_round_label(k)}: {n_fail} model evaluation(s) failed (model x config) -- left out.")
        if ref_phonon is not None:
            try:
                with capture_library_noise(library_noise, f"phonons {_round_label(k)}"):
                    res["phonons"] = phonon_errors(ref_phonon, committee_phonons(calcs, ref_phonon, to_ang),
                                                   args.phonon_mesh)
            except Exception as e:
                warn(f"Phonons: {_round_label(k)}: comparison failed ({e}).")
        res["seconds"] = time.monotonic() - t0
        seeds = {mlff_workflow.resolve(m["model"], wf_base): m.get("seed")
                 for m in (entry or {}).get("committee", []) if m.get("model")}
        res["seeds"] = [seeds.get(p) if seeds.get(p) is not None else lab
                        for p, lab in zip(models, _model_labels(models))]
        scored.append((k, res))
        calcs_by_round[k] = calcs
    finish_progress_line()
    if not scored:
        fail("No round could be scored.")
    if detail_k not in dict(scored):
        detail_k = scored[-1][0]
    detail = dict(scored)[detail_k]

    found = None
    if not args.skip_foundation_comparison:
        with capture_library_noise(library_noise, "foundation model"):
            fcalc = foundation_calc(size, args.device)
            fp = predict_all(fcalc, test_atoms)
        all_idx = list(range(len(test_atoms)))
        found = {"overall": metrics(fp, test_atoms, all_idx),
                 "groups": {g: metrics(fp, test_atoms, [i for i, kk in enumerate(kinds) if kk == g])
                            for g in sorted(set(kinds))},
                 "_preds": fp, "elements": element_errors(test_atoms, fp)}
        if md_atoms:
            with capture_library_noise(library_noise, "foundation model, MD check"):
                found["_md_preds"] = predict_all(fcalc, md_atoms)
        if found["overall"] is None:
            warn("Baseline: the foundation model failed on every test config.")
            found = None

    # ---- [2] models scored
    print_section("[2] MODELS SCORED", f_out)
    rows = []
    for k, models, entry in committees:
        if k not in dict(scored):
            continue
        e = entry or {}
        st = e.get("settings", {})
        oos = ((e.get("training_errors") or {}).get("out_of_sample") or {}).get("F_rmse_meV")
        cv = (e.get("training_errors") or {}).get("cv_folds")
        rows.append(([_round_label(k), str(len(models)), str(len(e.get("configs", [])) or "-"),
                      str(st.get("epochs", "-")), st.get("foundation_model", "-"),
                      e.get("init_from") or "foundation" if entry else "-",
                      f"{cv}-fold" if cv else ("split" if entry else "-"), _fmt(oos),
                      (e.get("created") or "-")[:16].replace("T", " ")], None))
    print_table(["Round", "Models", "Train cfgs", "Epochs", "Foundation", "Start", "Valid.", "Stage-2 OOS F",
                 "Trained"], rows, f_out)
    print_dual("Stage-2 OOS F = out-of-sample force RMSE on the TRAINING configs (meV/Ang), from 4.17.2.", f_out)
    found_sizes = {(e or {}).get("settings", {}).get("foundation_model") for _, _, e in committees if e}
    if len(found_sizes - {None}) > 1:
        warn(f"Setup: the scored rounds start from different foundation models ({', '.join(sorted(found_sizes - {None}))}) "
             "-- part of the change between rounds is the model size, not the data.")

    # ---- [3] test errors, every round
    print_section("[3] TEST ERRORS (committee mean, every round)", f_out)
    best_k = min(scored, key=lambda kv: kv[1]["overall"]["F_rmse_meV"])[0]
    f_rms_ref = max(1000 * float(np.sqrt(np.mean(np.concatenate(
        [a.arrays["REF_forces"].ravel() for a in test_atoms]) ** 2))), 1e-9)
    rows = []
    for k, v in scored:
        o = v["overall"]
        gain = (100 * (found["overall"]["F_rmse_meV"] - o["F_rmse_meV"]) / found["overall"]["F_rmse_meV"]
                if found else None)
        rows.append(([_round_label(k), _fmt(o["E_rmse_meV"], "{:.2f}"), _fmt(o["E_rmse_rel_meV"], "{:.2f}"),
                      f"{o['F_rmse_meV']:.1f}", _fmt(o.get("F_rmse_sigma_meV")), f"{o['F_mae_meV']:.1f}",
                      f"{o['F_maxerr_meV']:.0f}", f"{o['F_rel_pct']:.1f}", _fmt(o["S_rmse_GPa"], "{:.3f}"),
                      _fmt(gain, "{:+.0f} %")], 'green' if k == best_k and len(scored) > 1 else None))
    if found:
        fo = found["overall"]
        rows.append((["foundation", "(scale)", _fmt(fo["E_rmse_rel_meV"], "{:.2f}"), f"{fo['F_rmse_meV']:.1f}",
                      "-", f"{fo['F_mae_meV']:.1f}", f"{fo['F_maxerr_meV']:.0f}",
                      f"{100 * fo['F_rmse_meV'] / f_rms_ref:.1f}",
                      _fmt(fo["S_rmse_GPa"], "{:.3f}"), "-"], 'yellow'))
    print_table(["Model", "E abs", "E rel", "F RMSE", "+/-", "F MAE", "F max", "F rel %", "S RMSE",
                 "vs found."], rows, f_out)
    print_dual("Units: E meV/atom, F meV/Ang, S GPa. E abs includes a constant offset; E rel removes it (what "
               "matters for", f_out)
    print_dual("energy differences). +/- = bootstrap std over test configs. F rel = F RMSE / RMS SIESTA force. "
               "F max = largest", f_out)
    print_dual("component error. vs found. = force RMSE gain over the raw foundation model."
               + (f" Best round: {best_k} (green)." if len(scored) > 1 else ""), f_out)

    # ---- [4] detailed round
    print_section(f"[4] {_round_label(detail_k).upper()} IN DETAIL", f_out)
    print_dual(color_text("Per kind:", 'bold'), f_out)
    rows = []
    for g, v in detail["groups"].items():
        if v is None:
            continue
        fg = (found or {}).get("groups", {}).get(g) if found else None
        gain = (100 * (fg["F_rmse_meV"] - v["F_rmse_meV"]) / fg["F_rmse_meV"]
                if fg and fg["F_rmse_meV"] >= _MIN_SPREAD_MEV else None)   # no forces to compare (equilibrium)
        rows.append(([g, str(v["n"]), f"{v['F_rmse_meV']:.1f}", f"{v['F_maxerr_meV']:.0f}",
                      _fmt(fg["F_rmse_meV"] if fg else None), _fmt(gain, "{:+.0f} %"),
                      _fmt(v["E_rmse_rel_meV"], "{:.2f}"), _fmt(v["S_rmse_GPa"], "{:.3f}")], None))
    print_table(["Kind", "N", "F RMSE", "F max", "Found. F", "Gain", "E rel", "S RMSE"], rows, f_out)
    worst_g = max(((g, v) for g, v in detail["groups"].items() if v), key=lambda gv: gv[1]["F_rmse_meV"])
    if len(detail["groups"]) > 1:
        print_dual(f"Hardest kind      : {worst_g[0]} (F {worst_g[1]['F_rmse_meV']:.1f} meV/Ang)", f_out)

    print_dual("\n" + color_text("Per element (forces on the atoms of each element):", 'bold'), f_out)
    rows = []
    fel = (found or {}).get("elements") or {}
    for el, v in detail["elements"].items():
        fe = fel.get(el)
        gain = (100 * (fe["F_rmse_meV"] - v["F_rmse_meV"]) / fe["F_rmse_meV"]
                if fe and fe["F_rmse_meV"] >= _MIN_SPREAD_MEV else None)
        rows.append(([el, str(v["n_atoms"]), f"{v['F_rmse_meV']:.1f}", f"{v['F_mae_meV']:.1f}",
                      f"{v['F_maxatom_meV']:.0f}", f"{v['F_ref_rms_meV']:.0f}",
                      _fmt(100 * v["F_rmse_meV"] / v["F_ref_rms_meV"] if v["F_ref_rms_meV"] >= _MIN_SPREAD_MEV else None),
                      _fmt(fe["F_rmse_meV"] if fe else None), _fmt(gain, "{:+.0f} %")], None))
    print_table(["Element", "Atoms", "F RMSE", "F MAE", "Worst atom", "Ref RMS", "F rel %", "Found. F", "Gain"],
                rows, f_out)
    if len(detail["elements"]) > 1:
        els = detail["elements"]
        worst_el = max(els, key=lambda e: els[e]["F_rmse_meV"])
        best_el = min(els, key=lambda e: els[e]["F_rmse_meV"])
        ratio = els[worst_el]["F_rmse_meV"] / max(els[best_el]["F_rmse_meV"], 1e-9)
        print_dual(f"Hardest element   : {worst_el} ({ratio:.1f}x the error on {best_el})", f_out)
        if ratio > 2:
            warn(f"Elements: forces on {worst_el} are predicted {ratio:.1f}x worse than on {best_el} -- if that is "
                 f"not just larger forces (compare 'F rel %' in [4]), add configs that move {worst_el} "
                 "(its environments are under-sampled).")

    print_dual("\n" + color_text("Per config (worst first):", 'bold'), f_out)
    order = sorted((i for i, r in enumerate(detail["per_config"]) if r is not None),
                   key=lambda i: -detail["per_config"][i]["F_rmse_meV"])
    fp_ = found["_preds"] if found else None
    has_nov = any(x is not None for x in novelty)
    rows = []
    show_n = order if len(order) <= 30 else order[:20]
    for rank, i in enumerate(show_n):
        r = detail["per_config"][i]
        a = test_atoms[i]
        ratio = (r["err_max_meV"] / r["spread_max_meV"]
                 if r["spread_max_meV"] is not None and r["spread_max_meV"] >= _MIN_SPREAD_MEV else None)
        cells = [os.path.basename(a.info["folder"]), kinds[i], str(len(a)), f"{ds_test['fmax'][i]:.2f}",
                 f"{r['F_rmse_meV']:.1f}", _fmt(_f_rmse(fp_[i], a) if fp_ and fp_[i] else None),
                 f"{a.get_chemical_symbols()[r['worst_atom']]}{r['worst_atom']} {r['err_max_meV']:.0f}",
                 f"{r['dE_meV']:+.2f}", _fmt(r["spread_max_meV"]), _fmt(ratio)]
        if has_nov:
            cells.append(_fmt(novelty[i], "{:.2f}"))
        cells.append(outside[i] or "in")
        rows.append((cells, 'yellow' if rank < 3 and len(order) > 3 else None))
    heads = ["Config", "Kind", "Atoms", "Fmax ref", "F RMSE", "Found. F", "Worst atom", "dE", "Spread",
             "Err/spr"] + (["Novelty"] if has_nov else []) + ["Range"]
    print_table(heads, rows, f_out)
    if len(order) > len(show_n):
        print_dual(f"... {len(order) - len(show_n)} more (all in --save-data).", f_out)
    print_dual("Fmax ref eV/Ang; F RMSE, worst atom |dF|, spread (max over atoms) meV/Ang; dE meV/atom with the "
               "round's offset removed;", f_out)
    print_dual("Err/spr = worst-atom error / max spread; Range: in, or outside the training range in V/F/E"
               + ("; novelty from Stage 1 (in typical spacings)." if has_nov else "."), f_out)

    print_dual("\n" + color_text("Analysis:", 'bold'), f_out)
    cfg_err = np.array([detail["per_config"][i]["F_rmse_meV"] if detail["per_config"][i] else np.nan
                        for i in range(len(test_atoms))])
    c_f = _corr(ds_test["fmax"], cfg_err)
    if c_f is not None:
        print_dual(f"Error vs force    : corr(config F RMSE, SIESTA max|F|) = {c_f:+.2f}"
                   + (" -- the error grows with how far the config is from equilibrium" if c_f > 0.5 else ""), f_out)
    if has_nov:
        c_n = _corr([x if x is not None else np.nan for x in novelty], cfg_err)
        if c_n is not None:
            print_dual(f"Error vs novelty  : corr(config F RMSE, Stage-1 novelty) = {c_n:+.2f}"
                       + (" -- configs unlike the training set are predicted worse" if c_n > 0.5 else ""), f_out)
    if any(outside):
        ins = [cfg_err[i] for i, o in enumerate(outside) if not o and np.isfinite(cfg_err[i])]
        outs = [cfg_err[i] for i, o in enumerate(outside) if o and np.isfinite(cfg_err[i])]
        if ins and outs:
            print_dual(f"Inside vs outside : mean config F RMSE {np.mean(ins):.1f} inside the training range, "
                       f"{np.mean(outs):.1f} outside", f_out)
    oos2 = (((detail_entry or {}).get("training_errors") or {}).get("out_of_sample") or {}).get("F_rmse_meV")
    if oos2:
        ratio = detail["overall"]["F_rmse_meV"] / oos2
        print_dual(f"Generalization    : test F {detail['overall']['F_rmse_meV']:.1f} vs Stage-2 out-of-sample "
                   f"{oos2:.1f} meV/Ang ({ratio:.2f}x)", f_out)
        if ratio > _GENERALIZATION_WARN_RATIO:
            warn(f"Generalization: the test error is {ratio:.1f}x the Stage-2 out-of-sample error -- the test set "
                 "probes configurations the training set does not cover (see Range/Novelty in [4]); sample "
                 "them in training (4.17.4, or stb-mlff with the test set's settings).")
        elif ratio < 0.5:
            print_dual("                    (the test set is easier than the training configs)", f_out)
    o = detail["overall"]
    if o["E_rmse_meV"] is not None and o["E_rmse_rel_meV"] is not None:
        print_dual(f"Energy offset     : {o['E_offset_meV']:+.2f} meV/atom (constant; E rel {o['E_rmse_rel_meV']:.2f} "
                   "is the error on energy differences)", f_out)

    # ---- [5] committee
    if len(detail["members"]) > 1:
        print_section(f"[5] COMMITTEE ({_round_label(detail_k)})", f_out)
        fr = [m["F_rmse_meV"] for m in detail["members"] if m]
        med = float(np.median(fr))
        rows = []
        for seed, m in zip(detail["seeds"], detail["members"]):
            if not m:
                rows.append(([str(seed), "-", "-", "-", "-", "failed"], 'red'))
                continue
            flag = m["F_rmse_meV"] > _OUTLIER_MEMBER_RATIO * med
            rows.append(([str(seed), _fmt(m["E_rmse_rel_meV"], "{:.2f}"), f"{m['F_rmse_meV']:.1f}",
                          f"{m['F_mae_meV']:.1f}", _fmt(m["S_rmse_GPa"], "{:.3f}"),
                          "outlier" if flag else ("best" if m["F_rmse_meV"] == min(fr) else "")],
                         'yellow' if flag else None))
            if flag:
                warn(f"Committee: member {seed} of {_round_label(detail_k)} has F {m['F_rmse_meV']:.1f} meV/Ang, "
                     f"{m['F_rmse_meV'] / med:.1f}x the median member -- it drags the committee mean; check its "
                     "training (4.17.2 table [3]) or retrain it.")
        rows.append((["mean", _fmt(o["E_rmse_rel_meV"], "{:.2f}"), f"{o['F_rmse_meV']:.1f}", f"{o['F_mae_meV']:.1f}",
                      _fmt(o["S_rmse_GPa"], "{:.3f}"), "committee"], 'green'))
        print_table(["Seed", "E rel", "F RMSE", "F MAE", "S RMSE", ""], rows, f_out)
        gain = 100 * (np.mean(fr) - o["F_rmse_meV"]) / np.mean(fr)
        print_dual(f"Averaging         : committee mean {o['F_rmse_meV']:.1f} vs mean member {np.mean(fr):.1f}, "
                   f"best member {min(fr):.1f} meV/Ang ({gain:+.0f} % from averaging)", f_out)
        if o["F_rmse_meV"] > min(fr) * 1.05:
            print_dual("                    (the best single member beats the mean: the members are not "
                       "equally good)", f_out)
        cal = detail["calibration"]
        if cal:
            print_dual(f"Spread            : mean {cal['spread_mean_meV']:.1f}, max {cal['spread_max_meV']:.1f} meV/Ang "
                       "(max over atoms, per config)", f_out)
            if cal["error_over_spread"] is not None:
                print_dual(f"True error/spread : {cal['error_over_spread']:.1f}x (median over configs with spread >= "
                           f"{_MIN_SPREAD_MEV:g} meV/Ang) -- an error bar from the spread is ~"
                           f"{cal['error_over_spread']:.1f} x spread", f_out)
            cc = (f"{cal['config_corr']:+.2f}" if cal["config_corr"] is not None
                  else f"n/a (needs >= {_MIN_CONFIGS_CORR} configs)")
            ca = (f"{cal['atom_corr']:+.2f}" if cal["atom_corr"] is not None
                  else f"n/a (needs >= {_MIN_ATOMS_CORR} atoms)")
            print_dual(f"Correlation       : per config {cc}; per atom {ca}", f_out)
            print_dual("                    (per config: does the spread pick the worst configs -- what 4.17.4 "
                       "needs; per atom: the worst atom)", f_out)
            if (cal["error_over_spread"] or 0) > _CALIBRATION_WARN_RATIO:
                warn(f"Committee: the spread underestimates the error ~{cal['error_over_spread']:.0f}x on the test set "
                     "-- use it to RANK active-learning candidates, not as an error bar (members fine-tuned from "
                     "the same foundation model are nearly the same model).")
            if cal["config_corr"] is not None and cal["config_corr"] < 0.3:
                warn(f"Committee: its spread barely tracks the test error per config (corr {cal['config_corr']:+.2f}) "
                     "-- active learning by spread will pick close to at random here.")
    else:
        print_section("[5] COMMITTEE", f_out)
        print_dual(f"{_round_label(detail_k)} is a single model: no spread, no calibration. Train a committee "
                   "(4.17.2 --committee 3 or --cv 5) for both.", f_out)

    # ---- [6] phonons
    if args.phonon_ref:
        print_section("[6] PHONONS vs SIESTA", f_out)
        st = ph_info.get("residual_status")
        if st == "subtracted":
            print_dual(f"Residual forces  : subtracted (undisplaced supercell, max {ph_info['residual_max']:.4f} eV/Ang "
                       "per atom) -- SIESTA's and each model's own", f_out)
        else:
            print_dual(color_text(f"Residual forces  : NOT subtracted from the SIESTA forces ({st})", 'yellow'), f_out)
            warn(f"Phonons: {args.phonon_ref}/perfect_supercell/ " + ("has no finished SIESTA run" if st == "not finished"
                 else "is missing" if st == "absent" else "is unreadable") + " -- the reference's residual forces stay in "
                 "the SIESTA force constants (each model's are subtracted), which biases Gamma-point modes; run SIESTA "
                 "in the undisplaced supercell (stb-phononsCreate writes it).")
        rows = []
        for k, v in scored:
            ph = v.get("phonons")
            if not ph:
                continue
            rows.append(([_round_label(k), f"{ph['rmse_THz']:.2f}", f"{ph['acoustic_rmse_THz']:.2f}",
                          _fmt(ph["optical_rmse_THz"], "{:.2f}"), f"{ph['mean_shift_THz']:+.2f}",
                          f"{ph['max_abs_THz']:.2f}", f"{ph['model_max_THz']:.2f} / {ph['ref_max_THz']:.2f}",
                          f"{ph['model_min_THz']:+.2f} / {ph['ref_min_THz']:+.2f}", str(ph["model_imaginary_q"])],
                         None))
            if ph["model_imaginary_q"] > ph["ref_imaginary_q"]:
                warn(f"Phonons: {_round_label(k)} has imaginary modes at {ph['model_imaginary_q']} of {ph['n_q']} "
                     f"q-points (min {ph['model_min_THz']:.2f} THz) that SIESTA does not have -- the model is "
                     "dynamically unstable there; add near-equilibrium configs (small --stdev) to training.")
        print_table(["Round", "RMSE", "Acoustic", "Optical", "Shift", "Max |dev|", "Top (ML/ref)",
                     "Min (ML/ref)", "Imag. q"], rows, f_out)
        print_dual(f"THz on a {'x'.join(map(str, args.phonon_mesh))} q-mesh, the SAME displacements as SIESTA; "
                   "shift = mean ML - ref (systematic softening < 0).", f_out)
    else:
        print_section("[6] PHONONS", f_out)
        print_dual("Not compared (no --phonon-ref). A finished stb-phononsCreate folder of the same system gives a "
                   "q-point by q-point check.", f_out)

    # ---- [7] equation of state
    print_section("[7] EQUATION OF STATE (test configs spanning a volume range)", f_out)
    eos_plot, eos_store = None, {}
    ek = eos_kinds(test_atoms, kinds)
    for g, ii in ek.items():
        vol = np.array([test_atoms[i].get_volume() / len(test_atoms[i]) for i in ii])
        fmx = max(ds_test["fmax"][i] for i in ii)
        print_dual(color_text(f"Kind '{g}': {len(ii)} configs, V/atom {vol.min():.3f}-{vol.max():.3f} Ang^3 "
                              f"({100 * (vol.max() - vol.min()) / vol.mean():.1f} %), max|F| up to {fmx:.2f} eV/Ang",
                              'bold'), f_out)
        sf = eos_fit_set(test_atoms, ii, [test_atoms[i].info["REF_energy"] / len(test_atoms[i]) for i in ii])
        if sf is None:
            warn(f"EOS: the SIESTA E(V) of the '{g}' test configs could not be fitted -- skipped.")
            continue
        fits = {"siesta": sf}
        per = [(_round_label(k), k, eos_fit_set(test_atoms, ii, [v["_mean_preds"][i]["e"] if v["_mean_preds"][i]
                                                                 else None for i in ii])) for k, v in scored]
        if found:
            fits["foundation"] = eos_fit_set(test_atoms, ii, [found["_preds"][i]["e"] if found["_preds"][i]
                                                               else None for i in ii])
            per.append(("foundation", None, fits["foundation"]))
        fits["detail"] = dict((k, f) for _, k, f in per if k is not None).get(detail_k)
        rows = [(["SIESTA", f"{sf['v0']:.3f}", "-", f"{sf['b0_gpa']:.1f}", "-", _fmt(sf["bprime"], "{:.2f}"),
                  f"{sf['r2']:.4f}"], 'blue')]
        for lab, k, fit in per:
            if fit is None:
                rows.append(([lab, "-", "-", "-", "-", "-", "fit failed"], None))
                continue
            dv = 100 * (fit["v0"] - sf["v0"]) / sf["v0"]
            db = 100 * (fit["b0_gpa"] - sf["b0_gpa"]) / sf["b0_gpa"]
            rows.append(([lab, f"{fit['v0']:.3f}", f"{dv:+.2f}", f"{fit['b0_gpa']:.1f}", f"{db:+.1f}",
                          _fmt(fit["bprime"], "{:.2f}"), f"{fit['r2']:.4f}"],
                         'green' if k == detail_k else ('yellow' if k is None else None)))
            if k is not None:
                eos_store.setdefault(k, {})[g] = {"v0": round(fit["v0"], 4), "b0_gpa": round(fit["b0_gpa"], 2),
                                                 "bprime": None if fit["bprime"] is None else round(fit["bprime"], 3),
                                                 "dV0_pct": round(dv, 3), "dB0_pct": round(db, 2),
                                                 "siesta_v0": round(sf["v0"], 4), "siesta_b0_gpa": round(sf["b0_gpa"], 2)}
            if k == detail_k and (abs(db) > 10 or abs(dv) > 1):
                warn(f"EOS: {_round_label(k)} gives B0 {fit['b0_gpa']:.1f} GPa ({db:+.0f} %) and V0 {dv:+.2f} % against "
                     f"SIESTA on the '{g}' test configs -- the volume response is off; train with stress "
                     "(4.17.2 default) on volume-strained configs (stb-mlff --volume-strains).")
        print_table(["Model", "V0 (Ang^3/at)", "dV0 %", "B0 (GPa)", "dB0 %", "B0'", "R2"], rows, f_out)
        if sf["r2"] < 0.99:
            print_dual(f"(SIESTA R2 {sf['r2']:.3f}: these configs are not pure strains -- the rattle adds energy; "
                       "compare models with each other, the fit is equally noisy for all.)", f_out)
        if sf["outside"]:
            warn(f"EOS: the SIESTA V0 of the '{g}' test configs lies outside their volume range -- B0 is an extrapolation.")
        if eos_plot is None:
            eos_plot = (g, fits)
    if not ek:
        print_dual("No test kind spans a volume range (needs >= 5 distinct volumes over >= 1 %: the fit has 4 "
                   "parameters). Stage 1 cycles the test kinds one config each, so give it 5+ --volume-strains "
                   "values and --n-test >= (number of --stdev values + number of --volume-strains).", f_out)
    if ek:
        print_dual("Birch-Murnaghan fit of E/atom vs V/atom, the same configs for every model; dV0/dB0 vs SIESTA.", f_out)
        from stb.core.kspace import detect_vacuum_axes
        a0 = test_atoms[next(iter(ek.values()))[0]]
        if any(detect_vacuum_axes(a0.get_scaled_positions(), np.asarray(a0.cell), 5.0)):
            print_dual("(Not a bulk cell: V includes the vacuum, so B0 is not the material's bulk modulus -- read the "
                       "table as model vs SIESTA on the same E(V) curve only.)", f_out)

    # ---- [8] MD check
    print_section("[8] MD CHECK (error along a trajectory run with the model)", f_out)
    md_plot, md_store = None, {}
    scored_d = dict(scored)
    if md_atoms:
        md_plot = {}
        test_f = detail["overall"]["F_rmse_meV"]
        fpm = (found or {}).get("_md_preds")
        for run in sorted({i.get("run", "?") for i in md_info}):
            ii = sorted((j for j, inf in enumerate(md_info) if inf.get("run", "?") == run),
                        key=lambda j: md_info[j].get("time_fs", 0))
            inf0 = md_info[ii[0]]
            drv = inf0.get("round")
            print_dual(color_text(f"Run {run}: NVT {inf0.get('temperature', '?')} K with {_round_label(drv)}'s model "
                                  f"({inf0.get('model', '?')}), {len(ii)} frame(s) of {len(md_atoms[ii[0]])} atoms",
                                  'bold'), f_out)
            rows, pts = [], []
            for j in ii:
                a, inf = md_atoms[j], md_info[j]
                rd = scored_d[detail_k]["md"]["per_config"][j] if scored_d[detail_k].get("md") else None
                rv = (scored_d[drv]["md"]["per_config"][j]
                      if drv in scored_d and drv != detail_k and scored_d[drv].get("md") else None)
                ff = _f_rmse(fpm[j], a) if fpm and fpm[j] else None
                rows.append(([str(inf.get("time_fs", "?")), _fmt(inf.get("t_inst"), "{:.0f}"),
                              f"{np.linalg.norm(a.arrays['REF_forces'], axis=1).max():.2f}",
                              _fmt(rd["F_rmse_meV"] if rd else None)]
                             + ([_fmt(rv["F_rmse_meV"] if rv else None)] if drv in scored_d and drv != detail_k else [])
                             + [_fmt(ff), _fmt(rd["dE_meV"] if rd else None, "{:+.2f}"),
                                _fmt(rd["spread_max_meV"] if rd else None)], None))
                if rd:
                    pts.append((inf.get("time_fs", 0), rd["F_rmse_meV"]))
            heads = (["t (fs)", "T (K)", "Fmax ref", f"F {_round_label(detail_k)}"]
                     + ([f"F {_round_label(drv)} (ran MD)"] if drv in scored_d and drv != detail_k else [])
                     + ["Found. F", "dE", "Spread"])
            print_table(heads, rows, f_out)
            md_plot[run] = pts
            for k, v in scored:
                if v.get("md"):
                    md_store.setdefault(k, {})[run] = [
                        [md_info[j].get("time_fs"), round(v["md"]["per_config"][j]["F_rmse_meV"], 2)]
                        for j in ii if v["md"]["per_config"][j]]
            if len(pts) >= 3:
                n3 = max(1, len(pts) // 3)
                early = float(np.mean([p[1] for p in pts[:n3]]))
                late = float(np.mean([p[1] for p in pts[-n3:]]))
                c = _corr([p[0] for p in pts], [p[1] for p in pts])
                print_dual(f"Trend             : first third {early:.1f} -> last third {late:.1f} meV/Ang "
                           f"({late / max(early, 1e-9):.1f}x), corr(t, F) {_fmt(c, '{:+.2f}')}; test set "
                           f"{test_f:.1f}", f_out)
                if late > _MD_GROWTH_WARN * early:
                    warn(f"MD check: along run {run} the error grows {late / early:.1f}x -- the trajectory drifts into "
                         "configurations the model has not learned; add these frames to training "
                         "(--md-check-to-train) and retrain.")
            if pts and np.mean([p[1] for p in pts]) > 2 * test_f:
                warn(f"MD check: run {run} frames are predicted {np.mean([p[1] for p in pts]) / test_f:.1f}x worse than "
                     "the test set -- the MD visits configurations the test set does not contain, so the test score "
                     "is optimistic for this use.")
        if args.md_check_to_train:
            done = {mlff_workflow.rel(a.info["folder"], wf_base) for a in md_atoms}
            moved = 0
            for c in state["configs"]:
                if c.get("set") == MD_CHECK_DIR and c["folder"] in done:
                    c["set"], c["kind"] = "train", "md_check"
                    moved += 1
            print_dual(color_text(f"[OK] {moved} finished MD-check frame(s) moved to the training set -- the next "
                                  "stb-mlffAnalysis (4.17.2) trains on them.", 'green'), f_out)
    if md_pending:
        print_dual(f"{len(md_pending)} frame(s) still need SIESTA:", f_out)
        for run_dir in sorted({os.path.dirname(d) for d in md_pending}):
            first = next(d for d in md_pending if os.path.dirname(d) == run_dir)
            print_dual(f"  {_siesta_loop(run_dir, first)}", f_out)
            n_here = sum(1 for d in md_pending if os.path.dirname(d) == run_dir)
            warn(f"SIESTA: run it in {os.path.relpath(run_dir)}/ ({n_here} MD-check frame(s)) -- they are scored "
                 "on the next stb-mlffValidate.")
    if args.md_check:
        dims = list(args.md_check_supercell)
        models_b = dict((k, m) for k, m, _ in committees)[best_k]
        run = (f"r{best_k:02d}_{args.md_check_temperature:g}K"
               + ("" if dims == [1, 1, 1] else "_" + "x".join(map(str, dims))))
        run_dir = os.path.join(wf_base, MD_CHECK_DIR, run)
        if any((c.get("md") or {}).get("run") == run for c in state["configs"]):
            print_dual(f"--md-check: run {run} is already in the workflow file -- not rewritten (delete its entries "
                       "and folder to redo it, or change --md-check-temperature).", f_out)
        else:
            ref = None
            try:
                ref = md_check_reference(state, wf_base, dims)
            except Exception as e:
                warn(f"MD check: not written ({e}).")
            if ref:
                t0 = time.monotonic()
                try:
                    with capture_library_noise(library_noise, "MD check"):
                        frames, stopped = run_with_spinner(
                            run_md_check, ref["atoms"], models_b[0], args.device, args.md_check_temperature,
                            args.md_check_steps, args.md_check_frames, 7, label=f"MD check ({run})")
                except Exception as e:
                    frames, stopped = [], str(e)
                new_entries = []
                for n, (step, fr, t_inst) in enumerate(frames):
                    pos = fr.get_scaled_positions() if ref["fractional"] else fr.get_positions()
                    folder = write_config_folder(run_dir, n, ref["template"], ref["basename"], pos, None,
                                                 ref["pseudos"], extra_block=ref["extra"], prefix="frame",
                                                 calc_template=ref["calc"])
                    new_entries.append({"folder": folder, "set": MD_CHECK_DIR, "kind": "md_check",
                                        "origin": f"MD {args.md_check_temperature:g} K t={step} fs, round {best_k} model",
                                        "stage": "validate", "round": best_k, "n_atoms": len(fr),
                                        "md": {"run": run, "time_fs": step, "temperature": args.md_check_temperature,
                                               "t_inst": round(t_inst, 1), "round": best_k,
                                               "model": mlff_workflow.rel(models_b[0], wf_base), "supercell": dims}})
                if ref["tmp"]:
                    shutil.rmtree(ref["tmp"], ignore_errors=True)
                mlff_workflow.add_configs(state, new_entries, wf_base)
                if stopped:
                    warn(f"MD check: the MD with {_round_label(best_k)}'s model stopped early ({stopped}) -- it is "
                         f"unstable at {args.md_check_temperature:g} K; the {len(frames)} frame(s) before that were written.")
                if new_entries:
                    print_dual(color_text(f"[OK] --md-check: {len(new_entries)} frame(s) of a {args.md_check_steps} fs NVT "
                                          f"run at {args.md_check_temperature:g} K with {_round_label(best_k)}'s model "
                                          f"({len(ref['atoms'])} atoms, {time.monotonic() - t0:.0f} s) written to "
                                          f"{os.path.relpath(run_dir)}/", 'green'), f_out)
                    print_dual("  Run SIESTA in them, then run stb-mlffValidate again to see the error along the run:", f_out)
                    print_dual(f"  {_siesta_loop(run_dir, new_entries[0]['folder'])}", f_out)
                    warn(f"SIESTA: run it in {os.path.relpath(run_dir)}/ ({len(new_entries)} MD-check frame(s)) -- "
                         "they are scored on the next stb-mlffValidate.")
    if not md_atoms and not md_pending and not args.md_check:
        print_dual("Not run. --md-check runs an MD with the best model and writes frames for SIESTA; the next "
                   "validation then shows whether the error grows along the trajectory.", f_out)

    # ---- [9] round comparison
    print_section("[9] ROUND COMPARISON", f_out)
    cmp_plot, cmp_store = None, None
    if args.compare:
        ka, kb = args.compare
        ra, rb = scored_d[ka], scored_d[kb]
        rows_c, better, worse, same = compare_rounds(ra, rb, test_atoms, kinds)
        oa, ob = ra["overall"], rb["overall"]
        print_dual(f"Overall           : F {oa['F_rmse_meV']:.1f} -> {ob['F_rmse_meV']:.1f} meV/Ang "
                   f"({100 * (ob['F_rmse_meV'] - oa['F_rmse_meV']) / oa['F_rmse_meV']:+.1f} %), E rel "
                   f"{_fmt(oa['E_rmse_rel_meV'], '{:.2f}')} -> {_fmt(ob['E_rmse_rel_meV'], '{:.2f}')} meV/atom, S "
                   f"{_fmt(oa['S_rmse_GPa'], '{:.3f}')} -> {_fmt(ob['S_rmse_GPa'], '{:.3f}')} GPa", f_out)
        sig = ob.get("F_rmse_sigma_meV")
        if sig is not None and abs(ob["F_rmse_meV"] - oa["F_rmse_meV"]) < sig:
            print_dual(f"                    (within the test set's noise, +/- {sig:.1f})", f_out)
        print_dual(f"Configs           : {better} better, {worse} worse, {same} within +/-{_COMPARE_BAND_PCT:g} % "
                   f"(of {len(rows_c)})", f_out)
        rows = []
        for g in sorted(set(kinds)):
            va, vb = ra["groups"].get(g), rb["groups"].get(g)
            if not va or not vb:
                continue
            sel = [r for r in rows_c if kinds[r[0]] == g and r[3] is not None]
            ch = 100 * (vb["F_rmse_meV"] - va["F_rmse_meV"]) / va["F_rmse_meV"] if va["F_rmse_meV"] >= _MIN_SPREAD_MEV else None
            rows.append(([g, str(va["n"]), f"{va['F_rmse_meV']:.1f}", f"{vb['F_rmse_meV']:.1f}", _fmt(ch, "{:+.1f} %"),
                          str(sum(r[3] < -_COMPARE_BAND_PCT for r in sel)), str(sum(r[3] > _COMPARE_BAND_PCT for r in sel))],
                         'green' if ch is not None and ch < -_COMPARE_BAND_PCT else
                         ('yellow' if ch is not None and ch > _COMPARE_BAND_PCT else None)))
        print_table(["Kind", "N", f"F round {ka}", f"F round {kb}", "Change", "Better", "Worse"], rows, f_out)
        ranked = sorted(rows_c, key=lambda r: -(r[3] if r[3] is not None else -1e9))
        shown = ranked if len(ranked) <= 20 else ranked[:10] + ranked[-10:]
        rows = [([os.path.basename(test_atoms[i].info["folder"]), kinds[i], f"{fa:.1f}", f"{fb:.1f}",
                  _fmt(ch, "{:+.1f} %")],
                 'yellow' if ch is not None and ch > _COMPARE_BAND_PCT else
                 ('green' if ch is not None and ch < -_COMPARE_BAND_PCT else None)) for i, fa, fb, ch in shown]
        print_dual("\n" + color_text("Per config (largest regression first):", 'bold'), f_out)
        print_table(["Config", "Kind", f"F round {ka}", f"F round {kb}", "Change"], rows, f_out)
        if len(ranked) > len(shown):
            print_dual(f"... {len(ranked) - len(shown)} in between not shown.", f_out)
        if worse > better:
            warn(f"Compare: round {kb} predicts more test configs worse ({worse}) than better ({better}) than round {ka}.")
        cmp_plot = (ka, kb, rows_c)
        cmp_store = {"rounds": [ka, kb], "better": better, "worse": worse, "unchanged": same}
    else:
        print_dual("Not requested. --compare A B shows, config by config and per kind, what got better or worse "
                   "between two rounds.", f_out)

    # ---- [10] history and verdict
    print_section("[10] ROUND HISTORY, LEARNING CURVE & VERDICT", f_out)
    rounds_only = [(k, v) for k, v in scored if k != "models"]
    code, lc = None, None
    best_info = None
    if rounds_only:
        rows = []
        best_f = None
        prev = None
        for k, v in rounds_only:
            f = v["overall"]["F_rmse_meV"]
            entry = dict((kk, e) for kk, _, e in committees)[k] or {}
            d_prev = _fmt(100 * (f - prev) / prev if prev else None, "{:+.1f} %")
            best_f = f if best_f is None else min(best_f, f)
            rows.append(([str(k), _fmt(_n_train(entry), "{}"), f"{f:.1f}",
                          _fmt(v["overall"].get("F_rmse_sigma_meV")), d_prev,
                          _fmt(v["overall"]["E_rmse_rel_meV"], "{:.2f}"), _fmt(v["overall"]["S_rmse_GPa"], "{:.3f}"),
                          _fmt((v.get("phonons") or {}).get("rmse_THz"), "{:.2f}"),
                          _fmt((v.get("calibration") or {}).get("spread_mean_meV"))],
                         'green' if k == best_k and len(rounds_only) > 1 else None))
            prev = f
        print_table(["Round", "Train cfgs", "F RMSE", "+/-", "vs prev", "E rel", "S RMSE", "Phonons",
                     "Spread"], rows, f_out)
        code, message = verdict([(k, v["overall"]["F_rmse_meV"]) for k, v in rounds_only],
                                args.target_force, args.plateau,
                                rounds_only[-1][1]["overall"].get("F_rmse_sigma_meV"))
        color = {'stop': 'green', 'regressed': 'red'}.get(code, 'yellow')
        print_dual(color_text(f"[{code.upper()}] {message}", color), f_out)
        if args.target_force is None:
            print_dual("(No --target-force given: the verdict can only say continue/plateau/regressed.)", f_out)
        entries_by_k = dict((kk, e) for kk, _, e in committees)
        lc = learning_curve([(_n_train(entries_by_k[k]), v["overall"]["F_rmse_meV"]) for k, v in rounds_only],
                            args.target_force)
        if "b" in lc:
            n_now = max(n for n, _ in lc["points"])
            print_dual(f"Learning curve    : F = {_power_law(lc)} over {len(lc['points'])} training-set "
                       f"sizes (R2 {lc['r2']:.2f}; best round per N, assumes the same recipe)", f_out)
            if lc["b"] < -0.02:
                print_dual("                    the test error GROWS with the training-set size -- the added configs "
                           "pull the model away from the test region", f_out)
                warn("Learning curve: the test error grows as configs are added -- the newer configs are unlike the "
                     "test set (see 'Far from the test' in [1]: usually an active-learning round sampled too hard). "
                     "Go back to the best round and sample closer to the region of use (4.17.4: smaller --stdev, "
                     "--md-temperature near the use temperature, --supercell for longer wavelengths).")
            elif lc["b"] <= 0.02:
                print_dual("                    more configs have not lowered the test error -- the sampling, not the "
                           "amount of data, limits the model", f_out)
                warn("Learning curve: the test error does not fall with the training-set size -- add DIFFERENT "
                     "configurations (4.17.4 --supercell / higher --md-temperature), not more of the same.")
            elif lc["n_target"]:
                need = lc["n_target"]
                print_dual(f"                    ~{need:.0f} training configs for the {args.target_force:g} meV/Ang target "
                           f"(now {n_now}" + (", already there)" if need <= n_now else f", ~{need - n_now:.0f} more)"),
                           f_out)
            if lc["r2"] < 0.5:
                print_dual("                    (poor fit: the rounds differ in more than their size -- read it as a "
                           "trend only)", f_out)
        else:
            print_dual(f"Learning curve    : needs rounds with 3+ different training-set sizes (have "
                       f"{len(lc['points'])}: {', '.join(str(n) for n, _ in lc['points'])} configs)", f_out)
    else:
        print_dual("(--model committee: no round history to compare.)", f_out)
        if args.target_force is not None:
            f = detail["overall"]["F_rmse_meV"]
            met = f <= args.target_force
            print_dual(color_text(f"[{'TARGET MET' if met else 'TARGET NOT MET'}] test F RMSE {f:.1f} "
                                  f"{'<=' if met else '>'} {args.target_force:g} meV/Ang.",
                                  'green' if met else 'yellow'), f_out)
    best_models = dict((k, m) for k, m, _ in committees).get(best_k, [])

    # ---- outputs
    png = near_wf(args.png)
    show_requested = args.show
    if args.show and not can_show():
        warn("Plots: --show needs a display (no interactive matplotlib backend here) -- the figure is only saved.")
        args.show = False
    try:
        plot_validation(scored, detail_k, found, test_atoms, kinds, args.target_force, png, args.show,
                        extras={"learning": lc if lc and len(lc["points"]) >= 2 else None, "eos": eos_plot,
                                "md": md_plot if md_plot and any(md_plot.values()) else None,
                                "compare": cmp_plot})
    except Exception as e:
        warn(f"Plots: summary figure not written ({e}).")
        png = None

    data_path = None
    if args.save_data:
        data_path = near_wf(DATA_FILE)
        with open(data_path, "w") as f:
            f.write(f"# Per test config errors of {_round_label(detail_k)}, stb-mlffValidate\n"
                    "# folder kind n_atoms E_ref(eV/atom) Fmax_ref(eV/Ang) F_rmse(meV/Ang) "
                    "F_rmse_found(meV/Ang) dE(meV/atom) err_max(meV/Ang) spread_max(meV/Ang) novelty range\n")
            for i, a in enumerate(test_atoms):
                r = detail["per_config"][i]
                ff = _f_rmse(fp_[i], a) if fp_ and fp_[i] else np.nan
                nan = float("nan")
                f.write(f"{os.path.basename(a.info['folder'])} {kinds[i]} {len(a)} "
                        f"{a.info['REF_energy'] / len(a):.6f} {ds_test['fmax'][i]:.4f} "
                        f"{r['F_rmse_meV'] if r else nan:.2f} {ff:.2f} {r['dE_meV'] if r else nan:.3f} "
                        f"{r['err_max_meV'] if r else nan:.2f} "
                        f"{(r['spread_max_meV'] if r and r['spread_max_meV'] is not None else nan):.2f} "
                        f"{novelty[i] if novelty[i] is not None else nan:.3f} {outside[i] or 'in'}\n")

    worst_path = None
    if args.view_worst > 0:
        from ase.io import write as ase_write
        frames = []
        for i in order[:args.view_worst]:
            r = detail["per_config"][i]
            b = test_atoms[i].copy()
            b.calc = None
            b.new_array("F_error", r["err_atom"])
            b.set_initial_charges(r["err_atom"])
            b.info = {"folder": os.path.basename(test_atoms[i].info["folder"]), "kind": kinds[i],
                      "F_rmse_meV": round(r["F_rmse_meV"], 1), "round": str(detail_k)}
            frames.append(b)
        worst_path = near_wf(WORST_FILE)
        ase_write(worst_path, frames, format="extxyz")

    if state and rounds_only:
        for k, v in rounds_only:
            entry = mlff_workflow.get_round(state, k)
            pc = [{"folder": mlff_workflow.rel(a.info["folder"], wf_base),
                   "F_rmse_meV": round(r["F_rmse_meV"], 2), "dE_meV": round(r["dE_meV"], 3),
                   "spread_max_meV": None if r["spread_max_meV"] is None else round(r["spread_max_meV"], 2)}
                  for a, r in zip(test_atoms, v["per_config"]) if r is not None]
            entry["validation"] = {
                "time": datetime.now().isoformat(timespec="seconds"),
                "n_test": len(test_atoms),
                "overall": v["overall"], "groups": v["groups"], "members": v["members"],
                "calibration": v["calibration"],
                "phonons": {kk: vv for kk, vv in (v.get("phonons") or {}).items() if not kk.startswith("_")} or None,
                "per_config": pc,
                "elements": v.get("elements"),
                "eos": eos_store.get(k),
                "md_check": md_store.get(k),
                "foundation": ({"model": size, "overall": found["overall"], "groups": found["groups"]}
                               if found else None),
                "verdict": code if k == rounds_only[-1][0] else None}
            if k == rounds_only[-1][0]:
                if lc:
                    entry["validation"]["learning_curve"] = {kk: vv for kk, vv in lc.items()}
                if cmp_store:
                    entry["validation"]["compare"] = cmp_store
        for k, v in rounds_only:
            try:
                refresh_card_validation(state, k, wf_base, mlff_workflow.get_round(state, k)["validation"])
            except Exception as e:
                warn(f"Model card: validation not attached to round {k}'s cards ({e}).")
        try:
            best_info = export_best_model(state, rounds_only, min(
                rounds_only, key=lambda kv: kv[1]["overall"]["F_rmse_meV"])[0], wf_base, args.workflow,
                code, len(test_atoms), warn)
        except Exception as e:
            warn(f"Best model: not exported ({type(e).__name__}: {e}).")
        mlff_workflow.record(state, "validate", "stb-mlffValidate", sys.argv[1:])
        mlff_workflow.save(state, args.workflow)

    # ---- [8] files
    print_section("[11] SUMMARY & FILES", f_out)
    print_dual(f"Status            : OK, {len(scored)} committee(s) on {len(test_atoms)} test config(s) in "
               f"{time.monotonic() - t_start:.0f} s" + (f" ({len(warns)} warning(s), see [13])" if warns else ""),
               f_out)
    if code:
        print_dual(f"Verdict           : {code.upper()}", f_out)
    if best_info:
        print_dual(f"Best model        : round {best_info['round']}, seed {best_info['seed']} (test F RMSE "
                   f"{best_info['test_F_rmse_meV']:.1f} meV/Ang; its committee mean {best_info['committee_test_F_rmse_meV']:.1f})"
                   f" -> {mlff_workflow.resolve(best_info['model'], wf_base)}", f_out)
        print_dual(f"                    + {mlff_workflow.resolve(best_info['card'], wf_base)}  ({mlff_workflow.resolve(BEST_DIR, wf_base)}/ "
                   "is rewritten by every Stage 3 run)", f_out)
    elif best_models:
        print_dual(f"Best model(s)     : {_round_label(best_k)}: {best_models[0]}"
                   + (f" (+{len(best_models) - 1} more)" if len(best_models) > 1 else ""), f_out)
    if png:
        print_dual(f"Figure            : {png}", f_out)
    if data_path:
        print_dual(f"Per-config errors : {data_path}", f_out)
    if worst_path:
        print_dual(f"Worst test configs: {worst_path} (per-atom F_error; color by initial charge in ASE)", f_out)
    if state and rounds_only:
        print_dual(f"Workflow file     : {args.workflow} (scores stored per round)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    # ---- [9] next steps
    print_section("[12] NEXT STEPS", f_out)
    if best_info:
        use = (f"--custom-model {mlff_workflow.resolve(best_info['model'], wf_base)} "
               f"--custom-model-card {mlff_workflow.resolve(best_info['card'], wf_base)} in stb-mlmd, stb-mlphonons, stb-mlelastic, ...")
    else:
        use = (f"--custom-model {best_models[0]} --custom-model-card "
               f"{os.path.join(os.path.dirname(best_models[0]), 'mlff_model_card.json')} in stb-mlmd, stb-mlphonons, stb-mlelastic, ..."
               if best_models else "")
    steps = {
        "stop": [f"Use the model: {use}"],
        "regressed": [f"Use {_round_label(best_k)} ({use})",
                      "Retrain the latest data from the foundation model, with more --epochs (4.17.2) -- "
                      "without the configs far above the test set in energy, if [1] lists any."],
        "plateau": ["Change the sampling: stb-mlffActiveLearning (4.17.4) with --supercell, a higher "
                    "--md-temperature or a larger --stdev, then SIESTA, 4.17.2, 4.17.3.",
                    "Or stop and use the best round: " + use],
        "continue": ["Propose new configs: stb-mlffActiveLearning (4.17.4); run SIESTA in them.",
                     "Retrain: stb-mlffAnalysis (4.17.2, --init-from previous is faster), then validate again."],
    }.get(code, [f"Use the model: {use}"] if use else [])
    for n, s_ in enumerate(steps, 1):
        print_dual(f"{n}. {s_}", f_out)
    extras = [t for t, used in (("--show (display the figure)", show_requested),
                                ("--save-data (per-config errors)", args.save_data),
                                ("--view-worst N (worst test configs in ASE's viewer)", args.view_worst),
                                ("--phonon-ref DIR (phonons vs SIESTA)", args.phonon_ref),
                                ("--target-force F (a stop criterion)", args.target_force is not None),
                                ("--compare A B (two rounds config by config)", args.compare),
                                ("--md-check (error along an MD run)", args.md_check or md_atoms or md_pending))
              if not used]
    if extras:
        print_dual(f"Also available    : {'; '.join(extras)}.", f_out)

    # ---- [10] warnings
    print_section("[13] WARNINGS", f_out)
    if warns:
        for k, w in enumerate(warns, 1):
            print_dual(color_text(f"{k:>2}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)
    if library_noise:
        print_dual(f"Library messages  : {len(library_noise)} informational message(s) from mace/torch/phonopy "
                   "captured (model loading, dtype).", f_out)

    if f_out:
        f_out.close()
    print("\n" + "-" * 60)
    print(color_text("Validation done.\n", 'bold'))
    if worst_path:
        from stb.core.ase_view import view_structure_interactive
        from ase.io import read as ase_read
        print(color_text(f"[INFO] Opening the {args.view_worst} worst test config(s) in ASE's viewer "
                         "(View > Colors > By initial charge = per-atom force error, eV/Ang).", 'cyan'))
        view_structure_interactive(ase_read(worst_path, ":"))


if __name__ == "__main__":
    main()
