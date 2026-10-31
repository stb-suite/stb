#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import glob
import math
import time
import shutil
import argparse
import tempfile
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, mlff_workflow
from stb.core.fdf_directives import parse_kgrid, parse_directives
from stb.core.pseudopotentials import BANKS, get_required_pseudos, resolve_pseudo_source
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                          print_progress_line, finish_progress_line, capture_library_noise)
from stb.mlff import build_rattled_configs, write_config_folder, resolve_reference
from stb.core.atomic_env import EnvironmentSet, FP_NEIGHBOURS
from stb.mlmd import build_dynamics

REPORT_FILE = "stb_mlffActiveLearning_report.txt"


DEFAULT_MIN_NOVELTY = 0.5
# Seeds: rattle candidates use their own offset (Stage 1 uses seed + i), and the
# base seed moves by ROUND_SEED_STEP per round, so every round screens a fresh,
# still reproducible pool instead of the same geometries again.
RATTLE_SEED_OFFSET = 30000
ROUND_SEED_STEP = 1000
# A pick whose shortest interatomic distance is below this fraction of the
# training set's shortest one is a close contact SIESTA may not converge on.
_CLOSE_CONTACT_FRAC = 0.85
# Below this share of the pool being new (novelty >= --min-novelty) the
# sampling is close to saturated.
_SATURATION_PCT = 20.0
_DMIN_CUTOFF = 4.0


def next_config_index(output_dir):
    """Returns 1 + the highest NNN already used by an mlff_config_NNN/ folder
    under `output_dir` -- in its configs_round_*/ subfolders or (older layout)
    directly in it -- so a batch continues the project-wide numbering instead
    of overwriting a folder."""
    return mlff_workflow.next_config_index(output_dir)


# --------------------------------------------------------------------------
# Scoring: committee disagreement (2+ models) or max predicted force (1 model)
# --------------------------------------------------------------------------

def evaluate_candidate(atoms, calcs):
    """Single-point forces of every model in `calcs` on `atoms`. Returns
    {'max_f': max per-atom |F| of the committee-mean force (eV/Ang),
     'sigma': max over atoms of the committee force spread (eV/Ang) -- the
              RMS over models of |F_model - F_mean| per atom -- or None for a
              single model,
     'score': sigma with a committee, max_f with one model,
     'spread_atom' / 'f_atom': the per-atom spread (or None) and |F_mean|,
     'energy': committee-mean energy per atom (eV)}
    or None if any model fails (e.g. a NaN force).

    The committee spread is the actual uncertainty signal: models trained
    on the same data with different seeds agree where the data constrain
    them and drift apart where they do not. One model's own max force only
    says the geometry is far from equilibrium, which selects the most
    distorted candidates rather than the ones the model gets wrong.
    """
    forces, energies = [], []
    for calc in calcs:
        a = atoms.copy()
        a.calc = calc
        try:
            f = a.get_forces()
            e = a.get_potential_energy()
        except Exception:
            return None
        if not np.all(np.isfinite(f)):
            return None
        forces.append(f)
        energies.append(e / len(a))
    forces = np.array(forces)
    mean = forces.mean(axis=0)
    f_atom = np.linalg.norm(mean, axis=1)
    sigma = spread = None
    if len(calcs) > 1:
        spread = np.sqrt(np.mean(np.sum((forces - mean) ** 2, axis=2), axis=0))
        sigma = float(spread.max())
    max_f = float(f_atom.max())
    return {"max_f": max_f, "sigma": sigma, "score": sigma if sigma is not None else max_f,
            "sigma_e": float(np.std(energies)) if len(calcs) > 1 else None, "e_models": list(energies),
            "spread_atom": spread, "f_atom": f_atom, "energy": float(np.mean(energies))}


def rattle_candidates(base_atoms, n_candidates, stdevs, seed):
    """Independent Gaussian-displaced copies of `base_atoms` (the reference
    cell, or its --supercell), stdevs cycled round-robin. Seeded with
    seed + RATTLE_SEED_OFFSET + i: Stage 1 rattles config i with seed + i
    through the same generator, and ASE draws stdev x one fixed normal vector
    per seed -- without the offset candidate i would be training config i
    displaced along the very same direction, only scaled."""
    positions = build_rattled_configs(base_atoms, "cartesian", n_candidates, stdevs,
                                      seed + RATTLE_SEED_OFFSET)
    out = []
    for i, pos in enumerate(positions):
        a = base_atoms.copy()
        a.positions = pos
        out.append((a, f"rattle {stdevs[i % len(stdevs)]:g} Ang"))
    return out


def strain_candidates(base_atoms, n_candidates, strain_max_pct, stdevs, seed):
    """Candidates with a random cell deformation: a symmetric strain tensor
    with every component uniform in +/- strain_max (off-diagonal halved),
    only along periodic axes (a vacuum axis keeps its length and angles),
    atoms carried with the cell (same fractional coordinates) and then
    rattled with `stdevs` cycled. Teaches volume and shear response (stress)
    without going back to Stage 1."""
    from stb.core.kspace import detect_vacuum_axes
    vac = detect_vacuum_axes(base_atoms.get_scaled_positions(), np.asarray(base_atoms.cell), 5.0)
    per = np.array([not v for v in vac], float)
    mask = np.outer(per, per)
    rng = np.random.default_rng(seed)
    scaled0 = base_atoms.get_scaled_positions(wrap=False)
    v0 = base_atoms.get_volume()
    out = []
    for i in range(n_candidates):
        e = rng.uniform(-1, 1, (3, 3)) * strain_max_pct / 100.0
        e = 0.5 * (e + e.T)
        e[~np.eye(3, dtype=bool)] *= 0.5
        e *= mask
        a = base_atoms.copy()
        a.set_cell(np.asarray(base_atoms.cell) @ (np.eye(3) + e), scale_atoms=False)
        a.set_scaled_positions(scaled0)
        stdev = stdevs[i % len(stdevs)]
        a.rattle(stdev=stdev, seed=seed + 20000 + i)
        dv = 100 * (a.get_volume() / v0 - 1)
        out.append((a, f"strain {strain_max_pct:g} % (rattle {stdev:g} Ang) dV {dv:+.1f} %"))
    return out


def center_energy_spread(evaluations):
    """Energy spread without each member's constant offset: the members are
    fitted with their own E0s, so their energies differ by a constant that
    says nothing about this candidate. Each model's mean deviation from the
    committee mean over the pool is removed before the spread is taken."""
    ok = [ev for ev in evaluations if ev.get("e_models") and len(ev["e_models"]) > 1]
    if len(ok) < 2:
        return
    e = np.array([ev["e_models"] for ev in ok])          # (n_cand, n_models)
    offsets = (e - e.mean(axis=1, keepdims=True)).mean(axis=0)
    for ev, row in zip(ok, e - offsets):
        ev["sigma_e"] = float(np.std(row))


def apply_score_mode(evaluations, mode, committee):
    """Sets each evaluation's 'score' for --score: 'force' (max per-atom force
    spread; max|F| with one model), 'energy' (spread of the energy per atom
    across the committee, eV/atom) or 'both' (each spread divided by its pool
    median, summed -- a candidate counts when either is unusual). Returns the
    mode actually used (energy/both need a committee)."""
    if not committee:
        for ev in evaluations:
            ev["score"] = ev["max_f"]
        return "force"
    center_energy_spread(evaluations)
    sf = np.array([ev["sigma"] for ev in evaluations])
    se = np.array([ev["sigma_e"] for ev in evaluations])
    if mode == "energy":
        for ev in evaluations:
            ev["score"] = ev["sigma_e"]
    elif mode == "both":
        mf, me = max(np.median(sf), 1e-12), max(np.median(se), 1e-12)
        for ev in evaluations:
            ev["score"] = ev["sigma"] / mf + ev["sigma_e"] / me
    else:
        for ev in evaluations:
            ev["score"] = ev["sigma"]
    return mode


def select_diverse(candidates, evaluations, top_k, reference, scale, min_novelty, max_force, diversity,
                   pool_factor=20):
    """Batch selection that trades score for diversity: at every step the
    candidate maximizing  score / max_score * min(1, novelty / diversity)  is
    picked, novelty measured against the training set AND the picks so far --
    so a candidate next to an earlier pick loses value instead of being
    taken as a near-duplicate. Only the `pool_factor` x top_k highest-score
    candidates are considered. Same filters and return as select_candidates."""
    order = sorted(range(len(candidates)), key=lambda i: evaluations[i]["score"], reverse=True)
    eligible, n_force = [], 0
    for i in order:
        if max_force and evaluations[i]["max_f"] > max_force:
            n_force += 1
        elif not evaluations[i].get("too_high"):
            eligible.append(i)
    consider = eligible[:max(pool_factor * top_k, 50)]
    smax = max((evaluations[i]["score"] for i in consider), default=1.0) or 1.0
    picked, rejected = [], set()
    while len(picked) < top_k:
        best, best_u, best_nov = None, -1.0, None
        taken = {p[0] for p in picked}
        for i in consider:
            if i in taken or i in rejected:
                continue
            nov = reference.novelty(candidates[i][0], scale) if scale else math.inf
            if min_novelty > 0 and nov < min_novelty:
                rejected.add(i)
                continue
            u = evaluations[i]["score"] / smax * (min(1.0, nov / diversity) if scale and diversity > 0 else 1.0)
            if u > best_u:
                best, best_u, best_nov = i, u, nov
        if best is None:
            break
        picked.append((best, evaluations[best], best_nov))
        reference.add(candidates[best][0])
    return picked, n_force, len(rejected)


def scan_values(method, given):
    """Settings tried by --auto-sampling (stdev Ang / temperature K / strain %)."""
    if given:
        return list(given)
    return {"rattle": [0.03, 0.05, 0.1, 0.15, 0.2, 0.3],
            "md": [300.0, 600.0, 900.0, 1200.0, 1500.0],
            "strain": [1.0, 2.0, 3.0, 5.0, 8.0]}[method]


KB_MEV = 0.08617333   # Boltzmann constant, meV/K


def equipartition_t(de_mev):
    """Temperature (K) whose harmonic potential energy, 3/2 kT per atom, is
    `de_mev` meV/atom above the minimum -- a physical reading of how hard a
    geometry is distorted (graphene melts near 4500 K = ~580 meV/atom)."""
    return de_mev / (1.5 * KB_MEV)


def choose_scanned(rows, max_force, max_close_pct=10.0, max_energy=None, train_de=None, frontier=1.5):
    """The scanned setting to sample with. 'Physical' settings: MD did not
    blow up, at most `max_close_pct` % close contacts, median max|F| within
    --max-force and median energy above the reference within --max-energy (if
    set). The committee spread grows with distortion, so the highest score
    alone always sits at the edge of the physical range -- a 0.15 Ang rattle of
    graphene, 1-2 eV/atom picks that made the model worse. So, when the
    training set's own energy range is known (`train_de`, its 90th-percentile
    energy above the reference, meV/atom -- robust to a few outliers), the choice is the highest-score setting whose
    median energy stays within `frontier` x that: the training range grows
    step by step, round after round, instead of jumping far from it. None
    physical -> the mildest setting. Returns (value, reason)."""
    ok = [r for r in rows if not r["stopped"] and r["close_pct"] <= max_close_pct
          and (not max_force or r["maxf_med"] <= max_force)
          and (not max_energy or r.get("de_med") is None or r["de_med"] <= max_energy) and r["n"] > 0]
    crit = (f"<= {max_close_pct:g} % close contacts, no MD blow-up"
            + (f", median max|F| <= {max_force:g} eV/Ang" if max_force else "")
            + (f", median E <= {max_energy:g} meV/atom above the reference" if max_energy else ""))
    if not ok:
        return rows[0]["value"], "no setting is physical by these criteria -- the mildest one"
    if train_de and train_de > 0 and frontier and all(r.get("de_med") is not None for r in ok):
        limit = frontier * train_de
        within = [r for r in ok if r["de_med"] <= limit]
        if within:
            best = max(within, key=lambda r: r["score_med"])
            return best["value"], (f"highest median score among the settings whose median energy stays within "
                                   f"--frontier {frontier:g} x the training set's 90th-percentile energy ({train_de:.0f} -> "
                                   f"{limit:.0f} meV/atom above the reference), with {crit}")
        best = min(ok, key=lambda r: r["de_med"])
        return best["value"], (f"every physical setting lies beyond --frontier {frontier:g} x the training set's "
                               f"90th-percentile energy ({limit:.0f} meV/atom) -- the mildest one")
    best = max(ok, key=lambda r: r["score_med"])
    return best["value"], f"highest median score among the settings with {crit}"


def md_candidates(base_atoms, calc, n_steps, temperatures_k, stride, seed, friction_fs=0.01, info=None,
                  ensemble="nvt", pressure_gpa=0.0):
    """Frames from one short NVT trajectory per temperature, driven by
    `calc` (the first model) through stb-mlmd's own build_dynamics -- same
    verified fs/GPa unit conversions, not a re-implementation. Every
    `stride`-th frame is kept as a candidate; scoring is done afterwards by
    every committee member.

    Temperatures are deliberately meant to be higher than a production
    temperature, purely to cover more configuration space per step (a
    common MLIP training-set trick, not a claim about the material's real
    operating conditions). Several temperatures widen the sampled amplitudes.
    A poorly-fitted model can make this MD numerically unstable -- a step
    that raises (e.g. a NaN force) ends that trajectory, keeping the frames
    already collected, rather than aborting the whole scan. With `info` (a
    dict), each temperature's steps run, mean instantaneous temperature and
    the reason it stopped early (or None) are stored there.
    """
    out = []
    for t_index, temperature in enumerate(temperatures_k):
        atoms = base_atoms.copy()
        atoms.calc = calc
        dyn = build_dynamics(atoms, ensemble, timestep_fs=1.0, temperature_k=temperature,
                             friction_fs=friction_fs, pressure_gpa=pressure_gpa, seed=seed + t_index)
        temps, stopped, done = [], None, 0
        for step in range(1, n_steps + 1):
            try:
                dyn.run(1)
            except Exception as e:
                stopped = f"step {step}: {type(e).__name__}"
                break
            done = step
            if step % stride == 0:
                if not np.all(np.isfinite(atoms.get_positions())):
                    stopped = f"step {step}: non-finite positions"
                    break
                temps.append(float(atoms.get_temperature()))
                frame = atoms.copy()
                frame.calc = None
                frame.wrap()
                out.append((frame, f"md{'-npt' if ensemble == 'npt' else ''} {temperature:g} K step {step}"))
        if info is not None:
            info[temperature] = {"steps": done, "frames": len(temps),
                                 "T_mean": float(np.mean(temps)) if temps else None,
                                 "T_std": float(np.std(temps)) if temps else None, "stopped": stopped}
    return out


# --------------------------------------------------------------------------
# Novelty / diversity filter
# --------------------------------------------------------------------------

def read_training_configs(dirs, ref_basename):
    """ASE Atoms of every folder in `dirs` (the configs already labelled, or
    written and waiting to be), read from their own .fdf.
    Returns (atoms_list, n_unreadable)."""
    atoms_list, n_bad = [], 0
    for d in dirs:
        if not os.path.isdir(d):
            continue
        fdf = os.path.join(d, ref_basename)
        if not os.path.isfile(fdf):
            others = [f for f in sorted(glob.glob(os.path.join(d, "*.fdf")))
                      if os.path.basename(f) != "config_extra.fdf"]
            fdf = others[0] if others else None
        try:
            atoms_list.append(AseAtomsAdaptor.get_atoms(
                structure_io.to_pymatgen(structure_io.read_fdf(fdf))))
        except Exception:
            n_bad += 1
    return atoms_list, n_bad


def select_candidates(candidates, evaluations, top_k, reference, scale, min_novelty, max_force):
    """Greedy selection by decreasing score: a candidate is skipped when its
    predicted max|F| exceeds `max_force` (if set), or when its novelty --
    relative to the training set AND to the candidates already picked -- is
    below `min_novelty` (0 disables). Every pick is added to `reference`, so
    the next one must also differ from it. Returns (picked, n_force, n_redundant),
    picked = list of (candidate index, evaluation, novelty)."""
    order = sorted(range(len(candidates)), key=lambda i: evaluations[i]["score"], reverse=True)
    picked, n_force, n_redundant = [], 0, 0
    for i in order:
        if len(picked) >= top_k:
            break
        ev = evaluations[i]
        if max_force and ev["max_f"] > max_force:
            n_force += 1
            continue
        if ev.get("too_high"):          # above --max-energy, counted by the caller
            continue
        atoms = candidates[i][0]
        nov = reference.novelty(atoms, scale) if scale else math.inf
        if min_novelty > 0 and nov < min_novelty:
            n_redundant += 1
            continue
        picked.append((i, ev, nov))
        reference.add(atoms)
    return picked, n_force, n_redundant


def supercell_kgrid_block(calc_text, dims):
    """(config_extra block, note) for a supercell of the reference: the
    template's own Monkhorst-Pack divisions divided by `dims` (rounded up),
    so the supercell samples the Brillouin zone like the reference cell
    does and its energies stay on the same footing as the training set."""
    kgrid = parse_kgrid(calc_text)
    if kgrid is None:
        if "kgridcutoff" in parse_directives(calc_text):
            return "", ("template sets kgrid.Cutoff -- SIESTA sizes the supercell grid itself "
                        "(same real-space cutoff), no override needed.")
        return "", ("[WARNING] template has no k-grid (Gamma only): the supercell then samples "
                    "the Brillouin zone more densely than the reference cell, so its energies "
                    "are not on exactly the same footing as the training set.")
    new = [max(1, math.ceil(k / n)) for k, n in zip(kgrid, dims)]
    exact = all(k % n == 0 for k, n in zip(kgrid, dims))
    block = ("# Supercell k-grid (stb-mlffActiveLearning --supercell): the reference's\n"
             f"# {kgrid[0]}x{kgrid[1]}x{kgrid[2]} divided by {dims[0]}x{dims[1]}x{dims[2]}.\n"
             f"kgrid.MonkhorstPack [{new[0]} {new[1]} {new[2]}]\n\n")
    note = (f"k-grid {kgrid[0]}x{kgrid[1]}x{kgrid[2]} -> {new[0]}x{new[1]}x{new[2]} "
            + ("(equivalent sampling)" if exact else
               "([WARNING] not divisible -- rounded up, slightly denser than the reference)"))
    return block, note


# --------------------------------------------------------------------------
# Analysis helpers
# --------------------------------------------------------------------------

def shortest_distance(atoms, cutoff=_DMIN_CUTOFF):
    from ase.neighborlist import neighbor_list
    d = neighbor_list("d", atoms, cutoff)
    return float(d.min()) if len(d) else cutoff


def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return None
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _fmt(v, f="{:.3f}"):
    return "-" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f.format(v)


def _group(origin):
    """Pool group of a candidate: its sampling setting ('rattle 0.15 Ang', 'md 800 K')."""
    return origin.split(" step")[0].split(" dV")[0]


def _group_key(g):
    """Sort 'md 300 K' before 'md 1500 K'."""
    nums = re.findall(r"[0-9.]+", g)
    return (g.split()[0], float(nums[0]) if nums else 0.0)


def training_summary(train_dirs):
    """Labelled training configs (SIESTA finished): max|F| (eV/Ang), shortest
    distance, wall time per atom^3 (s), number of atoms -- the yardstick the
    picks are compared with. Empty lists when nothing is labelled."""
    from stb.core import siesta_log
    out = {"fmax": [], "dmin": [], "n_atoms": [], "wall": [], "n_labelled": 0}
    for d in train_dirs:
        label = structure_io.find_fdf_system_label(d)
        if not label:
            continue
        out_file = siesta_log.find_out_file(d, label)
        fa = os.path.join(d, f"{label}.FA")
        if not (out_file and os.path.isfile(fa)):
            continue
        forces = siesta_log.read_fa_forces(fa)
        if forces is None:
            continue
        out["n_labelled"] += 1
        out["fmax"].append(float(np.linalg.norm(forces, axis=1).max()))
        out["n_atoms"].append(len(forces))
        try:
            with open(out_file, errors="replace") as f:
                m = re.findall(r"Elapsed wall time \(sec\)\s*=\s*([0-9.]+)", f.read())
            if m:
                out["wall"].append(float(m[-1]))
        except OSError:
            pass
    return out


def plot_al(out_path, candidates, evaluations, pool_nov, picked, committee, min_novelty, max_force,
            md_info, show, score_mode="force", extras=None):
    """Six panels on the candidate pool: score distribution (picks marked),
    score vs predicted max|F|, score vs novelty, score per sampling setting,
    the per-atom uncertainty of the picks, and the score along each MD run
    (md) or the max|F| distribution (rattle)."""
    score = np.array([e["score"] for e in evaluations])
    maxf = np.array([e["max_f"] for e in evaluations])
    nov = np.array([np.nan if (v is None or math.isinf(v)) else v for v in pool_nov])
    pk = np.array([i for i, _, _ in picked], int)
    lab_s = ({"force": "Committee force spread, max over atoms (eV/Ang)",
              "energy": "Committee energy spread (eV/atom)",
              "both": "Force + energy spread (x pool medians)"}[score_mode] if committee
             else "Predicted max|F| (eV/Ang)")
    extras = {k: v for k, v in (extras or {}).items() if v}
    nrows = 2 + (1 if extras else 0)
    fig, axes = plt.subplots(nrows, 3, figsize=(16, 4.5 * nrows), squeeze=False)

    ax = axes[0, 0]
    bins = np.linspace(score.min(), score.max() * 1.001 + 1e-9, 30)
    ax.hist(score, bins=bins, color="0.7", label="pool")
    if len(pk):
        ax.hist(score[pk], bins=bins, color="tab:red", label="picked")
        ax.axvline(score[pk].min(), color="tab:red", ls=":", lw=1)
    ax.set_xlabel(lab_s); ax.set_ylabel("Candidates")
    ax.set_title("Score distribution"); ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.scatter(maxf, score, s=12, color="0.6", label="pool")
    if len(pk):
        ax.scatter(maxf[pk], score[pk], s=40, color="tab:red", label="picked", zorder=3)
    if max_force:
        ax.axvline(max_force, color="k", ls="--", lw=0.8, label="--max-force")
    ax.set_xlabel("Predicted max|F| (eV/Ang)"); ax.set_ylabel(lab_s)
    ax.set_title("Score vs force magnitude"); ax.legend(fontsize=8)

    ax = axes[0, 2]
    if np.isfinite(nov).any():
        ax.scatter(nov, score, s=12, color="0.6", label="pool")
        if len(pk):
            ax.scatter(nov[pk], score[pk], s=40, color="tab:red", label="picked", zorder=3)
        if min_novelty > 0:
            ax.axvline(min_novelty, color="k", ls="--", lw=0.8, label="--min-novelty")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "novelty filter off", ha="center", transform=ax.transAxes)
    ax.set_xlabel("Novelty vs training set (typical spacings)"); ax.set_ylabel(lab_s)
    ax.set_title("Score vs novelty")

    ax = axes[1, 0]
    groups = sorted({_group(c[1]) for c in candidates}, key=_group_key)
    data = [score[[i for i, c in enumerate(candidates) if _group(c[1]) == g]] for g in groups]
    ax.boxplot(data, showfliers=False)
    for x, (g, d) in enumerate(zip(groups, data), 1):
        ax.scatter(np.full(len(d), x) + np.random.default_rng(x).uniform(-0.15, 0.15, len(d)), d, s=8, color="0.4")
    for i in pk:
        ax.scatter(groups.index(_group(candidates[i][1])) + 1, score[i], s=40, color="tab:red", zorder=3)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels(groups, rotation=15, fontsize=8)
    ax.set_ylabel(lab_s); ax.set_title("Score per sampling setting (red = picked)")

    ax = axes[1, 1]
    if len(pk) and committee:
        for n, (i, ev, _) in enumerate(picked):
            ax.plot(np.sort(ev["spread_atom"])[::-1], "o-", ms=3, label=f"pick {n + 1}")
        ax.set_xlabel("Atoms, sorted"); ax.set_ylabel("Committee spread (eV/Ang)")
        ax.set_title("Per-atom uncertainty of the picks")
        ax.legend(fontsize=7)
    elif len(pk):
        for n, (i, ev, _) in enumerate(picked):
            ax.plot(np.sort(ev["f_atom"])[::-1], "o-", ms=3, label=f"pick {n + 1}")
        ax.set_xlabel("Atoms, sorted"); ax.set_ylabel("Predicted |F| (eV/Ang)")
        ax.set_title("Per-atom force of the picks (single model)")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "no pick", ha="center", transform=ax.transAxes)

    ax = axes[1, 2]
    if md_info:
        for g in groups:
            ii = [i for i, c in enumerate(candidates) if _group(c[1]) == g]
            steps = [int(candidates[i][1].split("step ")[1]) for i in ii]
            line, = ax.plot(steps, score[ii], "o-", ms=3, label=g)
            for i in pk:
                if i in ii:
                    ax.plot(int(candidates[i][1].split("step ")[1]), score[i], "o", ms=8, color="tab:red")
        ax.set_xlabel("MD step (fs)"); ax.set_ylabel(lab_s)
        ax.set_title("Score along each MD run (red = picked)")
        ax.legend(fontsize=8)
    else:
        bins = np.linspace(maxf.min(), maxf.max() * 1.001 + 1e-9, 30)
        ax.hist(maxf, bins=bins, color="0.7", label="pool")
        if len(pk):
            ax.hist(maxf[pk], bins=bins, color="tab:red", label="picked")
        ax.set_xlabel("Predicted max|F| (eV/Ang)"); ax.set_ylabel("Candidates")
        ax.set_title("Force magnitude of the pool"); ax.legend(fontsize=8)

    slots = list(axes[2]) if nrows > 2 else []
    for ax, (name, data) in zip(slots, extras.items()):
        if name == "scan":
            rows, chosen, unit = data
            v = [r["value"] for r in rows]
            ax.plot(v, [r["score_med"] for r in rows], "o-", color="k", label="median score")
            ax.set_xlabel(f"Setting ({unit})"); ax.set_ylabel(lab_s)
            ax2 = ax.twinx()
            ax2.plot(v, [r["close_pct"] for r in rows], "s--", color="tab:red", label="close contacts %")
            if any(r["new_pct"] is not None for r in rows):
                ax2.plot(v, [r["new_pct"] if r["new_pct"] is not None else np.nan for r in rows], "^--",
                         color="tab:blue", label="new territory %")
            ax2.set_ylabel("%"); ax2.set_ylim(-5, 105)
            ax.axvline(chosen, color="tab:green", ls=":", lw=1.5)
            h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, fontsize=7)
            ax.set_title(f"Auto-sampling scan (chosen {chosen:g})")
        elif name == "energy":
            sf, se = data
            ax.scatter(sf, 1000 * se, s=12, color="0.6", label="pool")
            if len(pk):
                ax.scatter(sf[pk], 1000 * se[pk], s=40, color="tab:red", label="picked", zorder=3)
            ax.set_xlabel("Force spread, max over atoms (eV/Ang)"); ax.set_ylabel("Energy spread (meV/atom)")
            ax.set_title("Energy vs force disagreement"); ax.legend(fontsize=8)
        elif name == "history":
            lab = [h["label"].split()[0] for h in data]
            x = np.arange(len(data))
            ax.plot(x, [h["max"] if h["max"] is not None else np.nan for h in data], "o-", label="max score")
            ax.plot(x, [h["med"] if h["med"] is not None else np.nan for h in data], "s--", label="median score")
            ax.set_xticks(x); ax.set_xticklabels(lab, rotation=20, fontsize=8)
            ax.set_xlabel("Round.run"); ax.set_ylabel("Score (as recorded)")
            ax.set_title("Active-learning history"); ax.legend(fontsize=8)
    for ax in slots[len(extras):]:
        ax.axis("off")

    fig.suptitle("stb-mlffActiveLearning -- candidate screening", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 1 - 0.035 * 2 / nrows])
    fig.savefig(out_path, dpi=140)
    if show:
        plt.show()
    plt.close(fig)


def can_show():
    import matplotlib
    return matplotlib.get_backend().lower() not in ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template")


def _versions():
    out = {}
    for mod in ("mace", "torch", "ase"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = "?"
    return out


# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Stage 4 (Active Learning) of the Custom ML Force Field workflow: uses "
                     "already fine-tuned model(s) (stb-mlffAnalysis) to screen a fresh batch "
                     "of candidate configurations (independent rattled snapshots, or frames "
                     "from short NVT trajectories driven by the model itself via "
                     "--sampling-method md, reusing stb-mlmd's own dynamics), optionally in a "
                     "--supercell of the reference, and writes the most informative ones as new "
                     "mlff_config_NNN/ folders, ready for the user to label with real SIESTA and "
                     "feed back into another round of stb-mlffAnalysis. With 2+ --model (a "
                     "committee, e.g. the same data trained with different --seed) the score is "
                     "the committee's force disagreement -- an actual uncertainty signal; with "
                     "one model it falls back to the largest predicted force, a cheap heuristic "
                     "that favours the most distorted geometries. Picks are filtered for "
                     "novelty: one whose atomic environments are all close to the training set "
                     "or to an earlier pick is skipped. The report describes the pool, the "
                     "selection funnel, each pick against the training set (extrapolation, close "
                     "contacts, estimated error from Validate's calibration) and the SIESTA "
                     "cost; every warning goes to one final section.",
        epilog="Example usage:\n"
               "  stb-mlffActiveLearning                    # latest round's committee from mlff_workflow.json\n"
               "  stb-mlffActiveLearning --dry-run --show   # screen and look, write nothing\n"
               "  stb-mlffActiveLearning --model m_s1.model m_s2.model m_s3.model \\\n"
               "                          -s structure.fdf -c calc.fdf --n-candidates 50 --stdev 0.15 --top-k 5\n"
               "  stb-mlffActiveLearning --model m_s1.model m_s2.model m_s3.model --file si.fdf \\\n"
               "                          --sampling-method md --md-temperature 300 600 900 \\\n"
               "                          --supercell 2 2 2",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--workflow", default=None, metavar="JSON",
                        help=f"Workflow state file (default: {mlff_workflow.RUN_DIR}/"
                             f"{mlff_workflow.WORKFLOW_FILE}, or ./{mlff_workflow.WORKFLOW_FILE} "
                             "inside that folder): gives the "
                             "committee, the reference .fdf, the training set for the novelty "
                             "filter and the output directory when those options are omitted; the "
                             "picked folders are recorded in it.")
    parser.add_argument("--round", type=int, default=None, metavar="K",
                        help="Take the committee of round K from --workflow (default: the latest).")
    parser.add_argument("--model", default=None, nargs="+", metavar="PATH",
                        help="Fine-tuned MACE model file(s), overriding --workflow's committee. Two "
                             "or more form a committee: the score becomes their force disagreement "
                             "(stb-mlffAnalysis --committee N trains one). With one model the score "
                             "is its largest predicted force.")
    parser.add_argument("-s", "--structure", default=None,
                        help="Reference structure .fdf (default: the one in --workflow), rewritten "
                             "with each pick's coordinates; with -c/--calc, the calc.fdf template "
                             "(default: calc.fdf) that %%includes it, copied into every folder.")
    parser.add_argument("-c", "--calc", default=None,
                        help="Calculation template that %%includes --structure (default: calc.fdf "
                             "with --structure, else the one in --workflow).")
    parser.add_argument("--file", default=None,
                        help="Instead of -s/-c: one self-contained reference .fdf -- basis/"
                             "pseudos/SCF settings are reused verbatim in every candidate written out.")
    parser.add_argument("--pseudo-dir", default=None,
                        help="Directory (or bundled bank) to search for <element>.psf/.psml "
                             "pseudopotentials (default: the one stb-mlff recorded in --workflow, "
                             "else the current directory).")
    parser.add_argument("--supercell", type=int, nargs=3, default=[1, 1, 1], metavar=("NA", "NB", "NC"),
                        help="Sample candidates in an NA x NB x NC supercell of --file (default: "
                             "1 1 1, the reference cell itself). A small training cell only "
                             "constrains the interactions its own periodicity allows -- e.g. an "
                             "8-atom diamond cell sees the Gamma and X phonons but not L -- and a "
                             "supercell adds the rest. Each folder gets the reference k-grid "
                             "divided by the supercell size in its config_extra.fdf. SIESTA cost "
                             "grows steeply with the atom count.")
    parser.add_argument("--sampling-method", choices=["rattle", "md", "strain"], default="rattle",
                        help="How to generate candidate configurations (default: rattle):\n"
                             "  rattle - independent, static Gaussian-displaced snapshots "
                             "(cheap: 1 evaluation per candidate and model).\n"
                             "  md     - short NVT trajectories with the (first) model itself "
                             "(stb-mlmd's own dynamics, reused directly); more physically "
                             "representative of the configuration space the system actually "
                             "explores, at the cost of running real MD steps (--md-ensemble npt lets the "
                             "cell breathe too).\n"
                             "  strain - random cell deformations (volume and shear, up to --strain-max, "
                             "periodic axes only), rattled: the volume/stress response.")
    parser.add_argument("--strain-max", type=float, default=3.0, metavar="PCT",
                        help="Largest strain component, %% (--sampling-method strain only, default: 3).")
    parser.add_argument("--md-ensemble", choices=["nvt", "npt"], default="nvt",
                        help="MD ensemble for --sampling-method md (default: nvt); npt (Berendsen, bulk only, "
                             "--md-pressure) samples the cell as well, so the picks carry their own lattice.")
    parser.add_argument("--md-pressure", type=float, default=0.0, metavar="GPA",
                        help="Pressure for --md-ensemble npt, GPa (default: 0).")
    parser.add_argument("--score", choices=["force", "energy", "both"], default="force",
                        help="What the committee's disagreement is measured on (default: force):\n"
                             "  force  - max over atoms of the per-atom force spread\n"
                             "  energy - spread of the energy per atom (catches configs whose forces agree but "
                             "whose energy does not, e.g. another volume)\n"
                             "  both   - force and energy spreads, each divided by its pool median, summed")
    parser.add_argument("--selection", choices=["diverse", "greedy"], default="diverse",
                        help="How the top-k batch is chosen (default: diverse):\n"
                             "  diverse - score x min(1, novelty/--diversity), novelty vs the training set "
                             "and the picks so far: no near-duplicate picks\n"
                             "  greedy  - highest score first (only --min-novelty filters)")
    parser.add_argument("--diversity", type=float, default=1.0, metavar="D",
                        help="Novelty (typical training spacings) at which a candidate keeps its full score "
                             "in --selection diverse (default: 1).")
    parser.add_argument("--auto-sampling", action="store_true",
                        help="First scan several settings of the sampling method (stdev / temperature / "
                             "strain, --scan-values) with a few candidates each, then sample with the one "
                             "of highest median score that is still physical (few close contacts, stable MD, "
                             "--max-force, --max-energy) and whose median energy stays within --frontier x "
                             "the training set's -- one step beyond the labelled range, not its edge.")
    parser.add_argument("--scan-values", type=float, nargs="+", default=None, metavar="V",
                        help="Settings for --auto-sampling (default: stdev 0.03-0.3 Ang, 300-1500 K, or "
                             "strain 1-8 %%).")
    parser.add_argument("--scan-n", type=int, default=12,
                        help="Candidates per scanned setting (default: 12).")
    parser.add_argument("--target-force", type=float, default=None, metavar="MEV_ANG",
                        help="Error that is good enough (meV/Ang, as in stb-mlffValidate): the report says "
                             "whether any candidate is still expected above it -- if none is, another round "
                             "on this sampling adds little.")
    parser.add_argument("--n-candidates", type=int, default=50,
                        help="Number of candidate configurations to screen (--sampling-method "
                             "rattle only, default: 50).")
    parser.add_argument("--stdev", type=float, nargs="+", default=[0.15],
                        help="Rattling standard deviation(s) in Angstrom for the candidate "
                             "batch (--sampling-method rattle only, default: 0.15 -- typically "
                             "larger than the amplitude(s) used for the original training set, "
                             "to probe configuration space not yet covered).")
    parser.add_argument("--md-steps", type=int, default=500,
                        help="NVT steps per temperature (--sampling-method md only, default: 500).")
    parser.add_argument("--md-temperature", type=float, nargs="+", default=[800.0],
                        help="NVT temperature(s), K -- one trajectory each (--sampling-method md "
                             "only, default: 800 -- deliberately higher than a typical production "
                             "temperature, purely to cover more configuration space per step; "
                             "not a claim about the material's real operating conditions).")
    parser.add_argument("--stride", type=int, default=10,
                        help="Keep one frame every Nth MD step as a candidate (--sampling-method md "
                             "only, default: 10).")
    parser.add_argument("--top-k", type=int, default=5,
                        help="Number of candidates to write out (default: 5).")
    parser.add_argument("--min-novelty", type=float, default=DEFAULT_MIN_NOVELTY,
                        help="Skip a candidate whose atomic environments are all closer than this "
                             "to the training set or to an earlier pick (default: "
                             f"{DEFAULT_MIN_NOVELTY:g}; 0 disables). Measured on a model-independent "
                             f"fingerprint (sorted distances to the {FP_NEIGHBOURS} nearest neighbours) "
                             "in units of the training set's own typical spacing -- 1 means as far "
                             "as two independent training configs are from each other.")
    parser.add_argument("--train-path", default=None, metavar="GLOB",
                        help="Folders the novelty filter compares against (default: the training "
                             "configs in --workflow, else every mlff_config_* under <output-dir>, in "
                             "configs_round_*/ or directly in it -- labelled or not). Quote the pattern.")
    parser.add_argument("--max-force", type=float, default=mlff_workflow.UNPHYSICAL_FMAX, metavar="EV_ANG",
                        help="Skip candidates whose predicted max|F| exceeds this, eV/Ang -- "
                             "geometries too distorted to matter for the intended use (default: "
                             f"{mlff_workflow.UNPHYSICAL_FMAX:g}, the workflow's 'unphysical' line, the same "
                             "Stage 1 and 2 warn at; 0 = no limit).")
    parser.add_argument("--d3", action="store_true",
                        help="Add DFT-D3 (DFTD3 .true.) to every pick's config_extra.fdf. Not needed with a "
                             "workflow file: the picks follow Stage 1's --d3 recorded there; for explicit "
                             "--model/-s/-c runs without one.")
    parser.add_argument("--max-energy", type=float, default=mlff_workflow.UNPHYSICAL_DE, metavar="MEV_ATOM",
                        help="Skip candidates whose predicted energy lies more than this above the "
                             "reference's, meV/atom (committee mean; default: "
                             f"{mlff_workflow.UNPHYSICAL_DE:g}, the workflow's 'unphysical' line; 0 = no "
                             "limit). The report also gives each energy as a temperature (3/2 kT per atom).")
    parser.add_argument("--frontier", type=float, default=1.5, metavar="FACTOR",
                        help="--auto-sampling: choose the highest-score setting whose median predicted "
                             "energy stays within FACTOR x the training set's 90th-percentile energy above "
                             "the reference (default: 1.5), so each round extends the training range a step "
                             "instead of jumping to the most distorted physical setting. 0 = off (the "
                             "highest score among the physical settings).")
    parser.add_argument("--seed", type=int, default=42,
                        help="Base random seed for rattling / straining / initial MD velocities (default: 42). "
                             f"The seed used is --seed + {ROUND_SEED_STEP} x the round the picks will train, so "
                             "each round screens new candidates; rattles also get their own offset, so they never "
                             "repeat Stage 1's displacement directions.")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"],
                        help="Device to run the model(s) on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Run folder (default: --workflow's directory, else the current one). With a "
                             "workflow file the picks go to <output-dir>/configs_round_NN/ -- NN = the round "
                             "they will train, the one after the latest --, else to <output-dir> itself; the "
                             "mlff_config_NNN numbering continues project-wide.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Screen and report (figure and data files too) but write no folder and "
                             "leave the workflow file untouched -- to tune --stdev/--md-temperature/"
                             "--supercell before committing to SIESTA runs.")
    parser.add_argument("--show", action="store_true",
                        help="Also display the screening figure on screen (it is always saved).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the whole pool, one line per candidate (origin, score, spread, "
                             "max|F|, novelty, shortest distance, picked), to mlff_al_<tag>_pool.dat. "
                             "Off by default.")
    parser.add_argument("--view-picks", action="store_true",
                        help="Write the picks to mlff_al_<tag>_picks.xyz (per-atom committee spread as "
                             "'spread' and as the initial charge, so ASE's viewer colors by it: View > "
                             "Colors > By initial charge) and open them.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE} (in the output directory). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlffActiveLearning {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Active Learning)",
            "Suggest new SIESTA configurations to label, using the fine-tuned model(s)",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CUSTOM ML FORCE FIELD -- ACTIVE LEARNING:", 'bold'))
    print("-" * 60)
    t_start = time.monotonic()

    warns = []           # every warning, printed together in the last section
    library_noise = []   # third-party prints/warnings, counted there too

    def warn(msg):
        warns.append(msg)

    f_out = None
    tmp_dir = None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if warns:
            print_section("WARNINGS SO FAR", f_out)
            for k, w in enumerate(warns, 1):
                print_dual(color_text(f"  {k}. {w}", 'yellow'), f_out)
        if f_out:
            f_out.close()
        if tmp_dir:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        sys.exit(1)

    args.workflow = mlff_workflow.default_path(args.workflow)
    try:
        state = mlff_workflow.load(args.workflow)
    except ValueError as e:
        fail(str(e))
    wf_base = mlff_workflow.base_dir(args.workflow)
    round_entry = None
    models_from_round = args.model is None   # the committee IS that round's (its quality numbers apply)
    explicit_ref = bool(args.file or args.structure or args.calc)
    if args.model is None or not explicit_ref:
        if state is None:
            fail(f"'{args.workflow}' not found -- pass --model and -s/-c (or --file), or run the earlier "
                 "stages (stb-mlff, stb-mlffAnalysis) to create it.")
    if args.model is None:
        round_entry = mlff_workflow.get_round(state, args.round)
        if round_entry is None:
            avail = ", ".join(str(r["round"]) for r in (state or {}).get("rounds", [])) or "none"
            fail(f"No training round{' ' + str(args.round) if args.round else ''} in '{args.workflow}' "
                 f"(it has: {avail}) -- run stb-mlffAnalysis (7.1.2) first, or pass --model.")
        args.model = mlff_workflow.committee_models(round_entry, wf_base)
        if not args.model:
            fail(f"Round {round_entry['round']} in '{args.workflow}' has no model file.")
    elif state is not None and args.round is not None:
        round_entry = mlff_workflow.get_round(state, args.round)
    if args.file and (args.structure or args.calc):
        fail("Pass either -s/--structure + -c/--calc or --file, not both.")
    if args.file:
        ref_file, calc_template, error = resolve_reference(None, None, args.file)
    elif explicit_ref:
        ref_file, calc_template, error = resolve_reference(
            args.structure or "structure.fdf", args.calc or "calc.fdf", None)
    else:
        recorded = state.get("reference", {})
        if not recorded.get("file"):
            fail(f"'{args.workflow}' records no reference .fdf -- pass -s/-c or --file.")
        if recorded.get("calc"):
            ref_file, calc_template, error = resolve_reference(
                mlff_workflow.resolve(recorded["file"], wf_base),
                mlff_workflow.resolve(recorded["calc"], wf_base), None)
        else:
            ref_file, calc_template, error = resolve_reference(
                None, None, mlff_workflow.resolve(recorded["file"], wf_base))
    if args.output_dir is None:
        args.output_dir = os.path.relpath(wf_base) if state is not None else "."
    if args.pseudo_dir is None:
        recorded = (state or {}).get("reference", {}).get("pseudo_dir")
        if recorded:
            candidate = mlff_workflow.resolve(recorded, wf_base)
            args.pseudo_dir = recorded if recorded in BANKS else candidate
        else:
            args.pseudo_dir = "."
    proposer_round = round_entry["round"] if round_entry else (
        mlff_workflow.get_round(state)["round"] if state and state["rounds"] else None)
    if round_entry is None and state is not None and proposer_round is not None:
        round_entry = mlff_workflow.get_round(state, proposer_round)
    n_al_runs = len((round_entry or {}).get("active_learning", []))
    tag = f"r{proposer_round:02d}_{n_al_runs + 1}" if proposer_round is not None else "models"
    if not args.dry_run:   # the report and figure go to the run folder, created here if new
        os.makedirs(args.output_dir, exist_ok=True)
    out_dir_files = args.output_dir if os.path.isdir(args.output_dir) else "."
    near_out = lambda name: os.path.join(out_dir_files, name)
    # the picks train the next round: <run>/configs_round_NN/ (no workflow file: the output dir itself)
    target_round = mlff_workflow.next_round_number(state) if state is not None else None
    seed_eff = args.seed + ROUND_SEED_STEP * (target_round or 0)
    picks_root = (mlff_workflow.configs_dir(args.output_dir, target_round) if target_round
                  else args.output_dir)

    committee = len(args.model) > 1
    dims = list(args.supercell)
    is_supercell = dims != [1, 1, 1]
    if args.train_path:
        train_dirs = sorted(d for d in glob.glob(args.train_path) if os.path.isdir(d))
        train_pattern = args.train_path
    elif state is not None and mlff_workflow.folders(state, wf_base, "train"):
        train_dirs = [d for d in mlff_workflow.folders(state, wf_base, "train") if os.path.isdir(d)]
        train_pattern = f"{args.workflow} (training configs)"
    else:
        train_pattern = f"{args.output_dir}/[configs_round_*/]mlff_config_*"
        train_dirs = mlff_workflow.find_config_dirs(args.output_dir)

    report_path = near_out(REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    # seeds of the committee, from the workflow when it lists them
    seeds = {mlff_workflow.resolve(m["model"], wf_base): m.get("seed")
             for m in (round_entry or {}).get("committee", []) if m.get("model")} if state else {}
    model_labels = [f"seed {seeds[m]}" if seeds.get(m) is not None else os.path.basename(m) for m in args.model]

    print_dual(color_text("===== STB-MLFFACTIVELEARNING REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    vers = _versions()
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Versions          : stb {VERSION}, mace {vers['mace']}, torch {vers['torch']}, ase {vers['ase']}", f_out)
    print_dual(f"Device            : {args.device}", f_out)
    print_dual(f"Workflow file     : {args.workflow}" + ("" if state else " (not found)"), f_out)
    if round_entry is not None and models_from_round:
        print_dual(f"Committee from    : round {round_entry['round']} ({len(args.model)} model(s): "
                   f"{', '.join(model_labels)})", f_out)
    elif round_entry is not None:
        print_dual(f"Committee from    : --model ({len(args.model)} model(s)); picks recorded as proposed by round "
                   f"{round_entry['round']}", f_out)
    print_dual(f"Model(s)          : {args.model[0]}" + (f" (+{len(args.model) - 1} more)" if committee else ""), f_out)
    if calc_template:
        print_dual(f"Reference         : {ref_file} + {calc_template}", f_out)
    else:
        print_dual(f"Reference         : {ref_file}", f_out)
    print_dual(f"Sampling          : {args.sampling_method}" + (
        f", {args.n_candidates} candidates, stdev {' '.join(f'{s:g}' for s in args.stdev)} Ang"
        if args.sampling_method == "rattle" else
        f", {args.n_candidates} candidates, strain up to {args.strain_max:g} %, rattle "
        f"{' '.join(f'{s:g}' for s in args.stdev)} Ang" if args.sampling_method == "strain" else
        f", {args.md_steps} {args.md_ensemble.upper()} steps at {', '.join(f'{t:g}' for t in args.md_temperature)} K"
        + (f" and {args.md_pressure:g} GPa" if args.md_ensemble == "npt" else "")
        + f", one frame every {args.stride}") + f", supercell {dims[0]}x{dims[1]}x{dims[2]}, seed {seed_eff}"
        + (f" (--seed {args.seed} + {ROUND_SEED_STEP} x round {target_round})" if target_round else "")
        + (" -- the setting is chosen by --auto-sampling" if args.auto_sampling else ""), f_out)
    print_dual(f"Score             : {args.score} disagreement" if committee else
               "Score             : largest predicted force (single model)", f_out)
    print_dual(f"Selection         : {args.selection}" + (f" (diversity {args.diversity:g})" if args.selection == "diverse"
                                                         else "") + f", top {args.top_k}; novelty " + (
        f">= {args.min_novelty:g} vs the training set and earlier picks" if args.min_novelty > 0 else "off")
        + "; max|F| " + (f"<= {args.max_force:g} eV/Ang" if args.max_force else "no limit"), f_out)
    print_dual(f"Output            : {picks_root}/" + (f" (the configs that will train round {target_round})"
                                                     if target_round else "")
               + (" -- DRY RUN, nothing is written" if args.dry_run else ""), f_out)
    if args.target_force is not None:
        print_dual(f"Target            : {args.target_force:g} meV/Ang (is another round worth it? -- [5])", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    # ---- [1] reference, committee, training set
    print_section("[1] REFERENCE, COMMITTEE & TRAINING SET", f_out)
    if error:
        fail(error)
    for m in args.model:
        if not os.path.isfile(m):
            fail(f"Model file '{m}' not found.")
    for name, value in (("--top-k", args.top_k), ("--n-candidates", args.n_candidates),
                        ("--md-steps", args.md_steps), ("--stride", args.stride), ("--scan-n", args.scan_n)):
        if value < 1:
            fail(f"{name} must be >= 1 (got {value}).")
    if min(dims) < 1:
        fail(f"--supercell dimensions must be >= 1 (got {dims}).")
    try:
        ref_structure = structure_io.read_fdf(ref_file)
    except Exception as e:
        fail(f"Could not read '{ref_file}': {e}")
    unique_elements = sorted(ref_structure.species)

    try:
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):   # its bank banner; the source is printed below
            pseudo_path = resolve_pseudo_source(args.pseudo_dir)
    except ValueError as e:
        fail(str(e))
    pseudos_to_copy, missing = get_required_pseudos(unique_elements, pseudo_path)
    if missing:
        fail(f"Missing pseudopotential(s) for: {', '.join(missing)} (searched '{args.pseudo_dir}').")

    ase_atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(ref_structure))
    ref_basename = os.path.basename(ref_file)
    template = ref_file
    extra_block = ""
    print_dual(f"Reference cell    : {ase_atoms.get_chemical_formula()} ({len(ase_atoms)} atoms), cell "
               f"{' x '.join(f'{l:.3f}' for l in ase_atoms.cell.lengths())} Ang, shortest distance "
               f"{shortest_distance(ase_atoms):.3f} Ang", f_out)
    print_dual(f"Pseudopotentials  : {', '.join(unique_elements)} from '{args.pseudo_dir}'", f_out)
    if is_supercell:
        tmp_dir = tempfile.mkdtemp(prefix="stb_mlffal_")
        template = os.path.join(tmp_dir, ref_basename)
        try:
            n_super = structure_io.replicate_fdf_cell(ref_file, dims, template)
        except Exception as e:
            fail(f"Could not build the {dims} supercell of '{ref_file}': {e}")
        # SIESTA reads the calc file first, then the structure it %includes.
        settings_text = ""
        for path in ([calc_template] if calc_template else []) + [ref_file]:
            with open(path) as f:
                settings_text += f.read() + "\n"
        extra_block, k_note = supercell_kgrid_block(settings_text, dims)
        if calc_template:
            # The structure file's own NumberOfAtoms was updated, a calc.fdf one was not.
            extra_block += f"NumberOfAtoms {n_super}\n\n"
        print_dual(f"Supercell         : {dims[0]}x{dims[1]}x{dims[2]} = {n_super} atoms; "
                   f"{k_note.replace('[WARNING] ', '')}", f_out)
        if "[WARNING]" in k_note:
            warn(f"Supercell: {k_note.replace('[WARNING] ', '')}")
    # dispersion: the same as every other folder of the workflow (Stage 1 --d3)
    d3 = mlff_workflow.dftd3_info(state)
    if d3.get("forced") or args.d3:
        extra_block += mlff_workflow.DFTD3_BLOCK
        print_dual("Dispersion        : DFT-D3 in every pick's config_extra.fdf ("
                   + ("as Stage 1 --d3 recorded in the workflow file" if d3.get("forced") else "--d3") + ")", f_out)
        if args.d3 and state is not None and d3 and not d3.get("on"):
            warn("Dispersion: --d3 given, but the workflow's labels have NO DFT-D3 (Stage 1 ran without it) -- "
                 "these picks would sit on another energy surface; drop --d3, or redo Stage 1 with --d3.")
    elif d3.get("on"):
        print_dual("Dispersion        : DFT-D3 from the calc template (as in every folder of the workflow)", f_out)
    base_atoms = ase_atoms.repeat(dims)
    if args.sampling_method == "md" and args.md_ensemble == "npt":
        from stb.core.kspace import detect_vacuum_axes
        if any(detect_vacuum_axes(base_atoms.get_scaled_positions(), np.asarray(base_atoms.cell), 5.0)):
            fail("--md-ensemble npt needs a bulk cell (this one has a vacuum axis) -- use nvt, or "
                 "--sampling-method strain, which leaves vacuum axes alone.")
    template_lc = structure_io.read_fdf(template).lattice_constant
    if not committee and args.score != "force":
        warn(f"Score: --score {args.score} needs a committee -- the single model's max|F| is used.")

    try:
        with capture_library_noise(library_noise, "import mace"):
            from mace.calculators import MACECalculator
    except ImportError:
        fail("mace-torch not found. This tool needs the optional 'ml' extra: pip install stb_suite[ml]")
    with capture_library_noise(library_noise, "load models"):
        calcs = [MACECalculator(model_paths=[m], device=args.device) for m in args.model]

    # what the earlier stages say about this committee
    calib = None
    if round_entry is not None and not models_from_round:
        print_dual(f"Committee quality : --model given -- the scores stored for round {round_entry['round']} do not "
                   "apply to these models (no calibration)", f_out)
    if round_entry is not None and models_from_round:
        te = round_entry.get("training_errors") or {}
        val = round_entry.get("validation") or {}
        oos = (te.get("out_of_sample") or {}).get("F_rmse_meV")
        test_f = (val.get("overall") or {}).get("F_rmse_meV")
        calib = (val.get("calibration") or {}).get("error_over_spread")
        print_dual(f"Committee quality : out-of-sample F {_fmt(oos, '{:.1f}')} meV/Ang (7.1.2), test F "
                   f"{_fmt(test_f, '{:.1f}')} meV/Ang (7.1.3)" + (f", verdict {val['verdict'].upper()}"
                                                                    if val.get("verdict") else ""), f_out)
        if calib:
            print_dual(f"Calibration       : true error ~ {calib:.1f} x spread on the test set (7.1.3) -- used "
                       "below to turn a spread into an error estimate", f_out)
        elif committee:
            print_dual("Calibration       : not known yet (run stb-mlffValidate on this round to turn "
                       "spreads into error estimates)", f_out)
        if (val.get("verdict") or "") == "stop":
            warn("Committee: stb-mlffValidate already said STOP for this round (target met) -- more data is "
                 "optional.")
        if (val.get("verdict") or "") == "regressed":
            warn(f"Committee: stb-mlffValidate marked round {round_entry['round']} as REGRESSED -- picking with "
                 "a worse committee; consider --round <best round>.")

    train_atoms, n_bad = read_training_configs(train_dirs, ref_basename)
    tr = training_summary(train_dirs)
    print_dual(f"Training set      : {len(train_atoms)} config(s) from '{train_pattern}'"
               + (f", {tr['n_labelled']} labelled by SIESTA" if train_dirs else "")
               + (f" ({n_bad} unreadable)" if n_bad else ""), f_out)
    tr_dmin = [shortest_distance(a) for a in train_atoms]
    if tr["fmax"]:
        print_dual(f"Training range    : SIESTA max|F| {min(tr['fmax']):.2f}-{max(tr['fmax']):.2f} eV/Ang, shortest "
                   f"distance {min(tr_dmin):.3f}-{max(tr_dmin):.3f} Ang, "
                   f"{min(tr['n_atoms'])}-{max(tr['n_atoms'])} atoms", f_out)
    if len(train_atoms) - tr["n_labelled"] > 0 and train_dirs:
        print_dual(f"                    ({len(train_atoms) - tr['n_labelled']} written but not labelled yet -- "
                   "they count for novelty)", f_out)
    reference = EnvironmentSet()
    for a in train_atoms:
        reference.add(a)

    # energy scale: the committee-mean energy of the (undistorted) reference
    # cell, and the training configs' predicted energies above it
    def _comp(a):
        sym = a.get_chemical_symbols()
        return tuple(sorted((el, round(sym.count(el) / len(sym), 4)) for el in set(sym)))
    e_ref, train_de, train_de90 = None, None, None
    with capture_library_noise(library_noise, "reference energy"):
        ev_ref = evaluate_candidate(base_atoms, calcs)
        if ev_ref is not None:
            e_ref = ev_ref["energy"]
            same = [a for a in train_atoms if _comp(a) == _comp(base_atoms)]
            if len(same) > 200:
                same = [same[i] for i in np.random.default_rng(0).choice(len(same), 200, replace=False)]
            des = [1000 * (ev["energy"] - e_ref) for ev in (evaluate_candidate(a, calcs) for a in same)
                   if ev is not None]
            train_de = max(des) if des else None
            train_de90 = float(np.percentile(des, 90)) if des else None
    if train_de is not None:
        print_dual(f"Training energies : 90 % within {train_de90:.0f}, all within {train_de:.0f} meV/atom above the "
                   f"reference (model-predicted; ~{equipartition_t(max(train_de90, 0)):.0f} / "
                   f"~{equipartition_t(max(train_de, 0)):.0f} K as 3/2 kT per atom)", f_out)
    # the energy cap of every pick: --max-energy, and --frontier x the training
    # set's 90th-percentile energy (each round extends the labelled range a step)
    e_cap, e_cap_why = (args.max_energy or None), (f"--max-energy {args.max_energy:g}" if args.max_energy else "")
    if args.frontier and train_de90 and train_de90 > 0 and (e_cap is None or args.frontier * train_de90 < e_cap):
        e_cap = args.frontier * train_de90
        e_cap_why = f"--frontier {args.frontier:g} x the training set's 90th percentile {train_de90:.0f}"
    if e_cap is not None and e_ref is not None:
        print_dual(f"Energy cap        : picks at most {e_cap:.0f} meV/atom above the reference "
                   f"(~{equipartition_t(e_cap):.0f} K as 3/2 kT; {e_cap_why})", f_out)
    if not committee:
        warn("Score: single model -- the score is its largest predicted force, which favours the most "
             "distorted geometries rather than the ones it gets wrong. Train a committee (7.1.2 "
             "--committee 3 or --cv 5) for an uncertainty score.")

    # ---- [2] candidate pool
    print_section("[2] CANDIDATE POOL", f_out)

    def make_pool(value=None, n=None, steps=None, seed=None):
        """(candidates, md_info) for the sampling method; `value` replaces the
        stdev / temperature / strain setting (the --auto-sampling scan)."""
        info, sd, n = {}, seed_eff if seed is None else seed, n or args.n_candidates
        with capture_library_noise(library_noise, "candidate generation"):
            if args.sampling_method == "rattle":
                c = rattle_candidates(base_atoms, n, [value] if value is not None else args.stdev, sd)
            elif args.sampling_method == "strain":
                c = strain_candidates(base_atoms, n, value if value is not None else args.strain_max,
                                      args.stdev, sd)
            else:
                temps = [value] if value is not None else args.md_temperature
                st = steps or args.md_steps
                stride = args.stride if steps is None else max(1, st // n)
                print_progress_line(f"Running {len(temps)} MD trajectory(ies) with the model...", 0, 1)
                c = md_candidates(base_atoms, calcs[0], st, temps, stride, sd, info=info,
                                  ensemble=args.md_ensemble, pressure_gpa=args.md_pressure)
                finish_progress_line()
        return c, info

    def screen(cands, label):
        out = []
        with capture_library_noise(library_noise, "screening"):
            for n, cand in enumerate(cands, 1):
                print_progress_line(f"{label} {n}/{len(cands)} with {len(calcs)} model(s)...", n, len(cands))
                ev = evaluate_candidate(cand[0], calcs)
                if ev is not None:
                    ev["de"] = 1000 * (ev["energy"] - e_ref) if e_ref is not None else None
                    ev["too_high"] = bool(e_cap is not None and ev["de"] is not None and ev["de"] > e_cap)
                    out.append((cand, ev))
        finish_progress_line()
        return out

    def raw_score(ev):
        """Comparable across pools (the 'both' score is pool-normalized)."""
        if not committee:
            return ev["max_f"]
        return ev["sigma_e"] if args.score == "energy" else ev["sigma"]

    scale = reference.typical_spacing()
    scale_source = "training set"
    ref_dmin = min(tr_dmin) if tr_dmin else shortest_distance(ase_atoms)
    scan_rows, scan_choice = None, None
    if args.auto_sampling:
        unit = {"rattle": "Ang (stdev)", "md": "K", "strain": "% (strain)"}[args.sampling_method]
        scan_rows = []
        for v in scan_values(args.sampling_method, args.scan_values):
            cands, info = make_pool(v, n=args.scan_n, steps=min(args.md_steps, 20 * args.scan_n),
                                    seed=seed_eff + 7919)
            ev_s = screen(cands, f"Scanning {v:g} {unit.split()[0]}:")
            if committee:
                center_energy_spread([e for _, e in ev_s])
            sc = np.array([raw_score(e) for _, e in ev_s]) if ev_s else np.array([np.nan])
            mf = np.array([e["max_f"] for _, e in ev_s]) if ev_s else np.array([np.nan])
            dm = np.array([shortest_distance(c[0]) for c, _ in ev_s]) if ev_s else np.array([np.nan])
            nv = np.array([reference.novelty(c[0], scale) if scale else np.nan for c, _ in ev_s]) if ev_s else np.array([np.nan])
            nv = np.where(np.isinf(nv), np.nan, nv)
            de_s = [e["de"] for _, e in ev_s if e.get("de") is not None]
            scan_rows.append({"value": v, "n": len(ev_s), "score_med": float(np.nanmedian(sc)),
                              "de_med": float(np.median(de_s)) if de_s else None,
                              "score_max": float(np.nanmax(sc)), "maxf_med": float(np.nanmedian(mf)),
                              "new_pct": (float(100 * np.mean(nv >= args.min_novelty))
                                          if args.min_novelty > 0 and np.isfinite(nv).any() else None),
                              "close_pct": float(100 * np.mean(dm < _CLOSE_CONTACT_FRAC * ref_dmin)),
                              "stopped": any(i_["stopped"] for i_ in info.values())})
        scan_choice, why = choose_scanned(scan_rows, args.max_force, max_energy=args.max_energy,
                                          train_de=train_de90, frontier=args.frontier)
        print_dual(color_text(f"Auto-sampling scan ({args.scan_n} candidates per setting, a fresh seed):", 'bold'), f_out)
        rows = [([f"{r['value']:g}", str(r["n"]), f"{r['score_med']:.3f} / {r['score_max']:.3f}",
                  f"{r['maxf_med']:.2f}",
                  (f"{r['de_med']:.0f} (~{equipartition_t(max(r['de_med'], 0)):.0f} K)" if r["de_med"] is not None
                   else "-"),
                  _fmt(r["new_pct"], "{:.0f}"), f"{r['close_pct']:.0f}",
                  "blew up" if r["stopped"] else "ok"],
                 'green' if r["value"] == scan_choice else ('yellow' if r["stopped"] or r["close_pct"] > 10 else None))
                for r in scan_rows]
        print_table([f"Setting {unit}", "N", "Score med/max", "max|F| med", "E med (meV/atom)", "New %", "Close %",
                     "MD"], rows, f_out)
        print_dual("E med = median predicted energy above the reference; the temperature is that energy as "
                   "3/2 kT per atom (harmonic).", f_out)
        print_dual(f"Chosen            : {scan_choice:g} {unit.split()[0]} -- {why}", f_out)
        if scan_choice == max(r["value"] for r in scan_rows) and len(scan_rows) > 1:
            warn(f"Auto-sampling: the chosen setting ({scan_choice:g}) is the largest one scanned and still physical -- "
                 "larger values may be even more informative; extend the range with --scan-values.")
        if args.sampling_method == "rattle":
            args.stdev = [scan_choice]
        elif args.sampling_method == "md":
            args.md_temperature = [scan_choice]
        else:
            args.strain_max = scan_choice
        print_dual("", f_out)

    t0 = time.monotonic()
    candidates, md_info = make_pool()
    t_gen = time.monotonic() - t0
    n_attempted = (len(args.md_temperature) * (args.md_steps // args.stride) if args.sampling_method == "md"
                   else args.n_candidates)
    t0 = time.monotonic()
    evaluated = screen(candidates, "Screening candidate")
    t_eval = time.monotonic() - t0
    if not evaluated:
        fail("None of the candidate configurations could be evaluated by the model(s).")
    n_failed = len(candidates) - len(evaluated)
    candidates = [c for c, _ in evaluated]
    evaluations = [e for _, e in evaluated]
    n_atoms_c = len(base_atoms)
    print_dual(f"Screened          : {len(candidates)} of {n_attempted} candidate(s) of {n_atoms_c} atoms "
               f"({t_gen:.0f} s generating, {t_eval:.0f} s scoring, "
               f"{1000 * t_eval / max(len(candidates) * len(calcs), 1):.0f} ms per model evaluation)", f_out)
    if n_failed:
        warn(f"Pool: {n_failed} candidate(s) could not be evaluated (a model failed or gave NaN forces) -- left out.")
    if md_info:
        rows = []
        for t, inf in md_info.items():
            rows.append(([f"{t:g}", f"{inf['steps']}/{args.md_steps}", str(inf["frames"]),
                          _fmt(inf["T_mean"], "{:.0f}") + (f" +/- {inf['T_std']:.0f}" if inf["T_std"] else ""),
                          inf["stopped"] or "completed"], 'yellow' if inf["stopped"] else None))
            if inf["stopped"]:
                warn(f"MD: the {t:g} K trajectory with the first model stopped early ({inf['stopped']}) -- the "
                     "model is unstable there; its frames up to that point are in the pool (they are exactly "
                     "what the training set is missing).")
        print_table(["T target (K)", "Steps run", "Frames", "T inst (K)", "Status"], rows, f_out)

    score_mode = apply_score_mode(evaluations, args.score, committee)
    score = np.array([e["score"] for e in evaluations])
    sf_arr = np.array([e["sigma"] for e in evaluations]) if committee else None
    maxf = np.array([e["max_f"] for e in evaluations])
    dmin = np.array([shortest_distance(c[0]) for c in candidates])
    # novelty against the training set alone (before any pick), the measure of new territory
    # (computed even with --min-novelty 0: it is reported either way)
    if scale is None:
        pool = EnvironmentSet()
        for a, _ in candidates:
            pool.add(a)
        scale = pool.typical_spacing()
        scale_source = "candidate pool (training set has < 2 configs)"
    pool_nov = [reference.novelty(c[0], scale) if scale else math.inf for c in candidates]
    nov_arr = np.array([np.nan if math.isinf(v) else v for v in pool_nov])

    groups = sorted({_group(c[1]) for c in candidates}, key=_group_key)
    rows = []
    for g in groups:
        ii = [i for i, c in enumerate(candidates) if _group(c[1]) == g]
        sp = sf_arr if committee else score
        rows.append(([g, str(len(ii)), f"{np.median(maxf[ii]):.2f} / {maxf[ii].max():.2f}",
                      f"{np.median(sp[ii]):.3f} / {sp[ii].max():.3f}",
                      _fmt(np.nanmedian(nov_arr[ii]) if np.isfinite(nov_arr[ii]).any() else None, "{:.2f}"),
                      f"{dmin[ii].min():.3f}"], None))
    print_table(["Setting", "N", "max|F| med/max", ("F spread" if committee else "Score") + " med/max",
                 "Novelty med", "d_min"], rows, f_out)
    print_dual("Units: eV/Ang, novelty in typical training spacings, d_min = shortest interatomic distance (Ang).",
               f_out)
    s_unit = {"force": "eV/Ang", "energy": "eV/atom", "both": "(x pool median)"}[score_mode]
    pct = np.percentile(score, [50, 90])
    print_dual(f"{'Score (' + score_mode + ')':<18}: median {pct[0]:.3f}, 90th percentile {pct[1]:.3f}, max {score.max():.3f} {s_unit}",
               f_out)
    if committee and calib:
        p_f = np.percentile(sf_arr, [50, 90])
        print_dual(f"Estimated error   : force spread x {calib:.1f} (7.1.3 calibration) = median {1000 * calib * p_f[0]:.0f}, "
                   f"90th percentile {1000 * calib * p_f[1]:.0f}, max {1000 * calib * sf_arr.max():.0f} meV/Ang", f_out)
    if committee:
        se_arr = np.array([e["sigma_e"] for e in evaluations])
        ce = _corr(se_arr, sf_arr)
        e_only = int(np.sum((se_arr >= np.percentile(se_arr, 90)) & (sf_arr < np.median(sf_arr))))
        print_dual(f"Energy spread     : median {1000 * np.median(se_arr):.2f}, max {1000 * se_arr.max():.2f} meV/atom; "
                   f"corr with the force spread {_fmt(ce, '{:+.2f}')}; {e_only} candidate(s) in the top 10 % by energy "
                   "spread but below the median force spread"
                   + (" -- --score energy/both would add them" if e_only and score_mode == "force" else ""), f_out)
        c = _corr(sf_arr, maxf)
        print_dual(f"Spread vs force   : corr {_fmt(c, '{:+.2f}')}"
                   + (" -- the spread mostly tracks how distorted a candidate is; it ranks about like max|F|"
                      if c is not None and c > 0.9 else
                      " -- the spread carries information beyond the force magnitude" if c is not None else ""),
                   f_out)
        if args.sampling_method == "rattle" and len(args.stdev) == 1 and c is not None and c > 0.9 and score_mode == "force":
            warn("Pool: the committee spread follows the force magnitude almost exactly (corr > 0.9) -- the "
                 "picks will just be the most displaced candidates; sample several --stdev values or use "
                 "--sampling-method md for a more varied pool.")
    new_share = None
    if scale:
        print_dual(f"Novelty unit      : {scale:.4f} Ang (typical spacing of the {scale_source})", f_out)
        new_share = (100 * float(np.mean(nov_arr >= args.min_novelty))
                     if args.min_novelty > 0 and np.isfinite(nov_arr).any() else None)
        if new_share is not None:
            print_dual(f"New territory     : {new_share:.0f} % of the pool has novelty >= {args.min_novelty:g} "
                       "(environments unlike any training config)", f_out)
            if new_share < _SATURATION_PCT:
                warn(f"Pool: only {new_share:.0f} % of the candidates are new territory -- this sampling is close to "
                     "saturated by the training set; raise --stdev / --md-temperature / --strain-max, use a "
                     "--supercell, or let --auto-sampling choose.")
    else:
        warn("Novelty: could not size the novelty unit -- filter disabled.")
    cell_changed = [not np.allclose(np.asarray(c[0].cell), np.asarray(base_atoms.cell), atol=1e-6) for c in candidates]
    if any(cell_changed):
        dv = np.array([100 * (c[0].get_volume() / base_atoms.get_volume() - 1) for c in candidates])
        print_dual(f"Cell              : {sum(cell_changed)} candidate(s) with their own cell, volume change "
                   f"{dv.min():+.1f} to {dv.max():+.1f} % (written with that lattice)", f_out)
    if tr["fmax"]:
        beyond = 100 * float(np.mean(maxf > max(tr["fmax"]) * 1.1))
        print_dual(f"Beyond training   : {beyond:.0f} % of the pool has a predicted max|F| above the training "
                   f"set's largest SIESTA force ({max(tr['fmax']):.2f} eV/Ang)", f_out)
    de_arr = np.array([e["de"] for e in evaluations if e.get("de") is not None])
    if de_arr.size:
        print_dual(f"Energy above ref. : median {np.median(de_arr):.0f}, max {de_arr.max():.0f} meV/atom (model-predicted; "
                   f"median ~{equipartition_t(max(np.median(de_arr), 0)):.0f} K as 3/2 kT)"
                   + (f"; {int(np.sum(de_arr > e_cap))} above the cap ({e_cap:.0f})" if e_cap is not None else ""),
                   f_out)

    # ---- [3] selection
    print_section("[3] SELECTION", f_out)
    def fresh_reference():
        ref = EnvironmentSet()
        for a in train_atoms:
            ref.add(a)
        return ref

    def mutual_novelty(pk):
        out = []
        for a_i, (i, _, _) in enumerate(pk):
            others = EnvironmentSet()
            for b_i, (j, _, _) in enumerate(pk):
                if b_i != a_i:
                    others.add(candidates[j][0])
            out.append(others.novelty(candidates[i][0], scale))
        return out

    n_energy = sum(1 for e in evaluations if e.get("too_high")
                   and not (args.max_force and e["max_f"] > args.max_force))
    nov_scale = scale if args.min_novelty > 0 or args.selection == "diverse" else None
    g_picked, g_force, g_red = select_candidates(
        candidates, evaluations, args.top_k, fresh_reference(), scale if args.min_novelty > 0 else None,
        args.min_novelty, args.max_force)
    if args.selection == "diverse":
        picked, n_force, n_redundant = select_diverse(
            candidates, evaluations, args.top_k, fresh_reference(), nov_scale, args.min_novelty,
            args.max_force, args.diversity)
        examined = min(len(candidates) - n_force - n_energy, max(20 * args.top_k, 50))
        ex_label = "Considered (highest scores)"
    else:
        picked, n_force, n_redundant = g_picked, g_force, g_red
        examined = len(picked) + n_force + n_energy + n_redundant
        ex_label = "Examined (highest score first)"
    rows = [(["Screened", str(len(candidates))], None),
            ([ex_label, str(examined)], None),
            (["Skipped: max|F| > --max-force", str(n_force)], 'yellow' if n_force else None),
            (["Skipped: E above the energy cap", str(n_energy)], 'yellow' if n_energy else None),
            (["Skipped: novelty < --min-novelty", str(n_redundant)], 'yellow' if n_redundant else None),
            (["Picked", f"{len(picked)} of {args.top_k} requested"], 'green')]
    print_table(["Selection funnel", "Count"], rows, f_out)
    if n_energy and len(picked) < args.top_k:
        warn(f"Selection: the energy cap ({e_cap:.0f} meV/atom above the reference, {e_cap_why}) left only "
             f"{len(picked)} of {args.top_k} picks -- the pool is sampled harder than the training set reaches; "
             "use a smaller --stdev / --md-temperature (or let --auto-sampling choose), or raise --frontier / "
             "--max-energy (0 = off) if such geometries matter.")
    sel_compare = None
    if picked:
        ps = [ev["score"] for _, ev, _ in picked]
        rank_last = int(np.sum(score >= min(ps)))
        print_dual(f"Score cut         : picks have score {min(ps):.3f}-{max(ps):.3f} {s_unit}, i.e. within the top "
                   f"{rank_last} of {len(candidates)} by score ({100 * rank_last / len(candidates):.0f} %)", f_out)
        if args.selection == "diverse" and g_picked and scale:
            same = len({i for i, _, _ in picked} & {i for i, _, _ in g_picked})
            md_div = min(mutual_novelty(picked)) if len(picked) > 1 else None
            md_gr = min(mutual_novelty(g_picked)) if len(g_picked) > 1 else None
            ratio = np.mean([ev["score"] for _, ev, _ in picked]) / np.mean([ev["score"] for _, ev, _ in g_picked])
            sel_compare = {"same": same, "score_ratio": float(ratio), "min_mutual_diverse": md_div,
                           "min_mutual_greedy": md_gr}
            print_dual(f"vs greedy         : {len(picked) - same} of {len(picked)} pick(s) differ from score-only; mean "
                       f"score {100 * ratio:.0f} % of greedy's; closest pair novelty {_fmt(md_div, '{:.2f}')} (greedy "
                       f"{_fmt(md_gr, '{:.2f}')})", f_out)
    if len(picked) < args.top_k:
        warn(f"Selection: only {len(picked)} of the requested {args.top_k} candidate(s) passed the filters -- the "
             "training set already covers this sampling; try larger --stdev, higher --md-temperature, a "
             "--supercell, or a lower --min-novelty.")
    if n_redundant and n_redundant > 2 * max(len(picked), 1):
        print_dual(f"(The novelty filter rejected {n_redundant} high-score candidates: they resemble the training set "
                   "or each other.)", f_out)

    # ---- [4] the picks
    print_section("[4] PICKED CONFIGURATIONS", f_out)
    order = list(np.argsort(-score))
    rows, pick_rows, close = [], [], []
    for n, (i, ev, nov) in enumerate(picked, 1):
        atoms, origin = candidates[i]
        flags = []
        if tr["fmax"] and ev["max_f"] > max(tr["fmax"]) * 1.1:
            flags.append("F")
        if tr_dmin and dmin[i] < min(tr_dmin) * 0.97:
            flags.append("d")
        if train_de is not None and ev.get("de") is not None and ev["de"] > max(train_de, 0) * 1.1 + 1.0:
            flags.append("E")
        if tr_dmin and dmin[i] < min(tr_dmin) * _CLOSE_CONTACT_FRAC:
            close.append(f"{n} ({dmin[i]:.2f} Ang)")
        worst = (int(np.argmax(ev["spread_atom"])) if ev["spread_atom"] is not None else int(np.argmax(ev["f_atom"])))
        est = f"{1000 * calib * ev['sigma']:.0f}" if calib and committee else "-"
        pick_rows.append((i, ev, nov, flags, worst))
        rows.append(([str(n), origin, f"{ev['score']:.3f}", est,
                      _fmt(1000 * ev["sigma_e"] if ev["sigma_e"] is not None else None, "{:.2f}"), f"{ev['max_f']:.2f}",
                      _fmt(None if math.isinf(nov) else nov, "{:.2f}"), f"{dmin[i]:.3f}",
                      f"{order.index(i) + 1}", f"{atoms.get_chemical_symbols()[worst]}{worst}",
                      "".join(flags) or "in"], 'yellow' if flags else None))
    if rows:
        print_table(["#", "Origin", "Score", "Est. err", "E spr", "max|F|", "Novelty", "d_min", "Rank", "Worst atom",
                     "Range"], rows, f_out)
        print_dual(f"Score {s_unit}, max|F| eV/Ang; Est. err = calibration x force spread (meV/Ang, from 7.1.3); E spr = "
                   "energy spread (meV/atom); novelty vs training set + earlier picks;", f_out)
        print_dual("Rank by score in the pool; worst atom = largest " + ("spread" if committee else "force")
                   + "; Range: in, or beyond the training set in F (max|F|) / d (shortest distance) / E (energy).", f_out)
        pde = [ev["de"] for _, ev, _ in picked if ev.get("de") is not None]
        if pde:
            print_dual(f"Energy above ref. : picks {min(pde):.0f} .. {max(pde):.0f} meV/atom (~{equipartition_t(max(min(pde), 0)):.0f}"
                       f" .. {equipartition_t(max(max(pde), 0)):.0f} K as 3/2 kT)"
                       + (f"; training set up to {train_de:.0f}" if train_de is not None else ""), f_out)
            if train_de is not None and train_de > 0 and np.median(pde) > 2 * train_de:
                warn(f"Picks: their median energy ({np.median(pde):.0f} meV/atom above the reference, "
                     f"~{equipartition_t(np.median(pde)):.0f} K) is more than twice the training set's highest "
                     f"({train_de:.0f}) -- a jump far from the labelled region; such configs can make the model WORSE "
                     "near equilibrium (check with 7.1.3). A smaller --stdev / --md-temperature, or --max-energy, "
                     "keeps the next round closer.")
        if close:
            warn(f"Picks: {len(close)} pick(s) have atoms closer than {100 * _CLOSE_CONTACT_FRAC:.0f} % of the training "
                 f"set's shortest distance ({min(tr_dmin):.2f} Ang): pick {', '.join(close)} -- SIESTA may struggle to "
                 "converge them and they may be irrelevant for the intended use; --max-force or a smaller --stdev "
                 "keeps such geometries out.")
        n_out = sum(1 for r in pick_rows if r[3])
        if n_out:
            print_dual(f"Extrapolation     : {n_out} of {len(picked)} pick(s) lie beyond the training set -- that is "
                       "usually the point of active learning", f_out)
        if len(picked) > 1 and scale:
            pair = mutual_novelty(picked)
            print_dual(f"Diversity         : each pick vs the other picks, novelty min {min(pair):.2f}, median "
                       f"{np.median(pair):.2f}" + (" -- some picks are near-duplicates of each other"
                                                   if min(pair) < 0.5 * max(args.min_novelty, 0.1) else ""), f_out)
        by_group = {}
        for i, _, _ in picked:
            by_group[_group(candidates[i][1])] = by_group.get(_group(candidates[i][1]), 0) + 1
        print_dual("From              : " + ", ".join(f"{g} {n}" for g, n in sorted(by_group.items(), key=lambda kv: _group_key(kv[0]))), f_out)
    else:
        print_dual("No pick.", f_out)

    # ---- [5] is another round worth it?
    print_section("[5] IS ANOTHER ROUND WORTH IT?", f_out)
    stop_info = None
    if committee:
        est = 1000 * sf_arr * (calib or 1.0)
        kind = "estimated error (calibrated)" if calib else "force spread (a LOWER bound of the error: no calibration)"
        p90, top = float(np.percentile(est, 90)), float(est.max())
        print_dual(f"Pool uncertainty  : {kind}: 90th percentile {p90:.0f}, max {top:.0f} meV/Ang", f_out)
        val = (round_entry or {}).get("validation") or {} if models_from_round else {}
        if val.get("overall"):
            print_dual(f"Test set (7.1.3) : F RMSE {val['overall']['F_rmse_meV']:.1f} meV/Ang"
                       + (f", verdict {val['verdict'].upper()}" if val.get("verdict") else ""), f_out)
        if args.target_force is not None:
            n_above = int(np.sum(est > args.target_force))
            if top <= args.target_force and calib:
                code, msg = "stop", (f"even the most uncertain candidate is expected within {args.target_force:g} meV/Ang "
                                     "-- another round on this sampling adds little; change the sampling or stop.")
            elif top <= args.target_force:
                code, msg = "unsure", (f"every spread is below {args.target_force:g} meV/Ang, but without a calibration the "
                                       "spread underestimates the error -- run stb-mlffValidate on this round first.")
            elif p90 <= args.target_force:
                code, msg = "small", (f"only the tail is above the target: {n_above} candidate(s) expected above "
                                      f"{args.target_force:g} meV/Ang -- a small round (--top-k {max(1, min(n_above, args.top_k))}) "
                                      "is enough.")
            else:
                code, msg = "continue", (f"{n_above} of {len(est)} candidates ({100 * n_above / len(est):.0f} %) are expected "
                                         f"above {args.target_force:g} meV/Ang -- another round is worth it.")
            print_dual(color_text(f"[{code.upper()}] {msg}", 'green' if code == "stop" else 'yellow'), f_out)
            stop_info = {"target_meV": args.target_force, "verdict": code, "n_above": n_above,
                         "p90_meV": round(p90, 1), "max_meV": round(top, 1), "calibrated": bool(calib)}
        else:
            print_dual("(Give --target-force, the error good enough for your use, for a verdict.)", f_out)
            stop_info = {"p90_meV": round(p90, 1), "max_meV": round(top, 1), "calibrated": bool(calib)}
    else:
        print_dual("Needs a committee: a single model gives no uncertainty to compare with a target.", f_out)

    # ---- [6] history across rounds
    print_section("[6] ACTIVE-LEARNING HISTORY", f_out)
    hist = []
    for r in (state or {}).get("rounds", []):
        c_r = ((r.get("validation") or {}).get("calibration") or {}).get("error_over_spread")
        for n_run, al in enumerate(r.get("active_learning", []), 1):
            pool = al.get("pool") or {}
            hist.append({"label": f"{r['round']}.{n_run}", "round": r["round"], "sampling": al.get("sampling", "?"),
                         "settings": al.get("settings") or {}, "score": al.get("score", "force"),
                         "n": pool.get("n"), "med": pool.get("score_median"), "max": pool.get("score_max"),
                         "err_max": (1000 * c_r * pool["force_spread_max"] if c_r and pool.get("force_spread_max")
                                     else None),
                         "new": pool.get("new_share_pct"), "picked": len(al.get("picked", [])),
                         "now": False})
    cur_settings = ({"n_candidates": args.n_candidates, "stdev": args.stdev} if args.sampling_method == "rattle" else
                    {"n_candidates": args.n_candidates, "strain_max": args.strain_max, "stdev": args.stdev}
                    if args.sampling_method == "strain" else
                    {"md_steps": args.md_steps, "md_temperature": args.md_temperature, "stride": args.stride,
                     "ensemble": args.md_ensemble})
    hist.append({"label": f"{proposer_round}.{n_al_runs + 1}" + (" (now)" if not args.dry_run else " (dry run)"),
                 "round": proposer_round, "sampling": args.sampling_method, "score": score_mode,
                 "settings": cur_settings, "n": len(candidates), "med": float(np.median(score)), "max": float(score.max()),
                 "err_max": 1000 * calib * float(sf_arr.max()) if (calib and committee) else None,
                 "new": new_share, "picked": len(picked), "now": True})
    if len(hist) > 1:
        rows = [([h["label"], h["sampling"], h["score"], _fmt(h["n"], "{}"), _fmt(h["med"]), _fmt(h["max"]),
                  _fmt(h["err_max"], "{:.0f}"), _fmt(h["new"], "{:.0f}"), str(h["picked"])],
                 'green' if h["now"] else None) for h in hist]
        print_table(["Round.run", "Sampling", "Score", "Pool", "Score med", "Score max", "Est. err max", "New %",
                     "Picked"], rows, f_out)
        def _key(st):   # what makes two pools comparable (the pool size aside)
            return {k: v for k, v in (st or {}).items() if k not in ("n_candidates",)}
        same = [h for h in hist[:-1] if h["sampling"] == args.sampling_method and h["score"] == score_mode
                and h["max"] is not None and _key(h["settings"]) == _key(cur_settings)]
        if same:
            prev = same[-1]
            ch = 100 * (hist[-1]["max"] - prev["max"]) / prev["max"] if prev["max"] else None
            print_dual(f"Trend             : max score {prev['max']:.3f} (run {prev['label']}) -> {hist[-1]['max']:.3f} now "
                       f"({_fmt(ch, '{:+.0f} %')})" + (" -- the committee got more confident on this sampling"
                                                       if ch is not None and ch < -20 else
                                                       " -- not more confident yet" if ch is not None else ""), f_out)
        if not same:
            print_dual("Trend             : no earlier run with the same sampling, settings and score mode to compare with",
                       f_out)
        print_dual("Scores compare only within the same sampling, settings and score mode; Est. err uses each "
                   "round's own 7.1.3 calibration.", f_out)
    else:
        print_dual("First active-learning run of this workflow -- the next runs will be compared with it here.", f_out)

    # ---- [7] cost estimate
    print_section("[7] SIESTA COST ESTIMATE", f_out)
    if picked:
        print_dual(f"New calculations  : {len(picked)} single point(s) of {n_atoms_c} atoms", f_out)
        if tr["wall"] and tr["n_atoms"]:
            w = float(np.median(tr["wall"]))
            n_ref = int(np.median(tr["n_atoms"]))
            each = w * (n_atoms_c / n_ref) ** 3
            print_dual(f"Estimate          : ~{each:.0f} s each, ~{each * len(picked) / 60:.1f} min in total (median "
                       f"training run {w:.0f} s for {n_ref} atoms, scaled as N^3; same cores)", f_out)
            if each > 3600:
                warn(f"Cost: each new config takes an estimated {each / 3600:.1f} h -- a smaller --supercell or "
                     "fewer --top-k keeps the round affordable.")
        else:
            print_dual("Estimate          : no finished training run with a timing to scale from.", f_out)
    else:
        print_dual("Nothing to run.", f_out)

    # ---- write
    folders, picked_entries = [], []
    if picked and not args.dry_run:
        os.makedirs(picks_root, exist_ok=True)
        start_index = next_config_index(args.output_dir)
        for offset, (i, ev, nov) in enumerate(picked):
            atoms, origin = candidates[i]
            atoms = atoms.copy()
            atoms.wrap()
            positions = (atoms.get_scaled_positions() if ref_structure.coord_format == "fractional"
                         else atoms.get_positions())
            raw_lat = np.asarray(atoms.cell) / template_lc if cell_changed[i] else None
            folder = write_config_folder(
                picks_root, start_index + offset, template, ref_basename, positions, raw_lat,
                pseudos_to_copy, extra_block=extra_block, calc_template=calc_template)
            folders.append(folder)
            picked_entries.append({"folder": folder, "set": "train",
                                   "kind": f"al-{args.sampling_method}",
                                   "origin": origin + (f", supercell {dims[0]}x{dims[1]}x{dims[2]}" if is_supercell else ""),
                                   "stage": "active-learning", "round": proposer_round,
                                   "n_atoms": len(atoms), "score": ev["score"], "max_f": ev["max_f"],
                                   "novelty": None if math.isinf(nov) else nov,
                                   "d_min": round(float(dmin[i]), 4)})
    if tmp_dir:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        tmp_dir = None

    if state is not None and not args.dry_run:
        mlff_workflow.add_configs(state, picked_entries, wf_base)
        target = mlff_workflow.get_round(state, proposer_round) if proposer_round else None
        if target is not None:
            target.setdefault("active_learning", []).append({
                "time": datetime.now().isoformat(timespec="seconds"),
                "sampling": args.sampling_method, "supercell": dims,
                "settings": ({"n_candidates": args.n_candidates, "stdev": args.stdev}
                             if args.sampling_method == "rattle" else
                             {"n_candidates": args.n_candidates, "strain_max": args.strain_max, "stdev": args.stdev}
                             if args.sampling_method == "strain" else
                             {"md_steps": args.md_steps, "md_temperature": args.md_temperature,
                              "stride": args.stride, "ensemble": args.md_ensemble}),
                "score": score_mode, "selection": args.selection, "seed": seed_eff, "seed_base": args.seed,
                "auto_sampling": ({"rows": scan_rows, "chosen": scan_choice} if scan_rows else None),
                "selection_vs_greedy": sel_compare, "another_round": stop_info,
                "committee_size": len(args.model),
                "pool": {"n": len(candidates), "score_median": float(np.median(score)),
                         "score_max": float(score.max()),
                         "force_spread_max": float(sf_arr.max()) if committee else None,
                         "energy_spread_max": float(max(e["sigma_e"] for e in evaluations)) if committee else None,
                         "new_share_pct": (float(100 * np.mean(nov_arr >= args.min_novelty))
                                           if args.min_novelty > 0 and np.isfinite(nov_arr).any() else None)},
                "skipped": {"max_force": n_force, "max_energy": n_energy, "novelty": n_redundant},
                "picked": [mlff_workflow.rel(f, wf_base) for f in folders]})
        mlff_workflow.record(state, "active-learning", "stb-mlffActiveLearning", sys.argv[1:])
        mlff_workflow.save(state, args.workflow)

    show_requested = args.show
    if args.show and not can_show():
        warn("Plots: --show needs a display (no interactive matplotlib backend here) -- the figure is only saved.")
        args.show = False
    png = near_out(f"mlff_al_{tag}.png")
    try:
        plot_al(png, candidates, evaluations, pool_nov, picked, committee, args.min_novelty, args.max_force,
                md_info, args.show, score_mode,
                extras={"scan": (scan_rows, scan_choice, {"rattle": "stdev, Ang", "md": "K",
                                                          "strain": "strain, %"}[args.sampling_method])
                        if scan_rows else None,
                        "energy": (sf_arr, np.array([e["sigma_e"] for e in evaluations])) if committee else None,
                        "history": hist if len(hist) > 1 else None})
    except Exception as e:
        warn(f"Plots: screening figure not written ({e}).")
        png = None

    data_path = None
    if args.save_data:
        data_path = near_out(f"mlff_al_{tag}_pool.dat")
        picked_set = {i for i, _, _ in picked}
        with open(data_path, "w") as f:
            f.write("# Candidate pool, stb-mlffActiveLearning\n# index origin score(eV/Ang) spread(eV/Ang) "
                    "max_F(eV/Ang) novelty d_min(Ang) E(eV/atom) picked\n")
            for i, ((_, origin), ev) in enumerate(zip(candidates, evaluations)):
                nv = pool_nov[i]
                f.write(f"{i} \"{origin}\" {ev['score']:.5f} {_fmt(ev['sigma'], '{:.5f}')} {ev['max_f']:.5f} "
                        f"{'inf' if math.isinf(nv) else f'{nv:.4f}'} {dmin[i]:.4f} {ev['energy']:.6f} "
                        f"{int(i in picked_set)}\n")

    picks_path = None
    if args.view_picks and picked:
        from ase.io import write as ase_write
        frames = []
        for n, (i, ev, nov) in enumerate(picked, 1):
            b = candidates[i][0].copy()
            b.calc = None
            per_atom = ev["spread_atom"] if ev["spread_atom"] is not None else ev["f_atom"]
            b.new_array("spread" if committee else "max_F", per_atom)
            b.set_initial_charges(per_atom)
            b.info = {"pick": n, "origin": candidates[i][1], "score": round(ev["score"], 4),
                      "folder": os.path.basename(folders[n - 1]) if folders else ""}
            frames.append(b)
        picks_path = near_out(f"mlff_al_{tag}_picks.xyz")
        ase_write(picks_path, frames, format="extxyz")

    # ---- [8] files
    if folders:   # before the status line counts the warnings
        warn(f"SIESTA: run it in {os.path.relpath(picks_root)}/ ({len(folders)} new folder(s)) before the next "
             "stb-mlffAnalysis -- unlabelled picks are left out of training.")
    print_section("[8] SUMMARY & FILES", f_out)
    print_dual(f"Status            : OK, {len(picked)} pick(s) from {len(candidates)} candidate(s) in "
               f"{time.monotonic() - t_start:.0f} s" + (f" ({len(warns)} warning(s), see [10])" if warns else ""), f_out)
    if args.dry_run:
        print_dual(color_text("Dry run           : no folder written, workflow file untouched.", 'cyan'), f_out)
    for n, fd in enumerate(folders, 1):
        print_dual(f"Pick {n:<13}: {fd}", f_out)
    if state is not None and not args.dry_run:
        print_dual(f"Workflow file     : {args.workflow} ({len(folders)} config(s) added to the training set, "
                   f"proposed by round {proposer_round})", f_out)
    if png:
        print_dual(f"Figure            : {png}", f_out)
    if data_path:
        print_dual(f"Pool data         : {data_path}", f_out)
    if picks_path:
        print_dual(f"Picks (ASE)       : {picks_path} (color by initial charge = per-atom "
                   f"{'spread' if committee else 'force'}, eV/Ang)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    # ---- [7] next steps
    print_section("[9] NEXT STEPS", f_out)
    if folders:
        main_fdf = os.path.basename(calc_template) if calc_template else ref_basename
        label = structure_io.find_fdf_system_label(folders[0]) or "siesta"
        names = " ".join(os.path.relpath(f) for f in folders)
        print_dual(color_text(f"Run SIESTA in     : {os.path.relpath(picks_root)}/ ({len(folders)} new folder(s): "
                              f"{', '.join(os.path.basename(f) for f in folders)})", 'cyan'), f_out)
        print_dual("1. Run SIESTA in each new folder, e.g.:", f_out)
        print_dual(f'     for d in {names}; do (cd "$d" && siesta < {main_fdf} > {label}.out); done', f_out)
        print_dual("2. Train the next round: stb-mlffAnalysis (7.1.2) -- it takes the configs from the workflow "
                   "file; --init-from previous is faster.", f_out)
        print_dual("3. Score it: stb-mlffValidate (7.1.3); repeat 7.1.4 while the verdict says CONTINUE.", f_out)
    elif args.dry_run:
        print_dual("1. Happy with the picks? Run again without --dry-run to write them.", f_out)
    else:
        print_dual("1. Nothing picked: change the sampling (--stdev, --md-temperature, --supercell) or lower "
                   "--min-novelty.", f_out)
    extras = [t for t, used in (("--dry-run (screen only)", args.dry_run),
                                ("--show (display the figure)", show_requested),
                                ("--save-data (the whole pool)", args.save_data),
                                ("--view-picks (picks in ASE's viewer, colored by uncertainty)", args.view_picks),
                                ("--supercell NA NB NC (longer-range environments)", is_supercell),
                                ("--sampling-method md", args.sampling_method == "md"))
              if not used]
    if extras:
        print_dual(f"Also available    : {'; '.join(extras)}.", f_out)

    # ---- [10] warnings
    print_section("[10] WARNINGS", f_out)
    if warns:
        for k, w in enumerate(warns, 1):
            print_dual(color_text(f"{k:>2}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)
    if library_noise:
        print_dual(f"Library messages  : {len(library_noise)} informational message(s) from mace/torch captured "
                   "(model loading, dtype).", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text(("Screening done (dry run).\n" if args.dry_run else
                      "Candidates suggested -- label them, then retrain.\n"), 'bold'))
    if picks_path:
        from stb.core.ase_view import view_structure_interactive
        from ase.io import read as ase_read
        print(color_text(f"[INFO] Opening the {len(picked)} pick(s) in ASE's viewer (View > Colors > By initial "
                         f"charge = per-atom {'committee spread' if committee else 'force'}, eV/Ang).", 'cyan'))
        view_structure_interactive(ase_read(picks_path, ":"))


if __name__ == "__main__":
    main()
