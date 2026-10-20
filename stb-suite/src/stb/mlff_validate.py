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
import argparse
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from stb.core import mlff_workflow
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.mlff_analysis import collect_training_configs, _EV_A3_TO_GPA

REPORT_FILE = "stb_mlffValidate_report.txt"
DEFAULT_PLATEAU_PCT = 5.0
# Above this ratio of true force error to committee spread the spread is
# reported as uncalibrated: still usable to RANK candidates, not as an error bar.
_CALIBRATION_WARN_RATIO = 3.0


def predict_committee(calcs, atoms):
    """Per-model energy (eV), forces (n, 3) and stress (3x3 or None)."""
    out = []
    for calc in calcs:
        a = atoms.copy()
        a.calc = calc
        e = a.get_potential_energy()
        f = a.get_forces()
        try:
            s = a.get_stress(voigt=False) if "REF_stress" in atoms.info else None
        except Exception:
            s = None
        out.append((e, f, s))
    return out


def evaluate_round(calcs, test_atoms, groups):
    """Committee-mean and per-member errors on the test set, overall and per
    group, plus the spread/error calibration. Returns a JSON-ready dict."""
    n_models = len(calcs)
    rows = []   # per config
    member_err = [{"de": [], "df": [], "ds": []} for _ in range(n_models)]
    for atoms, group in zip(test_atoms, groups):
        preds = predict_committee(calcs, atoms)
        n = len(atoms)
        e_ref = atoms.info["REF_energy"] / n
        f_ref = atoms.arrays["REF_forces"]
        s_ref = np.asarray(atoms.info["REF_stress"]) if "REF_stress" in atoms.info else None
        es = np.array([p[0] / n for p in preds])
        fs = np.array([p[1] for p in preds])
        f_mean = fs.mean(axis=0)
        s_mean = (np.mean([p[2] for p in preds], axis=0)
                  if s_ref is not None and all(p[2] is not None for p in preds) else None)
        for m, (e, f, s) in enumerate(preds):
            member_err[m]["de"].append(e / n - e_ref)
            member_err[m]["df"].append((f - f_ref).ravel())
            if s is not None and s_ref is not None:
                member_err[m]["ds"].append((np.asarray(s) - s_ref).ravel())
        err_atom = np.linalg.norm(f_mean - f_ref, axis=1)
        spread_atom = (np.sqrt(np.mean(np.sum((fs - f_mean) ** 2, axis=2), axis=0))
                       if n_models > 1 else None)
        rows.append({"group": group, "de": float(es.mean() - e_ref),
                     "df": (f_mean - f_ref).ravel(), "f_ref": f_ref.ravel(), "f_pred": f_mean.ravel(),
                     "ds": (s_mean - s_ref).ravel() if s_mean is not None else None,
                     "err_max": float(err_atom.max()),
                     "spread_max": float(spread_atom.max()) if spread_atom is not None else None,
                     "err_atom": err_atom, "spread_atom": spread_atom})

    def summarize(sel):
        if not sel:
            return None
        de = np.array([r["de"] for r in sel])
        df = np.concatenate([r["df"] for r in sel])
        ds = [r["ds"] for r in sel if r["ds"] is not None]
        return {"n": len(sel),
                "E_rmse_meV": float(1000 * np.sqrt(np.mean(de ** 2))),
                "E_rmse_rel_meV": float(1000 * np.sqrt(np.mean((de - de.mean()) ** 2))),
                "F_rmse_meV": float(1000 * np.sqrt(np.mean(df ** 2))),
                "F_mae_meV": float(1000 * np.mean(np.abs(df))),
                "S_rmse_GPa": (float(np.sqrt(np.mean(np.concatenate(ds) ** 2)) * _EV_A3_TO_GPA)
                               if ds else None)}

    result = {"overall": summarize(rows),
              "groups": {g: summarize([r for r in rows if r["group"] == g])
                         for g in sorted(set(groups))},
              "members": []}
    for m in range(n_models):
        de = np.array(member_err[m]["de"])
        df = np.concatenate(member_err[m]["df"])
        ds = member_err[m]["ds"]
        result["members"].append({
            "E_rmse_meV": float(1000 * np.sqrt(np.mean(de ** 2))),
            "F_rmse_meV": float(1000 * np.sqrt(np.mean(df ** 2))),
            "S_rmse_GPa": (float(np.sqrt(np.mean(np.concatenate(ds) ** 2)) * _EV_A3_TO_GPA)
                           if ds else None)})
    if n_models > 1:
        spread = np.array([r["spread_max"] for r in rows])
        err = np.array([r["err_max"] for r in rows])
        sa = np.concatenate([r["spread_atom"] for r in rows])
        ea = np.concatenate([r["err_atom"] for r in rows])
        result["calibration"] = {
            "config_corr": _corr(spread, err),
            "atom_corr": _corr(sa, ea),
            "error_over_spread": float(np.median(err / np.maximum(spread, 1e-12))),
            "spread_mean_meV": float(1000 * spread.mean())}
    result["_parity"] = (np.concatenate([r["f_ref"] for r in rows]),
                         np.concatenate([r["f_pred"] for r in rows]))
    return result


def _corr(a, b):
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def foundation_baseline(model_size, device, test_atoms):
    """Raw foundation model on the test set: offset-removed E RMSE (its absolute
    energies are on another scale), F RMSE, S RMSE."""
    from stb.core import mace_relax
    calc = mace_relax.get_calculator(model=model_size, device=device)
    de, df, ds = [], [], []
    for atoms in test_atoms:
        a = atoms.copy()
        a.calc = calc
        de.append(a.get_potential_energy() / len(a) - atoms.info["REF_energy"] / len(a))
        df.append((a.get_forces() - atoms.arrays["REF_forces"]).ravel())
        if "REF_stress" in atoms.info:
            ds.append((a.get_stress(voigt=False) - np.asarray(atoms.info["REF_stress"])).ravel())
    de = np.array(de)
    return {"E_rmse_rel_meV": float(1000 * np.sqrt(np.mean((de - de.mean()) ** 2))),
            "F_rmse_meV": float(1000 * np.sqrt(np.mean(np.concatenate(df) ** 2))),
            "S_rmse_GPa": float(np.sqrt(np.mean(np.concatenate(ds) ** 2)) * _EV_A3_TO_GPA) if ds else None}


# --------------------------------------------------------------------------
# Phonons against a finished stb-phononsCreate run
# --------------------------------------------------------------------------

def load_siesta_phonons(phonon_dir, f_out):
    """(SIESTA phonopy object with force constants, internal_to_angstrom),
    from a finished stb-phononsCreate folder -- the same loader stb-phononsPos
    uses. Restores the working directory."""
    from stb.core.phonon_workflow import load_phonon_with_force_constants, detect_system_label
    label = detect_system_label(phonon_dir) or "siesta"
    phonon, to_ang, original = load_phonon_with_force_constants(phonon_dir, label, False, f_out)
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
    uc = phonopy_atoms_to_ase(ref_phonon.unitcell, to_ang)
    unitcell = PhonopyAtoms(symbols=uc.get_chemical_symbols(), cell=uc.cell[:],
                            scaled_positions=uc.get_scaled_positions())
    ml = Phonopy(unitcell, supercell_matrix=ref_phonon.supercell_matrix,
                 primitive_matrix=ref_phonon.primitive_matrix)
    dataset = {"natom": ref_phonon.dataset["natom"],
               "first_atoms": [{"number": d["number"],
                                "displacement": (np.array(d["displacement"]) * to_ang).tolist()}
                               for d in ref_phonon.dataset["first_atoms"]]}
    ml.dataset = dataset
    forces = []
    for sc in ml.supercells_with_displacements:
        from ase import Atoms
        atoms = Atoms(symbols=sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        fs = []
        for calc in calcs:
            a = atoms.copy()
            a.calc = calc
            fs.append(a.get_forces())
        forces.append(np.mean(fs, axis=0))
    ml.forces = np.array(forces)
    ml.produce_force_constants()
    ml.symmetrize_force_constants(show_drift=False)
    return ml


def phonon_errors(ref_phonon, ml_phonon, mesh):
    ref_phonon.run_mesh(mesh)
    q = ref_phonon.mesh.qpoints
    f_ref = ref_phonon.mesh.frequencies
    ml_phonon.run_qpoints(q)
    f_ml = ml_phonon.get_qpoints_dict()["frequencies"]
    err = np.sort(f_ml, axis=1) - np.sort(f_ref, axis=1)
    n_ac = 3
    return {"mesh": list(mesh), "rmse_THz": float(np.sqrt(np.mean(err ** 2))),
            "acoustic_rmse_THz": float(np.sqrt(np.mean(err[:, :n_ac] ** 2))),
            "optical_rmse_THz": (float(np.sqrt(np.mean(err[:, n_ac:] ** 2)))
                                 if err.shape[1] > n_ac else None),
            "max_abs_THz": float(np.abs(err).max()),
            "ref_max_THz": float(f_ref.max()), "model_max_THz": float(f_ml.max())}


# --------------------------------------------------------------------------

def verdict(history, target_force, plateau_pct):
    """(code, message) from the F RMSE of the scored rounds, oldest first."""
    latest = history[-1][1]
    if target_force is not None and latest <= target_force:
        return "stop", (f"Target met: test F RMSE {latest:.1f} <= {target_force:g} meV/Ang -- "
                        "stop here and use the model.")
    if len(history) < 2:
        return "continue", ("First scored round -- run an active-learning round (4.17.4), "
                            "label it, train (4.17.2) and validate again to see the trend.")
    prev = history[-2][1]
    gain = 100 * (prev - latest) / prev
    if gain < plateau_pct:
        return "plateau", (f"Round {history[-1][0]} changed the test F RMSE by {gain:+.1f} % "
                           f"(< {plateau_pct:g} %): more configs of the SAME kind are unlikely to help. "
                           "Change the sampling (larger --stdev, higher --md-temperature, "
                           "--supercell) or stop.")
    return "continue", (f"Round {history[-1][0]} lowered the test F RMSE by {gain:.1f} % -- "
                        "another active-learning round is likely to help.")


def plot_validation(scored, out_path):
    """F RMSE per round (overall + per group) and force parity of the latest round."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    rounds = [r for r, _ in scored]
    ax = axes[0]
    ax.plot(rounds, [v["overall"]["F_rmse_meV"] for _, v in scored], "o-", color="k", lw=2, label="all")
    groups = sorted({g for _, v in scored for g in v["groups"]})
    for g in groups:
        ys = [v["groups"].get(g, {}) or {} for _, v in scored]
        ax.plot(rounds, [y.get("F_rmse_meV", np.nan) for y in ys], "o--", label=g)
    ax.set_xlabel("Round")
    ax.set_ylabel("Test F RMSE (meV/Ang)")
    ax.set_xticks(rounds)
    ax.set_title("Force error on the test set")
    ax.legend(fontsize=8)
    f_ref, f_pred = scored[-1][1]["_parity"]
    ax = axes[1]
    ax.scatter(f_ref, f_pred, s=6, alpha=0.5)
    lo, hi = min(f_ref.min(), f_pred.min()), max(f_ref.max(), f_pred.max())
    ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=0.8)
    ax.set_xlabel("SIESTA force (eV/Ang)")
    ax.set_ylabel("Committee-mean force (eV/Ang)")
    ax.set_title(f"Force parity, round {scored[-1][0]} (test set)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Stage 3 (Validate) of the Custom ML Force Field workflow: scores every "
                     "training round recorded in mlff_workflow.json (its committee mean and "
                     "each member) against the held-out SIESTA test set from stb-mlff "
                     "--n-test -- energy, force and stress errors, overall and per kind of "
                     "configuration --, checks whether the committee's disagreement tracks "
                     "the real error, optionally compares phonons against a finished "
                     "stb-phononsCreate run, and says whether another active-learning round "
                     "is worth it. Reads finished calculations only; writes no new folders.",
        epilog="Example usage:\n"
               "  stb-mlffValidate                         # every round in mlff_workflow.json\n"
               "  stb-mlffValidate --target-force 30 --phonon-ref phonon_runs\n"
               "  stb-mlffValidate --model a.model b.model --test-path 'mlff_test/test_*'",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--workflow", default=mlff_workflow.WORKFLOW_FILE, metavar="JSON",
                        help=f"Workflow state file (default: {mlff_workflow.WORKFLOW_FILE}). Rounds "
                             "and test set are read from it; each round's metrics are stored back.")
    parser.add_argument("--rounds", type=int, nargs="+", default=None, metavar="K",
                        help="Rounds to score (default: all rounds in --workflow).")
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
                             "round was trained from, else small).")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="Do not score the raw foundation model on the test set.")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"],
                        help="Device to run the models on (default: cpu).")
    parser.add_argument("--png", default="mlff_validation.png", metavar="FILE",
                        help="Summary figure (default: mlff_validation.png, next to --workflow).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
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

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    try:
        state = mlff_workflow.load(args.workflow)
    except ValueError as e:
        fail(str(e))
    wf_base = mlff_workflow.base_dir(args.workflow)

    print_dual(color_text("===== STB-MLFFVALIDATE REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Workflow file     : {args.workflow}" + ("" if state else " (not found)"), f_out)

    # ---- what to score
    committees = []   # (round number or label, [model paths], round entry or None)
    if args.model:
        missing = [m for m in args.model if not os.path.isfile(m)]
        if missing:
            fail(f"Model file(s) not found: {', '.join(missing)}")
        committees.append(("models", args.model, None))
    else:
        if not state or not state["rounds"]:
            fail(f"No training rounds in '{args.workflow}' -- run stb-mlffAnalysis (4.17.2) first, "
                 "or pass --model.")
        numbers = args.rounds or [r["round"] for r in state["rounds"]]
        for k in numbers:
            entry = mlff_workflow.get_round(state, k)
            if entry is None:
                fail(f"Round {k} is not in '{args.workflow}'.")
            models = [m for m in mlff_workflow.committee_models(entry, wf_base) if os.path.isfile(m)]
            if not models:
                print_dual(color_text(f"[WARNING] Round {k}: no model file found -- skipped.", 'yellow'), f_out)
                continue
            committees.append((k, models, entry))
        if not committees:
            fail("None of the requested rounds has a model file.")

    if args.test_path:
        test_dirs = sorted(d for d in glob.glob(args.test_path) if os.path.isdir(d))
        test_source = f"--test-path '{args.test_path}'"
    elif state and mlff_workflow.folders(state, wf_base, "test"):
        test_dirs = [d for d in mlff_workflow.folders(state, wf_base, "test") if os.path.isdir(d)]
        test_source = f"{args.workflow} (test set)"
    else:
        test_dirs = sorted(d for d in glob.glob(os.path.join(wf_base, "mlff_test", "test_*")) if os.path.isdir(d))
        test_source = "'mlff_test/test_*'"
    print_dual(f"Test set from     : {test_source}", f_out)
    print_dual(f"Scoring           : " + ", ".join(
        f"round {k} ({len(m)} model(s))" if k != "models" else f"{len(m)} model(s) from --model"
        for k, m, _ in committees), f_out)
    if args.phonon_ref:
        print_dual(f"Phonon reference  : {args.phonon_ref} (mesh {' '.join(map(str, args.phonon_mesh))})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] TEST SET", f_out)
    if not test_dirs:
        fail("No test configs found -- generate them with stb-mlff --n-test N (and run SIESTA "
             "in them), or pass --test-path.")
    # One folder at a time, so each config keeps its folder (and so its kind).
    test_atoms, groups, used_dirs, n_skipped = [], [], [], 0
    for d in test_dirs:
        found, skipped, _ = collect_training_configs([d], f_out)
        n_skipped += skipped
        for a in found:
            entry = mlff_workflow.config_entry(state, d, wf_base) if state else None
            test_atoms.append(a)
            groups.append(entry.get("kind", "test") if entry else "test")
            used_dirs.append(d)
    if not test_atoms:
        fail("No finished test calculation -- run SIESTA in the test folders first.")
    print_dual(f"[OK] {len(test_atoms)} finished test config(s) ({n_skipped} skipped): " +
               ", ".join(f"{g} {groups.count(g)}" for g in sorted(set(groups))), f_out)
    if len(test_atoms) < 5:
        print_dual(color_text("[WARNING] Fewer than 5 test configs -- differences between rounds "
                              "of a few % are within the noise.", 'yellow'), f_out)
    if state:
        trained = {c for r in state["rounds"] for c in r.get("configs", [])}
        overlap = [d for d in used_dirs if mlff_workflow.rel(d, wf_base) in trained]
        if overlap:
            print_dual(color_text(f"[WARNING] {len(overlap)} test folder(s) were also trained on -- "
                                  "the scores are optimistic.", 'yellow'), f_out)

    try:
        from mace.calculators import MACECalculator
    except ImportError:
        fail("mace-torch not found. This tool needs the optional 'ml' extra: pip install stb_suite[ml]")

    ref_phonon = None
    if args.phonon_ref:
        if not os.path.isfile(os.path.join(args.phonon_ref, "phonopy_disp.yaml")):
            fail(f"'{args.phonon_ref}' has no phonopy_disp.yaml (expected a stb-phononsCreate folder).")
        try:
            ref_phonon, to_ang = load_siesta_phonons(args.phonon_ref, f_out)
        except SystemExit:
            raise
        except Exception as e:
            fail(f"Could not load the SIESTA phonons from '{args.phonon_ref}': {e}")

    print_section("[2] SCORES", f_out)
    scored = []
    for k, models, entry in committees:
        calcs = [MACECalculator(model_paths=[m], device=args.device, default_dtype="float64")
                 for m in models]
        res = evaluate_round(calcs, test_atoms, groups)
        if ref_phonon is not None:
            res["phonons"] = phonon_errors(ref_phonon, committee_phonons(calcs, ref_phonon, to_ang),
                                           args.phonon_mesh)
        scored.append((k, res))
        o = res["overall"]
        title = f"Round {k}" if k != "models" else "--model committee"
        print_dual(color_text(f"{title} ({len(models)} model(s), committee mean):", 'bold'), f_out)
        print_dual(f"  E RMSE            : {o['E_rmse_meV']:8.2f} meV/atom "
                   f"(offset removed {o['E_rmse_rel_meV']:.2f})", f_out)
        print_dual(f"  F RMSE / MAE      : {o['F_rmse_meV']:8.1f} / {o['F_mae_meV']:.1f} meV/Ang", f_out)
        if o["S_rmse_GPa"] is not None:
            print_dual(f"  Stress RMSE       : {o['S_rmse_GPa']:8.3f} GPa", f_out)
        for g, v in res["groups"].items():
            if len(res["groups"]) > 1:
                print_dual(f"    {g:<14}: F {v['F_rmse_meV']:7.1f} meV/Ang, E {v['E_rmse_meV']:6.2f} "
                           f"meV/atom ({v['n']} config(s))", f_out)
        if len(res["members"]) > 1:
            fr = [m["F_rmse_meV"] for m in res["members"]]
            print_dual(f"  Members F RMSE    : {', '.join(f'{x:.1f}' for x in fr)} meV/Ang "
                       f"(mean of members {np.mean(fr):.1f}, committee mean {o['F_rmse_meV']:.1f})", f_out)
        cal = res.get("calibration")
        if cal:
            corr_a = "n/a" if cal["atom_corr"] is None else f"{cal['atom_corr']:.2f}"
            corr_c = "n/a" if cal["config_corr"] is None else f"{cal['config_corr']:.2f}"
            print_dual(f"  Committee spread  : {cal['spread_mean_meV']:.1f} meV/Ang (mean max); true error / "
                       f"spread = {cal['error_over_spread']:.1f}x; corr per config {corr_c}, per atom {corr_a}", f_out)
            if cal["error_over_spread"] > _CALIBRATION_WARN_RATIO:
                print_dual(color_text(
                    f"  [NOTE] The spread underestimates the error ~{cal['error_over_spread']:.0f}x: use it "
                    "to RANK active-learning candidates, not as an error bar.", 'yellow'), f_out)
        ph = res.get("phonons")
        if ph:
            opt = "" if ph["optical_rmse_THz"] is None else f", optical {ph['optical_rmse_THz']:.2f}"
            print_dual(f"  Phonons vs SIESTA : RMSE {ph['rmse_THz']:.2f} THz (acoustic "
                       f"{ph['acoustic_rmse_THz']:.2f}{opt}), max |dev| {ph['max_abs_THz']:.2f}; "
                       f"top {ph['model_max_THz']:.2f} vs {ph['ref_max_THz']:.2f} THz", f_out)

    baseline = None
    if not args.skip_foundation_comparison:
        size = args.foundation_model
        if size is None:
            last = committees[-1][2]
            size = (last or {}).get("settings", {}).get("foundation_model", "small")
        baseline = foundation_baseline(size, args.device, test_atoms)
        s_txt = "" if baseline["S_rmse_GPa"] is None else f", stress {baseline['S_rmse_GPa']:.3f} GPa"
        print_dual(f"Foundation ({size}, no fine-tuning): F {baseline['F_rmse_meV']:.1f} meV/Ang, "
                   f"E (offset removed) {baseline['E_rmse_rel_meV']:.2f} meV/atom{s_txt}", f_out)

    print_section("[3] VERDICT", f_out)
    rounds_only = [(k, v) for k, v in scored if k != "models"]
    code = None
    if rounds_only:
        print_dual(f"{'round':>5}  {'F RMSE':>8}  {'E RMSE':>8}  {'S RMSE':>7}  {'phonons':>8}", f_out)
        for k, v in rounds_only:
            o = v["overall"]
            s_txt = f"{o['S_rmse_GPa']:7.3f}" if o["S_rmse_GPa"] is not None else f"{'-':>7}"
            p_txt = f"{v['phonons']['rmse_THz']:8.2f}" if v.get("phonons") else f"{'-':>8}"
            print_dual(f"{k:>5}  {o['F_rmse_meV']:8.1f}  {o['E_rmse_meV']:8.2f}  {s_txt}  {p_txt}", f_out)
        code, message = verdict([(k, v["overall"]["F_rmse_meV"]) for k, v in rounds_only],
                                args.target_force, args.plateau)
        print_dual(color_text(f"[{code.upper()}] {message}", 'green' if code == "stop" else 'yellow'), f_out)
    else:
        print_dual("(--model committee: no round history to compare.)", f_out)

    png = args.png if os.path.isabs(args.png) else os.path.relpath(os.path.join(wf_base, args.png))
    plot_validation(scored, png)
    print_dual(f"[OK] Saved {png}", f_out)

    if state and rounds_only:
        for k, v in rounds_only:
            entry = mlff_workflow.get_round(state, k)
            entry["validation"] = {
                "time": datetime.now().isoformat(timespec="seconds"),
                "n_test": len(test_atoms),
                **{key: val for key, val in v.items() if not key.startswith("_")},
                "foundation": baseline, "verdict": code if k == rounds_only[-1][0] else None}
        mlff_workflow.record(state, "validate", "stb-mlffValidate", sys.argv[1:])
        mlff_workflow.save(state, args.workflow)
        print_dual(f"[OK] Scores stored in {args.workflow}", f_out)

    print_section("[4] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Figure            : {png}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    if f_out:
        f_out.close()
    print("\n" + "-" * 60)
    print(color_text("Validation done.\n", 'bold'))


if __name__ == "__main__":
    main()
