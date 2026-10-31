#!/usr/bin/env python3
"""The SIESTA check of the 7.6 example (README section 7): the four adatom minima of the main case (T4, B2, H3, T1)
and the clean slab, each relaxed with SIESTA at the 7.1 training level of theory, starting from the model's own
geometry (stb-mladsorb's rank*.fdf), the bottom bilayer fixed as in the ML runs.

From each relaxation it reads the energy and the largest force at every CG step -- step 0 IS the SIESTA single point on
the model's geometry --, the final geometry (where the adatom ended, its height and bonds), and puts next to them:

    E_form(SIESTA)   = E(slab + adatom) - E(slab) - mu_Si,   mu_Si = E/N of the 7.1 bulk relaxation (same settings)
    E_form(SIESTA on the model's geometries) -- the same with step 0 of every run
    E_form of the models (reference_results/<case>/mladsorb_minima.dat)
    the round-3 model evaluated on SIESTA's final geometries (single points; needs mace)

Writes reference_results/siesta/compare.txt and figures/siesta_vs_ml.png, figures/siesta_relaxation.png. Folders still
running are read as far as they got (marked 'running').

    python3 siesta_compare.py [--siesta-dir siesta_runs] [--no-model] [--device cuda]
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference_results")
FIG = os.path.join(HERE, "figures")
SITES = ["T4", "B2", "H3", "T1"]                       # README names
SLOT = {"T1": 1, "T4": 2, "H3": 3, "B2": 5}           # stb-mladsorb site numbers in the main case
MU_SI = -920.141126 / 8                                 # eV/atom: 7.1 bulk relaxation (a = 5.480238 Ang), same settings
MODELS = [("main", "round 3 (fine-tuned)"), ("round2", "round 2 (fine-tuned)"),
          ("foundation_small", "MACE-MP-0 small"), ("foundation_medium", "MACE-MP-0 medium"),
          ("foundation_large", "MACE-MP-0 large")]


def read_run(folder):
    """Energies (E_KS, eV) and largest constrained force (eV/Ang) per CG step, finished or not, final geometry."""
    out = os.path.join(folder, "siesta.out")
    if not os.path.isfile(out):
        return None
    text = open(out, errors="replace").read()
    e = [float(x) for x in re.findall(r"siesta: E_KS\(eV\) =\s+(-?\d+\.\d+)", text)]
    f = [float(x) for x in re.findall(r"^\s+Max\s+(\d+\.\d+)\s+constrained", text, re.M)]
    n = min(len(e), len(f))
    res = {"E": np.array(e[:n]), "F": np.array(f[:n]), "done": "End of run" in text, "atoms": None}
    try:
        from ase.io import read
        res["atoms"] = read(os.path.join(folder, "siesta.xyz"), index=-1) if os.path.isfile(
            os.path.join(folder, "siesta.xyz")) else None
        if os.path.isfile(os.path.join(folder, "siesta.STRUCT_OUT")):
            from ase.io import read as rd
            res["atoms"] = rd(os.path.join(folder, "siesta.STRUCT_OUT"), format="struct_out")
    except Exception:
        pass
    return res


def ml_eform(case):
    """{site name: E_form} of one stb-mladsorb case from its minima table (slot -> README name)."""
    p = os.path.join(REF, case, "mladsorb_minima.dat")
    if not os.path.isfile(p):
        return {}
    name = {v: k for k, v in SLOT.items()}
    out = {}
    for line in open(p):
        if line.startswith("#"):
            continue
        c = line.split()
        slot = int(c[2])
        if slot in name and name[slot] not in out:
            out[name[slot]] = float(c[6])
    return out


def adatom_geometry(at):
    """Height of the adatom (the last atom) above layer 1, its bonds (< 2.9 Ang) and the nearest site of the cell. The
    sites are found on the main case's own relaxed slab (mladsorb_slab.xyz), so their numbers are the report's -- the
    order of pymatgen's sites changes from one slab geometry to another."""
    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor
    from stb import mladsorb as m
    ref = read(os.path.join(REF, "main", "mladsorb_slab.xyz"))
    ref.set_constraint()
    fi = m.SurfaceFrame(ref.cell)
    lay_i = m.layer_index(m.find_layers(fi.height(ref.positions)), len(ref))
    _, _, all_sites, _ = m.find_sites(AseAtomsAdaptor.get_structure(ref), fi, ["ontop", "bridge", "hollow"], 2.0,
                                      0.01, lay_i, ref.get_chemical_symbols())
    fr = m.SurfaceFrame(at.cell)
    layers = m.find_layers(fr.height(at.positions[:-1]))
    cls = m.classify(fr, at.positions[-1], all_sites)
    d = at.get_distances(len(at) - 1, list(range(len(at) - 1)), mic=True)
    return {"height": float(fr.height(at.positions[-1]) - layers[0]["h"]), "bonds": np.sort(d[d < 2.9]),
            "slot": cls["slot"], "dist": cls["dist"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--siesta-dir", default=os.path.join(REF, "siesta"),
                    help="Folder with slab/, T4/, B2/, H3/, T1/ (default: reference_results/siesta)")
    ap.add_argument("--no-model", action="store_true", help="Skip the round-3 model on SIESTA's final geometries")
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    runs = {k: read_run(os.path.join(args.siesta_dir, k)) for k in ["slab"] + SITES}
    if runs["slab"] is None or not len(runs["slab"]["E"]):
        sys.exit(f"no clean-slab run in {args.siesta_dir}/slab")
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    slab = runs["slab"]
    say("SIESTA check of the 7.6 example: relaxations at the 7.1 training level of theory (PBE, DZP, 300 Ry, k 4x4x1),")
    say("each started from the round-3 model's geometry, bottom bilayer fixed, CG to 0.05 eV/Ang.")
    say(f"mu_Si = {MU_SI:.6f} eV/atom (7.1 bulk relaxation, a = 5.480238 Ang)")
    say(f"clean slab: {len(slab['E'])} CG steps ({'done' if slab['done'] else 'RUNNING'}), max|F| {slab['F'][0]:.3f} -> "
        f"{slab['F'][-1]:.3f} eV/Ang, E {slab['E'][0]:.4f} -> {slab['E'][-1]:.4f} eV "
        f"(relaxation from the model's geometry {slab['E'][-1] - slab['E'][0]:+.4f} eV)")
    ideal = read_run(os.path.join(args.siesta_dir, "slab_ideal"))
    if ideal is not None and len(ideal["E"]):
        say(f"ideal slab (as cut, single point): E {ideal['E'][-1]:.4f} eV -> SIESTA's own relaxation of the clean slab "
            f"gains {slab['E'][-1] - ideal['E'][-1]:+.4f} eV; the model's relaxed slab lies "
            f"{slab['E'][0] - ideal['E'][-1]:+.4f} eV from the ideal one in SIESTA")
    # the isolated atom: SIESTA's cohesive energy against the model's E_ads - E_form
    for tag, lab in (("atom", "spin-unpolarized (the training setting)"), ("atom_spin", "spin-polarized")):
        r = read_run(os.path.join(args.siesta_dir, tag))
        if r is not None and len(r["E"]):
            say(f"isolated Si atom, {lab}: E {r['E'][-1]:.4f} eV -> cohesive energy E(atom) - mu = "
                f"{r['E'][-1] - MU_SI:.3f} eV")
    rep = os.path.join(REF, "main", "stb_mladsorb_report.txt")
    if os.path.isfile(rep):
        t = open(rep).read()
        m_iso = re.search(r"^Si\s+\|\s+1\s+\|[^|]+\|[^|]+\|\s+(-?\d+\.\d+)", t, re.M)
        m_mu = re.search(r"si_bulk\.fdf\s*\|\s*\d+\s*\|[^|]+\|\s*(-?\d+\.\d+)", t)
        if m_iso and m_mu:
            e_iso, mu_m = float(m_iso.group(1)), float(m_mu.group(1))
            say(f"round-3 model: isolated atom {e_iso:.4f} eV, mu {mu_m:.4f} eV/atom -> cohesive energy {e_iso - mu_m:.3f} "
                "eV: its isolated atom is the reference it never saw")
    say()
    ml = {c: ml_eform(c) for c, _ in MODELS}
    rows = {}
    say(f"{'site':<5}{'steps':>6}  {'state':<8}{'F0':>7}{'Fend':>7}  {'E_form SIESTA':>14}{'@ML geom':>10}"
        f"{'relax':>8}  " + "".join(f"{lab.split(' (')[0]:>18}" for _, lab in MODELS))
    for s in SITES:
        r = runs[s]
        if r is None or not len(r["E"]):
            say(f"{s:<5}  not run")
            continue
        ef = r["E"][-1] - slab["E"][-1] - MU_SI
        ef0 = r["E"][0] - slab["E"][0] - MU_SI
        rows[s] = {"ef": ef, "ef0": ef0, "run": r}
        say(f"{s:<5}{len(r['E']):>6}  {'done' if r['done'] else 'RUNNING':<8}{r['F'][0]:>7.3f}{r['F'][-1]:>7.3f}  "
            f"{ef:>14.4f}{ef0:>10.4f}{ef - ef0:>8.3f}  "
            + "".join(f"{ml[c].get(s, np.nan):>18.4f}" for c, _ in MODELS))
    say("E_form = E(slab+adatom) - E(slab) - mu_Si (eV); '@ML geom': step 0, SIESTA on the model's own geometry; "
        "'relax': what SIESTA's relaxation gained on top of it; F0/Fend: SIESTA's largest force at the model's geometry "
        "/ at the end (eV/Ang).")
    say()
    if "T4" in rows:
        say("Energy relative to T4 (eV):")
        say(f"{'site':<5}{'SIESTA':>10}{'@ML geom':>10}  " + "".join(f"{lab.split(' (')[0]:>18}" for _, lab in MODELS))
        for s in SITES:
            if s not in rows:
                continue
            say(f"{s:<5}{rows[s]['ef'] - rows['T4']['ef']:>10.4f}{rows[s]['ef0'] - rows['T4']['ef0']:>10.4f}  "
                + "".join(f"{ml[c].get(s, np.nan) - ml[c].get('T4', np.nan):>18.4f}" for c, _ in MODELS))
        say()

    # geometry at the end of each SIESTA relaxation
    if slab["atoms"] is not None:
        say("Where SIESTA's relaxation ended (adatom = the last atom):")
        for s in SITES:
            if s in rows and rows[s]["run"]["atoms"] is not None:
                try:
                    g = adatom_geometry(rows[s]["run"]["atoms"])
                    rows[s]["geom"] = g
                    end = {v: k for k, v in SLOT.items()}.get(g["slot"], f"#{g['slot']}")
                    say(f"  {s}: ended at {end} ({g['dist']:.2f} Ang from it), {g['height']:.3f} Ang above layer 1, "
                        f"bonds {' '.join(f'{b:.3f}' for b in g['bonds'])}")
                except Exception as e:
                    say(f"  {s}: geometry not analysed ({e})")
        say()

    # the round-3 model on SIESTA's final geometries
    if not args.no_model and slab["atoms"] is not None:
        try:
            from stb.core import model_card, mace_relax
            r3 = os.path.join(HERE, "..", "7.5-mlneb", "model_round3")
            mp, cp = os.path.join(r3, "mlff_model_compiled.model"), os.path.join(r3, "mlff_model_card.json")
            model_card.register(mp, cp, echo=False)
            model_card.mark_announced(mp)
            calc = mace_relax.get_calculator(model=mp, device=args.device)
            a = slab["atoms"].copy()
            a.calc = calc
            e_slab = a.get_potential_energy()
            from stb.core import structure_io
            from pymatgen.io.ase import AseAtomsAdaptor
            bulk = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(os.path.join(HERE,
                                                                                                          "si_bulk.fdf"))))
            bulk.calc = calc
            mu = bulk.get_potential_energy() / len(bulk)
            say("Round-3 model on SIESTA's final geometries (single points, its own mu_Si at a = 5.4802 Ang):")
            for s in SITES:
                if s in rows and rows[s]["run"]["atoms"] is not None:
                    b = rows[s]["run"]["atoms"].copy()
                    b.calc = calc
                    ef_m = b.get_potential_energy() - e_slab - mu
                    fm = np.linalg.norm(b.get_forces()[8:], axis=1).max()
                    rows[s]["ef_model_on_siesta"] = ef_m
                    say(f"  {s}: E_form {ef_m:.4f} eV (SIESTA {rows[s]['ef']:.4f}), model's largest force on the free "
                        f"atoms {fm:.3f} eV/Ang")
            say()
        except Exception as e:
            say(f"(round-3 model on SIESTA's geometries skipped: {e})")

    os.makedirs(os.path.join(REF, "siesta"), exist_ok=True)
    with open(os.path.join(REF, "siesta", "compare.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    names = {v: k for k, v in SLOT.items()}
    with open(os.path.join(REF, "siesta", "siesta_eform.dat"), "w") as f:
        f.write("# site E_form_relaxed(eV) E_form_at_model_geometry(eV) ended_at steps finished\n")
        for s_ in SITES:
            if s_ in rows:
                g = rows[s_].get("geom") or {}
                end = names.get(g.get("slot"), "off-site") if g else "?"
                f.write(f"{s_} {rows[s_]['ef']:.6f} {rows[s_]['ef0']:.6f} {end} {len(rows[s_]['run']['E'])} "
                        f"{int(rows[s_]['run']['done'])}\n")

    # figures
    have = [s for s in SITES if s in rows]
    if not have:
        return
    fig, ax = plt.subplots(figsize=(10, 4.6))
    names = {v: k for k, v in SLOT.items()}
    moved = {s_: names.get((rows[s_].get("geom") or {}).get("slot")) for s_ in have}
    moved = {s_: e for s_, e in moved.items() if e and e != s_}
    series = [("SIESTA, relaxed", [np.nan if s in moved else rows[s]["ef"] for s in have], "k"),
              ("SIESTA on the model's geometry", [rows[s]["ef0"] for s in have], "0.6")]
    cols = ["#1f77b4", "#17becf", "#ff7f0e", "#d62728", "#9467bd"]
    for (c, lab), col in zip(MODELS, cols):
        if ml[c]:
            series.append((lab, [ml[c].get(s, np.nan) for s in have], col))
    w = 0.8 / len(series)
    x = np.arange(len(have))
    for k, (lab, vals, col) in enumerate(series):
        ax.bar(x + (k - (len(series) - 1) / 2) * w, vals, w, color=col, label=lab, edgecolor="k", lw=0.4)
    for k, s_ in enumerate(have):
        if s_ in moved:
            ax.annotate(f"SIESTA: -> {moved[s_]}", (k - (len(series) - 1) / 2 * w, 0), xytext=(0, -14),
                        textcoords="offset points", fontsize=8, rotation=90, va="top", ha="center")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(have)
    ax.set_ylabel("E_form = E(slab+Si) - E(slab) - mu_Si (eV)")
    ax.set_title("Si adatom on Si(111) 2x2: formation energy per site, SIESTA vs the models")
    ax.legend(fontsize=7.5, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "siesta_vs_ml.png"), dpi=150)
    plt.close(fig)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for s, col in zip(have, cols):
        r = rows[s]["run"]
        a1.plot(np.arange(len(r["E"])), r["E"] - r["E"][0], "o-", ms=3, color=col, label=s)
        a2.semilogy(np.arange(len(r["F"])), r["F"], "o-", ms=3, color=col, label=s)
    r = slab
    a1.plot(np.arange(len(r["E"])), r["E"] - r["E"][0], "k--", label="clean slab")
    a2.semilogy(np.arange(len(r["F"])), r["F"], "k--", label="clean slab")
    a2.axhline(0.05, color="0.5", ls=":", label="0.05 eV/Ang")
    a1.set_xlabel("SIESTA CG step (0 = the model's geometry)")
    a1.set_ylabel("E - E(step 0) (eV)")
    a1.set_title("(a) energy SIESTA gains from the model's geometry")
    a2.set_xlabel("SIESTA CG step")
    a2.set_ylabel("largest force on the free atoms (eV/Ang)")
    a2.set_title("(b) SIESTA's forces")
    a1.legend(fontsize=8)
    a2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "siesta_relaxation.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
