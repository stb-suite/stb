#!/usr/bin/env python3
"""Silicon's elastic constants from SIESTA at the level of theory the fine-tuned model was trained on, against the model
(reference_results/base/) and experiment (README section 5.3). Needs only numpy and ase; the SIESTA runs themselves are
made by siesta_elastic.sh.

    python3 siesta_compare.py                  compare (reads reference_results/siesta_elastic/siesta_runs.dat)
    python3 siesta_compare.py --collect DIR    extract the finished runs of DIR (siesta_elastic.sh's work folder) into
                                               reference_results/siesta_elastic/siesta_runs.dat, then compare

The SIESTA runs (siesta_elastic.sh, ~6 min on 6 cores): the 8-atom cell si8.fdf -- SIESTA's own minimum at these
settings -- with the training calc.fdf (PBE, DZP, PAO.EnergyShift 0.01 Ry, Mesh.Cutoff 300 Ry, k 6x6x6, Si.psml):
  clamped/  xx and xy strains, -1 .. +1 %, 7 points, single points            -> clamped-ion constants
  relaxed/  the same cells, atoms relaxed by CG to 0.001 eV/Ang (fixed cell) -> relaxed-ion constants (as the model)
  hydro/    isotropic strains, -2 .. +2 % in lattice, 9 points               -> B from E(V) and from the pressure

What it computes, into reference_results/siesta_elastic/siesta_compare.json and figures/siesta_vs_model.png:
  1. C11, C12, C44 (relaxed and clamped) from the stress -- the slopes of sigma_xx, (sigma_yy + sigma_zz)/2 under the
     xx strain and sigma_xy under the xy shear, the same fit as stb-elasticAnalysis;
  2. the same constants from the ENERGY: C11 and C44 from the curvature of E(strain), B from a Birch-Murnaghan E(V),
     C12 = (3B - C11)/2. Two routes through the same runs; their difference is SIESTA's numerical noise (mesh, SCF);
  3. the moduli and derived quantities (Hill B, G, E, nu, Pugh, Cauchy, Zener, sound velocities, Debye temperature)
     with the same functions stb-mlelastic uses, for SIESTA, the model and experiment;
  4. the internal relaxation under shear: how far SIESTA's atoms move per % strain, against the model's.
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "reference_results", "siesta_elastic")
DAT = os.path.join(OUT_DIR, "siesta_runs.dat")
BASE = os.path.join(HERE, "reference_results", "base")
FIG = os.path.join(HERE, "figures", "siesta_vs_model.png")
sys.path.insert(0, os.path.join(HERE, "..", "..", "stb-suite", "src"))

EV_A3_GPA = 160.21766208
BOHR = 0.529177210903
# McSkimin & Andreatch 1964 (ultrasound, 298 K): C11, C12, C44 in GPa; README reference [15]
EXPERIMENT = {"C11": 165.6, "C12": 63.9, "C44": 79.5}
AMU = 1.66053906660e-27
A_SIESTA, A_EXPERIMENT = 5.48024, 5.4310   # Ang: SIESTA's (and the model's, -0.01 %) cell; experiment at 300 K


# ---------------------------------------------------------------------------------------------- collect
def read_lattice_fdf(path):
    txt = open(path).read()
    scale = float(re.search(r"LatticeConstant\s+([-\d.eE+]+)", txt).group(1))
    block = re.search(r"%block LatticeVectors(.*?)%endblock", txt, re.S | re.I).group(1)
    cell = np.array([[float(x) for x in l.split()[:3]] for l in block.strip().splitlines()]) * scale
    coords = re.search(r"%block AtomicCoordinatesAndAtomicSpecies(.*?)%endblock", txt, re.S | re.I).group(1)
    frac = np.array([[float(x) for x in l.split()[:3]] for l in coords.strip().splitlines()])
    return cell, frac


def read_xv(path):
    lines = open(path).read().split("\n")
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(3)]) * BOHR
    n = int(lines[3].split()[0])
    pos = np.array([[float(x) for x in lines[4 + i].split()[2:5]] for i in range(n)]) * BOHR
    return cell, pos


def read_out(path):
    txt = open(path, errors="ignore").read()
    if "End of run" not in txt:
        return None
    energy = float(re.findall(r"siesta: FreeEng\s*=\s*([-\d.]+)", txt)[-1])
    block = re.findall(r"siesta: Stress tensor \(static\) \(eV/Ang\*\*3\):\s*\n"
                       r"siesta:(.*)\nsiesta:(.*)\nsiesta:(.*)\n", txt)[-1]
    s = np.array([[float(x) for x in row.split()] for row in block]) * EV_A3_GPA
    fmax = [float(x) for x in re.findall(r"^\s+Max\s+([\d.]+)\s*$", txt, re.M)]
    steps = len(re.findall(r"Begin CG opt. move", txt))
    return energy, s, (fmax[-1] if fmax else float("nan")), max(steps - 1, 0)


def collect(work):
    ref_cell, ref_frac = read_lattice_fdf(os.path.join(work, "structure.fdf"))
    rows = []
    for folder in sorted(glob.glob(os.path.join(work, "*", "*", "strain_*"))):
        res = read_out(os.path.join(folder, "siesta.out"))
        if res is None:
            print(f"[skip] {folder}: unfinished")
            continue
        energy, s, fmax, steps = res
        cell, frac = read_lattice_fdf(os.path.join(folder, "structure.fdf"))
        F = np.linalg.solve(ref_cell, cell).T                 # cell rows: cell = ref_cell @ F^T
        eps = 0.5 * (F + F.T) - np.eye(3)
        voigt_e = [eps[0, 0], eps[1, 1], eps[2, 2], 2 * eps[1, 2], 2 * eps[0, 2], 2 * eps[0, 1]]
        shift = 0.0
        xv = os.path.join(folder, "siesta.XV")
        if os.path.exists(xv):
            _, pos = read_xv(xv)
            d = pos - frac @ cell                             # vs the homogeneously strained positions
            d -= np.round(d @ np.linalg.inv(cell)) @ cell
            d -= d.mean(axis=0)                               # a rigid translation is not a relaxation
            shift = float(np.linalg.norm(d, axis=1).max())
        rel = os.path.relpath(folder, work).split(os.sep)
        rows.append([rel[0], rel[1], *voigt_e, abs(np.linalg.det(cell)), energy,
                     s[0, 0], s[1, 1], s[2, 2], s[1, 2], s[0, 2], s[0, 1], fmax, steps, shift])
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(DAT, "w") as f:
        f.write("# SIESTA runs of siesta_elastic.sh (8-atom Si, training level of theory). Strain from the cells "
                "(Voigt, engineering shear);\n# E = FreeEng (eV); stress = 'Stress tensor (static)' in GPa, positive "
                "= tension; fmax = last max |F| (eV/Ang); shift = largest internal displacement (Ang)\n")
        f.write("# set dir e_xx e_yy e_zz g_yz g_xz g_xy V(Ang^3) E(eV) s_xx s_yy s_zz s_yz s_xz s_xy fmax "
                "cg_steps shift\n")
        for r in rows:
            f.write(f"{r[0]:8s} {r[1]:6s}" + "".join(f" {x:13.7f}" for x in r[2:17])
                    + f" {r[17]:4d} {r[18]:9.5f}\n")
    print(f"{len(rows)} runs -> {os.path.relpath(DAT, HERE)}")


# ---------------------------------------------------------------------------------------------- analyse
def load():
    data = []
    for line in open(DAT):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        data.append({"set": p[0], "dir": p[1], "e": np.array(p[2:8], float), "V": float(p[8]), "E": float(p[9]),
                     "s": np.array(p[10:16], float), "fmax": float(p[16]), "steps": int(p[17]),
                     "shift": float(p[18])})
    return data


def slope(x, y):
    a, b = np.polyfit(x, y, 1)
    r2 = 1 - np.sum((y - (a * x + b)) ** 2) / np.sum((y - y.mean()) ** 2)
    return a, b, r2


def stress_constants(data, which):
    xx = [d for d in data if d["set"] == which and d["dir"] == "xx"]
    xy = [d for d in data if d["set"] == which and d["dir"] == "xy"]
    e = np.array([d["e"][0] for d in xx])
    c11 = slope(e, np.array([d["s"][0] for d in xx]))
    c12 = slope(e, np.array([(d["s"][1] + d["s"][2]) / 2 for d in xx]))
    g = np.array([d["e"][5] for d in xy])
    c44 = slope(g, np.array([d["s"][5] for d in xy]))
    return {"C11": c11[0], "C12": c12[0], "C44": c44[0], "R2": {"C11": c11[2], "C12": c12[2], "C44": c44[2]},
            "intercept_GPa": max(abs(c11[1]), abs(c44[1]))}


def energy_curvature(points, V0):
    """E = E0 + V0/2 C x^2 (+ cubic): returns C in GPa."""
    x = np.array([p[0] for p in points])
    E = np.array([p[1] for p in points])
    c = np.polyfit(x, E, 3)
    return 2 * c[1] / V0 * EV_A3_GPA


def energy_constants(data, which):
    xx = [d for d in data if d["set"] == which and d["dir"] == "xx"]
    xy = [d for d in data if d["set"] == which and d["dir"] == "xy"]
    V0 = [d["V"] for d in xx if abs(d["e"][0]) < 1e-9][0]
    c11 = energy_curvature([(d["e"][0], d["E"]) for d in xx], V0)
    c44 = energy_curvature([(d["e"][5], d["E"]) for d in xy], V0)
    return {"C11": c11, "C44": c44, "V0": V0}


def bulk_from_hydro(data):
    from ase.eos import EquationOfState
    h = sorted([d for d in data if d["set"] == "hydro"], key=lambda d: d["V"])
    V = np.array([d["V"] for d in h])
    E = np.array([d["E"] for d in h])
    P = np.array([-(d["s"][0] + d["s"][1] + d["s"][2]) / 3 for d in h])
    v0, _, b = EquationOfState(V, E, eos="birchmurnaghan").fit()
    Vr = [d["V"] for d in h if abs(d["e"][0]) < 1e-9][0]
    c = np.polyfit(V, P, 3)
    b_stress = -Vr * np.polyval(np.polyder(c), Vr)
    # the local slope of P(V) between symmetric pairs: how much the stress alone scatters
    n = len(V) // 2
    local = {f"+-{abs(h[n + k]['e'][0]) * 100:.1f} %": -Vr * (P[n + k] - P[n - k]) / (V[n + k] - V[n - k])
             for k in range(1, n + 1)}
    return {"B_EV": b * EV_A3_GPA, "V0_EV": v0, "V_ref": Vr, "B_stress_cubic": b_stress,
            "B_stress_local": local, "P_ref": float(P[n])}


def model_constants():
    rows = [l for l in open(os.path.join(BASE, "elastic_tensor.dat")).read().split("\n")]
    i = rows.index("# projected on the point group")
    C = np.array([[float(x) for x in rows[i + 1 + k].split()] for k in range(6)])
    txt = open(os.path.join(BASE, "stb_mlelastic_report.txt")).read()
    clamped = float(re.search(r"C44 = C55 = C66\s*\|\s*([\d.]+)\s*\|\s*[\d.]+\s*\|", txt).group(1))
    shift = float(re.search(r"\nyz\s*\|.*?\|.*?\|.*?\|.*?\|.*?\|\s*[\d.]+\s*\|\s*([\d.]+)", txt).group(1))
    return {"C11": C[0, 0], "C12": C[0, 1], "C44": C[3, 3], "C44_clamped": clamped, "shift_per_pct": shift}


def cubic(c11, c12, c44):
    C = np.zeros((6, 6))
    C[:3, :3] = c12
    np.fill_diagonal(C[:3, :3], c11)
    C[3, 3] = C[4, 4] = C[5, 5] = c44
    return C


def derived(c, a):
    """Moduli and derived quantities of a cubic tensor, at the density of the conventional cell of side a (Ang)."""
    from stb.mlelastic import vrh, christoffel, debye_from_vrh
    C = cubic(c["C11"], c["C12"], c["C44"])
    h = vrh(C)["H"]
    B, G, E, nu = h["B"], h["G"], h["E"], h["nu"]
    n_density = 8 / (a * 1e-10) ** 3
    rho = n_density * 28.0855 * AMU
    theta = debye_from_vrh(B, G, rho, n_density)[3]
    v100 = christoffel(C, rho, [1, 0, 0])[0]
    return {"B": B, "G": G, "E": E, "nu": nu, "Pugh B/G": B / G, "Cauchy C12-C44": c["C12"] - c["C44"],
            "Zener 2C44/(C11-C12)": 2 * c["C44"] / (c["C11"] - c["C12"]), "Debye K": theta,
            "v[100] L/T (m/s)": [float(v100[2]), float(v100[0])]}


def figure(data, sets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    # (a) stress-strain: SIESTA points, model lines
    model = {}
    for line in open(os.path.join(BASE, "stress_strain.dat")):
        if line.startswith("#"):
            if "foundation" in line or "MACE-MP" in line:
                break
            continue
        p = line.split()
        if p:
            model.setdefault(p[0], []).append([float(x) for x in p[1:]])
    rx = [d for d in data if d["set"] == "relaxed" and d["dir"] == "xx"]
    ry = [d for d in data if d["set"] == "relaxed" and d["dir"] == "xy"]
    cy = [d for d in data if d["set"] == "clamped" and d["dir"] == "xy"]
    mx, my = np.array(model["xx"]), np.array(model["xy"])
    curves = [("σxx, strain xx (C11)", "C0", rx, lambda d: d["s"][0], mx[:, 1], mx[:, 2], lambda d: d["e"][0]),
              ("σyy, strain xx (C12)", "C1", rx, lambda d: d["s"][1], mx[:, 1], mx[:, 3], lambda d: d["e"][0]),
              ("σxy, shear xy (C44)", "C2", ry, lambda d: d["s"][5], my[:, 1], my[:, 7], lambda d: d["e"][5])]
    for label, col, pts, sf, xm, ym, ef in curves:
        ax1.plot(np.array([ef(d) for d in pts]) * 100, [sf(d) for d in pts], "o", color=col, label=f"SIESTA {label}")
        ax1.plot(xm * 100, ym, "-", color=col, lw=1.2, alpha=0.8)
    ax1.plot(np.array([d["e"][5] for d in cy]) * 100, [d["s"][5] for d in cy], "^", color="C2", mfc="none",
             label="SIESTA σxy, clamped ions")
    ax1.plot([], [], "-", color="0.4", label="fine-tuned model (lines)")
    ax1.set_xlabel("strain (%) — engineering for the shear")
    ax1.set_ylabel("stress (GPa)")
    ax1.set_title("(a) stress against strain, relaxed ions")
    ax1.axhline(0, color="0.8", lw=0.6)
    ax1.legend(fontsize=7.5, loc="upper left")
    # (b) constants
    keys = ["C11", "C12", "C44", "C44⁰", "B", "G"]
    groups = [("SIESTA (stress)", "C0"), ("SIESTA (energy)", "C9"), ("fine-tuned model", "C3"), ("experiment", "0.3")]
    w = 0.2
    for k, (name, col) in enumerate(groups):
        vals = [sets[name].get(key, np.nan) for key in keys]
        ax2.bar(np.arange(len(keys)) + (k - 1.5) * w, vals, w, color=col, label=name)
    ax2.set_xticks(np.arange(len(keys)))
    ax2.set_xticklabels(keys)
    ax2.set_ylabel("GPa")
    ax2.set_title("(b) constants and Hill moduli (C44⁰: clamped ions)")
    ax2.legend(fontsize=8, loc="upper right")
    ax2.set_ylim(0, 190)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130)
    print(f"figure -> {os.path.relpath(FIG, HERE)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--collect", metavar="DIR", help="siesta_elastic.sh's work folder to extract first")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    if args.collect:
        collect(args.collect)
    data = load()

    rel, cla = stress_constants(data, "relaxed"), stress_constants(data, "clamped")
    en_rel, en_cla = energy_constants(data, "relaxed"), energy_constants(data, "clamped")
    hyd = bulk_from_hydro(data)
    mod = model_constants()

    siesta_s = {"C11": rel["C11"], "C12": rel["C12"], "C44": rel["C44"], "C44⁰": cla["C44"]}
    siesta_e = {"C11": en_rel["C11"], "C12": (3 * hyd["B_EV"] - en_rel["C11"]) / 2, "C44": en_rel["C44"],
                "C44⁰": en_cla["C44"]}
    model = {"C11": mod["C11"], "C12": mod["C12"], "C44": mod["C44"], "C44⁰": mod["C44_clamped"]}
    sets = {"SIESTA (stress)": siesta_s, "SIESTA (energy)": siesta_e, "fine-tuned model": model,
            "experiment": dict(EXPERIMENT)}
    for name, v in sets.items():
        v.update(derived(v, A_EXPERIMENT if name == "experiment" else A_SIESTA))

    shear = sorted([d for d in data if d["set"] == "relaxed" and d["dir"] == "xy"], key=lambda d: d["e"][5])
    shift_per_pct = np.polyfit([abs(d["e"][5]) * 50 for d in shear], [d["shift"] for d in shear], 1)[0]

    # ------------------------------------------------------------------ print
    print("\nSIESTA RUNS")
    print(f"  {len(data)} runs; relaxed: max |F| after CG {max(d['fmax'] for d in data if d['set'] == 'relaxed'):.4f}"
          f" eV/Ang, {sum(d['steps'] for d in data if d['set'] == 'relaxed')} CG steps in all")
    print(f"  fit quality (relaxed): R^2 C11 {rel['R2']['C11']:.6f}  C12 {rel['R2']['C12']:.6f}  "
          f"C44 {rel['R2']['C44']:.6f}; largest intercept {rel['intercept_GPa']:.3f} GPa")
    print(f"  pressure of the reference cell: {hyd['P_ref']:+.3f} GPa")
    print("\nBULK MODULUS, three routes through SIESTA (GPa)")
    print(f"  (C11 + 2 C12)/3, stress fits      {(rel['C11'] + 2 * rel['C12']) / 3:7.2f}")
    print(f"  Birch-Murnaghan E(V), +-2 %       {hyd['B_EV']:7.2f}   (V0 {hyd['V0_EV']:.3f} vs reference "
          f"{hyd['V_ref']:.3f} Ang^3)")
    print(f"  -V dP/dV, cubic fit of P(V)       {hyd['B_stress_cubic']:7.2f}")
    print("  -V dP/dV between symmetric pairs  " + "  ".join(f"{k}: {v:.1f}" for k, v in hyd["B_stress_local"].items()))
    print("\nCONSTANTS (GPa)                 C11      C12      C44     C44⁰       B        G      E      nu")
    for name, v in sets.items():
        print(f"  {name:22s} " + " ".join(f"{v.get(k, float('nan')):8.2f}" for k in ("C11", "C12", "C44", "C44⁰", "B", "G", "E"))
              + f"  {v['nu']:.3f}")
    ref = sets["SIESTA (stress)"]
    print("\nMODEL vs SIESTA (stress route) and SIESTA vs EXPERIMENT, %")
    for k in ("C11", "C12", "C44", "C44⁰", "B", "G"):
        m = 100 * (model[k] - ref[k]) / ref[k]
        exp = sets["experiment"]
        x = 100 * (ref[k] - exp[k]) / exp[k] if k in exp else float("nan")
        print(f"  {k:5s} model {m:+6.1f} %   SIESTA vs experiment " + ("    —" if np.isnan(x) else f"{x:+6.1f} %"))
    print("\nINTERNAL RELAXATION UNDER SHEAR")
    print(f"  C44 clamped -> relaxed: SIESTA {siesta_s['C44⁰']:.1f} -> {siesta_s['C44']:.1f} "
          f"({100 * (siesta_s['C44'] / siesta_s['C44⁰'] - 1):+.1f} %); model {model['C44⁰']:.1f} -> {model['C44']:.1f} "
          f"({100 * (model['C44'] / model['C44⁰'] - 1):+.1f} %)")
    print(f"  largest atomic shift per % strain: SIESTA {shift_per_pct:.4f} Ang, model {mod['shift_per_pct']:.4f} Ang")
    print("\nDERIVED (SIESTA by stress and by energy, model, experiment)")
    for k in ("Pugh B/G", "Cauchy C12-C44", "Zener 2C44/(C11-C12)", "Debye K"):
        print(f"  {k:22s} " + "  ".join(f"{short} {sets[n][k]:7.2f}" for n, short in
                                        zip(sets, ("stress", "energy", "model", "experiment"))))

    out = {"siesta_stress": siesta_s, "siesta_energy": siesta_e, "model": model, "experiment": sets["experiment"],
           "fit": rel["R2"], "bulk_routes": hyd, "shift_per_pct": {"siesta": shift_per_pct,
                                                                   "model": mod["shift_per_pct"]}}
    with open(os.path.join(OUT_DIR, "siesta_compare.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"\n-> {os.path.relpath(os.path.join(OUT_DIR, 'siesta_compare.json'), HERE)}")
    if not args.no_figure:
        figure(data, sets)


if __name__ == "__main__":
    main()
