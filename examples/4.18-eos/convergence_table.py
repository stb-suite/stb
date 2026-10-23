#!/usr/bin/env python3
"""Side-by-side summary of the Si equation-of-state variants of README section 5: for each folder of
reference_results/si/ it reads the 11 finished SIESTA runs, fits a Birch-Murnaghan EOS with the same
code as stb-eosAnalysis, and prints V0, a0, B0, B0', the bulk modulus at the all-electron volume, the
total-energy shift and the Delta-test distance to the all-electron (WIEN2k) reference. With --png it
also draws figures/convergence.png.

    python3 convergence_table.py [--png] [--dir reference_results/si]

The all-electron numbers are the Si row of WIEN2k.txt of the Delta-test repository
(https://github.com/molmod/DeltaCodesDFT): V0 = 20.4530 A^3/atom, B0 = 88.545 GPa, B1 = 4.31.
"""
import argparse
import glob
import os

import numpy as np

from stb.core import siesta_log
from stb.core.eos_fit import fit_eos
from stb.eos_analysis import curve_tools

GPA = 160.21766208
REF = dict(v0=20.4530, b0=88.545, b1=4.31)          # all-electron PBE, WIEN2k 13.1
GROUPS = [
    ("k-points (Mesh 320 Ry, DZP, 0.02 Ry)", ["k04", "k06", "k08", "base", "k14"], ["4³", "6³", "8³", "10³", "14³"]),
    ("Mesh.Cutoff (k 10³, DZP, 0.02 Ry)", ["m150", "m200", "base", "m450", "m600"], ["150", "200", "320", "450", "600"]),
    ("Basis (k 10³, Mesh 320 Ry)", ["base", "dzp_s010", "dzp_s005", "tzp_s020", "tzp_s010", "tzp_s005"],
     ["DZP\n.02", "DZP\n.01", "DZP\n.005", "TZP\n.02", "TZP\n.01", "TZP\n.005"]),
]


def delta_mev(v0a, b0, b1):
    from ase.utils.deltacodesdft import delta
    return 1000 * delta(v0a, b0 / GPA, b1, REF["v0"], REF["b0"] / GPA, REF["b1"])


def analyse(folder):
    V, E, P = [], [], []
    for d in sorted(glob.glob(os.path.join(folder, "vol_*"))):
        out = os.path.join(d, "calc.out")
        e, cell, s = siesta_log.get_free_energy(out), siesta_log.get_outcell(out), siesta_log.get_stress_tensor(out)
        if e is None or cell is None:
            continue
        V.append(abs(np.linalg.det(cell)))
        E.append(e)
        P.append(np.nan if s is None else -np.trace(s) / 3 * GPA)
    V, E, P = np.array(V), np.array(E), np.array(P)
    fit = fit_eos(V, E, "birchmurnaghan")
    e_of_v, p_of_v, b_of_v = curve_tools(fit)
    v0a = fit["v0"] / 2
    return {"n": len(V), "v0": v0a, "a0": (8 * v0a) ** (1 / 3), "b0": fit["b0_gpa"], "b1": fit["bprime"],
            "e0": fit["e0"] / 2, "b_ref": float(b_of_v(2 * REF["v0"])), "delta": delta_mev(v0a, fit["b0_gpa"], fit["bprime"]),
            "dp": float(np.nanmax(np.abs(P - p_of_v(V)))), "r2": fit["r_squared"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="reference_results/si")
    ap.add_argument("--png", action="store_true", help="also write figures/convergence.png")
    args = ap.parse_args()
    results = {}
    for name in sorted({n for _, names, _ in GROUPS for n in names}):
        folder = os.path.join(args.dir, name)
        if os.path.isdir(folder):
            results[name] = analyse(folder)
    ref = results.get("k14") or results.get("base")
    for title, names, _ in GROUPS:
        print(f"\n{title}")
        print(f"{'run':10s} {'V0/atom':>8s} {'a0 (A)':>8s} {'B0':>7s} {'B0p':>6s} {'B(V_AE)':>8s} {'E0 shift':>10s} {'Delta':>8s} {'max dP':>7s}")
        print(f"{'':10s} {'(A^3)':>8s} {'':>8s} {'(GPa)':>7s} {'':>6s} {'(GPa)':>8s} {'(meV/at)':>10s} {'(meV/at)':>8s} {'(GPa)':>7s}")
        for n in names:
            r = results.get(n)
            if not r:
                continue
            print(f"{n:10s} {r['v0']:8.3f} {r['a0']:8.4f} {r['b0']:7.2f} {r['b1']:6.3f} {r['b_ref']:8.2f} "
                  f"{1000 * (r['e0'] - results['base']['e0']):10.2f} {r['delta']:8.2f} {r['dp']:7.3f}")
    print(f"\nAll-electron WIEN2k (Delta test): V0 {REF['v0']:.3f} A^3/atom, a0 {(8 * REF['v0']) ** (1 / 3):.4f} A, "
          f"B0 {REF['b0']:.3f} GPa, B0' {REF['b1']:.2f}")
    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 3, figsize=(15, 8))
        blue, orange = "#0072B2", "#D55E00"
        for col, (title, names, labels) in enumerate(GROUPS):
            have = [(n, l) for n, l in zip(names, labels) if n in results]
            x = np.arange(len(have))
            v0s = [results[n]["v0"] for n, _ in have]
            b0s = [results[n]["b0"] for n, _ in have]
            a = ax[0, col]
            b = a.twinx()
            segments = [(0, len(have))] if col < 2 else [(0, 3), (3, len(have))]     # DZP and TZP are separate series
            for lo, hi in segments:
                a.plot(x[lo:hi], v0s[lo:hi], "o-", color=blue)
                b.plot(x[lo:hi], b0s[lo:hi], "s--", color=orange)
            a.axhline(REF["v0"], color=blue, ls=":", lw=1)
            b.axhline(REF["b0"], color=orange, ls=":", lw=1)
            a.set_ylabel("V0 (Å³/atom)", color=blue)
            b.set_ylabel("B0 (GPa)", color=orange)
            a.set_xticks(x)
            a.set_xticklabels([l for _, l in have], fontsize=9)
            a.set_title(title, fontsize=10)
            a.grid(alpha=0.3)
            a = ax[1, col]
            if col == 0:
                e = [1000 * (results[n]["e0"] - results[have[-1][0]]["e0"]) for n, _ in have]
                a.plot(x, e, "o-", color="#009E73")
                a.set_ylabel("E0 − E0(densest) (meV/atom)")
                a.set_yscale("symlog", linthresh=0.1)
            elif col == 1:
                dp = [results[n]["dp"] for n, _ in have]
                a.plot(x, dp, "o-", color="#009E73")
                a.set_ylabel("max |P_SIESTA − P_fit| (GPa)")
                a.set_ylim(0, max(dp) * 1.2)
            else:
                d = [results[n]["delta"] for n, _ in have]
                a.bar(x, d, color="#AA3377")
                for xi, v in zip(x, d):
                    a.text(xi, v + 0.05, f"{v:.2f}", ha="center", fontsize=8)
                a.set_ylabel("Δ vs all-electron (meV/atom)")
            a.set_xticks(x)
            a.set_xticklabels([l for _, l in have], fontsize=9)
            a.grid(alpha=0.3)
        fig.suptitle("Si equation of state: what converges and what does not (dotted: all-electron WIEN2k reference)", fontsize=13)
        fig.tight_layout()
        os.makedirs("figures", exist_ok=True)
        fig.savefig("figures/convergence.png", dpi=130)
        print("figures/convergence.png written")


if __name__ == "__main__":
    main()
