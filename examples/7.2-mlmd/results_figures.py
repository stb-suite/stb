#!/usr/bin/env python3
"""Reads the GPU runs in reference_results/ (written by make_reference.sh) and prints the tables of README section 5,
redrawing its figures in figures/. Needs only numpy and matplotlib (no MACE, no GPU).

    python3 results_figures.py            tables + figures
    python3 results_figures.py --no-plots tables only

Cases read (each a folder of reference_results/ holding stb_mlmd_report.txt and the --save-gnuplot data files):
nvt300, nvt300_foundation, nve_dt{0.5,1,2,3,4,5}, nve_equipartition, size_N{8,64,512}, npt_T{300,...,1500},
npt_T300_{water,si}beta, melt_ft, melt_foundation; and siesta_check/*/forces_compare.json (siesta_check.py).
A missing case is skipped with a note.
"""
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference_results")
FIG = os.path.join(HERE, "figures")
KB = 8.617333262e-5
A_EXP_273 = 5.430741          # Ang, Okada & Tokumaru (1984) at 273.2 K
PLOTS = "--no-plots" not in sys.argv


# ---- reading ---------------------------------------------------------------------------------------

def case_dir(name):
    d = os.path.join(REF, name)
    return d if os.path.isfile(os.path.join(d, "stb_mlmd_report.txt")) else None


def report(name):
    d = case_dir(name)
    if d is None:
        return None
    return re.sub(r"\x1b\[[0-9;]*m", "", open(os.path.join(d, "stb_mlmd_report.txt")).read())


def data(name, kind):
    """kind = thermo | rdf | msd -> array, or None."""
    d = case_dir(name)
    files = glob.glob(os.path.join(d, f"*_md_{kind}.dat")) if d else []
    return np.loadtxt(files[0]) if files else None


def grab(text, pattern, group=1, cast=float):
    m = re.search(pattern, text or "", re.M)
    if not m:
        return None
    try:
        return cast(m.group(group))
    except (TypeError, ValueError):
        return None


def mean_pm(text, label):
    """'Label | mean +/- se | std | range' row of the [6] table -> (mean, se, std)."""
    m = re.search(rf"^{re.escape(label)}\s*\|\s*([-+0-9.]+)(?:\s*\+/-\s*([0-9.]+))?\s*\|\s*([0-9.]+)", text or "", re.M)
    if not m:
        return None, None, None
    return float(m.group(1)), (float(m.group(2)) if m.group(2) else None), float(m.group(3))


def fmt(x, f="{:.3f}", none="-"):
    return none if x is None else f.format(x)


def missing(names):
    gone = [n for n in names if case_dir(n) is None]
    if gone:
        print(f"  (not in reference_results/: {', '.join(gone)})")
    return gone


def table(rows, headers):
    widths = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    line = " | ".join(str(h).ljust(w) for h, w in zip(headers, widths))
    print("  " + line)
    print("  " + "-" * len(line))
    for r in rows:
        print("  " + " | ".join(str(c).ljust(w) for c, w in zip(r, widths)))


def okada_alpha(t):
    return (3.725 * (1 - np.exp(-5.88e-3 * (t - 124))) + 5.548e-4 * t) * 1e-6


def okada_strain(t1, t2):
    """ln(a(t2)/a(t1)) from the Okada-Tokumaru alpha(T)."""
    ts = np.linspace(t1, t2, 2001)
    return float(np.trapezoid(okada_alpha(ts), ts) if hasattr(np, "trapezoid") else np.trapz(okada_alpha(ts), ts))


def block_se(x, n):
    means = [b.mean() for b in np.array_split(np.asarray(x), n)]
    return float(np.std(means, ddof=1) / np.sqrt(n))


# ---- sections --------------------------------------------------------------------------------------

def s_reference():
    print("\n[5.2] Reference run: 216 Si atoms, NVT 300 K, 10 ps, fine-tuned model -- vs the raw foundation model")
    if missing(["nvt300", "nvt300_foundation"]):
        return
    rows = []
    for name, label in [("nvt300", "fine-tuned"), ("nvt300_foundation", "MACE-MP medium")]:
        r = report(name)
        t, t_se, t_sd = mean_pm(r, "Temperature (K)")
        p, p_se, _ = mean_pm(r, "Pressure (GPa)")
        rows.append([label, f"{t:.1f} +/- {t_se:.1f}", f"{t_sd:.1f}", fmt(grab(r, r"canonical T sqrt\(2/N_dof\) = ([0-9.]+)"), "{:.1f}"),
                     fmt(grab(r, r"^Heat capacity\s*:\s*([0-9.]+)"), "{:.2f}"), f"{p:+.3f} +/- {p_se:.3f}",
                     fmt(grab(r, r"first peak ([0-9.]+) Ang"), "{:.3f}"), fmt(grab(r, r"coordination ([0-9.]+)"), "{:.2f}"),
                     fmt(grab(r, r"^Lindemann ratio\s*:\s*([0-9.]+)")), fmt(grab(r, r"^Speed\s*:\s*([0-9.]+) ms/step"), "{:.1f}")])
    table(rows, ["model", "<T> (K)", "sigma_T", "canonical", "C_V (kB)", "<P> (GPa)", "g(r) peak", "coord.", "Lindemann", "ms/step"])

    th = data("nvt300", "thermo")
    if th is not None:
        e = th[:, 2] * 1000.0
        naive = e.std(ddof=1) / np.sqrt(len(e))
        print("\n  Block averaging of E_pot (meV/atom), nvt300, 1000 samples every 10 fs:")
        rows = [[n, len(e) // n * 10, f"{block_se(e, n):.4f}"] for n in (100, 50, 20, 10, 5, 2) if len(e) // n >= 1]
        table(rows, ["blocks", "block length (fs)", "SE of the mean"])
        print(f"  naive sigma/sqrt(M) = {naive:.4f} meV/atom -- too small: consecutive samples are correlated")
        if PLOTS:
            fig_block(e)


def s_siesta():
    files = sorted(glob.glob(os.path.join(REF, "siesta_check", "*", "forces_compare.json")))
    print("\n[5.3] Forces against SIESTA on MD frames (siesta_check.py)")
    if not files:
        print("  (no siesta_check/*/forces_compare.json)")
        return
    rows = []
    for f in files:
        d = json.load(open(f))
        for m in d["models"]:
            rows.append([d["frame"], d["n_atoms"], f"{d['siesta_rms_meV']:.0f}", m["name"], f"{m['F_rmse_meV']:.1f}",
                         f"{m['F_mae_meV']:.1f}", f"{m['F_max_meV']:.0f}", fmt(m.get("dE_meV_atom"), "{:+.1f}")])
    table(rows, ["frame", "N", "SIESTA RMS F", "model", "F RMSE", "F MAE", "max |dF|", "dE (meV/atom)"])
    print("  forces in meV/Ang; dE = model - SIESTA energy per atom after removing the offset measured on the 300 K frame")


def s_timestep():
    names = [f"nve_dt{dt}" for dt in ("0.5", "1", "2", "3", "4", "5")]
    print("\n[5.4] Energy conservation vs timestep: NVE, 216 atoms, 5 ps each, from the equilibrated 300 K structure")
    missing(names)
    rows, curves = [], []
    for n in names:
        r = report(n)
        if r is None:
            continue
        dt = float(n[6:])
        drift = grab(r, r"^Energy drift\s*:\s*([-+0-9.]+)")
        fl = grab(r, r"^E_tot fluctuation\s*:\s*([0-9.]+)")
        ratio = grab(r, r"meV/atom \(std\), ([0-9.]+) x the E_pot")
        t = mean_pm(r, "Temperature (K)")[0]
        rows.append([f"{dt:g}", fmt(drift, "{:+.4f}"), fmt(fl, "{:.4f}"), fmt(ratio, "{:.4f}"), fmt(t, "{:.1f}"),
                     fmt(grab(r, r"^Speed\s*:\s*([0-9.]+) ms/step"), "{:.1f}")])
        th = data(n, "thermo")
        if th is not None:
            curves.append((dt, th[:, 0], (th[:, 4] - th[0, 4]) * 1000.0, fl))
    if rows:
        table(rows, ["dt (fs)", "drift (meV/atom/ps)", "sigma(E_tot) (meV/atom)", "/ sigma(E_pot)", "<T> (K)", "ms/step"])
    if PLOTS and curves:
        fig_timestep(curves)


def s_equipartition():
    print("\n[5.5] Equipartition: perfect 216-atom crystal, velocities at 600 K, NVE")
    th = data("nve_equipartition", "thermo")
    if th is None:
        missing(["nve_equipartition"])
        return
    t, temp = th[:, 0], th[:, 1]
    late = t > t[-1] / 2
    rows = [[f"{tt:.0f}", f"{TT:.1f}"] for tt, TT in zip(t, temp) if tt in (2, 10, 20, 30, 40, 60, 80, 100, 200, 500, 1000)]
    table(rows, ["t (fs)", "T (K)"])
    print(f"  mean T over the second half: {temp[late].mean():.1f} K (600/2 = 300 for a harmonic solid); "
          f"E_tot constant to {np.ptp(th[:, 4]) * 1000:.3f} meV/atom")
    if PLOTS:
        fig_equipartition(th)


def s_size():
    names = ["size_N8", "size_N64", "nvt300", "size_N512"]
    print("\n[5.6] Finite size: NVT 300 K with 8, 64, 216 and 512 atoms")
    missing(names)
    rows, pts = [], []
    for n in names:
        r = report(n)
        if r is None:
            continue
        nat = grab(r, r"Read (\d+) atom", cast=int)
        t, t_se, t_sd = mean_pm(r, "Temperature (K)")
        theory = 300 * np.sqrt(2.0 / (3 * nat - 3))
        p, p_se, p_sd = mean_pm(r, "Pressure (GPa)")
        rows.append([nat, f"{t:.1f} +/- {t_se:.1f}", f"{t_sd:.1f}", f"{theory:.1f}", f"{t_sd / theory:.2f}",
                     fmt(grab(r, r"^Heat capacity\s*:\s*([0-9.]+)"), "{:.2f}"), f"{p:+.3f} +/- {p_se:.3f}", f"{p_sd:.3f}",
                     fmt(grab(r, r"^Lindemann ratio\s*:\s*([0-9.]+)")), fmt(grab(r, r"^Speed\s*:\s*([0-9.]+) ms/step"), "{:.1f}"),
                     fmt(grab(r, r"= ([0-9.]+) us/atom/step"), "{:.0f}")])
        pts.append((nat, t_sd / 300.0, data(n, "thermo")))
    if rows:
        table(rows, ["N", "<T> (K)", "sigma_T", "sqrt(2/N_dof) T", "ratio", "C_V (kB)", "<P> (GPa)", "sigma_P", "Lindemann",
                     "ms/step", "us/atom/step"])
    if PLOTS and pts:
        fig_size(pts)


def s_expansion():
    temps = [300, 600, 900, 1200, 1500]
    names = [f"npt_T{t}" for t in temps]
    print("\n[5.7] Thermal expansion: NPT at 0 GPa, 216 atoms, compressibility 1.0e-6 1/bar, 2 + 6 ps per temperature")
    missing(names)
    pts = []
    for t, n in zip(temps, names):
        th = data(n, "thermo")
        r = report(n)
        if th is None or r is None:
            continue
        a = th[:, 7] / 3.0                       # 3x3x3 conventional cells
        pts.append((t, a.mean(), block_se(a, 5), mean_pm(r, "Pressure (GPa)")[0], mean_pm(r, "Density (g/cm^3)")[0],
                    grab(r, r"^Lindemann ratio\s*:\s*([0-9.]+)"), "DRIFTING" in r))
    if not pts:
        return
    a300 = next((p[1] for p in pts if p[0] == 300), pts[0][1])
    rows = []
    for t, a, se, p, rho, lind, drift in pts:
        md = (a / a300 - 1) * 100
        exp = (np.exp(okada_strain(300, t)) - 1) * 100
        a_exp = A_EXP_273 * np.exp(okada_strain(273.2, t))
        rows.append([t, f"{a:.5f} +/- {se:.5f}", f"{a_exp:.5f}", f"{md:+.3f}", f"{exp:+.3f}", fmt(p, "{:+.3f}"),
                     fmt(rho, "{:.4f}"), fmt(lind), "yes" if drift else ""])
    table(rows, ["T (K)", "a MD (Ang)", "a exp (Ang)", "da/a MD (%)", "da/a exp (%)", "<P> (GPa)", "density", "Lindemann",
                 "V drifting"])
    ts = np.array([p[0] for p in pts], float)
    aa = np.array([p[1] for p in pts])
    if len(ts) >= 3:
        print("  alpha from neighbouring temperatures (1e-6 1/K), MD vs Okada-Tokumaru at the midpoint:")
        for i in range(len(ts) - 1):
            am = np.log(aa[i + 1] / aa[i]) / (ts[i + 1] - ts[i]) * 1e6
            tm = 0.5 * (ts[i] + ts[i + 1])
            print(f"    {ts[i]:.0f}-{ts[i + 1]:.0f} K: MD {am:.2f}, experiment {okada_alpha(tm) * 1e6:.2f}")
    if PLOTS:
        fig_expansion(pts, a300)

    print("\n  The barostat's compressibility: NPT 300 K from the static cell, 4 ps, sampled every 5 fs")
    rows, curves = [], []
    for n, label in [("npt_T300_waterbeta", "4.57e-5 (water, default)"), ("npt_T300_sibeta", "1.0e-6 (silicon)")]:
        th = data(n, "thermo")
        r = report(n)
        if th is None or r is None:
            missing([n])
            continue
        v = th[:, 6]
        late = th[:, 0] > 1000
        rows.append([label, f"{v[late].mean():.4f}", f"{v[late].std():.4f}", f"{mean_pm(r, 'Pressure (GPa)')[2]:.3f}",
                     fmt(grab(r, r"barostat response of B ~ ([0-9.]+) GPa"), "{:.1f}")])
        curves.append((label, th[:, 0], v))
    if rows:
        table(rows, ["compressibility (1/bar)", "<V> after 1 ps (Ang^3/atom)", "sigma_V", "sigma_P (GPa)",
                     "B the barostat assumes (GPa)"])
    if PLOTS and curves:
        fig_barostat(curves)


def s_melt():
    print("\n[5.8] 2500 K at 0 GPa: the melt, fine-tuned model vs the raw foundation model (and the 300/1500 K crystal)")
    names = ["npt_T300", "npt_T1500", "melt_ft", "melt_foundation"]
    missing(names)
    rows, curves = [], []
    for n in names:
        r = report(n)
        if r is None:
            continue
        outside = re.findall(r"^(.+?)\s*\|[^|]*\|[^|]*\|\s*OUTSIDE", r, re.M)
        rows.append([n, fmt(mean_pm(r, "Density (g/cm^3)")[0], "{:.3f}"), fmt(grab(r, r"first peak ([0-9.]+) Ang"), "{:.3f}"),
                     fmt(grab(r, r"coordination ([0-9.]+)"), "{:.2f}"), fmt(grab(r, r"D = ([0-9.eE+-]+) cm\^2/s"), "{:.2e}", "-"),
                     fmt(grab(r, r"^Lindemann ratio\s*:\s*([0-9.]+)")), fmt(grab(r, r"^Shortest distance : min ([0-9.]+)"), "{:.3f}"),
                     ", ".join(o.strip() for o in outside) or "none"])
        curves.append((n, data(n, "rdf"), data(n, "msd")))
    if rows:
        table(rows, ["run", "density", "g(r) peak", "coord.", "D (cm^2/s)", "Lindemann", "d_min", "OUTSIDE the trained range"])
    print("  experiment: crystal 2.33 g/cm^3 at 300 K; liquid ~2.57 g/cm^3 at 1687 K, coordination ~5-6.4, D ~ 2e-4 cm^2/s")
    if PLOTS and curves:
        fig_structure(curves)


# ---- figures ---------------------------------------------------------------------------------------

def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#555555"]


def fig_block(e):
    plt = _plt()
    ns = np.unique(np.logspace(np.log10(2), np.log10(len(e) // 2), 25).astype(int))
    se = [block_se(e, n) for n in ns]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.semilogx(len(e) / ns * 10, se, "o-", color=COLORS[0], label="block average")
    ax.axhline(e.std(ddof=1) / np.sqrt(len(e)), color=COLORS[1], ls="--", label="naive sigma/sqrt(M)")
    ax.axvline(len(e) / 5 * 10, color="#999999", ls=":", label="5 blocks (stb-mlmd)")
    ax.set_xlabel("block length (fs)")
    ax.set_ylabel("standard error of <E_pot> (meV/atom)")
    ax.set_title("Block averaging, NVT 300 K, 216 atoms")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "block_average.png"), dpi=130)
    plt.close(fig)


def fig_timestep(curves):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for (dt, t, de, fl), c in zip(curves, COLORS):
        if dt in (1.0, 3.0, 5.0):
            ax[0].plot(t / 1000, de, color=c, lw=1, label=f"dt = {dt:g} fs", zorder=10 - dt)
    ax[0].set_xlabel("t (ps)")
    ax[0].set_ylabel("E_tot(t) − E_tot(0) (meV/atom)")
    ax[0].set_title("NVE, 216 atoms, from 300 K")
    ax[0].legend(fontsize=8)
    dts = np.array([c[0] for c in curves])
    fls = np.array([c[3] for c in curves], float)
    ax[1].loglog(dts, fls, "o-", color=COLORS[0], label="sigma(E_tot), measured")
    ax[1].loglog(dts, fls[np.argmin(abs(dts - 1))] * dts ** 2, "--", color=COLORS[1], label="∝ dt²")
    ax[1].set_xlabel("dt (fs)")
    ax[1].set_ylabel("sigma(E_tot) (meV/atom)")
    ax[1].set_title("The Verlet energy oscillation grows as dt²")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "timestep.png"), dpi=130)
    plt.close(fig)


def fig_equipartition(th):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    t = th[:, 0]
    ax[0].plot(t, th[:, 1], color=COLORS[0], lw=1.2)
    ax[0].axhline(600, color="#999999", ls=":", label="600 K (initial velocities)")
    ax[0].axhline(300, color=COLORS[1], ls="--", label="300 K = 600/2")
    ax[0].set_xlim(0, min(t[-1], 600))
    ax[0].set_xlabel("t (fs)")
    ax[0].set_ylabel("T (K)")
    ax[0].set_title("NVE from the perfect lattice: the temperature halves")
    ax[0].legend(fontsize=8)
    e0 = th[0, 4]
    ax[1].plot(t, (th[:, 2] - th[0, 2]) * 1000, color=COLORS[1], lw=1.2, label="E_pot − E_pot(0)")
    ax[1].plot(t, th[:, 3] * 1000, color=COLORS[0], lw=1.2, label="E_kin")
    ax[1].plot(t, (th[:, 4] - e0) * 1000 + th[0, 3] * 1000, color="#333333", lw=1.2, label="E_tot (shifted)")
    ax[1].set_xlim(0, min(t[-1], 600))
    ax[1].set_xlabel("t (fs)")
    ax[1].set_ylabel("meV/atom")
    ax[1].set_title("Kinetic energy flows into potential energy")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "equipartition.png"), dpi=130)
    plt.close(fig)


def fig_size(pts):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    ns = np.logspace(np.log10(6), np.log10(700), 100)
    ax[0].loglog(ns, np.sqrt(2 / (3 * ns - 3)), color=COLORS[1], ls="--", label="sqrt(2 / (3N − 3))")
    ax[0].loglog([p[0] for p in pts], [p[1] for p in pts], "o", color=COLORS[0], ms=8, label="measured sigma_T / T")
    ax[0].set_xlabel("N (atoms)")
    ax[0].set_ylabel("sigma_T / T")
    ax[0].set_title("Canonical temperature fluctuation")
    ax[0].legend(fontsize=8)
    for (n, _, th), c in zip(pts, COLORS):
        if th is not None:
            ax[1].hist(th[:, 1], bins=40, density=True, histtype="step", lw=1.5, color=c, label=f"N = {n}")
    ax[1].set_xlabel("T (K)")
    ax[1].set_ylabel("probability density")
    ax[1].set_title("Distribution of the instantaneous temperature at 300 K")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "finite_size.png"), dpi=130)
    plt.close(fig)


def fig_expansion(pts, a300):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    ts = np.array([p[0] for p in pts], float)
    aa = np.array([p[1] for p in pts])
    se = np.array([p[2] for p in pts])
    tt = np.linspace(150, 1600, 200)
    ax[0].errorbar(ts, aa, yerr=se, fmt="o-", color=COLORS[0], label="MD, fine-tuned model (PBE)")
    ax[0].plot(tt, [A_EXP_273 * np.exp(okada_strain(273.2, x)) for x in tt], color=COLORS[1], label="experiment [17]")
    ax[0].axhline(5.48024, color="#999999", ls=":", label="SIESTA static (0 K, no zero-point)")
    ax[0].set_xlabel("T (K)")
    ax[0].set_ylabel("a (Å)")
    ax[0].set_title("Lattice constant: this PBE setup is ~1 % too large at every T")
    ax[0].legend(fontsize=8)
    ax[1].plot(ts, (aa / a300 - 1) * 100, "o-", color=COLORS[0], label="MD")
    ax[1].plot(tt, [(np.exp(okada_strain(300, x)) - 1) * 100 for x in tt], color=COLORS[1], label="experiment [17]")
    ax[1].set_xlabel("T (K)")
    ax[1].set_ylabel("a(T)/a(300 K) − 1 (%)")
    ax[1].set_title("Thermal expansion relative to 300 K")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "thermal_expansion.png"), dpi=130)
    plt.close(fig)


def fig_barostat(curves):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for (label, t, v), c in zip(curves, COLORS):
        ax.plot(t / 1000, v, color=c, lw=1, label=f"compressibility {label}")
    ax.set_xlabel("t (ps)")
    ax.set_ylabel("V (Å³/atom)")
    ax.set_title("Berendsen barostat, 216 Si at 300 K, 0 GPa")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "barostat.png"), dpi=130)
    plt.close(fig)


def fig_structure(curves):
    plt = _plt()
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.6))
    for (n, rdf, msd), c in zip(curves, COLORS):
        if rdf is not None:
            ax[0, 0].plot(rdf[:, 0], rdf[:, 1], color=c, lw=1.3, label=n)
        if msd is not None:
            ax[0, 1].plot(msd[:, 0] / 1000, msd[:, 1], color=c, lw=1.3, label=n)
    ax[0, 0].set_xlim(1.5, 8)
    ax[0, 0].set_xlabel("r (Å)")
    ax[0, 0].set_ylabel("g(r)")
    ax[0, 0].set_title("Radial distribution function")
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].set_xlabel("lag (ps)")
    ax[0, 1].set_ylabel("MSD (Å²)")
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_title("Mean-square displacement: plateau (solid) vs growth (liquid)")
    ax[0, 1].legend(fontsize=8)
    for n, c, label in [("melt_ft", COLORS[2], "fine-tuned"), ("melt_foundation", COLORS[3], "MACE-MP medium")]:
        th = data(n, "thermo")
        if th is None:
            continue
        ax[1, 0].plot(th[:, 0] / 1000, th[:, 6], color=c, lw=1.2, label=label)
        if n == "melt_ft":
            ax[1, 1].plot(th[:, 0] / 1000, th[:, 2], color=c, lw=1.2, label="fine-tuned, 2500 K")
    ax[1, 0].axhspan(18.722, 22.425, color="#cccccc", alpha=0.5, label="trained range (card)")
    ax[1, 0].set_xlabel("t (ps, production)")
    ax[1, 0].set_ylabel("V (Å³/atom)")
    ax[1, 0].set_title("Volume at 2500 K, 0 GPa")
    ax[1, 0].legend(fontsize=8)
    ax[1, 1].axhline(-115.0069, color="#555555", ls="--", label="lowest training energy (card)")
    th = data("nvt300", "thermo")
    if th is not None:
        ax[1, 1].axhline(th[:, 2].mean(), color=COLORS[0], ls=":", label="crystal at 300 K (nvt300)")
    ax[1, 1].set_xlabel("t (ps, production)")
    ax[1, 1].set_ylabel("E_pot (eV/atom, SIESTA scale)")
    ax[1, 1].set_title("The fine-tuned model sinks below its own crystal")
    ax[1, 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "structure_melt.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    s_reference()
    s_siesta()
    s_timestep()
    s_equipartition()
    s_size()
    s_expansion()
    s_melt()
    if PLOTS:
        print("\nfigures/ updated: block_average, timestep, equipartition, finite_size, thermal_expansion, barostat, structure_melt (.png)")
