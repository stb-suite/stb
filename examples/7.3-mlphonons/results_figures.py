#!/usr/bin/env python3
"""Reads reference_results/ (written by make_reference.sh) and prints the tables of README section 5, redrawing its
figures in figures/. Needs only numpy and matplotlib (no MACE, no phonopy).

    python3 results_figures.py            tables + figures
    python3 results_figures.py --no-plots tables only

Cases read (each a folder of reference_results/ holding stb_mlphonons_report.txt and the --save-data files): base,
base_444, foundation_{small,medium,large}, conv, qha, qha_foundation, qha_foundation_v0, dist_*, rattled_{norelax,
relax,tight}, elastic (stb_mlelastic_report.txt); and siesta_compare.json (siesta_compare.py). A missing case is
skipped with a note.
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
PLOTS = "--no-plots" not in sys.argv
KB_MEV = 0.08617333262          # meV/K
THZ_MEV = 4.135667696
R_GAS = 8.314462618
AMU = 1.66053906660e-27

# experiment
EXP_POINTS = {"Γ": [0, 0, 0, 15.53, 15.53, 15.53], "X": [4.49, 4.49, 12.32, 12.32, 13.90, 13.90],
              "L": [3.43, 3.43, 11.35, 12.60, 14.68, 14.68]}          # Nilsson & Nelin, PRB 6, 3777 (1972), 296 K
EXP_ELASTIC = {"C11": 165.64, "C12": 63.94, "C44": 79.51}            # GPa, Hall, Phys. Rev. 161, 756 (1967), 298 K
EXP_RHO = 2329.0                                                       # kg/m^3 at 300 K


def alpha_exp(T):
    """Okada & Tokumaru, J. Appl. Phys. 56, 314 (1984), 120-1500 K: linear expansion, 1/K."""
    T = np.asarray(T, float)
    return (3.725 * (1 - np.exp(-5.88e-3 * (T - 124))) + 5.548e-4 * T) * 1e-6


def note(msg):
    print(f"  ({msg})")


def report(name, fname="stb_mlphonons_report.txt"):
    p = os.path.join(REF, name, fname)
    return re.sub(r"\x1b\[[0-9;]*m", "", open(p).read()) if os.path.isfile(p) else None


def grab(text, pattern, cast=float, group=1):
    m = re.search(pattern, text or "", re.M)
    if not m:
        return None
    try:
        return cast(m.group(group))
    except (TypeError, ValueError):
        return None


def points(text, label="[^|]*"):
    """{q label: [frequencies]} from section [4] of a single-volume report (the first model listed, or `label`)."""
    out, cur = {}, None
    for line in (text or "").splitlines():
        m = re.match(r"^(\S*)\s*\|\s*(\[[^\]]*\]|)\s*\|\s*([^|]+?)\s*\|\s*(.+)$", line)
        if not m or "..." in line:
            continue
        q, model, freqs = m.group(1), m.group(3), m.group(4)
        if q:
            cur = q
        if cur and re.fullmatch(label, model) and cur not in out:
            try:
                out[cur] = [float(x) for x in freqs.split()]
            except ValueError:
                pass
    return out


def velocities(text, label):
    """{direction: (TA1, TA2, LA)} of the sound-velocity table that follows '[label] Sound velocities'."""
    out = {}
    i = (text or "").find(f"[{label}] Sound velocities")
    if i < 0:
        return out
    for line in text[i:].splitlines()[3:]:
        m = re.match(r"^Γ→(\S+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)", line)
        if not m:
            break
        out[m.group(1)] = tuple(float(m.group(k)) for k in (2, 3, 4))
    return out


def table(header, rows, fmt=None):
    w = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(header)]
    print("  " + " | ".join(str(h).ljust(w[i]) for i, h in enumerate(header)))
    print("  " + "-+-".join("-" * x for x in w))
    for r in rows:
        print("  " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(r)))


# ---- 5.2 the reference run --------------------------------------------------------------------------

def section_bands():
    print("\n=== 5.2 Frequencies at the special points vs experiment (THz) ===")
    base = report("base")
    if base is None:
        return note("base missing")
    ft = points(base, "Fine-tuned")
    fd = points(base, r"Foundation.*")
    rows = []
    for q in ("Γ", "X", "L"):
        e = EXP_POINTS[q]
        f = ft.get(q, [])
        g = fd.get(q, [])
        rows.append([q, " ".join(f"{x:5.2f}" for x in e[3 if q == "Γ" else 0:]),
                     " ".join(f"{x:5.2f}" for x in f[3 if q == "Γ" else 0:]),
                     " ".join(f"{x:5.2f}" for x in g[3 if q == "Γ" else 0:])])
    table(["q", "experiment", "fine-tuned", "MACE-MP-0 medium"], rows)
    pairs = []
    for q in ("Γ", "X", "L"):
        e, f = np.array(EXP_POINTS[q]), np.array(ft[q])
        k = e > 0
        pairs += list(zip(f[k], e[k]))
    pairs = np.array(pairs)
    d = (pairs[:, 0] - pairs[:, 1]) / pairs[:, 1] * 100
    print(f"  fine-tuned vs experiment over {len(d)} modes: mean {d.mean():+.1f} %, RMS {np.sqrt(np.mean(d ** 2)):.1f} %, "
          f"range {d.min():+.1f} .. {d.max():+.1f} %")
    return ft, fd


def bands_figure(siesta):
    base = os.path.join(REF, "base")
    gp = os.path.join(base, "plot", "bands.gplot")
    if not (os.path.isfile(os.path.join(base, "bands.dat")) and os.path.isfile(gp)):
        return note("base/bands.dat or plot/bands.gplot missing: no band figure")
    import matplotlib.pyplot as plt
    ticks = re.findall(r'"([^"]+)" ([0-9.]+)', re.search(r"set xtics \((.*)\)", open(gp).read()).group(1))
    ticks = [(lab.replace("{/Symbol G}", "Γ"), float(x)) for lab, x in ticks]
    fig, ax = plt.subplots(figsize=(9, 5.8))
    for fname, color, lab in (("bands.dat", "#2a6fdb", "fine-tuned (SIESTA PBE/DZP level)"),
                              ("bands_foundation.dat", "#d9480f", "MACE-MP-0 medium")):
        blocks = open(os.path.join(base, fname)).read().strip().split("\n\n")
        first = True
        for b in blocks:
            a = np.loadtxt([l for l in b.splitlines() if not l.startswith("#")])
            for col in range(1, a.shape[1]):
                ax.plot(a[:, 0], a[:, col], color=color, lw=1.4, label=lab if first else None)
                first = False
    xpos = {}
    for lab, x in ticks:
        xpos.setdefault(lab, []).append(x)
    first = True
    for q, freqs in EXP_POINTS.items():
        for x in xpos.get(q, []):
            f = np.array(freqs)
            ax.plot([x] * len(f[f > 0]), f[f > 0], "ko", ms=6, mfc="white", mew=1.5, zorder=5,
                    label="experiment (neutrons, 296 K)" if first else None)
            first = False
    if siesta:
        for key_s, mfc, lab in (("SIESTA", "#2f9e44", "SIESTA (7.1 reference, 2x2x2)"),
                                ("SIESTA, one atom repaired", "white", "SIESTA, one bad force repaired (5.9)")):
            if key_s not in siesta["phonons"]:
                continue
            sp = siesta["phonons"][key_s]["points"]
            first = True
            for q, key in (("Γ", "G"), ("X", "X"), ("L", "L")):
                for x in xpos.get(q, []):
                    f = np.array(sp[key])
                    ax.plot([x + 0.006] * len(f[f > 0.5]), f[f > 0.5], "s", color="#2f9e44", mfc=mfc, ms=5, zorder=4,
                            label=lab if first else None)
                    first = False
    for lab, x in ticks:
        ax.axvline(x, color="grey", lw=0.5)
    ax.set_xticks([x for _, x in ticks])
    ax.set_xticklabels([lab for lab, _ in ticks])
    ax.set_xlim(ticks[0][1], ticks[-1][1])
    ax.set_ylim(0, 17)
    ax.set_ylabel("frequency (THz)")
    ax.set_title("Silicon phonons: fine-tuned model, foundation model, SIESTA and experiment")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "bands_vs_experiment.png"), dpi=150)
    plt.close(fig)


# ---- 5.3 sound velocities and elastic constants -----------------------------------------------------

def section_elastic():
    print("\n=== 5.3 Sound velocities (m/s) and elastic constants (GPa) ===")
    r2, r4 = report("base"), report("base_444")
    rho = None
    if r2:
        rho = grab(r2, r"density ([0-9.]+) g/cm\^3") * 1000
    rows = []
    vel = {}
    for name, txt in (("2x2x2", r2), ("4x4x4", r4)):
        v = velocities(txt, "Fine-tuned")
        if v:
            vel[name] = v
            rows.append([name] + [f"{v[d][0]:.0f} / {v[d][2]:.0f}" for d in ("X", "K", "L")])
    vf = velocities(r2, "Foundation (medium)")
    if vf:
        rows.append(["foundation, 2x2x2"] + [f"{vf[d][0]:.0f} / {vf[d][2]:.0f}" for d in ("X", "K", "L")])
    c11, c12, c44 = (EXP_ELASTIC[k] * 1e9 for k in ("C11", "C12", "C44"))
    ve = {"X": (np.sqrt(c44 / EXP_RHO), np.sqrt(c11 / EXP_RHO)),
          "K": (np.sqrt((c11 - c12) / 2 / EXP_RHO), np.sqrt((c11 + c12 + 2 * c44) / 2 / EXP_RHO)),
          "L": (np.sqrt((c11 - c12 + c44) / 3 / EXP_RHO), np.sqrt((c11 + 2 * c12 + 4 * c44) / 3 / EXP_RHO))}
    rows.append(["experiment (from C_ij)"] + [f"{ve[d][0]:.0f} / {ve[d][1]:.0f}" for d in ("X", "K", "L")])
    table(["", "Γ→X [100] TA / LA", "Γ→K [110] TA1 / LA", "Γ→L [111] TA / LA"], rows)
    if "4x4x4" in vel and rho:
        v = vel["4x4x4"]
        C11 = rho * v["X"][2] ** 2 / 1e9
        C44 = rho * v["X"][0] ** 2 / 1e9
        C12 = 2 * rho * v["K"][2] ** 2 / 1e9 - C11 - 2 * C44
        c11mc12 = 2 * rho * v["K"][0] ** 2 / 1e9
        B = (C11 + 2 * C12) / 3
        print(f"  from the 4x4x4 velocities (rho = {rho:.1f} kg/m^3): C11 {C11:.1f}, C12 {C12:.1f}, C44 {C44:.1f} GPa; "
              f"check C11-C12 = {C11 - C12:.1f} vs 2 rho v_TA1[110]^2 = {c11mc12:.1f}; B = (C11+2C12)/3 = {B:.1f} GPa")
    el = report("elastic", "stb_mlelastic_report.txt")
    if el:
        vals = {}
        for k in ("C11", "C12", "C44"):
            v = grab(el, rf"\b{k}\s*[:=|]\s*([0-9.]+)")
            vals[k] = v
        b = grab(el, r"Bulk modulus.*?\b([0-9]+\.[0-9]+)")
        print(f"  stb-mlelastic (stress-strain, same model): " + ", ".join(f"{k} {v}" for k, v in vals.items())
              + (f"; B {b}" if b else ""))
    else:
        note("elastic missing")
    print(f"  experiment (Hall 1967): C11 {EXP_ELASTIC['C11']}, C12 {EXP_ELASTIC['C12']}, C44 {EXP_ELASTIC['C44']} GPa; "
          f"B = {(EXP_ELASTIC['C11'] + 2 * EXP_ELASTIC['C12']) / 3:.1f} GPa")


# ---- 5.4 DOS and thermal properties -----------------------------------------------------------------

def section_thermal():
    print("\n=== 5.4 DOS and thermal properties ===")
    base = report("base")
    if base is None:
        return note("base missing")
    for lab in ("Fine-tuned", "Foundation (medium)"):
        line = grab(base, rf"^\[{re.escape(lab)}\] Max frequency (.*)$", str)
        print(f"  [{lab}] {line}")
    th = os.path.join(REF, "base", "thermal.dat")
    if os.path.isfile(th):
        a = np.loadtxt(th)
        rows = []
        for T in (100, 300, 500, 1000):
            i = int(np.argmin(np.abs(a[:, 0] - T)))
            rows.append([f"{a[i, 0]:.0f}", f"{a[i, 7]:.3f}", f"{a[i, 7] * R_GAS:.2f}", f"{a[i, 6]:.3f}",
                         f"{a[i, 6] * R_GAS:.2f}", f"{a[i, 5]:.2f}"])
        table(["T (K)", "Cv kB/atom", "Cv J/K/mol", "S kB/atom", "S J/K/mol", "F meV/atom"], rows)
        print("  experiment at 298 K: Cp 20.00 J/K/mol (Cv ~ 19.98), S 18.81 J/K/mol; Debye temperature 645 K")
    if PLOTS and os.path.isfile(th):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        for path, color, lab in ((os.path.join(REF, "base", "dos.dat"), "#2a6fdb", "fine-tuned"),
                                 (os.path.join(REF, "base", "plot", "dos_1.dat"), "#d9480f", "MACE-MP-0 medium")):
            if os.path.isfile(path):
                d = np.loadtxt(path)
                axes[0].plot(d[:, 0], d[:, 1], color=color, lw=1.6, label=lab)
                axes[0].fill_between(d[:, 0], d[:, 1], color=color, alpha=0.12)
        axes[0].set_xlabel("frequency (THz)")
        axes[0].set_ylabel("DOS (states/THz per 2-atom cell)")
        axes[0].set_xlim(0, 17)
        axes[0].set_ylim(0, None)
        axes[0].legend(fontsize=8)
        axes[0].set_title("phonon density of states")
        a = np.loadtxt(th)
        axes[1].plot(a[:, 0], a[:, 7], color="#2a6fdb", lw=2, label="fine-tuned, harmonic (quantum)")
        p1 = os.path.join(REF, "base", "plot", "thermal_1.dat")
        if os.path.isfile(p1):
            b = np.loadtxt(p1)
            axes[1].plot(b[:, 0], b[:, 7], color="#d9480f", lw=1.6, label="MACE-MP-0 medium")
        T = np.linspace(1, a[-1, 0], 300)
        xd = 645.0 / T
        debye = []
        for x in xd:
            t = np.linspace(1e-6, min(x, 60.0), 400)        # the integrand is ~1e-20 beyond 60
            debye.append(9 / x ** 3 * np.trapezoid(t ** 4 * np.exp(-t) / (-np.expm1(-t)) ** 2, t))
        axes[1].plot(T, debye, "k:", lw=1.3, label="Debye model, θ_D = 645 K")
        axes[1].axhline(3, color="grey", ls="--", lw=1, label="classical (Dulong-Petit, MD)")
        axes[1].plot([298.15], [19.98 / R_GAS], "ks", ms=7, label="experiment, 298 K (Cv)")
        axes[1].set_xlabel("T (K)")
        axes[1].set_ylabel("C_V (k_B per atom)")
        axes[1].set_title("heat capacity")
        axes[1].legend(fontsize=8, loc="lower right")
        axes[1].set_ylim(0, 3.2)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "dos_heat_capacity.png"), dpi=150)
        plt.close(fig)


# ---- 5.5 foundation models --------------------------------------------------------------------------

def section_foundation(ft):
    print("\n=== 5.5 The foundation models (THz) ===")
    rows, data = [], {}
    for m in ("small", "medium", "large"):
        r = report(f"foundation_{m}")
        if r is None:
            note(f"foundation_{m} missing")
            continue
        p = points(r)
        data[m] = p
        v = velocities(r, f"MACE-MP-0 ({m})")
        rows.append([f"MACE-MP-0 {m}", f"{p['Γ'][3]:.2f}", f"{p['X'][0]:.2f}", f"{p['X'][4]:.2f}", f"{p['L'][0]:.2f}",
                     f"{p['L'][4]:.2f}", f"{v['X'][2]:.0f}" if v else "-",
                     f"{grab(r, r'Debye temperature ([0-9]+) K'):.0f}"])
    if ft:
        b = report("base")
        v = velocities(b, "Fine-tuned")
        rows.append(["fine-tuned", f"{ft['Γ'][3]:.2f}", f"{ft['X'][0]:.2f}", f"{ft['X'][4]:.2f}", f"{ft['L'][0]:.2f}",
                     f"{ft['L'][4]:.2f}", f"{v['X'][2]:.0f}", f"{grab(b, r'Debye temperature ([0-9]+) K'):.0f}"])
    rows.append(["experiment", "15.53", "4.49", "13.90", "3.43", "14.68", f"{np.sqrt(EXP_ELASTIC['C11'] * 1e9 / EXP_RHO):.0f}",
                 "645"])
    table(["model", "Γ optical", "X TA", "X TO", "L TA", "L TO", "v_LA[100] m/s", "θ_D K"], rows)
    if PLOTS and data:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 4.2))
        labels = ["Γ opt", "X TA", "X LA/LO", "X TO", "L TA", "L LA", "L LO", "L TO"]
        idx = [("Γ", 3), ("X", 0), ("X", 2), ("X", 4), ("L", 0), ("L", 2), ("L", 3), ("L", 4)]
        exp = np.array([EXP_POINTS[q][i] for q, i in idx])
        series = [(f"MACE-MP-0 {m}", data[m], c) for m, c in (("small", "#f59f00"), ("medium", "#d9480f"),
                                                                 ("large", "#862e9c")) if m in data]
        if ft:
            series.append(("fine-tuned", ft, "#2a6fdb"))
        x = np.arange(len(idx))
        for k, (lab, p, c) in enumerate(series):
            f = np.array([p[q][i] for q, i in idx])
            ax.bar(x + (k - (len(series) - 1) / 2) * 0.2, (f - exp) / exp * 100, width=0.2, color=c, label=lab)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("frequency vs experiment (%)")
        ax.set_title("Foundation models: optical modes ~30 % soft, zone-boundary TA modes nearly right", fontsize=10)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "foundation_softening.png"), dpi=150)
        plt.close(fig)


# ---- 5.6 convergence --------------------------------------------------------------------------------

def section_convergence():
    print("\n=== 5.6 Supercell convergence ===")
    r = report("conv")
    p = os.path.join(REF, "conv", "convergence.dat")
    if r is None or not os.path.isfile(p):
        return note("conv missing")
    rows = []
    for line in open(p):
        if line.startswith("#"):
            continue
        f = line.split()
        rows.append([f[0], f[1], f"{float(f[6]):.3f}", f"{float(f[7]):.3f}", f"{float(f[10]):.2f}", f"{float(f[8]):.4f}",
                     f"{float(f[9]):.4f}"])
    table(["supercell", "atoms", "RMSE vs 4x4x4 THz", "max |d| THz", "ZPE meV/at", "Cv(300) kB/at", "S(300) kB/at"], rows)
    print(grab(r, r"^(Converged .*)$", str))
    pts = {}
    for line in r.splitlines():
        m = re.match(r"^(\S*)\s*\|\s*(\dx\dx\d)\s*\|\s*(.+?)\s*\|\s*(\S+)\s*$", line)
        if m:
            q = m.group(1) or q
            pts.setdefault(q, {})[m.group(2)] = [float(x) for x in m.group(3).split()]
    if PLOTS and pts:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        a = np.array([[float(v) for v in l.split()[1:]] for l in open(p) if not l.startswith("#")])
        n = a[:, 0]
        axes[0].semilogy(n[:-1], a[:-1, 5], "o-", color="#2a6fdb", label="band-path RMSE vs 4x4x4")
        axes[0].semilogy(n[:-1], a[:-1, 6], "s--", color="#d9480f", label="largest difference")
        axes[0].axhline(0.1, color="grey", ls=":", lw=1)
        axes[0].text(n[0], 0.11, "tolerance 0.1 THz", fontsize=8, color="grey")
        axes[0].set_xscale("log")
        axes[0].set_xticks(n[:-1])
        axes[0].set_xticklabels([f"{int(x)}" for x in n[:-1]])
        axes[0].set_xlabel("atoms in the supercell")
        axes[0].set_ylabel("THz")
        axes[0].legend(fontsize=8)
        axes[0].set_title("convergence with the supercell")
        dims = ["1x1x1", "2x2x2", "3x3x3", "4x4x4"]
        for q, c in (("L", "#2a6fdb"), ("K", "#2f9e44"), ("X", "#868e96")):
            if q not in pts:
                continue
            f = np.array([pts[q][d] for d in dims])
            for b in range(f.shape[1]):
                axes[1].plot(range(4), f[:, b], "o-", color=c, lw=1, ms=4, label=q if b == 0 else None)
        axes[1].set_xticks(range(4))
        axes[1].set_xticklabels(dims)
        axes[1].set_ylabel("frequency (THz)")
        axes[1].set_title("X never changes; L is exact from 2x2x2; K converges")
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "convergence.png"), dpi=150)
        plt.close(fig)


# ---- 5.7 displacement -------------------------------------------------------------------------------

def section_displacement():
    print("\n=== 5.7 Displacement amplitude (THz) ===")
    rows, data = [], []
    for p in sorted(glob.glob(os.path.join(REF, "dist_*")), key=lambda s: float(s.split("_")[-1])):
        u = float(p.split("_")[-1])
        pt = points(report(os.path.basename(p)))
        if not pt:
            continue
        data.append((u, pt["Γ"][3], pt["X"][0], pt["L"][0], pt["X"][4]))
        rows.append([f"{u:g}", f"{pt['Γ'][3]:.2f}", f"{pt['X'][0]:.2f}", f"{pt['L'][0]:.2f}", f"{pt['X'][4]:.2f}"])
    if not rows:
        return note("dist_* missing")
    table(["u (Ang)", "Γ optical", "X TA", "L TA", "X TO"], rows)
    if PLOTS:
        import matplotlib.pyplot as plt
        d = np.array(data)
        ref = d[np.argmin(np.abs(d[:, 0] - 0.01))]
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for k, (lab, c) in enumerate((("Γ optical", "#d9480f"), ("X TA", "#2a6fdb"), ("L TA", "#2f9e44"),
                                      ("X TO", "#862e9c")), start=1):
            ax.semilogx(d[:, 0], (d[:, k] / ref[k] - 1) * 100, "o-", color=c, label=lab)
        ax.axvline(0.01, color="grey", ls=":", lw=1)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("displacement (Å)")
        ax.set_ylabel("change vs 0.01 Å (%)")
        ax.set_title("Anharmonicity shows from ~0.05 Å; the soft TA modes first")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "displacement.png"), dpi=150)
        plt.close(fig)


# ---- 5.8 rattled ------------------------------------------------------------------------------------

def section_rattled():
    print("\n=== 5.8 A structure that is not at a minimum ===")
    rows = []
    for c in ("rattled_norelax", "rattled_relax", "rattled_tight"):
        r = report(c)
        if r is None:
            note(f"{c} missing")
            continue
        relax = grab(r, r"^Pre-relax\s*: (.*)$", str) or ""
        f1 = re.search(r"max \|F\| ([0-9.]+) eV/Ang after (\d+)", relax)
        shift = grab(relax, r"shift ([0-9.]+) Ang")
        g = [float(x) for x in (re.search(r"^Results.*Gamma optical ([0-9., ]+) THz", r, re.M).group(1)).split(",")]
        g_top = [x for x in g if x > 14.0]
        rows.append([c.replace("rattled_", ""), f"{grab(r, r'^Input forces\s*: max \|F\| ([0-9.]+)'):.3f}",
                     f1.group(1) if f1 else "-", f1.group(2) if f1 else "-", f"{shift:.4f}" if shift else "-",
                     grab(r, r"^Symmetry\s*: (\S+)", str), grab(r, r"^Displacements\s*: (\d+)", str),
                     f"{min(g_top):.2f} .. {max(g_top):.2f}" if g_top else "-",
                     str(grab(r, r"imaginary q (\d+)/", int))])
    table(["pre-relax", "max|F| in", "max|F| out", "steps", "max shift Ang", "space group", "displ.",
           "Γ optical spread THz", "imag. q"], rows)


# ---- 5.9 SIESTA -------------------------------------------------------------------------------------

def section_siesta():
    print("\n=== 5.9 Against SIESTA ===")
    p = os.path.join(REF, "siesta_compare.json")
    if not os.path.isfile(p):
        note("siesta_compare.json missing (python3 siesta_compare.py)")
        return None
    s = json.load(open(p))
    rows = []
    for name, r in s["phonons"].items():
        pt = r["points"]
        e = r["rmse_vs_siesta"]
        rows.append([name, f"{pt['G'][3]:.2f}", f"{pt['X'][0]:.2f}", f"{pt['X'][4]:.2f}", f"{pt['L'][0]:.2f}",
                     f"{r['min']:.2f}", f"{r['imaginary_q']}/{r['n_q']}", f"{e['acoustic']:.2f} / {e['optical']:.2f}"])
    rows.append(["experiment", "15.53", "4.49", "13.90", "3.43", "-", "0", "-"])
    table(["", "Γ opt", "X TA", "X TO", "L TA", "min f", "imaginary q", "RMSE vs SIESTA ac./opt."], rows)
    print("  forces on the displaced supercell, RMS per component by neighbour shell (meV/Ang):")
    table(["r (Ang)", "atoms", "SIESTA", "model", "difference"],
          [[f"{x['r']:.3f}", x["n"], f"{x['siesta_rms']:.2f}", f"{x['model_rms']:.2f}", f"{x['diff_rms']:.2f}"]
           for x in s["force_shells_meV_per_ang"][:8]])
    sc = s["symmetry_scatter_meV_per_ang"]
    print(f"  symmetry scatter beyond 3 Ang ({sc['n_site_operations']} site operations): SIESTA {sc['SIESTA']:.2f}, "
          f"fine-tuned {sc['fine-tuned']:.4f}, foundation {sc['MACE-MP-0 medium']:.4f} meV/Ang; true force there "
          f"{sc['model_force_rms_beyond_3A']:.2f} (model) / {sc['SIESTA_true_force_rms_beyond_3A']:.2f} (SIESTA)")
    o = s.get("outlier")
    if o:
        print(f"  net force in SIESTA's .FA (meV/Ang): displaced supercell {o['net_force_meV_per_ang']['displaced']}, "
              f"perfect supercell {o['net_force_meV_per_ang']['perfect']}")
        print(f"  outlier: atom {o['atom']} at {o['distance_ang']:.2f} Ang from the displaced atom, force "
              f"{o['force_meV_per_ang']} meV/Ang; median of its {11} symmetry partners {o['partners_median_meV_per_ang']:.2f}")
    print("  noise experiment (model forces + Gaussian noise, 20 draws):")
    table(["u (Ang)", "sigma meV/Ang", "share with imag.", "median min f", "worst min f", "X TA"],
          [[f"{x['u_ang']:.2f}", f"{x['sigma_meV']:.1f}", f"{x['share_with_imaginary']:.2f}", f"{x['min_f_median']:.2f}",
            f"{x['min_f_worst']:.2f}", f"{x['X_TA_mean']:.2f} ± {x['X_TA_std']:.2f}"] for x in s["noise_experiment"]])
    if PLOTS:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        sh = s["force_shells_meV_per_ang"][1:12]
        r = [x["r"] for x in sh]
        axes[0].semilogy(r, [x["siesta_rms"] for x in sh], "s-", color="#2f9e44", label="SIESTA force")
        axes[0].semilogy(r, [x["model_rms"] for x in sh], "o-", color="#2a6fdb", label="model force")
        axes[0].semilogy(r, [x["diff_rms"] for x in sh], "x--", color="#d9480f", label="SIESTA − model")
        axes[0].axhline(sc["SIESTA"], color="#2f9e44", ls=":", lw=1.2, label="SIESTA's symmetry scatter (noise)")
        axes[0].set_xlabel("distance from the displaced atom (Å)")
        axes[0].set_ylabel("RMS force component (meV/Å)")
        axes[0].set_title("SIESTA vs model forces: one shell (5.97 Å) stands out", fontsize=10)
        axes[0].legend(fontsize=8)
        for u, c in ((0.01, "#d9480f"), (0.03, "#2a6fdb")):
            xs = [x for x in s["noise_experiment"] if x["u_ang"] == u]
            axes[1].plot([x["sigma_meV"] for x in xs], [x["share_with_imaginary"] for x in xs], "o-", color=c,
                         label=f"displacement {u} Å")
        axes[1].axvline(sc["SIESTA"], color="#2f9e44", ls=":", lw=1.2)
        axes[1].text(sc["SIESTA"] + 0.05, 0.75, "SIESTA's symmetry\nscatter", fontsize=8, color="#2f9e44")
        axes[1].set_xlabel("noise added to the model's forces (meV/Å)")
        axes[1].set_ylabel("share of draws with imaginary modes")
        axes[1].set_title("Random force errors / displacement -> imaginary modes", fontsize=10)
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "siesta_noise.png"), dpi=150)
        plt.close(fig)
    return s


# ---- 5.10 QHA ---------------------------------------------------------------------------------------

MD_ALPHA = [(450, 4.87), (750, 4.46), (1050, 4.59), (1350, 3.48)]     # 7.2 example, NPT MD, 10^-6/K (classical)


def qha_data(name):
    p = os.path.join(REF, name, "qha_thermal_expansion.dat")
    return np.loadtxt(p) if os.path.isfile(p) else None


def section_qha():
    print("\n=== 5.10 Quasi-harmonic approximation ===")
    rows = []
    for name, lab in (("qha", "fine-tuned"), ("qha_foundation", "MACE-MP-0 medium, input volumes"),
                      ("qha_foundation_v0", "MACE-MP-0 medium, own V0")):
        r = report(name)
        a = qha_data(name)
        if r is None or a is None:
            note(f"{name} missing")
            continue
        def at(T, col):
            return a[int(np.argmin(np.abs(a[:, 0] - T))), col]
        rows.append([lab, f"{grab(r, r'= ([0-9.]+) Ang\^3/atom \('):.4f}", f"{grab(r, r'^B0 \(bulk modulus\) : ([0-9.]+)'):.1f}",
                     f"{grab(r, r'\(\+([0-9.]+) %\)'):.2f}", f"{at(300, 4) * 1e6:.2f}", f"{at(1000, 4) * 1e6:.2f}",
                     f"{at(300, 5):.1f} / {at(1000, 5):.1f}", f"{at(300, 8):.2f}"])
    rows.append(["experiment", "20.02 (300 K)", "97.8", "-", f"{alpha_exp(300) * 1e6:.2f}", f"{alpha_exp(1000) * 1e6:.2f}",
                 "97.8 (298 K)", "0.46"])
    table(["", "V0 Ang^3/atom", "B0 GPa", "zero-point ΔV %", "α_L(300) 1e-6/K", "α_L(1000)", "B_T 300/1000 K GPa",
           "γ(300 K)"], rows)
    print("  experiment: V0 = a^3/8 with a = 5.4310 Ang at 300 K; B = (C11 + 2 C12)/3 from Hall; gamma = 3 alpha_L B V / C_V "
          "= 3 x 2.57e-6 x 97.8 GPa x 20.02 Ang^3 / 2.40 kB = 0.46")
    q = report("qha")
    if q:
        print("  mode Grueneisen parameters (fine-tuned):")
        for line in q.splitlines():
            if re.match(r"^(Γ|X|L) \|", line) and "+" in line:
                print("   ", line)
        print("  experiment: Γ-TO ~ +1.0; X-TA ~ -1.4; L-TA ~ -1.3 (Raman under pressure, Weinstein & Piermarini 1975)")
    if PLOTS:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
        T = np.linspace(125, 1500, 200)
        axes[0].plot(T, alpha_exp(T) * 1e6, "k-", lw=2, label="experiment (Okada & Tokumaru)")
        for name, c, lab in (("qha", "#2a6fdb", "QHA, fine-tuned"), ("qha_foundation_v0", "#d9480f", "QHA, MACE-MP-0 medium")):
            a = qha_data(name)
            if a is not None:
                axes[0].plot(a[:, 0], a[:, 4] * 1e6, color=c, lw=1.8, label=lab)
                axes[1].plot(a[:, 0], a[:, 5], color=c, lw=1.8, label=lab)
                axes[2].plot(a[:, 0], (a[:, 2] / a[0, 2]) ** (1 / 3) * 100 - 100, color=c, lw=1.8, label=lab)
        axes[0].plot([t for t, _ in MD_ALPHA], [v for _, v in MD_ALPHA], "s", color="#2f9e44", ms=7,
                     label="NPT MD, fine-tuned (7.2, classical)")
        axes[0].axhline(0, color="grey", lw=0.6)
        axes[0].set_xlabel("T (K)")
        axes[0].set_ylabel("α_L (10⁻⁶ K⁻¹)")
        axes[0].set_title("linear thermal expansion")
        axes[0].legend(fontsize=7)
        axes[0].set_ylim(-1, 12)
        axes[1].set_xlabel("T (K)")
        axes[1].set_ylabel("B_T (GPa)")
        axes[1].set_title("isothermal bulk modulus")
        axes[1].legend(fontsize=7)
        a = qha_data("qha")
        if a is not None:
            ins = axes[2].inset_axes([0.12, 0.55, 0.4, 0.38])
            k = a[:, 0] <= 200
            ins.plot(a[k, 0], a[k, 4] * 1e6, color="#2a6fdb")
            ins.axhline(0, color="grey", lw=0.5)
            ins.set_title("α_L below 200 K", fontsize=7)
            ins.tick_params(labelsize=6)
        axes[2].set_xlabel("T (K)")
        axes[2].set_ylabel("a(T)/a(0 K) − 1 (%)")
        axes[2].set_title("lattice constant")
        axes[2].legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "qha.png"), dpi=150)
        plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    siesta = section_siesta_quiet()
    ft_fd = section_bands()
    if PLOTS:
        bands_figure(siesta)
    section_elastic()
    section_thermal()
    section_foundation(ft_fd[0] if ft_fd else None)
    section_convergence()
    section_displacement()
    section_rattled()
    section_siesta()
    section_qha()
    if PLOTS:
        print("\nfigures/: bands_vs_experiment, dos_heat_capacity, foundation_softening, convergence, displacement, "
              "siesta_noise, qha (.png)")


def section_siesta_quiet():
    p = os.path.join(REF, "siesta_compare.json")
    return json.load(open(p)) if os.path.isfile(p) else None


if __name__ == "__main__":
    main()
