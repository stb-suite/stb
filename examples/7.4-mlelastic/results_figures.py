#!/usr/bin/env python3
"""Reads reference_results/ (written by make_reference.sh) and prints the tables of README section 5, redrawing its
summary figures in figures/. Needs only numpy and matplotlib (no MACE).

    python3 results_figures.py            tables + figures
    python3 results_figures.py --no-plots tables only

Cases read (each a folder of reference_results/ holding stb_mlelastic_report.txt and the --save-data files): base,
conv, fmax_*, foundation_{small,medium,large}, si216, rattled_noprerelax, rattled, mg_hcp, graphene. A missing case is
skipped with a note.
"""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference_results")
FIG = os.path.join(HERE, "figures")
PLOTS = "--no-plots" not in sys.argv

# experiment and other references
SI_EXP = {"C11": 165.64, "C12": 63.94, "C44": 79.51}        # GPa, Hall, Phys. Rev. 161, 756 (1967), 298 K
SI_RHO_EXP = 2329.0
SI_PHONON_52 = {"C11": 148.0, "C12": 70.2, "C44": 65.1}     # 7.3 example: from the 4x4x4 acoustic slopes
SI_PHONON_V = {"[100]": (8081, 5361, 5361), "[110]": (8768, 5365, 4122), "[111]": (8994, 4566, 4566)}
MG_EXP = {"C11": 59.4, "C12": 25.6, "C13": 21.4, "C33": 61.6, "C44": 16.4}   # room temperature (Simmons & Wang 1971)


def note(msg):
    print(f"  ({msg})")


def report(name):
    p = os.path.join(REF, name, "stb_mlelastic_report.txt")
    return re.sub(r"\x1b\[[0-9;]*m", "", open(p).read()) if os.path.isfile(p) else None


def section(text, title):
    out, on = [], False
    for line in (text or "").splitlines():
        if re.match(r"^\[\d+\] [A-Z]", line):     # '[100] | ...' rows of the sound table are not headers
            if on:
                break
            on = title in line
        if on:
            out.append(line)
    return "\n".join(out)


def grab(text, pattern, cast=float, group=1):
    m = re.search(pattern, text or "", re.M)
    if not m:
        return None
    try:
        return cast(m.group(group))
    except (TypeError, ValueError):
        return None


def tensor(text):
    """The first 6x6 stiffness matrix of the report (rows '  i=1 | ...')."""
    rows = re.findall(r"^  i=\d \|(.*)$", section(text, "STIFFNESS TENSOR"), re.M)[:6]
    return np.array([[float(x) for x in r.split()] for r in rows]) if len(rows) == 6 else None


def independent(text):
    """[(names, value)] from the 'Independent' list."""
    out = []
    for m in re.finditer(r"^  (C\d\d(?: = C\d\d)*)\s+(-?[\d.]+) ", section(text, "STIFFNESS TENSOR"), re.M):
        out.append((m.group(1), float(m.group(2))))
    return out


def hill(text):
    m = re.search(r"^Hill \(average\)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)", text or "", re.M)
    return tuple(float(x) for x in m.groups()) if m else None


def fit_rows(text, label):
    """{dir: (C, R2, intercept, max_shift, steps)} of one model's fit table."""
    sec = section(text, "STRESS-STRAIN FITS")
    i = sec.find(f"--- {label}")
    out = {}
    for line in sec[i:].splitlines()[3:]:
        m = re.match(r"^(\w\w)\s*\|\s*C\d\d\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([-+\d.]+)\s*\|\s*[\d.]+\s*\|\s*([\d.]+)"
                     r"\s*\|\s*[\d.]+\s*\|\s*(\d+)", line)
        if not m:
            break
        out[m.group(1)] = tuple(float(x) for x in m.groups()[1:])
    return out


def table(header, rows):
    w = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(header)]
    print("  " + " | ".join(str(h).ljust(w[i]) for i, h in enumerate(header)))
    print("  " + "-+-".join("-" * x for x in w))
    for r in rows:
        print("  " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(r)))


def vrh_cubic(c11, c12, c44):
    B = (c11 + 2 * c12) / 3
    gv = (c11 - c12 + 3 * c44) / 5
    gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
    G = (gv + gr) / 2
    return B, G, 9 * B * G / (3 * B + G), (3 * B - 2 * G) / (2 * (3 * B + G))


# ---- 5.2 the reference run ------------------------------------------------------------------

def section_base():
    print("\n=== 5.2 Silicon: constants and moduli (GPa) ===")
    t = report("base")
    if t is None:
        return note("base missing")
    C = tensor(t)
    print("  independent constants:", ", ".join(f"{n} {v:.2f}" for n, v in independent(t)))
    fr = fit_rows(t, "fine-tuned")
    fd = fit_rows(t, "MACE-MP-0 medium")
    rows = []
    cmp = section(t, "FINE-TUNED vs FOUNDATION")
    fc = {k: grab(cmp, rf"^{re.escape(k)}.*?\|\s*[\d.]+\s*\|\s*([\d.]+)") for k in ("C11 = C22 = C33", "C12 = C13 = C23",
                                                                                   "C44 = C55 = C66")}
    for key, name, exp, ph in (("C11 = C22 = C33", "C11", SI_EXP["C11"], SI_PHONON_52["C11"]),
                               ("C12 = C13 = C23", "C12", SI_EXP["C12"], SI_PHONON_52["C12"]),
                               ("C44 = C55 = C66", "C44", SI_EXP["C44"], SI_PHONON_52["C44"])):
        m, n = int(name[1]) - 1, int(name[2]) - 1
        v = C[m, n]
        rows.append([name, f"{v:.2f}", f"{ph:.1f}", f"{fc[key]:.2f}" if fc[key] else "-", f"{exp:.1f}",
                     f"{100 * (v / exp - 1):+.1f} %"])
    h = hill(t)
    hf = grab(cmp, r"^B Hill \(GPa\)\s*\|\s*[\d.]+\s*\|\s*([\d.]+)")
    gf = grab(cmp, r"^G Hill \(GPa\)\s*\|\s*[\d.]+\s*\|\s*([\d.]+)")
    eB, eG, eE, enu = vrh_cubic(SI_EXP["C11"], SI_EXP["C12"], SI_EXP["C44"])
    rows.append(["B (Hill)", f"{h[0]:.2f}", "93.5 (QHA B0)", fmt(hf), f"{eB:.1f}", f"{100 * (h[0] / eB - 1):+.1f} %"])
    rows.append(["G (Hill)", f"{h[1]:.2f}", "-", fmt(gf), f"{eG:.1f}", f"{100 * (h[1] / eG - 1):+.1f} %"])
    rows.append(["E (Hill)", f"{h[2]:.2f}", "-", "-", f"{eE:.1f}", f"{100 * (h[2] / eE - 1):+.1f} %"])
    rows.append(["Poisson ν", f"{h[3]:.3f}", "-", "-", f"{enu:.3f}", "-"])
    table(["", "fine-tuned", "7.3 phonons", "MACE-MP-0 medium", "experiment", "fine-tuned vs exp."], rows)
    print("  clamped-ion C44:", grab(section(t, "INTERNAL RELAXATION"), r"^C44 = C55 = C66\s*\|\s*([\d.]+)"),
          "GPa; fit R^2 (shear):", fr.get("yz", (None, None))[1], "; internal shift at 1 %:", fr.get("yz", [0] * 4)[3], "Ang")
    print("  reference state:", grab(t, r"^Lattice a, b, c\s*:.*-> ([\d.]+)", str), "Ang (fine-tuned),",
          grab(t, r"^--- MACE-MP-0 medium\n.*\nLattice a, b, c\s*:.*-> ([\d.]+)", str), "Ang (foundation)")


def fmt(v, spec="{:.2f}"):
    return "-" if v is None else spec.format(v)


def section_sound():
    print("\n=== 5.3 Sound velocities (m/s): Christoffel (this tool), 7.3 phonon slopes, experiment ===")
    t = report("base")
    if t is None:
        return note("base missing")
    sec = section(t, "SOUND VELOCITIES")
    c11, c12, c44 = SI_EXP["C11"] * 1e9, SI_EXP["C12"] * 1e9, SI_EXP["C44"] * 1e9
    rho = SI_RHO_EXP
    exp = {"[100]": (np.sqrt(c11 / rho), np.sqrt(c44 / rho), np.sqrt(c44 / rho)),
           "[110]": (np.sqrt((c11 + c12 + 2 * c44) / 2 / rho), np.sqrt(c44 / rho), np.sqrt((c11 - c12) / 2 / rho)),
           "[111]": (np.sqrt((c11 + 2 * c12 + 4 * c44) / 3 / rho),) + (np.sqrt((c11 - c12 + c44) / 3 / rho),) * 2}
    rows = []
    for d in ("[100]", "[110]", "[111]"):
        m = re.search(rf"^{re.escape(d)}\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)", sec, re.M)
        v = tuple(int(x) for x in m.groups())
        p = SI_PHONON_V[d]
        rows.append([d, " / ".join(map(str, v)), " / ".join(map(str, p)), " / ".join(f"{x:.0f}" for x in exp[d]),
                     f"{100 * max(abs(a / b - 1) for a, b in zip(v, p)):.1f} %"])
    table(["direction", "Christoffel L / T1 / T2", "7.3 phonons", "experiment", "largest diff. vs 7.3"], rows)
    print("  Debye temperature (Hill, isotropic):", grab(sec, r"^Debye temperature\s*:\s*([\d.]+)"), "K")


def section_fmax():
    print("\n=== 5.4 Force threshold of the internal relaxation ===")
    rows, data = [], []
    for p in sorted(glob.glob(os.path.join(REF, "fmax_*")), key=lambda s: -float(s.split("_")[-1])):
        f = float(p.split("_")[-1])
        t = report(os.path.basename(p))
        r = fit_rows(t, "fine-tuned")
        if "yz" not in r:
            continue
        data.append((f, r["yz"][0], r["yz"][1], r["yz"][4], r["xx"][0]))
        rows.append([f"{f:g}", f"{r['yz'][0]:.2f}", f"{r['yz'][1]:.6f}", f"{int(r['yz'][4])}", f"{r['xx'][0]:.2f}"])
    if not rows:
        return note("fmax_* missing")
    table(["--fmax (eV/Ang)", "C44 (GPa)", "R^2 (yz)", "relax steps (yz)", "C11 (GPa)"], rows)
    if PLOTS and data:
        import matplotlib.pyplot as plt
        d = np.array(data)
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
        axes[0].semilogx(d[:, 0], d[:, 1], "o-", color="#2a6fdb", label="C44 (relaxed ions)")
        axes[0].axhline(d[-1, 1], color="grey", ls=":", lw=1)
        axes[0].invert_xaxis()
        axes[0].set_xlabel("--fmax (eV/Å)")
        axes[0].set_ylabel("C44 (GPa)")
        axes[0].set_title("C44 converges once the atoms are relaxed", fontsize=10)
        axes[1].semilogx(d[:, 0], 1 - d[:, 2], "s-", color="#d9480f")
        axes[1].set_yscale("log")
        axes[1].invert_xaxis()
        axes[1].set_xlabel("--fmax (eV/Å)")
        axes[1].set_ylabel("1 − R² of the yz fit")
        axes[1].set_title("and the stress-strain points line up", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fmax.png"), dpi=150)
        plt.close(fig)


def section_amplitude():
    print("\n=== 5.5 Strain amplitude ===")
    t = report("conv")
    if t is None:
        return note("conv missing")
    for line in section(t, "STRAIN-AMPLITUDE").splitlines()[2:12]:
        print("  " + line)


def section_foundation():
    print("\n=== 5.6 Foundation models (GPa) ===")
    rows = []
    for m in ("small", "medium", "large"):
        t = report(f"foundation_{m}")
        if t is None:
            note(f"foundation_{m} missing")
            continue
        C = tensor(t)
        h = hill(t)
        a = grab(t, r"^Lattice a, b, c\s*:.*-> ([\d.]+)")
        rows.append([f"MACE-MP-0 {m}", f"{a:.4f}", f"{C[0, 0]:.1f}", f"{C[0, 1]:.1f}", f"{C[3, 3]:.1f}", f"{h[0]:.1f}",
                     f"{h[1]:.1f}", f"{h[3]:.3f}"])
    t = report("base")
    if t:
        C, h = tensor(t), hill(t)
        rows.append(["fine-tuned", f"{grab(t, r'^Lattice a, b, c\s*:.*-> ([\d.]+)'):.4f}", f"{C[0, 0]:.1f}",
                     f"{C[0, 1]:.1f}", f"{C[3, 3]:.1f}", f"{h[0]:.1f}", f"{h[1]:.1f}", f"{h[3]:.3f}"])
    eB, eG, _e, enu = vrh_cubic(SI_EXP["C11"], SI_EXP["C12"], SI_EXP["C44"])
    rows.append(["experiment", "5.4310", f"{SI_EXP['C11']}", f"{SI_EXP['C12']}", f"{SI_EXP['C44']}", f"{eB:.1f}",
                 f"{eG:.1f}", f"{enu:.3f}"])
    table(["model", "a (Ang)", "C11", "C12", "C44", "B", "G", "ν"], rows)
    if PLOTS and rows:
        import matplotlib.pyplot as plt
        names = ["C11", "C12", "C44", "B", "G"]
        exp = np.array([float(rows[-1][i]) for i in (2, 3, 4, 5, 6)])
        fig, ax = plt.subplots(figsize=(8.5, 4))
        cols = {"MACE-MP-0 small": "#f59f00", "MACE-MP-0 medium": "#d9480f", "MACE-MP-0 large": "#862e9c",
                "fine-tuned": "#2a6fdb"}
        series = [r for r in rows[:-1]]
        x = np.arange(len(names))
        for k, r in enumerate(series):
            v = np.array([float(r[i]) for i in (2, 3, 4, 5, 6)])
            ax.bar(x + (k - (len(series) - 1) / 2) * 0.2, 100 * (v / exp - 1), width=0.2, color=cols.get(r[0], "grey"),
                   label=r[0])
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(names)
        ax.set_ylabel("vs experiment (%)")
        ax.set_title("Silicon: the foundation models are soft, C44 most", fontsize=10)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "foundation_vs_experiment.png"), dpi=150)
        plt.close(fig)


def section_size_rattled():
    print("\n=== 5.7 Cell size, 5.8 a structure not at a minimum ===")
    rows = []
    for name, lab in (("base", "8 atoms (7 strains)"), ("si216", "216 atoms (5 strains)"),
                      ("rattled_noprerelax", "rattled, --no-pre-relax"), ("rattled", "rattled, pre-relaxed")):
        t = report(name)
        if t is None:
            note(f"{name} missing")
            continue
        C = tensor(t)
        sym = grab(t, r"^Crystal symmetry\s*:\s*(\S+)", str)
        dev = grab(t, r"^Projection\s*:.*\(([\d.]+) % of the largest", float)
        nind = grab(t, r"^Independent\s*:\s*(\d+)", int)
        f0 = grab(t, r"^Pre-relax\s*:.*max \|F\| ([\d.]+)")
        st = grab(t, r"largest component after ([\d.]+) GPa")
        p1 = sym == "P1"
        rows.append([lab, f"{C[0, 0]:.2f}", f"{C[0, 1]:.2f}", f"{C[3, 3]:.2f}", sym,
                     "-" if p1 or dev is None else f"{dev:.2f} %", "21 (no symmetry)" if p1 else str(nind), fmt(st),
                     grab(t, r"^Wall time\s*:\s*([\d.]+)", str) + " s"])
    table(["case", "C11", "C12", "C44", "symmetry", "deviation", "indep.", "residual stress GPa", "time"], rows)


def section_mg():
    print("\n=== 5.9 hcp magnesium (MACE-MP-0 medium) ===")
    t = report("mg_hcp")
    if t is None:
        return note("mg_hcp missing")
    C = tensor(t)
    print("  symmetry:", grab(t, r"^Crystal symmetry\s*:\s*(.*)$", str))
    for n, v in independent(t):
        print(f"  {n:<26} {v:8.2f}")
    rows = []
    for k, (m, n) in (("C11", (0, 0)), ("C12", (0, 1)), ("C13", (0, 2)), ("C33", (2, 2)), ("C44", (3, 3))):
        rows.append([k, f"{C[m, n]:.1f}", f"{MG_EXP[k]:.1f}", f"{100 * (C[m, n] / MG_EXP[k] - 1):+.0f} %"])
    rows.append(["C66 = (C11 - C12)/2", f"{C[5, 5]:.1f} (check {(C[0, 0] - C[0, 1]) / 2:.1f})",
                 f"{(MG_EXP['C11'] - MG_EXP['C12']) / 2:.1f}", "-"])
    table(["constant", "MACE-MP-0 medium", "experiment", "difference"], rows)
    print("  lattice:", grab(t, r"^Lattice a, b, c\s*:\s*(.*)$", str))
    print("  Born rows:")
    for line in section(t, "MECHANICAL STABILITY").splitlines()[4:9]:
        print("   " + line)


def section_graphene():
    print("\n=== 5.10 graphene (MACE-MP-0 medium) ===")
    t = report("graphene")
    if t is None:
        return note("graphene missing")
    for line in section(t, "IN-PLANE MODULI").splitlines()[2:12]:
        print("  " + line)
    print("  independent:", ", ".join(f"{n} {v:.1f}" for n, v in independent(t)))
    print("  lattice:", grab(t, r"^Lattice a, b, c\s*:\s*(.*)$", str))


def main():
    os.makedirs(FIG, exist_ok=True)
    section_base()
    section_sound()
    section_fmax()
    section_amplitude()
    section_foundation()
    section_size_rattled()
    section_mg()
    section_graphene()
    if PLOTS:
        print("\nfigures/: fmax.png, foundation_vs_experiment.png (+ the tool's own figures copied by make_reference.sh)")


if __name__ == "__main__":
    main()
