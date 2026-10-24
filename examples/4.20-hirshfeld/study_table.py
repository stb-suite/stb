#!/usr/bin/env python3
"""Collects the reference runs of the 4.20 example into the tables of README section 5, and draws its figures.

    python3 study_table.py [--root reference_results] [--png]

[A] the four systems: mean charge per element by every partition (Voronoi and Hirshfeld from SIESTA, simple
    Hirshfeld, Hirshfeld-I, Bader), the contraction ratio and the estimated remaining error
[B] the dipole of the molecules: point charges of every partition vs SIESTA's dipole of the density
[C] convergence of the iteration (H2O and NaCl): the step, the error estimated by Eq. 9 and the true error
    (against the --tol 1e-6 run in analysis_tight/)
[D] H2O, the mesh;  [E] H2O, the isolated-atom box;  [F] H2O, the basis
With --png: figures/h2o_studies.png (convergence, mesh/box, basis) and figures/stockholder_weights.png (the
Hirshfeld weight of H along the O-H bond, simple Hirshfeld vs Hirshfeld-I; needs the .RHO.xz of reference_results/h2o).
"""
import argparse
import glob
import lzma
import os
import re
import shutil
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SYSTEMS = ["h2o", "n2o", "nacl", "mgo"]
LABEL = {"h2o": "H2O", "n2o": "N2O", "nacl": "NaCl", "mgo": "MgO"}
OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]


def strip(s):
    return re.sub(r"\x1b\[[0-9;]*m", "", s)


def read_dat(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        rows.append({"idx": int(p[0]), "el": p[1], "zval": float(p[5]), "voronoi": float(p[6]),
                     "hsiesta": float(p[7]), "h0": float(p[8]), "hi": float(p[9]), "frac": float(p[11]),
                     "ion": p[12]})
    return rows


def read_iterations(path):
    head, data = None, []
    for line in open(path):
        if line.startswith("# iter"):
            head = line[1:].split()
        elif not line.startswith("#") and line.strip():
            data.append([float(x) for x in line.split()])
    return head, np.array(data)


def read_report(path):
    txt = strip(open(path).read())
    out = {}
    m = re.search(r"Contraction ratio\s+: ([0-9.]+)", txt)
    out["ratio"] = float(m.group(1)) if m else None
    m = re.search(r"Est\. remaining error\s+: ([0-9.eE+-]+)", txt)
    out["remaining"] = float(m.group(1)) if m else None
    m = re.search(r"Converged after (\d+) iteration", txt)
    out["iterations"] = int(m.group(1)) if m else None
    m = re.search(r"^\|mu\|\s*\|(.*)$", txt, re.M)
    if m:
        vals = [v.strip() for v in m.group(1).split("|")]
        out["mu"] = dict(zip(["voronoi", "hsiesta", "h0", "hi", "siesta"], [float(v) if v not in ("", "-", "N/A") else None for v in vals]))
    m = re.search(r"Density grid\s+: (\d+ x \d+ x \d+) \(spacing ([0-9.]+) Ang\)", txt)
    out["grid"] = m.group(1) if m else None
    out["spacing"] = float(m.group(2)) if m else None
    out["refs"] = {}
    for m in re.finditer(r"^(\w+)\s+\| (neutral|cation|anion)\s+\| (\d+)\s+\| ([0-9.]+|0 \(exact\))\s+\| ([0-9.]+|0 \(exact\))\s+\|"
                         r" [^|]*\| ([0-9.]+|-)\s+\|", txt, re.M):
        out["refs"][(m.group(1), m.group(2))] = {"r999": None if m.group(6) == "-" else float(m.group(6))}
    out["capped"] = "(capped)" in txt or "capped at the +-1 ion" in txt
    return out


def read_bader(path):
    if not os.path.isfile(path):
        return None
    q = {}
    for line in strip(open(path).read()).splitlines():
        p = [x.strip() for x in line.split("|")]
        if len(p) >= 5 and p[0].isdigit():
            try:
                q[int(p[0])] = float(p[4])
            except ValueError:
                pass
    return q or None


def mesh_used(study):
    for line in open(os.path.join(study, "combined", "calc.out")):
        if "Mesh cutoff (required, used)" in line:
            return float(line.split()[-2])
    return None


def species_mean(rows, key, q=None):
    out = {}
    for el in dict.fromkeys(r["el"] for r in rows):
        vals = [(q[r["idx"]] if q is not None else r[key]) for r in rows if r["el"] == el]
        out[el] = float(np.mean(vals))
    return out


def table(title, head, rows):
    print(f"\n{title}")
    w = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(head)]
    print(" | ".join(str(h).ljust(w[i]) for i, h in enumerate(head)))
    print("-+-".join("-" * x for x in w))
    for r in rows:
        print(" | ".join(str(c).ljust(w[i]) for i, c in enumerate(r)))


def f(x, fmt="{:+.3f}"):
    return "-" if x is None else fmt.format(x)


def convergence(root, name):
    """(n, step, estimated error, true error) for the slowest atom of a default-tol run, against the tight run."""
    head, it = read_iterations(os.path.join(root, name, "analysis_tight", "hirshfeld_iterations.dat"))
    q_conv = it[-1, 1:-1]
    rows = []
    for n in range(1, len(it)):
        step = it[n, -1]
        true = float(np.max(np.abs(it[n, 1:-1] - q_conv)))
        est = None
        if n >= 3:
            steps = it[n - 2:n + 1, -1]
            r = float(np.median(steps[1:] / steps[:-1]))
            est = step * r / (1 - r) if 0 < r < 1 else None
        rows.append((n, step, est, true))
    return rows, q_conv, head


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.path.join(HERE, "reference_results"))
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args()
    R = a.root

    # [A] the four systems
    rows, dip = [], []
    for s in SYSTEMS:
        d = os.path.join(R, s)
        if not os.path.isdir(d):
            continue
        dat = read_dat(os.path.join(d, "analysis", "hirshfeld_charges.dat"))
        rep = read_report(os.path.join(d, "analysis", "hirshfeld_analysis_report.txt"))
        tight = os.path.join(d, "analysis_tight", "hirshfeld_charges.dat")
        hi_conv = species_mean(read_dat(tight), "hi") if os.path.isfile(tight) else {}
        bader = read_bader(os.path.join(d, "bader", "stb_bader_report.txt"))
        means = {k: species_mean(dat, k) for k in ("voronoi", "hsiesta", "h0", "hi")}
        mb = species_mean(dat, None, bader) if bader else {}
        for el in means["hi"]:
            rows.append((LABEL[s], el, f(means["voronoi"][el]), f(means["hsiesta"][el]), f(means["h0"][el]),
                         f(means["hi"][el]), f(hi_conv.get(el)), f(mb.get(el)),
                         f(rep["ratio"], "{:.3f}"), f(rep["remaining"], "{:.4f}")))
        if "mu" in rep and rep["mu"].get("siesta") is not None:
            m = rep["mu"]
            dip.append((LABEL[s], f(m["voronoi"], "{:.3f}"), f(m["hsiesta"], "{:.3f}"), f(m["h0"], "{:.3f}"),
                        f(m["hi"], "{:.3f}"), f(m["siesta"], "{:.3f}")))
    table("[A] Mean charge per element (e), by partition; HI at the default --tol 0.005 and converged (--tol 1e-6)",
          ["System", "El", "Voronoi", "Hirsh. SIESTA", "Hirsh. (pass 0)", "HI", "HI conv.", "Bader", "ratio r",
           "est. error"], rows)
    if dip:
        table("[B] Dipole of the molecules, |mu| (D): point charges of each partition vs SIESTA's dipole of the density",
              ["System", "Voronoi", "Hirsh. SIESTA", "Hirsh. (pass 0)", "HI", "SIESTA"], dip)

    # [C] convergence
    conv = {}
    for s in SYSTEMS:
        if os.path.isfile(os.path.join(R, s, "analysis_tight", "hirshfeld_iterations.dat")):
            conv[s] = convergence(R, s)
            crow = [(n, f"{st:.6f}", f(es, "{:.6f}"), f"{tr:.6f}") for n, st, es, tr in conv[s][0] if n <= 14]
            table(f"[C] {LABEL[s]}: step, estimated remaining error (Eq. 9) and true error vs the converged charges (e)",
                  ["round", "step", "estimate", "true error"], crow)

    # [D], [E], [F] H2O studies
    def h2o_row(name, tag):
        d = os.path.join(R, name)
        if not os.path.isfile(os.path.join(d, "analysis", "hirshfeld_charges.dat")):
            return None
        dat = read_dat(os.path.join(d, "analysis", "hirshfeld_charges.dat"))
        rep = read_report(os.path.join(d, "analysis", "hirshfeld_analysis_report.txt"))
        o = [r for r in dat if r["el"] == "O"][0]
        bader = read_bader(os.path.join(d, "bader", "stb_bader_report.txt"))
        return {"tag": tag, "mesh": mesh_used(d), "spacing": rep["spacing"], "voronoi": o["voronoi"],
                "hsiesta": o["hsiesta"], "h0": o["h0"], "hi": o["hi"], "ratio": rep["ratio"],
                "rem": rep["remaining"], "mu_hi": rep.get("mu", {}).get("hi"), "mu": rep.get("mu", {}).get("siesta"),
                "r999_anion": rep["refs"].get(("O", "anion"), {}).get("r999"),
                "bader": bader.get(o["idx"]) if bader else None}

    mesh = [h2o_row(f"h2o_mesh{m}", f"{m}") for m in (200, 300)] + [h2o_row("h2o", "400")] + [h2o_row("h2o_mesh600", "600")]
    mesh = [x for x in mesh if x]
    table("[D] H2O, MeshCutoff (Ry): q(O) in e",
          ["requested", "used", "grid spacing (Ang)", "Voronoi", "Hirsh. (pass 0)", "HI", "ratio r", "|mu| HI (D)"],
          [(x["tag"], f(x["mesh"], "{:.1f}"), f(x["spacing"], "{:.4f}"), f(x["voronoi"]), f(x["h0"], "{:+.4f}"),
            f(x["hi"], "{:+.4f}"), f(x["ratio"], "{:.3f}"), f(x["mu_hi"], "{:.3f}")) for x in mesh])
    box = [h2o_row(f"h2o_box{v}", f"{v}") for v in (10, 14)] + [h2o_row("h2o", "20")] + [h2o_row("h2o_box28", "28")]
    box = [x for x in box if x]
    table("[E] H2O, isolated-atom box (Ang): q(O) in e",
          ["box", "Hirsh. (pass 0)", "HI", "ratio r"],
          [(x["tag"], f(x["h0"], "{:+.4f}"), f(x["hi"], "{:+.4f}"), f(x["ratio"], "{:.3f}")) for x in box])
    basis = [h2o_row(n, t) for n, t in (("h2o_szp_s001", "SZP 0.01"), ("h2o_dzp_s002", "DZP 0.02"), ("h2o", "DZP 0.01"),
                                         ("h2o_dzp_s0005", "DZP 0.005"), ("h2o_tzp_s001", "TZP 0.01"),
                                         ("h2o_tzp_s0005", "TZP 0.005"))]
    basis = [x for x in basis if x]
    table("[F] H2O, basis (PAO.BasisSize, PAO.EnergyShift in Ry): q(O) in e, dipoles in D",
          ["basis", "Voronoi", "Hirsh. SIESTA", "Hirsh. (pass 0)", "HI", "Bader", "r99.9% O- (Ang)", "|mu| HI", "|mu| SIESTA"],
          [(x["tag"], f(x["voronoi"]), f(x["hsiesta"]), f(x["h0"]), f(x["hi"]), f(x["bader"]),
            f(x["r999_anion"], "{:.2f}"), f(x["mu_hi"], "{:.2f}"), f(x["mu"], "{:.2f}")) for x in basis])

    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
        for k, s in enumerate(c for c in ("h2o", "n2o", "nacl", "mgo") if c in conv):
            rr = conv[s][0][:20]
            ax[0].semilogy([x[0] for x in rr], [max(x[3], 1e-7) for x in rr], "-o", ms=3.5, color=OKABE[k],
                           label=f"{LABEL[s]}")
            ax[0].semilogy([x[0] for x in rr if x[2]], [x[2] for x in rr if x[2]], "--", color=OKABE[k], alpha=0.9)
        ax[0].plot([], [], "k-o", ms=3.5, label="true error (vs --tol 1e-6)")
        ax[0].plot([], [], "k--", label="estimate, Eq. 9")
        ax[0].axhline(0.005, color="0.4", lw=1, ls=":")
        ax[0].text(19.5, 0.0058, "default --tol", color="0.3", fontsize=8, ha="right")
        ax[0].set_ylim(1e-6, 1)
        ax[0].set_xlabel("round")
        ax[0].set_ylabel("max |q - q_converged| (e)")
        ax[0].set_title("How far from the fixed point")
        ax[0].legend(fontsize=8, ncol=2)
        ax[0].grid(alpha=0.3)
        base = [x for x in mesh if x["tag"] == "400"]
        if base:
            q0 = base[0]["hi"]
            ax[1].plot([x["mesh"] for x in mesh], [1e3 * (x["hi"] - q0) for x in mesh], "-o", color=OKABE[0],
                       label="MeshCutoff used (Ry), bottom axis")
            if box:
                t = ax[1].twiny()
                t.plot([float(x["tag"]) for x in box], [1e3 * (x["hi"] - q0) for x in box], "--^", color=OKABE[1],
                       label="isolated-atom box (Ang), top axis")
                t.set_xlabel("isolated-atom box (Ang)", color=OKABE[1])
                t.legend(loc="lower right", fontsize=8)
            ax[1].set_ylim(-1.0, 1.0)
            ax[1].set_xlabel("MeshCutoff used (Ry)")
            ax[1].set_ylabel("q(O) - q(O) at 400 Ry, 20 Ang  (10^-3 e)")
            ax[1].set_title("H2O, Hirshfeld-I: mesh and box", pad=4, y=1.13)
            ax[1].legend(loc="upper right", fontsize=8)
            ax[1].grid(alpha=0.3)
        if basis:
            xs = np.arange(len(basis))
            for j, (key, lab, col) in enumerate((("voronoi", "Voronoi (SIESTA)", "0.75"), ("h0", "simple Hirshfeld", "0.45"),
                                                 ("hi", "Hirshfeld-I", "0.1"))):
                ax[2].bar(xs + (j - 1) * 0.26, [x[key] for x in basis], 0.24, color=col, label=lab)
            t2 = ax[2].twinx()
            t2.plot(xs, [x["r999_anion"] for x in basis], "-D", color=OKABE[1], ms=5, label="r(99.9 %) of O- (right)")
            t2.set_ylabel("r(99.9 %) of the O- reference (Ang)", color=OKABE[1])
            t2.set_ylim(1.6, 2.3)
            t2.legend(loc="lower right", fontsize=8)
            ax[2].set_xticks(xs, [x["tag"] for x in basis], rotation=20, fontsize=8)
            ax[2].set_ylim(-0.9, 0)
            ax[2].set_ylabel("q(O) (e)")
            ax[2].set_title("H2O: basis (size, PAO.EnergyShift in Ry)")
            ax[2].legend(loc="lower left", fontsize=8)
            ax[2].grid(alpha=0.3, axis="y")
        fig.tight_layout()
        out = os.path.join(HERE, "figures", "h2o_studies.png")
        fig.savefig(out, dpi=130)
        print(f"\nwrote {out}")
        weights_figure(R, os.path.join(HERE, "figures", "stockholder_weights.png"))


def weights_figure(R, out):
    """w_H along the O-H bond of water, simple Hirshfeld vs Hirshfeld-I, from the real reference atoms."""
    import sisl
    import matplotlib.pyplot as plt
    from stb.core import hirshfeld as hf, structure_io
    from stb.core.rho_io import read_total

    d = os.path.join(R, "h2o")
    tmp = tempfile.mkdtemp()
    try:
        def profile(sub):
            src = os.path.join(d, sub, "siesta.RHO.xz")
            dst = os.path.join(tmp, sub.replace("/", "_") + ".RHO")
            with lzma.open(src) as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            g = read_total(sisl.get_sile(dst))
            return hf.radial_density_profile(np.asarray(g.grid), g.lattice.cell, (0.5, 0.5, 0.5))

        prof = {k: profile(k) for k in ("neutral/O", "neutral/H", "ions/O/anion")}
        dat = read_dat(os.path.join(d, "analysis", "hirshfeld_charges.dat"))
        qO = [r for r in dat if r["el"] == "O"][0]["hi"]
        qH = [r for r in dat if r["el"] == "H"][0]["hi"]
        st = structure_io.read_fdf(os.path.join(d, "combined", "structure.fdf"))
        pm = structure_io.to_pymatgen(st)
        dOH = pm.get_distance(0, 1)
        x = np.linspace(-1.0, dOH + 1.0, 600)        # along O -> H, O at 0
        rO, rH = np.abs(x), np.abs(x - dOH)

        def ev(p, r):
            return np.interp(r, p[0], p[1], left=0, right=0)
        rhoO0, rhoH0, rhoOa = ev(prof["neutral/O"], rO), ev(prof["neutral/H"], rH), ev(prof["ions/O/anion"], rO)
        fO = min(abs(qO), 1)
        rhoO_hi = (1 - fO) * rhoO0 + fO * rhoOa
        rhoH_hi = (1 - min(abs(qH), 1)) * rhoH0           # H+ has zero density
        wH0 = rhoH0 / np.maximum(rhoO0 + rhoH0, 1e-30)
        wHi = rhoH_hi / np.maximum(rhoO_hi + rhoH_hi, 1e-30)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].semilogy(x, rhoO0, color=OKABE[0], label="O neutral")
        ax[0].semilogy(x, rhoO_hi, "--", color=OKABE[0], label=f"O pro-atom, q = {qO:+.2f}")
        ax[0].semilogy(x, rhoH0, color=OKABE[1], label="H neutral")
        ax[0].semilogy(x, np.maximum(rhoH_hi, 1e-6), "--", color=OKABE[1], label=f"H pro-atom, q = {qH:+.2f}")
        ax[0].set_ylim(1e-3, 3e2)
        ax[0].set_xlabel("position along O -> H (Ang)")
        ax[0].set_ylabel("density (e/Ang^3)")
        ax[0].set_title("Pro-atoms along the O-H bond")
        ax[0].legend(fontsize=8)
        ax[0].grid(alpha=0.3)
        ax[1].plot(x, wH0, color="0.45", label="w_H, simple Hirshfeld (neutral pro-atoms)")
        ax[1].plot(x, wHi, color="0.05", label="w_H, Hirshfeld-I (converged pro-atoms)")
        for k, (w, c, lab) in enumerate(((wH0, "0.45", "simple"), (wHi, "0.05", "Hirshfeld-I"))):
            i = np.nonzero((x > 0) & (x < dOH) & (w >= 0.5))[0]
            if i.size:
                ax[1].axvline(x[i[0]], color=c, ls=":", lw=1)
                ax[1].annotate(f"{lab}: w_H = 1/2 at {x[i[0]]:.2f} Ang from O", (x[i[0]], 0.5),
                               xytext=(-0.95, 0.30 - 0.08 * k), color=c, fontsize=8,
                               arrowprops=dict(arrowstyle="->", color=c, lw=0.8))
        for xp, lab in ((0, "O"), (dOH, "H")):
            ax[1].axvline(xp, color="0.8", lw=1)
            ax[1].text(xp, 1.03, lab, ha="center")
        ax[1].set_ylim(0, 1.1)
        ax[1].set_xlabel("position along O -> H (Ang)")
        ax[1].set_ylabel("share of H")
        ax[1].set_title("Hydrogen's share of the density")
        ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 0.88))
        ax[1].grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(out, dpi=130)
        print(f"wrote {out}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
