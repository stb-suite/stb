#!/usr/bin/env python3
"""Tables of section 6 of the 7.5 README and its summary figures, from reference_results/ only.

    python3 results_figures.py         # numpy + ase + matplotlib; no MACE, no SIESTA

Reads the stb-mlneb reports and band files of every reference case, and the SIESTA single points of the real
loop (reference_results/loop_round*/siesta_check/, foundation_siesta/), and writes figures/loop.png,
figures/vacancy_hop.png and figures/images.png.
"""
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "reference_results")
FIG = os.path.join(HERE, "figures")


def path(*p):
    return os.path.join(RES, *p)


def report(case):
    f = path(case, "stb_mlneb_report.txt")
    return open(f).read() if os.path.isfile(f) else None


def barrier(text, which="fine-tuned"):
    """(spline, highest image) of one model's BARRIER block."""
    if text is None:
        return None
    blk = text.split("BARRIER", 1)[1]
    if f"--- {which}" in blk:
        blk = blk.split(f"--- {which}", 1)[1]
    m = re.search(r"forward\s+\(initial -> final\)\s+\|\s+([-\d.]+)\s+\|\s+([-\d.]+)", blk)
    return (float(m.group(1)), float(m.group(2))) if m else None


def neb_line(text, which=None):
    """steps, residual force and time of the NEB BAND block."""
    pat = r"--- %s: ([^;\n]*); (\d+) shaping \+ (\d+) climbing step\(s\).*?residual force ([\d.]+).*?; ([\d.]+) s" % (
        re.escape(which) if which else r"[^:]+")
    m = re.search(pat, text)
    if not m:
        return None
    return {"status": m.group(1), "steps": int(m.group(2)) + int(m.group(3)), "res": float(m.group(4)),
            "time": float(m.group(5))}


def member_barriers(text):
    m = re.search(r"Barrier per member: ([\d., ]+) eV", text or "")
    return [float(x) for x in m.group(1).split(",")] if m else []


def siesta_energy(folder):
    """SIESTA free energy (eV) of a finished single point."""
    out = os.path.join(folder, "siesta.out")
    if not os.path.isfile(out):
        return None
    e = None
    with open(out, errors="replace") as f:
        for line in f:
            if line.strip().startswith("siesta:") and "Free energy" in line and "=" in line:
                try:
                    e = float(line.split("=")[1].split()[0])
                except ValueError:
                    pass
            elif re.match(r"\s*siesta:\s+FreeEng\s*=", line):
                e = float(line.split("=")[1].split()[0])
    return e


def siesta_profile(loop):
    d = path(loop, "siesta_check")
    out = {}
    for name in sorted(os.listdir(d)):
        e = siesta_energy(os.path.join(d, name))
        if e is not None:
            out[int(name.split("_")[1])] = e
    e0 = out[0]
    return {i: e - e0 for i, e in out.items()}


def band(case, foundation=False):
    f = path(case, "neb_band_foundation.dat" if foundation else "neb_band.dat")
    return np.loadtxt(f) if os.path.isfile(f) else None


def spline(case, k):
    f = path(case, "plot", f"spline_{k}.dat")
    return np.loadtxt(f) if os.path.isfile(f) else None


def table(rows, head):
    w = [max(len(str(r[j])) for r in rows + [head]) for j in range(len(head))]
    print(" | ".join(str(h).ljust(x) for h, x in zip(head, w)))
    print("-+-".join("-" * x for x in w))
    for r in rows:
        print(" | ".join(str(c).ljust(x) for c, x in zip(r, w)))
    print()


def main():
    os.makedirs(FIG, exist_ok=True)
    r2, r3 = report("loop_round2"), report("loop_round3")
    s2, s3 = siesta_profile("loop_round2"), siesta_profile("loop_round3")
    fd = {k: siesta_energy(path("foundation_siesta", k)) for k in ("fd_00", "fd_03", "ft_00", "ft_03")}
    b2, b3 = band("loop_round2"), band("loop_round3")
    bf = band("loop_round3", foundation=True)

    # ---- the loop ------------------------------------------------------------------------------------
    print("== Section 6: the loop (real runs, reference_results/loop_round*/)\n")
    rows = []
    for lab, b, s in (("round 2", b2, s2), ("round 3", b3, s3)):
        for i in range(len(b)):
            rows.append([lab, i, f"{b[i, 1]:.3f}", f"{b[i, 2]:+.4f}",
                         f"{s[i]:+.4f}" if i in s else "-", f"{b[i, 2] - s[i]:+.4f}" if i in s else "-"])
    table(rows, ["model", "image", "s (Ang)", "model E-E0 (eV)", "SIESTA E-E0 (eV)", "model - SIESTA"])
    sf = fd["fd_03"] - fd["fd_00"]
    m2, m3 = member_barriers(r2), member_barriers(r3)
    summary = [["round 2", f"{b2[3, 2]:.4f}", f"{min(m2):.3f}-{max(m2):.3f}", f"{s2[3]:.4f}",
                f"{100 * (b2[3, 2] - s2[3]) / s2[3]:+.0f} %"],
               ["round 3", f"{b3[3, 2]:.4f}", f"{min(m3):.3f}-{max(m3):.3f}", f"{s3[3]:.4f}",
                f"{100 * (b3[3, 2] - s3[3]) / s3[3]:+.0f} %"],
               ["MACE-MP-0 medium", f"{bf[3, 2]:.4f}", "-", f"{sf:.4f}", f"{100 * (bf[3, 2] - sf) / sf:+.0f} %"]]
    table(summary, ["model", "barrier (eV)", "committee", "SIESTA on its path", "model vs SIESTA"])
    print(f"SIESTA, foundation endpoint vs round-2 endpoint: {fd['fd_00'] - fd['ft_00']:+.4f} eV; "
          f"foundation saddle vs round-2 saddle: {fd['fd_03'] - fd['ft_03']:+.4f} eV")
    print(f"lowest SIESTA saddle - lowest SIESTA endpoint (both paths): "
          f"{min(fd['fd_03'], fd['ft_03']) - min(fd['fd_00'], fd['ft_00']):.4f} eV\n")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    cols = {"r2": "#c92a2a", "r3": "#2a6fdb", "fd": "#e8590c"}
    for key, case, k, b, lab in (("r2", "round2", 0, b2, "round 2"), ("r3", "round3", 0, b3, "round 3"),
                                 ("fd", "round3", 1, bf, "MACE-MP-0 medium")):
        sp = spline(case, k)
        if sp is not None:
            ax.plot(sp[:, 0], sp[:, 1], "-", color=cols[key], lw=1.6)
        ax.plot(b[:, 1], b[:, 2], "o", color=cols[key], ms=6, label=f"{lab}: {b[:, 2].max():.3f} eV")
    for key, s, b in (("r2", s2, b2), ("r3", s3, b3)):
        ax.plot([b[i, 1] for i in s], [s[i] for i in s], "s", mfc="white", mec=cols[key], mew=2, ms=8)
    ax.plot([bf[0, 1], bf[3, 1]], [0, sf], "s", mfc="white", mec=cols["fd"], mew=2, ms=8)
    ax.plot([], [], "s", mfc="white", mec="k", mew=2, ms=8, label="SIESTA on each model's path")
    ax.set_xlabel("reaction coordinate s (Ang)")
    ax.set_ylabel("E - E(initial) (eV)")
    ax.set_title("(a) The band of each model, and SIESTA on it")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    i2 = sorted(s2)
    ax.bar(np.array(i2) - 0.2, [b2[i, 2] - s2[i] for i in i2], 0.4, color=cols["r2"], label="round 2")
    i3 = sorted(s3)
    ax.bar(np.array(i3) + 0.2, [b3[i, 2] - s3[i] for i in i3], 0.4, color=cols["r3"], label="round 3")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("image")
    ax.set_ylabel("model - SIESTA, E - E0 (eV)")
    ax.set_title("(b) The model's error image by image")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[2]
    labels = ["round 2", "round 3", "MACE-MP-0\nmedium"]
    model = [b2[3, 2], b3[3, 2], bf[3, 2]]
    sies = [s2[3], s3[3], sf]
    x = np.arange(3)
    ax.bar(x - 0.2, model, 0.4, color=[cols["r2"], cols["r3"], cols["fd"]], label="model")
    ax.bar(x + 0.2, sies, 0.4, color="#adb5bd", label="SIESTA on the model's path")
    for xi, m in ((0, m2), (1, m3)):
        ax.errorbar(xi - 0.2, np.mean(m), yerr=[[np.mean(m) - min(m)], [max(m) - np.mean(m)]], fmt="none",
                    ecolor="k", capsize=5, lw=1.5)
    ax.plot([], [], "|", color="k", ms=10, mew=1.5, label="committee range")
    for y, lab in ((0.36, "GGA, 215-999 atoms, extrapolated [15]"), (0.45, "experiment, neutral V [11]")):
        ax.axhline(y, color="#495057", ls=":", lw=1)
        ax.text(2.45, y + 0.01, lab, fontsize=7, ha="right", color="#495057")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("migration barrier (eV)")
    ax.set_title("(c) Barriers: model, committee, SIESTA")
    ax.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 0.88))
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "loop.png"), dpi=150)
    plt.close(fig)

    # ---- the hop -------------------------------------------------------------------------------------
    from ase.io import read
    frames = read(path("loop_round3", "neb_path.xyz"), ":")
    a0 = frames[0]
    d0 = a0.get_distances(0, range(1, len(a0)), mic=True)
    near = 1 + np.argsort(d0)[:3]                               # the hopping atom's three bonds at the start
    df_ = frames[-1].get_distances(0, range(1, len(a0)), mic=True)
    near_f = 1 + np.argsort(df_)[:3]                            # its three bonds at the end
    s = b3[:, 1]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    for j, col, lab in ((near, "#2a6fdb", "its 3 neighbours at the start"),
                        (near_f, "#c92a2a", "its 3 neighbours at the end")):
        for n_, idx in enumerate(j):
            d = [fr.get_distance(0, idx, mic=True) for fr in frames]
            ax.plot(s, d, "o-", color=col, ms=4, lw=1.2, label=lab if n_ == 0 else None)
    ax.axhline(2.3730, color="#868e96", ls=":", lw=1)
    ax.text(0.75, 2.38, "bulk bond 2.373 Ang", fontsize=8, color="#495057")
    ax.set_xlabel("reaction coordinate s (Ang)")
    ax.set_ylabel("distance from the hopping atom (Ang)")
    ax.set_title("(a) Bonds broken and formed along the band (round 3)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    # projection onto the plane spanned by the hop direction [111] and [1-10]
    p = np.array([fr.positions[0] for fr in frames])
    u = np.array([1, 1, 1]) / np.sqrt(3)
    v = np.array([1, -1, 0]) / np.sqrt(2)
    centre = p[3]
    pos0 = frames[0].positions
    for k, a in enumerate(pos0):
        dvec = a - centre
        dvec -= frames[0].cell.T @ np.round(np.linalg.solve(frames[0].cell.T, dvec))
        if np.linalg.norm(dvec) < 4.2 and k != 0:
            ax.plot(dvec @ u, dvec @ v, "o", color="#adb5bd", ms=11, mec="#495057", zorder=2)
    for lab, site in (("vacancy (start)", np.zeros(3)), ("vacancy (end)", p[0])):
        dvec = site - centre
        dvec -= frames[0].cell.T @ np.round(np.linalg.solve(frames[0].cell.T, dvec))
        ax.plot(dvec @ u, dvec @ v, "s", mfc="none", mec="#c92a2a", ms=13, mew=1.5, zorder=3)
        ax.annotate(lab, (dvec @ u, dvec @ v), xytext=(0, -20), textcoords="offset points", ha="center",
                    fontsize=8, color="#c92a2a")
    q = np.array([(x - centre) for x in p])
    ax.plot(q @ u, q @ v, "o-", color="#2a6fdb", ms=6, lw=1.5, zorder=4, label="the hopping atom, image by image")
    ax.plot(q[3] @ u, q[3] @ v, "*", color="#2a6fdb", ms=16, mec="k", zorder=5, label="saddle (split vacancy)")
    ax.set_aspect("equal")
    ax.set_xlabel("along [111] (Ang)")
    ax.set_ylabel("along [1-10] (Ang)")
    ax.set_title("(b) The hop, projected (atoms within 4.2 Ang of the saddle)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "vacancy_hop.png"), dpi=150)
    plt.close(fig)
    print("Saddle: distance of the hopping atom to its six nearest atoms (round 3): "
          + " ".join(f"{x:.3f}" for x in np.sort(frames[3].get_distances(0, range(1, len(a0)), mic=True))[:6]) + "\n")

    # ---- the reference cases ---------------------------------------------------------------------------
    print("== Reference cases (reference_results/<case>/, make_reference.sh)\n")
    rows = []
    cases = [("foundation_small", "MACE-MP-0 small", None), ("foundation_medium", "MACE-MP-0 medium", None),
             ("foundation_large", "MACE-MP-0 large", None), ("round2", "round 2", "fine-tuned"),
             ("round3", "round 3, 7 images", "fine-tuned"), ("images_5", "round 3, 5 images", "fine-tuned"),
             ("images_9", "round 3, 9 images", "fine-tuned"), ("images_13", "round 3, 13 images", "fine-tuned"),
             ("no_idpp", "round 3, --no-idpp", "fine-tuned"),
             ("fmax_0.01", "round 3, 7 images, --fmax 0.01", "fine-tuned"),
             ("fmax_0.01_images_9", "round 3, 9 images, --fmax 0.01", "fine-tuned"),
             ("no_prerelax", "round 3, --no-pre-relax-endpoints", "fine-tuned")]
    img, img01 = [], []
    for case, lab, which in cases:
        t = report(case)
        if t is None:
            rows.append([lab, "(not run)", "", "", "", ""])
            continue
        bb = barrier(t, which or "MACE")
        nl = neb_line(t, which)
        rows.append([lab, f"{bb[0]:.4f}", f"{bb[1]:.4f}", nl["steps"] if nl else "?",
                     nl["status"] if nl else "?", f"{nl['time']:.0f}" if nl else "?"])
        m = re.match(r"round 3, (\d+) images(, --fmax 0.01)?$", lab)
        if m:
            (img01 if m.group(2) else img).append((int(m.group(1)), bb[0], nl["res"] if nl else np.nan))
    table(rows, ["case", "barrier, spline (eV)", "highest image (eV)", "NEB steps", "NEB", "NEB time (s)"])
    if img:
        img.sort()
        img01.sort()
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        n_, b_, r_ = zip(*img)
        ax.plot(n_, b_, "o-", color="#2a6fdb", ms=7, label="--fmax 0.05 (default)")
        for n, b, r in img:
            ax.annotate(f"{r:.3f}", (n, b), xytext=(4, 6), textcoords="offset points", fontsize=7, color="#2a6fdb")
        if img01:
            n1, b1, r1 = zip(*img01)
            ax.plot(n1, b1, "s-", color="#c92a2a", ms=7, label="--fmax 0.01")
        if s3:
            ax.axhline(s3[3], color="#2b8a3e", ls=":", lw=1.2, label="SIESTA on the 7-image path")
        ax.set_xlabel("--n-images (endpoints included)")
        ax.set_ylabel("barrier (eV)")
        ax.set_ylim(0.24, 0.32)
        ax.set_title("Round 3: barrier against images and --fmax\n(labels: residual NEB force, eV/Ang)", fontsize=10)
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "images.png"), dpi=150)
        plt.close(fig)
    print("figures/loop.png, figures/vacancy_hop.png" + (", figures/images.png" if img else ""))


if __name__ == "__main__":
    main()
