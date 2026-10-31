"""Figure of section 7 of the 5.11 README (figures/results.png), from runs/ and reference_results/ only:
(a) the Si spectrum, (b) how the Si activity moves with the displacement and the Optical.Mesh, (c) why graphene's
eps1(0) is not defined: its eps2(E) does not vanish at E -> 0, so the Kramers-Kronig integrand eps2/E diverges.

    python3 results_figures.py          (numpy + matplotlib; run from this folder)
"""
import gzip
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference_results")
SI_EXP = 520.0          # cm^-1, Si optical phonon at room temperature, about [README ref. 4]


def epsimg(path):
    with gzip.open(path, "rt") as f:
        rows = [line.split() for line in f if line.strip() and not line.startswith("#")]
    a = np.array(rows, dtype=float)
    return a[:, 0], a[:, 1]


def activity(report):
    for line in open(report):
        if re.match(r"^  1 ", line):
            return float(line.split()[9])
    return None


fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

# (a) Si spectrum
d = np.loadtxt(os.path.join(REF, "si", "raman_spectrum.dat"))
ax[0].plot(d[:, 0], d[:, 1] / d[:, 1].max(), color="#1f6f8b", lw=1.8, label="SIESTA PBE (T$_{2g}$, 512 cm$^{-1}$)")
ax[0].axvline(SI_EXP, color="#c0392b", ls="--", lw=1.2, label="experiment (≈ 520 cm$^{-1}$)")
ax[0].set_xlabel("Raman shift (cm$^{-1}$)")
ax[0].set_ylabel("intensity (normalized)")
ax[0].set_title("(a) bulk Si: one line")
ax[0].legend(frameon=False, fontsize=9, loc="upper right")

# (b) convergence of the activity
cases = [("si_d001", "δ = 0.01 Å"), ("si", "δ = 0.02 Å"), ("si_d004", "δ = 0.04 Å"), ("si_mesh16", "Optical.Mesh 16³")]
vals = [activity(os.path.join(REF, c, "stage3.txt")) for c, _ in cases]
ref = vals[1]
colors = ["#1f6f8b", "#1f6f8b", "#1f6f8b", "#e67e22"]
bars = ax[1].bar(range(len(cases)), [v / ref for v in vals], color=colors)
for b, v in zip(bars, vals):
    ax[1].text(b.get_x() + b.get_width() / 2, v / ref + 0.004, f"{v:.0f}", ha="center", fontsize=9)
ax[1].set_xticks(range(len(cases)), [label for _, label in cases], fontsize=9)
ax[1].set_ylim(0.95, 1.08)
ax[1].axhline(1.0, color="0.5", lw=0.8)
ax[1].set_ylabel("activity / activity(δ = 0.02 Å, mesh 32³)")
ax[1].set_title("(b) Si: displacement and optical mesh")

# (c) eps2 near E = 0: Si vs graphene
for path, name, color in (
        ("runs/si/raman_study/optical_disp/mode_01_plus_x/siesta.EPSIMG.gz", "Si (gap)", "#1f6f8b"),
        ("runs/graphene/raman_study/optical_disp/mode_02_plus_x/siesta.EPSIMG.gz", "graphene, +δ", "#c0392b"),
        ("runs/graphene/raman_study/optical_disp/mode_02_minus_x/siesta.EPSIMG.gz", "graphene, −δ", "#e67e22")):
    e, eps2 = epsimg(os.path.join(HERE, path))
    m = e <= 3.0
    ax[2].plot(e[m], eps2[m], color=color, lw=1.5, label=name)
ax[2].set_xlabel("E (eV)")
ax[2].set_ylabel("ε$_2$(E) along x")
ax[2].set_yscale("log")
ax[2].set_title("(c) ε$_2$ at low energy: ε$_1$(0) needs ε$_2$/E")
ax[2].legend(frameon=False, fontsize=9)

fig.tight_layout()
out = os.path.join(HERE, "figures", "results.png")
fig.savefig(out, dpi=130)
print("wrote", out)
for (c, label), v in zip(cases, vals):
    print(f"  {label:18s} activity {v:9.1f}  ({100 * (v / ref - 1):+.1f} %)")
