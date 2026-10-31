"""Figures of the 2.19 example (stb-nativecharges), from the charges SIESTA printed in the logs of runs/ (read with the
same parsers stb-nativecharges uses):

  figures/basis_dependence.png   water: the charge of O by method and basis; the dipole of the point charges vs SIESTA's
  figures/methods.png            NaCl and water (DZP) side by side: the three methods

    python3 theory_figures.py        # seconds
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stb.core import siesta_log

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "runs")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
METHODS = [("Mulliken", siesta_log.get_mulliken_charges, "#c0392b"),
           ("Hirshfeld", siesta_log.get_hirshfeld_charges, "#1f6f8b"),
           ("Voronoi", siesta_log.get_voronoi_charges, "#2e7d32")]
BASES = ["SZ", "DZ", "DZP", "TZP"]
DEBYE = 1 / 0.20819434
dz_h = 0.585882                       # height of the H atoms above O along z (Ang), from the input geometry

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
x = np.arange(len(BASES))
for k, (name, reader, color) in enumerate(METHODS):
    qo = [reader(os.path.join(R, "h2o_basis", b, "calc.out"))[0]["charge"] for b in BASES]
    ax[0].plot(x, qo, "o-", color=color, label=name, lw=1.6, ms=7)
    qh = [reader(os.path.join(R, "h2o_basis", b, "calc.out"))[1]["charge"] for b in BASES]
    ax[1].plot(x, [2 * q * dz_h * DEBYE for q in qh], "o-", color=color, label=f"{name} point charges", lw=1.6, ms=7)
siesta = []
for b in BASES:
    d = None
    for line in open(os.path.join(R, "h2o_basis", b, "calc.out")):
        if "Electric dipole (Debye)" in line:
            d = float(line.split("=")[1].split()[2])
    siesta.append(d)
ax[1].plot(x, siesta, "s--", color="k", label="SIESTA (the density)", ms=7)
ax[0].axhline(0, color="0.5", lw=0.8)
for a in ax:
    a.set_xticks(x)
    a.set_xticklabels(BASES)
    a.set_xlabel("basis")
    a.grid(alpha=0.3)
    a.legend(frameon=False, fontsize=8)
ax[0].set_ylabel("charge of O (e)")
ax[0].set_title("(a) water: q(O) by method and basis")
ax[1].set_ylabel("dipole along the C2 axis (D)")
ax[1].set_title("(b) the point charges' dipole vs SIESTA's")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "basis_dependence.png"), dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.4, 3.8))
labels = ["Na in NaCl", "Cl in NaCl", "O in H$_2$O (DZP)", "H in H$_2$O (DZP)"]
picks = [("nacl", 0), ("nacl", 1), ("h2o_basis/DZP", 0), ("h2o_basis/DZP", 1)]
w = 0.26
for k, (name, reader, color) in enumerate(METHODS):
    vals = [reader(os.path.join(R, d, "calc.out"))[i]["charge"] for d, i in picks]
    ax.bar(np.arange(4) + (k - 1) * w, vals, width=w, color=color, label=name)
ax.axhline(0, color="0.4", lw=0.8)
ax.set_xticks(np.arange(4))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("charge (e)")
ax.set_title("One density, three partitions")
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "methods.png"), dpi=150)
plt.close(fig)
print("figures/basis_dependence.png, methods.png")
