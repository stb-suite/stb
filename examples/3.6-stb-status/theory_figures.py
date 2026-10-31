"""Figures of the 3.6 example (stb-status), from the logs of project/ read with stb-status's own parser:

  figures/relax_progress.png   the converged and the cut-short relaxation: max force, stress and energy per move
  figures/strain_energy.png    the three single points of project/strain/: E(a) and the parabola through them

    python3 theory_figures.py        # seconds
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stb.status import parse_log

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "project")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

full = parse_log(os.path.join(P, "relax", "calc.out"))
short = parse_log(os.path.join(P, "relax_short", "calc.out"))
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for info, color, name in ((full, "#1f6f8b", "relax/ (converged)"), (short, "#c0392b", "relax_short/ (3 moves)")):
    mv = [s["step"] for s in info["steps"]]
    ax[0].semilogy(mv, [s["max_force"] for s in info["steps"]], "o-", color=color, label=name)
    ax[1].semilogy(mv, [s["max_stress"] / 10 for s in info["steps"]], "o-", color=color, label=name)
    e0 = full["final"]["Total"]
    ax[2].plot(mv, [1000 * (s["E_KS"] - e0) / 2 for s in info["steps"]], "o-", color=color, label=name)
ax[0].axhline(full["force_tol"], color="k", ls="--", lw=0.8, label="MD.MaxForceTol")
ax[1].axhline(full["stress_tol"], color="k", ls="--", lw=0.8, label="MD.MaxStressTol")
ax[0].set_ylabel("max force (eV/Å)")
ax[1].set_ylabel("largest stress component (GPa)")
ax[2].set_ylabel("E_KS − E_final (meV/atom)")
for a, t in zip(ax, ("(a) forces", "(b) stress", "(c) energy")):
    a.set_xlabel("CG move")
    a.set_title(t)
    a.legend(frameon=False, fontsize=8)
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "relax_progress.png"), dpi=150)
plt.close(fig)

a_vals, e_vals = [], []
for d in sorted(os.listdir(os.path.join(P, "strain"))):
    info = parse_log(os.path.join(P, "strain", d, "calc.out"))
    a_vals.append(float(d.split("_")[1]))
    e_vals.append(info["final"]["Total"] / 2)
a_vals, e_vals = np.array(a_vals), np.array(e_vals)
c = np.polyfit(a_vals, e_vals, 2)
a0 = -c[1] / (2 * c[0])
x = np.linspace(a_vals.min() - 0.03, a_vals.max() + 0.03, 100)
fig, ax = plt.subplots(figsize=(5.4, 3.8))
ax.plot(x, 1000 * (np.polyval(c, x) - np.polyval(c, a0)), color="#1f6f8b", lw=1.4, label="parabola")
ax.plot(a_vals, 1000 * (e_vals - np.polyval(c, a0)), "o", color="#c0392b", ms=7, label="strain/a_*")
ax.axvline(a0, color="0.5", ls=":", lw=0.9, label=f"minimum a = {a0:.3f} Å")
ax.set_xlabel("lattice constant a (Å)")
ax.set_ylabel("E − E_min (meV/atom)")
ax.set_title("Si: the series of strain/")
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "strain_energy.png"), dpi=150)
plt.close(fig)
print(f"parabola minimum a = {a0:.4f} Ang")
print("figures/relax_progress.png, strain_energy.png")
