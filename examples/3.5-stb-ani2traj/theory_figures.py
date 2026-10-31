"""Figures of the 3.5 example (stb-ani2traj), drawn from the extended-XYZ trajectories the tool writes -- every number
plotted (time, temperature, energies, volume, cell) is a per-frame value stb-ani2traj embedded in the file:

  figures/thermo.png    the three runs: NVE energies, NVT temperature, NPT volume and cell

    python3 theory_figures.py        # seconds
"""
import os
import shutil
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")


def convert(work, run, label, *extra):
    out = f"{label}_{'_'.join(x.strip('-') for x in extra) or 'plain'}.xyz"
    subprocess.run(["stb-ani2traj", "--no-intro", "-l", label, "-of", "xyz", "-o", out, *extra],
                   cwd=os.path.join(work, run), check=True, stdout=subprocess.DEVNULL)
    return read(os.path.join(work, run, out), index=":")


def series(frames, key):
    return np.array([f.info.get(key, np.nan) for f in frames], dtype=float)


os.makedirs(FIG, exist_ok=True)
work = tempfile.mkdtemp()
try:
    for r in ("water_nvt", "si_nve", "si_npt"):
        shutil.copytree(os.path.join(HERE, "siesta_runs", r), os.path.join(work, r))
    nve = convert(work, "si_nve", "si")
    nvt = convert(work, "water_nvt", "water")
    npt = convert(work, "si_npt", "si")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.0))
    t = series(nve, "Time")
    n = len(nve[0])
    e_ks, e_tot = series(nve, "E_KS"), series(nve, "E_tot")
    ax[0].plot(t, 1000 * (e_ks - np.nanmean(e_tot)) / n, color="#1f6f8b", lw=1.2, label="E_KS (potential)")
    ax[0].plot(t, 1000 * (e_tot - np.nanmean(e_tot)) / n, color="#c0392b", lw=1.6, label="E_tot (conserved)")
    ax[0].set_xlabel("time (fs)")
    ax[0].set_ylabel("energy − ⟨E_tot⟩ (meV/atom)")
    ax[0].set_title("(a) Si, NVE: kinetic ↔ potential, E_tot flat")
    ax[0].legend(frameon=False, fontsize=8)

    t = series(nvt, "Time")
    temp = series(nvt, "Temp_ion")
    ax[1].plot(t, temp, color="#c0392b", lw=1.2, label="T (ions)")
    ax[1].axhline(350, color="k", ls="--", lw=0.8, label="target 350 K")
    nw = len(nvt[0])
    sig = 350 * np.sqrt(2 / (3 * nw - 3))
    ax[1].axhspan(350 - sig, 350 + sig, color="0.85", zorder=0, label="±σ canonical")
    ax[1].set_xlabel("time (fs)")
    ax[1].set_ylabel("temperature (K)")
    ax[1].set_title("(b) 8 H$_2$O, Nosé NVT")
    ax[1].legend(frameon=False, fontsize=8)

    t = series(npt, "Time")
    vol = np.array([f.get_volume() for f in npt])
    abc = np.array([f.cell.lengths() for f in npt])
    ax[2].plot(t, vol, color="#1f6f8b", lw=1.6, label="V (from each frame's cell)")
    ax[2].set_xlabel("time (fs)")
    ax[2].set_ylabel("cell volume (Å$^3$)", color="#1f6f8b")
    a2 = ax[2].twinx()
    for k, (name, c) in enumerate(zip("abc", ("#c0392b", "#e67e22", "#8e44ad"))):
        a2.plot(t, abc[:, k], color=c, lw=0.9, label=name)
    a2.set_ylabel("cell lengths (Å)")
    a2.legend(frameon=False, fontsize=8, loc="lower right")
    ax[2].set_title("(c) Si, Nosé–Parrinello–Rahman NPT")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "thermo.png"), dpi=150)
    plt.close(fig)

finally:
    shutil.rmtree(work)
print("figures/thermo.png")
