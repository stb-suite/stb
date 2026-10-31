#!/usr/bin/env python3
"""Draws the theory figures of README section 1 (no model, no MACE run needed):

  figures/cubic_anisotropy.png    the directional Young's modulus E(n) of cubic crystals with the same C11 and C12 as
                                  the silicon model and four Zener ratios A = 2 C44 / (C11 - C12): 0.5, 1 (isotropic),
                                  1.71 (the model's silicon) and 3, in the (001) and (1-10) planes.
  figures/vrh_bounds.png          the Voigt (uniform strain) and Reuss (uniform stress) shear moduli of the same family
                                  against A: they coincide only for an isotropic crystal; Hill is their mean.
  figures/internal_relaxation.png a crystal with an internal coordinate u (diamond under shear): at fixed strain the
                                  energy E(u) has its minimum away from the homogeneous position u = 0, and the relaxed
                                  constant C = C0 - D^2/K is smaller than the clamped one. Silicon model numbers.

    python3 theory_figures.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
C11, C12 = 144.9, 68.2          # GPa, the silicon model of the example


def cubic(c11, c12, c44):
    C = np.zeros((6, 6))
    C[:3, :3] = c12
    np.fill_diagonal(C[:3, :3], c11)
    C[3, 3] = C[4, 4] = C[5, 5] = c44
    return C


def youngs(C, n):
    S = np.linalg.inv(C)
    s = np.zeros((3, 3, 3, 3))
    for m, (i, j) in enumerate(VOIGT):
        for k_, (k, l) in enumerate(VOIGT):
            v = S[m, k_] / ((2.0 if m >= 3 else 1.0) * (2.0 if k_ >= 3 else 1.0))
            for a, b in ((i, j), (j, i)):
                for p, q in ((k, l), (l, k)):
                    s[a, b, p, q] = v
    return 1.0 / np.einsum("ijkl,...i,...j,...k,...l->...", s, n, n, n, n)


def cubic_anisotropy():
    t = np.linspace(0, 2 * np.pi, 361)
    planes = (("(001) plane: [100] at 0°, [110] at 45°", np.stack([np.cos(t), np.sin(t), 0 * t], 1)),
              ("(1-10) plane: [001] at 0°, [111] at 54.7°, [110] at 90°",
               np.stack([np.sin(t) / np.sqrt(2), np.sin(t) / np.sqrt(2), np.cos(t)], 1)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), subplot_kw={"projection": "polar"})
    for A, c in ((0.5, "#f59f00"), (1.0, "#868e96"), (1.71, "#2a6fdb"), (3.0, "#d9480f")):
        C = cubic(C11, C12, A * (C11 - C12) / 2)
        for ax, (_t, N) in zip(axes, planes):
            ax.plot(t, youngs(C, N), color=c, lw=2 if A == 1.71 else 1.4,
                    label=f"A = {A:g}" + (" (silicon model)" if A == 1.71 else " (isotropic)" if A == 1 else ""))
    for ax, (title, _N) in zip(axes, planes):
        ax.set_title(f"Young's modulus E (GPa)\n{title}", fontsize=9, pad=14)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(os.path.join(FIG, "cubic_anisotropy.png"), dpi=150)
    plt.close(fig)


def vrh_bounds():
    A = np.linspace(0.3, 4, 200)
    GV, GR, AU = [], [], []
    for a in A:
        C = cubic(C11, C12, a * (C11 - C12) / 2)
        S = np.linalg.inv(C)
        gv = ((C[0, 0] * 3) - (C[0, 1] * 3) + 3 * 3 * C[3, 3]) / 15
        gr = 15 / (4 * 3 * S[0, 0] - 4 * 3 * S[0, 1] + 3 * 3 * S[3, 3])
        GV.append(gv)
        GR.append(gr)
        AU.append(5 * gv / gr - 6 + 1)          # B_V = B_R for a cubic crystal
    GV, GR = np.array(GV), np.array(GR)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(A, GV, color="#d9480f", lw=2, label="Voigt (uniform strain): upper bound")
    axes[0].plot(A, GR, color="#2a6fdb", lw=2, label="Reuss (uniform stress): lower bound")
    axes[0].plot(A, (GV + GR) / 2, "k--", lw=1.3, label="Hill: the mean")
    axes[0].axvline(1, color="grey", lw=0.6)
    axes[0].axvline(1.71, color="#2f9e44", lw=1, ls=":")
    axes[0].text(1.75, GV.max() - 8, "silicon\nmodel", fontsize=8, color="#2f9e44")
    axes[0].set_xlabel("Zener ratio A = 2 C44 / (C11 − C12)")
    axes[0].set_ylabel("polycrystalline shear modulus G (GPa)")
    axes[0].set_title(f"Cubic crystals with C11 = {C11}, C12 = {C12} GPa", fontsize=10)
    axes[0].legend(fontsize=8, loc="lower right")
    axes[1].plot(A, AU, color="#862e9c", lw=2)
    axes[1].axvline(1, color="grey", lw=0.6)
    axes[1].axvline(1.71, color="#2f9e44", lw=1, ls=":")
    axes[1].set_xlabel("Zener ratio A")
    axes[1].set_ylabel("universal anisotropy index A^U")
    axes[1].set_title("A^U = 5 G_V / G_R + B_V / B_R − 6: zero only when isotropic", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "vrh_bounds.png"), dpi=150)
    plt.close(fig)


def internal_relaxation():
    """E(eps, u) = 1/2 C0 eps^2 - D eps u + 1/2 K u^2 per unit volume; u* = D eps / K; C = C0 - D^2 / K."""
    C0, Crel = 97.4, 65.5                       # GPa: clamped and relaxed C44 of the silicon model
    K = 1.0
    D = np.sqrt((C0 - Crel) * K)
    eps = 0.01
    u = np.linspace(-0.6, 1.6, 300) * D * eps / K
    E = 0.5 * C0 * eps ** 2 - D * eps * u + 0.5 * K * u ** 2
    ustar = D * eps / K
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.plot(u / ustar, E * 1e3, color="#2a6fdb", lw=2)
    ax.plot([0], [0.5 * C0 * eps ** 2 * 1e3], "o", color="#d9480f", ms=9, label="clamped: u = 0 (homogeneous strain)")
    ax.plot([1], [0.5 * Crel * eps ** 2 * 1e3], "o", color="#2f9e44", ms=9, label="relaxed: u = u*, the minimum")
    ax.set_xlabel("internal coordinate u / u*")
    ax.set_ylabel("elastic energy at fixed shear strain (arb.)")
    ax.set_title("At fixed strain the atoms find a lower energy inside the cell", fontsize=10)
    ax.legend(fontsize=8)
    ax = axes[1]
    e = np.linspace(-1, 1, 50)
    ax.plot(e, C0 * e / 100, color="#d9480f", lw=2, label=f"clamped ions: C44⁰ = {C0} GPa")
    ax.plot(e, Crel * e / 100, color="#2f9e44", lw=2, label=f"relaxed ions: C44 = {Crel} GPa")
    ax.set_xlabel("engineering shear strain (%)")
    ax.set_ylabel("shear stress (GPa)")
    ax.set_title("C44 = C44⁰ − D²/K: the relaxation softens it by 33 %", fontsize=10)
    ax.legend(fontsize=8)
    ax.axhline(0, color="grey", lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "internal_relaxation.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    cubic_anisotropy()
    vrh_bounds()
    internal_relaxation()
    print("figures/cubic_anisotropy.png, vrh_bounds.png, internal_relaxation.png")
