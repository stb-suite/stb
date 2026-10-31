#!/usr/bin/env python3
"""Draws the theory figures of README section 1 (no model, no phonon run needed):

  figures/chain_dispersion.png   the 1D diatomic chain (two masses, one spring constant): an acoustic branch that
                                 starts linearly at q = 0 (a sound wave) and an optical branch in which the two
                                 atoms of a cell move against each other; equal masses close the gap.
  figures/single_mode.png        the thermodynamics of ONE quantum harmonic oscillator of frequency nu, against
                                 T / (h nu / k_B): energy (with its zero-point part), heat capacity, entropy, free
                                 energy, each against the classical value. A crystal is 3N such modes (Eq. 9-12).
  figures/supercell_folding.png  what a finite supercell sees: a 1D chain whose force constants reach 6 neighbours,
                                 computed in supercells of 2, 3 and 4 cells. Each supercell gives the exact frequency
                                 at the q-points commensurate with it (q = m/N) and an interpolation in between.
  figures/displacement_error.png the finite-difference trade-off: the error of a force constant from a displacement
                                 u, with an anharmonic term growing as u (one-sided) or u^2 (central difference) and
                                 a noise term sigma/u. Silicon numbers.

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


def chain_dispersion():
    """omega^2 = K (1/m1 + 1/m2) -+ K sqrt((1/m1 + 1/m2)^2 - 4 sin^2(qa/2) / (m1 m2)), nearest-neighbour spring K."""
    q = np.linspace(-np.pi, np.pi, 401)          # q a
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, (m1, m2, title) in zip(axes, [(1.0, 1.0, "equal masses (like Si)"), (1.0, 2.5, "unequal masses (like GaAs, AlP)")]):
        K = 1.0
        s = 1 / m1 + 1 / m2
        root = np.sqrt(s ** 2 - 4 * np.sin(q / 2) ** 2 / (m1 * m2))
        ac, op = np.sqrt(K * (s - root)), np.sqrt(K * (s + root))
        ax.plot(q / np.pi, ac, color="#2a6fdb", lw=2, label="acoustic: the cell moves as a whole")
        ax.plot(q / np.pi, op, color="#d9480f", lw=2, label="optical: the two atoms move against each other")
        v = np.sqrt(K / (2 * (m1 + m2)))             # slope d omega / d(qa) at q -> 0: the sound velocity (a = 1)
        qs = np.linspace(0, 0.35, 10)
        ax.plot(qs, v * qs * np.pi, "k:", lw=1.2, label="sound velocity: slope at q = 0")
        ax.set_xlabel("q (units of π/a)")
        ax.set_ylabel("ω (units of √(K/m₁))")
        ax.set_title(title)
        ax.set_xlim(-1, 1)
        ax.set_ylim(0, None)
        ax.axvline(0, color="grey", lw=0.5)
        if m1 != m2:
            ax.annotate("gap", xy=(1, 0.5 * (np.sqrt(2 * K / m2) + np.sqrt(2 * K / m1))), xytext=(0.55, 1.2),
                        arrowprops=dict(arrowstyle="->"))
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(FIG, "chain_dispersion.png"), dpi=150)
    plt.close(fig)


def single_mode():
    x = np.linspace(0.02, 3.0, 400)               # T / theta, theta = h nu / k_B
    y = 1 / x                                    # h nu / k_B T
    n = 1 / np.expm1(y)
    E = x * y * (0.5 + n)                        # in units of k_B theta: (1/2 + n) h nu
    Cv = y ** 2 * np.exp(y) / np.expm1(y) ** 2
    S = y * n - np.log(-np.expm1(-y))
    F = x * y * 0.5 + x * np.log(-np.expm1(-y))
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]
    ax.plot(x, E, color="#2a6fdb", lw=2, label="quantum: (½ + n) hν")
    ax.plot(x, x, "k--", lw=1, label="classical: k_B T")
    ax.axhline(0.5, color="grey", lw=0.8, ls=":")
    ax.text(0.05, 0.53, "zero-point energy ½hν", fontsize=8)
    ax.set_title("energy of the mode")
    ax.set_ylabel("E / (k_B θ)")
    ax = axes[1]
    ax.plot(x, Cv, color="#2a6fdb", lw=2, label="quantum")
    ax.axhline(1, color="k", ls="--", lw=1, label="classical (Dulong-Petit): k_B")
    for xv in (300 / 745, 1000 / 745):
        ax.axvline(xv, color="#d9480f", lw=0.8, ls=":")
    ax.text(300 / 745 + 0.03, 0.15, "Si optical mode\nat 300 K", fontsize=7, color="#d9480f")
    ax.set_title("heat capacity")
    ax.set_ylabel("C_V / k_B")
    ax = axes[2]
    ax.plot(x, S, color="#2a6fdb", lw=2, label="entropy S / k_B")
    ax.plot(x, F, color="#7048e8", lw=2, label="free energy F / (k_B θ)")
    ax.plot(x, -x * np.log(x), "k--", lw=1, label="classical F = k_B T ln(hν / k_B T)")
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_title("entropy and free energy")
    for a in axes:
        a.set_xlabel("T / θ,   θ = hν / k_B")
        a.legend(fontsize=8)
        a.set_xlim(0, 3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "single_mode.png"), dpi=150)
    plt.close(fig)


def supercell_folding():
    """A 1D monatomic chain with force constants to 6 neighbours (decaying, alternating sign). In an N-cell supercell
    the finite-displacement method sees the SUM of each force constant over its images, Phi_N(j) = sum_k Phi(j + kN),
    and phonopy interpolates the dispersion from it. Exact at q = m/N."""
    K = {1: 1.0, 2: -0.18, 3: 0.09, 4: -0.05, 5: 0.03, 6: -0.015}

    def omega2(q, fc):                           # q in units of 2 pi / a
        return sum(2 * k * (1 - np.cos(2 * np.pi * q * j)) for j, k in fc.items())

    q = np.linspace(0, 0.5, 400)
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.plot(q, np.sqrt(omega2(q, K)), color="k", lw=2.5, label="exact (force constants to 6 neighbours)")
    colors = {2: "#d9480f", 3: "#2f9e44", 4: "#2a6fdb"}
    for N, c in colors.items():
        # folded force constants: the pair of neighbours at +-j lands on the supercell atoms j mod N and -j mod N;
        # phonopy places each at its minimum image, so the pair coefficient at distance d = min(j mod N, N - j mod N)
        # collects every K_j with j = +-d (mod N). At d = N/2 one atom receives both +j and -j, split equally between
        # its two equidistant images -- the same pair coefficient again.
        fcN = {}
        for j, k in K.items():
            jj = j % N
            if jj == 0:
                continue                         # an image of the atom itself: (1 - cos) vanishes at every q it can see
            jj = min(jj, N - jj)
            fcN[jj] = fcN.get(jj, 0.0) + k
        w = np.sqrt(np.clip(omega2(q, fcN), 0, None))
        ax.plot(q, w, color=c, lw=1.4, ls="--", label=f"{N}-cell supercell")
        qc = np.arange(0, N // 2 + 1) / N
        ax.plot(qc, np.sqrt(omega2(qc, K)), "o", color=c, ms=7, mfc="none", mew=1.8)
    ax.set_xlabel("q (units of 2π/a)")
    ax.set_ylabel("ω (arb. units)")
    ax.set_title("A supercell is exact only at the q-points commensurate with it (circles)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "supercell_folding.png"), dpi=150)
    plt.close(fig)


def displacement_error():
    """Relative error of a force constant k from a displacement u: anharmonic (one-sided ~ c3 u / k, central ~ c4 u^2 / k)
    plus force noise sigma / (k u). Silicon-like numbers: k ~ 13 eV/Ang^2 (the on-site constant), a cubic term giving a
    1 % shift at u = 0.1 Ang, noise 3 meV/Ang (SIESTA in section 5.9) or 0.1 meV/Ang (float64 MACE)."""
    u = np.logspace(-3.3, -0.5, 300)
    k = 13.0
    a3 = 0.10                                   # relative error per Ang, one-sided
    a4 = 0.6                                    # relative error per Ang^2, central
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for sigma, ls, lab in ((3e-3, "-", "noise 3 meV/Å (SIESTA, section 5.9)"), (1e-4, "--", "noise 0.1 meV/Å (MACE, float64)")):
        noise = sigma / (k * u)
        ax.loglog(u, a3 * u + noise, color="#d9480f", ls=ls, lw=1.8, label=f"one-sided, {lab}")
        ax.loglog(u, a4 * u ** 2 + noise, color="#2a6fdb", ls=ls, lw=1.8, label=f"central (±u), {lab}")
    ax.axvline(0.01, color="grey", lw=0.8, ls=":")
    ax.text(0.0105, 6e-4, "default 0.01 Å", fontsize=8, color="grey")
    ax.set_xlabel("displacement u (Å)")
    ax.set_ylabel("relative error of the force constant")
    ax.set_title("Finite differences: anharmonicity grows with u, noise with 1/u")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "displacement_error.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    chain_dispersion()
    single_mode()
    supercell_folding()
    displacement_error()
    print("figures/chain_dispersion.png, single_mode.png, supercell_folding.png, displacement_error.png")
