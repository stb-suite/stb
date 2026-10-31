#!/usr/bin/env python3
"""Draws the two theory figures of README section 1 (no MD run, no model needed):

  figures/verlet.png       velocity Verlet on a harmonic oscillator with silicon's highest phonon
                           frequency (15.5 THz, period 64.5 fs): the energy error is bounded and
                           oscillates (it conserves a nearby "shadow" energy), grows as (omega dt)^2,
                           and the integration becomes unstable at omega dt = 2 (dt ~ 20.5 fs here).
  figures/classical_vs_quantum.png
                           what a classical MD misses: the Debye heat capacity of silicon (theta_D = 645 K)
                           against the classical 3 kB per atom, and the mean-square vibration amplitude,
                           classical (proportional to T) against quantum (zero-point motion included).

    python3 theory_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import quad

HBAR = 1.054571817e-34      # J s
KB = 1.380649e-23           # J/K
AMU = 1.66053906660e-27     # kg
M_SI = 28.0855 * AMU
THETA_D = 645.0             # K, Debye temperature of Si (Kittel)
NU_MAX = 15.5e12            # Hz, Si optical phonon at Gamma (~15.5 THz)


def verlet_energy(omega, dt, n_steps, x0=1.0, v0=0.0):
    """Velocity Verlet for x'' = -omega^2 x; returns t and E(t)/E(0)."""
    x, v = x0, v0
    e0 = 0.5 * v0 ** 2 + 0.5 * omega ** 2 * x0 ** 2
    t, e = [0.0], [1.0]
    a = -omega ** 2 * x
    for i in range(n_steps):
        v_half = v + 0.5 * dt * a
        x = x + dt * v_half
        a = -omega ** 2 * x
        v = v_half + 0.5 * dt * a
        t.append((i + 1) * dt)
        e.append((0.5 * v ** 2 + 0.5 * omega ** 2 * x ** 2) / e0)
        if abs(e[-1]) > 1e6:
            break
    return np.array(t), np.array(e)


def figure_verlet():
    omega = 2 * np.pi * NU_MAX * 1e-15          # rad/fs
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    colors = ["#0072B2", "#009E73", "#E69F00", "#D55E00"]
    for dt, c in zip([1.0, 3.0, 8.0, 16.0], colors):
        t, e = verlet_energy(omega, dt, int(400 / dt))
        ax[0].plot(t, e - 1, color=c, lw=1.2, marker="o" if dt > 4 else None, ms=3,
                   label=f"dt = {dt:g} fs (ω·dt = {omega * dt:.2f})")
    t, e = verlet_energy(omega, 21.0, 25)
    ax[0].plot(t, e - 1, color="#555555", ls="--", lw=1.4, label=f"dt = 21 fs (ω·dt = {omega * 21:.2f}): unstable")
    ax[0].set_ylim(-0.6, 1.4)
    ax[0].set_xlim(0, 400)
    ax[0].set_xlabel("t (fs)")
    ax[0].set_ylabel("E(t)/E(0) − 1")
    ax[0].set_title("Velocity Verlet, harmonic mode of 15.5 THz (period 64.5 fs)")
    ax[0].legend(fontsize=8, loc="upper left")
    ax[0].axhline(0, color="#999999", lw=0.8)

    wdt = np.linspace(0.02, 1.98, 200)
    amp = []
    for w in wdt:
        _, e = verlet_energy(1.0, w, int(200 / w) + 50)
        amp.append(e.max() - e.min())
    ax[1].loglog(wdt, amp, color="#0072B2", lw=2, label="numerical: max − min of E(t)/E(0)")
    ax[1].loglog(wdt, wdt ** 2 / 4, color="#D55E00", ls="--", lw=1.4, label="(ω·dt)² / 4")
    for dt in [0.5, 1, 2, 5]:
        w = omega * dt
        ax[1].axvline(w, color="#999999", lw=0.8, ls=":")
        ax[1].text(w, 2e-4, f" {dt:g} fs", rotation=90, fontsize=8, va="bottom", color="#555555")
    ax[1].axvline(2.0, color="#555555", lw=1.2)
    ax[1].text(1.95, 3e-3, "stability limit ω·dt = 2  ", ha="right", fontsize=8)
    ax[1].set_xlabel("ω·dt")
    ax[1].set_ylabel("relative energy oscillation")
    ax[1].set_title("Error grows as dt²; silicon's fastest mode at 1, 2, 5 fs")
    ax[1].legend(fontsize=8, loc="upper left")
    ax[1].set_ylim(1e-4, 10)
    fig.tight_layout()
    fig.savefig("figures/verlet.png", dpi=130)
    plt.close(fig)


def debye_cv(t):
    """Debye C_V / (3 N kB)."""
    x_d = THETA_D / t
    integral = quad(lambda x: x ** 4 * np.exp(x) / np.expm1(x) ** 2, 0, x_d)[0]
    return 3.0 / x_d ** 3 * integral


def debye_u2(t, quantum=True):
    """Mean-square displacement along one axis (Ang^2) in the Debye model."""
    pref = 3 * HBAR ** 2 / (M_SI * KB * THETA_D)
    if not quantum:
        return pref * t / THETA_D * 1e20          # high-T limit: 3 hbar^2 T / (m kB theta^2)
    x_d = THETA_D / t
    integral = quad(lambda x: x / np.expm1(x), 0, x_d)[0]
    return pref * ((t / THETA_D) ** 2 * integral + 0.25) * 1e20


def figure_classical_vs_quantum():
    temps = np.linspace(10, 1600, 300)
    cv = np.array([debye_cv(t) for t in temps])
    u2_q = np.array([debye_u2(t) for t in temps])
    u2_c = np.array([debye_u2(t, quantum=False) for t in temps])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    ax[0].plot(temps, 3 * cv, color="#0072B2", lw=2, label=f"quantum (Debye, θ_D = {THETA_D:.0f} K)")
    ax[0].axhline(3, color="#D55E00", ls="--", lw=1.6, label="classical MD, harmonic: 3 kB (Dulong–Petit)")
    for t in [300, 1200]:
        ax[0].plot([t], [3 * debye_cv(t)], "o", color="#0072B2")
        ax[0].annotate(f"{3 * debye_cv(t):.2f} kB at {t} K", (t, 3 * debye_cv(t)), textcoords="offset points",
                       xytext=(8, -14), fontsize=8)
    ax[0].set_xlabel("T (K)")
    ax[0].set_ylabel("C_V per atom (kB)")
    ax[0].set_ylim(0, 3.4)
    ax[0].set_title("Heat capacity: classical MD is right only for T ≳ θ_D")
    ax[0].legend(fontsize=8, loc="lower right")

    ax[1].plot(temps, u2_q, color="#0072B2", lw=2, label="quantum (with zero-point motion)")
    ax[1].plot(temps, u2_c, color="#D55E00", ls="--", lw=1.6, label="classical (∝ T)")
    for t in [300]:
        ax[1].annotate(f"300 K: quantum/classical = {debye_u2(t) / debye_u2(t, False):.2f}",
                       (t, debye_u2(t)), textcoords="offset points", xytext=(40, -40), fontsize=8,
                     arrowprops=dict(arrowstyle="->", color="#555555"))
    ax[1].set_xlabel("T (K)")
    ax[1].set_ylabel("<u_x²> per atom (Å²)")
    ax[1].set_title("Vibration amplitude of a Si atom (Debye model)")
    ax[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig("figures/classical_vs_quantum.png", dpi=130)
    plt.close(fig)
    for t in [100, 300, 600, 1200]:
        print(f"T = {t:5d} K: C_V = {3 * debye_cv(t):.3f} kB/atom, <u_x^2> quantum {debye_u2(t):.5f}, "
              f"classical {debye_u2(t, False):.5f} Ang^2 (ratio {debye_u2(t) / debye_u2(t, False):.3f})")


if __name__ == "__main__":
    figure_verlet()
    figure_classical_vs_quantum()
    print("figures/verlet.png, figures/classical_vs_quantum.png written")
