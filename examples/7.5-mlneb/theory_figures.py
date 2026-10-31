#!/usr/bin/env python3
"""Figures of section 1 of the 7.5 README (numpy + matplotlib only).

Every NEB here is computed, not drawn: a small chain-of-states optimizer on the Mueller-Brown
surface, a standard two-dimensional test potential with three minima and two saddle points, written
with the same force projections stb-mlneb uses through ASE (improved tangent, climbing image).

    python3 theory_figures.py          # writes figures/neb_*.png, figures/rate_sensitivity.png
"""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")

# ---- Mueller-Brown potential (K. Mueller, L. D. Brown, Theor. Chim. Acta 53, 75 (1979)) ----------
A = np.array([-200.0, -100.0, -170.0, 15.0])
a = np.array([-1.0, -1.0, -6.5, 0.7])
b = np.array([0.0, 0.0, 11.0, 0.6])
c = np.array([-10.0, -10.0, -6.5, 0.7])
X0 = np.array([1.0, 0.0, -0.5, -1.0])
Y0 = np.array([0.0, 0.5, 1.5, 1.0])
SCALE = 0.01   # energies in "eV": the original units / 100, so the barriers look like atomic ones


def V(x, y):
    x, y = np.asarray(x)[..., None], np.asarray(y)[..., None]
    dx, dy = x - X0, y - Y0
    return SCALE * np.sum(A * np.exp(a * dx ** 2 + b * dx * dy + c * dy ** 2), axis=-1)


def grad(p):
    x, y = p[..., 0:1], p[..., 1:2]
    dx, dy = x - X0, y - Y0
    e = A * np.exp(a * dx ** 2 + b * dx * dy + c * dy ** 2)
    gx = np.sum(e * (2 * a * dx + b * dy), axis=-1)
    gy = np.sum(e * (b * dx + 2 * c * dy), axis=-1)
    return SCALE * np.stack([gx, gy], axis=-1)


MIN_A = np.array([-0.5582, 1.4417])
MIN_B = np.array([0.6235, 0.0280])
MIN_C = np.array([-0.0500, 0.4667])
SAD_1 = np.array([-0.8220, 0.6243])
SAD_2 = np.array([0.2123, 0.2930])


def refine(p, saddle=False):
    """Newton on the gradient: exact stationary points of the surface."""
    for _ in range(50):
        h = 1e-5
        H = np.zeros((2, 2))
        for j in range(2):
            d = np.zeros(2)
            d[j] = h
            H[:, j] = (grad(p + d) - grad(p - d)) / (2 * h)
        p = p - np.linalg.solve(H, grad(p))
    return p


MIN_A, MIN_B, MIN_C = refine(MIN_A), refine(MIN_B), refine(MIN_C)
SAD_1, SAD_2 = refine(SAD_1, True), refine(SAD_2, True)


# ---- a chain-of-states optimizer ---------------------------------------------------------------

def tangents(R, E):
    """Improved tangent (Henkelman & Jonsson 2000): toward the higher neighbour; at an extremum, a
    bisection weighted by the energy differences."""
    n = len(R)
    T = np.zeros_like(R)
    for i in range(1, n - 1):
        tp, tm = R[i + 1] - R[i], R[i] - R[i - 1]
        if E[i + 1] > E[i] > E[i - 1]:
            t = tp
        elif E[i + 1] < E[i] < E[i - 1]:
            t = tm
        else:
            dmax = max(abs(E[i + 1] - E[i]), abs(E[i - 1] - E[i]))
            dmin = min(abs(E[i + 1] - E[i]), abs(E[i - 1] - E[i]))
            t = tp * dmax + tm * dmin if E[i + 1] > E[i - 1] else tp * dmin + tm * dmax
        T[i] = t / np.linalg.norm(t)
    return T


def band_forces(R, kind, k, climb=False):
    E = V(R[:, 0], R[:, 1])
    G = grad(R)
    F = np.zeros_like(R)
    T = tangents(R, E)
    top = 1 + int(np.argmax(E[1:-1]))
    for i in range(1, len(R) - 1):
        if kind == "plain":   # elastic band: full true force + full spring force
            F[i] = -G[i] + k * (R[i + 1] - 2 * R[i] + R[i - 1])
            continue
        t = T[i]
        f_perp = -G[i] + np.dot(G[i], t) * t
        f_spring = k * (np.linalg.norm(R[i + 1] - R[i]) - np.linalg.norm(R[i] - R[i - 1])) * t
        if climb and i == top:
            F[i] = -G[i] + 2 * np.dot(G[i], t) * t
        else:
            F[i] = f_perp + f_spring
    return F, E, T


def relax_band(R, kind="neb", k=1.0, climb=False, steps=20000, fmax=1e-3, dt=0.002):
    """FIRE (Bitzek et al. 2006) on the band; returns the band and the residual-force history."""
    R = R.copy()
    v = np.zeros_like(R)
    alpha, n_pos, dtmax = 0.1, 0, 0.02
    hist = []
    for _ in range(steps):
        F, E, _ = band_forces(R, kind, k, climb)
        res = np.abs(F[1:-1]).max()
        hist.append(res)
        if res < fmax:
            break
        P = np.sum(F * v)
        if P > 0:
            v = (1 - alpha) * v + alpha * F * np.linalg.norm(v) / max(np.linalg.norm(F), 1e-12)
            n_pos += 1
            if n_pos > 5:
                dt, alpha = min(dt * 1.1, dtmax), alpha * 0.99
        else:
            v[:] = 0.0
            dt, alpha, n_pos = dt * 0.5, 0.1, 0
        v += dt * F
        R[1:-1] += dt * v[1:-1]
    return R, np.array(hist)


def path_s(R):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(R, axis=0), axis=1))])


def hermite_spline(R, n=60):
    """Energy along the band as a piecewise cubic with the image energies AND the projected forces
    (dE/ds = -F . tangent) at each image -- the construction of ASE's fit_images."""
    E = V(R[:, 0], R[:, 1])
    G = grad(R)
    s = path_s(R)
    d = np.diff(R, axis=0)
    # tangent at each image: central difference of the path (ends: one-sided)
    T = np.zeros_like(R)
    T[0], T[-1] = d[0], d[-1]
    T[1:-1] = R[2:] - R[:-2]
    T /= np.linalg.norm(T, axis=1)[:, None]
    dEds = np.sum(G * T, axis=1)
    ss, ee = [], []
    for i in range(len(R) - 1):
        h = s[i + 1] - s[i]
        u = np.linspace(0, 1, n, endpoint=(i == len(R) - 2))
        h00, h10, h01, h11 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u, -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
        ee.append(h00 * E[i] + h10 * h * dEds[i] + h01 * E[i + 1] + h11 * h * dEds[i + 1])
        ss.append(s[i] + u * h)
    return np.concatenate(ss), np.concatenate(ee), s, E


def linear_band(p, q, n):
    return np.array([p + (q - p) * t for t in np.linspace(0, 1, n)])


# ---- figures -----------------------------------------------------------------------------------

def contour(ax, lim=(-1.5, 1.2, -0.3, 2.0)):
    x = np.linspace(lim[0], lim[1], 300)
    y = np.linspace(lim[2], lim[3], 300)
    Xg, Yg = np.meshgrid(x, y)
    Z = V(Xg, Yg)
    Z = np.minimum(Z, 1.0)
    ax.contourf(Xg, Yg, Z, levels=np.linspace(-1.5, 1.0, 26), cmap="Greys_r", alpha=0.55)
    ax.contour(Xg, Yg, Z, levels=np.linspace(-1.5, 1.0, 26), colors="k", linewidths=0.3, alpha=0.5)
    ax.set_xlim(lim[0], lim[1])
    ax.set_ylim(lim[2], lim[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def mark_points(ax):
    for p, lab in ((MIN_A, "A"), (MIN_C, "C"), (MIN_B, "B")):
        ax.plot(*p, "o", color="#2a6fdb", ms=6, zorder=6)
        ax.annotate(lab, p, xytext=(6, 4), textcoords="offset points", fontsize=10, color="#2a6fdb", weight="bold")
    for p, lab in ((SAD_1, "S1"), (SAD_2, "S2")):
        ax.plot(*p, "x", color="#c92a2a", ms=9, mew=2.2, zorder=6)
        ax.annotate(lab, p, xytext=(6, -12), textcoords="offset points", fontsize=10, color="#c92a2a", weight="bold")


def fig_bands():
    n = 13
    R0 = linear_band(MIN_A, MIN_B, n)
    plain_stiff, _ = relax_band(R0, "plain", k=8.0, steps=40000, fmax=2e-3)
    plain_soft, _ = relax_band(R0, "plain", k=0.05, steps=40000, fmax=2e-3)
    neb, h1 = relax_band(R0, "neb", k=1.0, fmax=1e-3)
    ci, h2 = relax_band(neb, "neb", k=1.0, climb=True, fmax=1e-3)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0))
    ax = axes[0]
    contour(ax)
    ax.plot(R0[:, 0], R0[:, 1], "o--", color="#868e96", ms=4, lw=1, label="linear interpolation (start)")
    ax.plot(plain_stiff[:, 0], plain_stiff[:, 1], "s-", color="#e8590c", ms=4, lw=1.2,
            label="plain elastic band, stiff springs")
    ax.plot(plain_soft[:, 0], plain_soft[:, 1], "^-", color="#ae3ec9", ms=4, lw=1.2,
            label="plain elastic band, soft springs")
    mark_points(ax)
    ax.set_title("(a) Without nudging: corner cutting or sliding down")
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.9)

    ax = axes[1]
    contour(ax)
    ax.plot(neb[:, 0], neb[:, 1], "o-", color="#2a6fdb", ms=4, lw=1.2, label="NEB")
    ax.plot(ci[:, 0], ci[:, 1], "o-", color="#c92a2a", ms=3, lw=0.8, label="climbing-image NEB")
    top = 1 + int(np.argmax(V(ci[1:-1, 0], ci[1:-1, 1])))
    ax.plot(*ci[top], "*", color="#c92a2a", ms=15, mec="k", zorder=7, label="climbing image = S1")
    mark_points(ax)
    ax.set_title("(b) NEB follows the minimum-energy path")
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.9)

    ax = axes[2]
    for R, col, lab, ms in ((neb, "#2a6fdb", "NEB", 5), (ci, "#c92a2a", "climbing-image NEB", 4)):
        ss, ee, s, E = hermite_spline(R)
        ax.plot(ss, ee - E[0], "-", color=col, lw=1.4)
        ax.plot(s, E - E[0], "o", color=col, ms=ms, label=f"{lab}: highest image {np.max(E) - E[0]:.3f}")
    eS1 = V(*SAD_1) - V(*MIN_A)
    eS2 = V(*SAD_2) - V(*MIN_A)
    eC = V(*MIN_C) - V(*MIN_A)
    ax.axhline(eS1, color="#c92a2a", ls=":", lw=1)
    ax.text(0.02, eS1 + 0.02, f"exact saddle S1: {eS1:.3f}", color="#c92a2a", fontsize=8, transform=ax.get_yaxis_transform())
    ax.axhline(eS2, color="#868e96", ls=":", lw=1)
    ax.text(0.62, eS2 + 0.02, f"S2: {eS2:.3f}", color="#495057", fontsize=8, transform=ax.get_yaxis_transform())
    ax.axhline(eC, color="#868e96", ls=":", lw=1)
    ax.text(0.40, eC - 0.07, f"intermediate minimum C: {eC:.3f}", color="#495057", fontsize=8,
            transform=ax.get_yaxis_transform())
    ax.set_xlabel("reaction coordinate s (path length)")
    ax.set_ylabel("E - E(A)")
    ax.set_title("(c) Energy along the band, force-consistent spline")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "neb_mueller_brown.png"), dpi=150)
    plt.close(fig)
    print(f"neb_mueller_brown.png: exact S1 {eS1:.4f}; NEB highest image "
          f"{V(neb[:, 0], neb[:, 1]).max() - V(*MIN_A):.4f}; CI {V(ci[:, 0], ci[:, 1]).max() - V(*MIN_A):.4f}; "
          f"CI image at {ci[top]}, S1 at {SAD_1}; steps {len(h1)} + {len(h2)}")


def fig_forces():
    """The force decomposition on one image of a band that is not converged yet."""
    n = 9
    R0 = linear_band(MIN_A, MIN_B, n)
    R, _ = relax_band(R0, "neb", k=1.0, steps=25)
    E = V(R[:, 0], R[:, 1])
    T = tangents(R, E)
    G = grad(R)
    i = 2
    t = T[i]
    ftrue = -G[i]
    fperp = ftrue - np.dot(ftrue, t) * t
    fpar_true = np.dot(ftrue, t) * t
    fspring = 1.0 * (np.linalg.norm(R[i + 1] - R[i]) - np.linalg.norm(R[i] - R[i - 1])) * t
    sc = 0.35 / np.linalg.norm(ftrue)

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    contour(ax, lim=(-1.15, -0.05, 0.2, 1.35))
    ax.plot(R[:, 0], R[:, 1], "o-", color="#495057", ms=6, lw=1.2, label="band (images)")
    p = R[i]
    kw = dict(angles="xy", scale_units="xy", scale=1, width=0.006, zorder=8)
    ax.quiver(*p, *(ftrue * sc), color="#212529", **kw, label=r"true force $-\nabla E$")
    ax.quiver(*p, *(fperp * sc), color="#2a6fdb", **kw, label=r"kept: $-\nabla E_\perp$ (moves the band sideways)")
    ax.quiver(*p, *(fpar_true * sc), color="#adb5bd", **kw, label=r"discarded: $-\nabla E_\parallel$ (would slide down)")
    ax.quiver(*p, *(t * 0.18), color="#2b8a3e", **kw, label=r"tangent $\hat\tau$")
    if np.linalg.norm(fspring) > 1e-6:
        ax.quiver(*p, *(fspring / np.linalg.norm(fspring) * 0.12), color="#e8590c", **kw,
                  label=r"kept: spring force $\parallel\hat\tau$ (length not to scale)")
    ax.annotate(f"image {i}", p, xytext=(-55, 12), textcoords="offset points", fontsize=9)
    mark_points(ax)
    ax.set_title("Nudging: the two projections of the NEB force")
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "neb_forces.png"), dpi=150)
    plt.close(fig)
    print("neb_forces.png")


def fig_spline():
    """Few images: the highest image misses the saddle, the force-consistent spline finds it."""
    n = 5
    R0 = linear_band(MIN_A, MIN_C, n)
    R, _ = relax_band(R0, "neb", k=1.0, fmax=1e-3)
    ss, ee, s, E = hermite_spline(R, n=200)
    eS1 = V(*SAD_1) - E[0]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    s_lin = np.linspace(0, 1, 2)
    for j in range(len(s) - 1):
        ax.plot(s[j:j + 2], E[j:j + 2] - E[0], "--", color="#adb5bd", lw=1, label="straight segments" if j == 0 else None)
    ax.plot(ss, ee - E[0], "-", color="#2a6fdb", lw=1.6, label=f"force-consistent spline: max {ee.max() - E[0]:.3f}")
    ax.plot(s, E - E[0], "o", color="#2a6fdb", ms=7, label=f"images (no climbing): highest {E.max() - E[0]:.3f}")
    G = grad(R)
    d = np.diff(R, axis=0)
    T = np.zeros_like(R)
    T[0], T[-1] = d[0], d[-1]
    T[1:-1] = R[2:] - R[:-2]
    T /= np.linalg.norm(T, axis=1)[:, None]
    for j in range(len(R)):
        slope = np.dot(G[j], T[j])
        ds = 0.06
        ax.plot([s[j] - ds, s[j] + ds], [E[j] - E[0] - slope * ds, E[j] - E[0] + slope * ds], "-", color="#e8590c",
                lw=2.5, label=r"slope $dE/ds = -F\cdot\hat\tau$ at each image" if j == 0 else None)
    ax.axhline(eS1, color="#c92a2a", ls=":", lw=1)
    ax.text(0.02, eS1 + 0.015, f"exact saddle: {eS1:.3f}", color="#c92a2a", fontsize=8,
            transform=ax.get_yaxis_transform())
    ax.set_xlabel("reaction coordinate s")
    ax.set_ylabel("E - E(A)")
    ax.set_title("5 images from A to C: the spline uses the forces too")
    ax.legend(fontsize=7.5, loc="lower center")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "neb_spline.png"), dpi=150)
    plt.close(fig)
    print(f"neb_spline.png: highest image {E.max() - E[0]:.4f}, spline {ee.max() - E[0]:.4f}, exact {eS1:.4f}")


def fig_rate():
    kB = 8.617333e-5
    nu = 1e13
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axes[0]
    invT = np.linspace(0.5, 4.0, 200)   # 1000/T
    T = 1000 / invT
    for Ea, col, lab in ((0.25, "#2b8a3e", "SIESTA on the path: 0.25 eV"), (0.31, "#2a6fdb", "round-3 model: 0.31 eV"),
                         (0.72, "#c92a2a", "round-2 model: 0.72 eV")):
        ax.semilogy(invT, nu * np.exp(-Ea / (kB * T)), color=col, lw=2, label=lab)
    ax.axvline(1000 / 300, color="#868e96", ls=":", lw=1)
    ax.text(1000 / 300 + 0.03, 1e11, "300 K", color="#495057", fontsize=8)
    ax.set_xlabel("1000 / T (1/K)")
    ax.set_ylabel(r"hop rate $\Gamma = \nu^* e^{-E_a/k_BT}$ (1/s), $\nu^* = 10^{13}$ Hz")
    ax.set_title("(a) Arrhenius: the barrier is the slope")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    dE = np.linspace(0, 0.5, 200)
    for T_, col in ((300, "#c92a2a"), (600, "#e8590c"), (1000, "#2a6fdb"), (1500, "#2b8a3e")):
        ax.semilogy(dE, np.exp(dE / (kB * T_)), color=col, lw=2, label=f"{T_} K")
    for d, lab in ((0.057, "round 3: +0.057 eV"), (0.468, "round 2: +0.468 eV")):
        ax.axvline(d, color="#495057", ls=":", lw=1)
        ax.text(d + 0.005, 3e7, lab, rotation=90, fontsize=8, color="#495057", va="center")
    ax.set_xlabel("error of the barrier (eV)")
    ax.set_ylabel(r"factor on the rate, $e^{\delta E_a/k_BT}$")
    ax.set_title("(b) What a barrier error does to the rate")
    ax.legend(fontsize=8, title="temperature")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "rate_sensitivity.png"), dpi=150)
    plt.close(fig)
    for d in (0.057, 0.1, 0.468):
        print(f"rate factor for {d} eV: 300 K {np.exp(d / (kB * 300)):.3g}, 1000 K {np.exp(d / (kB * 1000)):.3g}")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    print(f"Mueller-Brown (x0.01): A {MIN_A} {V(*MIN_A):.4f}; C {MIN_C} {V(*MIN_C):.4f}; B {MIN_B} {V(*MIN_B):.4f}; "
          f"S1 {SAD_1} {V(*SAD_1):.4f}; S2 {SAD_2} {V(*SAD_2):.4f}")
    fig_bands()
    fig_forces()
    fig_spline()
    fig_rate()
