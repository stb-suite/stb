#!/usr/bin/env python3
"""Figures of the theory part of the 3.1 README (figures/):

    coordinates.png   (a) one atom in the oblique surface cell of si111_2x2.fdf: fractional (u, v) along the cell
                      vectors vs Cartesian (x, y); (b) periodic images, wrapping into [0, 1) and the minimum-image
                      distance between two atoms on opposite edges of the cell
    tolerance.png     how many symmetry operations spglib finds in bulk Si (si_bulk.fdf) whose atoms are displaced at
                      random by up to A Ang, for the three tolerances the report uses -- why a relaxed structure needs a
                      tolerance, and why too large a one invents symmetry

    python3 theory_figures.py       (numpy, spglib, matplotlib; seconds)
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spglib
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")


def read(name):
    return AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(os.path.join(HERE, name))))


def fig_coordinates():
    slab = read("si111_2x2.fdf")
    a, b = np.array(slab.cell[0][:2]), np.array(slab.cell[1][:2])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.4))

    # (a) fractional vs Cartesian
    for t in np.linspace(0, 1, 5):
        p0, p1 = t * a, t * a + b
        a1.plot([p0[0], p1[0]], [p0[1], p1[1]], color="0.85", lw=0.8)
        p0, p1 = t * b, t * b + a
        a1.plot([p0[0], p1[0]], [p0[1], p1[1]], color="0.85", lw=0.8)
    cell = np.array([[0, 0], a, a + b, b, [0, 0]])
    a1.plot(cell[:, 0], cell[:, 1], "k-", lw=1.5)
    a1.annotate("", xy=a, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=2))
    a1.annotate("", xy=b, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=2))
    a1.text(*(a * 0.55 + [0, -0.45]), "a", color="#1f77b4", fontsize=13, weight="bold")
    a1.text(*(b * 0.5 + [-0.6, 0]), "b", color="#1f77b4", fontsize=13, weight="bold")
    u, v = 0.625, 0.375
    r = u * a + v * b
    a1.plot([0, u * a[0]], [0, u * a[1]], color="#ff7f0e", lw=3)
    a1.plot([u * a[0], r[0]], [u * a[1], r[1]], color="#ff7f0e", lw=3, label=f"u a + v b  (u = {u}, v = {v})")
    a1.plot([0, r[0]], [r[1], r[1]], "--", color="#2ca02c", lw=1.4)
    a1.plot([r[0], r[0]], [0, r[1]], "--", color="#2ca02c", lw=1.4,
            label=f"x = {r[0]:.3f} Ang, y = {r[1]:.3f} Ang")
    a1.scatter(*r, s=180, c="#d62728", edgecolors="k", zorder=5)
    a1.set_aspect("equal")
    a1.set_xlabel("x (Ang)")
    a1.set_ylabel("y (Ang)")
    a1.set_title("(a) one position, two coordinates: r = u a + v b")
    a1.legend(loc="lower right", fontsize=9)

    # (b) periodic images and the minimum image
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            off = i * a + j * b
            c = cell + off
            a2.plot(c[:, 0], c[:, 1], "-", color="k" if (i, j) == (0, 0) else "0.75",
                          lw=1.5 if (i, j) == (0, 0) else 0.8)
    p = 0.05 * a + 0.5 * b          # near the left edge
    q = 0.95 * a + 0.5 * b          # near the right edge
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            off = i * a + j * b
            for pt, col in ((p, "#1f77b4"), (q, "#d62728")):
                a2.scatter(*(pt + off), s=60 if (i, j) == (0, 0) else 25, c=col,
                           alpha=1 if (i, j) == (0, 0) else 0.35, edgecolors="k", lw=0.4)
    a2.plot([p[0], q[0]], [p[1], q[1]], ":", color="0.4", lw=1.5)
    a2.text(*((p + q) / 2 + [-1.6, 0.3]), f"inside the cell: {np.linalg.norm(q - p):.2f} Ang", fontsize=9)
    qi = q - a
    a2.plot([p[0], qi[0]], [p[1], qi[1]], "-", color="#2ca02c", lw=2.5)
    a2.text(*(qi + [-1.6, -0.9]), f"minimum image: {np.linalg.norm(qi - p):.2f} Ang", color="#2ca02c", fontsize=9)
    w = 1.04 * a + 0.2 * b
    wr = w - a
    a2.scatter(*w, s=60, marker="s", c="#9467bd", edgecolors="k")
    a2.scatter(*wr, s=60, marker="s", c="#9467bd", edgecolors="k")
    a2.annotate("", xy=wr, xytext=w, arrowprops=dict(arrowstyle="-|>", color="#9467bd", lw=1.5))
    a2.text(*(w + [0.2, -0.5]), "u = 1.04 -> wrapped to 0.04", color="#9467bd", fontsize=9)
    a2.set_aspect("equal")
    a2.set_xlim(-6, 19)
    a2.set_ylim(-5, 11)
    a2.set_xlabel("x (Ang)")
    a2.set_ylabel("y (Ang)")
    a2.set_title("(b) periodic images, wrapping and the minimum image")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "coordinates.png"), dpi=150)
    plt.close(fig)


def fig_tolerance():
    si = read("si_bulk.fdf")
    amps = np.logspace(-4, -0.5, 36)
    tols = (1e-3, 1e-2, 1e-1)
    seeds = range(8)
    res = {t: [] for t in tols}
    for amp in amps:
        for t in tols:
            ops = []
            for sd in seeds:
                b = si.copy()
                b.positions += np.random.default_rng(sd).uniform(-amp, amp, b.positions.shape)
                ds = spglib.get_symmetry_dataset((np.array(b.cell), b.get_scaled_positions(), b.numbers), symprec=t)
                ops.append(len(ds.rotations) if ds is not None else 1)
            res[t].append(np.median(ops))
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for t, c in zip(tols, ("#1f77b4", "#ff7f0e", "#2ca02c")):
        ax.semilogx(amps, res[t], "o-", color=c, ms=4, label=f"symprec = {t:g} Ang")
        ax.axvline(t, color=c, ls=":", lw=1)
    ax.set_yscale("log", base=2)
    ax.set_yticks([1, 2, 4, 8, 16, 48, 192])
    ax.set_yticklabels(["1 (P1)", "2", "4", "8", "16", "48", "192 (Fd-3m)"])
    ax.set_xlabel("random displacement of every atom, up to A (Ang)")
    ax.set_ylabel("symmetry operations found (median of 8)")
    ax.set_title("Bulk Si (8-atom cell): what spglib finds, against noise and tolerance")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "tolerance.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_coordinates()
    fig_tolerance()
    print("figures written to", FIG)
