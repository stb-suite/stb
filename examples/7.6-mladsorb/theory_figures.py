#!/usr/bin/env python3
"""Figures of the theory part of the 7.6 README (figures/):

    si111_sites.png   the bulk-terminated Si(111) slab of the example: side view (bilayers, where it is cut, the
                      dangling bonds) and top view with the four high-symmetry adatom sites T1, T4, H3, B2 -- taken
                      from stb-mladsorb's own site finder, so the labels are those of the report
    populations.png   Boltzmann population of the most stable of two sites against their energy difference, and how
                      a multiplicity ratio shifts it

    python3 theory_figures.py      (needs stb_suite; seconds, no model)
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io
from stb import mladsorb as m

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
KB = 8.617333262e-5
NAMES = {1: "T1", 2: "T4", 3: "H3", 5: "B2"}   # stb-mladsorb site # -> the name used in the literature


def fig_sites():
    st = structure_io.read_fdf(os.path.join(HERE, "si111_2x2.fdf"))
    pmg = structure_io.to_pymatgen(st)
    at = AseAtomsAdaptor.get_atoms(pmg)
    frame = m.SurfaceFrame(at.cell)
    h = frame.height(at.positions)
    layers = m.find_layers(h)
    lay = m.layer_index(layers, len(at))
    _, slots, all_sites, _ = m.find_sites(pmg, frame, ["ontop", "bridge", "hollow"], 2.0, 0.01, lay,
                                          at.get_chemical_symbols())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})

    # side view: projected along a (the [1-10] direction), so the horizontal axis is the in-plane direction
    # perpendicular to it and the bilayers show their zig-zag; three periodic copies along b
    top = layers[0]["h"]
    shades = {1: "0.1", 2: "0.35", 3: "0.55", 4: "0.65", 5: "0.8", 6: "0.85"}
    y = at.positions @ frame.e2
    period = frame.b @ frame.e2
    from ase.neighborlist import neighbor_list
    I, J, D = neighbor_list("ijD", at, 2.5)
    for c in range(-1, 2):
        for i, d in zip(I, D):
            a1.plot([y[i] + c * period, y[i] + d @ frame.e2 + c * period],
                    [h[i] - top, h[i] + d @ frame.n - top], "-", color="0.55", lw=1.4, zorder=1)
        for i in range(len(at)):
            a1.scatter(y[i] + c * period, h[i] - top, s=180, c=shades[lay[i]], edgecolors="k", zorder=3)
        for i in layers[0]["idx"]:
            a1.plot([y[i] + c * period] * 2, [0.15, 1.0], ":", color="#c0392b", lw=2.2, zorder=2)
    for k, L in enumerate(layers, 1):
        a1.text(-period - 1.0, L["h"] - top, f"L{k}", va="center", ha="right", fontsize=9)
    a1.axhline((layers[3]["h"] + layers[4]["h"]) / 2 - top, color="0.3", ls="--", lw=0.8)
    a1.text(-period + 0.2, 1.3, "one dangling bond per surface atom (dotted)", color="#c0392b", fontsize=9)
    a1.text(-period + 0.2, (layers[3]["h"] + layers[4]["h"]) / 2 - top + 0.2,
            "--relax-depth 4: L1-L4 free, bottom bilayer fixed", fontsize=8.5)
    a1.set_xlim(-period - 2.2, 2 * period - 0.5)
    a1.set_ylim(-8.2, 1.9)
    a1.set_xlabel("in-plane position along [11-2] (Ang)")
    a1.set_ylabel("height above layer 1 (Ang)")
    a1.set_title("(a) Si(111) 2x2, three bilayers, cut between bilayers")

    # top view
    out = m._cell_outline(frame)
    lo, hi = out.min(axis=0) - 0.8, out.max(axis=0) + 0.8
    shifts = [frame.xy(i * frame.a + j * frame.b)[0] for i in (-1, 0, 1) for j in (-1, 0, 1)]

    def rep(p):
        pts = np.array([p + s for s in shifts])
        return pts[np.all((pts >= lo) & (pts <= hi), axis=1)]

    sizes = {1: 260, 2: 150, 3: 70, 4: 45}
    for L, lab in ((4, "L4"), (3, "L3"), (2, "L2 (below T4)"), (1, "L1 (surface)")):
        pts = np.vstack([rep(m._in_cell_xy(frame, at.positions[i])) for i in range(len(at)) if lay[i] == L])
        a2.scatter(pts[:, 0], pts[:, 1], s=sizes[L], c=shades[L], edgecolors="k", linewidths=0.4, label=lab,
                   zorder=1 + (5 - L))
    a2.plot(out[:, 0], out[:, 1], "k-", lw=0.8)
    colors = {2: "#1f77b4", 1: "#ff7f0e", 5: "#2ca02c", 3: "#9467bd"}
    for sl in slots:
        if sl["slot"] not in NAMES:
            continue
        pts = np.vstack([rep(m._in_cell_xy(frame, x["coord"])) for x in all_sites if x["slot"] == sl["slot"]])
        a2.scatter(pts[:, 0], pts[:, 1], marker="X", s=110, c=colors[sl["slot"]], edgecolors="k", linewidths=0.6,
                   zorder=9, label=f"{NAMES[sl['slot']]}: stb-mladsorb #{sl['slot']} {sl['type']} ({sl['below']})")
    a2.set_aspect("equal")
    a2.set_xlim(lo[0], hi[0])
    a2.set_ylim(lo[1], hi[1])
    a2.set_xlabel("x (Ang)")
    a2.set_ylabel("y (Ang)")
    a2.set_title("(b) top view: the four high-symmetry adatom sites")
    a2.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "si111_sites.png"), dpi=150)
    plt.close(fig)


def fig_populations():
    de = np.linspace(0, 0.3, 301)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for T, c in ((300, "#1f77b4"), (600, "#ff7f0e"), (1000, "#d62728")):
        p = 1 / (1 + np.exp(-de / (KB * T)))
        a1.plot(1000 * de, 100 * p, color=c, lw=2, label=f"{T} K")
    a1.axhline(90, color="0.6", ls=":")
    a1.set_xlabel("dE = E(site 2) - E(site 1) (meV)")
    a1.set_ylabel("population of site 1 (%)")
    a1.set_title("(a) two sites, equal multiplicity")
    a1.legend()
    T = 300
    for r, c in ((1, "#1f77b4"), (3, "#2ca02c"), (1 / 3, "#9467bd")):
        p = 1 / (1 + r * np.exp(-de / (KB * T)))
        a2.plot(1000 * de, 100 * p, color=c, lw=2, label=f"m2/m1 = {r:.2g}")
    a2.set_xlabel("dE (meV)")
    a2.set_ylabel("population of site 1 (%)")
    a2.set_title("(b) 300 K: the multiplicity ratio shifts the balance")
    a2.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "populations.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_sites()
    fig_populations()
    print("figures written to", FIG)
