#!/usr/bin/env python3
"""Summary of the 7.6 reference set (README section 6): reads reference_results/<case>/mladsorb_minima.dat (and the
binding-curve data), prints one table -- E_form and the energy above T4 of the four Si(111) adatom sites in every case --
and draws:

    figures/relax_depth.png     E_form of each site against --relax-depth (0, 1, 4, 7 Ang)
    figures/models.png          E_form of each site with every model, and SIESTA's when reference_results/siesta/ exists
    figures/binding_curves.png  E_ads against the height of the adatom held fixed, per site

    python3 results_figures.py       (numpy + matplotlib, seconds; no model)
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference_results")
FIG = os.path.join(HERE, "figures")
SITES = ["T4", "H3", "B2", "T1"]
SLOT = {1: "T1", 2: "T4", 3: "H3", 5: "B2"}            # stb-mladsorb site numbers on these Si(111) slabs
COLORS = {"T4": "#1f77b4", "H3": "#9467bd", "B2": "#2ca02c", "T1": "#ff7f0e"}
CASES = [("main", "round 3, top 4 A relaxed (main)"), ("depth_0", "round 3, substrate fixed"),
         ("depth_1", "round 3, top 1 A relaxed"), ("depth_7", "round 3, top 7 A relaxed"),
         ("start_heights", "round 3, 4 start heights"), ("thick_4bl", "round 3, 4 bilayers"),
         ("cell_3x3", "round 3, 3x3 cell (1/9 ML)"), ("round2", "round 2"),
         ("foundation_small", "MACE-MP-0 small"), ("foundation_medium", "MACE-MP-0 medium"),
         ("foundation_large", "MACE-MP-0 large")]


def minima(case):
    """{site: (E_form, E_ads, height, CN)} of the lowest minimum at each named site of one case."""
    p = os.path.join(REF, case, "mladsorb_minima.dat")
    if not os.path.isfile(p):
        return None
    out = {}
    for line in open(p):
        if line.startswith("#"):
            continue
        c = line.split()
        name = SLOT.get(int(c[2]))
        if name and name not in out:
            out[name] = (float(c[6]), float(c[4]), float(c[9]), int(c[10]))
    return out


def siesta():
    """{site: E_form} of the SIESTA relaxations that stayed at their site (reference_results/siesta/siesta_eform.dat,
    written by siesta_compare.py); a relaxation that ended at another site is no minimum of SIESTA there."""
    p = os.path.join(REF, "siesta", "siesta_eform.dat")
    if not os.path.isfile(p):
        return {}
    out = {}
    for line in open(p):
        if line.startswith("#"):
            continue
        c = line.split()
        if c[3] == c[0]:
            out[c[0]] = float(c[1])
    return out


def main():
    data = {c: minima(c) for c, _ in CASES}
    si = siesta()
    head = f"{'case':<34}" + "".join(f"{s:>9}" for s in SITES) + "   |" + "".join(f"{s:>8}" for s in SITES[1:])
    print("E_form = E(slab+Si) - E(slab) - mu_Si (eV), lowest minimum at each site   |   energy above T4 (eV)")
    print(head)
    print("-" * len(head))
    rows = [(lab, data[c]) for c, lab in CASES if data[c] is not None]
    if si:
        rows.append(("SIESTA, relaxed (section 7)", {k: (v, None, None, None) for k, v in si.items()}))
    for lab, d in rows:
        ef = [d[s][0] if s in d else np.nan for s in SITES]
        t4 = d.get("T4", (np.nan,))[0]
        print(f"{lab:<34}" + "".join(f"{x:>9.3f}" if np.isfinite(x) else f"{'-':>9}" for x in ef) + "   |"
              + "".join(f"{x - t4:>8.3f}" if np.isfinite(x) else f"{'-':>8}" for x in ef[1:]))
    print("'-': no run ended at that site (not a minimum of that model there; for SIESTA: its relaxation from the "
          "model's geometry ended at another site).")
    os.makedirs(FIG, exist_ok=True)

    # relaxation depth
    depths = [(0, "depth_0"), (1, "depth_1"), (4, "main"), (7, "depth_7")]
    if all(data[c] for _, c in depths):
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
        for s in SITES:
            x = [d for d, c in depths if s in data[c]]
            y = [data[c][s][0] for d, c in depths if s in data[c]]
            ax.plot(x, y, "o-", color=COLORS[s], label=s, ms=7)
        if si:
            for s in SITES:
                if s in si:
                    ax.axhline(si[s], color=COLORS[s], ls=":", lw=1)
            ax.plot([], [], "k:", label="SIESTA (top 4 A relaxed)")
        ax.set_xlabel("--relax-depth (Ang): substrate atoms this close to the top relax")
        ax.set_ylabel("E_form (eV)")
        ax.set_title("Round-3 model: how much of the slab relaxes")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "relax_depth.png"), dpi=150)
        plt.close(fig)

    # models
    mods = [("main", "round 3"), ("round2", "round 2"), ("foundation_small", "MP-0 small"),
            ("foundation_medium", "MP-0 medium"), ("foundation_large", "MP-0 large")]
    mods = [(c, lab) for c, lab in mods if data[c]]
    if mods:
        fig, ax = plt.subplots(figsize=(9, 4.4))
        series = [(lab, data[c]) for c, lab in mods]
        if si:
            series = [("SIESTA", {k: (v,) for k, v in si.items()})] + series
        w = 0.8 / len(series)
        x = np.arange(len(SITES))
        cols = ["k", "#1f77b4", "#17becf", "#ff7f0e", "#d62728", "#9467bd"][-len(series):] if not si else \
            ["k", "#1f77b4", "#17becf", "#ff7f0e", "#d62728", "#9467bd"]
        for k, (lab, d) in enumerate(series):
            ax.bar(x + (k - (len(series) - 1) / 2) * w, [d[s][0] if s in d else np.nan for s in SITES], w,
                   color=cols[k], label=lab, edgecolor="k", lw=0.4)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(SITES)
        ax.set_ylabel("E_form (eV)")
        ax.set_title("Si adatom on Si(111) 2x2, top 4 Ang relaxed: E_form per site (missing bar: not a minimum)")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "models.png"), dpi=150)
        plt.close(fig)

    # binding curves
    p = os.path.join(REF, "binding_curve", "mladsorb_binding_curve.dat")
    if os.path.isfile(p):
        blocks, cur, name = {}, [], None
        for line in open(p):
            if line.startswith("# Si #"):
                if name and cur:
                    blocks[name] = np.array(cur)
                slot = int(line.split("#")[2].split()[0])
                name, cur = SLOT.get(slot), []
            elif line.strip() and not line.startswith("#"):
                c = line.split()
                if len(c) < 5 or c[4] == "0":          # column 5: 1 = too close to an atom, left out by the tool
                    cur.append([float(c[0]), float(c[1])])
        if name and cur:
            blocks[name] = np.array(cur)
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
        for s in SITES:
            if s in blocks:
                b = blocks[s]
                ax.plot(b[:, 0], b[:, 1], "o-", color=COLORS[s], label=s, ms=5)
        ax.set_xlabel("adatom height above layer 1, held fixed (Ang)")
        ax.set_ylabel("E_ads (eV)")
        ax.set_title("Round-3 model: binding curves at fixed height")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "binding_curves.png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
