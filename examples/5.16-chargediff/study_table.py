#!/usr/bin/env python3
"""Helper of example_5.16.sh (README section 5): collects every finished run of reference_results/ and prints

  [A] the three cases: counterpoise and uncorrected interaction energies, BSSE, deformation and binding energies,
      charge moved and the net charge of each fragment;
  [B] the water dimer with five basis sets (same geometry);
  [C] the water dimer with five MeshCutoff values;
  [D] the water dimer along the O...O distance;

with --png it also draws figures/water_studies.png. The numbers come from the analysis JSON that
stb-chargediffAnalysis wrote for every set (--save-data) and from the calc.out of the runs without ghost atoms
(fragmentN_alone/, the fragment in its own basis) and of the relaxed monomers.

    python3 study_table.py [--png] [--ref reference_results]
"""
import argparse
import glob
import json
import os

import numpy as np

from stb.core import siesta_log

HERE = os.path.dirname(os.path.abspath(__file__))
EV_TO_KCAL = 23.0605
MONOMERS = {"water": ["h2o", "h2o"], "nh3bh3": ["nh3", "bh3"], "li_graphene": [None, None]}


def energy(path):
    return siesta_log.get_free_energy(path) if os.path.isfile(path) else None


def load(ref, name):
    """everything about one set: the analysis results, the uncorrected fragment energies, the mesh used"""
    d = os.path.join(ref, name)
    res = json.load(open(glob.glob(os.path.join(d, "analysis", "*results.json"))[0]))
    frags = sorted(k for k in res["energies_ev"] if k != "combined")
    alone = [energy(os.path.join(d, f + "_alone", "calc.out")) for f in frags]
    e = res["energies_ev"]
    out = {"res": res, "frags": frags, "e_cp": res["e_int_cp_ev"], "alone": alone}
    if None not in alone:
        out["e_raw"] = e["combined"] - sum(alone)
        out["bsse"] = out["e_cp"] - out["e_raw"]          # > 0: the ghost orbitals lowered each fragment
    return out


def deformation(ref, name, s):
    """E_def = sum_i [E_i(complex geometry, own basis) - E_i(relaxed alone)] (eV)"""
    if name == "li_graphene":
        flat = energy(os.path.join(ref, "graphene_flat", "calc.out"))
        return None if flat is None or s["alone"][0] is None else s["alone"][0] - flat   # the Li atom: nothing to relax
    mons = MONOMERS[name]
    rel = [energy(os.path.join(ref, "relax", m, "calc.out")) for m in mons]
    if None in rel or None in s["alone"]:
        return None
    return sum(a - r for a, r in zip(s["alone"], rel))


def oo_distance(ref, name):
    from stb.core import structure_io
    st = structure_io.read_fdf(os.path.join(ref, name, "combined", "structure.fdf"))
    lat = np.asarray(st.lattice)
    p = np.array([q for _, q in st.atoms], float)
    if st.coord_format != "cartesian":
        p = p @ lat
    d = (p[3] - p[0]) @ np.linalg.inv(lat)
    d -= np.round(d)
    return float(np.linalg.norm(d @ lat))


def kcal(x):
    return "-" if x is None else f"{x * EV_TO_KCAL:+.2f}"


def table(head, rows):
    w = [max(len(str(r[k])) for r in [head] + rows) for k in range(len(head))]
    print("  " + " | ".join(str(h).ljust(w[k]) for k, h in enumerate(head)))
    print("  " + "-+-".join("-" * x for x in w))
    for r in rows:
        print("  " + " | ".join(str(c).ljust(w[k]) for k, c in enumerate(r)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default=os.path.join(HERE, "reference_results"))
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args()
    ref = a.ref

    print("[A] THE THREE CASES (energies in kcal/mol; 1 eV = 23.06 kcal/mol)")
    rows = []
    for name, label in (("water", "(H2O)2"), ("nh3bh3", "H3N-BH3"), ("li_graphene", "Li / graphene 3x3")):
        if not os.path.isdir(os.path.join(ref, name)):
            continue
        s = load(ref, name)
        edef = deformation(ref, name, s)
        ebind = None if edef is None else s["e_cp"] + edef
        pf = s["res"]["per_fragment_e"]
        rows.append([label, kcal(s["e_cp"]), kcal(s.get("e_raw")), kcal(s.get("bsse")), kcal(edef), kcal(ebind),
                     f"{s['res']['charge_moved_e']:.3f}", " / ".join(f"{pf[f]:+.3f}" for f in s["frags"])])
    table(["system", "E_int CP", "E_int raw", "BSSE", "E_def", "E_bind", "moved (e)", "fragments (e)"], rows)
    print("  E_int raw = E(AB) - sum E(fragment, own basis); BSSE = E_int CP - E_int raw; E_bind = E_int CP + E_def.")

    print("\n[B] WATER DIMER: BASIS SET (same geometry, MeshCutoff 320 Ry)")
    basis = [("water_dzp_s002", "DZP, 0.02 Ry"), ("water", "DZP, 0.01 Ry (base)"), ("water_dzp_s0005", "DZP, 0.005 Ry"),
             ("water_tzp_s001", "TZP, 0.01 Ry"), ("water_tzp_s0005", "TZP, 0.005 Ry")]
    rows, bstud = [], []
    for name, label in basis:
        if not os.path.isdir(os.path.join(ref, name)):
            continue
        s = load(ref, name)
        pf = s["res"]["per_fragment_e"]
        rows.append([label, kcal(s["e_cp"]), kcal(s.get("e_raw")), kcal(s.get("bsse")),
                     f"{s['res']['charge_moved_e']:.4f}", f"{pf['fragment2']:+.4f}"])
        bstud.append((label, s["e_cp"] * EV_TO_KCAL, s.get("e_raw", np.nan) * EV_TO_KCAL, pf["fragment2"]))
    table(["basis", "E_int CP", "E_int raw", "BSSE", "moved (e)", "donor (e)"], rows)

    print("\n[C] WATER DIMER: MESH (DZP, 0.01 Ry)")
    rows = []
    for m in (150, 200, 320, 450, 600):
        name = "water" if m == 320 else f"water_mesh{m}"
        if not os.path.isdir(os.path.join(ref, name)):
            continue
        s = load(ref, name)
        r = s["res"]
        used = [c["detail"] for c in r["checks"] if c["name"].startswith("Same mesh")][0]
        rows.append([f"{m}", used, kcal(s["e_cp"]), f"{r['integral_e']:+.1e}", f"{r['charge_moved_e']:.4f}",
                     f"{r['per_fragment_e']['fragment2']:+.4f}", f"{r['maximum']['value']:.4f}"])
    table(["MeshCutoff (Ry)", "used", "E_int CP", "integral (e)", "moved (e)", "donor (e)", "max (e/Ang^3)"], rows)

    print("\n[D] WATER DIMER: O...O DISTANCE (donor moved rigidly; DZP, 0.01 Ry, 320 Ry)")
    rows, scan = [], []
    for d in sorted(glob.glob(os.path.join(ref, "water_r*"))) + [os.path.join(ref, "water")]:
        name = os.path.basename(d)
        s = load(ref, name)
        r = oo_distance(ref, name)
        scan.append((r, s["e_cp"] * EV_TO_KCAL, s.get("e_raw", np.nan) * EV_TO_KCAL,
                     s["res"]["per_fragment_e"]["fragment2"], s["res"]["charge_moved_e"]))
    scan.sort()
    for r, ecp, eraw, q, mv in scan:
        rows.append([f"{r:.3f}", f"{ecp:+.2f}", f"{eraw:+.2f}", f"{eraw - ecp:+.2f}", f"{q:+.4f}", f"{mv:.4f}"])
    table(["d(O...O) (Ang)", "E_int CP", "E_int raw", "-BSSE", "donor (e)", "moved (e)"], rows)
    far = [(r, q) for r, _, _, q, _ in scan if r >= 2.9 and q > 0]
    if len(far) >= 3:
        k, lnq0 = np.polyfit([r for r, _ in far], np.log([q for _, q in far]), 1)
        print(f"  Net charge of the donor beyond 2.9 Ang: q(R) ~ exp(-R / {-1 / k:.2f} Ang) (fit of {len(far)} points)")
        i = int(np.argmin([e for _, e, _, _, _ in scan]))
        print(f"  Minimum of E_int CP along the scan: {scan[i][1]:+.2f} kcal/mol at {scan[i][0]:.2f} Ang "
              "(rigid scan: the monomers keep the geometry of the relaxed dimer)")

    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        r = np.array([x[0] for x in scan])
        ax[0].plot(r, [x[1] for x in scan], "o-", color="#0072B2", lw=2, label="counterpoise (ghost basis)")
        ax[0].plot(r, [x[2] for x in scan], "s--", color="#D55E00", lw=2, label="uncorrected (own basis)")
        ax[0].axhline(0, color="k", lw=0.6)
        ax[0].axhline(-5.02, color="#009E73", ls=":", lw=1.5, label="CCSD(T)/CBS $D_e$ = 5.02 (Klopper 2000)")
        ax[0].axvline(2.912, color="#009E73", ls=":", lw=1.0)
        ax[0].set_xlabel("d(O···O) (Å)")
        ax[0].set_ylabel("E_int (kcal/mol)")
        ax[0].set_title("Interaction energy along a rigid scan")
        ax[0].legend(fontsize=8)
        ax[1].semilogy(r, [x[3] for x in scan], "o-", color="#0072B2", lw=2, label="net charge of the donor")
        ax[1].semilogy(r, [x[4] for x in scan], "s--", color="#CC79A7", lw=2, label="charge moved ½∫|Δρ|")
        ax[1].set_xlabel("d(O···O) (Å)")
        ax[1].set_ylabel("charge (e)")
        ax[1].set_title("Charge transfer vs polarization")
        ax[1].legend(fontsize=8)
        x = np.arange(len(bstud))
        ax[2].bar(x - 0.2, [b[1] for b in bstud], 0.4, color="#0072B2", label="counterpoise")
        ax[2].bar(x + 0.2, [b[2] for b in bstud], 0.4, color="#D55E00", label="uncorrected")
        ax[2].axhline(-5.02, color="#009E73", ls=":", lw=1.5)
        ax[2].set_xticks(x)
        ax[2].set_xticklabels([b[0].replace(" (base)", "\n(base)").replace(", ", "\n") for b in bstud], fontsize=8)
        ax[2].set_ylabel("E_int at the relaxed geometry (kcal/mol)")
        ax[2].set_title("Basis set: BSSE shrinks, CP converges from above")
        ax[2].legend(fontsize=8)
        fig.tight_layout()
        os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
        out = os.path.join(HERE, "figures", "water_studies.png")
        fig.savefig(out, dpi=140)
        print(f"\n  figure: {os.path.relpath(out, HERE)}")


if __name__ == "__main__":
    main()
