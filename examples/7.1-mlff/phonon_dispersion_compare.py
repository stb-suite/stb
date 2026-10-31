#!/usr/bin/env python3
"""Graphene phonon dispersion along Gamma-M-K-Gamma: SIESTA vs the fine-tuned committees of
an MLFF workflow vs the raw MACE-MP model -- all on the SAME displaced supercell of a finished
stb-phononsCreate folder, each with its OWN forces on the undisplaced supercell subtracted
(perfect_supercell/), symprec 0.01 like the suite.

Helper of README section 8.5 (figures/phonon_dispersion.png). Needs the trained models, which
the example does not ship -- run it on your own run folder:

    python phonon_dispersion_compare.py --phonon-dir phonons/phonon_runs \\
        --workflow mlff_runs/mlff_workflow.json --rounds 1 2 \\
        [--extra-ref OTHER_PHONON_DIR "label"] -o phonon_dispersion.png
"""
import argparse
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
B2A = 0.529177210903          # phonopy's SIESTA interface works in Bohr
SYMPREC = 0.01
THZ_TO_CM = 33.35641


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phonon-dir", required=True, help="finished stb-phononsCreate folder (phonon_runs/)")
    ap.add_argument("--workflow", default="mlff_runs/mlff_workflow.json")
    ap.add_argument("--rounds", type=int, nargs="+", default=[1])
    ap.add_argument("--foundation", default="medium", help="MACE-MP size to compare ('' = none)")
    ap.add_argument("--extra-ref", nargs=2, metavar=("DIR", "LABEL"), default=None,
                    help="another finished SIESTA phonon folder of the same supercell (e.g. another k-grid)")
    ap.add_argument("-o", "--output", default="phonon_dispersion.png")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import phonopy
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    from stb.core import mlff_workflow, siesta_log, mace_relax
    from stb.core.cli import capture_library_noise
    from mace.calculators import MACECalculator
    noise = []

    ph0 = phonopy.load(os.path.join(a.phonon_dir, "phonopy_disp.yaml"), symprec=SYMPREC,
                       produce_fc=False, log_level=0)
    uc = PhonopyAtoms(symbols=ph0.unitcell.symbols, cell=np.asarray(ph0.unitcell.cell) * B2A,
                      scaled_positions=ph0.unitcell.scaled_positions)
    disps = [{"number": d["number"], "displacement": np.array(d["displacement"]) * B2A}
             for d in ph0.dataset["first_atoms"]]

    def build(forces):
        ph = Phonopy(uc, supercell_matrix=ph0.supercell_matrix, primitive_matrix=ph0.primitive_matrix,
                     symprec=SYMPREC)
        ph.dataset = {"natom": len(ph0.supercell.symbols), "first_atoms": [dict(d) for d in disps]}
        ph.forces = np.array(forces)
        ph.produce_force_constants()
        ph.symmetrize_force_constants(show_drift=False)
        return ph

    def siesta_model(directory):
        f0 = siesta_log.read_fa_forces(os.path.join(directory, "perfect_supercell", "siesta.FA"))
        fd = [siesta_log.read_fa_forces(os.path.join(directory, f"disp-{k + 1:03d}", "siesta.FA"))
              for k in range(len(disps))]
        f0 = f0 - f0.mean(axis=0)                 # phonopy removes each set's drift
        return build([f - f0 for f in fd])

    models = {"SIESTA": siesta_model(a.phonon_dir)}
    if a.extra_ref:
        models[a.extra_ref[1]] = siesta_model(a.extra_ref[0])

    template = build([np.zeros((len(ph0.supercell.symbols), 3)) for _ in disps])
    as_ase = lambda c: Atoms(symbols=c.symbols, cell=c.cell, scaled_positions=c.scaled_positions, pbc=True)
    perfect = as_ase(template.supercell)
    displaced = [as_ase(sc) for sc in template.supercells_with_displacements]

    def mean_forces(calcs, atoms):
        out = []
        for c in calcs:
            x = atoms.copy()
            x.calc = c
            out.append(x.get_forces())
        return np.mean(out, axis=0)

    state = mlff_workflow.load(a.workflow)
    base = mlff_workflow.base_dir(a.workflow)
    sets = []
    for k in a.rounds:
        paths = mlff_workflow.committee_models(mlff_workflow.get_round(state, k), base)
        with capture_library_noise(noise, f"round {k}"):
            sets.append((f"fine-tuned, round {k} ({len(paths)} models)",
                         [MACECalculator(model_paths=[m], device="cpu", default_dtype="float64") for m in paths]))
    if a.foundation:
        with capture_library_noise(noise, "foundation"):
            sets.append((f"MACE-MP {a.foundation} (no fine-tuning)",
                         [mace_relax.get_calculator(model=a.foundation, device="cpu")]))
    for name, calcs in sets:
        with capture_library_noise(noise, name):
            f0 = mean_forces(calcs, perfect)
            models[name] = build([mean_forces(calcs, x) - f0 for x in displaced])

    # Gamma-M-K-Gamma in the reduced coordinates of the phonopy primitive cell
    prim = np.asarray(models["SIESTA"].primitive.cell)
    rec = np.linalg.inv(prim).T
    b1, b2 = rec[0], rec[1]
    cands = [(b1 + b2) / 3, (2 * b1 + b2) / 3, (b1 + 2 * b2) / 3, (2 * b1 - b2) / 3, (b1 - b2) / 3]
    k_cart = min(cands, key=lambda q: abs(np.linalg.norm(q) - np.linalg.norm(b1) / np.sqrt(3)))
    to_frac = lambda c: np.linalg.solve(rec.T, c)
    nodes, labels, npts = [np.zeros(3), to_frac(b1 / 2), to_frac(k_cart), np.zeros(3)], \
        ["$\\Gamma$", "M", "K", "$\\Gamma$"], 120
    qs, dist, ticks, d0 = [], [], [0.0], 0.0
    for p, q in zip(nodes[:-1], nodes[1:]):
        step = np.linalg.norm((q - p) @ rec) / (npts - 1)
        for j, t in enumerate(np.linspace(0, 1, npts)):
            if qs and j == 0:
                continue
            qs.append(p + (q - p) * t)
            dist.append(d0 + j * step)
        d0 = dist[-1]
        ticks.append(d0)
    qs, dist = np.array(qs), np.array(dist)
    curves = {}
    for name, ph in models.items():
        ph.run_qpoints(qs)
        curves[name] = np.sort(ph.get_qpoints_dict()["frequencies"], axis=1)

    palette = {"SIESTA": dict(color="k", lw=2.4, ls="-")}
    others = [dict(color="tab:blue", lw=1.6, ls="-"), dict(color="tab:green", lw=1.3, ls="--"),
              dict(color="tab:purple", lw=1.3, ls="--")]
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    i_ml = 0
    for name, f in curves.items():
        if name == "SIESTA":
            st = palette["SIESTA"]
        elif a.extra_ref and name == a.extra_ref[1]:
            st = dict(color="0.55", lw=1.2, ls="-.")
        elif name.startswith("MACE-MP"):
            st = dict(color="tab:orange", lw=1.5, ls=":")
        else:
            st = others[i_ml % len(others)]
            i_ml += 1
        for j in range(f.shape[1]):
            ax.plot(dist, f[:, j], label=name if j == 0 else None, **st)
    for t in ticks:
        ax.axvline(t, color="0.75", lw=0.6)
    ax.axhline(0, color="0.75", lw=0.6)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlim(dist[0], dist[-1])
    ax.set_ylim(bottom=min(-1.0, min(f.min() for f in curves.values()) - 1))
    ax.set_ylabel("Frequency (THz)")
    ax2 = ax.twinx()
    ax2.set_ylim(np.array(ax.get_ylim()) * THZ_TO_CM)
    ax2.set_ylabel("cm$^{-1}$")
    n_sc = int(round(abs(np.linalg.det(ph0.supercell_matrix)) ** 0.5))
    ax.set_title(f"Graphene phonons: SIESTA vs fine-tuned MACE vs MACE-MP\n"
                 f"{n_sc}x{n_sc}x1 supercell, the same displacement, residual forces subtracted")
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(a.output, dpi=150)

    ref = curves["SIESTA"]
    print(f"{'model':42} {'RMSE':>6} {'max|d|':>7}  Gamma E2g   K top   M top  (THz)")
    for name, f in curves.items():
        d = f - ref
        print(f"{name:42} {np.sqrt(np.mean(d ** 2)):6.2f} {np.abs(d).max():7.2f}  "
              f"{f[0, -1]:9.2f} {f[2 * (npts - 1), -1]:7.2f} {f[npts - 1, -1]:7.2f}")
    print(f"figure: {a.output}")


if __name__ == "__main__":
    main()
