#!/usr/bin/env python3
"""The SIESTA phonons of the 7.1 example against the fine-tuned model and the foundation model, on the SAME displaced
supercell, plus two experiments on the SIESTA forces (README section 5.9). Needs the 'ml' extra; runs on the CPU in
about a minute.

    python3 siesta_compare.py [--device cpu|cuda]

Input (reference_results/siesta_phonons/, copied from the 7.1 silicon run):
  phonopy_disp.yaml   stb-phononsCreate's phonopy file: 2x2x2 conventional supercell (64 atoms), ONE symmetry-reduced
                      displacement of 0.01 Ang (lengths in Bohr, phonopy's SIESTA convention)
  FORCE_SETS          SIESTA's forces on that supercell (eV/Ang), the undisplaced supercell's residual subtracted
  FORCE_SETS.raw      the same before the subtraction
  calc.fdf, config_extra.fdf, structure.fdf   the SIESTA input of the displaced supercell

What it computes, into reference_results/siesta_compare.json:
  1. frequencies at Gamma, X, L and on a 12^3 mesh for SIESTA, the fine-tuned model and MACE-MP-0 medium, each with
     the same displacement and the same phonopy settings, against experiment (Nilsson & Nelin 1972);
  2. the force on every atom of the displaced supercell, SIESTA vs model, by neighbour shell;
  3. the SYMMETRY SCATTER of each force set: the displacement along x leaves the site symmetry of the displaced atom
     partly intact, so symmetry-equivalent atoms must feel symmetry-equivalent forces. Averaging the force set over
     those operations and subtracting it measures the noise of the force set by itself, with no model involved;
  4. a noise experiment: Gaussian noise of a given size added to the model's clean forces, 20 random draws per size,
     at the 0.01 Ang displacement SIESTA used and at 0.03 Ang -- how much noise it takes to produce SIESTA's
     imaginary acoustic modes, and how a larger displacement suppresses it.
"""
import argparse
import contextlib
import io
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
SIESTA = os.path.join(HERE, "reference_results", "siesta_phonons")
MODEL_DIR = os.path.join(HERE, "..", "7.2-mlmd", "model")
OUT = os.path.join(HERE, "reference_results", "siesta_compare.json")
BOHR = 0.529177210903
MESH = [12, 12, 12]
# Si phonons at 296 K, inelastic neutron scattering: G. Nilsson, G. Nelin, Phys. Rev. B 6, 3777 (1972), THz
EXPERIMENT = {"G": [0, 0, 0, 15.53, 15.53, 15.53], "X": [4.49, 4.49, 12.32, 12.32, 13.90, 13.90],
              "L": [3.43, 3.43, 11.35, 12.60, 14.68, 14.68]}
QPOINTS = {"G": [0, 0, 0], "X": [0.5, 0.0, 0.5], "L": [0.5, 0.5, 0.5]}     # primitive-cell reduced coordinates


def quiet():
    return contextlib.ExitStack() if os.environ.get("STB_VERBOSE") else _Quiet()


class _Quiet(contextlib.ExitStack):
    def __enter__(self):
        super().__enter__()
        self.enter_context(contextlib.redirect_stdout(io.StringIO()))
        self.enter_context(contextlib.redirect_stderr(io.StringIO()))
        return self


def analyse(phonon):
    """Frequencies at G/X/L, mesh statistics. Acoustic modes at Gamma excluded from the minimum."""
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()
    pts = {k: np.sort(phonon.get_frequencies(q)).round(3).tolist() for k, q in QPOINTS.items()}
    phonon.run_mesh(MESH)
    f = np.array(phonon.mesh.frequencies)
    g = np.all(np.abs(phonon.mesh.qpoints) < 1e-9, axis=1)
    f_off = f.copy()
    f_off[g, :3] = np.nan
    return {"points": pts, "mesh": f, "min": float(np.nanmin(f_off)),
            "imaginary_q": int(np.sum(np.nanmin(f_off, axis=1) < -0.1)), "n_q": int(len(f))}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    args = ap.parse_args()

    import phonopy
    import spglib
    from ase import Atoms
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    sys.path.insert(0, os.path.join(HERE, "..", "..", "stb-suite", "src"))
    from stb.core import model_card, mace_relax

    ref = phonopy.load(os.path.join(SIESTA, "phonopy_disp.yaml"), force_sets_filename=os.path.join(SIESTA, "FORCE_SETS"),
                       symprec=0.01, log_level=0)
    raw = phonopy.load(os.path.join(SIESTA, "phonopy_disp.yaml"),
                       force_sets_filename=os.path.join(SIESTA, "FORCE_SETS.raw"), symprec=0.01, log_level=0)
    disp = ref.dataset["first_atoms"][0]
    i0 = disp["number"]
    u_bohr = np.array(disp["displacement"])
    F_siesta = np.array(disp["forces"])                       # phonopy's SIESTA interface: eV/Ang
    F_siesta_raw = np.array(raw.dataset["first_atoms"][0]["forces"])

    # the same structures in Angstrom, for the ML models
    uc = PhonopyAtoms(symbols=ref.unitcell.symbols, cell=ref.unitcell.cell * BOHR,
                      scaled_positions=ref.unitcell.scaled_positions)
    sc0 = ref.supercell
    sc_pos, sc_cell = sc0.positions * BOHR, sc0.cell * BOHR

    def ml_phonon(forces, u_ang):
        p = Phonopy(uc, supercell_matrix=ref.supercell_matrix, primitive_matrix=ref.primitive_matrix, symprec=0.01)
        p.dataset = {"natom": len(sc0), "first_atoms": [{"number": i0, "displacement": u_ang, "forces": forces}]}
        return p

    with quiet():
        model_card.setup_from_args(os.path.join(MODEL_DIR, "mlff_model_compiled.model"),
                                   os.path.join(MODEL_DIR, "mlff_model_card.json"), "--custom-model")
        calcs = {"fine-tuned": mace_relax.get_calculator(os.path.join(MODEL_DIR, "mlff_model_compiled.model"),
                                                         device=args.device),
                 "MACE-MP-0 medium": mace_relax.get_calculator("medium", device=args.device)}

    def forces(calc, pos):
        a = Atoms(numbers=sc0.numbers, positions=pos, cell=sc_cell, pbc=True)
        a.calc = calc
        with quiet():
            return a.get_forces()

    u_dir = u_bohr / np.linalg.norm(u_bohr)
    out = {"experiment": EXPERIMENT, "mesh": MESH, "displacement_ang": float(np.linalg.norm(u_bohr) * BOHR)}
    results, model_forces = {}, {}
    results["SIESTA"] = analyse(ref)
    results["SIESTA (raw forces)"] = analyse(raw)
    for name, calc in calcs.items():
        f_perfect = forces(calc, sc_pos)
        fsets = {}
        for u in (0.01, 0.03):
            pos = sc_pos.copy()
            pos[i0] += u * u_dir
            fsets[u] = forces(calc, pos) - f_perfect
        model_forces[name] = fsets
        results[name] = analyse(ml_phonon(fsets[0.01], 0.01 * u_dir))

    # forces by neighbour shell
    dv = sc_pos - sc_pos[i0]
    fr = dv @ np.linalg.inv(sc_cell)
    fr -= np.round(fr)
    r = np.linalg.norm(fr @ sc_cell, axis=1)
    shells = sorted(set(np.round(r, 2)))
    Fm = model_forces["fine-tuned"][0.01]
    rows = []
    for s in shells:
        idx = np.abs(r - s) < 0.01
        rows.append({"r": float(s), "n": int(idx.sum()),
                     "siesta_rms": float(1000 * np.sqrt(np.mean(F_siesta[idx] ** 2))),
                     "model_rms": float(1000 * np.sqrt(np.mean(Fm[idx] ** 2))),
                     "diff_rms": float(1000 * np.sqrt(np.mean((F_siesta[idx] - Fm[idx]) ** 2)))})
    out["force_shells_meV_per_ang"] = rows

    # symmetry scatter: operations of the supercell that keep the displaced atom and its displacement in place
    dataset = spglib.get_symmetry((sc_cell, sc0.scaled_positions, sc0.numbers), symprec=1e-3)
    rots_c = [sc_cell.T @ R @ np.linalg.inv(sc_cell.T) for R in dataset["rotations"]]
    ops = []
    for R, Rc, t in zip(dataset["rotations"], rots_c, dataset["translations"]):
        if not np.allclose(Rc @ u_dir, u_dir, atol=1e-6):
            continue
        new = (sc0.scaled_positions @ R.T + t) % 1.0
        perm = []
        for p in new:
            dd = sc0.scaled_positions - p
            dd -= np.round(dd)
            perm.append(int(np.argmin(np.linalg.norm(dd @ sc_cell, axis=1))))
        if perm[i0] == i0:
            ops.append((Rc, np.array(perm)))

    def symmetrize(F):
        acc = np.zeros_like(F)
        for Rc, perm in ops:
            acc[perm] += F @ Rc.T
        return acc / len(ops)

    def scatter(F, beyond=3.0):
        d = F - symmetrize(F)
        return float(1000 * np.sqrt(np.mean(d[r > beyond] ** 2)))

    out["symmetry_scatter_meV_per_ang"] = {
        "n_site_operations": len(ops),
        "SIESTA": scatter(F_siesta), "SIESTA (raw)": scatter(F_siesta_raw),
        "fine-tuned": scatter(Fm), "MACE-MP-0 medium": scatter(model_forces["MACE-MP-0 medium"][0.01]),
        "SIESTA_true_force_rms_beyond_3A": float(1000 * np.sqrt(np.mean(F_siesta[r > 3] ** 2))),
        "model_force_rms_beyond_3A": float(1000 * np.sqrt(np.mean(Fm[r > 3] ** 2)))}

    # the outlier: SIESTA's own .FA files keep the net force that phonopy removes from FORCE_SETS (it subtracts the
    # mean force per atom, the "drift"). The atom whose force departs most from its symmetry partners is then
    # replaced by the value its partners imply (symmetrization iterated so that its own force drops out).
    def read_fa(path):
        lines = open(path).read().split("\n")
        n = int(lines[0])
        return np.array([[float(x) for x in s.split()[1:4]] for s in lines[1:n + 1]])

    fa_disp, fa_perf = read_fa(os.path.join(SIESTA, "disp-001.FA")), read_fa(os.path.join(SIESTA, "perfect_supercell.FA"))
    dev = np.linalg.norm(F_siesta - symmetrize(F_siesta), axis=1)
    dev[r < 3] = 0.0
    k_out = int(np.argmax(dev))
    partners = [j for j in range(len(r)) if abs(r[j] - r[k_out]) < 0.01 and j != k_out]
    repaired = F_siesta.copy()
    for _ in range(30):
        repaired[k_out] = symmetrize(repaired)[k_out]
    repaired -= repaired.mean(axis=0)
    rep = Phonopy(ref.unitcell, supercell_matrix=ref.supercell_matrix, primitive_matrix=ref.primitive_matrix,
                  symprec=0.01, calculator="siesta")
    rep.dataset = {"natom": len(sc0), "first_atoms": [{"number": i0, "displacement": u_bohr, "forces": repaired}]}
    with quiet():
        results["SIESTA, one atom repaired"] = analyse(rep)
    out["outlier"] = {
        "net_force_meV_per_ang": {"displaced": (1000 * fa_disp.sum(axis=0)).round(2).tolist(),
                                  "perfect": (1000 * fa_perf.sum(axis=0)).round(2).tolist()},
        "atom": k_out + 1, "distance_ang": float(r[k_out]),
        "force_meV_per_ang": (1000 * F_siesta[k_out]).round(2).tolist(),
        "partners_median_meV_per_ang": float(1000 * np.median(np.linalg.norm(F_siesta[partners], axis=1))),
        "deviation_from_symmetry_meV_per_ang": float(1000 * dev[k_out]),
        "next_largest_deviation_meV_per_ang": float(1000 * np.sort(dev)[-2])}

    ref_mesh = results["SIESTA"]["mesh"]
    for name, res in results.items():
        d = res["mesh"] - ref_mesh
        res["rmse_vs_siesta"] = {"all": float(np.sqrt(np.mean(d ** 2))), "acoustic": float(np.sqrt(np.mean(d[:, :3] ** 2))),
                               "optical": float(np.sqrt(np.mean(d[:, 3:] ** 2))), "mean_shift": float(d.mean())}
        res.pop("mesh")
    out["phonons"] = results

    # noise experiment
    rng = np.random.default_rng(0)
    noise_rows = []
    for u in (0.01, 0.03):
        clean = model_forces["fine-tuned"][u]
        for sigma in (0.0, 0.5, 1.0, 2.0, 3.0, 5.0):
            mins, xta, imag = [], [], []
            for _ in range(1 if sigma == 0 else 20):
                n = rng.normal(0, sigma / 1000, clean.shape)
                n -= n.mean(axis=0)
                with quiet():
                    a = analyse(ml_phonon(clean + n, u * u_dir))
                mins.append(a["min"])
                imag.append(a["imaginary_q"])
                xta.append(a["points"]["X"][0])
            noise_rows.append({"u_ang": u, "sigma_meV": sigma, "draws": len(mins),
                               "share_with_imaginary": float(np.mean(np.array(imag) > 0)),
                               "min_f_median": float(np.median(mins)), "min_f_worst": float(np.min(mins)),
                               "X_TA_mean": float(np.mean(xta)), "X_TA_std": float(np.std(xta))})
    out["noise_experiment"] = noise_rows

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)

    # a readable summary
    print("Frequencies (THz), sorted")
    for k in QPOINTS:
        print(f"  {k}  experiment " + " ".join(f"{x:6.2f}" for x in EXPERIMENT[k]))
        for name, res in results.items():
            print(f"     {name:26s} " + " ".join(f"{x:6.2f}" for x in res["points"][k]))
    print("\n12^3 mesh: min frequency, q-points with an imaginary mode, RMSE vs SIESTA (acoustic / optical)")
    for name, res in results.items():
        e = res["rmse_vs_siesta"]
        print(f"  {name:26s} min {res['min']:6.2f}  imaginary {res['imaginary_q']:3d}/{res['n_q']}  "
              f"RMSE {e['all']:.2f} ({e['acoustic']:.2f} / {e['optical']:.2f})")
    print("\nForces by neighbour shell (meV/Ang RMS per component): r, n, SIESTA, model, difference")
    for row in rows[:8]:
        print(f"  {row['r']:6.3f} {row['n']:3d} {row['siesta_rms']:8.2f} {row['model_rms']:8.2f} {row['diff_rms']:8.2f}")
    print("\nSymmetry scatter beyond 3 Ang (meV/Ang):", json.dumps(out["symmetry_scatter_meV_per_ang"]))
    o = out["outlier"]
    print(f"\nNet force (meV/Ang): displaced {o['net_force_meV_per_ang']['displaced']}, perfect {o['net_force_meV_per_ang']['perfect']}")
    print(f"Outlier: atom {o['atom']} at {o['distance_ang']:.2f} Ang, force {o['force_meV_per_ang']} meV/Ang, its partners "
          f"{o['partners_median_meV_per_ang']:.2f}; deviation from symmetry {o['deviation_from_symmetry_meV_per_ang']:.1f} "
          f"(next largest {o['next_largest_deviation_meV_per_ang']:.1f})")
    print("\nNoise experiment: u, sigma, share with imaginary modes, median/worst min f, X-TA mean +- std")
    for row in noise_rows:
        print(f"  {row['u_ang']:.2f} {row['sigma_meV']:4.1f}  {row['share_with_imaginary']:4.2f}  "
              f"{row['min_f_median']:6.2f} {row['min_f_worst']:6.2f}  {row['X_TA_mean']:5.2f} +- {row['X_TA_std']:.2f}")
    print(f"\n[OK] {OUT}")


if __name__ == "__main__":
    main()
