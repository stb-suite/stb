#!/usr/bin/env python3
"""Compares the models with SIESTA on MD frames: for every folder reference_results/siesta_check/<frame>/ holding a
finished SIESTA single point (structure.fdf, siesta.out, siesta.FA -- same level of theory as the model's training
data, k-grid divided by the 3x3x3 supercell), evaluates the fine-tuned model of model/ and the MACE-MP medium
foundation model on the same structure (float64) and writes <frame>/forces_compare.json: force RMSE / MAE / max
error, and the energy error per atom after removing each model's constant offset measured on the frame named
REFERENCE_FRAME (energies of different models live on different scales; only differences are meaningful).
results_figures.py prints the table.

    python3 siesta_check.py [--device cuda]        (needs the 'ml' extra)
"""
import argparse
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FRAMES = os.path.join(HERE, "reference_results", "siesta_check")
MODEL = os.path.join(HERE, "model", "mlff_model_compiled.model")
CARD = os.path.join(HERE, "model", "mlff_model_card.json")
REFERENCE_FRAME = "nvt300_final"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    args = ap.parse_args()

    from pymatgen.io.ase import AseAtomsAdaptor
    from stb.core import mace_relax, model_card, siesta_log, structure_io

    model_card.register(MODEL, CARD, echo=False)
    calcs = [("fine-tuned", mace_relax.get_calculator(model=MODEL, device=args.device, dtype="float64")),
             ("MACE-MP medium", mace_relax.get_calculator(model="medium", device=args.device, dtype="float64"))]

    frames = {}
    for d in sorted(glob.glob(os.path.join(FRAMES, "*"))):
        if not os.path.isfile(os.path.join(d, "siesta.FA")):
            continue
        atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(os.path.join(d, "structure.fdf"))))
        f_ref = np.asarray(siesta_log.read_fa_forces(os.path.join(d, "siesta.FA")), float)
        e_ref = siesta_log.get_free_energy(os.path.join(d, "siesta.out"))
        res = {}
        for name, calc in calcs:
            a = atoms.copy()
            a.calc = calc
            f = a.get_forces()
            df = f - f_ref
            res[name] = {"F_rmse_meV": float(np.sqrt((df ** 2).mean()) * 1000), "F_mae_meV": float(np.abs(df).mean() * 1000),
                         "F_max_meV": float(np.abs(df).max() * 1000),
                         "E_minus_siesta_eV_atom": float((a.get_potential_energy() - e_ref) / len(a))}
        frames[os.path.basename(d)] = (d, len(atoms), float(np.sqrt((f_ref ** 2).mean()) * 1000), res)

    if not frames:
        print(f"no finished SIESTA frame under {FRAMES}")
        return
    offsets = {name: frames[REFERENCE_FRAME][3][name]["E_minus_siesta_eV_atom"] for name, _ in calcs} \
        if REFERENCE_FRAME in frames else {}
    for frame, (d, n, rms, res) in frames.items():
        models = []
        for name, r in res.items():
            de = (r["E_minus_siesta_eV_atom"] - offsets[name]) * 1000 if offsets else None
            models.append({"name": name, **r, "dE_meV_atom": de})
            print(f"{frame:18s} {name:15s} F RMSE {r['F_rmse_meV']:7.1f}  MAE {r['F_mae_meV']:7.1f}  max {r['F_max_meV']:6.0f} meV/Ang"
                  + (f"  dE {de:+7.2f} meV/atom" if de is not None else ""))
        json.dump({"frame": frame, "n_atoms": n, "siesta_rms_meV": rms, "energy_offset_reference": REFERENCE_FRAME,
                   "models": models}, open(os.path.join(d, "forces_compare.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
