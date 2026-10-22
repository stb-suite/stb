#!/usr/bin/env python3
"""Gamma-point E2g frequency of graphene by frozen phonon, from finished SIESTA runs.

Helper of example_4.17.sh (Section 8.1 of the README): the two atoms of the primitive
cell are displaced by +u and -u along x (the E2g optical mode at Gamma), SIESTA gives the
restoring forces, and for this two-atom mode

    omega^2 = -(F1x - F2x) / (m * (u1 - u2)),   u1 - u2 = 2u,

so f = omega / 2 pi. Run it on reference_results/kgrid_e2g/ (folders k<N>_u<u>/ with
siesta.FA), or with --write to create the input folders yourself:

    python e2g_frozen_phonon.py reference_results/kgrid_e2g
    python e2g_frozen_phonon.py --write my_e2g --structure graphene.fdf --calc calc_primitive.fdf \\
        --kgrids 15 16 30 45 --amplitudes 0.005 0.01
"""
import argparse
import glob
import os
import re
import shutil
import sys

import numpy as np

AMU = 1.66053906660e-27          # kg
EV_PER_A2 = 1.602176634e-19 / 1e-20   # eV/Ang^2 -> J/m^2
M_C = 12.011 * AMU


def frequency_thz(fa_path, u):
    """E2g frequency (THz) from the .FA of a +u/-u displaced two-atom cell."""
    rows = np.loadtxt(fa_path, skiprows=1)
    f = rows[:, 1:4]
    w2 = -((f[0, 0] - f[1, 0]) / M_C) / (2.0 * u) * EV_PER_A2
    return np.sqrt(w2) / (2.0 * np.pi) / 1e12


def write_folders(out, structure, calc, kgrids, amplitudes, pseudo):
    from stb.core import structure_io
    text = open(structure).read()
    s = structure_io.to_pymatgen(structure_io.read_fdf(structure))
    calc_text = open(calc).read()
    for k in kgrids:
        for u in amplitudes:
            d = os.path.join(out, f"k{k}_u{u:g}")
            os.makedirs(d, exist_ok=True)
            cart = s.cart_coords.copy()
            cart[0, 0] += u
            cart[1, 0] -= u
            frac = s.lattice.get_fractional_coords(cart)
            block = "\n".join(["%block AtomicCoordinatesAndAtomicSpecies"]
                              + [f"  {x:.10f} {y:.10f} {z:.10f} 1" for x, y, z in frac]
                              + ["%endblock AtomicCoordinatesAndAtomicSpecies"])
            new = re.sub(r"%block AtomicCoordinatesAndAtomicSpecies.*?%endblock AtomicCoordinatesAndAtomicSpecies",
                         block, text, flags=re.S | re.I)
            with open(os.path.join(d, "structure.fdf"), "w") as fh:
                fh.write(new)
            with open(os.path.join(d, "calc.fdf"), "w") as fh:
                fh.write("%include config_extra.fdf\n" + calc_text)
            with open(os.path.join(d, "config_extra.fdf"), "w") as fh:
                fh.write("# single point, k-grid of this test (first occurrence wins)\n"
                         "MD.TypeOfRun CG\nMD.Steps 0\nMD.NumCGsteps 0\nMD.VariableCell false\n"
                         f"kgrid.MonkhorstPack [{k} {k} 1]\n")
            if pseudo:
                shutil.copy(pseudo, d)
            print(f"wrote {d}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory", nargs="?", default="reference_results/kgrid_e2g")
    ap.add_argument("--write", metavar="DIR", help="write the input folders instead of reading results")
    ap.add_argument("--structure", default="graphene.fdf")
    ap.add_argument("--calc", default="calc_primitive.fdf")
    ap.add_argument("--pseudo", default=None, help="C.psf to copy into each folder")
    ap.add_argument("--kgrids", type=int, nargs="+", default=[15, 16, 30, 45])
    ap.add_argument("--amplitudes", type=float, nargs="+", default=[0.005, 0.01])
    a = ap.parse_args()
    if a.write:
        write_folders(a.write, a.structure, a.calc, a.kgrids, a.amplitudes, a.pseudo)
        return 0
    rows = []
    for d in glob.glob(os.path.join(a.directory, "k*_u*")):
        m = re.match(r"k(\d+)_u([0-9.]+)$", os.path.basename(d))
        fa = glob.glob(os.path.join(d, "*.FA"))
        if m and fa:
            rows.append((int(m.group(1)), float(m.group(2)), frequency_thz(fa[0], float(m.group(2)))))
    if not rows:
        print(f"No finished k<N>_u<u>/ folders in {a.directory}")
        return 1
    print(f"{'k-grid (primitive)':>18} | {'u (Ang)':>7} | {'E2g (THz)':>9} | {'(cm^-1)':>7} | contains K?")
    print("-" * 64)
    for k, u, f in sorted(rows):
        print(f"{str(k) + ' x ' + str(k):>18} | {u:7.3f} | {f:9.2f} | {f * 33.35641:7.0f} | "
              f"{'yes' if k % 3 == 0 else 'no'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
