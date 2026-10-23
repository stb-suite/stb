#!/usr/bin/env python3
"""Builds a finished-looking stb-chargediffPrep run for stb-chargediffAnalysis's test, with a KNOWN answer.

usage: make_synthetic_run.py <run dir>   (the run dir must already hold the folders written by stb-chargediffPrep)

Every atom of the 4-atom chain C N C N (3 Ang apart along c, listed interleaved so the fragment folders, which group
their atoms by species label, list them in a different order from the combined one) carries one electron as a
Gaussian (sigma 0.5 Ang). The combined density is the sum of the fragment densities plus a planted transfer:
0.05 e from the first C to the first N and 0.02 e from the second C to the second N. The analysis must therefore find

    integral of Delta rho = 0;  charge moved = 0.07 e;  C fragment -0.07 e, N fragment +0.07 e;
    per atom  -0.05, +0.05, -0.02, +0.02 e;  E_int = -0.150 eV (energies written below)

and a cube whose positive part, summed over its own voxel volume in bohr^3, is 0.07 e (the unit of the cube file)."""
import os
import sys

import numpy as np
import sisl

from stb.core import structure_io

SIGMA = 0.5                                       # Ang
TRANSFER = {0: 0.05, 2: 0.02}                     # index of the donor atom (C) -> electrons given to the next atom (N)
ENERGIES = {"combined": -100.000, "fragment1": -49.900, "fragment2": -49.950}   # eV; E_int = -0.150 eV
MESH = 400.0

OUT_TEMPLATE = """Executable      : siesta
Version         : 5.4.2

InitMesh: MESH = {n[0]} x {n[1]} x {n[2]} = {total}
InitMesh: Mesh cutoff (required, used) =   {mesh:.3f}   {mesh_used:.3f} Ry
Total number of electrons:     {electrons:.6f}
siesta: FreeEng =  {energy:.6f}
SCF cycle converged after 11 iterations
"""


def gaussian(shape, lattice, center, electrons):
    """periodic Gaussian (minimum image) normalised to `electrons` on this grid; (e/Ang^3)"""
    frac = [(np.arange(n) + 0.0) / n for n in shape]
    inv = np.linalg.inv(lattice)
    c = np.asarray(center) @ inv
    d = []
    for k in range(3):
        f = frac[k] - c[k]
        f -= np.round(f)
        d.append(f * np.linalg.norm(lattice[k]))     # orthogonal cell: fraction x length
    r2 = d[0][:, None, None] ** 2 + d[1][None, :, None] ** 2 + d[2][None, None, :] ** 2
    g = np.exp(-r2 / (2 * SIGMA ** 2))
    dv = abs(np.linalg.det(lattice)) / g.size
    return g * (electrons / (g.sum() * dv))


def write_rho(path, shape, lattice, positions, grid):
    geom = sisl.Geometry(positions, atoms=sisl.Atom(6), lattice=sisl.Lattice(lattice))
    g = sisl.Grid(shape, geometry=geom)
    g.grid[...] = grid.astype(np.float32)
    sisl.get_sile(path, mode="w").write_grid(g)


def main(run):
    comb = structure_io.read_fdf(os.path.join(run, "combined", "structure.fdf"))
    lattice = np.asarray(comb.lattice)
    pos = np.array([p for _, p in comb.atoms], float)
    if comb.coord_format != "cartesian":
        pos = pos @ lattice
    shape = (32, 32, 72)
    symbols = [s.split("_")[0] for s, _ in comb.atoms]
    parts = [gaussian(shape, lattice, p, 1.0) for p in pos]
    frag = {"fragment1": sum(parts[i] for i in range(len(pos)) if symbols[i] == "C"),
            "fragment2": sum(parts[i] for i in range(len(pos)) if symbols[i] == "N")}
    combined = frag["fragment1"] + frag["fragment2"]
    for donor, e in TRANSFER.items():
        combined = combined - gaussian(shape, lattice, pos[donor], e) + gaussian(shape, lattice, pos[donor + 1], e)
    grids = {"combined": combined, **frag}
    for name, grid in grids.items():
        folder = os.path.join(run, name)
        write_rho(os.path.join(folder, "siesta.RHO"), shape, lattice, pos, grid)
        dv = abs(np.linalg.det(lattice)) / grid.size
        with open(os.path.join(folder, "calc.out"), "w") as fh:
            fh.write(OUT_TEMPLATE.format(n=shape, total=int(np.prod(shape)), mesh=MESH, mesh_used=MESH,
                                         electrons=float(grid.sum() * dv), energy=ENERGIES[name]))
    print("synthetic finished run written in", run)


if __name__ == "__main__":
    main(sys.argv[1])
