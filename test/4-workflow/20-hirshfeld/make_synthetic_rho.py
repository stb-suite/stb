#!/usr/bin/env python3
"""Generates a tiny synthetic (physically meaningless) SIESTA .RHO for a
given structure.fdf, via sisl -- exercises the Hirshfeld-I workflow's
PLUMBING (file I/O, manifest chaining, sign decision, convergence loop)
without a real SIESTA run, the same "commit a generator script, not a
multi-KB binary" convention as
test/4-workflow/19-chargediff/analysis/make_synthetic_rho.py. Density =
sum of one Gaussian bump per atom actually present in the given
structure.fdf (works unmodified for both a single-isolated-atom folder
and the multi-atom combined folder), minimum-image aware -- same
convention core/hirshfeld.py itself uses.

Usage: make_synthetic_rho.py <structure.fdf> <output.RHO> <n_grid> <amp1,amp2,...> <sigma1,sigma2,...>
"""
import sys

import numpy as np
import sisl

from stb.core import structure_io, kspace


def main():
    structure_path, out_path, n_grid_str, amps_str, sigmas_str = sys.argv[1:6]
    n_grid = int(n_grid_str)
    amps = [float(x) for x in amps_str.split(",")]
    sigmas = [float(x) for x in sigmas_str.split(",")]

    structure = structure_io.read_fdf(structure_path)
    positions = np.array([pos for _, pos in structure.atoms])
    frac = kspace.to_fractional(positions, structure.lattice,
                                 structure.coord_format == "cartesian")
    assert len(frac) == len(amps) == len(sigmas), \
        f"structure has {len(frac)} atom(s) but got {len(amps)} amplitude(s)/{len(sigmas)} sigma(s)"

    shape = (n_grid, n_grid, n_grid)
    fx, fy, fz = np.meshgrid(np.arange(n_grid) / n_grid, np.arange(n_grid) / n_grid,
                              np.arange(n_grid) / n_grid, indexing='ij')
    grid_frac = np.stack([fx, fy, fz], axis=-1)

    total = np.zeros(shape, dtype=np.float32)
    for i in range(len(frac)):
        disp = grid_frac - frac[i]
        disp -= np.round(disp)
        cart = disp @ structure.lattice
        r2 = np.sum(cart**2, axis=-1)
        total += (amps[i] * np.exp(-r2 / (2 * sigmas[i]**2))).astype(np.float32)

    lattice = sisl.Lattice(structure.lattice)
    cart_positions = frac @ structure.lattice
    geom = sisl.Geometry(cart_positions, atoms=sisl.Atom(1), lattice=lattice)
    grid = sisl.Grid(shape, geometry=geom)
    grid.grid[...] = total
    sisl.get_sile(out_path, mode="w").write_grid(grid)
    print(f"Wrote synthetic RHO '{out_path}' from '{structure_path}' "
          f"({len(frac)} atom(s), grid {shape}).")


if __name__ == "__main__":
    main()
