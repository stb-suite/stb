#!/usr/bin/env python3
"""Generates tiny synthetic SIESTA .RHO grid files for stb-chargediffAnalysis's
smoke test, via sisl -- a real (if physically meaningless) binary .RHO file
a few KB in size, instead of committing a multi-MB real SIESTA output to the
repo (see test/6-utils/3-cube/o2.RHO for the "commit a real one" convention
used elsewhere; a tiny synthetic grid is preferred here since this test only
needs to exercise the subtraction/shape-check/cube-write PLUMBING, not real
physics). Also writes a matching combined.fdf (for the --cube geometry read)
via stb's own structure_io.write_fdf, so the fixture is self-consistent.

Run with no arguments; writes into its own directory.
"""
import os
import numpy as np
import sisl

from stb.core import structure_io

HERE = os.path.dirname(os.path.abspath(__file__))


def write_rho(path, shape, seed):
    lattice = sisl.Lattice(np.diag([5.0, 5.0, 20.0]))
    geom = sisl.Geometry([[0, 0, 2], [0, 0, 4]], atoms=sisl.Atom(6), lattice=lattice)
    grid = sisl.Grid(shape, geometry=geom)
    rng = np.random.default_rng(seed)
    grid.grid[...] = rng.random(shape).astype(np.float32)
    sisl.get_sile(path, mode="w").write_grid(grid)


def write_geometry_fdf(path):
    structure = structure_io.FdfStructure(
        lattice=np.diag([5.0, 5.0, 20.0]),
        lattice_constant=1.0,
        species=["C"],
        species_meta={"C": {"id": "1", "Z": 6}},
        atoms=[("C", np.array([0.0, 0.0, 0.1])), ("C", np.array([0.0, 0.0, 0.2]))],
        coord_format="fractional",
    )
    structure_io.write_fdf(structure, path)


if __name__ == "__main__":
    write_rho(os.path.join(HERE, "combined.RHO"), (6, 6, 24), seed=0)
    write_rho(os.path.join(HERE, "frag1.RHO"), (6, 6, 24), seed=1)
    write_rho(os.path.join(HERE, "frag2.RHO"), (6, 6, 24), seed=2)
    write_rho(os.path.join(HERE, "badshape.RHO"), (4, 4, 16), seed=3)
    write_geometry_fdf(os.path.join(HERE, "combined.fdf"))
    print("Synthetic .RHO/.fdf fixtures written.")
