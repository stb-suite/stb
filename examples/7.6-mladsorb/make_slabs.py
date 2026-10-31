#!/usr/bin/env python3
"""Writes the Si(111) slabs of the 7.6 example: bulk-terminated (1x1, unreconstructed) Si(111) cut between two
bilayers, so every surface atom keeps ONE dangling bond (a cut through the short bond would leave three), at the
lattice constant of the 7.1 SIESTA training set (5.480238 Ang), 7.5 Ang of vacuum on each side.

    si111_2x2.fdf      2x2 surface cell, 3 bilayers (24 atoms) -- the main system (one adatom = 1/4 monolayer)
    si111_2x2_4bl.fdf  2x2, 4 bilayers (32 atoms)                -- thickness check
    si111_3x3.fdf      3x3, 3 bilayers (54 atoms)                -- coverage check (1/9 monolayer)
"""
from ase.build import diamond111
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io

A_SI = 5.480238   # Ang, 7.1 training reference (SIESTA PBE/DZP)


def slab(n, bilayers):
    s = diamond111("Si", size=(n, n, 2 * bilayers + 2), a=A_SI, vacuum=7.5, orthogonal=False)
    z = s.positions[:, 2]
    s = s[(z > z.min() + 0.1) & (z < z.max() - 0.1)]     # drop the two single-atom layers: whole bilayers left
    s.center(vacuum=7.5, axis=2)
    s.wrap()
    return s


for name, n, bl in (("si111_2x2.fdf", 2, 3), ("si111_2x2_4bl.fdf", 2, 4), ("si111_3x3.fdf", 3, 3)):
    s = slab(n, bl)
    structure_io.write_fdf(structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(s), coord_format="fractional"),
                           name)
    print(f"{name}: {len(s)} atoms, a = {s.cell.lengths()[0]:.4f} Ang, c = {s.cell.lengths()[2]:.3f} Ang")
