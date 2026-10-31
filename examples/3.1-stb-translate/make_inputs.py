#!/usr/bin/env python3
"""Writes the hand-made inputs of the 3.1 example (the real ones -- si_bulk.fdf, si111_2x2.fdf,
si111_adatom_relaxed.STRUCT_OUT -- are copies of the 7.1/5.10 runs):

    nacl_asym.cif         rock salt as a database gives it: space group Fm-3m and its 192 operations, 2 atoms in the
                          asymmetric unit (pymatgen's CifWriter with symmetry)
    si_bohr.fdf           the 2-atom primitive cell of Si with LatticeConstant 1.0 Bohr (the convention of phonopy's
                          SIESTA writer): the same crystal as si_bulk.fdf in other units
    water.xyz             an H2O molecule (no cell: an XYZ file has none -- --lattice-vectors gives it a 15 Ang box)
    si_noisy_<A>.fdf      si_bulk.fdf with every atom displaced at random by up to A Ang (0.002, 0.02, 0.2): the
                          symmetry-tolerance study
    si_mixup.vasp         POSCAR of si_bulk.fdf with Cartesian coordinates under the word 'Direct' (a classic mix-up)
    si_wrongscale.vasp    POSCAR with the scale factor 0.5 instead of 1.0 (cell and atoms shrunk 2x)

    python3 make_inputs.py      (seconds)
"""
import os

import numpy as np
from ase.build import molecule
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

HERE = os.path.dirname(os.path.abspath(__file__))


def p(name):
    return os.path.join(HERE, name)


def main():
    s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    CifWriter(s, symprec=0.01).write_file(p("nacl_asym.cif"))

    bohr = 0.529177210903
    a = 5.480238 / bohr / 2                 # the fcc primitive vectors in Bohr
    with open(p("si_bohr.fdf"), "w") as f:
        f.write("# Si, 2-atom primitive cell, lattice in Bohr (as phonopy's SIESTA writer does)\n"
                "NumberOfSpecies 1\nNumberOfAtoms 2\n%block ChemicalSpeciesLabel\n 1 14 Si\n"
                "%endblock ChemicalSpeciesLabel\nLatticeConstant 1.0 Bohr\n%block LatticeVectors\n"
                f"  0.000000 {a:.6f} {a:.6f}\n  {a:.6f} 0.000000 {a:.6f}\n  {a:.6f} {a:.6f} 0.000000\n"
                "%endblock LatticeVectors\nAtomicCoordinatesFormat Fractional\n"
                "%block AtomicCoordinatesAndAtomicSpecies\n 0.00 0.00 0.00 1\n 0.25 0.25 0.25 1\n"
                "%endblock AtomicCoordinatesAndAtomicSpecies\n")

    w = molecule("H2O")
    w.positions += 7.5
    with open(p("water.xyz"), "w") as f:
        f.write(f"{len(w)}\nH2O, Cartesian Ang, centred in a 15 Ang box (the box comes from --lattice-vectors)\n")
        for sym, pos in zip(w.get_chemical_symbols(), w.positions):
            f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")

    si = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(p("si_bulk.fdf"))))
    rng = np.random.default_rng(7)
    for amp in (0.002, 0.02, 0.2):
        b = si.copy()
        b.positions += rng.uniform(-amp, amp, b.positions.shape)
        frac = b.get_scaled_positions(wrap=True)
        with open(p(f"si_noisy_{amp:g}.fdf"), "w") as f:
            f.write(f"# si_bulk.fdf with every atom displaced by up to {amp:g} Ang (seed 7)\n"
                    "NumberOfSpecies 1\nNumberOfAtoms 8\n%block ChemicalSpeciesLabel\n 1 14 Si\n"
                    "%endblock ChemicalSpeciesLabel\nLatticeConstant 1.0 Ang\n%block LatticeVectors\n")
            for v in b.cell:
                f.write(f"  {v[0]:.8f} {v[1]:.8f} {v[2]:.8f}\n")
            f.write("%endblock LatticeVectors\nAtomicCoordinatesFormat Fractional\n"
                    "%block AtomicCoordinatesAndAtomicSpecies\n")
            for q in frac:
                f.write(f"  {q[0]:.8f} {q[1]:.8f} {q[2]:.8f} 1\n")
            f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

    def poscar(name, scale, label, rows):
        with open(p(name), "w") as f:
            f.write(f"{label}\n{scale}\n")
            for v in si.cell:
                f.write(f"  {v[0]:.8f} {v[1]:.8f} {v[2]:.8f}\n")
            f.write(f"Si\n{len(si)}\nDirect\n")
            for q in rows:
                f.write(f"  {q[0]:.8f} {q[1]:.8f} {q[2]:.8f}\n")

    poscar("si_mixup.vasp", 1.0, "Si8 -- Cartesian coordinates under 'Direct' (a mix-up)", si.positions)
    poscar("si_wrongscale.vasp", 0.5, "Si8 -- scale factor 0.5 instead of 1.0", si.get_scaled_positions())
    print("inputs written in", HERE)


if __name__ == "__main__":
    main()
