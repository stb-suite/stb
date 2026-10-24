#!/usr/bin/env python3
"""Helper of make_reference.sh: writes the starting structures of the 4.20 example.

    build_inputs.py molecule NAME OUT.fdf [BOX]   a molecule of ASE's G2 set centred in a cubic box (default 14 Ang)
    build_inputs.py rocksalt El1 El2 A OUT.fdf    2-atom primitive cell of a rock-salt crystal, cubic lattice constant A
"""
import sys

import numpy as np


def write_structure(path, cell, symbols, cart, comment):
    from ase.data import atomic_numbers
    labels = list(dict.fromkeys(symbols))
    ids = {lab: k for k, lab in enumerate(labels, 1)}
    frac = (np.asarray(cart, float) @ np.linalg.inv(np.asarray(cell, float))) % 1.0
    with open(path, "w") as f:
        f.write(f"# {comment}\n\nNumberOfSpecies {len(labels)}\nNumberOfAtoms {len(symbols)}\n\n")
        f.write("%block ChemicalSpeciesLabel\n")
        for lab in labels:
            f.write(f" {ids[lab]} {atomic_numbers[lab]} {lab}\n")
        f.write("%endblock ChemicalSpeciesLabel\n\nLatticeConstant 1.0 Ang\n\n%block LatticeVectors\n")
        for row in cell:
            f.write(" %.8f %.8f %.8f\n" % tuple(row))
        f.write("%endblock LatticeVectors\n\nAtomicCoordinatesFormat Fractional\n")
        f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for s, fr in zip(symbols, frac):
            f.write("  %.10f %.10f %.10f %d\n" % (*fr, ids[s]))
        f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")


def molecule(name, out, box=14.0):
    from ase.build import molecule as g2
    m = g2(name)
    pos = m.get_positions()
    pos += box / 2 - (pos.max(0) + pos.min(0)) / 2
    write_structure(out, np.diag([box] * 3), m.get_chemical_symbols(), pos,
                    f"{name} from ASE's G2 set (MP2 geometry), centred in a {box:g} Ang cubic box; "
                    "make_reference.sh relaxes it with SIESTA")


def rocksalt(el1, el2, a, out):
    cell = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    write_structure(out, cell, [el1, el2], [[0, 0, 0], [a / 2, 0, 0]],
                    f"{el1}{el2} rock salt, primitive fcc cell, a = {a:g} Ang (start; make_reference.sh relaxes the cell)")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "molecule":
        molecule(sys.argv[2], sys.argv[3], *(float(x) for x in sys.argv[4:5]))
    elif cmd == "rocksalt":
        rocksalt(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5])
    else:
        sys.exit(__doc__)
