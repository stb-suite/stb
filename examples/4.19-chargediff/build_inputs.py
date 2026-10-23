#!/usr/bin/env python3
"""Helper of make_reference.sh: writes the starting structures of the 4.19 example and the derived folders
that stb-chargediffPrep does not write (the distance scan, the fragments without ghost atoms, the isolated
monomers). Every structure keeps its atom order, so --fragment-sizes splits it as intended.

    build_inputs.py nh3bh3   OUT.fdf               ammonia borane, staggered C3v start (B-N 1.66 Ang), 12 Ang box
    build_inputs.py graphene OUT.fdf               3x3 graphene (18 C) + Li above a hollow site, 20 Ang along c
    build_inputs.py scan     DIMER.fdf R OUT.fdf   the dimer with the second molecule moved rigidly along O1->O2
                                                   so that d(O...O) = R (Ang)
    build_inputs.py alone    FRAGMENT_DIR OUT_DIR  copy of a fragment folder WITHOUT its ghost atoms: the fragment
                                                   in its own basis only (the uncorrected, BSSE-affected energy)
    build_inputs.py monomer  STRUCT.fdf FIRST N OUT.fdf   atoms FIRST..FIRST+N-1 (1-based) as a structure of their own
"""
import os
import re
import shutil
import sys

import numpy as np

from stb.core import structure_io


def write_structure(path, cell, labels_z, atoms, comment):
    """atoms: [(label, cartesian xyz)] in the order to keep; labels_z: {label: Z} in species-id order"""
    ids = {lab: k for k, lab in enumerate(labels_z, 1)}
    inv = np.linalg.inv(np.asarray(cell, float))
    with open(path, "w") as f:
        f.write(f"# {comment}\n\nNumberOfSpecies {len(labels_z)}\nNumberOfAtoms {len(atoms)}\n\n")
        f.write("%block ChemicalSpeciesLabel\n")
        for lab, z in labels_z.items():
            f.write(f" {ids[lab]} {z} {lab}\n")
        f.write("%endblock ChemicalSpeciesLabel\n\nLatticeConstant 1.0 Ang\n\n%block LatticeVectors\n")
        for row in cell:
            f.write(" %.8f %.8f %.8f\n" % tuple(row))
        f.write("%endblock LatticeVectors\n\nAtomicCoordinatesFormat Fractional\n")
        f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for lab, xyz in atoms:
            fr = np.asarray(xyz, float) @ inv
            f.write("  %.10f %.10f %.10f %d\n" % (*fr, ids[lab]))
        f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")


def nh3bh3(out):
    box = 12.0
    c = np.array([box / 2] * 3)
    b, n = c + [0, 0, -0.83], c + [0, 0, 0.83]
    atoms = [("N", n), ("H", None), ("H", None), ("H", None), ("B", b), ("H", None), ("H", None), ("H", None)]
    # N-H 1.02 Ang at 110 deg from the N->B... direction; B-H 1.21 Ang at 104 deg from the B->N direction; staggered
    for k in range(3):
        phi = np.radians(120 * k)
        t = np.radians(180 - 110)
        atoms[1 + k] = ("H", n + 1.02 * np.array([np.sin(t) * np.cos(phi), np.sin(t) * np.sin(phi), np.cos(t)]))
        phi = np.radians(120 * k + 60)
        t = np.radians(180 - 104)
        atoms[5 + k] = ("H", b + 1.21 * np.array([np.sin(t) * np.cos(phi), np.sin(t) * np.sin(phi), -np.cos(t)]))
    write_structure(out, np.diag([box] * 3), {"N": 7, "B": 5, "H": 1}, atoms,
                    "Ammonia borane H3N-BH3, staggered C3v starting geometry (make_reference.sh relaxes it). "
                    "Atom order: N H H H (donor NH3), then B H H H (acceptor BH3) -- --fragment-sizes 4,4.")


def graphene(out):
    a, vac, n = 2.467, 20.0, 3
    a1, a2 = np.array([a, 0, 0]), np.array([a / 2, a * np.sqrt(3) / 2, 0])     # 60 deg: (a1 + a2)/3 is a C-C bond
    cell = np.array([n * a1, n * a2, [0, 0, vac]])
    z0 = vac / 2
    atoms = []
    for i in range(n):
        for j in range(n):
            o = i * a1 + j * a2
            atoms.append(("C", o + (a1 + a2) / 3 + [0, 0, z0]))
            atoms.append(("C", o + 2 * (a1 + a2) / 3 + [0, 0, z0]))
    # hollow site: the centre of a hexagon (the origin of a unit cell is a hexagon centre here)
    atoms.append(("Li", (a1 + a2) + [0, 0, z0 + 1.75]))
    write_structure(out, cell, {"C": 6, "Li": 3}, atoms,
                    "Li on 3x3 graphene (18 C), hollow site, 1.75 Ang above the sheet; a = 2.467 Ang; 20 Ang along c. "
                    "Atom order: 18 C, then Li -- --fragment-sizes 18,1.")


def scan(dimer, r, out):
    st = structure_io.read_fdf(dimer)
    lat = np.asarray(st.lattice)
    pos = np.array([p for _, p in st.atoms], float)
    if st.coord_format != "cartesian":
        pos = pos @ lat
    d = pos[3] - pos[0]
    shift = (float(r) - np.linalg.norm(d)) * d / np.linalg.norm(d)
    pos[3:] += shift
    labels = {s: int(st.species_meta[s]["Z"]) for s in st.species}
    write_structure(out, lat, labels, [(st.atoms[i][0], pos[i]) for i in range(len(pos))],
                    f"Water dimer with the donor moved rigidly along O1->O2: d(O...O) = {float(r):.3f} Ang")


def alone(frag_dir, out_dir):
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    for name in os.listdir(frag_dir):          # the inputs only; the ghost pseudopotentials are not needed
        if name in ("structure.fdf", "calc.fdf", "config_extra.fdf") or (name.endswith(".psml") and "_ghost" not in name):
            shutil.copy(os.path.join(frag_dir, name), out_dir)
    path = os.path.join(out_dir, "structure.fdf")
    text = open(path).read()
    m = re.search(r"%block ChemicalSpeciesLabel\n(.*?)%endblock", text, re.S)
    keep_species, new_id = [], {}
    for line in m.group(1).strip().split("\n"):
        k, z, lab = line.split()
        if not lab.endswith("_ghost"):
            new_id[k] = str(len(new_id) + 1)            # SIESTA wants the ids 1..N without gaps
            keep_species.append(f" {new_id[k]}   {z}   {lab}")
    m2 = re.search(r"%block AtomicCoordinatesAndAtomicSpecies\n(.*?)%endblock", text, re.S)
    keep_atoms = []
    for line in m2.group(1).strip().split("\n"):
        x, y, z, k = line.split()[:4]
        if k in new_id:
            keep_atoms.append(f"  {x}   {y}   {z}   {new_id[k]}")
    text = text.replace(m.group(1), "\n".join(keep_species) + "\n").replace(m2.group(1), "\n".join(keep_atoms) + "\n")
    text = re.sub(r"NumberOfSpecies\s+\d+", f"NumberOfSpecies    {len(keep_species)}", text, flags=re.I)
    text = re.sub(r"NumberofAtoms\s+\d+", f"NumberofAtoms      {len(keep_atoms)}", text, flags=re.I)
    open(path, "w").write("# the fragment WITHOUT ghost atoms (its own basis only), from " + frag_dir + "\n" + text)


def monomer(src, first, n, out):
    st = structure_io.read_fdf(src)
    lat = np.asarray(st.lattice)
    pos = np.array([p for _, p in st.atoms], float)
    if st.coord_format != "cartesian":
        pos = pos @ lat
    idx = range(int(first) - 1, int(first) - 1 + int(n))
    used = [st.atoms[i][0] for i in idx]
    labels = {s: int(st.species_meta[s]["Z"]) for s in st.species if s in used}
    write_structure(out, lat, labels, [(st.atoms[i][0], pos[i]) for i in idx],
                    f"atoms {first}..{int(first) + int(n) - 1} of {os.path.basename(src)}, alone")


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"nh3bh3": nh3bh3, "graphene": graphene, "scan": scan, "alone": alone, "monomer": monomer}[cmd](*args)
