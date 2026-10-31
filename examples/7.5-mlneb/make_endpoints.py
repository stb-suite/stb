#!/usr/bin/env python3
"""Builds the two endpoints of the vacancy hop of the 7.5 example from structure.fdf.

    python3 make_endpoints.py [OUTDIR]        # writes OUTDIR/vac_initial.fdf, OUTDIR/vac_final.fdf (default .)

1. a 2x2x2 supercell of the 8-atom diamond cell (64 Si, a = 10.96 Ang), atoms in ase's repeat() order
   (stb.core.structure_io.replicate_fdf_cell, as stb-supercell writes it);
2. initial: the atom at the origin removed -- a vacancy at (0, 0, 0), 63 atoms;
3. final: the same 63 atoms, with its nearest neighbour at (1/8, 1/8, 1/8) of the supercell (one bond,
   a sqrt(3)/4 = 2.373 Ang away) moved into the empty site. The vacancy has hopped to (1/8, 1/8, 1/8).

Both files keep the SAME atom order: atom 1 is the one that hops, every other atom stays put. NEB needs
that correspondence -- image i of the band interpolates atom j of the initial state to atom j of the final.
"""
import os
import sys
import tempfile

from stb.core import structure_io

HERE = os.path.dirname(os.path.abspath(__file__))


def coordinate_block(lines):
    i0 = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("%block atomiccoordinatesandatomicspecies"))
    i1 = next(i for i in range(i0, len(lines)) if l_end(lines[i]))
    return i0, i1


def l_end(line):
    return line.strip().lower().startswith("%endblock atomiccoordinatesandatomicspecies")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        sc = os.path.join(tmp, "si64.fdf")
        structure_io.replicate_fdf_cell(os.path.join(HERE, "structure.fdf"), (2, 2, 2), sc)
        lines = open(sc).read().splitlines()
    i0, i1 = coordinate_block(lines)
    atoms = [l.split() for l in lines[i0 + 1:i1]]
    assert [float(x) for x in atoms[0][:3]] == [0.0, 0.0, 0.0], "atom 0 of the supercell should sit at the origin"
    assert [float(x) for x in atoms[1][:3]] == [0.125, 0.125, 0.125], "atom 1 should be its neighbour at 1/8"
    a = float(lines[next(i for i, l in enumerate(lines) if "%block LatticeVectors" in l) + 1].split()[0])
    header = ["# Si63: 2x2x2 diamond supercell (a = %.5f Ang) with one vacancy -- %s",
              "# written by make_endpoints.py (examples/7.5-mlneb)"]
    for name, first, what in (("vac_initial.fdf", None, "vacancy at (0,0,0); atom 1 at (1/8,1/8,1/8)"),
                              ("vac_final.fdf", ["0.00000000", "0.00000000", "0.00000000"],
                               "atom 1 hopped into (0,0,0); vacancy now at (1/8,1/8,1/8)")):
        rows = atoms[1:]
        if first:
            rows = [first + rows[0][3:]] + rows[1:]
        text = [header[0] % (a, what), header[1], ""]
        for l in lines[:i0 + 1]:
            if l.startswith("#"):
                continue
            text.append(f"NumberOfAtoms    {len(rows)}" if l.strip().lower().startswith("numberofatoms") else l)
        text += [f"  {float(r[0]):.8f}   {float(r[1]):.8f}   {float(r[2]):.8f}   {r[3]}" for r in rows]
        text += lines[i1:]
        with open(os.path.join(out, name), "w") as f:
            f.write("\n".join(text) + "\n")
        print(f"{os.path.join(out, name)}: 63 Si -- {what}")
    print("hop length: sqrt(3)/8 x %.5f = %.4f Ang" % (a, 3 ** 0.5 / 8 * a))


if __name__ == "__main__":
    main()
