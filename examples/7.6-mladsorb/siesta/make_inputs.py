#!/usr/bin/env python3
"""SIESTA inputs of the 7.6 example's check (README section 7): the clean slab and the four adatom minima of the main
case, each to be RELAXED with SIESTA at the 7.1 training level of theory, starting from the model's geometry.

    calc.fdf           the 7.1 training calc.fdf (PBE, DZP, 300 Ry, %include structure.fdf), with
                       %include config_extra.fdf on top: k-grid 4x4x1 (k x a ~ 31 Ang in-plane, as the training's 33),
                       CG to 0.05 eV/Ang (the model's --fmax), fixed cell, the bottom bilayer fixed (Geometry.Constraints,
                       the 8 lowest atoms -- the atoms --relax-depth 4 fixes)
    structure.fdf      the model's relaxed geometry (reference_run/main/rank*.fdf, mladsorb_slab.xyz)
    Si.psml            the 4.17 pseudopotential (PseudoDojo PBE standard)

    python3 siesta/make_inputs.py <output folder> [<pseudopotential Si.psml>]
    bash siesta/run_siesta.sh <output folder>/{T4,H3,T1,B2,slab}      # 6 MPI processes each, ~1 h each
"""
import os
import shutil
import sys

import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(HERE, "reference_run", "main")
CALC = os.path.join(HERE, "..", "7.5-mlneb", "calc.fdf")      # the 7.1 training calc.fdf, shipped with 7.5
JOBS = {"slab": None, "T1": "rank*_ontop_site1_Si.fdf", "T4": "rank*_ontop_site2_Si.fdf",
        "H3": "rank*_bridge_site3_Si.fdf", "B2": "rank*_bridge_site5_Si.fdf"}
EXTRA = ("# stb 7.6 example: relaxation at the 7.1 level of theory, bottom bilayer fixed\n"
         "kgrid.MonkhorstPack   [4 4 1]\nMD.TypeOfRun          CG\nMD.NumCGsteps         300\n"
         "MD.MaxForceTol        0.05 eV/Ang\nMD.VariableCell       false\nWriteForces           true\n"
         "WriteCoorXmol         true\n%block Geometry.Constraints\n  atom {fixed}\n%endblock Geometry.Constraints\n")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    pseudo = sys.argv[2] if len(sys.argv) > 2 else None
    import glob
    for name, pattern in JOBS.items():
        d = os.path.join(out, name)
        os.makedirs(d, exist_ok=True)
        if pattern is None:
            at = read(os.path.join(MAIN, "mladsorb_slab.xyz"))
            st = structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(at), coord_format="fractional")
        else:
            st = structure_io.read_fdf(glob.glob(os.path.join(MAIN, pattern))[0])
        structure_io.write_fdf(st, os.path.join(d, "structure.fdf"))
        at = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(st))
        fixed = sorted(int(i) + 1 for i in np.argsort(at.positions[:24, 2])[:8])
        with open(os.path.join(d, "calc.fdf"), "w") as f:
            f.write("%include config_extra.fdf\n" + open(CALC).read())
        with open(os.path.join(d, "config_extra.fdf"), "w") as f:
            f.write(EXTRA.format(fixed=" ".join(map(str, fixed))))
        if pseudo:
            shutil.copy(pseudo, d)
        print(f"{d}: {len(at)} atoms, fixed {fixed}")


if __name__ == "__main__":
    main()
