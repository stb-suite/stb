#!/usr/bin/env python
"""Writes the SIESTA regression of the bundled PseudoDojo bank: one X-C molecule per element of
every table (C-C for carbon), in a 12 Ang box.

    python make_cases.py -o cases [--tables dojo dojo-fr ...]

Per case folder <out>/<bank>/<case>/: calc.fdf (self-contained) and the pseudopotentials, with the
blocks the suite writes itself (core.pseudopotentials.ensure_species_blocks: synthetic-atoms for the
sr3plus lanthanides, non-perturbative polarization for Li/Be). Cases:
  sr, sr3plus tables  X-C without spin-orbit (sr3plus: C from dojo)
  fr-pbe tables       X-C without and with spin-orbit
  fr-pbesol tables    X-C with spin-orbit; without it only X-X for Ge, Sb, Te (the others have no
                      scalar-relativistic part and are refused by the suite)
Settings: DZP, PAO.EnergyShift 0.02 Ry, Mesh.Cutoff 420 Ry, no FilterCutoff, 3 SCF iterations,
SCF.MustConverge false. It is a 'does SIESTA run it' test, not an accuracy test.
"""
import argparse
import os
import shutil

from ase.data import atomic_numbers, covalent_radii
from stb.core import pseudopotentials as pp

XC = {"pbe": ("GGA", "PBE"), "pbesol": ("GGA", "PBEsol"), "pw": ("LDA", "PW92")}
BOX = 12.0


def calc_text(xc, species, coords, soc):
    func, auth = XC[xc]
    rows = "".join(f"  {i} {z} {lb}\n" for i, (z, lb) in enumerate(species, 1))
    atoms = "".join(f"  {x:.6f} {y:.6f} {z:.6f} {k}\n" for (x, y, z), k in coords)
    return (f"SystemLabel siesta\nNumberOfAtoms {len(coords)}\nNumberOfSpecies {len(species)}\n"
            f"%block ChemicalSpeciesLabel\n{rows}%endblock ChemicalSpeciesLabel\n"
            f"LatticeConstant 1.0 Ang\n%block LatticeVectors\n  {BOX} 0 0\n  0 {BOX} 0\n  0 0 {BOX}\n"
            "%endblock LatticeVectors\nAtomicCoordinatesFormat Ang\n"
            f"%block AtomicCoordinatesAndAtomicSpecies\n{atoms}%endblock AtomicCoordinatesAndAtomicSpecies\n"
            f"XC.Functional {func}\nXC.Authors {auth}\n"
            f"Spin {'spin-orbit' if soc else 'non-polarized'}\n"
            "PAO.BasisSize DZP\nPAO.EnergyShift 0.02 Ry\nMesh.Cutoff 420 Ry\n"
            "MaxSCFIterations 3\nSCF.MustConverge false\nElectronicTemperature 300 K\n")


def write_case(folder, xc, first, second, src_first, src_second, soc):
    os.makedirs(folder, exist_ok=True)
    za, zb = atomic_numbers[first], atomic_numbers[second]
    d = covalent_radii[za] + covalent_radii[zb]
    c = BOX / 2
    if first == second:
        species = [(za, first)]
        coords = [((c - d / 2, c, c), 1), ((c + d / 2, c, c), 1)]
    else:
        species = [(za, first), (zb, second)]
        coords = [((c - d / 2, c, c), 1), ((c + d / 2, c, c), 2)]
    shutil.copy(os.path.join(src_first, f"{first}.psml"), folder)
    if first != second:
        shutil.copy(os.path.join(src_second, f"{second}.psml"), folder)
    with open(os.path.join(folder, "calc.fdf"), "w") as f:
        f.write(calc_text(xc, species, coords, soc))
    pp.ensure_species_blocks(folder)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="cases")
    ap.add_argument("--tables", nargs="*", default=[n for n in pp.bank_names() if n != "virtual_vault"])
    args = ap.parse_args()
    total = 0
    for bank in args.tables:
        rel, xc, _acc = pp.parse_dojo_name(pp.canonical_bank(bank))
        src = pp.bank_path(bank)
        carbon_src = pp.bank_path("dojo") if rel == "sr3plus" else src
        manifest = pp.bank_manifest(bank)["files"]
        for name in sorted(manifest):
            el = name[:-5]
            base = os.path.join(args.out, pp.canonical_bank(bank))
            plans = []
            if rel in ("sr", "sr3plus"):
                plans.append(("nosoc", False))
            elif xc == "pbe":
                plans += [("nosoc", False), ("soc", True)]
            else:
                plans.append(("soc", True))
                if manifest[name]["scalar_part"]:
                    plans.append(("nosoc_dimer", False))
            for tag, soc in plans:
                partner = el if tag == "nosoc_dimer" or el == "C" else "C"
                write_case(os.path.join(base, f"{el}_{tag}"), xc, el, partner, src, carbon_src, soc)
                total += 1
    print(f"{total} cases written under {args.out}/")


if __name__ == "__main__":
    main()
