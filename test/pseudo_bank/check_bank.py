#!/usr/bin/env python
"""Bundled PseudoDojo bank, checks that need no SIESTA (run by test.sh).

1. every bank folder matches its manifest.json (files, md5) and no file has ekb = 0 projectors;
2. per table, the projector sets SIESTA 5.4.2 needs: every sr/sr3plus file has a scalar-relativistic
   set, every fr file an lj set; without spin-orbit, fr_pbe runs everywhere and fr_pbesol only for
   Ge, Sb and Te;
3. the sr3plus files all have the valence the synthetic-atoms block is written for;
4. the guards of core/pseudopotentials on small fixtures: an ekb = 0 file is refused, fr files with
   no scalar part are an error without spin-orbit, the synthetic-atoms block goes to each folder's
   own species index (ghosts too), Li/Be get the non-perturbative scheme.
Prints one line per check and exits 1 when any fails.
"""
import os
import shutil
import sys
import tempfile

from stb.core import pseudopotentials as pp

FAILED = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}{(': ' + detail) if detail and not ok else ''}")
    if not ok:
        FAILED.append(name)


def fdf(species, extra=""):
    """A minimal run folder input: ChemicalSpeciesLabel from [(Z, label)] plus `extra`."""
    rows = "".join(f"  {i} {z} {lb}\n" for i, (z, lb) in enumerate(species, 1))
    return ("SystemLabel t\nXC.Functional GGA\nXC.Authors PBE\nPAO.BasisSize DZP\n" + extra
            + f"NumberOfSpecies {len(species)}\n%block ChemicalSpeciesLabel\n{rows}%endblock ChemicalSpeciesLabel\n")


def main():
    # 1-3: the bank itself
    for name in pp.bank_names():
        problems = pp.verify_bank(name)
        check(f"{name}: files match manifest.json, no ekb = 0", not problems, "; ".join(problems[:5]))
        if name == "virtual_vault":
            continue
        manifest = pp.bank_manifest(name)
        files = manifest["files"]
        relativity = pp.parse_dojo_name(name)[0]
        if relativity in ("sr", "sr3plus"):
            bad = [f for f, d in files.items() if not d["scalar_part"] or d["lj"]]
            check(f"{name}: scalar-relativistic set in every file, no lj", not bad, ", ".join(bad))
        else:
            no_lj = [f for f, d in files.items() if not d["lj"]]
            check(f"{name}: lj set in every file", not no_lj, ", ".join(no_lj))
            scalar = sorted(f[:-5] for f, d in files.items() if d["scalar_part"])
            if "pbesol" in name:
                check(f"{name}: scalar part only in Ge, Sb, Te", scalar == ["Ge", "Sb", "Te"], ", ".join(scalar))
            else:
                check(f"{name}: scalar part in every file", len(scalar) == len(files))
        if relativity == "sr3plus":
            folder = pp.bank_path(name)
            odd = [f for f in files if {(n, l, o) for n, l, o in
                                        pp.read_pseudo_header(os.path.join(folder, f))["valence"]}
                   != pp.SR3PLUS_VALENCE]
            check(f"{name}: every file has valence 5s2 5p6 5d1 6s2", not odd, ", ".join(odd))

    # 4: guards on fixtures
    sr = pp.bank_path("dojo")
    work = tempfile.mkdtemp(prefix="stb_pseudo_bank_")
    try:
        # ekb = 0: a copy of a bank file with one projector's ekb set to zero
        bad_dir = os.path.join(work, "bad")
        os.makedirs(bad_dir)
        text = open(os.path.join(sr, "Pb.psml")).read()
        import re
        text = re.sub(r'ekb="[^"]*"', 'ekb="0"', text, count=1)
        open(os.path.join(bad_dir, "Pb.psml"), "w").write(text)
        refused = False
        try:
            pp.copy_pseudo(bad_dir, "Pb", work)
        except SystemExit:
            refused = True
        check("ekb = 0 file refused by copy_pseudo, nothing copied",
              refused and not os.path.exists(os.path.join(work, "Pb.psml")))
        run = os.path.join(work, "run_ekb")
        os.makedirs(run)
        shutil.copy(os.path.join(bad_dir, "Pb.psml"), run)
        open(os.path.join(run, "calc.fdf"), "w").write(fdf([(82, "Pb")]))
        levels = [lvl for lvl, _m in pp.check_pseudo_level(run)]
        check("ekb = 0 file in a hand-made folder is an [ERROR]", "error" in levels)

        # fr_pbesol without spin-orbit
        run = os.path.join(work, "run_fr")
        os.makedirs(run)
        shutil.copy(os.path.join(pp.bank_path("dojo-fr-pbesol"), "Fe.psml"), run)
        open(os.path.join(run, "calc.fdf"), "w").write(fdf([(26, "Fe")]).replace("PBE\n", "PBEsol\n"))
        check("fr-pbesol Fe without spin-orbit is an [ERROR]",
              any(lvl == "error" and "scalar-relativistic" in m for lvl, m in pp.check_pseudo_level(run)))
        open(os.path.join(run, "calc.fdf"), "w").write(
            fdf([(26, "Fe")], "Spin spin-orbit\n").replace("PBE\n", "PBEsol\n"))
        check("fr-pbesol Fe with spin-orbit is accepted", not pp.check_pseudo_level(run))

        # sr3plus: block at each folder's own index, ghosts too
        run = os.path.join(work, "run_ln")
        os.makedirs(run)
        ln = pp.bank_path("dojo-sr3plus")
        shutil.copy(os.path.join(sr, "C.psml"), os.path.join(run, "C.psml"))
        shutil.copy(os.path.join(ln, "Ce.psml"), os.path.join(run, "Ce_ghost.psml"))
        shutil.copy(os.path.join(ln, "Gd.psml"), os.path.join(run, "Gd.psml"))
        open(os.path.join(run, "calc.fdf"), "w").write(fdf([(6, "C"), (-58, "Ce_ghost"), (64, "Gd")]))
        check("sr3plus without the block is an [ERROR]",
              any(lvl == "error" and "synthetic-atoms" in m for lvl, m in pp.check_pseudo_level(run)))
        pp.ensure_species_blocks(run)
        entries = pp._synthetic_block_entries(open(os.path.join(run, "calc.fdf")).read())
        check("synthetic-atoms written at species 2 (ghost Ce) and 3 (Gd)",
              entries == {2: pp.SR3PLUS_SYNTHETIC, 3: pp.SR3PLUS_SYNTHETIC}, str(entries))
        check("sr3plus with the block passes", not pp.check_pseudo_level(run))

        # Li and Be ghosts, no polarization scheme in the input
        run = os.path.join(work, "run_be")
        os.makedirs(run)
        for el, label in (("Be", "Be"), ("Li", "Li_ghost"), ("O", "O")):
            shutil.copy(os.path.join(sr, f"{el}.psml"), os.path.join(run, f"{label}.psml"))
        open(os.path.join(run, "calc.fdf"), "w").write(fdf([(4, "Be"), (-3, "Li_ghost"), (8, "O")]))
        pp.ensure_species_blocks(run)
        check("Be and Li_ghost get the non-perturbative scheme",
              pp.missing_nonperturbative(open(os.path.join(run, "calc.fdf")).read()) == [])
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print(f"\n{len(FAILED)} check(s) failed" if FAILED else "\nAll checks passed")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
