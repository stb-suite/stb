#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import os
import sys
import math
import difflib
import argparse
from datetime import datetime

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_pyxtal, require_mace

REPORT_FILE = "stb_crystalcast_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-unitcell/stb-supercell/stb-fetch/
# stb-kgrid/stb-mlrelax (core/kspace.py's other callers) -- used to detect
# which axes of a generated/transformed structure are vacuum-padded, for
# both the structure-validation density check (skipped on a vacuum axis)
# and to tell a --dim 0 boxed cluster (vacuum on all 3 axes) apart from a
# --dim 1/2 structure at report time.
VACUUM_GAP_ANG = 10.0

MOLECULE_FILE_EXTENSIONS = ("xyz", "gjf", "g03", "json")

DIM_LABELS = {
    3: "space group",
    2: "layer group",
    1: "rod group",
    0: "point group",
}


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as stb-unitcell/
    stb-crystalbuilder/stb-passivate/stb-amorphize's own _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def resolve_group(spec):
    """Returns spec as an int if it looks like one, else as-is (a symbol string).

    Covers all four --dim cases: space group (int or e.g. 'Fm-3m'), layer/rod
    group (int), or point group (int or Schoenflies symbol e.g. 'D3d') -- pyxtal
    accepts either form directly for its `group` argument in every dimension.
    """
    try:
        return int(spec)
    except ValueError:
        return spec


def check_lattice_type(a, b, c, alpha, beta, gamma, ltype, tol=1e-3):
    """Returns an error string if a/b/c/alpha/beta/gamma don't satisfy the
    geometric constraints of `ltype` (e.g. cubic needs a==b==c and all angles
    90), or None if they're compatible. pyxtal's own Lattice.from_para()
    doesn't validate this itself -- it silently builds whatever cell the raw
    numbers describe and just tags it with `ltype`, so an incompatible
    --lattice would otherwise generate a structure with the wrong symmetry
    and no error at all (only caught after the fact, and only for --dim 3,
    by the requested-vs-detected space-group check much later).
    """
    def close(x, y):
        return math.isclose(x, y, abs_tol=tol)

    satisfied = {
        "cubic": close(a, b) and close(b, c) and close(alpha, 90) and close(beta, 90) and close(gamma, 90),
        "tetragonal": close(a, b) and close(alpha, 90) and close(beta, 90) and close(gamma, 90),
        "orthorhombic": close(alpha, 90) and close(beta, 90) and close(gamma, 90),
        "hexagonal": close(a, b) and close(alpha, 90) and close(beta, 90) and close(gamma, 120),
        "trigonal": close(a, b) and close(alpha, 90) and close(beta, 90) and close(gamma, 120),
        "rhombohedral": close(a, b) and close(b, c) and close(alpha, beta) and close(beta, gamma),
        "monoclinic": close(alpha, 90) and close(gamma, 90),
        "triclinic": True,
        "spherical": True,
        "ellipsoidal": True,
    }.get(ltype)

    if satisfied is None:
        return None  # unrecognized ltype -- don't block on something we don't know how to check
    if not satisfied:
        return (f"--lattice {a} {b} {c} {alpha} {beta} {gamma} is not a valid {ltype} cell "
                f"(the crystal system required by this --group/--dim).")
    return None


def resolve_output_paths(output, count):
    """Returns a list of `count` output paths: just [output] for count == 1,
    else numbered '<stem>_<i><ext>' (1-indexed) -- the same convention used
    by --count > 1 in generation mode, reused here for --subgroup/--supergroup's
    multiple candidates. Raises ValueError if count > 1 and any numbered name
    already exists, to avoid silently overwriting unrelated files.
    """
    if count == 1:
        return [output]
    stem, ext = os.path.splitext(output)
    ext = ext or ".fdf"
    names = [f"{stem}_{i}{ext}" for i in range(1, count + 1)]
    preexisting = [name for name in names if os.path.exists(name)]
    if preexisting:
        raise ValueError(
            f"refusing to overwrite existing file(s): {', '.join(preexisting)} -- "
            "move them aside or choose a different -o.")
    return names


def molecule_to_boxed_structure(molecule, vacuum):
    """Places a non-periodic pyxtal dim=0 cluster (a pymatgen Molecule) into an
    orthorhombic box, vacuum Ang from the outermost atom on every side --
    same convention as stb-molecule's --vacuum, needed because dim=0 has no
    lattice of its own to write into an .fdf.
    """
    coords = molecule.cart_coords
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    lengths = (maxs - mins) + 2 * vacuum
    shifted = coords - mins + vacuum
    species = [str(site.specie) for site in molecule]
    lattice = Lattice.orthorhombic(*lengths)
    return Structure(lattice, species, shifted, coords_are_cartesian=True)


def structure_to_fdf(pmg_structure, species_order=None, coord_format="fractional"):
    """Builds an FdfStructure from a pymatgen Structure, assigning fresh
    species ids from scratch (no prior species_meta to preserve) -- shared by
    every write path in this file (generation, --substitute, --subgroup,
    --supergroup) since none of them start from an existing .fdf's own
    species numbering.

    `species_order`, if given, fixes the order species ids are assigned in
    (the plain-atomic generation case: the user's own --species order).
    Otherwise ids are assigned in the structure's own site order -- used for
    --molecular (species aren't known until after generation, so there's no
    user order to preserve) and for the transform modes (--substitute/
    --subgroup/--supergroup), which have no original --species list at all.
    """
    symbols = species_order if species_order is not None \
        else dict.fromkeys(site.specie.symbol for site in pmg_structure)
    species_meta = {}
    for symbol in symbols:
        species_meta = structure_io.ensure_species_id(species_meta, symbol)
    return structure_io.from_pymatgen(pmg_structure, species_meta=species_meta, coord_format=coord_format)


def validate_molecule_species(entries):
    """Checks each --species entry is usable in --molecular mode before the
    generation loop starts, rather than letting pyxtal fail mid-batch.

    A bare name (e.g. 'H2O') must match pyxtal's small bundled molecule
    collection -- checked here explicitly because pyxtal's own lookup
    (molecule_collection, a process-wide singleton) silently reuses the last
    successfully resolved molecule instead of raising when a later name in
    the same run doesn't match, e.g. --species H2O BadName would silently
    build BadName out of H2O's geometry. A dotted entry is a path (.xyz,
    .gjf, .g03, .json) or a SMILES string (.smi, needs the optional 'rdkit'
    package -- not pre-validated here, pyxtal raises its own clear error).
    """
    from pyxtal.database.collection import Collection
    known_names = list(Collection("molecules"))
    lower_names = {n.lower() for n in known_names}

    for entry in entries:
        parts = entry.split(".")
        if len(parts) > 1:
            ext = parts[-1].lower()
            if ext in MOLECULE_FILE_EXTENSIONS:
                if not os.path.exists(entry):
                    raise ValueError(f"--species '{entry}': file not found.")
            elif ext != "smi":
                raise ValueError(
                    f"--species '{entry}': unsupported file extension '.{parts[-1]}' "
                    f"(expected one of .{', .'.join(MOLECULE_FILE_EXTENSIONS)}, or .smi "
                    "for a SMILES string).")
        elif entry.lower() not in lower_names:
            msg = f"--species '{entry}' is not a known molecule name (case-insensitive)."
            suggestions = difflib.get_close_matches(entry, known_names, n=3)
            if suggestions:
                msg += f" Did you mean: {', '.join(suggestions)}?"
            msg += (" Run with --list-molecules to see all available names, or give a "
                     "path (.xyz/.gjf/.g03/.json) or SMILES string (name.smi) for a "
                     "custom molecule.")
            raise ValueError(msg)


def parse_substitutions(entries):
    """Parses ['Cl:F', 'Na:K'] into {'Cl': 'F', 'Na': 'K'}, validating both
    sides are real element symbols.
    """
    subs = {}
    for entry in entries:
        parts = [p.strip() for p in entry.split(":")]
        if len(parts) != 2 or not parts[0] or not parts[1]:
            raise ValueError(f"--substitute '{entry}': expected OLD:NEW, e.g. 'Cl:F'.")
        old, new = parts
        for symbol in (old, new):
            Element(symbol)  # raises ValueError with a clear message if invalid
        subs[old] = new
    return subs


def load_crystal_from_file(args, f_out):
    """Reads -f/--file and builds a pyxtal object via from_seed(). Shared by
    --substitute/--subgroup/--supergroup -- all three operate on an existing
    (always 3D-periodic) structure instead of generating one from scratch.
    """
    from pyxtal import pyxtal

    if not args.file:
        _fail("this mode requires -f/--file.", f_out)
    if not os.path.exists(args.file):
        _fail(f"File '{args.file}' not found.", f_out)

    try:
        structure = structure_io.read_fdf(args.file)
        pmg_structure = structure_io.to_pymatgen(structure)
    except (FileNotFoundError, ValueError) as e:
        _fail(str(e), f_out)

    crystal = pyxtal()
    try:
        crystal.from_seed(seed=pmg_structure, tol=args.symprec, a_tol=args.angle_tolerance)
    except Exception as e:
        _fail(f"could not determine the symmetry of '{args.file}' -- {e}", f_out)

    return crystal, pmg_structure


def describe_candidate_group(dim, pmg_structure, crystal_group_label, symprec, angle_tolerance, f_out):
    """Returns (label, detected_number) for the symmetry group of one
    generated candidate, appropriate to `dim` (`detected_number` is None
    except for dim 3, where it's used for the requested-vs-detected check):
      - dim 3: re-detects the SPACE GROUP of the actual atomic positions via
        SpacegroupAnalyzer (not just pyxtal's own build-time label) -- this is
        what catches a random placement landing on a higher-symmetry special
        position, or (a real failure) a lower one than requested.
      - dim 0: re-detects the POINT GROUP via core_symmetry.point_group_label
        (pymatgen's non-periodic PointGroupAnalyzer, same convention as
        stb-molecule) -- like dim 3, this works on the actual final geometry,
        so it still reflects reality after an --ml-rank relax.
      - dim 2/1: no independent layer-/rod-group re-detector exists anywhere
        in this codebase (spglib's get_layergroup is only wrapped for the
        vacuum-heuristic case core/symmetry.py already covers, and it has no
        rod-group equivalent at all) -- reports pyxtal's own as-built
        `crystal_group_label` instead. Documented limitation: for dim 2/1,
        this label is NOT re-verified after --ml-rank the way dim 3/0 are.
    """
    if dim == 3:
        try:
            sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
            number = sga.get_space_group_number()
            return f"{sga.get_space_group_symbol()} (No. {number})", number
        except Exception as e:
            print_dual(color_text(f"[WARNING] Space group re-detection failed: {e}", 'yellow'), f_out)
            return "N/A", None
    if dim == 0:
        label = core_symmetry.point_group_label(pmg_structure)
        return (label if label else "N/A (could not be determined)"), None
    return crystal_group_label, None


def run_substitute(args, crystal, pmg_structure, f_out):
    """Core logic for --substitute: swaps elements in the already-loaded
    `crystal`/`pmg_structure` (see load_crystal_from_file, called by main()
    under [1] INPUT STRUCTURE, before this runs) while preserving the
    symmetry framework, and returns (pmg_after, substitutions) for main() to
    validate/report/write. Detail lines are print_dual'd here; section
    headers/structure-validation/symmetry-table/writing are all main()'s
    responsibility, shared with every other mode.
    """
    try:
        substitutions = parse_substitutions(args.substitute)
    except ValueError as e:
        _fail(str(e), f_out)

    present = {site.specie.symbol for site in pmg_structure}
    missing = [old for old in substitutions if old not in present]
    if missing:
        _fail(f"{', '.join(missing)} not present in '{args.file}' "
              f"(species present: {', '.join(sorted(present))}).", f_out)

    print_dual(f"Original formula : {pmg_structure.composition.reduced_formula}", f_out)
    for old, new in substitutions.items():
        print_dual(f"Substituting     : {old} -> {new}", f_out)

    crystal.substitute(substitutions)
    new_pmg = crystal.to_pymatgen()

    if len(new_pmg) != len(pmg_structure):
        print_dual(color_text(
            f"[NOTE] Atom count changed from {len(pmg_structure)} to {len(new_pmg)} -- "
            "this substitution made two previously-distinct species identical, so the "
            "structure collapsed to a smaller true primitive cell. Still a valid structure, "
            "just not a direct atom-for-atom relabeling of the original.", 'yellow'), f_out)

    print_dual(f"Output formula   : {new_pmg.composition.reduced_formula}", f_out)
    return new_pmg, substitutions


def run_group_transform(args, mode, crystal, f_out):
    """Core logic shared by --subgroup/--supergroup: runs pyxtal's search on
    the already-loaded `crystal` (see load_crystal_from_file, called by
    main() under [1] INPUT STRUCTURE, before this runs). The two modes only
    differ in which pyxtal method they call, how they unwrap its result,
    and their user-facing messages. Returns (candidates, found_label) for
    main() to validate/report/write -- extracted here once --supergroup was
    added alongside the pre-existing --subgroup, rather than keeping two
    near-identical copies in sync by hand.
    """
    from pyxtal.msg import Error as PyxtalError

    if args.count < 1:
        _fail("--count must be at least 1.", f_out)

    if mode == "subgroup":
        eps = args.eps if args.eps is not None else 0.05
        group_type = args.group_type if args.group_type is not None else "t"
        if not (eps > 0):
            _fail("--eps must be a positive, finite number.", f_out)
        print_dual(f"Perturbation --eps    : {eps}", f_out)
        print_dual(f"Group type            : {group_type}", f_out)
        if args.target_group is not None:
            print_dual(f"Target space group    : {args.target_group}", f_out)
        try:
            candidates = crystal.subgroup(H=args.target_group, eps=eps, group_type=group_type)
        except (PyxtalError, ValueError, RuntimeError) as e:
            _fail(str(e), f_out)
        no_candidates_msg = ("no subgroup candidates found for the given constraints -- try a "
                              "different --target-group, --group-type, or --eps.")
        found_label = "subgroup candidates"
    else:
        print_dual(f"Target space group    : {args.target_group}", f_out)
        d_tol = args.d_tol if args.d_tol is not None else 1.0
        if not (d_tol > 0):
            _fail("--d-tol must be a positive, finite number.", f_out)
        print_dual(f"Displacement --d-tol  : {d_tol}", f_out)
        try:
            result = crystal.supergroup(G=args.target_group, d_tol=d_tol)
        except (PyxtalError, ValueError, RuntimeError) as e:
            _fail(str(e), f_out)
        candidates = result[0] if isinstance(result, tuple) else result
        no_candidates_msg = ("no supergroup structure found towards the given --target-group -- "
                              "try a larger --d-tol or a different --target-group.")
        found_label = "supergroup candidates"

    if not candidates:
        _fail(no_candidates_msg, f_out)

    print_dual(f"{found_label.capitalize()} found : {len(candidates)}", f_out)
    return candidates, found_label


def run_ml_rank(pmg_structures, model_arg, model_desc, f_out):
    """Relaxes (positions only) each pmg_structure in `pmg_structures` --
    a list of (out_name, pmg_structure, is_isolated) -- with the requested
    MACE model, prints them ranked by relaxed energy, and returns a NEW list
    of (out_name, relaxed_pmg_structure, is_isolated) in the SAME order (not
    re-sorted -- ranking is informational, file naming/order is untouched).
    Mirrors stb-defect --ml-rank's exact convention (a fast, relative
    pre-screen, not an absolute DFT formation energy), reusing
    core/mace_relax.py. Writing the relaxed structure to disk is main()'s
    job (the single, later WRITING OUTPUT FILE(S) section), not this
    function's -- so the header comment written there can already include
    the relax outcome.

    `is_isolated` (True for --dim 0's vacuum-boxed cluster) disables periodic
    boundary conditions before relaxing -- ASE's pymatgen adapter otherwise
    always sets pbc=True, which would let MACE see spurious interactions
    with the cluster's own periodic images across the vacuum box.

    A relax failure on one candidate (e.g. a pathological geometry) is
    caught and skipped rather than aborting the whole ranking -- consistent
    with how the generation loop itself treats a single failed attempt.
    """
    print_dual(f"Model : {model_desc}", f_out)
    calc = mace_relax.get_calculator(model_arg)
    for line in mace_relax.describe_model(model_arg, calc):
        print_dual(line, f_out)

    updated = []
    rankings = []
    for out_name, pmg, is_isolated in pmg_structures:
        atoms = AseAtomsAdaptor.get_atoms(pmg)
        if is_isolated:
            atoms.pbc = False
        try:
            mace_relax.relax(atoms, calc, fmax=0.05, max_steps=200)
            energy = atoms.get_potential_energy()
        except Exception as e:
            print_dual(color_text(f"[WARNING] Relaxing {out_name} failed ({e}) -- kept "
                                   "unrelaxed, excluded from ranking.", 'yellow'), f_out)
            updated.append((out_name, pmg, is_isolated))
            continue
        relaxed_pmg = AseAtomsAdaptor.get_structure(atoms)
        updated.append((out_name, relaxed_pmg, is_isolated))
        rankings.append((out_name, energy))

    if not rankings:
        print_dual(color_text("[ERROR] Every candidate failed to relax -- no ranking to show.", 'red'), f_out)
        return updated

    rankings.sort(key=lambda r: r[1])
    e_min = rankings[0][1]
    name_width = max(24, max(len(out_name) for out_name, _ in rankings) + 2)
    print_dual(f"\n{color_text('ML-ranked structures (relaxed energy, most stable first):', 'bold')}", f_out)
    print_dual(f"  {'Rank':<5}{'File':<{name_width}}{'Energy (eV)':<16}{'dE (eV)':<10}", f_out)
    for rank, (out_name, energy) in enumerate(rankings, start=1):
        print_dual(f"  {rank:<5}{out_name:<{name_width}}{energy:<16.4f}{energy - e_min:<10.4f}", f_out)
    print_dual(color_text(
        "\n[NOTE] A relative comparison from a fast ML potential, not an absolute DFT "
        "formation energy -- use it to prioritize which candidate(s) to relax with SIESTA, "
        "not as a final answer.", 'yellow'), f_out)
    return updated


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Casts one or more random structures compatible with a given symmetry group.", 'bold')}
Give a symmetry group and a composition (species + how many of each) -- pyxtal
places the atoms on randomly chosen, symmetry-compatible Wyckoff positions
for you. This is the inverse of stb-crystalbuilder: use crystalbuilder when
you already know the exact Wyckoff sites you want, use crystalcast when you
want valid candidate structures generated for you (e.g. as starting guesses
for structure prediction). --molecular packs whole rigid molecules (instead
of bare atoms) into the symmetry group. --ml-rank pre-screens generated
candidates by MACE relaxed energy. --substitute/--subgroup/--supergroup
all instead operate on an existing structure read via -f/--file:
--substitute swaps elements while preserving the symmetry framework,
--subgroup/--supergroup search for related lower/higher-symmetry
structures. (To inspect an existing structure's own symmetry/Wyckoff
decomposition instead of generating or transforming one, use
stb-symmetry or stb-unitcell.)""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --group 225 --species Ni O --num-ions 4 8\n"
               "  %(prog)s --group Fd-3m --species Fe O --num-ions 8 16 \\\n"
               "      --count 5 --seed 42 -o spinel.fdf\n"
               "  %(prog)s --group 225 --species Ni O --num-ions 4 8 --count 5 \\\n"
               "      --ml-rank --model medium\n"
               "  %(prog)s --group 225 --species Na Cl --num-ions 4 4 \\\n"
               "      --lattice 5.64 5.64 5.64 90 90 90\n"
               "  %(prog)s --group 225 --species Na Cl --num-ions 4 4 --sites 4a 4b\n"
               "  %(prog)s --dim 2 --group 65 --species C --num-ions 6 --thickness 3.4\n"
               "  %(prog)s --dim 0 --group D3d --species C --num-ions 6 --vacuum 12\n"
               "  %(prog)s --molecular --group 19 --species H2O --num-ions 4\n"
               "  %(prog)s --molecular --group 14 --species aspirin --num-ions 4 -o aspirin.fdf\n"
               "  %(prog)s --list-molecules\n"
               "  %(prog)s --substitute Cl:F -f rocksalt.fdf -o rocksalt_f.fdf\n"
               "  %(prog)s --subgroup -f spinel_1.fdf --count 5 -o distorted.fdf\n"
               "  %(prog)s --supergroup --target-group 225 -f distorted_1.fdf -o parent.fdf\n"
               "  %(prog)s --group 225 --species Ni O --num-ions 4 8 --save-report --view\n"
    )

    parser.add_argument("--substitute", nargs="+", default=None, metavar="OLD:NEW",
                        help="Substitute elements in an existing structure (read via -f/--file), "
                             "preserving its symmetry framework: one or more OLD:NEW pairs, e.g. "
                             "--substitute Cl:F Na:K. Mutually exclusive with --subgroup/"
                             "--supergroup and the generation options below.")
    parser.add_argument("--subgroup", action="store_true",
                        help="Generate lower-symmetry subgroup variant(s) of an existing "
                             "structure (read via -f/--file) -- a symmetry-breaking distortion. "
                             "Up to --count candidates are written (numbered like --count > 1 "
                             "in generation mode). Optionally narrow the search with "
                             "--target-group/--group-type/--eps. Mutually exclusive with "
                             "--substitute/--supergroup and the generation options.")
    parser.add_argument("--supergroup", action="store_true",
                        help="Search for higher-symmetry supergroup structure(s) of an existing "
                             "structure (read via -f/--file) towards --target-group (required -- "
                             "the installed pyxtal doesn't support auto-searching all possible "
                             "supergroups). Mutually exclusive with --substitute/--subgroup and "
                             "the generation options below.")
    parser.add_argument("--target-group", type=int, default=None,
                        help="Target space group number. Optional filter for --subgroup (only "
                             "search towards this specific subgroup); required destination for "
                             "--supergroup.")
    parser.add_argument("--group-type", choices=["t", "k", "t+k"], default=None,
                        help="--subgroup only: 't' (translationengleiche), 'k' (klassengleiche), "
                             "or 't+k' relation to search (default: t).")
    parser.add_argument("--eps", type=float, default=None,
                        help="--subgroup only: perturbation applied to atomic coordinates when "
                             "breaking symmetry (default: 0.05).")
    parser.add_argument("--d-tol", type=float, default=None,
                        help="--supergroup only: maximum atomic-displacement tolerance allowed "
                             "when searching for a higher-symmetry parent (default: 1.0).")
    parser.add_argument("-f", "--file", type=str, default=None,
                        help="Input .fdf structure file. Required with --substitute/--subgroup/"
                             "--supergroup, unused otherwise.")

    parser.add_argument("--dim", type=int, choices=[3, 2, 1, 0], default=3,
                        help="Structure dimensionality: 3 = bulk (space group, default), "
                             "2 = layer (layer group, periodic in a/b with vacuum along c), "
                             "1 = rod/wire (rod group, periodic along c with vacuum in a/b), "
                             "0 = isolated cluster (point group, no periodicity). "
                             "--molecular does not support --dim 0 (upstream pyxtal limitation). "
                             "Generation mode only.")
    parser.add_argument("--group", type=str, default=None,
                        help="Symmetry group identifier, meaning depends on --dim: space group "
                             "number/symbol (dim 3, e.g. 225 or 'Fm-3m'), layer group number "
                             "1-80 (dim 2), rod group number 1-75 (dim 1), or point group "
                             "number/Schoenflies symbol 1-32 (dim 0, e.g. 20 or 'D3d'). "
                             "Required in generation mode.")
    parser.add_argument("--molecular", action="store_true",
                        help="Treat each --species entry as a whole rigid molecule instead of "
                             "a bare element: a name from pyxtal's bundled collection (see "
                             "--list-molecules), a path to a .xyz/.gjf/.g03/.json file, or a "
                             "SMILES string as 'SMILES.smi' (needs the optional 'rdkit' package).")
    parser.add_argument("--list-molecules", action="store_true",
                        help="Print the molecule names bundled with pyxtal (usable directly as "
                             "--species with --molecular) and exit.")
    parser.add_argument("--species", nargs="+", default=None,
                        help="Element symbols (or, with --molecular, molecule identifiers) in "
                             "the structure, e.g. --species Ni O. Required in generation mode.")
    parser.add_argument("--num-ions", nargs="+", type=int, default=None,
                        help="Number of atoms (or, with --molecular, molecules) of each "
                             "--species, same order and count, e.g. --num-ions 4 8. Required "
                             "in generation mode.")
    parser.add_argument("--volume-factor", type=float, default=1.1,
                        help="Scales the estimated cell volume before placing atoms; raise "
                             "it if generation keeps failing to fit atoms without overlap "
                             "(default: 1.1, pyxtal's own default).")
    parser.add_argument("--thickness", type=float, default=None,
                        help="Layer thickness in Ang (--dim 2 only; default: chosen "
                             "automatically by pyxtal).")
    parser.add_argument("--area", type=float, default=None,
                        help="Rod cross-sectional area in Ang^2 (--dim 1 only; default: "
                             "chosen automatically by pyxtal).")
    parser.add_argument("--vacuum", type=float, default=None,
                        help="Vacuum padding in Ang around the cluster, on every side "
                             "(--dim 0 only, matches stb-molecule's --vacuum; default: 10.0).")
    parser.add_argument("--lattice", type=float, nargs=6, default=None,
                        metavar=("A", "B", "C", "ALPHA", "BETA", "GAMMA"),
                        help="Fix the cell instead of estimating it from --volume-factor, e.g. "
                             "for matching a known experimental cell: A B C ALPHA BETA GAMMA "
                             "(Ang/degrees). Not valid with --dim 0 (use --vacuum instead). For "
                             "--dim 2/1, only the periodic direction(s) are actually honored -- "
                             "the vacuum direction(s) are still sized by --thickness/--area/"
                             "pyxtal's own default, same as without --lattice.")
    parser.add_argument("--sites", nargs="+", default=None,
                        help="Pre-assign Wyckoff positions instead of leaving every site random: "
                             "one entry per --species, each a comma-separated list of Wyckoff "
                             "labels summing to that species' --num-ions, e.g. --species Na O "
                             "--num-ions 4 8 --sites 4a 4b,4c. Default: fully random assignment.")
    parser.add_argument("--max-attempts", type=int, default=10,
                        help="Internal retries pyxtal allows itself per structure to find a "
                             "non-overlapping placement before giving up on it (default: 10).")
    parser.add_argument("--count", type=int, default=1,
                        help="Number of independent random structures to generate in generation "
                             "mode, or number of candidates to keep in --subgroup/--supergroup "
                             "mode (default: 1).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed. Same seed + same inputs reproduces the same batch "
                             "of structures (default: not seeded, different every run). With "
                             "--molecular, only the lattice is reproduced this way -- molecule "
                             "orientations are not (an upstream pyxtal limitation). Unused "
                             "outside generation mode.")

    parser.add_argument("--ml-rank", action="store_true",
                        help="Generation mode only. After generating --count structures, "
                             "quickly relaxes each one's positions with a MACE potential, "
                             "prints them ranked by relaxed energy, and writes each output file "
                             "with its relaxed geometry. A relative comparison from a fast ML "
                             "potential, not an absolute DFT formation energy.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-rank (default: small). "
                             "Only valid together with --ml-rank.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-rank, instead "
                             "of a MACE-MP-0 foundation size. Only valid together with --ml-rank.")

    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision: for generation, used to re-detect each --dim 3 "
                             "candidate's actual space group; for --substitute/--subgroup/"
                             "--supergroup, the tolerance passed to pyxtal's own symmetry "
                             "detection (default: 1e-3, matches stb-symmetry/stb-unitcell).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Angle tolerance in degrees for symmetry detection (default: 5.0, "
                             "matches stb-symmetry/stb-unitcell).")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) of every structure written "
                             "by this run (for --substitute/--subgroup/--supergroup, the input "
                             "structure is shown first, as frame 0) -- page through frames in "
                             "ase-gui. Needs a display. Off by default.")

    parser.add_argument("-o", "--output", type=str, default="crystalcast.fdf",
                        help="Output .fdf file name (default: crystalcast.fdf). With "
                             "--count > 1 (generation, --subgroup, or --supergroup), each "
                             "structure is written as '<output>_<N>.fdf'; the tool refuses to "
                             "run if any of those numbered names already exist, to avoid "
                             "silently overwriting unrelated files.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-crystalcast {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.list_molecules:
        require_pyxtal()
        from pyxtal.database.collection import Collection
        names = list(Collection("molecules"))
        print(color_text(f"Molecules bundled with pyxtal ({len(names)}):", 'bold'))
        for name in names:
            print(f"  {name}")
        sys.exit(0)

    mode_flags = {"substitute": bool(args.substitute), "subgroup": args.subgroup, "supergroup": args.supergroup}
    active_modes = [name for name, on in mode_flags.items() if on]
    if len(active_modes) > 1:
        parser.error(f"--{' and --'.join(active_modes)} are mutually exclusive.")
    mode = active_modes[0] if active_modes else "generate"

    if mode == "generate":
        if args.molecular and args.dim == 0:
            parser.error("--molecular does not support --dim 0 (upstream pyxtal limitation). "
                         "Use --dim 3, 2, or 1.")
        if not args.group:
            parser.error("--group is required in generation mode.")
        if not args.species:
            parser.error("--species is required in generation mode.")
        if not args.num_ions:
            parser.error("--num-ions is required in generation mode.")
        if args.lattice is not None and args.dim == 0:
            parser.error("--lattice is not valid with --dim 0 -- use --vacuum instead.")
    else:
        if args.ml_rank:
            parser.error("--ml-rank is only valid in generation mode.")
        if mode == "supergroup" and args.target_group is None:
            parser.error("--supergroup requires --target-group (the installed pyxtal doesn't "
                         "support auto-searching all possible supergroups).")

    if (args.custom_model or args.model != "small") and not args.ml_rank:
        parser.error("--model/--custom-model are only valid together with --ml-rank.")

    require_pyxtal()
    if args.ml_rank:
        require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-CRYSTALCAST REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"
    model_arg = args.custom_model if args.custom_model else args.model

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Mode           : {mode}", f_out)
    if mode == "generate":
        print_dual(f"Dimension      : {args.dim} ({DIM_LABELS[args.dim]})", f_out)
        print_dual(f"Requested group: {args.group}", f_out)
        print_dual(f"Molecular      : {'yes' if args.molecular else 'no'}", f_out)
        print_dual(f"Structures req.: {args.count}", f_out)
        print_dual(f"Seed           : {args.seed if args.seed is not None else 'none (not reproducible)'}", f_out)
        print_dual(f"ML rank        : {model_desc if args.ml_rank else 'no'}", f_out)
    else:
        print_dual(f"Input file     : {args.file}", f_out)
        if mode == "subgroup":
            print_dual(f"Target group   : {args.target_group if args.target_group is not None else 'auto-search'}", f_out)
        elif mode == "supergroup":
            print_dual(f"Target group   : {args.target_group}", f_out)
        else:
            print_dual(f"Substitutions  : {', '.join(args.substitute)}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)

    for flag, value, valid_modes in (
        ("--eps", args.eps, ("subgroup",)),
        ("--group-type", args.group_type, ("subgroup",)),
        ("--d-tol", args.d_tol, ("supergroup",)),
        ("--target-group", args.target_group, ("subgroup", "supergroup")),
    ):
        if value is not None and mode not in valid_modes:
            print_dual(color_text(
                f"[NOTE] {flag} is ignored (only used with --{'/--'.join(valid_modes)}).", 'yellow'), f_out)

    # ------------------------------------------------------------------
    # --substitute / --subgroup / --supergroup: always a single 3D
    # -periodic input structure, so they share one report shape.
    # ------------------------------------------------------------------
    if mode in ("substitute", "subgroup", "supergroup"):
        print_section("[1] INPUT STRUCTURE", f_out)
        crystal, pmg_before = load_crystal_from_file(args, f_out)
        print_dual(f"Original formula     : {pmg_before.composition.reduced_formula}", f_out)
        print_dual(f"Atoms                : {len(pmg_before)}", f_out)
        print_dual(f"Original space group : {crystal.group.symbol} (No. {crystal.group.number})", f_out)

        print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
        try:
            frac_before = [site.frac_coords for site in pmg_before]
            vacuum_before = kspace.detect_vacuum_axes(frac_before, pmg_before.lattice.matrix, VACUUM_GAP_ANG)
            structure_checks.run_malformation_checks(pmg_before, vacuum_before, f_out)
        except Exception as e:
            print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

        section_title = {"substitute": "SUBSTITUTION", "subgroup": "SUBGROUP SEARCH",
                         "supergroup": "SUPERGROUP SEARCH"}[mode]
        print_section(f"[3] {section_title}", f_out)
        if mode == "substitute":
            pmg_after, substitutions = run_substitute(args, crystal, pmg_before, f_out)
            candidates_pmg = [pmg_after]
            # Substitution preserves the symmetry framework by construction --
            # crystal.group is unchanged, so it's the same authoritative label
            # for every (single) candidate here as pyxtal's own subgroup()/
            # supergroup() results are for their candidates below.
            candidate_group_labels = [f"{crystal.group.symbol} (No. {crystal.group.number})"]
        else:
            candidates, found_label = run_group_transform(args, mode, crystal, f_out)
            candidates_pmg = [c.to_pymatgen() for c in candidates]
            # pyxtal's own per-candidate group (what the search actually built
            # towards), NOT a pymatgen re-detection -- verified live that
            # pymatgen's SpacegroupAnalyzer can report an inconsistent
            # space-group/point-group pair for these slightly-distorted
            # candidates (e.g. still "Fm-3m (225)" alongside point group
            # "4/mmm", which is not a real m-3m subgroup pairing); pyxtal's
            # own label is authoritative since it's what subgroup()/
            # supergroup() actually constructed.
            candidate_group_labels = [f"{c.group.symbol} (No. {c.group.number})" for c in candidates]
        print_dual(f"{len(candidates_pmg)} candidate(s) produced.", f_out)
        if mode != "substitute":
            for i, label in enumerate(candidate_group_labels, start=1):
                print_dual(f"  #{i}: {label}", f_out)

        try:
            names = resolve_output_paths(args.output, args.count if mode != "substitute" else 1)
        except ValueError as e:
            _fail(str(e), f_out)

        if mode != "substitute" and len(candidates_pmg) > len(names):
            print_dual(color_text(
                f"[NOTE] {len(candidates_pmg) - len(names)} more {found_label} found but not "
                "written -- raise --count to keep more of them.", 'yellow'), f_out)
            candidates_pmg = candidates_pmg[:len(names)]
            candidate_group_labels = candidate_group_labels[:len(names)]

        print_section("[4] STRUCTURE VALIDATION (post-transform)", f_out)
        for i, pmg in enumerate(candidates_pmg, start=1):
            if len(candidates_pmg) > 1:
                print_dual(f"-- Candidate #{i} ({names[i - 1]}) --", f_out)
            try:
                frac_after = [site.frac_coords for site in pmg]
                vacuum_after = kspace.detect_vacuum_axes(frac_after, pmg.lattice.matrix, VACUUM_GAP_ANG)
                structure_checks.run_malformation_checks(pmg, vacuum_after, f_out)
            except Exception as e:
                print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

        print_section("[5] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
        before_info = core_symmetry.symmetry_summary(pmg_before, args.symprec, VACUUM_GAP_ANG)
        if len(candidates_pmg) == 1:
            after_info = core_symmetry.symmetry_summary(candidates_pmg[0], args.symprec, VACUUM_GAP_ANG)
            if "Error" in before_info or "Error" in after_info:
                print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
            else:
                properties = ["Crystal System", "Space Group", "Point Group", "Hall Symbol"]
                after_values = dict(after_info)
                # pyxtal's own label (see note above), reformatted to match
                # symmetry_summary's own "SYMBOL (number)" convention (no "No.")
                # so the Before/After table reads consistently.
                after_values["Space Group"] = candidate_group_labels[0].replace("(No. ", "(")
                rows = [([prop, str(before_info.get(prop, "N/A")), str(after_values.get(prop, "N/A"))], None)
                        for prop in properties]
                print_table(["Property", "Before", "After"], rows, f_out)
        else:
            print_dual(f"Before: Space Group {before_info.get('Space Group', 'N/A')}, "
                       f"Point Group {before_info.get('Point Group', 'N/A')}", f_out)
            for i, (pmg, glabel) in enumerate(zip(candidates_pmg, candidate_group_labels), start=1):
                print_dual(f"  #{i} ({names[i - 1]}): Space Group {glabel}", f_out)

        print_section("[6] WRITING OUTPUT FILE(S)", f_out)
        written = []
        for i, (pmg, out_name, glabel) in enumerate(zip(candidates_pmg, names, candidate_group_labels), start=1):
            if mode == "substitute":
                header_comment = [
                    "Structure produced by stb-crystalcast --substitute.",
                    f"Input file: {args.file}.",
                    f"Substitutions applied: {', '.join(f'{o}->{n}' for o, n in substitutions.items())}.",
                    f"Formula: {pmg_before.composition.reduced_formula} -> {pmg.composition.reduced_formula}.",
                    f"Space group: {glabel}.",
                ]
            else:
                header_comment = [
                    f"Structure produced by stb-crystalcast --{mode}.",
                    f"Input file: {args.file} (original space group {before_info.get('Space Group', 'N/A')}).",
                    f"Target group: {args.target_group if args.target_group is not None else 'auto-search'}.",
                    f"Detected space group: {glabel}.",
                ]
            new_structure = structure_to_fdf(pmg)
            structure_io.write_fdf(new_structure, out_name, header_comment=header_comment)
            written.append(out_name)
            print_dual(color_text(f"[OK] #{i}: {pmg.composition.reduced_formula} "
                                  f"({len(pmg)} atoms) written to '{out_name}'.", 'green'), f_out)

        if len(written) < len(names):
            print_dual(color_text(
                f"[WARNING] Partial success: {len(written)} of {len(names)} requested "
                "structure(s) written.", 'yellow'), f_out)

        print_section("[7] REFERENCES", f_out)
        bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.PYXTAL]
        citations.write_bib_file(BIB_FILE, bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)

        print_section("[8] SUMMARY & FILES", f_out)
        print_dual(f"Status         : {'OK' if len(written) == len(names) else 'PARTIAL'}", f_out)
        print_dual(f"Structures     : {len(written)} of {len(names)} written", f_out)
        print_dual(f"Input file     : {args.file}", f_out)
        print_dual(f"Output file(s) : {', '.join(written)}", f_out)
        print_dual(f"References     : {BIB_FILE}", f_out)
        if report_path:
            print_dual(f"Report         : {report_path}", f_out)

        if f_out:
            f_out.close()

        if args.view:
            atoms_list = [AseAtomsAdaptor.get_atoms(pmg_before)] + \
                [AseAtomsAdaptor.get_atoms(pmg) for pmg in candidates_pmg]
            view_structure_interactive(atoms_list)

        if len(written) < len(names):
            sys.exit(1)
        return

    # ------------------------------------------------------------------
    # generate (default mode)
    # ------------------------------------------------------------------
    from pyxtal import pyxtal
    from pyxtal.msg import Error as PyxtalError

    print_section("[1] BUILD PARAMETERS", f_out)

    if len(args.species) != len(args.num_ions):
        _fail(f"--species has {len(args.species)} entries but --num-ions has "
              f"{len(args.num_ions)} -- they must match one-to-one.", f_out)

    if args.molecular:
        try:
            validate_molecule_species(args.species)
        except ValueError as e:
            _fail(str(e), f_out)
    else:
        for symbol in args.species:
            try:
                Element(symbol)
            except ValueError as e:
                _fail(str(e), f_out)

    if args.count < 1:
        _fail("--count must be at least 1.", f_out)

    group = resolve_group(args.group)
    requested_number = None
    if args.dim == 3:
        try:
            requested_number = SpaceGroup(str(group)).int_number if isinstance(group, str) \
                else SpaceGroup.from_int_number(group).int_number
        except ValueError:
            requested_number = None

    species_label = "Molecule" if args.molecular else "Species"
    print_dual(f"Requested {DIM_LABELS[args.dim]} : {args.group}", f_out)
    rows = [([symbol, str(n)], None) for symbol, n in zip(args.species, args.num_ions)]
    print_table([species_label, "Count"], rows, f_out)

    if args.molecular and args.seed is not None:
        print_dual(color_text(
            "[NOTE] --seed only reproduces the lattice for --molecular structures -- "
            "molecule orientations still vary between runs (an upstream pyxtal limitation).",
            'yellow'), f_out)
    for flag, value, required_dim in (
        ("--thickness", args.thickness, 2), ("--area", args.area, 1), ("--vacuum", args.vacuum, 0)
    ):
        if value is not None and args.dim != required_dim:
            print_dual(color_text(f"[NOTE] {flag} is ignored (only used with --dim {required_dim}).", 'yellow'), f_out)
    vacuum = args.vacuum if args.vacuum is not None else 10.0

    fixed_lattice = None
    if args.lattice is not None:
        from pyxtal.symmetry import Group
        from pyxtal.lattice import Lattice as PyxtalLattice
        ltype = Group(group, dim=args.dim).lattice_type
        mismatch = check_lattice_type(*args.lattice, ltype=ltype)
        if mismatch is not None:
            _fail(mismatch, f_out)
        try:
            fixed_lattice = PyxtalLattice.from_para(*args.lattice, ltype=ltype)
        except ValueError as e:
            _fail(str(e), f_out)
        print_dual(f"Fixed lattice  : a={args.lattice[0]} b={args.lattice[1]} "
                   f"c={args.lattice[2]} alpha={args.lattice[3]} beta={args.lattice[4]} "
                   f"gamma={args.lattice[5]} ({ltype})", f_out)

    sites = None
    if args.sites is not None:
        if len(args.sites) != len(args.species):
            _fail(f"--sites has {len(args.sites)} entries but --species has "
                  f"{len(args.species)} -- one --sites entry is needed per --species.", f_out)
        sites = [[wp.strip() for wp in entry.split(",")] for entry in args.sites]

        from pyxtal.symmetry import Group as SitesGroup
        try:
            site_group = SitesGroup(group, dim=args.dim)
            all_zero_dof = True
            for symbol, wp_letters in zip(args.species, sites):
                print_dual(f"{species_label} sites : {symbol} -> {', '.join(wp_letters)}", f_out)
                for letter in wp_letters:
                    if site_group.get_wyckoff_position(letter).get_dof() > 0:
                        all_zero_dof = False
        except (IndexError, ValueError) as e:
            _fail(f"--sites contains a Wyckoff label not valid for this group -- {e}", f_out)

        if args.count > 1 and all_zero_dof:
            print_dual(color_text(
                f"[NOTE] Every assigned Wyckoff site has zero free parameters -- all "
                f"{args.count} requested structures will be identical (--sites leaves "
                "nothing left to randomize; --seed won't change that either).", 'yellow'), f_out)

    rng = np.random.default_rng(args.seed)

    try:
        names = resolve_output_paths(args.output, args.count)
    except ValueError as e:
        _fail(str(e), f_out)

    print_section("[2] GENERATION", f_out)
    results = []  # (out_name, pmg_structure, is_isolated, group_label)
    for i in range(1, args.count + 1):
        crystal = pyxtal(molecular=args.molecular)
        try:
            crystal.from_random(
                dim=args.dim,
                group=group,
                species=args.species,
                numIons=args.num_ions,
                factor=args.volume_factor,
                thickness=args.thickness if args.dim == 2 else None,
                area=args.area if args.dim == 1 else None,
                lattice=fixed_lattice,
                sites=sites,
                max_count=args.max_attempts,
                random_state=rng,
            )
        except (RuntimeError, PyxtalError, ValueError) as e:
            # Not all attempt-independent: pyxtal can fail here for a hard, deterministic
            # reason (invalid --group, incompatible composition -- fails identically every
            # time) or for a per-attempt reason tied to this random draw (placement/volume
            # retries exhausted, or -- for --molecular -- an orientation-dependent
            # conformer/connectivity error). The two aren't reliably distinguishable by
            # exception type alone, so every failure is treated the same way here: skip
            # this attempt and keep whatever was already generated, rather than discarding
            # earlier successes in the same --count batch.
            print_dual(color_text(f"[WARNING] Attempt #{i}: {e} -- skipped.", 'yellow'), f_out)
            continue

        if args.dim == 0:
            molecule = crystal.to_pymatgen()
            pmg_structure = molecule_to_boxed_structure(molecule, vacuum)
            coord_format = "cartesian"
            overlap_check_target = molecule
        else:
            pmg_structure = crystal.to_pymatgen()
            coord_format = "fractional"
            overlap_check_target = pmg_structure

        # Checked on overlap_check_target (the raw dim=0 Molecule, not the boxed
        # Structure) because the box's own periodicity would otherwise make a real
        # atom look "close" to its own periodic image across the artificial vacuum cell.
        min_dist = structure_io.min_pairwise_distance(overlap_check_target)
        if min_dist is not None and min_dist < 0.5:
            print_dual(color_text(
                f"[WARNING] Structure #{i} has atoms unusually close together "
                f"({min_dist:.3f} Ang) -- check --volume-factor.", 'yellow'), f_out)

        out_name = names[i - 1]
        results.append((out_name, pmg_structure, args.dim == 0, f"{crystal.group.symbol} (No. {crystal.group.number})"))
        print_dual(f"  -> #{i}: {pmg_structure.composition.reduced_formula} "
                   f"({len(pmg_structure)} atoms), coord_format={coord_format}, "
                   f"as-built {DIM_LABELS[args.dim]} {crystal.group.symbol} (No. {crystal.group.number})", f_out)

    if not results:
        _fail("none of the requested structures could be generated -- see the "
              "per-attempt messages above.", f_out)

    ml_rank_used = args.ml_rank
    if args.ml_rank:
        print_section("[3] ML RANKING (MACE)", f_out)
        ranked = run_ml_rank([(n, p, iso) for n, p, iso, _ in results], model_arg, model_desc, f_out)
        results = [(n, p, iso, glabel) for (n, p, iso), (_, _, _, glabel)
                   in zip(ranked, results)]

    print_section("[4] STRUCTURE VALIDATION", f_out)
    for i, (out_name, pmg, is_isolated, _) in enumerate(results, start=1):
        if len(results) > 1:
            print_dual(f"-- Candidate #{i} ({out_name}) --", f_out)
        try:
            frac_coords = [site.frac_coords for site in pmg]
            vacuum_axes = kspace.detect_vacuum_axes(frac_coords, pmg.lattice.matrix, VACUUM_GAP_ANG)
            structure_checks.run_malformation_checks(pmg, vacuum_axes, f_out)
        except Exception as e:
            print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[5] SYMMETRY ANALYSIS", f_out)
    print_dual(f"Requested {DIM_LABELS[args.dim]}: {args.group}", f_out)
    detected_labels = []
    for i, (out_name, pmg, is_isolated, group_label) in enumerate(results, start=1):
        detected, detected_number = describe_candidate_group(
            args.dim, pmg, group_label, args.symprec, args.angle_tolerance, f_out)
        detected_labels.append(detected)
        note = ""
        if args.dim == 3 and requested_number is not None and detected_number is not None \
                and detected_number != requested_number:
            if detected_number > requested_number:
                note = color_text(" [NOTE] higher-symmetry special position", 'yellow')
            else:
                note = color_text(" [WARNING] LOWER symmetry than requested", 'yellow')
        print_dual(f"  #{i} ({out_name}): {detected}{note}", f_out)

    print_section("[6] WRITING OUTPUT FILE(S)", f_out)
    written = []
    for i, ((out_name, pmg, is_isolated, group_label), detected) in enumerate(zip(results, detected_labels), start=1):
        species_order = None if args.molecular else args.species
        coord_format = "cartesian" if args.dim == 0 else "fractional"
        new_structure = structure_to_fdf(pmg, species_order=species_order, coord_format=coord_format)
        header_comment = [
            f"Structure cast by stb-crystalcast (mode: generate) from a random placement "
            f"in {DIM_LABELS[args.dim]} {args.group} (dim={args.dim}).",
            "Composition: " + ", ".join(f"{s} x {n}" for s, n in zip(args.species, args.num_ions)) + ".",
            f"Detected {DIM_LABELS[args.dim]}: {detected}.",
        ]
        if ml_rank_used:
            header_comment.append(f"ML-ranked with {model_desc} (positions-only relax).")
        structure_io.write_fdf(new_structure, out_name, header_comment=header_comment)
        written.append(out_name)
        print_dual(color_text(f"[OK] #{i}: {pmg.composition.reduced_formula} "
                              f"({len(pmg)} atoms) written to '{out_name}'.", 'green'), f_out)

    if len(written) < args.count:
        print_dual(color_text(
            f"[WARNING] Partial success: {len(written)} of {args.count} structure(s) written "
            f"({args.count - len(written)} skipped -- see warnings above).", 'yellow'), f_out)

    print_section("[7] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.PYXTAL]
    if ml_rank_used:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[8] SUMMARY & FILES", f_out)
    print_dual(f"Status         : {'OK' if len(written) == args.count else 'PARTIAL'}", f_out)
    print_dual(f"Structures     : {len(written)} of {args.count} written", f_out)
    print_dual(f"Output file(s) : {', '.join(written)}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    if args.view:
        atoms_list = [AseAtomsAdaptor.get_atoms(pmg) for _, pmg, _, _ in results]
        view_structure_interactive(atoms_list)

    if len(written) < args.count:
        sys.exit(1)


if __name__ == "__main__":
    main()
