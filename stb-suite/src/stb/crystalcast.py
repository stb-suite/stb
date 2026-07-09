#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import os
import sys
import math
import difflib
import argparse
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.core.periodic_table import Element
from stb.core import structure_io
from stb.core.cli import color_text, show_intro
from stb.core.deps import require_pyxtal, require_mace

MOLECULE_FILE_EXTENSIONS = ("xyz", "gjf", "g03", "json")

DIM_LABELS = {
    3: "space group",
    2: "layer group",
    1: "rod group",
    0: "point group",
}


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


def load_crystal_from_file(args):
    """Reads -f/--file and builds a pyxtal object via from_seed(). Shared by
    --analyze/--substitute/--subgroup/--supergroup -- all four operate on an
    existing structure instead of generating one from scratch.
    """
    from pyxtal import pyxtal

    if not args.file:
        print(color_text("Error: this mode requires -f/--file.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.file):
        print(color_text(f"Error: File '{args.file}' not found.", 'red'))
        sys.exit(1)

    print(f"  {color_text('Input file:', 'cyan')} {args.file}")

    try:
        structure = structure_io.read_fdf(args.file)
        pmg_structure = structure_io.to_pymatgen(structure)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    crystal = pyxtal()
    try:
        crystal.from_seed(seed=pmg_structure, tol=args.symprec, a_tol=args.angle_tolerance)
    except Exception as e:
        print(color_text(f"Error: could not determine the symmetry of '{args.file}' -- {e}", 'red'))
        sys.exit(1)

    return crystal, pmg_structure


def run_analyze(args):
    crystal, pmg_structure = load_crystal_from_file(args)

    print(f"\n  {color_text('Space group:', 'cyan')} {crystal.group.symbol} (No. {crystal.group.number})")
    print(f"  {color_text('Formula:', 'cyan')} {pmg_structure.composition.reduced_formula}")
    print(f"  {color_text('Symmetrically distinct sites:', 'cyan')} {len(crystal.atom_sites)}")

    print(f"\n{color_text('Wyckoff sites (paste each --site straight into stb-crystalbuilder):', 'bold')}")
    for site in crystal.atom_sites:
        x, y, z = site.position
        print(f"  --site {site.specie:<3} {x:.6f} {y:.6f} {z:.6f}"
              f"   {color_text(f'(Wyckoff {site.wp.multiplicity}{site.wp.letter})', 'yellow')}")


def run_substitute(args):
    try:
        substitutions = parse_substitutions(args.substitute)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    crystal, pmg_structure = load_crystal_from_file(args)

    present = {site.specie.symbol for site in pmg_structure}
    missing = [old for old in substitutions if old not in present]
    if missing:
        print(color_text(
            f"Error: {', '.join(missing)} not present in '{args.file}' "
            f"(species present: {', '.join(sorted(present))}).", 'red'))
        sys.exit(1)

    print(f"  {color_text('Original formula:', 'cyan')} {pmg_structure.composition.reduced_formula}")
    for old, new in substitutions.items():
        print(f"  {color_text('Substituting:', 'cyan')} {old} -> {new}")

    crystal.substitute(substitutions)
    new_pmg = crystal.to_pymatgen()

    if len(new_pmg) != len(pmg_structure):
        print(color_text(
            f"  Warning: atom count changed from {len(pmg_structure)} to {len(new_pmg)} -- "
            "this substitution made two previously-distinct species identical, so the "
            "structure collapsed to a smaller true primitive cell. Still a valid structure, "
            "just not a direct atom-for-atom relabeling of the original.", 'yellow'))

    new_structure = structure_to_fdf(new_pmg)
    structure_io.write_fdf(new_structure, args.output)

    print(f"\n  {color_text('Output formula:', 'cyan')} {new_pmg.composition.reduced_formula}")
    print(f"{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


def write_candidates(candidates, names, label="candidates"):
    """Writes up to len(names) pyxtal candidate structures (from --subgroup/
    --supergroup) to their corresponding output paths. Shared by both since
    they only differ in how `candidates` was produced.
    """
    written = []
    for i, (candidate, out_name) in enumerate(zip(candidates, names), start=1):
        pmg = candidate.to_pymatgen()
        new_structure = structure_to_fdf(pmg)
        structure_io.write_fdf(new_structure, out_name)
        written.append(out_name)
        print(f"  {color_text('->', 'green')} #{i}: {pmg.composition.reduced_formula} "
              f"({len(pmg)} atoms), space group {candidate.group.symbol} "
              f"(No. {candidate.group.number}): {out_name}")
    if len(candidates) > len(names):
        print(color_text(
            f"  ({len(candidates) - len(names)} more {label} found but not written -- "
            "raise --count to keep more of them.)", 'yellow'))
    return written


def run_group_transform(args, mode):
    """Shared body of --subgroup/--supergroup: load the structure, resolve
    output paths, run pyxtal's search, write the results, and report
    success/partial-success/failure. The two modes only differ in which
    pyxtal method they call, how they unwrap its result, and their
    user-facing messages -- extracted here once --supergroup was added
    alongside the pre-existing --subgroup, rather than keeping two
    near-identical copies of this flow in sync by hand.
    """
    from pyxtal.msg import Error as PyxtalError

    crystal, _ = load_crystal_from_file(args)
    print(f"  {color_text('Original space group:', 'cyan')} {crystal.group.symbol} (No. {crystal.group.number})")

    if args.count < 1:
        print(color_text("Error: --count must be at least 1.", 'red'))
        sys.exit(1)

    try:
        names = resolve_output_paths(args.output, args.count)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    if mode == "subgroup":
        eps = args.eps if args.eps is not None else 0.05
        group_type = args.group_type if args.group_type is not None else "t"
        if not (eps > 0):
            print(color_text("Error: --eps must be a positive, finite number.", 'red'))
            sys.exit(1)
        try:
            candidates = crystal.subgroup(H=args.target_group, eps=eps, group_type=group_type)
        except (PyxtalError, ValueError, RuntimeError) as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)
        no_candidates_msg = ("no subgroup candidates found for the given constraints -- try a "
                              "different --target-group, --group-type, or --eps.")
        found_label = "subgroup candidates"
    else:
        print(f"  {color_text('Target space group:', 'cyan')} {args.target_group}")
        d_tol = args.d_tol if args.d_tol is not None else 1.0
        if not (d_tol > 0):
            print(color_text("Error: --d-tol must be a positive, finite number.", 'red'))
            sys.exit(1)
        try:
            result = crystal.supergroup(G=args.target_group, d_tol=d_tol)
        except (PyxtalError, ValueError, RuntimeError) as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)
        candidates = result[0] if isinstance(result, tuple) else result
        no_candidates_msg = ("no supergroup structure found towards the given --target-group -- "
                              "try a larger --d-tol or a different --target-group.")
        found_label = "supergroup candidates"

    if not candidates:
        print(color_text(f"Error: {no_candidates_msg}", 'red'))
        sys.exit(1)

    print(f"  {color_text(found_label.capitalize() + ' found:', 'cyan')} {len(candidates)}")
    print()
    written = write_candidates(candidates, names, label=found_label)

    if len(written) < len(names):
        print(color_text(
            f"\nPartial success: {len(written)} of {len(names)} requested structure(s) written "
            f"(only {len(candidates)} candidate(s) were found).", 'yellow'))
        sys.exit(1)

    print(f"\n{color_text('Success:', 'green')} {len(written)} structure(s) written.")


def run_subgroup(args):
    run_group_transform(args, "subgroup")


def run_supergroup(args):
    run_group_transform(args, "supergroup")


def run_ml_rank(results):
    """Relaxes (positions only) each (out_name, pmg_structure, is_isolated)
    in `results` with MACE-MP-0, overwrites each file with its relaxed
    geometry, and prints them ranked by relaxed energy. Mirrors stb-defect
    --ml-rank's exact convention (same rationale: a fast, relative
    pre-screen, not an absolute DFT formation energy), reusing
    core/mace_relax.py.

    `is_isolated` (True for --dim 0's vacuum-boxed cluster) disables periodic
    boundary conditions before relaxing -- ASE's pymatgen adapter otherwise
    always sets pbc=True, which would let MACE see spurious interactions
    with the cluster's own periodic images across the vacuum box.

    A relax failure on one candidate (e.g. a pathological geometry) is
    caught and skipped rather than aborting the whole ranking -- consistent
    with how the generation loop itself treats a single failed attempt.
    """
    from stb.core import mace_relax
    from pymatgen.io.ase import AseAtomsAdaptor

    print(f"\n{color_text('ML ranking:', 'cyan')} relaxing each structure with MACE-MP-0 (positions only)...")
    calc = mace_relax.get_calculator()

    rankings = []
    for out_name, pmg, is_isolated in results:
        atoms = AseAtomsAdaptor.get_atoms(pmg)
        if is_isolated:
            atoms.pbc = False
        try:
            mace_relax.relax(atoms, calc, fmax=0.05, max_steps=200)
            energy = atoms.get_potential_energy()
        except Exception as e:
            print(color_text(f"  Warning: relaxing {out_name} failed ({e}) -- kept unrelaxed, "
                              "excluded from ranking.", 'yellow'))
            continue
        relaxed_pmg = AseAtomsAdaptor.get_structure(atoms)
        relaxed_structure = structure_to_fdf(relaxed_pmg)
        structure_io.write_fdf(relaxed_structure, out_name)
        rankings.append((out_name, energy))

    if not rankings:
        print(color_text("  Error: every candidate failed to relax -- no ranking to show.", 'red'))
        return

    rankings.sort(key=lambda r: r[1])
    e_min = rankings[0][1]
    print(f"\n{color_text('ML-ranked structures (MACE-MP-0, relaxed energy, most stable first):', 'bold')}")
    print(f"  {'Rank':<5}{'File':<24}{'Energy (eV)':<16}{'dE (eV)':<10}")
    for rank, (out_name, energy) in enumerate(rankings, start=1):
        print(f"  {rank:<5}{out_name:<24}{energy:<16.4f}{energy - e_min:<10.4f}")
    print(color_text(
        "\n  Note: a relative comparison from a fast ML potential, not an absolute DFT "
        "formation energy -- use it to prioritize which candidate(s) to relax with SIESTA, "
        "not as a final answer.", 'yellow'))


def main():
    require_pyxtal()
    from pyxtal import pyxtal
    from pyxtal.msg import Error as PyxtalError
    from pyxtal.database.collection import Collection

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Casts one or more random structures compatible with a given symmetry group.", 'bold')}
Give a symmetry group and a composition (species + how many of each) -- pyxtal
places the atoms on randomly chosen, symmetry-compatible Wyckoff positions
for you. This is the inverse of stb-crystalbuilder: use crystalbuilder when
you already know the exact Wyckoff sites you want, use crystalcast when you
want valid candidate structures generated for you (e.g. as starting guesses
for structure prediction). --molecular packs whole rigid molecules (instead
of bare atoms) into the symmetry group. --ml-rank pre-screens generated
candidates by MACE-MP-0 relaxed energy. --analyze/--substitute/--subgroup/
--supergroup all instead operate on an existing structure read via
-f/--file: --analyze prints its Wyckoff decomposition, --substitute swaps
elements while preserving the symmetry framework, --subgroup/--supergroup
search for related lower/higher-symmetry structures.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --group 225 --species Ni O --num-ions 4 8\n"
               "  %(prog)s --group Fd-3m --species Fe O --num-ions 8 16 \\\n"
               "      --count 5 --seed 42 -o spinel.fdf\n"
               "  %(prog)s --group 225 --species Ni O --num-ions 4 8 --count 5 --ml-rank\n"
               "  %(prog)s --group 225 --species Na Cl --num-ions 4 4 \\\n"
               "      --lattice 5.64 5.64 5.64 90 90 90\n"
               "  %(prog)s --group 225 --species Na Cl --num-ions 4 4 --sites 4a 4b\n"
               "  %(prog)s --dim 2 --group 65 --species C --num-ions 6 --thickness 3.4\n"
               "  %(prog)s --dim 0 --group D3d --species C --num-ions 6 --vacuum 12\n"
               "  %(prog)s --molecular --group 19 --species H2O --num-ions 4\n"
               "  %(prog)s --molecular --group 14 --species aspirin --num-ions 4 -o aspirin.fdf\n"
               "  %(prog)s --list-molecules\n"
               "  %(prog)s --analyze -f spinel_1.fdf\n"
               "  %(prog)s --substitute Cl:F -f rocksalt.fdf -o rocksalt_f.fdf\n"
               "  %(prog)s --subgroup -f spinel_1.fdf --count 5 -o distorted.fdf\n"
               "  %(prog)s --supergroup --target-group 225 -f distorted_1.fdf -o parent.fdf\n"
    )

    parser.add_argument("--analyze", action="store_true",
                        help="Analyze an existing structure instead of generating one: reads "
                             "-f/--file and prints its Wyckoff decomposition. Mutually exclusive "
                             "with --substitute/--subgroup/--supergroup and the generation "
                             "options below.")
    parser.add_argument("--substitute", nargs="+", default=None, metavar="OLD:NEW",
                        help="Substitute elements in an existing structure (read via -f/--file), "
                             "preserving its symmetry framework: one or more OLD:NEW pairs, e.g. "
                             "--substitute Cl:F Na:K. Mutually exclusive with --analyze/"
                             "--subgroup/--supergroup and the generation options below.")
    parser.add_argument("--subgroup", action="store_true",
                        help="Generate lower-symmetry subgroup variant(s) of an existing "
                             "structure (read via -f/--file) -- a symmetry-breaking distortion. "
                             "Up to --count candidates are written (numbered like --count > 1 "
                             "in generation mode). Optionally narrow the search with "
                             "--target-group/--group-type/--eps. Mutually exclusive with "
                             "--analyze/--substitute/--supergroup and the generation options.")
    parser.add_argument("--supergroup", action="store_true",
                        help="Search for higher-symmetry supergroup structure(s) of an existing "
                             "structure (read via -f/--file) towards --target-group (required -- "
                             "the installed pyxtal doesn't support auto-searching all possible "
                             "supergroups). Mutually exclusive with --analyze/--substitute/"
                             "--subgroup and the generation options below.")
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
                        help="Input .fdf structure file. Required with --analyze/--substitute/"
                             "--subgroup/--supergroup, unused otherwise.")

    parser.add_argument("--dim", type=int, choices=[3, 2, 1, 0], default=3,
                        help="Structure dimensionality: 3 = bulk (space group, default), "
                             "2 = layer (layer group, periodic in a/b with vacuum along c), "
                             "1 = rod/wire (rod group, periodic along c with vacuum in a/b), "
                             "0 = isolated cluster (point group, no periodicity). "
                             "--molecular does not support --dim 0 (upstream pyxtal limitation).")
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
                             "quickly relaxes each one's positions with the MACE-MP-0 potential "
                             "(needs the optional 'ml' extra: pip install stb_suite[ml]), prints "
                             "them ranked by relaxed energy, and overwrites each output file "
                             "with its relaxed geometry. A relative comparison from a fast ML "
                             "potential, not an absolute DFT formation energy.")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision: for generation, used in the --dim 3 post-build "
                             "verification step; for --analyze/--substitute/--subgroup/"
                             "--supergroup, the tolerance passed to pyxtal's own symmetry "
                             "detection (default: 1e-3, matches stb-symmetry/stb-unitcell).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Angle tolerance in degrees for reading an existing structure "
                             "(--analyze/--substitute/--subgroup/--supergroup only; default: "
                             "5.0, matches stb-symmetry/stb-unitcell).")
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
        names = list(Collection("molecules"))
        print(color_text(f"Molecules bundled with pyxtal ({len(names)}):", 'bold'))
        for name in names:
            print(f"  {name}")
        sys.exit(0)

    mode_flags = {"analyze": args.analyze, "substitute": bool(args.substitute),
                  "subgroup": args.subgroup, "supergroup": args.supergroup}
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

    if args.ml_rank:
        require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    for flag, value, valid_modes in (
        ("--eps", args.eps, ("subgroup",)),
        ("--group-type", args.group_type, ("subgroup",)),
        ("--d-tol", args.d_tol, ("supergroup",)),
        ("--target-group", args.target_group, ("subgroup", "supergroup")),
    ):
        if value is not None and mode not in valid_modes:
            print(color_text(f"Note: {flag} is ignored (only used with --{'/--'.join(valid_modes)}).", 'yellow'))

    if mode == "analyze":
        print("\n" + color_text("Analyze a structure's Wyckoff decomposition:", 'bold'))
        print("-" * 60)
        run_analyze(args)
        return
    if mode == "substitute":
        print("\n" + color_text("Substitute elements in an existing structure:", 'bold'))
        print("-" * 60)
        run_substitute(args)
        return
    if mode == "subgroup":
        print("\n" + color_text("Search for lower-symmetry subgroup structure(s):", 'bold'))
        print("-" * 60)
        run_subgroup(args)
        return
    if mode == "supergroup":
        print("\n" + color_text("Search for higher-symmetry supergroup structure(s):", 'bold'))
        print("-" * 60)
        run_supergroup(args)
        return

    print("\n" + color_text("Cast random structure(s) from a symmetry group:", 'bold'))
    print("-" * 60)

    if len(args.species) != len(args.num_ions):
        print(color_text(
            f"Error: --species has {len(args.species)} entries but --num-ions has "
            f"{len(args.num_ions)} -- they must match one-to-one.", 'red'))
        sys.exit(1)

    if args.molecular:
        try:
            validate_molecule_species(args.species)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)
    else:
        for symbol in args.species:
            try:
                Element(symbol)
            except ValueError as e:
                print(color_text(f"Error: {e}", 'red'))
                sys.exit(1)

    if args.count < 1:
        print(color_text("Error: --count must be at least 1.", 'red'))
        sys.exit(1)

    group = resolve_group(args.group)
    requested_number = None
    if args.dim == 3:
        try:
            requested_number = SpaceGroup(str(group)).int_number if isinstance(group, str) \
                else SpaceGroup.from_int_number(group).int_number
        except ValueError:
            requested_number = None

    species_label = "Molecule" if args.molecular else "Species"
    unit_label = "molecules" if args.molecular else "atoms"
    print(f"  {color_text('Dimension:', 'cyan')} {args.dim}")
    print(f"  {color_text(f'Requested {DIM_LABELS[args.dim]}:', 'cyan')} {args.group}")
    for symbol, n in zip(args.species, args.num_ions):
        print(f"  {color_text(f'{species_label}:', 'cyan')} {symbol} x {n}")
    print(f"  {color_text('Structures requested:', 'cyan')} {args.count}")
    if args.molecular and args.seed is not None:
        print(color_text(
            "  Note: --seed only reproduces the lattice for --molecular structures -- "
            "molecule orientations still vary between runs (an upstream pyxtal limitation).",
            'yellow'))
    for flag, value, required_dim in (
        ("--thickness", args.thickness, 2), ("--area", args.area, 1), ("--vacuum", args.vacuum, 0)
    ):
        if value is not None and args.dim != required_dim:
            print(color_text(f"  Note: {flag} is ignored (only used with --dim {required_dim}).", 'yellow'))
    vacuum = args.vacuum if args.vacuum is not None else 10.0

    fixed_lattice = None
    if args.lattice is not None:
        from pyxtal.symmetry import Group
        from pyxtal.lattice import Lattice as PyxtalLattice
        ltype = Group(group, dim=args.dim).lattice_type
        mismatch = check_lattice_type(*args.lattice, ltype=ltype)
        if mismatch is not None:
            print(color_text(f"Error: {mismatch}", 'red'))
            sys.exit(1)
        try:
            fixed_lattice = PyxtalLattice.from_para(*args.lattice, ltype=ltype)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)
        print(f"  {color_text('Fixed lattice:', 'cyan')} a={args.lattice[0]} b={args.lattice[1]} "
              f"c={args.lattice[2]} alpha={args.lattice[3]} beta={args.lattice[4]} "
              f"gamma={args.lattice[5]} ({ltype})")

    sites = None
    if args.sites is not None:
        if len(args.sites) != len(args.species):
            print(color_text(
                f"Error: --sites has {len(args.sites)} entries but --species has "
                f"{len(args.species)} -- one --sites entry is needed per --species.", 'red'))
            sys.exit(1)
        sites = [[wp.strip() for wp in entry.split(",")] for entry in args.sites]

        from pyxtal.symmetry import Group as SitesGroup
        try:
            site_group = SitesGroup(group, dim=args.dim)
            all_zero_dof = True
            for symbol, wp_letters in zip(args.species, sites):
                print(f"  {color_text(f'{species_label} sites:', 'cyan')} {symbol} -> {', '.join(wp_letters)}")
                for letter in wp_letters:
                    if site_group.get_wyckoff_position(letter).get_dof() > 0:
                        all_zero_dof = False
        except (IndexError, ValueError) as e:
            print(color_text(
                f"Error: --sites contains a Wyckoff label not valid for this group -- {e}", 'red'))
            sys.exit(1)

        if args.count > 1 and all_zero_dof:
            print(color_text(
                f"  Note: every assigned Wyckoff site has zero free parameters -- all "
                f"{args.count} requested structures will be identical (--sites leaves "
                "nothing left to randomize; --seed won't change that either).", 'yellow'))

    rng = np.random.default_rng(args.seed)

    try:
        names = resolve_output_paths(args.output, args.count)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    written = []
    results = []
    print()
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
            # this attempt and keep whatever was already written, rather than discarding
            # earlier successes in the same --count batch. A --group/composition that's
            # truly invalid will simply fail again on every remaining attempt and the
            # "none of the requested structures" error below will catch it.
            print(color_text(f"  Attempt #{i}: {e} -- skipped.", 'yellow'))
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

        species_order = None if args.molecular else args.species
        new_structure = structure_to_fdf(pmg_structure, species_order=species_order, coord_format=coord_format)

        # Checked on overlap_check_target (the raw dim=0 Molecule, not the boxed
        # Structure) because the box's own periodicity would otherwise make a real
        # atom look "close" to its own periodic image across the artificial vacuum cell.
        min_dist = structure_io.min_pairwise_distance(overlap_check_target)
        if min_dist is not None and min_dist < 0.5:
            print(color_text(
                f"  Warning: structure #{i} has atoms unusually close together "
                f"({min_dist:.3f} Ang) -- check --volume-factor.", 'yellow'))

        out_name = names[i - 1]
        structure_io.write_fdf(new_structure, out_name)
        written.append(out_name)
        results.append((out_name, pmg_structure, args.dim == 0))

        if args.dim == 3:
            sga = SpacegroupAnalyzer(pmg_structure, symprec=args.symprec)
            detected_number = sga.get_space_group_number()
            detected_symbol = sga.get_space_group_symbol()
            print(f"  {color_text('->', 'green')} #{i}: {pmg_structure.composition.reduced_formula} "
                  f"({len(pmg_structure)} atoms), space group {detected_symbol} "
                  f"(No. {detected_number}): {out_name}")
            if requested_number is not None and requested_number != detected_number:
                if detected_number > requested_number:
                    print(color_text(
                        f"     Note: ended up at a higher-symmetry space group than requested "
                        f"(No. {requested_number}) -- the random placement landed on special "
                        "positions. Still valid, just be aware.", 'yellow'))
                else:
                    print(color_text(
                        f"     Warning: detected space group (No. {detected_number}) differs "
                        f"from the requested one (No. {requested_number}) and is not "
                        "higher-symmetry -- the generated structure may not actually have the "
                        "intended symmetry.", 'yellow'))
        else:
            print(f"  {color_text('->', 'green')} #{i}: {pmg_structure.composition.reduced_formula} "
                  f"({len(pmg_structure)} atoms), {DIM_LABELS[args.dim]} "
                  f"{crystal.group.symbol} (No. {crystal.group.number}): {out_name}")

    if not written:
        print(color_text(
            "\nError: none of the requested structures could be generated -- "
            "see the per-attempt messages above.", 'red'))
        sys.exit(1)

    if args.ml_rank:
        run_ml_rank(results)

    if len(written) < args.count:
        print(color_text(
            f"\nPartial success: {len(written)} of {args.count} structure(s) written "
            f"({args.count - len(written)} skipped -- see warnings above).", 'yellow'))
        sys.exit(1)

    print(f"\n{color_text('Success:', 'green')} {len(written)} of {args.count} structure(s) written.")


if __name__ == "__main__":
    main()
