#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import os
import sys
import difflib
import argparse
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.core.periodic_table import Element
from stb.core import structure_io
from stb.core.cli import color_text, show_intro
from stb.core.deps import require_pyxtal

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


def run_analyze(args):
    from pyxtal import pyxtal

    if not args.file:
        print(color_text("Error: --analyze requires -f/--file.", 'red'))
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

    print(f"\n  {color_text('Space group:', 'cyan')} {crystal.group.symbol} (No. {crystal.group.number})")
    print(f"  {color_text('Formula:', 'cyan')} {pmg_structure.composition.reduced_formula}")
    print(f"  {color_text('Symmetrically distinct sites:', 'cyan')} {len(crystal.atom_sites)}")

    print(f"\n{color_text('Wyckoff sites (paste each --site straight into stb-crystalbuilder):', 'bold')}")
    for site in crystal.atom_sites:
        x, y, z = site.position
        print(f"  --site {site.specie:<3} {x:.6f} {y:.6f} {z:.6f}"
              f"   {color_text(f'(Wyckoff {site.wp.multiplicity}{site.wp.letter})', 'yellow')}")


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
for structure prediction). --analyze runs the reverse direction: given an
existing structure, print its Wyckoff decomposition. --molecular packs whole
rigid molecules (instead of bare atoms) into the symmetry group.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --group 225 --species Ni O --num-ions 4 8\n"
               "  %(prog)s --group Fd-3m --species Fe O --num-ions 8 16 \\\n"
               "      --count 5 --seed 42 -o spinel.fdf\n"
               "  %(prog)s --dim 2 --group 65 --species C --num-ions 6 --thickness 3.4\n"
               "  %(prog)s --dim 0 --group D3d --species C --num-ions 6 --vacuum 12\n"
               "  %(prog)s --molecular --group 19 --species H2O --num-ions 4\n"
               "  %(prog)s --molecular --group 14 --species aspirin --num-ions 4 -o aspirin.fdf\n"
               "  %(prog)s --list-molecules\n"
               "  %(prog)s --analyze -f spinel_1.fdf\n"
    )

    parser.add_argument("--analyze", action="store_true",
                        help="Analyze an existing structure instead of generating one: reads "
                             "-f/--file and prints its Wyckoff decomposition. All generation "
                             "options below are ignored in this mode.")
    parser.add_argument("-f", "--file", type=str, default=None,
                        help="Input .fdf structure file. Required with --analyze, unused otherwise.")

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
                             "Required unless --analyze is given.")
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
                             "the structure, e.g. --species Ni O. Required unless --analyze or "
                             "--list-molecules is given.")
    parser.add_argument("--num-ions", nargs="+", type=int, default=None,
                        help="Number of atoms (or, with --molecular, molecules) of each "
                             "--species, same order and count, e.g. --num-ions 4 8. Required "
                             "unless --analyze or --list-molecules is given.")
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
    parser.add_argument("--max-attempts", type=int, default=10,
                        help="Internal retries pyxtal allows itself per structure to find a "
                             "non-overlapping placement before giving up on it (default: 10).")
    parser.add_argument("--count", type=int, default=1,
                        help="Number of independent random structures to generate (default: 1).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed. Same seed + same inputs reproduces the same batch "
                             "of structures (default: not seeded, different every run). With "
                             "--molecular, only the lattice is reproduced this way -- molecule "
                             "orientations are not (an upstream pyxtal limitation).")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision: for generation, used in the --dim 3 post-build "
                             "verification step; for --analyze, the tolerance passed to pyxtal's "
                             "own symmetry detection (default: 1e-3, matches stb-symmetry/"
                             "stb-unitcell).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Angle tolerance in degrees, --analyze only (default: 5.0, matches "
                             "stb-symmetry/stb-unitcell).")
    parser.add_argument("-o", "--output", type=str, default="crystalcast.fdf",
                        help="Output .fdf file name (default: crystalcast.fdf). With "
                             "--count > 1, each structure is written as '<output>_<N>.fdf'; "
                             "the tool refuses to run if any of those numbered names already "
                             "exist, to avoid silently overwriting unrelated files.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-crystalcast {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.list_molecules:
        names = list(Collection("molecules"))
        print(color_text(f"Molecules bundled with pyxtal ({len(names)}):", 'bold'))
        for name in names:
            print(f"  {name}")
        sys.exit(0)

    if not args.analyze:
        if args.molecular and args.dim == 0:
            parser.error("--molecular does not support --dim 0 (upstream pyxtal limitation). "
                         "Use --dim 3, 2, or 1.")
        if not args.group:
            parser.error("--group is required unless --analyze or --list-molecules is given.")
        if not args.species:
            parser.error("--species is required unless --analyze or --list-molecules is given.")
        if not args.num_ions:
            parser.error("--num-ions is required unless --analyze or --list-molecules is given.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.analyze:
        print("\n" + color_text("Analyze a structure's Wyckoff decomposition:", 'bold'))
        print("-" * 60)
        run_analyze(args)
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

    rng = np.random.default_rng(args.seed)
    stem, ext = os.path.splitext(args.output)
    ext = ext or ".fdf"

    if args.count > 1:
        candidate_names = [f"{stem}_{i}{ext}" for i in range(1, args.count + 1)]
        preexisting = [name for name in candidate_names if os.path.exists(name)]
        if preexisting:
            print(color_text(
                f"Error: refusing to overwrite existing file(s): {', '.join(preexisting)} -- "
                "move them aside or choose a different -o.", 'red'))
            sys.exit(1)

    written = []
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

        species_meta = {}
        element_symbols = dict.fromkeys(site.specie.symbol for site in pmg_structure) \
            if args.molecular else args.species
        for symbol in element_symbols:
            species_meta = structure_io.ensure_species_id(species_meta, symbol)
        new_structure = structure_io.from_pymatgen(
            pmg_structure, species_meta=species_meta, coord_format=coord_format)

        # Checked on overlap_check_target (the raw dim=0 Molecule, not the boxed
        # Structure) because the box's own periodicity would otherwise make a real
        # atom look "close" to its own periodic image across the artificial vacuum cell.
        min_dist = structure_io.min_pairwise_distance(overlap_check_target)
        if min_dist is not None and min_dist < 0.5:
            print(color_text(
                f"  Warning: structure #{i} has atoms unusually close together "
                f"({min_dist:.3f} Ang) -- check --volume-factor.", 'yellow'))

        out_name = f"{stem}_{i}{ext}" if args.count > 1 else args.output
        structure_io.write_fdf(new_structure, out_name)
        written.append(out_name)

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

    if len(written) < args.count:
        print(color_text(
            f"\nPartial success: {len(written)} of {args.count} structure(s) written "
            f"({args.count - len(written)} skipped -- see warnings above).", 'yellow'))
        sys.exit(1)

    print(f"\n{color_text('Success:', 'green')} {len(written)} of {args.count} structure(s) written.")


if __name__ == "__main__":
    main()
