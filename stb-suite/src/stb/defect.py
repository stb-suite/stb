#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.10.0"

import sys
import os
import argparse
import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro


def parse_index_list(spec, n_atoms):
    """Parses '3,7' into 0-indexed [2, 6], validating 1-based range and no duplicates."""
    try:
        indices_1based = [int(x) for x in spec.split(",")]
    except ValueError:
        raise ValueError(f"--index must be a comma-separated list of integers, got '{spec}'.")

    if len(set(indices_1based)) != len(indices_1based):
        raise ValueError(f"--index contains duplicate entries: {spec}")

    for idx in indices_1based:
        if not (1 <= idx <= n_atoms):
            raise ValueError(f"--index {idx} is out of range (1 to {n_atoms}).")

    return [idx - 1 for idx in indices_1based]


def convert_position(position, given_format, structure):
    """Converts a user-given position into whatever coord_format `structure` uses."""
    position = np.array(position, dtype=float)
    if given_format == structure.coord_format:
        return position
    if given_format == "cartesian":
        return position @ np.linalg.inv(structure.lattice)
    return position @ structure.lattice


def resolve_nearest(structure, position, given_format, filter_species=None):
    """Index (0-based) of the atom closest to `position`, using pymatgen's
    periodic minimum-image distance (PeriodicSite.distance) so a target near
    a cell edge still resolves correctly through the periodic image.

    PeriodicSite always wants fractional coordinates relative to the given
    lattice, regardless of the source file's own coord_format -- this is
    independent of convert_position(), which instead targets whatever
    coord_format the .fdf itself uses.
    """
    pmg_structure = structure_io.to_pymatgen(structure)
    position = np.array(position, dtype=float)
    if given_format == "cartesian":
        frac_position = pmg_structure.lattice.get_fractional_coords(position)
    else:
        frac_position = position
    dummy = PeriodicSite("X", frac_position, pmg_structure.lattice)

    candidate_indices = list(range(len(pmg_structure)))
    if filter_species is not None:
        candidate_indices = [
            i for i in candidate_indices if pmg_structure[i].specie.symbol == filter_species
        ]
        if not candidate_indices:
            raise ValueError(f"No atoms of species '{filter_species}' found in the structure.")

    distances = [(dummy.distance(pmg_structure[i]), i) for i in candidate_indices]
    distances.sort(key=lambda d: d[0])
    return distances[0][1]


def find_inequivalent_sites(pmg_structure, symprec, filter_species=None):
    """Returns (sites, space_group_label) where `sites` is a list of
    (index, wyckoff_letter, multiplicity) -- one representative atom (0-based
    index) per symmetrically distinct site, via spglib's equivalent-atoms
    mapping. If filter_species is given, only representatives of that
    species are returned (their multiplicity still counts all symmetry
    -equivalent atoms of that site, filtered species or not, since
    spglib never groups atoms of different species together).
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
    dataset = sga.get_symmetry_dataset()

    groups = {}
    for i, rep in enumerate(dataset.equivalent_atoms):
        groups.setdefault(rep, []).append(i)

    sites = []
    for rep in sorted(groups):
        symbol = pmg_structure[rep].specie.symbol
        if filter_species is not None and symbol != filter_species:
            continue
        sites.append((rep, dataset.wyckoffs[rep], len(groups[rep])))

    space_group = f"{dataset.international} (No. {dataset.number})"
    return sites, space_group


def apply_vacancy_or_substitution(atoms, species_list, species_meta, indices_set, defect_type, new_species=None):
    """Returns (new_atoms, new_species_list, new_species_meta) for removing (vacancy)
    or replacing (substitution) the atoms at `indices_set`, without mutating the inputs.
    """
    if defect_type == "vacancy":
        new_atoms = [pair for i, pair in enumerate(atoms) if i not in indices_set]
        return new_atoms, species_list, species_meta

    new_atoms = [(new_species if i in indices_set else sym, pos) for i, (sym, pos) in enumerate(atoms)]
    new_species_list = list(species_list)
    if new_species not in new_species_list:
        new_species_list.append(new_species)
    new_species_meta = structure_io.ensure_species_id(species_meta, new_species)
    return new_atoms, new_species_list, new_species_meta


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Introduces a point defect (vacancy, substitution, interstitial) into a SIESTA FDF structure.", 'bold')}
Select the site by raw atom index, by the position closest to a target, or
automatically at every symmetrically distinct site.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf --type vacancy --index 5\n"
               "  %(prog)s -f structure.fdf --type substitution --index 5 --new-species Ge\n"
               "  %(prog)s -f structure.fdf --type vacancy --nearest 0.5 0.5 0.5 --filter-species O\n"
               "  %(prog)s -f structure.fdf --type interstitial --position 0.5 0.5 0.5 --species N\n"
               "  %(prog)s -f structure.fdf --type vacancy --all-inequivalent-sites --filter-species O\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf).")
    parser.add_argument("--type", choices=["vacancy", "substitution", "interstitial"], required=True,
                        help="Defect type to introduce.")

    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--index", type=str, default=None,
                           help="Comma-separated 1-indexed atom(s) to remove/substitute, e.g. '3,7'.")
    selection.add_argument("--nearest", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
                           help="Selects the atom closest to this position (vacancy/substitution).")
    selection.add_argument("--all-inequivalent-sites", action="store_true",
                           help="Generate one output structure per symmetrically distinct site "
                                "(via spglib), instead of a single combined structure -- e.g. to "
                                "compare defect formation energies across all inequivalent sites "
                                "rather than picking one by hand. Writes '<output>_site<N>.fdf' "
                                "per site instead of a single --output file. Narrow to one "
                                "species with --filter-species.")
    parser.add_argument("--nearest-format", choices=["fractional", "cartesian"], default="fractional",
                        help="How to interpret --nearest. Default: fractional.")
    parser.add_argument("--filter-species", type=str, default=None,
                        help="Restrict --nearest's or --all-inequivalent-sites' search to atoms "
                             "of this element.")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision for --all-inequivalent-sites "
                             "(default: 1e-3, matches stb-symmetry/stb-unitcell).")

    parser.add_argument("--new-species", type=str, default=None,
                        help="Replacement element symbol (required for --type substitution).")

    parser.add_argument("--position", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
                        help="Position of the new atom (required for --type interstitial).")
    parser.add_argument("--position-format", choices=["fractional", "cartesian"], default="fractional",
                        help="How to interpret --position. Default: fractional.")
    parser.add_argument("--species", type=str, default=None,
                        help="Element symbol of the new atom (required for --type interstitial).")

    parser.add_argument("-o", "--output", type=str, default="defect.fdf",
                        help="Output .fdf file name (default: defect.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-defect {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.type in ("vacancy", "substitution"):
        if args.index is None and args.nearest is None and not args.all_inequivalent_sites:
            parser.error("--type vacancy/substitution requires --index, --nearest, or --all-inequivalent-sites.")
        if args.position is not None or args.species is not None:
            parser.error("--position/--species are only valid with --type interstitial.")
        if args.type == "substitution" and args.new_species is None:
            parser.error("--type substitution requires --new-species.")
        if args.type == "vacancy" and args.new_species is not None:
            parser.error("--new-species is only valid with --type substitution.")
    else:
        if args.position is None or args.species is None:
            parser.error("--type interstitial requires --position and --species.")
        if args.index is not None or args.nearest is not None or args.new_species is not None:
            parser.error("--index/--nearest/--new-species are only valid with --type vacancy/substitution.")
        if args.all_inequivalent_sites:
            parser.error("--all-inequivalent-sites is only valid with --type vacancy/substitution.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Introduce a point defect into a structure:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    n_atoms = len(structure.atoms)
    pmg_before = structure_io.to_pymatgen(structure)
    print(f"  {color_text('Input formula:', 'cyan')} {pmg_before.composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {n_atoms}")
    print(f"  {color_text('Defect type:', 'cyan')} {args.type}")

    for symbol in filter(None, [args.new_species, args.species]):
        try:
            Element(symbol)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

    atoms = list(structure.atoms)
    species_list = list(structure.species)
    species_meta = dict(structure.species_meta)

    if args.all_inequivalent_sites:
        if args.filter_species is not None:
            try:
                Element(args.filter_species)
            except ValueError as e:
                print(color_text(f"Error: {e}", 'red'))
                sys.exit(1)

        sites, space_group = find_inequivalent_sites(pmg_before, args.symprec, args.filter_species)
        if not sites:
            print(color_text(
                f"Error: no atoms of species '{args.filter_species}' found in the structure.", 'red'))
            sys.exit(1)

        print(f"  {color_text('Space group:', 'cyan')} {space_group}")
        print(f"\n  {color_text('Symmetrically distinct sites found:', 'cyan')} {len(sites)}")
        for idx, wyckoff, multiplicity in sites:
            print(f"    #{idx + 1:<4} {atoms[idx][0]:<3} Wyckoff {wyckoff:<3} (multiplicity {multiplicity})")

        if args.type == "substitution":
            print(f"\n  {color_text('Substitution:', 'cyan')} -> {args.new_species}")

        stem, ext = os.path.splitext(args.output)
        ext = ext or ".fdf"
        written = []
        print()
        for idx, wyckoff, multiplicity in sites:
            one_indexed = idx + 1
            new_atoms, new_species_list, new_species_meta = apply_vacancy_or_substitution(
                atoms, species_list, species_meta, {idx}, args.type, args.new_species)

            new_structure = structure_io.FdfStructure(
                lattice=structure.lattice,
                lattice_constant=structure.lattice_constant,
                species=new_species_list,
                species_meta=new_species_meta,
                atoms=new_atoms,
                coord_format=structure.coord_format,
            )
            out_name = f"{stem}_site{one_indexed}{ext}"
            structure_io.write_fdf(new_structure, out_name)
            written.append(out_name)
            print(f"  {color_text('->', 'green')} site #{one_indexed} "
                  f"({atoms[idx][0]}, Wyckoff {wyckoff}): {out_name}")

        print(f"\n{color_text('Success:', 'green')} {len(written)} structure(s) written.")
        return

    if args.type in ("vacancy", "substitution"):
        try:
            if args.index is not None:
                indices = parse_index_list(args.index, n_atoms)
            else:
                if args.filter_species is not None:
                    try:
                        Element(args.filter_species)
                    except ValueError as e:
                        print(color_text(f"Error: {e}", 'red'))
                        sys.exit(1)
                indices = [resolve_nearest(structure, args.nearest, args.nearest_format, args.filter_species)]
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

        indices_set = set(indices)
        selected = [(i + 1, atoms[i][0], atoms[i][1]) for i in sorted(indices_set)]
        for one_indexed, symbol, position in selected:
            print(f"  {color_text('Selected site:', 'cyan')} #{one_indexed} ({symbol}) at {position}")

        if args.type == "vacancy":
            atoms = [pair for i, pair in enumerate(atoms) if i not in indices_set]
        else:
            atoms = [(args.new_species if i in indices_set else sym, pos) for i, (sym, pos) in enumerate(atoms)]
            print(f"  {color_text('Substitution:', 'cyan')} -> {args.new_species}")
            if args.new_species not in species_list:
                species_list.append(args.new_species)
            species_meta = structure_io.ensure_species_id(species_meta, args.new_species)

    else:
        position = convert_position(args.position, args.position_format, structure)
        atoms = atoms + [(args.species, position)]
        print(f"  {color_text('Interstitial:', 'cyan')} {args.species} at {args.position} ({args.position_format})")
        if args.species not in species_list:
            species_list.append(args.species)
        species_meta = structure_io.ensure_species_id(species_meta, args.species)

    new_structure = structure_io.FdfStructure(
        lattice=structure.lattice,
        lattice_constant=structure.lattice_constant,
        species=species_list,
        species_meta=species_meta,
        atoms=atoms,
        coord_format=structure.coord_format,
    )

    pmg_after = structure_io.to_pymatgen(new_structure)
    print(f"\n  {color_text('Output formula:', 'cyan')} {pmg_after.composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {len(atoms)}")

    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
