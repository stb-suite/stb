#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import sys
import os
import time
import argparse
from datetime import datetime

import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_mace
from stb.core.symmetry import find_inequivalent_sites

REPORT_FILE = "stb_defect_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-fetch/stb-kgrid/stb-slab
# (core/kspace.py's other callers), used to detect vacuum-padded axes on
# both the input structure and every output defect structure.
VACUUM_GAP_ANG = 10.0


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


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a
    space-group label -- same shape as the rest of the suite's own
    _validate_structure(), wrapped in try/except by the caller (a
    validation failure is reported, never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group      : {sg_label}", f_out)
    return sg_label


def _relax_structure(pmg_structure, calc, vacuum_axes, relax_cell, f_out):
    """Relaxes one structure with an already-loaded MACE calculator, returns
    (relaxed_pmg_structure, e0, e1, converged, steps_used). Shared by both
    the single-defect and --all-inequivalent-sites paths so a --ml-relax run
    behaves identically either way."""
    atoms_ase = AseAtomsAdaptor.get_atoms(pmg_structure)
    atoms_ase.calc = calc
    e0 = atoms_ase.get_potential_energy()
    cell_mask = mace_relax.build_cell_mask(vacuum_axes) if relax_cell else None
    converged, steps_used = mace_relax.relax(atoms_ase, calc, cell_mask=cell_mask, fmax=0.05, max_steps=200)
    e1 = atoms_ase.get_potential_energy()
    return AseAtomsAdaptor.get_structure(atoms_ase), e0, e1, converged, steps_used


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
               "  %(prog)s -f structure.fdf --type vacancy --index 5 --ml-relax --save-report --view\n"
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
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry precision for --all-inequivalent-sites and the before/after "
                             "symmetry table (default: 0.01, pymatgen's own default, matches "
                             "stb-symmetry/stb-unitcell).")
    parser.add_argument("--ml-rank", action="store_true",
                        help="Only valid with --all-inequivalent-sites. Quickly relaxes each "
                             "candidate site's local geometry (positions only) with a MACE "
                             "potential (needs the optional 'ml' extra: pip install "
                             "stb_suite[ml]), then prints them ranked by relaxed energy -- a "
                             "fast pre-screen for which site is most likely favorable, before "
                             "spending DFT time on all of them. Overwrites each site's output "
                             "file with its relaxed geometry. This is a relative comparison "
                             "from a fast ML potential, not an absolute DFT formation energy. "
                             "Mutually exclusive with --ml-relax (which already relaxes every "
                             "output structure, just without the ranking table).")

    parser.add_argument("--new-species", type=str, default=None,
                        help="Replacement element symbol (required for --type substitution).")

    parser.add_argument("--position", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
                        help="Position of the new atom (required for --type interstitial).")
    parser.add_argument("--position-format", choices=["fractional", "cartesian"], default="fractional",
                        help="How to interpret --position. Default: fractional.")
    parser.add_argument("--species", type=str, default=None,
                        help="Element symbol of the new atom (required for --type interstitial).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax every output structure with a MACE potential before "
                             "writing it out (needs the optional 'ml' extra: pip install "
                             "stb_suite[ml]) -- positions only by default. Off by default. "
                             "Mutually exclusive with --ml-rank.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the cell -- any vacuum-padded axis "
                             "always stays exactly fixed. Only valid together with --ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax/--ml-rank (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax/--ml-rank, "
                             "instead of a MACE-MP-0 foundation size.")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the MACE model on (default: cpu).")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the input "
                             "structure and every output structure (page through frames in "
                             "ase-gui) after writing the output file(s). Needs a display. "
                             "Off by default.")

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

    if args.ml_rank and not args.all_inequivalent_sites:
        parser.error("--ml-rank is only valid with --all-inequivalent-sites.")
    if args.ml_rank and args.ml_relax:
        parser.error("--ml-rank and --ml-relax are mutually exclusive (--ml-rank already relaxes "
                      "every candidate, just without --ml-relax-cell's cell relaxation).")
    if args.ml_relax_cell and not args.ml_relax:
        parser.error("--ml-relax-cell is only valid together with --ml-relax.")
    if (args.custom_model or args.model != "small") and not (args.ml_relax or args.ml_rank):
        parser.error("--model/--custom-model are only valid together with --ml-relax/--ml-rank.")

    if args.ml_rank or args.ml_relax:
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

    def fail(message):
        print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-DEFECT REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time        : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file       : {args.filename}", f_out)
    print_dual(f"Defect type      : {args.type}", f_out)
    print_dual(f"Output file      : {args.output}", f_out)
    print_dual(f"ML pre-relax     : {model_desc if (args.ml_relax or args.ml_rank) else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell       : {'yes (vacuum axes fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    if not os.path.exists(args.filename):
        fail(f"File '{args.filename}' not found.")

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        fail(str(e))

    n_atoms = len(structure.atoms)
    pmg_before = structure_io.to_pymatgen(structure)

    print_section("[1] INPUT STRUCTURE", f_out)
    frac_before = [site.frac_coords for site in pmg_before]
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_before, pmg_before.lattice.matrix, VACUUM_GAP_ANG)
    print_dual(f"Formula          : {pmg_before.composition.reduced_formula}", f_out)
    print_dual(f"Atoms            : {n_atoms}", f_out)
    print_dual(f"Dimensionality   : {kspace.dimensionality_label(vacuum_axes_before)}", f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    try:
        _validate_structure(pmg_before, vacuum_axes_before, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    for symbol in filter(None, [args.new_species, args.species]):
        try:
            Element(symbol)
        except ValueError as e:
            fail(str(e))

    atoms = list(structure.atoms)
    species_list = list(structure.species)
    species_meta = dict(structure.species_meta)

    print_section("[3] SITE SELECTION", f_out)

    if args.all_inequivalent_sites:
        if args.filter_species is not None:
            try:
                Element(args.filter_species)
            except ValueError as e:
                fail(str(e))

        sites, space_group = find_inequivalent_sites(pmg_before, args.symprec, args.filter_species)
        if not sites:
            fail(f"no atoms of species '{args.filter_species}' found in the structure.")

        print_dual(f"Space group      : {space_group}", f_out)
        print_dual(f"Symmetrically distinct sites found: {len(sites)}", f_out)
        rows = [([f"#{idx + 1}", atoms[idx][0], wyckoff, str(multiplicity)], None)
                for idx, wyckoff, multiplicity in sites]
        print_table(["Site", "Species", "Wyckoff", "Multiplicity"], rows, f_out)

        if args.type == "substitution":
            print_dual(f"Substitution     : -> {args.new_species}", f_out)

        model_arg = args.custom_model if args.custom_model else args.model
        calc = None
        if args.ml_rank or args.ml_relax:
            try:
                calc = mace_relax.get_calculator(model_arg, device=args.device)
            except ValueError as e:
                fail(str(e))

        stem, ext = os.path.splitext(args.output)
        ext = ext or ".fdf"
        results = []  # (one_indexed, symbol, wyckoff, pmg_site_structure, out_name)

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
            site_pmg = structure_io.to_pymatgen(new_structure)
            out_name = f"{stem}_site{one_indexed}{ext}"
            results.append([one_indexed, atoms[idx][0], wyckoff, site_pmg, new_species_meta, out_name, None])

        if args.ml_rank:
            print_section("[4] ML RANKING (MACE)", f_out)
            print_dual(f"Model            : {model_desc}", f_out)
            print_dual(f"Device           : {args.device}", f_out)
            print_dual("Relaxing each candidate site (positions only)...", f_out)
            for row in results:
                site_pmg = row[3]
                frac_site = [site.frac_coords for site in site_pmg]
                vacuum_axes_site = kspace.detect_vacuum_axes(frac_site, site_pmg.lattice.matrix, VACUUM_GAP_ANG)
                relaxed_pmg, e0, e1, converged, steps_used = _relax_structure(
                    site_pmg, calc, vacuum_axes_site, relax_cell=False, f_out=f_out)
                row[3] = relaxed_pmg
                row[6] = e1

            rankings = sorted(results, key=lambda r: r[6])
            e_min = rankings[0][6]
            print_dual(color_text("ML-ranked sites (relaxed energy, most stable first):", 'bold'), f_out)
            print_dual(f"  {'Rank':<5}{'Site':<7}{'Species':<9}{'Wyckoff':<9}{'Energy (eV)':<14}{'dE (eV)':<10}", f_out)
            for rank, (one_indexed, symbol, wyckoff, _, _, _, energy) in enumerate(rankings, start=1):
                print_dual(f"  {rank:<5}#{one_indexed:<6}{symbol:<9}{wyckoff:<9}{energy:<14.4f}{energy - e_min:<10.4f}", f_out)
            print_dual(color_text(
                "Note: a relative comparison from a fast ML potential, not an absolute "
                "DFT formation energy -- use it to prioritize which site(s) to relax with "
                "SIESTA, not as a final answer.", 'yellow'), f_out)
        elif args.ml_relax:
            print_section("[4] ML PRE-RELAXATION (MACE)", f_out)
            print_dual(f"Model            : {model_desc}", f_out)
            print_dual(f"Device           : {args.device}", f_out)
            print_dual(f"Cell relaxation  : {'yes (vacuum axes fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)
            for row in results:
                one_indexed, symbol, wyckoff, site_pmg = row[0], row[1], row[2], row[3]
                frac_site = [site.frac_coords for site in site_pmg]
                vacuum_axes_site = kspace.detect_vacuum_axes(frac_site, site_pmg.lattice.matrix, VACUUM_GAP_ANG)
                relaxed_pmg, e0, e1, converged, steps_used = _relax_structure(
                    site_pmg, calc, vacuum_axes_site, relax_cell=args.ml_relax_cell, f_out=f_out)
                row[3] = relaxed_pmg
                row[6] = (converged, steps_used, e0, e1)
                print_dual(f"  Site #{one_indexed} ({symbol}, Wyckoff {wyckoff}): "
                           f"{steps_used} step(s), {'converged' if converged else 'NOT converged'}, "
                           f"dE = {e1 - e0:+.6f} eV", f_out)

        print_section("[5] STRUCTURE VALIDATION, SYMMETRY & WRITING OUTPUT FILE(S)", f_out)
        output_files = []
        view_atoms = [AseAtomsAdaptor.get_atoms(pmg_before)]
        for one_indexed, symbol, wyckoff, site_pmg, site_species_meta, out_name, ml_info in results:
            print_dual(f"\n{color_text(f'--- [Site {one_indexed} ({symbol}, Wyckoff {wyckoff})] ---', 'bold')}", f_out)
            print_dual(f"Formula          : {site_pmg.composition.reduced_formula}", f_out)
            print_dual(f"Atoms            : {len(site_pmg)}", f_out)

            frac_site = [site.frac_coords for site in site_pmg]
            vacuum_axes_site = kspace.detect_vacuum_axes(frac_site, site_pmg.lattice.matrix, VACUUM_GAP_ANG)
            try:
                _validate_structure(site_pmg, vacuum_axes_site, f_out)
            except Exception as e:
                print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

            before_info = core_symmetry.symmetry_summary(pmg_before, args.symprec, VACUUM_GAP_ANG)
            after_info = core_symmetry.symmetry_summary(site_pmg, args.symprec, VACUUM_GAP_ANG)
            if "Error" in before_info or "Error" in after_info:
                print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
            else:
                properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
                rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                        for prop in properties]
                print_table(["Property", "Before", "After"], rows, f_out)

            header_comment = [
                f"{args.type.capitalize()} defect at site #{one_indexed} ({symbol}, Wyckoff {wyckoff}) "
                f"introduced by stb-defect from {args.filename}.",
            ]
            if args.type == "substitution":
                header_comment.append(f"Substituted with {args.new_species}.")
            if ml_info is not None and args.ml_relax:
                converged, steps_used, e0, e1 = ml_info
                header_comment.append(
                    f"ML pre-relaxed with {model_desc} "
                    f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
                    f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV).")
            elif args.ml_rank:
                header_comment.append(f"ML-relaxed with {model_desc} for ranking (E = {ml_info:.6f} eV).")

            final_structure = structure_io.from_pymatgen(site_pmg, species_meta=site_species_meta)
            structure_io.write_fdf(final_structure, out_name, header_comment=header_comment)
            print_dual(color_text(f"[OK] Structure written to '{out_name}'.", 'green'), f_out)
            output_files.append(out_name)
            view_atoms.append(AseAtomsAdaptor.get_atoms(site_pmg))

        print_section("[6] REFERENCES", f_out)
        bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
        if args.ml_rank or args.ml_relax:
            bib_entries.append(citations.MACE)
            if not args.custom_model:
                bib_entries.append(citations.MACE_MP)
        citations.write_bib_file(BIB_FILE, bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)

        print_section("[7] SUMMARY & FILES", f_out)
        print_dual("Status                    : OK", f_out)
        print_dual(f"{len(output_files)} structure(s) written.", f_out)
        shown_files = output_files[:5]
        for of in shown_files:
            print_dual(f"Output file               : {of}", f_out)
        if len(output_files) > len(shown_files):
            print_dual(f"                            ... and {len(output_files) - len(shown_files)} more", f_out)
        print_dual(f"References                : {BIB_FILE}", f_out)
        if report_path:
            print_dual(f"Report                    : {report_path}", f_out)

        if f_out:
            f_out.close()

        if args.view:
            view_structure_interactive(view_atoms)

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
                        fail(str(e))
                indices = [resolve_nearest(structure, args.nearest, args.nearest_format, args.filter_species)]
        except ValueError as e:
            fail(str(e))

        indices_set = set(indices)
        selected = [(i + 1, atoms[i][0], atoms[i][1]) for i in sorted(indices_set)]
        for one_indexed, symbol, position in selected:
            print_dual(f"Selected site    : #{one_indexed} ({symbol}) at {position}", f_out)

        if args.type == "vacancy":
            atoms = [pair for i, pair in enumerate(atoms) if i not in indices_set]
        else:
            atoms = [(args.new_species if i in indices_set else sym, pos) for i, (sym, pos) in enumerate(atoms)]
            print_dual(f"Substitution     : -> {args.new_species}", f_out)
            if args.new_species not in species_list:
                species_list.append(args.new_species)
            species_meta = structure_io.ensure_species_id(species_meta, args.new_species)

    else:
        position = convert_position(args.position, args.position_format, structure)
        atoms = atoms + [(args.species, position)]
        print_dual(f"Interstitial     : {args.species} at {args.position} ({args.position_format})", f_out)
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
    print_dual(f"Output formula   : {pmg_after.composition.reduced_formula}", f_out)
    print_dual(f"Output atoms     : {len(atoms)}", f_out)

    ml_relax_info = None
    if args.ml_relax:
        print_section("[4] ML PRE-RELAXATION (MACE)", f_out)
        model_arg = args.custom_model if args.custom_model else args.model
        print_dual(f"Model            : {model_desc}", f_out)
        print_dual(f"Device           : {args.device}", f_out)
        print_dual(f"Cell relaxation  : {'yes (vacuum axes fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)
        try:
            calc = mace_relax.get_calculator(model_arg, device=args.device)
        except ValueError as e:
            fail(str(e))
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        frac_after = [site.frac_coords for site in pmg_after]
        vacuum_axes_after = kspace.detect_vacuum_axes(frac_after, pmg_after.lattice.matrix, VACUUM_GAP_ANG)
        t0 = time.time()
        pmg_after, e0, e1, converged, steps_used = _relax_structure(
            pmg_after, calc, vacuum_axes_after, relax_cell=args.ml_relax_cell, f_out=f_out)
        wall_time = time.time() - t0
        ml_relax_info = (converged, steps_used, e0, e1)

        print_dual(f"Steps used : {steps_used} "
                   f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
        print_dual(f"Wall time  : {wall_time:.1f} s", f_out)
        n_atoms_after = len(pmg_after)
        rows = [
            (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
              f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms_after:+.6f}/atom)"], None),
        ]
        print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

    print_section("[5] STRUCTURE VALIDATION, SYMMETRY & WRITING OUTPUT FILE", f_out)
    frac_after = [site.frac_coords for site in pmg_after]
    vacuum_axes_after = kspace.detect_vacuum_axes(frac_after, pmg_after.lattice.matrix, VACUUM_GAP_ANG)
    try:
        _validate_structure(pmg_after, vacuum_axes_after, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    before_info = core_symmetry.symmetry_summary(pmg_before, args.symprec, VACUUM_GAP_ANG)
    after_info = core_symmetry.symmetry_summary(pmg_after, args.symprec, VACUUM_GAP_ANG)
    if "Error" in before_info or "Error" in after_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
    else:
        properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
        rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                for prop in properties]
        print_table(["Property", "Before", "After"], rows, f_out)

    header_comment = [f"{args.type.capitalize()} defect introduced by stb-defect from {args.filename}."]
    if args.type == "vacancy":
        header_comment.append(f"Removed {len(selected)} atom(s): "
                               + ", ".join(f"#{i} ({s})" for i, s, _ in selected) + ".")
    elif args.type == "substitution":
        header_comment.append(f"Substituted {len(selected)} atom(s) with {args.new_species}: "
                               + ", ".join(f"#{i} ({s})" for i, s, _ in selected) + ".")
    else:
        header_comment.append(f"Added interstitial {args.species} at {args.position} ({args.position_format}).")
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV).")

    if ml_relax_info is not None:
        new_structure = structure_io.from_pymatgen(pmg_after, species_meta=species_meta)
    structure_io.write_fdf(new_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[6] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    if args.view:
        input_atoms = AseAtomsAdaptor.get_atoms(pmg_before)
        final_atoms = AseAtomsAdaptor.get_atoms(pmg_after)
        view_structure_interactive([input_atoms, final_atoms])


if __name__ == "__main__":
    main()
