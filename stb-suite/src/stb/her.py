#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
import numpy as np
from pymatgen.core import Molecule
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from stb.core import structure_io, kspace
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                           capture_library_noise)
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.symmetry import symmetry_summary, find_inequivalent_sites
from stb.core.adsorption_sites import (
    center_slab_in_vacuum, write_site_plot, cluster_candidate_coords,
    min_adsorbate_image_distance, _MIN_LATERAL_IMAGE_SEPARATION_ANG,
    parse_positions_file, write_positions_file, label_fragments, write_fragment_manifest,
    CONFIG_EXTRA_FILE, FIXED_CELL_BLOCK, DIPOLE_CORRECTION_BLOCK, SPIN_POLARIZED_BLOCK,
    VDW_CORRECTION_BLOCK,
)
from stb.neb import wrap_into_cell

REPORT_FILE = "her_stage1.txt"


def reorient_vacuum_to_c(structure, vacuum_axes):
    """Same relabeling trick as stb-adsorb's own reorient_vacuum_to_c
    (adsorb.py) -- duplicated, not imported (HER is self-contained, does
    not depend on adsorb.py, per the approved plan). Permutes lattice
    vectors so the single vacuum-padded axis becomes c (the convention
    AdsorbateSiteFinder assumes), restoring right-handedness if the
    permutation would flip it. A pure relabeling, not a rotation -- the
    physical structure is unchanged.
    """
    vacuum_axis = vacuum_axes.index(True)
    order = [i for i in range(3) if i != vacuum_axis] + [vacuum_axis]
    new_lattice = structure.lattice[order]
    if np.linalg.det(new_lattice) < 0:
        order[0], order[1] = order[1], order[0]
        new_lattice = structure.lattice[order]

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    new_frac_coords = frac_coords[:, order]
    new_atoms = [(sym, new_frac_coords[i]) for i, (sym, _pos) in enumerate(structure.atoms)]

    new_structure = structure_io.FdfStructure(
        lattice=new_lattice,
        lattice_constant=structure.lattice_constant,
        species=structure.species,
        species_meta=structure.species_meta,
        atoms=new_atoms,
        coord_format="fractional",
        raw_lines=[],
    )
    return new_structure, [False, False, True]


def resolve_slab_orientation(structure, vacuum_gap):
    """Same validation as stb-adsorb's own resolve_slab_orientation
    (duplicated, see reorient_vacuum_to_c's docstring): requires exactly
    one vacuum-padded axis (a slab/2D material, not bulk 3D or an
    isolated molecule/wire), reorienting it to c if needed.
    """
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
    if sum(vacuum_axes) != 1:
        detected = ', '.join(axis for axis, is_vac in zip('abc', vacuum_axes) if is_vac) or 'none'
        print(color_text(
            f"[ERROR] stb-her needs a slab/2D material with vacuum along exactly one axis "
            f"(HER adsorbs H onto a surface); detected vacuum axis/axes: {detected}. A bulk 3D "
            "structure (no vacuum) or a wire/molecule (vacuum on 2-3 axes) doesn't have a "
            "single well-defined surface.", 'red'))
        sys.exit(1)

    if not vacuum_axes[2]:
        old_axis = 'abc'[vacuum_axes.index(True)]
        structure, vacuum_axes = reorient_vacuum_to_c(structure, vacuum_axes)
        print(color_text(
            f"[INFO] Input structure's vacuum axis was '{old_axis}', not c -- relabeled the "
            "lattice vectors so vacuum is on c. Pure relabeling, not a rotation: every "
            "structure.fdf written below uses this relabeled a/b/c order.", 'yellow'))

    return structure


def write_site_folder(out_dir, ads_structure, calc_text, species_meta, pp_path, n_substrate):
    """Writes structure.fdf + calc.fdf + config_extra.fdf + copied pseudos
    for one candidate site. Mirrors stb-adsorb's own write_reference_folder
    (adsorb.py, menu 4.8) -- and the same config_extra.fdf sidecar
    convention stb-raman/stb-ir (menu 4.11/4.12) use for their own forced
    directives -- rewritten locally (see module docstring for why the
    surrounding site-finding/geometry logic isn't imported), but reusing
    the SAME shared block constants/CONFIG_EXTRA_FILE name from
    core.adsorption_sites rather than defining HER-specific duplicates,
    since 'keep the cell fixed' and 'force Slab.DipoleCorrection' are not
    HER-specific concerns.

    FIXED_CELL_BLOCK, DIPOLE_CORRECTION_BLOCK, SPIN_POLARIZED_BLOCK and
    VDW_CORRECTION_BLOCK are all unconditional here (not opt-in, unlike
    stb-adsorb's own --force-spin-style flags): FIXED_CELL_BLOCK because
    this is one independent site relaxation, not something that should be
    moving the underlying slab's own lattice; DIPOLE_CORRECTION_BLOCK
    because adsorbing H on only one face of a slab breaks whatever
    inversion/mirror symmetry the clean slab had along the surface normal,
    giving the cell a net dipole moment along a PERIODIC direction --
    without the correction, that spurious periodic-image field
    contaminates the total energy (and therefore the site ranking Stage 2
    does); SPIN_POLARIZED_BLOCK because H adsorption intermediates can
    have radical/open-shell character depending on the surface (most
    simply, an odd total valence-electron count, which a spin-restricted
    calculation cannot represent correctly at all) -- mandatory, not
    merely advisory, since a spin-restricted run of a genuinely open-shell
    site silently produces a wrong energy with no warning; costs nothing
    for a genuinely closed-shell site (converges to zero moment);
    VDW_CORRECTION_BLOCK (Grimme DFT-D3) because standard GGA functionals
    miss dispersion forces entirely, which matter for physisorption-like
    distant sites -- mandatory here too (used to be only an advisory
    [NOTE] in the report), since a distant/weakly-bound site's energy
    ranking against a closer, more strongly-bound one is exactly the kind
    of comparison dispersion can flip. `%include config_extra.fdf` is
    prepended to the UNTOUCHED user --calc template text
    (structure_io.prepend_include) rather than editing it in place --
    SIESTA's fdf reader is first-occurrence-wins, so this overrides
    whatever the template itself might already set.

    Fragment labeling ('<symbol>_slab'/'<symbol>_ads') is also reused
    directly from core.adsorption_sites (label_fragments/
    write_fragment_manifest): every atom is relabeled depending on whether
    it came from the input slab (index < n_substrate) or is the adsorbed H
    (index >= n_substrate), so a slab that already carries its own H (e.g.
    a passivated dangling-bond termination from stb-passivate/stb-slab
    --passivate) is never confused with the freshly-adsorbed H after a
    write_fdf/read_fdf round trip (write_fdf groups atoms by species, not
    physical origin). A fragment_manifest.json is written alongside
    structure.fdf so stb-herRefs' BSSE ghost-triad construction can trust
    the label split instead of re-deriving it from a fragile atom count.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(ads_structure, species_meta=species_meta,
                                                coord_format="fractional")
    fdf_structure = label_fragments(fdf_structure, n_substrate)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
        f.write(FIXED_CELL_BLOCK + DIPOLE_CORRECTION_BLOCK + SPIN_POLARIZED_BLOCK
                 + VDW_CORRECTION_BLOCK)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(structure_io.prepend_include(calc_text, CONFIG_EXTRA_FILE))
    for label in fdf_structure.species:
        real_symbol = structure_io.real_element(label, fdf_structure.species_meta)
        copy_pseudo(pp_path, real_symbol, out_dir, dest_label=label)
    write_fragment_manifest(out_dir, fdf_structure, n_substrate)
    return fdf_structure


def min_h_slab_distance(pmg_structure, n_substrate):
    """Minimum periodic distance between any slab atom and the adsorbed H
    -- catches a --height that places H unphysically close to (or inside)
    the slab. Same narrower-than-min_pairwise_distance reasoning as
    stb-adsorb's own min_adsorbate_slab_distance.
    """
    n_total = len(pmg_structure)
    if n_substrate == 0 or n_substrate >= n_total:
        return None
    dm = pmg_structure.distance_matrix
    return float(dm[:n_substrate, n_substrate:].min())


def print_symmetry_section(pmg_structure, symprec, vacuum_gap, f_out, library_warnings):
    """Prints [1] SLAB SYMMETRY: space/point group, crystal system, Hall
    symbol, layer group (core.symmetry.symmetry_summary), plus a table of
    symmetrically distinct substrate atoms (Wyckoff letter + multiplicity,
    core.symmetry.find_inequivalent_sites) -- the user explicitly asked for
    symmetry AND sites to be reported in full, unlike stb-adsorb's own
    Stage 1, which only ever reports a raw-vs-reduced candidate COUNT, never
    the actual space/point group or Wyckoff labels. Fails soft (prints a
    [NOTE], never sys.exit) on any symmetry-detection failure -- same
    'return/print a note, caller decides' contract as core/siesta_log.py's
    parsers.
    """
    with capture_library_noise(library_warnings, "pymatgen SpacegroupAnalyzer (symmetry_summary)"):
        summary = symmetry_summary(pmg_structure, symprec=symprec, vacuum_gap_ang=vacuum_gap)
    if "Error" in summary:
        print_dual(color_text(f"  [NOTE] Symmetry analysis failed: {summary['Error']} -- "
                    "continuing without a symmetry report.", 'yellow'), f_out)
        return
    print_table(["Property", "Value"], [([k, v], None) for k, v in summary.items()], f_out)

    try:
        with capture_library_noise(library_warnings,
                                    "pymatgen SpacegroupAnalyzer (find_inequivalent_sites)"):
            sites, space_group_label = find_inequivalent_sites(pmg_structure, symprec)
    except Exception as e:
        print_dual(color_text(f"  [NOTE] Inequivalent-site analysis failed: {e}.", 'yellow'), f_out)
        return
    print_dual(f"\n  {len(sites)} symmetrically distinct substrate atom(s) (space group "
                f"{space_group_label}):", f_out)
    site_rows = [([str(idx), pmg_structure[idx].specie.symbol, wyckoff, str(mult)], None)
                 for idx, wyckoff, mult in sites]
    print_table(["Atom idx (0-based)", "Element", "Wyckoff", "Multiplicity"], site_rows, f_out)
    print_dual(color_text("  Note: these are the DISTINCT ATOMS of the clean slab itself "
        "(before H is added) -- not the H-adsorption candidate coordinates, found separately "
        "in [3] ADSORPTION SITES: FINDING & COUNT below.", 'yellow'), f_out)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 3: finds candidate H-adsorption sites on a slab/2D "
        "material and writes one relaxation folder per site.", 'bold')}
Part of the HER (Hydrogen Evolution Reaction) workflow -- computes the computational hydrogen
electrode (CHE, Norskov et al.) descriptor Delta-G_H*. Unlike the general-purpose stb-adsorb
(any adsorbate, multi-adsorbate, ML ranking, height sweeps), this tool is deliberately narrower:
the adsorbate is always a single H atom, and EVERY symmetrically distinct site of --site-type is
written (no --site-index/--all-sites choice) -- HER's whole point is finding the GLOBAL most
stable site, so there is no reason to write only some of them. Doesn't run SIESTA -- run each
folder's relaxation yourself, then use stb-herRefs.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s -s slab.fdf -c calc.fdf -p dojo --site-type all\n"
               "  %(prog)s -s slab.fdf -c calc.fdf -p dojo --site-type ontop --both-sides\n"
               "  %(prog)s -s slab.fdf -c calc.fdf -p dojo --position 1.23 0.71 --height 1.6\n"
               "  %(prog)s -s slab.fdf -c calc.fdf -p dojo --positions-file sites/site_positions.dat\n"
    )

    parser.add_argument("-s", "--structure", type=str, required=True,
                         help="Input slab/2D structure.fdf (vacuum along one axis).")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template for the site relaxations (kgrid, basis, XC, "
                              "%%include structure.fdf, MD.TypeOfRun CG, etc.) -- untouched and "
                              "copied into every 'sites/site_*/' folder, with a config_extra.fdf "
                              "sidecar %%include'd on top forcing a fixed cell, "
                              "Slab.DipoleCorrection, Spin polarized and DFTD3 (all mandatory, "
                              "not opt-in -- see write_site_folder's docstring).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--site-type", choices=["ontop", "bridge", "hollow", "all"], default="all",
                         help="Which family/families of adsorption sites to search (default: "
                              "all -- HER wants the global minimum across site types). Ignored "
                              "with --position/--positions-file.")
    parser.add_argument("--height", type=float, default=1.5,
                         help="H adsorption distance above the surface, in Ang (default: 1.5). "
                              "Also used as the offset along the surface normal for "
                              "--position/--positions-file.")
    parser.add_argument("--both-sides", action="store_true",
                         help="Adsorb on both exposed faces (for a free-standing 2D material "
                              "with vacuum on both sides). Needs a symmetry operation mapping "
                              "top to bottom (same requirement as stb-adsorb --both-sides). "
                              "Requires a concrete --site-type (not 'all'); mutually exclusive "
                              "with --position/--positions-file.")
    parser.add_argument("--symprec", type=float, default=0.01,
                         help="Symmetry-reduction tolerance for site-finding (default: 0.01). "
                              "Also used for the [1] SLAB SYMMETRY report. Ignored with "
                              "--position/--positions-file.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang (default: 10.0).")
    parser.add_argument("--position", type=float, nargs=2, metavar=("X", "Y"), default=None,
                         help="Manual Cartesian (X, Y) position in Ang for a single H site, "
                              "overriding automatic site-finding entirely (--site-type/--symprec "
                              "are ignored when this is given). Mutually exclusive with "
                              "--positions-file and --both-sides.")
    parser.add_argument("--positions-file", type=str, default=None, metavar="PATH",
                         help="'.dat' file with one or more H-adsorption site positions, as "
                              "in-plane fractional (a, b) coordinates relative to the file's own "
                              "LATTICE_A/LATTICE_B/LATTICE_C lines ('number type frac_a frac_b' "
                              "per site line; '#' comments and blank lines ignored) -- used "
                              "INSTEAD of automatic site-finding (--site-type/--symprec are "
                              "ignored). EVERY position in the file gets its own site_N/ folder "
                              "-- unlike stb-adsorb, stb-her has no --site-index/--all-sites "
                              "choice to make (HER always wants every distinct site). stb-her "
                              "always writes 'sites/site_positions.dat' with the lattice and "
                              "positions actually used this run (auto-found or manual alike), "
                              "in the same format, so it can be inspected, hand-edited, and fed "
                              "back here for a future run. Mutually exclusive with --position "
                              "and --both-sides.")
    parser.add_argument("-O", "--output-dir", type=str, default="her_study",
                         help="Root directory for the whole HER workflow (default: her_study).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-her {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.position is not None and args.positions_file:
        parser.error("--position and --positions-file are mutually exclusive.")
    if args.both_sides and args.position is not None:
        parser.error("--both-sides and --position are mutually exclusive.")
    if args.both_sides and args.positions_file:
        parser.error("--both-sides and --positions-file are not supported together yet.")
    if args.both_sides and args.site_type == "all":
        parser.error("--both-sides requires a concrete --site-type (not 'all').")
    if args.positions_file and not os.path.isfile(args.positions_file):
        parser.error(f"--positions-file file not found: {args.positions_file}")

    print("\n" + color_text("HER WORKFLOW -- STAGE 1: ADSORPTION SITES", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
        sys.exit(1)
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    structure = structure_io.read_fdf(args.structure)
    structure = resolve_slab_orientation(structure, args.vacuum_gap)
    pmg_structure = structure_io.to_pymatgen(structure)
    pmg_structure, slab_recentered, recenter_shift_ang = center_slab_in_vacuum(pmg_structure)
    slab_species_meta = structure_io.species_dict(structure)
    n_substrate = len(pmg_structure)

    with open(args.calc) as f:
        calc_text = f.read()
    calc_text, include_status, old_include_name = structure_io.fix_structure_include(calc_text)

    h_molecule = Molecule(["H"], [[0.0, 0.0, 0.0]])
    finder = AdsorbateSiteFinder(pmg_structure)
    site_types = ["ontop", "bridge", "hollow"] if args.site_type == "all" else [args.site_type]

    output_root = args.output_dir
    clean_slab_dir = os.path.join(output_root, "clean_slab_source")
    sites_root = os.path.join(output_root, "sites")
    os.makedirs(sites_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, 'w') if report_path else None
    library_warnings = []  # collected via capture_library_noise, reported in [6]

    print_dual(f"{color_text('===== HER STAGE 1 REPORT (ADSORPTION SITES) =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Structure       : {args.structure}", f_out)
    print_dual(f"Calc template   : {args.calc}", f_out)
    if include_status == "fixed":
        print_dual(color_text(
            f"[NOTE] --calc's structure include was '%include {old_include_name}' -- "
            "auto-corrected to '%include structure.fdf' in every generated site_*/calc.fdf "
            "(every folder's structure is always written under that exact name).", 'yellow'), f_out)
    elif include_status == "not_found":
        print_dual(color_text(
            "[NOTE] Could not find a '%include ...struct....fdf'-style line in --calc to "
            "auto-correct -- make sure it '%include structure.fdf' itself so SIESTA picks up "
            "the generated geometry.", 'yellow'), f_out)
    print_dual(f"Pseudo dir      : {args.pseudo_dir or '(none)'}", f_out)
    print_dual(f"Output root     : {output_root}", f_out)
    if args.position is not None:
        print_dual(f"Site selection  : manual --position ({args.position[0]:.4f}, "
                    f"{args.position[1]:.4f}) Ang", f_out)
    elif args.positions_file:
        print_dual(f"Site selection  : manual --positions-file ({args.positions_file})", f_out)
    else:
        print_dual(f"Site selection  : automatic, site type '{args.site_type}'"
                    + (", both faces" if args.both_sides else ""), f_out)
    print_dual(f"Site type       : {args.site_type}", f_out)
    print_dual(f"Height          : {args.height:.2f} Ang", f_out)
    print_dual(f"Both sides      : {'yes' if args.both_sides else 'no'}", f_out)
    print_dual(f"symprec         : {args.symprec}", f_out)
    print_dual(f"Vacuum gap      : {args.vacuum_gap:.1f} Ang", f_out)
    print_dual("Fixed/dipole/spin/vdW: yes -> MD.VariableCell false + Slab.DipoleCorrection + "
                "Spin polarized + DFTD3, via config_extra.fdf (forced in every site_*/)", f_out)
    print_dual(f"Save report     : {'yes -> ' + report_path if args.save_report else 'no'}", f_out)

    print_section('[1] SLAB SYMMETRY', f_out)
    print_symmetry_section(pmg_structure, args.symprec, args.vacuum_gap, f_out, library_warnings)

    print_section('[2] CLEAN SLAB REFERENCE', f_out)
    if slab_recentered:
        print_dual(color_text(
            f"  [INFO] Slab was near a cell boundary along c (uneven vacuum split) -- shifted by "
            f"{recenter_shift_ang:+.3f} Ang to center it in the vacuum gap. This avoids a real "
            "risk: an atom starting near frac z=0 can drift slightly negative during SIESTA "
            "relaxation and wrap around to frac z~1 (the top of the cell) instead, which reads "
            "as an enormous, spurious displacement to anything comparing this relaxed geometry "
            "against another snapshot of the same atom later.", 'yellow'), f_out)
    else:
        print_dual("  Slab already centered in the vacuum gap (no boundary-wrap risk).", f_out)
    os.makedirs(clean_slab_dir, exist_ok=True)
    clean_fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=slab_species_meta,
                                                       coord_format="fractional")
    structure_io.write_fdf(clean_fdf_structure, os.path.join(clean_slab_dir, "structure.fdf"))
    print_dual(f"  {color_text('[OK]', 'green')} {clean_slab_dir}/structure.fdf (pristine geometry, "
                "reflecting the recentering above if any -- stb-herRefs recomputes it single-point "
                "with the winning site's exact numerical settings, for consistency).", f_out)

    print_section('[3] ADSORPTION SITES: FINDING & COUNT', f_out)

    lateral_dist = min_adsorbate_image_distance(np.zeros((1, 3)), pmg_structure.lattice.matrix)
    if lateral_dist < _MIN_LATERAL_IMAGE_SEPARATION_ANG:
        print_dual(color_text(
            f"  [WARNING] This slab's own periodic images are only {lateral_dist:.2f} Ang apart "
            f"in-plane (recommended >= {_MIN_LATERAL_IMAGE_SEPARATION_ANG:.0f} Ang) -- a "
            "symmetry-equivalent H site can end up spuriously close to its own periodic copy. "
            "Use a bigger in-plane supercell (stb-supercell) if so.", 'yellow'), f_out)
    else:
        print_dual(f"  Lateral periodic-image separation: {lateral_dist:.2f} Ang (>= "
                    f"{_MIN_LATERAL_IMAGE_SEPARATION_ANG:.0f} Ang recommended) -- OK.", f_out)

    plot_path = os.path.join(sites_root, "adsorption_sites.png")
    positions_file_out = os.path.join(sites_root, "site_positions.dat")

    def _save_positions_file(cands):
        entries = [(i + 1, label, coord[0], coord[1]) for i, (label, coord) in enumerate(cands)]
        write_positions_file(positions_file_out, pmg_structure.lattice.matrix, entries)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {positions_file_out} ({len(entries)} "
                    "position(s), fractional -- feed it back via --positions-file to reuse/"
                    "hand-pick among exactly these sites in a future run)", f_out)

    def _save_site_plot(numbered_sites=None):
        with capture_library_noise(library_warnings, "pymatgen plot_slab"):
            write_site_plot(pmg_structure, plot_path, show=False, numbered_sites=numbered_sites)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {plot_path}"
                    + (" (each number is the site_N this run writes it as)"
                       if numbered_sites is not None else
                       " (sanity check only -- pymatgen's own defaults)"), f_out)

    candidates = None  # [(label_or_site_type, cart_coord), ...] -- one entry per H site to write
    both_structs = None

    if args.both_sides:
        from pymatgen.core.surface import Slab
        slab_obj = Slab(
            lattice=pmg_structure.lattice,
            species=[site.specie for site in pmg_structure],
            coords=pmg_structure.frac_coords,
            miller_index=(0, 0, 1),
            oriented_unit_cell=pmg_structure,
            shift=0.0,
            scale_factor=np.eye(3),
            reorient_lattice=False,
        )
        with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
            probe_raw = finder.find_adsorption_sites(distance=args.height, symm_reduce=0,
                                                       positions=site_types)
            probe_reduced = finder.find_adsorption_sites(distance=args.height,
                                                           symm_reduce=args.symprec,
                                                           positions=site_types)
        n_raw, n_reduced = len(probe_raw[site_types[0]]), len(probe_reduced[site_types[0]])
        print_table(["Site type", "Raw candidates", "After symmetry reduction", "Reduced away"],
                    [([site_types[0], str(n_raw), str(n_reduced), str(n_raw - n_reduced)], None)],
                    f_out)
        print_dual("  (top face only -- --both-sides mirrors each selected site onto the bottom "
                    "face too, doubling the final count.)", f_out)
        _save_site_plot()
        try:
            with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder (both-sides)"):
                both_structs = AdsorbateSiteFinder(slab_obj).adsorb_both_surfaces(
                    h_molecule, repeat=(1, 1, 1),
                    find_args={"distance": args.height, "symm_reduce": args.symprec,
                               "positions": site_types})
        except RuntimeError as e:
            print_dual(color_text(
                f"[ERROR] Could not place H on both faces: {e}. --both-sides needs a slab "
                "with a symmetry operation mapping top to bottom (a free-standing 2D "
                "material, or a centrosymmetric slab).", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        if not both_structs:
            print_dual(color_text(f"[ERROR] No {site_types[0]} sites found.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        print_dual(f"  {len(both_structs)} configuration(s) after mirroring onto both faces.", f_out)

    elif args.position is not None:
        print_dual("  Manual override (--position): no site-finding performed -- exactly 1 "
                    "site.", f_out)
        x, y = args.position
        z_ref = float(np.max(pmg_structure.cart_coords[:, 2]))
        mvec = np.asarray(finder.mvec, dtype=float)
        coord = np.array([x, y, z_ref]) + args.height * mvec
        if abs(mvec[2]) < 0.999:
            angle_deg = np.degrees(np.arccos(np.clip(abs(mvec[2]), -1.0, 1.0)))
            print_dual(color_text(
                f"  [INFO] Surface normal is not parallel to Cartesian z (tilted "
                f"{angle_deg:.2f} deg) -- height is measured along the true normal "
                f"{tuple(round(float(v), 4) for v in mvec)}, not a plain z-offset.", 'yellow'),
                f_out)
        candidates = [("manual", coord)]
        _save_site_plot([(1, coord)])
        _save_positions_file(candidates)

    elif args.positions_file:
        print_dual("  Manual override (--positions-file): reading candidate positions directly "
                    "from the file.", f_out)
        try:
            file_lattice, file_positions = parse_positions_file(args.positions_file)
        except ValueError as e:
            print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        print_dual(f"  {len(file_positions)} position(s) read from {args.positions_file}.", f_out)
        if not np.allclose(file_lattice, pmg_structure.lattice.matrix, atol=1e-3):
            print_dual(color_text(
                "  [WARNING] --positions-file's own LATTICE_A/LATTICE_B/LATTICE_C differ from "
                "this run's actual structure lattice -- positions were converted using the "
                "FILE's lattice, so they may not land where expected. Regenerate the file from "
                "this exact structure if unsure.", 'yellow'), f_out)
        z_ref = float(np.max(pmg_structure.cart_coords[:, 2]))
        mvec = np.asarray(finder.mvec, dtype=float)
        candidates = [(st, np.array([x, y, z_ref]) + args.height * mvec)
                      for _number, st, x, y in file_positions]
        if abs(mvec[2]) < 0.999:
            angle_deg = np.degrees(np.arccos(np.clip(abs(mvec[2]), -1.0, 1.0)))
            print_dual(color_text(
                f"  [INFO] Surface normal is not parallel to Cartesian z (tilted "
                f"{angle_deg:.2f} deg) -- height is measured along the true normal "
                f"{tuple(round(float(v), 4) for v in mvec)}, not a plain z-offset.", 'yellow'),
                f_out)
        if not candidates:
            print_dual(color_text("[ERROR] --positions-file has no usable positions.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        _save_site_plot([(i + 1, coord) for i, (_st, coord) in enumerate(candidates)])
        _save_positions_file(candidates)
        print_dual(f"  {len(candidates)} position(s) from --positions-file -- every one will be "
                    "written below.", f_out)

    else:
        rows, raw_counts, reduced_counts = [], {}, {}
        with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
            for st in site_types:
                raw = finder.find_adsorption_sites(distance=args.height, symm_reduce=0,
                                                    positions=[st])
                reduced = finder.find_adsorption_sites(distance=args.height,
                                                        symm_reduce=args.symprec, positions=[st])
                raw_counts[st], reduced_counts[st] = len(raw[st]), len(reduced[st])
        for st in site_types:
            rows.append(([st, str(raw_counts[st]), str(reduced_counts[st]),
                          str(raw_counts[st] - reduced_counts[st])], None))
        total_raw, total_reduced = sum(raw_counts.values()), sum(reduced_counts.values())
        if len(site_types) > 1:
            rows.append((["TOTAL", str(total_raw), str(total_reduced),
                          str(total_raw - total_reduced)], 'cyan'))
        print_table(["Site type", "Raw candidates", "After symm. reduction",
                     "Reduced away (equivalent)"], rows, f_out)

        with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
            found = finder.find_adsorption_sites(distance=args.height, symm_reduce=args.symprec,
                                                  positions=site_types)
            anchor = next((c for st in site_types for c in found[st]), None)
            if anchor is not None:
                for st in site_types:
                    found[st] = cluster_candidate_coords(found[st], pmg_structure.lattice.matrix,
                                                          reference=anchor)
            candidates = [(st, coord) for st in site_types for coord in found[st]]

        if not candidates:
            print_dual(color_text(f"[ERROR] No {'/'.join(site_types)} sites found.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        _save_site_plot([(i + 1, coord) for i, (_st, coord) in enumerate(candidates)])
        _save_positions_file(candidates)
        print_dual(f"  {len(candidates)} symmetrically distinct '{args.site_type}' site(s) "
                    "found -- every one will be written below (HER has no --site-index/"
                    "--all-sites choice: it wants the global minimum across ALL distinct "
                    "sites).", f_out)

    if candidates is not None:
        print_dual(f"  Configuration count : {len(candidates)} H site folder(s), plus "
                    "clean_slab_source/ (no separate H reference here -- see stb-herRefs "
                    "Stage 2 for the H2 gas-phase/BSSE references).", f_out)

    print_section('[4] WRITING SITE FOLDERS', f_out)
    print_dual("  Each site's structure.fdf labels atoms by fragment (species "
                "'<symbol>_slab'/'<symbol>_ads') so the adsorbed H can never be confused with "
                "a slab atom of the same element (e.g. a slab already passivated with H) -- see "
                "fragment_manifest.json in each folder. The forced fixed-cell/dipole-correction/"
                "spin-polarized/DFTD3 directives live in each folder's own config_extra.fdf, "
                "%include'd on top of your --calc template (never edited in place).", f_out)
    site_records = []  # (label, pmg_structure_with_H)
    if args.both_sides:
        for i, s in enumerate(both_structs, start=1):
            site_records.append((f"site_{i}_{site_types[0]}_bothsides", s))
    else:
        for i, (st, coord) in enumerate(candidates, start=1):
            ads_struct = finder.add_adsorbate(h_molecule, coord)
            site_records.append((f"site_{i}_{st}", ads_struct))

    report_rows = []  # (label, dir)
    for label, ads_struct in site_records:
        ads_struct = wrap_into_cell(ads_struct)
        min_dist = min_h_slab_distance(ads_struct, n_substrate)
        if min_dist is not None and min_dist < 0.7:
            print_dual(color_text(
                f"  [WARNING] {label}: closest slab-H distance is only {min_dist:.3f} Ang -- "
                "likely overlapping atoms (--height too small, or a bad --position/"
                "--positions-file entry). Check before running SIESTA.", 'yellow'), f_out)
        site_dir = os.path.join(sites_root, label)
        write_site_folder(site_dir, ads_struct, calc_text, slab_species_meta, args.pseudo_dir,
                           n_substrate)
        print_dual(f"  {color_text('[OK]', 'green')} {site_dir}", f_out)
        report_rows.append((label, site_dir))

    print_section('[5] SUMMARY & NEXT STEPS', f_out)
    print_dual(f"{len(site_records)} site folder(s) written under '{sites_root}'.", f_out)
    if os.path.isfile(plot_path):
        print_dual(f"Site plot            : {plot_path}", f_out)
    if os.path.isfile(positions_file_out):
        print_dual(f"Site positions       : {positions_file_out} (feed back via "
                    "--positions-file)", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{sites_root}/site_*/' folder (a full relaxation).", f_out)
    print_dual(f"  2. Once they're done, run: stb-herRefs --directory {output_root}", f_out)

    print_section('[6] LIBRARY WARNINGS', f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (pymatgen/spglib) during this run -- "
            "collected here instead of interleaved with the report above; harmless in almost "
            "every case, but worth a look if a section above looks suspicious.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.write("\n# SITE_TABLE -- for your own reference (NOT parsed by stb-herRefs, which "
                     "scans 'sites/site_*/' folders directly via os.scandir)\n")
        f_out.write(f"# {'label':<24}{'dir'}\n")
        for label, site_dir in report_rows:
            f_out.write(f"{label:<26}{site_dir}\n")
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Adsorption site folders ready for Stage 2 (stb-herRefs).\n", 'bold'))


if __name__ == "__main__":
    main()
