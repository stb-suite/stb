#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import shutil
import argparse
import numpy as np
from pymatgen.core import Molecule
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo

REPORT_FILE = "oer_stage1.txt"
_DIPOLE_CORR_RE = re.compile(r'Slab\.DipoleCorrection\s+\S+', re.IGNORECASE)
_OH_BOND_LENGTH_ANG = 0.970  # experimental gas-phase OH radical equilibrium bond length --
                              # CG relaxation refines it, same spirit as --height and HER's
                              # _H2_BOND_LENGTH_ANG.


def reorient_vacuum_to_c(structure, vacuum_axes):
    """Same relabeling trick as stb-her's own reorient_vacuum_to_c (her.py) --
    duplicated, not imported (OER is self-contained, same policy as HER not
    depending on adsorb.py). Permutes lattice vectors so the single
    vacuum-padded axis becomes c (the convention AdsorbateSiteFinder
    assumes), restoring right-handedness if the permutation would flip it.
    A pure relabeling, not a rotation -- the physical structure is
    unchanged.
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
    """Same validation as stb-her's own resolve_slab_orientation
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
            f"[ERROR] stb-oer needs a slab/2D material with vacuum along exactly one axis "
            f"(OER adsorbs OH onto a surface); detected vacuum axis/axes: {detected}. A bulk 3D "
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


def force_dipole_correction(calc_text):
    """Substitutes/appends 'Slab.DipoleCorrection T' -- structurally
    required here, not just recommended: adsorbing OH on only one face of
    a slab breaks whatever inversion/mirror symmetry the clean slab had
    along the surface normal, giving the cell a net dipole moment along a
    PERIODIC direction -- without the dipole correction, that spurious
    periodic-image field contaminates the total energy (and therefore the
    site ranking Stage 3 does). Same substitute-or-append convention as
    stb-her's own force_dipole_correction.
    """
    new_text, count = _DIPOLE_CORR_RE.subn('Slab.DipoleCorrection   T', calc_text)
    if count == 0:
        new_text += "\nSlab.DipoleCorrection   T\n"
    return new_text


def write_site_folder(out_dir, ads_structure, calc_text, species_meta, pp_path):
    """Writes structure.fdf + calc.fdf + copied pseudos for one candidate
    site (or the clean-slab copy). Mirrors stb-her's own write_site_folder;
    also reused by stb-oerIntermediates' --o-strategy/--ooh-strategy
    search paths (same workflow, sibling script, direct import -- not the
    "never import a sibling workflow's module" HER/adsorb.py policy, which
    only applies ACROSS workflows).
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(ads_structure, species_meta=species_meta,
                                                coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = {site.specie.symbol for site in ads_structure}
    for sym in sorted(symbols):
        copy_pseudo(pp_path, sym, out_dir)
    return fdf_structure


def min_ads_slab_distance(pmg_structure, n_substrate):
    """Minimum periodic distance between any slab atom and the adsorbed
    OH group -- catches a --height that places OH unphysically close to
    (or inside) the slab. Same narrower-than-min_pairwise_distance
    reasoning as stb-her's own min_h_slab_distance (renamed here since
    the adsorbate is no longer a single H).
    """
    n_total = len(pmg_structure)
    if n_substrate == 0 or n_substrate >= n_total:
        return None
    dm = pmg_structure.distance_matrix
    return float(dm[:n_substrate, n_substrate:].min())


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 4: finds candidate OH-adsorption sites on a slab/2D "
        "material and writes one relaxation folder per site.", 'bold')}
Part of the OER (Oxygen Evolution Reaction) workflow -- computes the computational hydrogen
electrode (CHE, Rossmeisl et al. 2007; Man et al. 2011) descriptors Delta-G1..4, the theoretical
overpotential eta, and the potential-determining step, via the 3 adsorbed intermediates OH*, O*,
OOH*. This is Stage 1: it only searches for the OH* site (the other two intermediates are derived
or independently searched from the winning OH* site by stb-oerIntermediates, Stage 2). Every
symmetrically distinct site of --site-type is written (no --site-index/--all-sites choice) --
OER's whole point is finding the GLOBAL most stable site. Doesn't run SIESTA -- run each folder's
relaxation yourself, then use stb-oerIntermediates.

[NOTE] No O2 gas-phase reference is ever generated by this workflow (stb-oer/stb-oerIntermediates/
stb-oerRefs/stb-oerAnalysis) -- intentional, not an oversight. DFT/GGA badly describes the O2
triplet ground state (systematic error of at least ~0.3 eV, Sargeant et al. 2021), so
stb-oerAnalysis derives G(O2) from the experimental total reaction free energy instead. See
stb-oerAnalysis --help for details.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s -s slab.fdf -c calc.fdf -p dojo --site-type all\n"
               "  %(prog)s -s slab.fdf -c calc.fdf -p dojo --site-type ontop --both-sides\n"
    )

    parser.add_argument("-s", "--structure", type=str, required=True,
                         help="Input slab/2D structure.fdf (vacuum along one axis).")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template for the site relaxations (kgrid, basis, XC, "
                              "%%include structure.fdf, MD.TypeOfRun CG, etc.) -- copied into "
                              "every 'sites/site_*/' folder, with Slab.DipoleCorrection forced on "
                              "(structurally required, see --help of force_dipole_correction).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--site-type", choices=["ontop", "bridge", "hollow", "all"], default="all",
                         help="Which family/families of adsorption sites to search (default: "
                              "all -- OER wants the global minimum across site types).")
    parser.add_argument("--height", type=float, default=1.8,
                         help="OH adsorption distance above the surface, in Ang (default: 1.8 -- "
                              "slightly higher than a bare H atom's default, since OH sits "
                              "further from the surface at equilibrium; still just a CG-"
                              "relaxation starting guess).")
    parser.add_argument("--oh-bond-length", type=float, default=_OH_BOND_LENGTH_ANG,
                         help="O-H bond length of the adsorbing OH group, in Ang (default: "
                              f"{_OH_BOND_LENGTH_ANG}, the experimental gas-phase OH radical "
                              "equilibrium value). CG relaxation refines it further.")
    parser.add_argument("--both-sides", action="store_true",
                         help="Adsorb on both exposed faces (for a free-standing 2D material "
                              "with vacuum on both sides). Needs a symmetry operation mapping "
                              "top to bottom (same requirement as stb-her --both-sides).")
    parser.add_argument("--symprec", type=float, default=0.01,
                         help="Symmetry-reduction tolerance for site-finding (default: 0.01).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang (default: 10.0).")
    parser.add_argument("-O", "--output-dir", type=str, default="oer_study",
                         help="Root directory for the whole OER workflow (default: oer_study).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-oer {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OER WORKFLOW -- STAGE 1: ADSORPTION SITES (OH*)", 'bold'))
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
    slab_species_meta = structure_io.species_dict(structure)
    n_substrate = len(pmg_structure)

    with open(args.calc) as f:
        calc_text = force_dipole_correction(f.read())

    # O anchored at local origin (the site coordinate), H further along
    # local +z -- AdsorbateSiteFinder.add_adsorbate anchors whichever atom
    # sits at the molecule's most-negative local z to the surface site, so
    # this places O on the surface and H pointing outward (the physically
    # correct OH* orientation). Built by hand rather than via ASE's G2
    # 'OH' entry: G2's own internal geometry has H at the MORE negative
    # local z (H=-0.870, O=+0.109), which would anchor H to the surface
    # and point O away -- backwards. Same "hand-roll the geometry" choice
    # as her_refs.py::build_h2_structure, not adsorb.py::resolve_adsorbate's
    # G2-lookup pattern.
    oh_molecule = Molecule(["O", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, args.oh_bond_length]])
    finder = AdsorbateSiteFinder(pmg_structure)
    site_types = ["ontop", "bridge", "hollow"] if args.site_type == "all" else [args.site_type]

    output_root = args.output_dir
    clean_slab_dir = os.path.join(output_root, "clean_slab_source")
    sites_root = os.path.join(output_root, "sites")
    os.makedirs(sites_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, 'w') if report_path else None

    print_dual(f"{color_text('===== OER STAGE 1 REPORT (ADSORPTION SITES) =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Structure       : {args.structure}", f_out)
    print_dual(f"Calc template   : {args.calc}", f_out)
    print_dual(f"Site type       : {args.site_type}", f_out)
    print_dual(f"Height          : {args.height:.2f} Ang", f_out)
    print_dual(f"O-H bond length : {args.oh_bond_length:.3f} Ang", f_out)
    print_dual(f"Both sides      : {'yes' if args.both_sides else 'no'}", f_out)
    print_dual(f"symprec         : {args.symprec}", f_out)
    print_dual(f"Output root     : {output_root}", f_out)

    print_section('[1] CLEAN SLAB REFERENCE', f_out)
    os.makedirs(clean_slab_dir, exist_ok=True)
    shutil.copy(args.structure, os.path.join(clean_slab_dir, "structure.fdf"))
    print_dual(f"  {color_text('[OK]', 'green')} {clean_slab_dir}/structure.fdf (pristine geometry, "
                "copied as-is -- stb-oerRefs recomputes it single-point with the winning site's "
                "exact numerical settings, for consistency).", f_out)

    print_section('[2] ADSORPTION SITES (OH)', f_out)

    site_records = []  # (label, pmg_structure_with_OH)
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
        if len(site_types) != 1:
            print_dual(color_text(
                "[ERROR] --both-sides requires a concrete --site-type (not 'all').", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        try:
            both_structs = AdsorbateSiteFinder(slab_obj).adsorb_both_surfaces(
                oh_molecule, repeat=(1, 1, 1),
                find_args={"distance": args.height, "symm_reduce": args.symprec,
                           "positions": site_types})
        except RuntimeError as e:
            print_dual(color_text(
                f"[ERROR] Could not place OH on both faces: {e}. --both-sides needs a slab "
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
        for i, s in enumerate(both_structs, start=1):
            site_records.append((f"site_{i}_{site_types[0]}_bothsides", s))
    else:
        found = finder.find_adsorption_sites(distance=args.height, symm_reduce=args.symprec,
                                              positions=site_types)
        candidates = [(st, coord) for st in site_types for coord in found[st]]
        if not candidates:
            print_dual(color_text(f"[ERROR] No {'/'.join(site_types)} sites found.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        for i, (st, coord) in enumerate(candidates, start=1):
            ads_struct = finder.add_adsorbate(oh_molecule, coord)
            site_records.append((f"site_{i}_{st}", ads_struct))

    report_rows = []  # (label, dir)
    for label, ads_struct in site_records:
        min_dist = min_ads_slab_distance(ads_struct, n_substrate)
        if min_dist is not None and min_dist < 0.5:
            print_dual(color_text(
                f"  [WARNING] {label}: closest slab-OH distance is only {min_dist:.3f} Ang -- "
                "likely overlapping atoms (--height too small). Check before running SIESTA.",
                'yellow'), f_out)
        site_dir = os.path.join(sites_root, label)
        write_site_folder(site_dir, ads_struct, calc_text, slab_species_meta, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {site_dir}", f_out)
        report_rows.append((label, site_dir))

    print_section('[3] SUMMARY & NEXT STEPS', f_out)
    print_dual(f"{len(site_records)} site folder(s) written under '{sites_root}'.", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{sites_root}/site_*/' folder (a full relaxation).", f_out)
    print_dual(f"  2. Once they're done, run: stb-oerIntermediates --directory {output_root}", f_out)
    print_dual(color_text(
        "\n[NOTE] Spin polarization: OH is an intrinsically open-shell (doublet) radical, "
        "unlike neutral H2/atomic H -- 'Spin polarized' in your calc.fdf is strongly "
        "recommended, not just worth considering.", 'yellow'), f_out)
    print_dual(color_text(
        "[NOTE] Van der Waals: standard GGA misses dispersion forces, which can matter for "
        "physisorption-like distant sites -- a vdW correction (e.g. DFTD3) is recommended.",
        'yellow'), f_out)
    print_dual(color_text(
        "[NOTE] No O2 gas-phase reference is generated anywhere in this workflow -- "
        "intentional, see stb-oerAnalysis --help.", 'cyan'), f_out)
    print_dual(color_text(
        "[NOTE] This workflow assumes the adsorbate evolution mechanism (AEM: OH*/O*/OOH* "
        "bound at a single site) -- it does NOT model the lattice oxygen evolution mechanism "
        "(LOER, where lattice oxygen atoms themselves participate, common on some oxide/"
        "oxyhydroxide catalysts and linked to catalyst instability; see Exner, ChemCatChem "
        "2021, 'On the Lattice Oxygen Evolution Mechanism: Avoiding Pitfalls'). If your "
        "material is expected to favor LOER (e.g. a perovskite or Ru/Ir oxide operating at "
        "high anodic potential), treat this workflow's Delta-G/eta as an AEM-only estimate, "
        "not the full mechanistic picture.", 'cyan'), f_out)

    if f_out:
        f_out.write("\n# SITE_TABLE -- parsed by stb-oerIntermediates, do not reorder the columns\n")
        f_out.write(f"# {'label':<24}{'dir'}\n")
        for label, site_dir in report_rows:
            f_out.write(f"{label:<26}{site_dir}\n")
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Adsorption site folders ready for Stage 2 (stb-oerIntermediates).\n", 'bold'))


if __name__ == "__main__":
    main()
