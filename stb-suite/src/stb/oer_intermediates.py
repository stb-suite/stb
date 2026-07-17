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
import argparse
import numpy as np
from pymatgen.core import Molecule
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.pseudopotentials import resolve_pseudo_source, link_pseudo
from stb.core.siesta_log import get_free_energy, get_outcell, report_quality_diagnostics
from stb.core.deps import require_mace

REPORT_FILE = "oer_stage2.txt"
_OOH_OO_BOND_ANG = 1.45   # illustrative peroxo-like O-O starting bond length -- NOT pinned to a
                          # specific catalyst's literature value, same "CG relaxation refines it"
                          # spirit as --height. See build_ooh_structure's docstring.
_OOH_OH_BOND_ANG = 0.970  # experimental gas-phase OH bond length, reused for the new terminal H.
_OOH_BEND_DEG = 100.0     # H2O2-like O-O-H angle -- illustrative starting geometry, not fitted.

_DIPOLE_CORR_RE = re.compile(r'Slab\.DipoleCorrection\s+\S+', re.IGNORECASE)
_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_NUMCGSTEPS_RE = re.compile(r'MD\.NumCGsteps\s+\d+', re.IGNORECASE)
_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)
_RELAXED_COORDS_RE = re.compile(r'outcoor:\s*Relaxed atomic coordinates\s*\(fractional\)', re.IGNORECASE)


def force_dipole_correction(calc_text):
    """Duplicated from oer.py (self-contained stage -- each stage in this
    suite re-reads persisted files/duplicates helpers rather than
    importing a sibling stage's module, same policy her_refs.py/
    her_analysis.py already follow for HER). Structurally required: ANY
    one-sided adsorbate (OH, bare O, or OOH alike) breaks the clean
    slab's inversion/mirror symmetry along the surface normal.
    """
    new_text, count = _DIPOLE_CORR_RE.subn('Slab.DipoleCorrection   T', calc_text)
    if count == 0:
        new_text += "\nSlab.DipoleCorrection   T\n"
    return new_text


def force_relaxation(calc_text, max_steps=200):
    """Duplicated from her_refs.py::force_relaxation. Forces
    'MD.TypeOfRun CG' + 'MD.NumCGsteps <max_steps>' -- every geometry
    this module writes (derived O*/OOH* starting guesses, or fresh
    --strategy search candidates) MUST reach its own relaxed minimum
    before any BSSE/ZPE/final-energy calculation trusts it; explicit is
    safer than silently inheriting whatever (or nothing) a borrowed
    calc.fdf template happens to have.
    """
    new_text, count = _RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _NUMCGSTEPS_RE.subn(f'MD.NumCGsteps         {max_steps}', new_text)
    if count == 0:
        new_text += f"MD.NumCGsteps         {max_steps}\n"
    return new_text


def force_system_label(calc_text, label):
    """Duplicated from her_refs.py::force_system_label -- every derived
    folder needs its own distinct SystemLabel so SIESTA's own per-run
    files never collide between folders run from the same working tree.
    """
    new_text, count = _LABEL_RE.subn(f'SystemLabel {label}', calc_text)
    if count == 0:
        new_text += f"\nSystemLabel {label}\n"
    return new_text


def find_winning_site(sites_root, out_file, f_out):
    """Duplicated from her_refs.py::find_winning_site. Scans every
    'site_*/' (or candidate) folder's FreeEng and returns (winning_dir,
    winning_energy, all_results).
    """
    site_dirs = sorted(
        d.path for d in os.scandir(sites_root) if d.is_dir() and d.name.startswith("site_")
    )
    if not site_dirs:
        print_dual(color_text(f"[ERROR] No 'site_*' folders found in '{sites_root}'.", 'red'), f_out)
        sys.exit(1)

    results = []
    best_dir, best_energy = None, float('inf')
    for d in site_dirs:
        out_path = os.path.join(d, out_file)
        energy = get_free_energy(out_path)
        results.append((os.path.basename(d), energy))
        if energy is not None and energy < best_energy:
            best_energy = energy
            best_dir = d
    if best_dir is None:
        print_dual(color_text("[ERROR] No site folder has a readable energy -- did SIESTA finish?",
                               'red'), f_out)
        sys.exit(1)
    return best_dir, best_energy, results


def read_relaxed_structure(out_path, template):
    """Duplicated from her_refs.py::read_relaxed_structure. Reads the
    LAST 'outcoor: Relaxed atomic coordinates (fractional)' block from a
    SIESTA .out file -- SIESTA never rewrites structure.fdf in place, so
    every downstream derived structure (O*, OOH*) needs to start from the
    geometry the OH* relaxation actually reached. Returns None if the
    block isn't found.
    """
    try:
        with open(out_path, errors='ignore') as f:
            lines = f.readlines()
    except OSError:
        return None

    start_idx = None
    for i in range(len(lines) - 1, -1, -1):
        if _RELAXED_COORDS_RE.search(lines[i]):
            start_idx = i + 1
            break
    if start_idx is None:
        return None

    id_to_symbol = {str(info['id']): sym for sym, info in template.species_meta.items()}
    n_atoms = len(template.atoms)
    new_atoms = []
    try:
        for i in range(n_atoms):
            parts = lines[start_idx + i].split()
            frac = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
            symbol = id_to_symbol[parts[3]]
            new_atoms.append((symbol, frac))
    except (IndexError, ValueError, KeyError):
        return None

    lattice = get_outcell(out_path)
    if lattice is None:
        lattice = template.lattice

    return structure_io.FdfStructure(
        lattice=lattice, lattice_constant=template.lattice_constant,
        species=template.species, species_meta=template.species_meta,
        atoms=new_atoms, coord_format="fractional", raw_lines=[],
    )


def _perpendicular_axis(v):
    """Unit vector perpendicular to v, used as the Rodrigues rotation
    axis in build_ooh_structure/build_ooh_molecule. Picks the global x
    axis as the cross-product reference instead of z whenever v is
    nearly parallel to z, so the cross product never degenerates.
    """
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(v, ref)) > 0.95:
        ref = np.array([1.0, 0.0, 0.0])
    perp = np.cross(v, ref)
    return perp / np.linalg.norm(perp)


def _rodrigues_rotate(v, axis, angle_rad):
    """Rotates vector v about the unit vector axis by angle_rad
    (Rodrigues' rotation formula)."""
    axis = axis / np.linalg.norm(axis)
    return (v * np.cos(angle_rad) + np.cross(axis, v) * np.sin(angle_rad)
            + axis * np.dot(axis, v) * (1.0 - np.cos(angle_rad)))


def build_o_structure(relaxed_oh, h_index):
    """O* = OH*'s relaxed geometry minus the H atom -- every substrate
    atom AND the adsorbed O keep exactly the position OH*'s own
    relaxation converged to. Same single-index removal her_refs.py's
    remove_atom already does for HER's '03_slab_deformed'.
    """
    new_atoms = [a for i, a in enumerate(relaxed_oh.atoms) if i != h_index]
    present = {sym for sym, _ in new_atoms}
    species_meta = {k: v for k, v in relaxed_oh.species_meta.items() if k in present}
    species = [s for s in relaxed_oh.species if s in present]
    return structure_io.FdfStructure(
        lattice=relaxed_oh.lattice, lattice_constant=relaxed_oh.lattice_constant,
        species=species, species_meta=species_meta, atoms=new_atoms,
        coord_format=relaxed_oh.coord_format, raw_lines=[],
    )


def build_ooh_structure(relaxed_oh, o_index, h_index,
                         oo_bond_ang=_OOH_OO_BOND_ANG, oh_bond_ang=_OOH_OH_BOND_ANG,
                         bend_deg=_OOH_BEND_DEG):
    """Derives a STARTING OOH* geometry from the relaxed OH* site: every
    substrate atom and O1 keep exactly the position OH*'s relaxation
    converged to; the original H is dropped (it re-forms as the new
    terminal group's H instead) and O2 is appended continuing outward
    along the O1->H1 direction (approximating the local "away from the
    surface" direction after relaxation) at `oo_bond_ang`, then a new H
    on O2 bent off the O2->O1 axis by `bend_deg` (Rodrigues rotation
    about an axis perpendicular to O1-O2, matching a peroxide-like O-O-H
    angle) at `oh_bond_ang`.

    This is a documented STARTING geometry for further CG relaxation
    (this module always writes it into a relaxation folder, never a
    single-point one -- see main()), NOT a transition-state guess, same
    spirit as stb-oer's --height default. Adsorbate atoms are always O1
    (unchanged), O2 (new), H (new) appended last, preserving the "slab
    first, adsorbate last" invariant later BSSE/ZPE code in stb-oerRefs
    assumes.
    """
    o_symbol, o_frac = relaxed_oh.atoms[o_index]
    h_symbol, h_frac = relaxed_oh.atoms[h_index]
    lattice = relaxed_oh.lattice
    inv_lattice = np.linalg.inv(lattice)
    o_cart = o_frac @ lattice
    h_cart = h_frac @ lattice

    axis1 = h_cart - o_cart
    axis1 = axis1 / np.linalg.norm(axis1)
    o2_cart = o_cart + axis1 * oo_bond_ang

    rot_axis = _perpendicular_axis(axis1)
    bend_dir = _rodrigues_rotate(-axis1, rot_axis, np.radians(bend_deg))
    hnew_cart = o2_cart + bend_dir * oh_bond_ang

    o2_frac = o2_cart @ inv_lattice
    hnew_frac = hnew_cart @ inv_lattice

    new_atoms = [a for i, a in enumerate(relaxed_oh.atoms) if i != h_index]
    new_atoms.append(("O", o2_frac))
    new_atoms.append(("H", hnew_frac))

    return structure_io.FdfStructure(
        lattice=lattice, lattice_constant=relaxed_oh.lattice_constant,
        species=list(relaxed_oh.species), species_meta=dict(relaxed_oh.species_meta),
        atoms=new_atoms, coord_format=relaxed_oh.coord_format, raw_lines=[],
    )


def build_ooh_molecule(oo_bond_ang=_OOH_OO_BOND_ANG, oh_bond_ang=_OOH_OH_BOND_ANG,
                        bend_deg=_OOH_BEND_DEG):
    """Standalone OOH group for --ooh-strategy search's fresh site-finding
    (as opposed to build_ooh_structure's derive-from-OH* path): O1
    anchored at the local origin (bonds to the surface, same anchor
    convention as oer.py's oh_molecule -- AdsorbateSiteFinder.add_adsorbate
    anchors the most-negative-local-z atom), O2 along local +z at
    `oo_bond_ang`, H bent off O2 by `bend_deg` -- same illustrative
    starting geometry as build_ooh_structure's, just built fresh instead
    of continuing an existing OH* orientation.
    """
    o1 = np.array([0.0, 0.0, 0.0])
    axis1 = np.array([0.0, 0.0, 1.0])
    o2 = o1 + axis1 * oo_bond_ang
    rot_axis = _perpendicular_axis(axis1)
    bend_dir = _rodrigues_rotate(-axis1, rot_axis, np.radians(bend_deg))
    hnew = o2 + bend_dir * oh_bond_ang
    return Molecule(["O", "O", "H"], [o1, o2, hnew])


def write_relax_folder(out_dir, fdf_structure, calc_text, pp_path):
    """Writes structure.fdf + calc.fdf + linked pseudos for one
    relaxation folder (derived single candidate, or one --strategy
    search candidate). Same shape as oer.py's write_site_folder /
    her_refs.py's write_folder, duplicated locally (self-contained
    stage).
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = sorted({sym for sym, _ in fdf_structure.atoms})
    for sym in symbols:
        link_pseudo(pp_path, sym, out_dir)


def ml_prerelax_adsorbate(fdf_structure, n_substrate, model, device, fmax):
    """Positions-only MACE-MP-0 relax of ONLY the adsorbate atoms (every
    index >= n_substrate), with the substrate held fixed via ase's
    FixAtoms -- same pattern as stb-adsorb's own --ml-rank
    (mace_relax.get_calculator + relax + FixAtoms(indices=range(n_substrate))),
    just operating on this module's FdfStructure representation instead
    of a pymatgen Structure directly. A fast classical-potential
    pre-screen to improve O*/OOH*'s starting geometry -- especially
    --strategy derived's hand-built guess (build_ooh_structure's O-O bond
    length/bend angle are explicitly illustrative, not literature-fitted)
    -- before the real, much more expensive SIESTA CG relaxation this
    module always still writes a folder for. NOT a substitute for that
    real relaxation, just a better-informed starting point for it.
    Callers must call core.deps.require_mace() themselves first (see
    main()), matching every other MACE consumer in this suite.
    """
    from ase.constraints import FixAtoms
    from pymatgen.io.ase import AseAtomsAdaptor
    from stb.core import mace_relax

    pmg_structure = structure_io.to_pymatgen(fdf_structure)
    ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    ase_atoms.set_constraint(FixAtoms(indices=list(range(n_substrate))))
    calc = mace_relax.get_calculator(model=model, device=device)
    mace_relax.relax(ase_atoms, calc, fmax=fmax, max_steps=200)
    relaxed_pmg = AseAtomsAdaptor.get_structure(ase_atoms)
    species_meta = structure_io.species_dict(fdf_structure)
    return structure_io.from_pymatgen(relaxed_pmg, species_meta=species_meta, coord_format="fractional")


def search_intermediate_sites(pmg_clean, clean_species_meta, ads_molecule, site_type, height,
                               symprec, label_prefix, out_root, base_calc_text, pp_path,
                               ml_prerelax, ml_model, ml_device, ml_fmax, f_out):
    """--strategy search: an independent AdsorbateSiteFinder sweep for
    `ads_molecule` on the CLEAN slab -- same site-finding machinery as
    Stage 1's own OH* search (oer.py's main()), duplicated here rather
    than imported (self-contained stage, see module docstring). Writes
    one CG-relaxation folder per symmetrically distinct site under
    `out_root/<label_prefix>_candidates/site_N_type/`. If `ml_prerelax`,
    every candidate's adsorbate geometry is MACE-MP-0 pre-relaxed
    (substrate fixed) before being written out, same idea as
    ml_prerelax_adsorbate but applied per-candidate here since each one
    starts from AdsorbateSiteFinder's own fixed-height, non-relaxed
    placement. Returns the number of folders written.
    """
    finder = AdsorbateSiteFinder(pmg_clean)
    site_types = ["ontop", "bridge", "hollow"] if site_type == "all" else [site_type]
    found = finder.find_adsorption_sites(distance=height, symm_reduce=symprec, positions=site_types)
    candidates = [(st, coord) for st in site_types for coord in found[st]]
    if not candidates:
        print_dual(color_text(
            f"[ERROR] No {'/'.join(site_types)} sites found for the {label_prefix} search.", 'red'), f_out)
        sys.exit(1)

    candidates_root = os.path.join(out_root, f"{label_prefix}_candidates")
    os.makedirs(candidates_root, exist_ok=True)
    n_substrate = len(pmg_clean)

    mace_calc = None
    if ml_prerelax:
        from stb.core import mace_relax
        from ase.constraints import FixAtoms
        from pymatgen.io.ase import AseAtomsAdaptor
        mace_calc = mace_relax.get_calculator(model=ml_model, device=ml_device)
        print_dual(f"  {color_text('ML pre-relax:', 'cyan')} relaxing each {label_prefix.upper()} "
                    "candidate's adsorbate atoms with MACE-MP-0 (substrate fixed) before writing "
                    "it out ...", f_out)

    n_written = 0
    for i, (st, coord) in enumerate(candidates, start=1):
        ads_struct = finder.add_adsorbate(ads_molecule, coord)
        if mace_calc is not None:
            ase_atoms = AseAtomsAdaptor.get_atoms(ads_struct)
            ase_atoms.set_constraint(FixAtoms(indices=list(range(n_substrate))))
            mace_relax.relax(ase_atoms, mace_calc, fmax=ml_fmax, max_steps=200)
            ads_struct = AseAtomsAdaptor.get_structure(ase_atoms)
        label = f"site_{i}_{st}"
        site_dir = os.path.join(candidates_root, label)
        fdf_structure = structure_io.from_pymatgen(
            ads_struct, species_meta=clean_species_meta, coord_format="fractional")
        calc_text = force_system_label(
            force_relaxation(force_dipole_correction(base_calc_text)), f"oer_{label_prefix}_{label}")
        write_relax_folder(site_dir, fdf_structure, calc_text, pp_path)
        print_dual(f"  {color_text('[OK]', 'green')} {site_dir}", f_out)
        n_written += 1
    return n_written


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 4: derives (or independently searches for) the O* and "
        "OOH* intermediates from the winning OH* site, and writes their own CG-relaxation "
        "folder(s).", 'bold')}
Picks the lowest-FreeEng OH* site from Stage 1's 'sites/site_*/', reads its RELAXED geometry, then
for O* and OOH* independently (--o-strategy/--ooh-strategy):

  - 'derived' (default): builds ONE starting geometry from the winning OH* site (O* = OH* minus
    H; OOH* = OH* plus a second O-H group in a chemically plausible but illustrative starting
    orientation -- see build_ooh_structure --help/docstring) and writes ONE CG-relaxation folder
    ('intermediates/o_star/', 'intermediates/ooh_star/'). Cheap: no extra site search.

  - 'search': an independent AdsorbateSiteFinder sweep for that one adsorbate on the CLEAN slab
    (same machinery as Stage 1's own OH* search), writing MULTIPLE CG-relaxation candidate
    folders ('intermediates/o_candidates/site_*/', 'intermediates/ooh_candidates/site_*/') --
    stb-oerRefs (Stage 3) picks the winner once you've relaxed them all. More rigorous (the true
    global-minimum O*/OOH* site need not be the same as OH*'s), at the cost of N extra SIESTA
    relaxations per intermediate searched this way.

Every folder written here MUST be relaxed via SIESTA (MD.TypeOfRun CG is forced on) before
running stb-oerRefs -- unlike a single-point evaluation, O*/OOH*'s starting geometry (especially
in 'derived' mode) is only a reasonable guess, not already at equilibrium.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory oer_study\n"
               "  %(prog)s --directory oer_study --ooh-strategy search --ooh-site-type ontop\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="oer_study",
                         help="Root directory written by stb-oer (default: oer_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each site folder (default: calc.out).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--o-strategy", choices=["derived", "search"], default="derived",
                         help="How to obtain O*'s starting geometry (default: derived).")
    parser.add_argument("--ooh-strategy", choices=["derived", "search"], default="derived",
                         help="How to obtain OOH*'s starting geometry (default: derived).")
    parser.add_argument("--oo-bond-length", type=float, default=_OOH_OO_BOND_ANG,
                         help="Illustrative starting O-O bond length for OOH*, in Ang (default: "
                              f"{_OOH_OO_BOND_ANG}). Refined by CG relaxation.")
    parser.add_argument("--ooh-oh-bond-length", type=float, default=_OOH_OH_BOND_ANG,
                         help="Illustrative starting O-H bond length for OOH*'s new terminal H, "
                              f"in Ang (default: {_OOH_OH_BOND_ANG}).")
    parser.add_argument("--ooh-bend-deg", type=float, default=_OOH_BEND_DEG,
                         help="Illustrative starting O-O-H bend angle for OOH*, in degrees "
                              f"(default: {_OOH_BEND_DEG}, H2O2-like).")
    parser.add_argument("--o-site-type", choices=["ontop", "bridge", "hollow", "all"], default="all",
                         help="--o-strategy search only: which site type(s) to search (default: all).")
    parser.add_argument("--o-height", type=float, default=1.5,
                         help="--o-strategy search only: O adsorption height in Ang (default: 1.5).")
    parser.add_argument("--o-symprec", type=float, default=0.01,
                         help="--o-strategy search only: symmetry-reduction tolerance (default: 0.01).")
    parser.add_argument("--ooh-site-type", choices=["ontop", "bridge", "hollow", "all"], default="all",
                         help="--ooh-strategy search only: which site type(s) to search (default: all).")
    parser.add_argument("--ooh-height", type=float, default=2.0,
                         help="--ooh-strategy search only: OOH adsorption height in Ang (default: 2.0).")
    parser.add_argument("--ooh-symprec", type=float, default=0.01,
                         help="--ooh-strategy search only: symmetry-reduction tolerance (default: 0.01).")
    parser.add_argument("--ml-prerelax", action="store_true",
                         help="Relax O*/OOH*'s adsorbate atoms (positions only, substrate fixed) "
                              "with MACE-MP-0 before writing the CG-relaxation folder(s) -- same "
                              "idea as stb-adsorb --ml-rank. Mainly useful for --strategy derived's "
                              "hand-built starting geometry (the O-O bond length/bend angle are "
                              "illustrative, not literature-fitted), but also applies to --strategy "
                              "search's AdsorbateSiteFinder placements. A better-informed starting "
                              "point for the real SIESTA relaxation, NOT a substitute for it. Needs "
                              "the optional 'ml' extra.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small",
                         help="MACE-MP-0 model size, with --ml-prerelax (default: small).")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu",
                         help="Device for --ml-prerelax (default: cpu).")
    parser.add_argument("--ml-fmax", type=float, default=0.05,
                         help="Force convergence threshold in eV/Ang, with --ml-prerelax "
                              "(default: 0.05).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-oerIntermediates {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OER WORKFLOW -- STAGE 2: O*/OOH* INTERMEDIATES", 'bold'))
    print("-" * 60)

    output_root = args.directory
    sites_root = os.path.join(output_root, "sites")
    clean_slab_source = os.path.join(output_root, "clean_slab_source", "structure.fdf")
    if not os.path.isdir(sites_root) or not os.path.isfile(clean_slab_source):
        print(color_text(f"[ERROR] '{sites_root}' or '{clean_slab_source}' not found -- run "
                          "stb-oer (Stage 1) first.", 'red'))
        sys.exit(1)

    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    if args.ml_prerelax:
        require_mace()

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== OER STAGE 2 REPORT (O*/OOH* INTERMEDIATES) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Directory       : {output_root}", f_out)
        print_dual(f"O* strategy     : {args.o_strategy}", f_out)
        print_dual(f"OOH* strategy   : {args.ooh_strategy}", f_out)
        print_dual(f"ML pre-relax    : {'yes' if args.ml_prerelax else 'no'}"
                    + (f" (model={args.ml_model}, device={args.ml_device})" if args.ml_prerelax else ""),
                    f_out)

        print_dual(f"\n{color_text('[1] WINNING OH* SITE', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        winning_dir, winning_energy, all_results = find_winning_site(sites_root, args.file, f_out)
        for label, energy in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(sites_root, label) == winning_dir else ""
            energy_str = f"{energy:.6f} eV" if energy is not None else "(no energy)"
            print_dual(f"  {label:<28}{energy_str}{marker}", f_out)
        print_dual(f"Winning OH* site : {os.path.basename(winning_dir)} ({winning_energy:.6f} eV)", f_out)
        report_quality_diagnostics(os.path.basename(winning_dir),
                                    os.path.join(winning_dir, args.file), 0.05, f_out)

        winning_template = structure_io.read_fdf(os.path.join(winning_dir, "structure.fdf"))
        relaxed_oh = read_relaxed_structure(os.path.join(winning_dir, args.file), winning_template)
        if relaxed_oh is None:
            print_dual(color_text(
                f"[ERROR] Could not read relaxed coordinates from '{winning_dir}/{args.file}' -- "
                "did the relaxation finish?", 'red'), f_out)
            sys.exit(1)
        n_total = len(relaxed_oh.atoms)
        o_index = n_total - 2  # OH* atoms are always appended [O, H] last by stb-oer
        h_index = n_total - 1

        with open(winning_dir + "/calc.fdf") as f:
            site_calc_text = f.read()

        clean_fdf_structure = structure_io.read_fdf(clean_slab_source)
        pmg_clean = structure_io.to_pymatgen(clean_fdf_structure)
        clean_species_meta = structure_io.species_dict(clean_fdf_structure)

        print_dual(f"\n{color_text('[2] O* GEOMETRY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        if args.o_strategy == "derived":
            o_dir = os.path.join(output_root, "intermediates", "o_star")
            o_structure = build_o_structure(relaxed_oh, h_index)
            prerelax_note = ""
            if args.ml_prerelax:
                print_dual(f"  {color_text('ML pre-relax:', 'cyan')} relaxing O*'s adsorbate atom "
                            "with MACE-MP-0 (substrate fixed) ...", f_out)
                o_structure = ml_prerelax_adsorbate(o_structure, o_index, args.ml_model,
                                                     args.ml_device, args.ml_fmax)
                prerelax_note = ", ML pre-relaxed"
            o_calc = force_system_label(force_relaxation(site_calc_text), "oer_o_star")
            write_relax_folder(o_dir, o_structure, o_calc, args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} {o_dir} (derived from the winning OH* "
                        f"site: H removed{prerelax_note})", f_out)
        else:
            o_molecule = Molecule(["O"], [[0.0, 0.0, 0.0]])
            n_written = search_intermediate_sites(
                pmg_clean, clean_species_meta, o_molecule, args.o_site_type, args.o_height,
                args.o_symprec, "o", os.path.join(output_root, "intermediates"), site_calc_text,
                args.pseudo_dir, args.ml_prerelax, args.ml_model, args.ml_device, args.ml_fmax, f_out)
            print_dual(f"{n_written} O* candidate folder(s) written -- run SIESTA in all of them, "
                        "stb-oerRefs will pick the winner.", f_out)

        print_dual(f"\n{color_text('[3] OOH* GEOMETRY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        if args.ooh_strategy == "derived":
            ooh_dir = os.path.join(output_root, "intermediates", "ooh_star")
            ooh_structure = build_ooh_structure(relaxed_oh, o_index, h_index,
                                                 args.oo_bond_length, args.ooh_oh_bond_length,
                                                 args.ooh_bend_deg)
            prerelax_note = ""
            if args.ml_prerelax:
                print_dual(f"  {color_text('ML pre-relax:', 'cyan')} relaxing OOH*'s adsorbate "
                            "atoms with MACE-MP-0 (substrate fixed) ...", f_out)
                ooh_structure = ml_prerelax_adsorbate(ooh_structure, o_index, args.ml_model,
                                                       args.ml_device, args.ml_fmax)
                prerelax_note = ", ML pre-relaxed"
            ooh_calc = force_system_label(force_relaxation(site_calc_text), "oer_ooh_star")
            write_relax_folder(ooh_dir, ooh_structure, ooh_calc, args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} {ooh_dir} (derived from the winning OH* "
                        f"site: +O at {args.oo_bond_length:.3f} Ang, +H bent "
                        f"{args.ooh_bend_deg:.1f} deg -- illustrative starting geometry, refined "
                        f"by CG relaxation{prerelax_note})", f_out)
        else:
            ooh_molecule = build_ooh_molecule(args.oo_bond_length, args.ooh_oh_bond_length,
                                               args.ooh_bend_deg)
            n_written = search_intermediate_sites(
                pmg_clean, clean_species_meta, ooh_molecule, args.ooh_site_type, args.ooh_height,
                args.ooh_symprec, "ooh", os.path.join(output_root, "intermediates"), site_calc_text,
                args.pseudo_dir, args.ml_prerelax, args.ml_model, args.ml_device, args.ml_fmax, f_out)
            print_dual(f"{n_written} OOH* candidate folder(s) written -- run SIESTA in all of "
                        "them, stb-oerRefs will pick the winner.", f_out)

        print_dual(f"\n{color_text('[4] SUMMARY & NEXT STEPS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual("  1. Run SIESTA in every folder written above.", f_out)
        print_dual(f"  2. Once they're done, run: stb-oerRefs --directory {output_root}", f_out)

        f_out.write("\nO* strategy     : " + args.o_strategy + "\n")
        f_out.write("OOH* strategy   : " + args.ooh_strategy + "\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Intermediate geometries ready for Stage 3 (stb-oerRefs).\n", 'bold'))


if __name__ == "__main__":
    main()
