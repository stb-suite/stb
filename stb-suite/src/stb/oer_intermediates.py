#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import argparse
import numpy as np
import ase.io as ase_io
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.siesta_log import get_free_energy, get_outcell, report_quality_diagnostics
from stb.core.deps import require_mace
from stb.core.adsorption_sites import (
    ADSORBATE_LABEL_SUFFIX, CONFIG_EXTRA_FILE, FIXED_CELL_BLOCK, DIPOLE_CORRECTION_BLOCK,
    SPIN_POLARIZED_BLOCK, VDW_CORRECTION_BLOCK, deduplicate_orientations,
)
from stb.core.bsse import strip_config_extra_include

REPORT_FILE = "oer_stage2.txt"
_OOH_OO_BOND_ANG = 1.45   # illustrative peroxo-like O-O starting bond length -- NOT pinned to a
                          # specific catalyst's literature value, same "CG relaxation refines it"
                          # spirit as --height. See build_ooh_structure's docstring.
_OOH_OH_BOND_ANG = 0.970  # experimental gas-phase OH bond length, reused for the new terminal H.
_OOH_BEND_DEG = 100.0     # H2O2-like O-O-H angle -- illustrative starting geometry, not fitted.

# config_extra.fdf content for the O*/OOH* folders this module writes: these
# are independent SITE RELAXATIONS, same physical situation as stb-oer's own Stage 1 sites
# (a fixed cell, forced Slab.DipoleCorrection/Spin polarized/DFTD3 -- see
# oer.py's write_site_folder docstring for why each is mandatory), not
# single-point references (contrast HER's Stage 2, her_refs.py, which
# writes mostly single-point folders and therefore needs its own
# SINGLE_POINT_BLOCK-bearing constant instead).
_CONFIG_EXTRA = FIXED_CELL_BLOCK + DIPOLE_CORRECTION_BLOCK + SPIN_POLARIZED_BLOCK + VDW_CORRECTION_BLOCK

_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_MD_STEPS_RE = re.compile(r'MD\.Steps\s+\d+', re.IGNORECASE)
_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)
_RELAXED_COORDS_RE = re.compile(r'outcoor:\s*Relaxed atomic coordinates\s*\(fractional\)', re.IGNORECASE)


def _bare_ase_atoms(fdf_structure):
    """Converts fdf_structure (which may carry a Stage-1/2 fragment label
    like 'O_ads', or a stack of them) to a bare-element ase.Atoms for
    extended-XYZ export -- OVITO/VMD need real chemical symbols, and
    structure_io.to_pymatgen (which AseAtomsAdaptor.get_atoms would need)
    cannot parse a fragment label as an element. Used for O*/OOH*'s
    single-frame trajectories, where no bare-label pymatgen
    Structure exists at any point (build_o_structure/build_ooh_structure
    build directly on relaxed_oh's own already-fragment-labeled atoms).
    """
    from ase import Atoms
    symbols = [structure_io.real_element(sym, fdf_structure.species_meta)
               for sym, _pos in fdf_structure.atoms]
    frac = np.array([pos for _, pos in fdf_structure.atoms])
    cart = frac @ fdf_structure.lattice if fdf_structure.coord_format == "fractional" else frac
    return Atoms(symbols=symbols, positions=cart, cell=fdf_structure.lattice, pbc=True)


def force_relaxation(calc_text, max_steps=200):
    """Duplicated from her_refs.py::force_relaxation. Forces
    'MD.TypeOfRun CG' + 'MD.Steps <max_steps>' -- every O*/OOH* starting
    guess this module writes MUST reach its own relaxed minimum
    before any BSSE/ZPE/final-energy calculation trusts it; explicit is
    safer than silently inheriting whatever (or nothing) a borrowed
    calc.fdf template happens to have. Uses MD.Steps, not the older
    MD.NumCGsteps spelling (deprecated, no longer honored by current
    SIESTA).
    """
    new_text, count = _RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _MD_STEPS_RE.subn(f'MD.Steps              {max_steps}', new_text)
    if count == 0:
        new_text += f"MD.Steps              {max_steps}\n"
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

    `relaxed_oh.atoms[o_index]`'s own label (e.g. 'O_ads', a Stage-1
    fragment label -- see oer.py's write_site_folder) is reused verbatim
    for O1; the two brand-new atoms (O2, H-new) are given that SAME
    '<element>_ads' fragment suffix (not a bare 'O'/'H') so every
    adsorbate atom in the returned structure is consistently
    fragment-labeled, matching label_fragments' own convention, with a
    fresh species id/Z registered in species_meta for each.
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

    species_meta = dict(relaxed_oh.species_meta)

    def _new_ads_label(element_symbol, z):
        label = f"{element_symbol}{ADSORBATE_LABEL_SUFFIX}"
        if label not in species_meta:
            used_ids = {str(info['id']) for info in species_meta.values()}
            next_id = 1
            while str(next_id) in used_ids:
                next_id += 1
            species_meta[label] = {'id': str(next_id), 'Z': z}
        return label

    o2_label = _new_ads_label("O", 8)
    hnew_label = _new_ads_label("H", 1)

    new_atoms = [a for i, a in enumerate(relaxed_oh.atoms) if i != h_index]
    new_atoms.append((o2_label, o2_frac))
    new_atoms.append((hnew_label, hnew_frac))

    species = list(dict.fromkeys(sym for sym, _ in new_atoms))
    return structure_io.FdfStructure(
        lattice=lattice, lattice_constant=relaxed_oh.lattice_constant,
        species=species, species_meta=species_meta,
        atoms=new_atoms, coord_format=relaxed_oh.coord_format, raw_lines=[],
    )


def write_relax_folder(out_dir, fdf_structure, calc_text, pp_path):
    """Writes structure.fdf + calc.fdf + config_extra.fdf + copied
    pseudos for one derived-intermediate relaxation folder (O* or OOH*,
    built from the winning OH* site). Same config_extra.fdf sidecar
    convention as oer.py's write_site_folder / her_refs.py's write_folder
    (4.8/4.11/4.12/4.13 model), duplicated locally (self-contained
    stage): `_CONFIG_EXTRA` (fixed cell + mandatory Slab.DipoleCorrection/
    Spin polarized/DFTD3, same combination as a Stage-1 site -- these ARE
    independent site relaxations) is written to config_extra.fdf, and
    `%include config_extra.fdf` is prepended to the UNTOUCHED calc_text
    (structure_io.prepend_include) rather than editing directives into it
    in place. Pseudopotential copying resolves the real element behind
    any Stage-1 fragment label ('<real>_slab'/'<real>_ads') via
    structure_io.real_element (Z-based, robust to any suffix), with
    dest_label preserving the fragment-suffixed filename SIESTA's own
    ChemicalSpeciesLabel block expects.
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
        f.write(_CONFIG_EXTRA)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(structure_io.prepend_include(calc_text, CONFIG_EXTRA_FILE))
    present_labels = sorted({sym for sym, _ in fdf_structure.atoms})
    for label in present_labels:
        real_symbol = structure_io.real_element(label, fdf_structure.species_meta)
        copy_pseudo(pp_path, real_symbol, out_dir, dest_label=label)


def _mace_relax_adsorbate_scored(fdf_structure, n_substrate, model, device, fmax):
    """Positions-only MACE-MP-0 relax of ONLY the adsorbate atoms (every
    index >= n_substrate), with the substrate held fixed via ase's
    FixAtoms -- same pattern as stb-adsorb's own --ml-rank
    (mace_relax.get_calculator + relax + FixAtoms(indices=range(n_substrate))),
    just operating on this module's FdfStructure representation instead
    of a pymatgen Structure directly. Returns (energy_eV, relaxed_fdf_structure)
    -- the energy is what lets sample_ooh_orientations rank candidates;
    ml_prerelax_adsorbate (below) is a thin wrapper for callers that only
    want the relaxed geometry, discarding the energy.

    `fdf_structure` may carry Stage-1 fragment labels ('<real>_slab'/
    '<real>_ads') -- structure_io.to_pymatgen/AseAtomsAdaptor cannot parse
    those as real elements (pymatgen's Structure constructor tries to
    resolve each species string as an actual element/species), so the
    MACE relax runs on a BARE-element copy instead, and only the relaxed
    POSITIONS (same atom count/order -- FixAtoms never adds, removes, or
    reorders atoms) are copied back onto the ORIGINAL fragment-labeled
    species/species_meta afterward. No fragment-identity information is
    lost: only Cartesian positions change during this relax. Callers must
    call core.deps.require_mace() themselves first (see main()), matching
    every other MACE consumer in this suite.
    """
    from ase.constraints import FixAtoms
    from stb.core import mace_relax

    bare_species_meta = {}
    for orig_sym, _pos in fdf_structure.atoms:
        real_sym = structure_io.real_element(orig_sym, fdf_structure.species_meta)
        if real_sym not in bare_species_meta:
            bare_species_meta[real_sym] = {'id': str(len(bare_species_meta) + 1),
                                            'Z': abs(fdf_structure.species_meta[orig_sym]['Z'])}
    bare_atoms = [(structure_io.real_element(sym, fdf_structure.species_meta), pos)
                  for sym, pos in fdf_structure.atoms]
    bare_structure = structure_io.FdfStructure(
        lattice=fdf_structure.lattice, lattice_constant=fdf_structure.lattice_constant,
        species=list(dict.fromkeys(sym for sym, _ in bare_atoms)), species_meta=bare_species_meta,
        atoms=bare_atoms, coord_format=fdf_structure.coord_format, raw_lines=[],
    )

    pmg_structure = structure_io.to_pymatgen(bare_structure)
    ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    ase_atoms.set_constraint(FixAtoms(indices=list(range(n_substrate))))
    # dispersion=True: level-of-theory-matching with the real SIESTA
    # relaxation this is a starting point for -- config_extra.fdf forces
    # DFTD3 unconditionally on every folder this module writes (_CONFIG_EXTRA),
    # so the MACE pre-relax should include dispersion too, same rationale as
    # stb-adsorb's own --ml-rank/--ml-prerelax calculators.
    calc = mace_relax.get_calculator(model=model, device=device, dispersion=True)
    mace_relax.relax(ase_atoms, calc, fmax=fmax, max_steps=200)
    energy = ase_atoms.get_potential_energy()
    ase_atoms.wrap()
    relaxed_frac = AseAtomsAdaptor.get_structure(ase_atoms).frac_coords

    new_atoms = [(orig_sym, relaxed_frac[i]) for i, (orig_sym, _pos) in enumerate(fdf_structure.atoms)]
    relaxed_structure = structure_io.FdfStructure(
        lattice=fdf_structure.lattice, lattice_constant=fdf_structure.lattice_constant,
        species=fdf_structure.species, species_meta=fdf_structure.species_meta,
        atoms=new_atoms, coord_format=fdf_structure.coord_format, raw_lines=[],
    )
    return energy, relaxed_structure


def ml_prerelax_adsorbate(fdf_structure, n_substrate, model, device, fmax):
    """Thin wrapper around _mace_relax_adsorbate_scored for the single-
    geometry (no orientation sampling) pre-relax path -- a fast classical
    -potential pre-screen to improve O*/OOH*'s starting geometry -- especially
    OOH*'s hand-built guess (build_ooh_structure's O-O bond length/bend
    angle are explicitly illustrative, not literature-fitted) -- before the
    real, much more expensive SIESTA CG relaxation this module always still
    writes a folder for. NOT a substitute for that real relaxation, just a
    better-informed starting point for it. Discards the MACE energy (only
    sample_ooh_orientations' ranking needs it).
    """
    _energy, relaxed_structure = _mace_relax_adsorbate_scored(fdf_structure, n_substrate, model, device, fmax)
    return relaxed_structure


_OOH_ORIENTATION_MAX_POLAR_DEG = 70.0  # illustrative cap, not literature-pinned (see
                                        # _sample_ooh_directions' docstring for why it can't be 180)


def _sample_ooh_directions(n_polar, n_azimuthal, max_polar_deg=_OOH_ORIENTATION_MAX_POLAR_DEG):
    """Systematically samples n_polar x n_azimuthal unit directions for the
    O1->O2 bond, restricted to the HEMISPHERE pointing away from the
    surface (local +z, the module's own vacuum-normal convention) -- polar
    angle from 0 (straight up, build_ooh_structure's own default) to
    `max_polar_deg`, never all the way to 180 deg.

    core.adsorption_sites.generate_systematic_orientations' FULL-sphere
    Fibonacci-lattice sampling is the wrong tool here: it's built for
    AdsorbateSiteFinder's TOUCHDOWN convention, where a translation step
    afterward guarantees the molecule's own most-negative-z atom lands
    exactly AT the surface regardless of which way the rest of it
    initially points. build_ooh_structure_oriented has no such
    translation -- a sampled direction is used DIRECTLY as an offset from
    O1's own already-fixed real position, so a downward-pointing sample
    (verified live: it drove O2 to within 0.32 Ang of a substrate atom,
    producing a ~-134800 eV MACE 'energy' from the resulting force
    explosion) would send O2 straight into the substrate. Restricting to
    this hemisphere is the fix.

    n_polar=1 (default) returns [(0,0,1)] alone (straight up), matching
    build_ooh_structure's own unsampled default exactly.
    """
    n_polar = max(1, n_polar)
    n_azimuthal = max(1, n_azimuthal)
    polar_degs = [0.0] if n_polar == 1 else np.linspace(0.0, max_polar_deg, n_polar)
    az_degs = [0.0] if n_azimuthal == 1 else np.linspace(0.0, 360.0, n_azimuthal, endpoint=False)
    directions = []
    for polar_deg in polar_degs:
        polar = np.radians(polar_deg)
        for az_deg in az_degs:
            az = np.radians(az_deg)
            directions.append(np.array([
                np.sin(polar) * np.cos(az),
                np.sin(polar) * np.sin(az),
                np.cos(polar),
            ]))
    return directions


def build_ooh_structure_oriented(relaxed_oh, o_index, h_index, direction,
                                  oo_bond_ang=_OOH_OO_BOND_ANG, oh_bond_ang=_OOH_OH_BOND_ANG,
                                  bend_deg=_OOH_BEND_DEG):
    """Like build_ooh_structure, but O2 continues from the REAL O1 (
    relaxed_oh.atoms[o_index]'s own position -- the winning OH* site's own
    relaxed oxygen, which never moves) along `direction` (a unit vector,
    see _sample_ooh_directions) instead of continuing OH*'s own O-H bond
    direction; H is then bent off that SAME axis by `bend_deg`, exactly
    the Rodrigues construction build_ooh_structure itself uses.
    """
    o_symbol, o_frac = relaxed_oh.atoms[o_index]
    h_symbol, _h_frac = relaxed_oh.atoms[h_index]
    lattice = relaxed_oh.lattice
    inv_lattice = np.linalg.inv(lattice)
    o1_cart = o_frac @ lattice

    axis1 = direction / np.linalg.norm(direction)
    o2_cart = o1_cart + axis1 * oo_bond_ang
    rot_axis = _perpendicular_axis(axis1)
    bend_dir = _rodrigues_rotate(-axis1, rot_axis, np.radians(bend_deg))
    hnew_cart = o2_cart + bend_dir * oh_bond_ang

    o2_frac = o2_cart @ inv_lattice
    hnew_frac = hnew_cart @ inv_lattice

    species_meta = dict(relaxed_oh.species_meta)

    def _new_ads_label(element_symbol, z):
        label = f"{element_symbol}{ADSORBATE_LABEL_SUFFIX}"
        if label not in species_meta:
            used_ids = {str(info['id']) for info in species_meta.values()}
            next_id = 1
            while str(next_id) in used_ids:
                next_id += 1
            species_meta[label] = {'id': str(next_id), 'Z': z}
        return label

    o2_label = _new_ads_label("O", 8)
    hnew_label = _new_ads_label("H", 1)

    new_atoms = [a for i, a in enumerate(relaxed_oh.atoms) if i != h_index]
    new_atoms.append((o2_label, o2_frac))
    new_atoms.append((hnew_label, hnew_frac))

    species = list(dict.fromkeys(sym for sym, _ in new_atoms))
    return structure_io.FdfStructure(
        lattice=lattice, lattice_constant=relaxed_oh.lattice_constant,
        species=species, species_meta=species_meta,
        atoms=new_atoms, coord_format=relaxed_oh.coord_format, raw_lines=[],
    )


def sample_ooh_orientations(relaxed_oh, o_index, h_index, n_orientations_polar, n_orientations_azimuthal,
                             oo_bond_ang, oh_bond_ang, bend_deg, ml_rank, ml_model, ml_device, ml_fmax,
                             orientation_top_k, orientation_rmsd_tol, f_out):
    """Samples n_polar x n_azimuthal OOH* starting-orientation candidates,
    all anchored at the SAME O1 (the winning OH* site's own relaxed
    oxygen -- see build_ooh_structure_oriented's docstring; this is
    orientation sampling WITHIN one fixed site, never the old, now-removed
    per-intermediate SITE search). Without `ml_rank`, every sampled
    orientation is returned unranked (energy=None). With `ml_rank`, each is
    MACE-MP-0 relaxed (substrate fixed, O1/O2/H free -- same convention as
    ml_prerelax_adsorbate) and scored, then ranked, deduplicated
    (`orientation_rmsd_tol`, RMSD over the O1/O2/H adsorbate atoms only),
    and optionally trimmed to the `orientation_top_k` best.

    Returns a list of (fdf_structure, energy_or_None) for the kept
    candidates, best-first when ml_rank is True.
    """
    directions = _sample_ooh_directions(n_orientations_polar, n_orientations_azimuthal)

    if not ml_rank:
        print_dual(color_text(
            f"  [NOTE] Orientation sampling without --ml-prerelax: all {len(directions)} OOH* "
            "orientation(s) at the winning OH* site are written as their own relaxation folder "
            "below, unscreened.", 'yellow'), f_out)
        return [(build_ooh_structure_oriented(relaxed_oh, o_index, h_index, d,
                                               oo_bond_ang, oh_bond_ang, bend_deg), None)
                for d in directions]

    n_substrate = o_index  # every index before O1 is substrate (O1/H were OH*'s own adsorbate atoms)
    scored = []  # (energy, ase_atoms) -- same shape deduplicate_orientations expects
    structures = []  # fdf_structure per candidate, same order as `scored`, pre-relax (unused if kept)
    for d in directions:
        candidate = build_ooh_structure_oriented(relaxed_oh, o_index, h_index, d,
                                                  oo_bond_ang, oh_bond_ang, bend_deg)
        energy, relaxed_structure = _mace_relax_adsorbate_scored(
            candidate, n_substrate, ml_model, ml_device, ml_fmax)
        ase_atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(
            _bare_copy_for_ase(relaxed_structure)))
        scored.append((energy, ase_atoms))
        structures.append(relaxed_structure)

    order = sorted(range(len(scored)), key=lambda i: scored[i][0])
    scored_sorted = [scored[i] for i in order]
    structures_sorted = [structures[i] for i in order]
    e_min = scored_sorted[0][0]

    kept = deduplicate_orientations(scored_sorted, n_substrate, rmsd_tol=orientation_rmsd_tol)
    if orientation_top_k is not None:
        kept = kept[:orientation_top_k]

    for rank, i in enumerate(kept, start=1):
        energy = scored_sorted[i][0]
        print_dual(f"    orientation {rank}/{len(kept)}: E_MACE = {energy:.4f} eV, "
                    f"dE = {energy - e_min:+.4f} eV vs. best", f_out)
    print_dual(f"  {len(directions)} OOH* orientation(s) sampled at the winning OH* site -> "
                f"{len(kept)} unique kept"
                + (f" (--orientation-top-k {orientation_top_k})" if orientation_top_k is not None else ""),
                f_out)
    return [(structures_sorted[i], scored_sorted[i][0]) for i in kept]


def _bare_copy_for_ase(fdf_structure):
    """Bare-element copy of fdf_structure (fragment labels resolved to
    real elements) suitable for structure_io.to_pymatgen/AseAtomsAdaptor,
    used only to build an ase.Atoms for deduplicate_orientations' RMSD
    check (which needs .get_positions(), not fragment identity).
    """
    bare_species_meta = {}
    for orig_sym, _pos in fdf_structure.atoms:
        real_sym = structure_io.real_element(orig_sym, fdf_structure.species_meta)
        if real_sym not in bare_species_meta:
            bare_species_meta[real_sym] = {'id': str(len(bare_species_meta) + 1),
                                            'Z': abs(fdf_structure.species_meta[orig_sym]['Z'])}
    bare_atoms = [(structure_io.real_element(sym, fdf_structure.species_meta), pos)
                  for sym, pos in fdf_structure.atoms]
    return structure_io.FdfStructure(
        lattice=fdf_structure.lattice, lattice_constant=fdf_structure.lattice_constant,
        species=list(dict.fromkeys(sym for sym, _ in bare_atoms)), species_meta=bare_species_meta,
        atoms=bare_atoms, coord_format=fdf_structure.coord_format, raw_lines=[],
    )


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 4: derives the O* and OOH* intermediates from the "
        "winning OH* site, and writes their own CG-relaxation folder.", 'bold')}
Picks the lowest-FreeEng OH* site from Stage 1's 'sites/site_*/', reads its RELAXED geometry, and
builds starting geometries for O* and OOH* from it (O* = OH* minus H, always ONE folder -- a bare
O atom has no orientation to sample; OOH* = OH* plus a second O-H group in a chemically plausible
but illustrative starting orientation -- see build_ooh_structure --help/docstring), writing
'intermediates/o_star/' and 'intermediates/ooh_star/'. With --ooh-n-orientations-polar/-azimuthal
(> 1), OOH*'s new O-H group is instead sampled over several orientations AT THAT SAME SITE (O1,
the site's own relaxed oxygen, never moves -- see --help), optionally MACE-MP-0 ranked
(--ml-prerelax) down to the --orientation-top-k best, writing 'intermediates/ooh_star_orientN/'
folders instead of a single one.

[IMPORTANT] O*/OOH* are ALWAYS derived from the SAME winning OH* site, never independently
site-searched: the computational hydrogen electrode (CHE) descriptor's overpotential/PDS is only
physically meaningful as a single active site progressing through OH*->O*->OOH* (Rossmeisl et al.
2007; Man et al. 2011) -- the well-known ~3.2 eV universal scaling relation between
Delta-G(OOH*) and Delta-G(OH*) is itself derived assuming exactly this shared M-O bond across all
three intermediates. Deriving each intermediate from a DIFFERENT site (an earlier, now-removed
--o-strategy/--ooh-strategy search mode) mixes three different local bonding environments into
one nominal 'pathway', silently breaking that assumption and producing an eta/PDS that does not
correspond to any single physically realizable active site.

Every folder written here MUST be relaxed via SIESTA (MD.TypeOfRun CG is forced on) before
running stb-oerRefs -- O*/OOH*'s starting geometry is only a reasonable guess, not already at
equilibrium. Every folder gets its own config_extra.fdf sidecar (fixed cell,
Slab.DipoleCorrection, Spin polarized, DFTD3, all mandatory -- same combination Stage 1's own
site folders get, %include'd on top of your untouched --calc template, never edited in place).
Each derived geometry is also saved as a 1-frame extended-XYZ trajectory, viewable in OVITO/VMD.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory oer_study\n"
               "  %(prog)s --directory oer_study --ml-prerelax\n"
               "  %(prog)s --directory oer_study --ooh-n-orientations-polar 4 "
               "--ooh-n-orientations-azimuthal 4 --ml-prerelax --ml-device cuda "
               "--orientation-top-k 2\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="oer_study",
                         help="Root directory written by stb-oer (default: oer_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each site folder (default: calc.out).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--oo-bond-length", type=float, default=_OOH_OO_BOND_ANG,
                         help="Illustrative starting O-O bond length for OOH*, in Ang (default: "
                              f"{_OOH_OO_BOND_ANG}). Refined by CG relaxation.")
    parser.add_argument("--ooh-oh-bond-length", type=float, default=_OOH_OH_BOND_ANG,
                         help="Illustrative starting O-H bond length for OOH*'s new terminal H, "
                              f"in Ang (default: {_OOH_OH_BOND_ANG}).")
    parser.add_argument("--ooh-bend-deg", type=float, default=_OOH_BEND_DEG,
                         help="Illustrative starting O-O-H bend angle for OOH*, in degrees "
                              f"(default: {_OOH_BEND_DEG}, H2O2-like).")
    parser.add_argument("--ooh-n-orientations-polar", type=int, default=1,
                         help="Systematically sample this many initial OOH* orientations AT THE "
                              "WINNING OH* SITE before relaxing (default: 1, the single-"
                              "orientation behavior). O1 (the site's own relaxed oxygen) NEVER "
                              "moves -- only where O2/H point FROM it varies (same Fibonacci-"
                              "sphere sampling as stb-oer's own --n-orientations-polar, applied "
                              "to the local O1-O2-H template). No equivalent flag exists for O* "
                              "(a bare O atom has no orientation to sample). Combined with "
                              "--ooh-n-orientations-azimuthal below (total = polar x azimuthal). "
                              "Without --ml-prerelax, EVERY sampled orientation is written as its "
                              "own 'ooh_star_orientN/' folder directly. With --ml-prerelax, each "
                              "is MACE-MP-0 relaxed and ranked first, and only the unique/"
                              "--orientation-top-k survivors become relaxation folders.")
    parser.add_argument("--ooh-n-orientations-azimuthal", type=int, default=1,
                         help="Evenly spaced in-plane rotations sampled per polar direction (see "
                              "--ooh-n-orientations-polar). Default 1.")
    parser.add_argument("--orientation-top-k", type=int, default=None,
                         help="With --ml-prerelax and OOH* orientation sampling: keep only the N "
                              "best-ranked unique (post-deduplication) orientations, instead of "
                              "every surviving one. Unset (default) keeps every unique orientation.")
    parser.add_argument("--orientation-rmsd-tol", type=float, default=0.3,
                         help="With --ml-prerelax and OOH* orientation sampling: two relaxed "
                              "orientations within this RMSD (Ang, O1/O2/H only) AND within 0.01 "
                              "eV of each other are treated as duplicates (default: 0.3).")
    parser.add_argument("--ml-prerelax", action="store_true",
                         help="Relax O*/OOH*'s adsorbate atoms (positions only, substrate fixed) "
                              "with MACE-MP-0 before writing the CG-relaxation folder -- same "
                              "idea as stb-adsorb --ml-rank. Useful because the hand-built OOH* "
                              "starting geometry (O-O bond length/bend angle) is illustrative, "
                              "not literature-fitted. With OOH* orientation sampling (see "
                              "--ooh-n-orientations-polar/-azimuthal above), also drives the "
                              "MACE-MP-0 ranking/deduplication of sampled orientations. A "
                              "better-informed starting point for the real SIESTA relaxation, NOT "
                              "a substitute for it. Needs the optional 'ml' extra.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="medium",
                         help="MACE-MP-0 model size, with --ml-prerelax (default: medium).")
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

    if args.ooh_n_orientations_polar < 1 or args.ooh_n_orientations_azimuthal < 1:
        parser.error("--ooh-n-orientations-polar/--ooh-n-orientations-azimuthal must be >= 1.")
    if args.orientation_top_k is not None and not args.ml_prerelax:
        parser.error("--orientation-top-k is only valid with --ml-prerelax.")

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

        print_section('[0] RUN METADATA', f_out)
        print_dual(f"Directory       : {output_root}", f_out)
        print_dual(f"ML pre-relax    : {'yes' if args.ml_prerelax else 'no'}"
                    + (f" (model={args.ml_model}, device={args.ml_device})" if args.ml_prerelax else ""),
                    f_out)

        print_section('[1] WINNING OH* SITE', f_out)
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
            # The winning OH* site's own calc.fdf (written by oer.py's
            # write_site_folder) is itself '%include config_extra.fdf' +
            # the untouched user template -- strip that include before
            # using this text as the base for THIS stage's own derived
            # folders, each of which gets its OWN config_extra.fdf (see
            # write_relax_folder). A no-op if the site folder predates
            # this convention (plain calc_text already).
            site_calc_text = strip_config_extra_include(f.read())

        print_section('[2] O* GEOMETRY', f_out)
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
        o_traj_path = os.path.join(output_root, "intermediates", "o_trajectory.xyz")
        o_frame = _bare_ase_atoms(o_structure)
        o_frame.info["site_label"] = "o_star"
        ase_io.write(o_traj_path, [o_frame], format="extxyz")
        print_dual(f"  {color_text('[Saved]', 'cyan')} {o_traj_path} (1 frame, OVITO/VMD-"
                    "viewable)", f_out)

        print_section('[3] OOH* GEOMETRY', f_out)
        ooh_orientation_sampling = (args.ooh_n_orientations_polar > 1
                                     or args.ooh_n_orientations_azimuthal > 1)
        if not ooh_orientation_sampling:
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
            ooh_traj_path = os.path.join(output_root, "intermediates", "ooh_trajectory.xyz")
            ooh_frame = _bare_ase_atoms(ooh_structure)
            ooh_frame.info["site_label"] = "ooh_star"
            ase_io.write(ooh_traj_path, [ooh_frame], format="extxyz")
            print_dual(f"  {color_text('[Saved]', 'cyan')} {ooh_traj_path} (1 frame, OVITO/VMD-"
                        "viewable)", f_out)
        else:
            if args.ml_prerelax:
                print_dual(f"  {color_text('ML pre-relax:', 'cyan')} sampling "
                            f"{args.ooh_n_orientations_polar}x{args.ooh_n_orientations_azimuthal} "
                            "OOH* orientations AT THE WINNING OH* SITE (O1 fixed), relaxing/"
                            "ranking each with MACE-MP-0 (substrate fixed) ...", f_out)
            kept = sample_ooh_orientations(
                relaxed_oh, o_index, h_index, args.ooh_n_orientations_polar,
                args.ooh_n_orientations_azimuthal, args.oo_bond_length, args.ooh_oh_bond_length,
                args.ooh_bend_deg, args.ml_prerelax, args.ml_model, args.ml_device, args.ml_fmax,
                args.orientation_top_k, args.orientation_rmsd_tol, f_out)
            trajectory_frames = []
            for rank, (ooh_structure, energy) in enumerate(kept, start=1):
                ooh_dir = os.path.join(output_root, "intermediates", f"ooh_star_orient{rank}")
                ooh_calc = force_system_label(force_relaxation(site_calc_text),
                                               f"oer_ooh_star_orient{rank}")
                write_relax_folder(ooh_dir, ooh_structure, ooh_calc, args.pseudo_dir)
                energy_note = f", E_MACE = {energy:.4f} eV" if energy is not None else ""
                print_dual(f"  {color_text('[OK]', 'green')} {ooh_dir} (orientation {rank}/"
                            f"{len(kept)}, anchored at the winning OH* site's own O1{energy_note})",
                            f_out)
                frame = _bare_ase_atoms(ooh_structure)
                frame.info["site_label"] = f"ooh_star_orient{rank}"
                if energy is not None:
                    frame.info["energy_eV"] = round(energy, 6)
                trajectory_frames.append(frame)
            ooh_traj_path = os.path.join(output_root, "intermediates", "ooh_trajectory.xyz")
            ase_io.write(ooh_traj_path, trajectory_frames, format="extxyz")
            print_dual(f"  {color_text('[Saved]', 'cyan')} {ooh_traj_path} "
                        f"({len(trajectory_frames)} frame(s), OVITO/VMD-viewable)", f_out)

        print_section('[4] SUMMARY & NEXT STEPS', f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual("  1. Run SIESTA in every folder written above.", f_out)
        print_dual(f"  2. Once they're done, run: stb-oerRefs --directory {output_root}", f_out)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Intermediate geometries ready for Stage 3 (stb-oerRefs).\n", 'bold'))


if __name__ == "__main__":
    main()
