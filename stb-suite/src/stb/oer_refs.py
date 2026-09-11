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
import json
import argparse
import numpy as np
from pymatgen.core.periodic_table import Element
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.calc_directives import force_single_point
from stb.core.siesta_log import get_free_energy, get_outcell, check_scf_and_force, report_quality_diagnostics
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders
from stb.core.adsorption_sites import (
    CONFIG_EXTRA_FILE, FIXED_CELL_BLOCK, DIPOLE_CORRECTION_BLOCK, SPIN_POLARIZED_BLOCK,
    VDW_CORRECTION_BLOCK, SINGLE_POINT_BLOCK,
)
from stb.core.bsse import strip_config_extra_include

REPORT_FILE = "oer_stage3.txt"
_DEFAULT_LOCAL_DISPLACEMENT_ANG = 0.015
_DEFAULT_VACUUM_BOX_ANG = 15.0
_H2_BOND_LENGTH_ANG = 0.741  # experimental equilibrium bond length -- CG relaxation refines it

# config_extra.fdf content for every SINGLE-POINT folder this module writes
# (00_clean_slab, 04_slab_deformed, the BSSE ghost quartets, every local/full
# ZPE displacement folder): fixed cell + the winning OH* site's own mandatory
# Slab.DipoleCorrection/Spin polarized/DFTD3 (oer.py's write_site_folder,
# where all three are mandatory) -- propagating the winning site's own
# numerical settings into every derived single-point folder is required for
# a physically meaningful energy difference (the "level-of-theory
# propagation" convention documented in CLAUDE.md). Exact mirror of
# her_refs.py's own _SINGLE_POINT_CONFIG_EXTRA.
_SINGLE_POINT_CONFIG_EXTRA = (FIXED_CELL_BLOCK + DIPOLE_CORRECTION_BLOCK + SPIN_POLARIZED_BLOCK
                               + VDW_CORRECTION_BLOCK + SINGLE_POINT_BLOCK)
# For 02_h2_molecule/03_h2o_molecule only: these RELAX (not single-point) and
# get 'Spin non-polarized' forced via the local force_spin_nonpolarized()
# TEXT helper below (not through config_extra.fdf, and NOT
# SPIN_POLARIZED_BLOCK -- both are closed-shell singlets, see that
# function's own docstring) -- so this omits SPIN_POLARIZED_BLOCK/
# SINGLE_POINT_BLOCK, keeping only the fixed-cell/dipole/vdW directives
# every derived folder in this module needs.
_RELAX_CONFIG_EXTRA = FIXED_CELL_BLOCK + DIPOLE_CORRECTION_BLOCK + VDW_CORRECTION_BLOCK

_KGRID_RE = re.compile(r'kgrid[._]MonkhorstPack\s+\[.*?\]', re.IGNORECASE)
_SPIN_RE = re.compile(r'(Spin\s+)(\S+)', re.IGNORECASE)
_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)
_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_MD_STEPS_RE = re.compile(r'MD\.Steps\s+\d+', re.IGNORECASE)
_RELAXED_COORDS_RE = re.compile(r'outcoor:\s*Relaxed atomic coordinates\s*\(fractional\)', re.IGNORECASE)


# --- duplicated from her_refs.py / oer_intermediates.py (self-contained stage) -------------

def force_gamma_kgrid(calc_text):
    new_text, count = _KGRID_RE.subn('kgrid.MonkhorstPack   [1  1  1]', calc_text)
    if count == 0:
        raise ValueError("Could not find a 'kgrid.MonkhorstPack' tag in the calc.fdf template.")
    return new_text


def force_spin_nonpolarized(calc_text):
    """Forces 'Spin non-polarized' -- used for BOTH H2 and H2O (the only
    two callers of this function in this module), never 'Spin polarized':
    both are unambiguous closed-shell singlets (H2: 2 electrons filling
    the bonding sigma orbital; H2O: 8 valence electrons, all paired, no
    radical character at all), unlike the winning OH*/O*/OOH* site (which
    CAN genuinely be open-shell, hence 'Spin polarized' there via
    config_extra.fdf, inherited through site_calc_text). Forcing
    'Spin polarized' on either gas-phase reference risks the SCF
    converging to a spurious nonzero magnetic moment depending on the
    initial guess, silently biasing 0.5*G(H2)/G(H2O) in every Delta-G
    computed against them -- same reasoning as HER's own her_refs.py
    (_SPIN_NONPOLARIZED_BLOCK), which this was a real, live divergence
    from until fixed here.
    """
    new_text, count = _SPIN_RE.subn(r'\g<1>non-polarized', calc_text)
    if count == 0:
        return calc_text + "\nSpin                non-polarized\n"
    return new_text


def force_relaxation(calc_text, max_steps=200):
    new_text, count = _RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _MD_STEPS_RE.subn(f'MD.Steps              {max_steps}', new_text)
    if count == 0:
        new_text += f"MD.Steps              {max_steps}\n"
    return new_text


def force_system_label(calc_text, label):
    new_text, count = _LABEL_RE.subn(f'SystemLabel {label}', calc_text)
    if count == 0:
        new_text += f"\nSystemLabel {label}\n"
    return new_text


def find_winning_site(sites_root, out_file, f_out, prefix="site_"):
    """Scans every folder under sites_root whose name starts with `prefix`
    for the lowest FreeEng -- used to find the winning OH* site under
    sites/ (default prefix), and reused by _locate_intermediate below
    (prefix='ooh_star_orient') to pick the winner among OOH* orientation
    candidates sampled AT THE SAME SITE (never a different one -- see
    oer_intermediates.py's module docstring).
    """
    site_dirs = sorted(
        d.path for d in os.scandir(sites_root) if d.is_dir() and d.name.startswith(prefix)
    )
    if not site_dirs:
        print_dual(color_text(f"[ERROR] No '{prefix}*' folders found in '{sites_root}'.", 'red'), f_out)
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
        print_dual(color_text(f"[ERROR] No folder under '{sites_root}' has a readable energy -- "
                               "did SIESTA finish?", 'red'), f_out)
        sys.exit(1)
    return best_dir, best_energy, results


def read_relaxed_structure(out_path, template):
    """Reads the LAST 'outcoor: Relaxed atomic coordinates (fractional)'
    block from a SIESTA .out file -- SIESTA never rewrites structure.fdf
    in place. Returns None if the block isn't found.
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


def build_h2_structure(vacuum_box, bond_length_ang=_H2_BOND_LENGTH_ANG):
    """H2 molecule centered in its own cubic vacuum-box periodic cell --
    the CHE gas-phase reference for every (H+ + e-) pair in the 4 OER
    steps, exactly as in HER. Duplicated verbatim from her_refs.py.
    """
    half = bond_length_ang / 2.0
    center = vacuum_box / 2.0
    cart = np.array([
        [center - half, center, center],
        [center + half, center, center],
    ])
    frac = cart / vacuum_box
    lattice = np.eye(3) * vacuum_box
    species_meta = {"H": {"id": "1", "Z": 1}}
    atoms = [("H", frac[0]), ("H", frac[1])]
    return structure_io.FdfStructure(
        lattice=lattice, lattice_constant=1.0, species=["H"], species_meta=species_meta,
        atoms=atoms, coord_format="fractional", raw_lines=[],
    )


def build_h2o_structure(vacuum_box):
    """Isolated H2O molecule centered in its own cubic vacuum-box periodic
    cell, geometry from ASE's G2 'H2O' entry (Angstrom, gas-phase
    equilibrium) -- unlike H2's single bond-length scalar, H2O's
    angle+bond-length combination is more naturally sourced from ASE's
    validated database than hand-rolled, same external-geometry-source
    technique as adsorb.py::isolated_adsorbate_structure (duplicated
    here, not imported -- self-contained workflow). CG-relaxed further by
    the caller (force_relaxation), same treatment as H2.
    """
    from ase.build import molecule as ase_build_molecule
    atoms = ase_build_molecule("H2O")
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    center = np.array([vacuum_box / 2.0] * 3)
    cart = positions - positions.mean(axis=0) + center
    lattice = np.eye(3) * vacuum_box
    frac = cart @ np.linalg.inv(lattice)

    species_meta = {}
    for sym in symbols:
        if sym not in species_meta:
            species_meta[sym] = {'id': str(len(species_meta) + 1), 'Z': Element(sym).Z}
    species = list(dict.fromkeys(symbols))
    atoms_list = [(sym, frac[i]) for i, sym in enumerate(symbols)]
    return structure_io.FdfStructure(
        lattice=lattice, lattice_constant=1.0, species=species, species_meta=species_meta,
        atoms=atoms_list, coord_format="fractional", raw_lines=[],
    )


def make_ghost_variant(base_structure, ghost_start, ghost_end):
    """Returns a copy of base_structure with atoms in [ghost_start,
    ghost_end) turned into ghost species ('<symbol>_ghost', negative Z,
    same real pseudopotential file via copy_pseudo's dest_label) --
    same Boys-Bernardi counterpoise convention as her_refs.py's own
    make_ghost_variant, duplicated here. Adsorbate atoms are always
    appended last (guaranteed by oer.py/oer_intermediates.py), so a
    contiguous [start, end) range always selects exactly "the slab" or
    "the adsorbate", regardless of how many adsorbate atoms there are.
    `symbol` here is already a Stage-1/2 fragment label
    ('<real>_slab'/'<real>_ads', see oer.py's write_site_folder) rather
    than a bare element symbol, so the real Z is read straight out of
    `species_meta` (already declared for every label present) instead of
    constructing a pymatgen Element from the label text -- Element(symbol)
    would raise on a non-bare label like 'O_ads'.
    """
    species_meta = dict(base_structure.species_meta)
    new_atoms = []
    for i, (symbol, pos) in enumerate(base_structure.atoms):
        if ghost_start <= i < ghost_end:
            label = f"{symbol}_ghost"
            if label not in species_meta:
                real_z = species_meta[symbol]['Z']
                used_ids = {str(info['id']) for info in species_meta.values()}
                next_id = 1
                while str(next_id) in used_ids:
                    next_id += 1
                species_meta[label] = {'id': str(next_id), 'Z': -abs(real_z)}
        else:
            label = symbol
        new_atoms.append((label, pos))

    species = list(dict.fromkeys(sym for sym, _ in new_atoms))
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=species, species_meta=species_meta, atoms=new_atoms,
        coord_format=base_structure.coord_format, raw_lines=[],
    )


def remove_atoms(base_structure, indices):
    """Generalizes a single-index atom removal (her_refs.py::remove_atom)
    to a list -- '04_slab_deformed' needs to remove BOTH of OH*'s atoms
    (O and H), not just one.
    """
    indices = set(indices)
    new_atoms = [a for i, a in enumerate(base_structure.atoms) if i not in indices]
    present = {sym for sym, _ in new_atoms}
    species_meta = {k: v for k, v in base_structure.species_meta.items() if k in present}
    species = [s for s in base_structure.species if s in present]
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=species, species_meta=species_meta, atoms=new_atoms,
        coord_format=base_structure.coord_format, raw_lines=[],
    )


def isolate_atoms(base_structure, indices):
    """Generalizes a single-index atom isolation (her_refs.py::isolate_atom)
    to a list -- OOH*'s isolated-adsorbate BSSE reference needs
    {O1,O2,Hnew} (3 atoms), not just 1. Reverts each isolated atom's
    Stage-1/2 fragment label ('<real>_ads') back to the bare real element
    symbol via structure_io.real_element: with no slab atoms left in this
    single-fragment folder, there is nothing left to disambiguate from
    (same convention as her_refs.py's own isolate_atom).
    """
    indices = sorted(indices)
    kept = [base_structure.atoms[i] for i in indices]
    species_meta = {}
    new_kept = []
    for symbol, pos in kept:
        real_symbol = structure_io.real_element(symbol, base_structure.species_meta)
        if real_symbol not in species_meta:
            real_z = abs(base_structure.species_meta[symbol]['Z'])
            species_meta[real_symbol] = {'id': str(len(species_meta) + 1), 'Z': real_z}
        new_kept.append((real_symbol, pos))
    species = list(dict.fromkeys(sym for sym, _ in new_kept))
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=species, species_meta=species_meta, atoms=new_kept,
        coord_format=base_structure.coord_format, raw_lines=[],
    )


def formula_summary(fdf_structure):
    """Returns a compact, human-readable formula string for a written
    folder's report row, e.g. 'B9N9' or 'B9N9 +H(ghost)' -- real atoms
    grouped by element (via structure_io.real_element, so a Stage-1/2
    fragment label like 'B_slab' or a ghost label like 'H_ads_ghost' both
    collapse to their real element), ghost atoms (negative Z) called out
    separately since they contribute zero electrons/charge despite
    sharing the real pseudopotential. Exact mirror of her_refs.py's own
    formula_summary.
    """
    from collections import Counter
    real_counts, ghost_counts = Counter(), Counter()
    for label, _ in fdf_structure.atoms:
        element = structure_io.real_element(label, fdf_structure.species_meta)
        if fdf_structure.species_meta[label]['Z'] < 0:
            ghost_counts[element] += 1
        else:
            real_counts[element] += 1
    formula = "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(real_counts.items())) or "-"
    if ghost_counts:
        formula += " +" + "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(ghost_counts.items()))
        formula += "(ghost)"
    return formula


def write_folder(out_dir, fdf_structure, calc_text, pp_path, config_extra_content):
    """Writes structure.fdf + calc.fdf + config_extra.fdf + copied pseudos
    for one derived reference folder, following the same config_extra.fdf
    sidecar convention as oer.py/oer_intermediates.py (4.8/4.11/4.12/4.13
    model): `config_extra_content` (this folder's own combination of
    forced directives -- single-point vs. relaxation, dipole correction,
    fixed cell) is written as-is to config_extra.fdf, and
    `%include config_extra.fdf` is prepended to the UNTOUCHED `calc_text`
    (structure_io.prepend_include) rather than editing directives into it
    in place. Handles both a Stage-1/2 fragment label
    ('<real>_slab'/'<real>_ads') and a ghost label stacked on top of one
    ('<real>_slab_ghost', from make_ghost_variant) transparently for
    pseudopotential copying: the real element behind ANY label is
    recovered via structure_io.real_element (Z-based, robust to any
    suffix or stack of suffixes -- naive string-slicing off '_ghost'
    alone, this function's previous approach, silently mis-resolved a
    plain fragment label like 'O_ads' to a nonexistent 'O_ads.psf' source
    pseudopotential). Duplicated from her_refs.py's own write_folder.
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
        f.write(config_extra_content)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(structure_io.prepend_include(calc_text, CONFIG_EXTRA_FILE))
    present_labels = sorted({symbol for symbol, _ in fdf_structure.atoms})
    for label in present_labels:
        real_symbol = structure_io.real_element(label, fdf_structure.species_meta)
        copy_pseudo(pp_path, real_symbol, out_dir, dest_label=label)


def write_bsse_triad(out_root, prefix, base_structure, n_substrate, calc_text, pp_path,
                      config_extra_content):
    """Writes a Boys-Bernardi counterpoise QUARTET at base_structure's
    geometry: <prefix>_slab_only/ (slab atoms only, adsorbate REMOVED --
    not ghosted), <prefix>_slab_ghost/ (real slab + ghost adsorbate),
    <prefix>_adsorbate_ghost_slab/ (ghost slab + real adsorbate),
    <prefix>_isolated/ (adsorbate alone). BSSE_slab = E(slab_only) -
    E(slab_ghost) and BSSE_adsorbate = E(isolated) - E(adsorbate_ghost_slab)
    -- both differences MUST compare the SAME geometry across basis sets,
    which is why slab_only is written HERE (at base_structure's own
    geometry) rather than reusing '04_slab_deformed' (always built at
    OH*'s geometry, which generally differs from O*'s/OOH*'s -- reusing
    the diagnostic-only deformed-slab folder there would silently compare
    two different geometries and corrupt the counterpoise correction).
    Called once per intermediate (OH*, O*, then OOH*), each at ITS OWN
    relaxed geometry -- a single triad shared across all 3 (built once,
    reused for the other two) was removed after it was found to
    systematically under-correct OH*/O* relative to their true per-species
    BSSE (see main()'s own [5] BSSE CORRECTION comment for the real numbers).

    Returns a list of (label, fdf_structure, run_type) for the caller's
    own [4]/[5] recap table (see formula_summary/print_table in main()).
    """
    n_total = len(base_structure.atoms)
    written = []

    slab_only = remove_atoms(base_structure, list(range(n_substrate, n_total)))
    write_folder(os.path.join(out_root, f"{prefix}_slab_only"), slab_only,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_slab_only"), pp_path,
                 config_extra_content)
    written.append((f"{prefix}_slab_only", slab_only, "single-point (BSSE)"))

    slab_ghost = make_ghost_variant(base_structure, n_substrate, n_total)
    write_folder(os.path.join(out_root, f"{prefix}_slab_ghost"), slab_ghost,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_slab_ghost"), pp_path,
                 config_extra_content)
    written.append((f"{prefix}_slab_ghost", slab_ghost, "single-point (BSSE)"))

    ads_ghost_slab = make_ghost_variant(base_structure, 0, n_substrate)
    write_folder(os.path.join(out_root, f"{prefix}_adsorbate_ghost_slab"), ads_ghost_slab,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_adsorbate_ghost_slab"),
                 pp_path, config_extra_content)
    written.append((f"{prefix}_adsorbate_ghost_slab", ads_ghost_slab, "single-point (BSSE)"))

    isolated = isolate_atoms(base_structure, list(range(n_substrate, n_total)))
    write_folder(os.path.join(out_root, f"{prefix}_isolated"), isolated,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_isolated"), pp_path,
                 config_extra_content)
    written.append((f"{prefix}_isolated", isolated, "single-point (BSSE)"))

    return written


def write_local_zpe_folders(zpe_dir, relaxed_structure, local_indices, displacement_ang, calc_text,
                             pp_path, system_label, config_extra_content):
    """Generalizes her_refs.py::write_local_zpe_folders from ONE local
    atom (3 DOF, 6 folders) to a LIST of local atom indices (n atoms, 3n
    DOF, 6n folders): for every local atom and axis, writes a +/-
    displacement folder (every other atom, including OTHER local atoms,
    held fixed at its relaxed position). For a single local atom this
    produces the exact same 6-folder disp_001..disp_006 shape HER's own
    H*-only case already does (axis-then-sign order matches HER's fixed
    _LOCAL_DISPLACEMENTS list). Writes a JSON sidecar (zpe_local_meta.json)
    with the local atom indices/symbols (per-atom mass lookup in Stage 4),
    displacement_ang, the exact (atom_index, axis, sign) order, and the
    SystemLabel used (so Stage 4 knows which '<label>.FA' filename to
    read in each disp_NNN/ folder without guessing).
    """
    os.makedirs(zpe_dir, exist_ok=True)
    inv_lattice = np.linalg.inv(relaxed_structure.lattice)
    # Resolved to the BARE real element (not relaxed_structure's own
    # Stage-1/2 fragment label, e.g. 'O_ads') before being written to
    # zpe_local_meta.json below -- stb-oerAnalysis looks up atomic mass by
    # this symbol directly via pymatgen's Element(), which cannot parse a
    # fragment-suffixed label.
    local_symbols = [structure_io.real_element(relaxed_structure.atoms[i][0], relaxed_structure.species_meta)
                      for i in local_indices]

    order = []
    for atom_index in local_indices:
        for axis in range(3):
            for sign in (1.0, -1.0):
                order.append({"atom_index": atom_index, "axis": axis, "sign": sign})

    for i, entry in enumerate(order, start=1):
        atom_index, axis, sign = entry["atom_index"], entry["axis"], entry["sign"]
        symbol, frac = relaxed_structure.atoms[atom_index]
        cart = frac @ relaxed_structure.lattice
        delta_cart = np.zeros(3)
        delta_cart[axis] = sign * displacement_ang
        new_frac = (cart + delta_cart) @ inv_lattice
        new_atoms = list(relaxed_structure.atoms)
        new_atoms[atom_index] = (symbol, new_frac)
        disp_structure = structure_io.FdfStructure(
            lattice=relaxed_structure.lattice, lattice_constant=relaxed_structure.lattice_constant,
            species=relaxed_structure.species, species_meta=relaxed_structure.species_meta,
            atoms=new_atoms, coord_format=relaxed_structure.coord_format, raw_lines=[],
        )
        disp_dir = os.path.join(zpe_dir, f"disp_{i:03d}")
        write_folder(disp_dir, disp_structure, calc_text, pp_path, config_extra_content)

    with open(os.path.join(zpe_dir, "zpe_local_meta.json"), "w") as f:
        json.dump({
            "local_indices": list(local_indices),
            "local_symbols": local_symbols,
            "displacement_ang": displacement_ang,
            "order": order,
            "system_label": system_label,
        }, f)


def _locate_intermediate(output_root, name, out_file, f_out):
    """Locates the final relaxed geometry+calc.fdf for O*/OOH* ('name' is
    'o' or 'ooh'): the single 'intermediates/<name>_star/' folder if it
    exists, else (OOH* orientation sampling only -- O* never has this,
    a bare O atom has no orientation to sample) the winner among
    'intermediates/<name>_star_orient*/' -- MULTIPLE ORIENTATIONS AT THE
    SAME SITE (O1 never moves between them, see oer_intermediates.py's
    sample_ooh_orientations), not the old, now-removed per-intermediate
    SITE search. Either way, stb-oerIntermediates (Stage 2) derives O*/OOH*
    from the SAME winning OH* site unconditionally. Returns
    (relaxed_structure, source_dir, calc_text).
    """
    label_upper = f"{name.upper()}*"
    intermediates_root = os.path.join(output_root, "intermediates")
    single_dir = os.path.join(intermediates_root, f"{name}_star")
    if os.path.isdir(single_dir):
        source_dir = single_dir
    else:
        orient_prefix = f"{name}_star_orient"
        if not os.path.isdir(intermediates_root) or not any(
                d.is_dir() and d.name.startswith(orient_prefix)
                for d in os.scandir(intermediates_root)):
            print_dual(color_text(f"[ERROR] Neither '{single_dir}' nor '{orient_prefix}*' found "
                                   "-- run stb-oerIntermediates (Stage 2) first.", 'red'), f_out)
            sys.exit(1)
        source_dir, energy, all_results = find_winning_site(
            intermediates_root, out_file, f_out, prefix=orient_prefix)
        n_readable = sum(1 for _l, e in all_results if e is not None)
        print_dual(f"  {len(all_results)} {label_upper} orientation(s) scanned at the winning OH* "
                    f"site, {n_readable} with a readable FreeEng.", f_out)
        for label, e in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(intermediates_root, label) == source_dir else ""
            e_str = f"{e:.6f} eV" if e is not None else "(no energy)"
            print_dual(f"  {label:<28}{e_str}{marker}", f_out)
        readable = [e for _l, e in all_results if e is not None]
        if len(readable) > 1:
            spread = max(readable) - min(readable)
            print_dual(f"  Energy spread across readable orientations: {spread:.4f} eV (max - min).",
                        f_out)

    template = structure_io.read_fdf(os.path.join(source_dir, "structure.fdf"))
    relaxed = read_relaxed_structure(os.path.join(source_dir, out_file), template)
    if relaxed is None:
        print_dual(color_text(
            f"[ERROR] Could not read relaxed coordinates from '{source_dir}/{out_file}' -- did "
            "the relaxation finish?", 'red'), f_out)
        sys.exit(1)
    winning_energy = get_free_energy(os.path.join(source_dir, out_file))
    energy_str = f"{winning_energy:.6f} eV" if winning_energy is not None else "(no energy)"
    print_dual(f"{label_upper} site : {os.path.relpath(source_dir, output_root)} "
                f"({energy_str}, derived from the winning OH* site)", f_out)
    report_quality_diagnostics(os.path.basename(source_dir), os.path.join(source_dir, out_file),
                                0.05, f_out)
    with open(os.path.join(source_dir, "calc.fdf")) as f:
        # source_dir's own calc.fdf (written by oer.py/oer_intermediates.py)
        # is itself '%include config_extra.fdf' + the untouched user
        # template -- strip that include before reusing this text as the
        # base for THIS stage's own derived folders (see write_folder).
        calc_text = strip_config_extra_include(f.read())
    return relaxed, source_dir, calc_text


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 3 of 4: builds the H2/H2O gas-phase references, the BSSE "
        "counterpoise correction, and the ZPE/entropy calculation folders for all 3 intermediates "
        "(OH*, O*, OOH*).", 'bold')}
Locates the final relaxed geometry of OH* (Stage 1's winning site), O*, and OOH* (Stage 2's
derived-and-relaxed intermediates, both always from that SAME winning OH* site), then writes (all
single-point unless noted): '00_clean_slab/' (pristine slab), '02_h2_molecule/' and
'03_h2o_molecule/' (RELAX, gas-phase CHE references), '04_slab_deformed/' (OH*'s geometry with
both adsorbate atoms removed -- diagnostic only), the BSSE counterpoise triads (one PER
INTERMEDIATE -- 05_bsse_OH_*/05_bsse_O_*/05_bsse_OOH_*, 9 folders total, each triad built at
THAT intermediate's own relaxed geometry -- a single triad shared across all 3 intermediates
was removed after a real system showed it systematically under-corrects OH*/O* relative to
their true per-species BSSE, changing eta by ~0.30 V), and the ZPE calculation folder(s)
(--zpe-mode local/full -- no 'standard' mode: OER's per-species thermal corrections are always
computed from your own DFT/phonon data, never a hardcoded literature default, see
stb-oerAnalysis --help). Doesn't run SIESTA -- run each folder yourself, then use
stb-oerAnalysis.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory oer_study --zpe-mode local\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="oer_study",
                         help="Root directory written by stb-oer/stb-oerIntermediates "
                              "(default: oer_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (default: reuse whatever the winning "
                              "OH* site's own folder already has).")
    parser.add_argument("--zpe-mode", choices=["local", "full"], default="local",
                         help="ZPE/entropy calculation mode (default: local). 'local': "
                              "finite-difference displacements of each intermediate's own "
                              "adsorbate atoms only (fast, decoupled-oscillator approximation). "
                              "'full': complete Phonopy displacement sets for the clean slab AND "
                              "all 3 intermediates (needed to subtract the clean slab's own "
                              "phonon ZPE/entropy -- see stb-oerAnalysis --help).")
    parser.add_argument("--displacement", type=float, default=_DEFAULT_LOCAL_DISPLACEMENT_ANG,
                         help="Finite-difference displacement in Ang, for --zpe-mode local's "
                              "adsorbate displacements or --zpe-mode full's Phonopy displacements "
                              f"(default: {_DEFAULT_LOCAL_DISPLACEMENT_ANG}).")
    parser.add_argument("--supercell", type=int, nargs=3, default=[1, 1, 1],
                         help="Supercell dimensions for --zpe-mode full's phonon calculations "
                              "(default: 1 1 1).")
    parser.add_argument("--vacuum-box", type=float, default=_DEFAULT_VACUUM_BOX_ANG,
                         help=f"Cubic box side (Ang) for the isolated H2/H2O references (default: "
                              f"{_DEFAULT_VACUUM_BOX_ANG}).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-oerRefs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OER WORKFLOW -- STAGE 3: REFERENCES, BSSE & ZPE PREP", 'bold'))
    print("-" * 60)

    output_root = args.directory
    sites_root = os.path.join(output_root, "sites")
    clean_slab_source = os.path.join(output_root, "clean_slab_source", "structure.fdf")
    stage2_report = os.path.join(output_root, "oer_stage2.txt")
    if not os.path.isfile(stage2_report):
        print(color_text(f"[ERROR] '{stage2_report}' not found -- run stb-oerIntermediates "
                          "(Stage 2) first.", 'red'))
        sys.exit(1)

    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== OER STAGE 3 REPORT (REFERENCES, BSSE & ZPE PREP) =====', 'magenta')}", f_out)

        print_section('[0] RUN METADATA', f_out)
        print_dual(f"Directory       : {output_root}", f_out)
        print_dual(f"SIESTA output   : {args.file} (scanned in every site_*/candidates/ and read "
                    "back below)", f_out)
        print_dual(f"ZPE mode        : {args.zpe_mode}", f_out)
        print_dual(f"Displacement    : {args.displacement} Ang", f_out)
        print_dual(f"Vacuum box      : {args.vacuum_box:.1f} Ang (H2/H2O gas-phase references)", f_out)
        if args.zpe_mode == "full":
            print_dual(f"Supercell       : {args.supercell[0]} {args.supercell[1]} "
                        f"{args.supercell[2]} (--zpe-mode full's phonon calculations)", f_out)
        print_dual(f"Pseudo dir      : {args.pseudo_dir or '(reuse the winning OH* site itself has)'}",
                    f_out)
        print_dual(f"Report          : {report_path}", f_out)

        print_section('[1] WINNING OH* SITE', f_out)
        winning_oh_dir, winning_oh_energy, all_results = find_winning_site(sites_root, args.file, f_out)
        n_readable = sum(1 for _l, e in all_results if e is not None)
        print_dual(f"  {len(all_results)} site(s) scanned, {n_readable} with a readable FreeEng.",
                    f_out)
        for label, energy in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(sites_root, label) == winning_oh_dir else ""
            energy_str = f"{energy:.6f} eV" if energy is not None else "(no energy)"
            print_dual(f"  {label:<28}{energy_str}{marker}", f_out)
        readable = [e for _l, e in all_results if e is not None]
        if len(readable) > 1:
            spread = max(readable) - min(readable)
            print_dual(f"  Energy spread across readable sites: {spread:.4f} eV (max - min).", f_out)
        print_dual(f"Winning OH* site : {os.path.basename(winning_oh_dir)} ({winning_oh_energy:.6f} eV)", f_out)
        report_quality_diagnostics(os.path.basename(winning_oh_dir),
                                    os.path.join(winning_oh_dir, args.file), 0.05, f_out)

        oh_template = structure_io.read_fdf(os.path.join(winning_oh_dir, "structure.fdf"))
        relaxed_oh = read_relaxed_structure(os.path.join(winning_oh_dir, args.file), oh_template)
        if relaxed_oh is None:
            print_dual(color_text(
                f"[ERROR] Could not read relaxed coordinates from '{winning_oh_dir}/{args.file}' "
                "-- did the relaxation finish?", 'red'), f_out)
            sys.exit(1)
        n_oh_total = len(relaxed_oh.atoms)
        n_substrate = n_oh_total - 2
        oh_o_index, oh_h_index = n_oh_total - 2, n_oh_total - 1
        with open(winning_oh_dir + "/calc.fdf") as f:
            # Same stale-include stripping as _locate_intermediate above --
            # the winning OH* site's own calc.fdf (oer.py's
            # write_site_folder) is '%include config_extra.fdf' + the
            # untouched user template.
            site_calc_text = strip_config_extra_include(f.read())

        print_section('[2] O* FINAL GEOMETRY', f_out)
        relaxed_o, o_source_dir, o_calc_text = _locate_intermediate(
            output_root, "o", args.file, f_out)

        print_section('[3] OOH* FINAL GEOMETRY', f_out)
        relaxed_ooh, ooh_source_dir, ooh_calc_text = _locate_intermediate(
            output_root, "ooh", args.file, f_out)
        n_ooh_total = len(relaxed_ooh.atoms)

        print_section('[4] REFERENCE FOLDERS', f_out)
        ref_rows = []  # (label, fdf_structure, run_type) -- table at the end

        clean_template = structure_io.read_fdf(clean_slab_source)
        clean_dir = os.path.join(output_root, "00_clean_slab")
        clean_calc = force_system_label(force_single_point(site_calc_text), "oer_clean_slab")
        write_folder(clean_dir, clean_template, clean_calc, args.pseudo_dir, _SINGLE_POINT_CONFIG_EXTRA)
        ref_rows.append(("00_clean_slab", clean_template, "single-point"))
        print_dual(f"  {color_text('[OK]', 'green')} {clean_dir}", f_out)

        h2_dir = os.path.join(output_root, "02_h2_molecule")
        h2_structure = build_h2_structure(args.vacuum_box)
        h2_calc = force_system_label(
            force_spin_nonpolarized(force_gamma_kgrid(force_relaxation(site_calc_text))),
            "oer_h2_molecule")
        write_folder(h2_dir, h2_structure, h2_calc, args.pseudo_dir, _RELAX_CONFIG_EXTRA)
        ref_rows.append(("02_h2_molecule", h2_structure, "CG relax (Gamma-only)"))
        print_dual(f"  {color_text('[OK]', 'green')} {h2_dir} (Gamma-only, spin-unpolarized "
                    "singlet, relaxes)", f_out)

        h2o_dir = os.path.join(output_root, "03_h2o_molecule")
        h2o_structure = build_h2o_structure(args.vacuum_box)
        h2o_calc = force_system_label(
            force_spin_nonpolarized(force_gamma_kgrid(force_relaxation(site_calc_text))),
            "oer_h2o_molecule")
        write_folder(h2o_dir, h2o_structure, h2o_calc, args.pseudo_dir, _RELAX_CONFIG_EXTRA)
        ref_rows.append(("03_h2o_molecule", h2o_structure, "CG relax (Gamma-only)"))
        print_dual(f"  {color_text('[OK]', 'green')} {h2o_dir} (Gamma-only, spin-unpolarized "
                    "singlet, relaxes)", f_out)

        deformed_dir = os.path.join(output_root, "04_slab_deformed")
        deformed_structure = remove_atoms(relaxed_oh, [oh_o_index, oh_h_index])
        deformed_calc = force_system_label(force_single_point(site_calc_text), "oer_slab_deformed")
        write_folder(deformed_dir, deformed_structure, deformed_calc, args.pseudo_dir,
                     _SINGLE_POINT_CONFIG_EXTRA)
        ref_rows.append(("04_slab_deformed", deformed_structure, "single-point (diagnostic)"))
        print_dual(f"  {color_text('[OK]', 'green')} {deformed_dir} (diagnostic, not used in "
                    "Delta-G directly)", f_out)

        print_dual("", f_out)
        print_table(["Folder", "Atoms", "Formula", "Run type"],
                    [([label, str(len(fdf.atoms)), formula_summary(fdf), run_type], None)
                     for label, fdf, run_type in ref_rows], f_out)

        print_section('[5] BSSE CORRECTION', f_out)
        # Always 3 separate triads, one per intermediate at ITS OWN relaxed
        # geometry -- a single triad shared across all 3 (built once at
        # OOH*'s geometry, the largest fragment, and reused for OH*/O* too)
        # used to be the default ("--bsse-mode shared", now removed): it
        # systematically UNDER-corrected OH*/O* relative to their true
        # per-species BSSE (verified live on a real B9N9/OH-O-OOH system:
        # OH* needed +0.7512 eV, not the shared triad's +0.4552 eV; O*
        # needed +1.0762 eV -- a bare atom's own small DZP basis has the
        # least of its own functions to fall back on, hence the largest
        # ghost-basis correction of the three), inflating eta by ~0.30 V.
        bsse_rows = []
        bsse_rows += write_bsse_triad(output_root, "05_bsse_OH", relaxed_oh, n_substrate,
                                       site_calc_text, args.pseudo_dir, _SINGLE_POINT_CONFIG_EXTRA)
        bsse_rows += write_bsse_triad(output_root, "05_bsse_O", relaxed_o, n_substrate,
                                       o_calc_text, args.pseudo_dir, _SINGLE_POINT_CONFIG_EXTRA)
        bsse_rows += write_bsse_triad(output_root, "05_bsse_OOH", relaxed_ooh, n_substrate,
                                       ooh_calc_text, args.pseudo_dir, _SINGLE_POINT_CONFIG_EXTRA)
        print_dual(f"  {color_text('[OK]', 'green')} 3 separate triads written "
                    "(05_bsse_OH_*, 05_bsse_O_*, 05_bsse_OOH_*), each at its own "
                    "intermediate's geometry", f_out)
        print_dual("", f_out)
        print_table(["Folder", "Atoms", "Formula", "Run type"],
                    [([label, str(len(fdf.atoms)), formula_summary(fdf), run_type], None)
                     for label, fdf, run_type in bsse_rows], f_out)

        print_section('[6] ZPE PREPARATION', f_out)
        print_dual(color_text(
            "[NOTE] H2O's own ZPE/entropy is computed from its G2 gas-phase starting geometry "
            "(build_h2o_structure), not a re-relaxed one -- unlike OH*/O*/OOH*'s hand-built "
            "guesses (which start far from equilibrium and genuinely need CG relaxation first), "
            "ASE's G2 'H2O' entry is already a validated near-equilibrium gas-phase geometry, so "
            "no extra relax-then-reread step is needed here.", 'cyan'), f_out)
        intermediates = [
            ("OH", relaxed_oh, [oh_o_index, oh_h_index], site_calc_text, "oer_zpe_oh"),
            ("O", relaxed_o, [n_substrate], o_calc_text, "oer_zpe_o"),
            ("OOH", relaxed_ooh, [n_substrate, n_substrate + 1, n_substrate + 2], ooh_calc_text,
             "oer_zpe_ooh"),
            ("H2O", h2o_structure, [0, 1, 2], site_calc_text, "oer_zpe_h2o"),
        ]
        if args.zpe_mode == "local":
            zpe_rows = []
            for name, relaxed, local_indices, calc_text, label in intermediates:
                zpe_dir = os.path.join(output_root, f"08_zpe_calc_{name}")
                local_calc = force_system_label(force_single_point(calc_text), label)
                write_local_zpe_folders(zpe_dir, relaxed, local_indices, args.displacement,
                                         local_calc, args.pseudo_dir, label, _SINGLE_POINT_CONFIG_EXTRA)
                n_folders = len(local_indices) * 6
                print_dual(f"  {color_text('[OK]', 'green')} {zpe_dir}/disp_001.."
                            f"disp_{n_folders:03d} ({len(local_indices)} local "
                            "atom(s), decoupled-oscillator approximation)", f_out)
                zpe_rows.append([name, str(len(local_indices)), f"disp_001..disp_{n_folders:03d}",
                                 str(n_folders), f"{args.displacement} Ang"])
            print_dual("", f_out)
            print_table(["Intermediate", "Local atom(s)", "Folder range", "Total folders",
                         "Displacement"], [(row, None) for row in zpe_rows], f_out)
        else:  # full
            print_dual(color_text(
                "[NOTE] 'full' mode needs a full phonon calculation of the clean slab AND all 3 "
                "intermediates (see stb-oerAnalysis --help for why the clean-slab subtraction is "
                "physically necessary) -- roughly 4x the cost of a single full phonon run.",
                'yellow'), f_out)
            from phonopy.interface.siesta import read_siesta
            supercell_matrix = [[args.supercell[0], 0, 0], [0, args.supercell[1], 0], [0, 0, args.supercell[2]]]

            full_targets = [("clean", clean_template, site_calc_text, "oer_zpe_clean")] + [
                (name, relaxed, calc_text, label) for name, relaxed, _idx, calc_text, label in intermediates
            ]
            full_rows = []
            for name, structure_for_phonons, calc_text, label in full_targets:
                ref_fdf_path = os.path.join(output_root, f"09_zpe_calc_{name}", "_reference.fdf")
                os.makedirs(os.path.dirname(ref_fdf_path), exist_ok=True)
                structure_io.write_fdf(structure_for_phonons, ref_fdf_path)
                unitcell = read_siesta(ref_fdf_path)
                phonon, supercells = build_phonon_displacements(unitcell, supercell_matrix, args.displacement)
                zpe_calc = force_system_label(force_single_point(calc_text), label)
                folders, yaml_path = write_displacement_folders(
                    os.path.join(output_root, f"09_zpe_calc_{name}"), phonon, supercells,
                    "structure.fdf", winning_oh_dir + "/calc.fdf", [])
                # phonopy's own siesta writer (write_displacement_folders ->
                # write_siesta) rebuilds each disp-NNN/structure.fdf from
                # PhonopyAtoms, which tracks only atomic numbers -- it always
                # declares bare real-element species labels, never
                # structure_for_phonons' own Stage-1/2 fragment labels
                # ('<real>_slab'/'<real>_ads'). Resolve to the real element
                # before copying, or a fragment-labeled pseudopotential name
                # (e.g. 'O_ads.psf', which doesn't exist) would silently
                # fail to copy.
                real_symbols = sorted({structure_io.real_element(sym, structure_for_phonons.species_meta)
                                        for sym, _ in structure_for_phonons.atoms})
                for d in folders:
                    for sym in real_symbols:
                        copy_pseudo(args.pseudo_dir, sym, d)
                    with open(os.path.join(d, CONFIG_EXTRA_FILE), "w") as f:
                        f.write(_SINGLE_POINT_CONFIG_EXTRA)
                    with open(os.path.join(d, "calc.fdf"), "w") as f:
                        f.write(structure_io.prepend_include(zpe_calc, CONFIG_EXTRA_FILE))
                print_dual(f"  {color_text('[OK]', 'green')} {len(folders)} displacement folder(s) "
                            f"under 09_zpe_calc_{name}/ (symmetry-reduced from "
                            f"{len(supercells)} raw candidate(s))", f_out)
                full_rows.append([name, f"{args.supercell[0]}x{args.supercell[1]}x{args.supercell[2]}",
                                   str(len(structure_for_phonons.atoms)), str(len(supercells)),
                                   str(len(folders))])
            print_dual("", f_out)
            print_table(["Target", "Supercell", "Atoms", "Raw displacements", "Folders (reduced)"],
                        [(row, None) for row in full_rows], f_out)

        print_section('[7] SUMMARY & NEXT STEPS', f_out)
        winning_o_energy = get_free_energy(os.path.join(o_source_dir, args.file))
        winning_ooh_energy = get_free_energy(os.path.join(ooh_source_dir, args.file))
        print_dual(f"Winning OH* site     : {os.path.relpath(winning_oh_dir, output_root)} "
                    f"({winning_oh_energy:.6f} eV)", f_out)
        print_dual(f"O* site (derived)    : {os.path.relpath(o_source_dir, output_root)}"
                    + (f" ({winning_o_energy:.6f} eV)" if winning_o_energy is not None else ""),
                    f_out)
        print_dual(f"OOH* site (derived)  : {os.path.relpath(ooh_source_dir, output_root)}"
                    + (f" ({winning_ooh_energy:.6f} eV)" if winning_ooh_energy is not None else ""),
                    f_out)
        print_dual(f"Reference folders    : {len(ref_rows)} ({', '.join(label for label, _f, _r in ref_rows)})",
                    f_out)
        print_dual(f"BSSE folders         : {len(bsse_rows)} (3 separate triads, one per intermediate)", f_out)
        n_zpe_folders = (sum(int(row[3]) for row in zpe_rows) if args.zpe_mode == "local"
                          else sum(int(row[4]) for row in full_rows))
        print_dual(f"ZPE folders          : {n_zpe_folders} ({args.zpe_mode} mode)", f_out)
        print_dual(f"Pseudo dir           : {args.pseudo_dir or '(reused per-folder from the winning site)'}",
                    f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual("  1. Run SIESTA in every folder written above.", f_out)
        print_dual(f"  2. Once they're done, run: stb-oerAnalysis --directory {output_root}", f_out)

        f_out.write("\nWinning OH* dir : " + os.path.relpath(winning_oh_dir, output_root) + "\n")
        f_out.write("Winning O* dir  : " + os.path.relpath(o_source_dir, output_root) + "\n")
        f_out.write("Winning OOH* dir: " + os.path.relpath(ooh_source_dir, output_root) + "\n")
        f_out.write("ZPE mode        : " + args.zpe_mode + "\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Reference, BSSE and ZPE folders ready for Stage 4 (stb-oerAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
