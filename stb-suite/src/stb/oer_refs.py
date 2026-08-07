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
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.calc_directives import force_single_point
from stb.core.siesta_log import get_free_energy, get_outcell, check_scf_and_force, report_quality_diagnostics
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders

REPORT_FILE = "oer_stage3.txt"
_DEFAULT_LOCAL_DISPLACEMENT_ANG = 0.015
_DEFAULT_VACUUM_BOX_ANG = 15.0
_H2_BOND_LENGTH_ANG = 0.741  # experimental equilibrium bond length -- CG relaxation refines it

_KGRID_RE = re.compile(r'kgrid[._]MonkhorstPack\s+\[.*?\]', re.IGNORECASE)
_SPIN_RE = re.compile(r'(Spin\s+)(\S+)', re.IGNORECASE)
_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)
_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_NUMCGSTEPS_RE = re.compile(r'MD\.NumCGsteps\s+\d+', re.IGNORECASE)
_RELAXED_COORDS_RE = re.compile(r'outcoor:\s*Relaxed atomic coordinates\s*\(fractional\)', re.IGNORECASE)


# --- duplicated from her_refs.py / oer_intermediates.py (self-contained stage) -------------

def force_gamma_kgrid(calc_text):
    new_text, count = _KGRID_RE.subn('kgrid.MonkhorstPack   [1  1  1]', calc_text)
    if count == 0:
        raise ValueError("Could not find a 'kgrid.MonkhorstPack' tag in the calc.fdf template.")
    return new_text


def force_spin_polarized(calc_text):
    new_text, count = _SPIN_RE.subn(r'\g<1>polarized', calc_text)
    if count == 0:
        return calc_text + "\nSpin                polarized\n"
    return new_text


def force_relaxation(calc_text, max_steps=200):
    new_text, count = _RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _NUMCGSTEPS_RE.subn(f'MD.NumCGsteps         {max_steps}', new_text)
    if count == 0:
        new_text += f"MD.NumCGsteps         {max_steps}\n"
    return new_text


def force_system_label(calc_text, label):
    new_text, count = _LABEL_RE.subn(f'SystemLabel {label}', calc_text)
    if count == 0:
        new_text += f"\nSystemLabel {label}\n"
    return new_text


def find_winning_site(sites_root, out_file, f_out):
    """Scans every 'site_*' folder under sites_root for the lowest FreeEng
    -- reused for the OH* winner (sites/), and (for --strategy search
    intermediates) the O*/OOH* winner among their own candidates/ folder.
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
    """
    species_meta = dict(base_structure.species_meta)
    new_atoms = []
    for i, (symbol, pos) in enumerate(base_structure.atoms):
        if ghost_start <= i < ghost_end:
            label = f"{symbol}_ghost"
            if label not in species_meta:
                real_z = Element(symbol).Z
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
    {O1,O2,Hnew} (3 atoms), not just 1.
    """
    indices = sorted(indices)
    kept = [base_structure.atoms[i] for i in indices]
    species_meta = {}
    for symbol, _ in kept:
        if symbol not in species_meta:
            species_meta[symbol] = {'id': str(len(species_meta) + 1), 'Z': Element(symbol).Z}
    species = list(dict.fromkeys(sym for sym, _ in kept))
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=species, species_meta=species_meta, atoms=kept,
        coord_format=base_structure.coord_format, raw_lines=[],
    )


def write_folder(out_dir, fdf_structure, calc_text, pp_path):
    """Writes structure.fdf + calc.fdf + copied pseudos for one derived
    reference folder. Handles ghost species ('<symbol>_ghost' resolves
    back to the real element's pseudopotential via copy_pseudo's
    dest_label) transparently. Duplicated from her_refs.py.
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    present_labels = sorted({symbol for symbol, _ in fdf_structure.atoms})
    for label in present_labels:
        real_symbol = label[:-len("_ghost")] if label.endswith("_ghost") else label
        copy_pseudo(pp_path, real_symbol, out_dir, dest_label=label)


def write_bsse_triad(out_root, prefix, base_structure, n_substrate, calc_text, pp_path):
    """Writes a Boys-Bernardi counterpoise QUARTET at base_structure's
    geometry: <prefix>_slab_only/ (slab atoms only, adsorbate REMOVED --
    not ghosted), <prefix>_slab_ghost/ (real slab + ghost adsorbate),
    <prefix>_adsorbate_ghost_slab/ (ghost slab + real adsorbate),
    <prefix>_isolated/ (adsorbate alone). BSSE_slab = E(slab_only) -
    E(slab_ghost) and BSSE_adsorbate = E(isolated) - E(adsorbate_ghost_slab)
    -- both differences MUST compare the SAME geometry across basis sets,
    which is why slab_only is written HERE (at base_structure's own
    geometry) rather than reusing '04_slab_deformed' (always built at
    OH*'s geometry): in --bsse-mode shared, the triad is built at OOH*'s
    geometry, which generally differs from OH*'s -- reusing the
    diagnostic-only deformed-slab folder there would silently compare two
    different geometries and corrupt the counterpoise correction.
    Factored into one function since OER calls it either once
    (--bsse-mode shared, at the OOH* geometry -- the largest fragment) or
    three times (--bsse-mode full, once per intermediate at ITS OWN
    geometry).
    """
    n_total = len(base_structure.atoms)
    slab_only = remove_atoms(base_structure, list(range(n_substrate, n_total)))
    write_folder(os.path.join(out_root, f"{prefix}_slab_only"), slab_only,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_slab_only"), pp_path)

    slab_ghost = make_ghost_variant(base_structure, n_substrate, n_total)
    write_folder(os.path.join(out_root, f"{prefix}_slab_ghost"), slab_ghost,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_slab_ghost"), pp_path)

    ads_ghost_slab = make_ghost_variant(base_structure, 0, n_substrate)
    write_folder(os.path.join(out_root, f"{prefix}_adsorbate_ghost_slab"), ads_ghost_slab,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_adsorbate_ghost_slab"),
                 pp_path)

    isolated = isolate_atoms(base_structure, list(range(n_substrate, n_total)))
    write_folder(os.path.join(out_root, f"{prefix}_isolated"), isolated,
                 force_system_label(force_single_point(calc_text), f"oer_{prefix}_isolated"), pp_path)


def write_local_zpe_folders(zpe_dir, relaxed_structure, local_indices, displacement_ang, calc_text,
                             pp_path, system_label):
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
    local_symbols = [relaxed_structure.atoms[i][0] for i in local_indices]

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
        write_folder(disp_dir, disp_structure, calc_text, pp_path)

    with open(os.path.join(zpe_dir, "zpe_local_meta.json"), "w") as f:
        json.dump({
            "local_indices": list(local_indices),
            "local_symbols": local_symbols,
            "displacement_ang": displacement_ang,
            "order": order,
            "system_label": system_label,
        }, f)


def _locate_intermediate(output_root, name, strategy, out_file, f_out):
    """Locates the final relaxed geometry+calc.fdf for O*/OOH* ('name' is
    'o' or 'ooh'), regardless of --o-strategy/--ooh-strategy: 'derived'
    reads the single 'intermediates/<name>_star/' folder directly;
    'search' picks the winner among 'intermediates/<name>_candidates/
    site_*/' via find_winning_site (same logic already used for OH*
    itself). Returns (relaxed_structure, source_dir, calc_text).
    """
    if strategy == "derived":
        source_dir = os.path.join(output_root, "intermediates", f"{name}_star")
        if not os.path.isdir(source_dir):
            print_dual(color_text(f"[ERROR] '{source_dir}' not found -- run stb-oerIntermediates "
                                   "(Stage 2) first.", 'red'), f_out)
            sys.exit(1)
    else:
        candidates_root = os.path.join(output_root, "intermediates", f"{name}_candidates")
        if not os.path.isdir(candidates_root):
            print_dual(color_text(f"[ERROR] '{candidates_root}' not found -- run "
                                   "stb-oerIntermediates (Stage 2) first.", 'red'), f_out)
            sys.exit(1)
        source_dir, energy, all_results = find_winning_site(candidates_root, out_file, f_out)
        for label, e in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(candidates_root, label) == source_dir else ""
            e_str = f"{e:.6f} eV" if e is not None else "(no energy)"
            print_dual(f"  {label:<28}{e_str}{marker}", f_out)

    template = structure_io.read_fdf(os.path.join(source_dir, "structure.fdf"))
    relaxed = read_relaxed_structure(os.path.join(source_dir, out_file), template)
    if relaxed is None:
        print_dual(color_text(
            f"[ERROR] Could not read relaxed coordinates from '{source_dir}/{out_file}' -- did "
            "the relaxation finish?", 'red'), f_out)
        sys.exit(1)
    report_quality_diagnostics(os.path.basename(source_dir), os.path.join(source_dir, out_file),
                                0.05, f_out)
    with open(os.path.join(source_dir, "calc.fdf")) as f:
        calc_text = f.read()
    return relaxed, source_dir, calc_text


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 3 of 4: builds the H2/H2O gas-phase references, the BSSE "
        "counterpoise correction, and the ZPE/entropy calculation folders for all 3 intermediates "
        "(OH*, O*, OOH*).", 'bold')}
Locates the final relaxed geometry of OH* (Stage 1's winning site), O*, and OOH* (Stage 2's
derived-and-relaxed or independently-searched-and-relaxed intermediates), then writes (all
single-point unless noted): '00_clean_slab/' (pristine slab), '02_h2_molecule/' and
'03_h2o_molecule/' (RELAX, gas-phase CHE references), '04_slab_deformed/' (OH*'s geometry with
both adsorbate atoms removed -- diagnostic only), the BSSE counterpoise triad(s)
(--bsse-mode shared: one triad at OOH*'s geometry, reused for all 3 intermediates; --bsse-mode
full: one triad per intermediate, 9 folders total), and the ZPE calculation folder(s)
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
    parser.add_argument("--bsse-mode", choices=["shared", "full"], default="shared",
                         help="BSSE counterpoise scope (default: shared). 'shared': one triad "
                              "built at OOH*'s geometry (the largest fragment), reused for all 3 "
                              "intermediates -- cheaper, and the BSSE bias mostly cancels in the "
                              "O*-OH* and OOH*-O* energy differences anyway. 'full': a separate "
                              "triad per intermediate (9 folders), for per-species rigor.")
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

    o_strategy, ooh_strategy = None, None
    with open(stage2_report) as f:
        for line in f:
            if line.startswith("O* strategy"):
                o_strategy = line.split(":", 1)[1].strip()
            elif line.startswith("OOH* strategy"):
                ooh_strategy = line.split(":", 1)[1].strip()
    if o_strategy is None or ooh_strategy is None:
        print(color_text(f"[ERROR] Could not recover O*/OOH* strategy from '{stage2_report}'.", 'red'))
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
        print_dual(f"BSSE mode       : {args.bsse_mode}", f_out)
        print_dual(f"ZPE mode        : {args.zpe_mode}", f_out)
        print_dual(f"Displacement    : {args.displacement} Ang", f_out)

        print_section('[1] WINNING OH* SITE', f_out)
        winning_oh_dir, winning_oh_energy, all_results = find_winning_site(sites_root, args.file, f_out)
        for label, energy in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(sites_root, label) == winning_oh_dir else ""
            energy_str = f"{energy:.6f} eV" if energy is not None else "(no energy)"
            print_dual(f"  {label:<28}{energy_str}{marker}", f_out)
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
            site_calc_text = f.read()

        print_section('[2] O* FINAL GEOMETRY', f_out)
        relaxed_o, o_source_dir, o_calc_text = _locate_intermediate(
            output_root, "o", o_strategy, args.file, f_out)
        print_dual(f"Source: {o_source_dir} ({o_strategy})", f_out)

        print_section('[3] OOH* FINAL GEOMETRY', f_out)
        relaxed_ooh, ooh_source_dir, ooh_calc_text = _locate_intermediate(
            output_root, "ooh", ooh_strategy, args.file, f_out)
        print_dual(f"Source: {ooh_source_dir} ({ooh_strategy})", f_out)
        n_ooh_total = len(relaxed_ooh.atoms)

        print_section('[4] REFERENCE FOLDERS', f_out)

        clean_template = structure_io.read_fdf(clean_slab_source)
        clean_dir = os.path.join(output_root, "00_clean_slab")
        clean_calc = force_system_label(force_single_point(site_calc_text), "oer_clean_slab")
        write_folder(clean_dir, clean_template, clean_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {clean_dir}", f_out)

        h2_dir = os.path.join(output_root, "02_h2_molecule")
        h2_structure = build_h2_structure(args.vacuum_box)
        h2_calc = force_system_label(
            force_spin_polarized(force_gamma_kgrid(force_relaxation(site_calc_text))),
            "oer_h2_molecule")
        write_folder(h2_dir, h2_structure, h2_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {h2_dir} (Gamma-only, spin-polarized, relaxes)", f_out)

        h2o_dir = os.path.join(output_root, "03_h2o_molecule")
        h2o_structure = build_h2o_structure(args.vacuum_box)
        h2o_calc = force_system_label(
            force_spin_polarized(force_gamma_kgrid(force_relaxation(site_calc_text))),
            "oer_h2o_molecule")
        write_folder(h2o_dir, h2o_structure, h2o_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {h2o_dir} (Gamma-only, spin-polarized, relaxes)", f_out)

        deformed_dir = os.path.join(output_root, "04_slab_deformed")
        deformed_structure = remove_atoms(relaxed_oh, [oh_o_index, oh_h_index])
        deformed_calc = force_system_label(force_single_point(site_calc_text), "oer_slab_deformed")
        write_folder(deformed_dir, deformed_structure, deformed_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {deformed_dir} (diagnostic, not used in "
                    "Delta-G directly)", f_out)

        print_section('[5] BSSE CORRECTION', f_out)
        if args.bsse_mode == "shared":
            write_bsse_triad(output_root, "05_bsse", relaxed_ooh, n_substrate, ooh_calc_text,
                              args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} 05_bsse_slab_only/, 05_bsse_slab_ghost/, "
                        "05_bsse_adsorbate_ghost_slab/, 05_bsse_isolated/ (built at OOH*'s "
                        "geometry, reused for OH*/O*/OOH* alike)", f_out)
        else:
            write_bsse_triad(output_root, "05_bsse_OH", relaxed_oh, n_substrate, site_calc_text,
                              args.pseudo_dir)
            write_bsse_triad(output_root, "05_bsse_O", relaxed_o, n_substrate, o_calc_text,
                              args.pseudo_dir)
            write_bsse_triad(output_root, "05_bsse_OOH", relaxed_ooh, n_substrate, ooh_calc_text,
                              args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} 3 separate triads written "
                        "(05_bsse_OH_*, 05_bsse_O_*, 05_bsse_OOH_*), each at its own "
                        "intermediate's geometry", f_out)

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
            for name, relaxed, local_indices, calc_text, label in intermediates:
                zpe_dir = os.path.join(output_root, f"08_zpe_calc_{name}")
                local_calc = force_system_label(force_single_point(calc_text), label)
                write_local_zpe_folders(zpe_dir, relaxed, local_indices, args.displacement,
                                         local_calc, args.pseudo_dir, label)
                print_dual(f"  {color_text('[OK]', 'green')} {zpe_dir}/disp_001.."
                            f"disp_{len(local_indices) * 6:03d} ({len(local_indices)} local "
                            "atom(s), decoupled-oscillator approximation)", f_out)
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
                symbols = sorted({sym for sym, _ in structure_for_phonons.atoms})
                for d in folders:
                    for sym in symbols:
                        copy_pseudo(args.pseudo_dir, sym, d)
                    with open(os.path.join(d, "calc.fdf"), "w") as f:
                        f.write(zpe_calc)
                print_dual(f"  {color_text('[OK]', 'green')} {len(folders)} displacement folder(s) "
                            f"under 09_zpe_calc_{name}/", f_out)

        print_section('[7] SUMMARY & NEXT STEPS', f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual("  1. Run SIESTA in every folder written above.", f_out)
        print_dual(f"  2. Once they're done, run: stb-oerAnalysis --directory {output_root}", f_out)

        f_out.write("\nWinning OH* dir : " + os.path.relpath(winning_oh_dir, output_root) + "\n")
        f_out.write("Winning O* dir  : " + os.path.relpath(o_source_dir, output_root) + "\n")
        f_out.write("Winning OOH* dir: " + os.path.relpath(ooh_source_dir, output_root) + "\n")
        f_out.write("BSSE mode       : " + args.bsse_mode + "\n")
        f_out.write("ZPE mode        : " + args.zpe_mode + "\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Reference, BSSE and ZPE folders ready for Stage 4 (stb-oerAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
