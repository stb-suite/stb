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

REPORT_FILE = "her_stage2.txt"
_DEFAULT_LOCAL_DISPLACEMENT_ANG = 0.015
_DEFAULT_VACUUM_BOX_ANG = 15.0
_H2_BOND_LENGTH_ANG = 0.741  # experimental equilibrium bond length -- CG relaxation refines it

_KGRID_RE = re.compile(r'kgrid[._]MonkhorstPack\s+\[.*?\]', re.IGNORECASE)
_SPIN_RE = re.compile(r'(Spin\s+)(\S+)', re.IGNORECASE)
_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)
_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_NUMCGSTEPS_RE = re.compile(r'MD\.NumCGsteps\s+\d+', re.IGNORECASE)
_RELAXED_COORDS_RE = re.compile(r'outcoor:\s*Relaxed atomic coordinates\s*\(fractional\)', re.IGNORECASE)

# Local-mode displacement order (axis, sign) -- matches the row order
# written to disp_NNN/ folders and the sidecar metadata Stage 3 reads back.
_LOCAL_DISPLACEMENTS = [(0, 1.0), (0, -1.0), (1, 1.0), (1, -1.0), (2, 1.0), (2, -1.0)]


def force_gamma_kgrid(calc_text):
    """Same substitution as stb-adsorb's own force_gamma_kgrid
    (duplicated, not imported -- HER is self-contained, see her.py's
    module docstring): Gamma-only k-grid, always correct for an isolated
    molecule in a vacuum box.
    """
    new_text, count = _KGRID_RE.subn('kgrid.MonkhorstPack   [1  1  1]', calc_text)
    if count == 0:
        raise ValueError("Could not find a 'kgrid.MonkhorstPack' tag in the calc.fdf template.")
    return new_text


def force_spin_polarized(calc_text):
    """Same as stb-adsorb's own force_spin_polarized (duplicated) --
    forces Spin polarized, appending the tag if absent rather than
    erroring (SIESTA defaults to non-polarized, so an absent tag is a
    normal template).
    """
    new_text, count = _SPIN_RE.subn(r'\g<1>polarized', calc_text)
    if count == 0:
        return calc_text + "\nSpin                polarized\n"
    return new_text


def force_relaxation(calc_text, max_steps=200):
    """Forces 'MD.TypeOfRun CG' + 'MD.NumCGsteps <max_steps>' -- H2's own
    bond length must reach ITS equilibrium (unlike every other reference
    folder here, which is forced single-point via
    core.calc_directives.force_single_point instead). The winning site's
    own calc.fdf template isn't guaranteed to already set MD.TypeOfRun
    (some users configure it externally, or rely on the site relaxation
    being driven some other way) -- explicit is safer than silently
    inheriting whatever (or nothing) the template happens to have.
    """
    new_text, count = _RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _NUMCGSTEPS_RE.subn(f'MD.NumCGsteps         {max_steps}', new_text)
    if count == 0:
        new_text += f"MD.NumCGsteps         {max_steps}\n"
    return new_text


def force_system_label(calc_text, label):
    """Substitutes/appends SystemLabel -- every derived reference folder
    needs its own distinct label so SIESTA's own per-run files (.DM,
    .ion, ...) never collide between folders run from the same working
    tree layout.
    """
    new_text, count = _LABEL_RE.subn(f'SystemLabel {label}', calc_text)
    if count == 0:
        new_text += f"\nSystemLabel {label}\n"
    return new_text


def find_winning_site(sites_root, out_file, f_out):
    """Scans every 'site_*/' folder's FreeEng (core.siesta_log.get_free_energy
    -- NOT a raw 'Total =' line, matching every other energy-difference
    workflow in this suite) and returns (winning_dir, winning_energy,
    all_results) where all_results is [(label, energy_or_None), ...] in
    scan order, for the report.
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
    """Reads the LAST 'outcoor: Relaxed atomic coordinates (fractional)'
    block from a SIESTA .out file and returns an updated FdfStructure
    (same lattice_constant/species/species_meta as `template`, new atomic
    positions from the actual relaxation result) -- SIESTA never rewrites
    structure.fdf in place, so every downstream derived structure (the
    deformed slab, the BSSE ghost triad, the local/full ZPE geometry)
    needs to start from the geometry the relaxation actually reached, not
    Stage 1's pre-relaxation input. Uses core.siesta_log.get_outcell for
    an updated lattice (only different from `template`'s if the
    calc.fdf used MD.VariableCell T -- not the default HER expects, but
    handled anyway rather than silently ignored). Returns None (same
    fail-soft contract as every core.siesta_log parser) if the block
    isn't found -- caller must treat this as "relaxation didn't finish".
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
    the CHE gas-phase reference state (NOT the same thing as an
    'isolated adsorbate' reference: the species adsorbed on the surface
    is ATOMIC H, but the reference reservoir the Delta-G_H* formula
    compares against is MOLECULAR H2 -- see her_analysis.py). Starts at
    the experimental equilibrium bond length; the caller's calc.fdf
    relaxes it further (MD.TypeOfRun CG, not single-point -- H2's own
    bond length should reach ITS equilibrium, unlike every other
    reference folder here which is evaluated at a fixed, already-known
    geometry).
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


def make_ghost_variant(base_structure, ghost_start, ghost_end):
    """Returns a copy of `base_structure` with atoms in [ghost_start,
    ghost_end) turned into ghost species ('<symbol>_ghost', negative Z,
    same real pseudopotential file via copy_pseudo's dest_label) --
    same SIESTA ghost-atom Boys-Bernardi counterpoise convention as
    stb-adsorb's own make_ghost_variant (adsorb.py), duplicated here
    (HER doesn't import from adsorb.py, see her.py's module docstring).
    `base_structure.atoms` must have the slab atoms first and H last
    (guaranteed here: Stage 1 always appends H via
    AdsorbateSiteFinder.add_adsorbate/adsorb_both_surfaces, which only
    ever appends).
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


def remove_atom(base_structure, index):
    """Returns a copy of `base_structure` with the atom at `index`
    (0-based) removed entirely (not ghosted) -- used for '03_slab_deformed'
    (winning site's relaxed geometry minus H).
    """
    new_atoms = [a for i, a in enumerate(base_structure.atoms) if i != index]
    present = {sym for sym, _ in new_atoms}
    species_meta = {k: v for k, v in base_structure.species_meta.items() if k in present}
    species = [s for s in base_structure.species if s in present]
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=species, species_meta=species_meta, atoms=new_atoms,
        coord_format=base_structure.coord_format, raw_lines=[],
    )


def isolate_atom(base_structure, index):
    """Returns a copy of `base_structure` keeping ONLY the atom at
    `index` (same cell/lattice, every other atom removed entirely, not
    ghosted) -- used for '07_h_isolated' (H alone, no slab atoms at all,
    real or ghost -- the reference the H-Ghost-Slab BSSE term needs).
    """
    symbol, pos = base_structure.atoms[index]
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=[symbol], species_meta={symbol: {'id': '1', 'Z': Element(symbol).Z}},
        atoms=[(symbol, pos)], coord_format=base_structure.coord_format, raw_lines=[],
    )


def write_folder(out_dir, fdf_structure, calc_text, pp_path):
    """Writes structure.fdf + calc.fdf + copied pseudos for one derived
    reference folder. Handles ghost species (a '<symbol>_ghost' label
    resolves back to the real element's pseudopotential via
    copy_pseudo's dest_label) transparently.
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    present_labels = sorted({symbol for symbol, _ in fdf_structure.atoms})
    for label in present_labels:
        real_symbol = label[:-len("_ghost")] if label.endswith("_ghost") else label
        copy_pseudo(pp_path, real_symbol, out_dir, dest_label=label)


def write_local_zpe_folders(zpe_dir, relaxed_structure, h_index, displacement_ang, calc_text, pp_path):
    """Writes 6 single-point folders (H displaced +/-x, +/-y, +/-z by
    `displacement_ang` from its relaxed position, every other atom held
    fixed at ITS relaxed position) for the partial-Hessian 'local' ZPE
    mode -- a common surface-science approximation treating the light H
    adsorbate as a decoupled 3-DOF oscillator against a rigid substrate
    (valid when the adsorbate is much lighter than the substrate; ignores
    adsorbate-substrate vibrational coupling, unlike --zpe-mode full).
    Writes a JSON sidecar (zpe_local_meta.json) recording which atom
    index is H and the exact axis/sign/displacement order, so
    stb-herAnalysis doesn't have to re-derive it from folder names.
    """
    os.makedirs(zpe_dir, exist_ok=True)
    symbol, h_frac = relaxed_structure.atoms[h_index]
    h_cart = h_frac @ relaxed_structure.lattice
    inv_lattice = np.linalg.inv(relaxed_structure.lattice)

    for i, (axis, sign) in enumerate(_LOCAL_DISPLACEMENTS, start=1):
        delta_cart = np.zeros(3)
        delta_cart[axis] = sign * displacement_ang
        new_cart = h_cart + delta_cart
        new_frac = new_cart @ inv_lattice
        new_atoms = list(relaxed_structure.atoms)
        new_atoms[h_index] = (symbol, new_frac)
        disp_structure = structure_io.FdfStructure(
            lattice=relaxed_structure.lattice, lattice_constant=relaxed_structure.lattice_constant,
            species=relaxed_structure.species, species_meta=relaxed_structure.species_meta,
            atoms=new_atoms, coord_format=relaxed_structure.coord_format, raw_lines=[],
        )
        disp_dir = os.path.join(zpe_dir, f"disp_{i:03d}")
        write_folder(disp_dir, disp_structure, calc_text, pp_path)

    with open(os.path.join(zpe_dir, "zpe_local_meta.json"), "w") as f:
        json.dump({
            "h_index": h_index,
            "displacement_ang": displacement_ang,
            "order": [{"axis": axis, "sign": sign} for axis, sign in _LOCAL_DISPLACEMENTS],
        }, f)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 3: picks the winning H-adsorption site, builds the H2 "
        "gas-phase reference, the BSSE ghost triad, and the ZPE/entropy calculation folders.", 'bold')}
Scans Stage 1's 'sites/site_*/' for the lowest-FreeEng relaxed site, then writes (all single-point
unless noted): '00_clean_slab/' (pristine slab, same numerical settings as the winning site),
'02_h2_molecule/' (isolated H2 molecule -- RELAXES, since its bond length must reach its own
equilibrium), '03_slab_deformed/' (winning site minus H -- diagnostic only), the BSSE
counterpoise triad ('04_slab_ghost/', '06_h_ghost_slab/', '07_h_isolated/'), and the ZPE
calculation folder(s) for --zpe-mode local/full (nothing for 'standard', which uses the fixed
Norskov offset in Stage 3 instead). Doesn't run SIESTA -- run each folder yourself, then use
stb-herAnalysis.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory her_study --zpe-mode local\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="her_study",
                         help="Root directory written by stb-her (default: her_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each site folder (default: calc.out).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (default: reuse whatever the winning "
                              "site's own folder already has via -- resolved fresh here since "
                              "the H2/ghost/isolated references need their own pseudo copies).")
    parser.add_argument("--zpe-mode", choices=["standard", "local", "full"], default="local",
                         help="ZPE/entropy calculation mode (default: local). 'standard': none "
                              "generated here, stb-herAnalysis uses the fixed Norskov offset "
                              "(+0.24 eV). 'local': 6 finite-difference displacements of H only "
                              "(fast, decoupled-oscillator approximation). 'full': complete "
                              "Phonopy displacement sets for BOTH the winning site and the clean "
                              "slab (needed so stb-herAnalysis can subtract the clean slab's own "
                              "phonon ZPE/entropy -- see stb-herAnalysis --help for why that "
                              "subtraction is physically necessary).")
    parser.add_argument("--displacement", type=float, default=_DEFAULT_LOCAL_DISPLACEMENT_ANG,
                         help="Finite-difference displacement in Ang, for --zpe-mode local's H "
                              f"displacements or --zpe-mode full's Phonopy displacements (default: "
                              f"{_DEFAULT_LOCAL_DISPLACEMENT_ANG}).")
    parser.add_argument("--supercell", type=int, nargs=3, default=[1, 1, 1],
                         help="Supercell dimensions for --zpe-mode full's phonon calculation "
                              "(default: 1 1 1 -- appropriate for a large-enough slab supercell "
                              "already; increase only if your slab's own lateral repeat is small).")
    parser.add_argument("--vacuum-box", type=float, default=_DEFAULT_VACUUM_BOX_ANG,
                         help=f"Cubic box side (Ang) for the isolated H2 reference (default: "
                              f"{_DEFAULT_VACUUM_BOX_ANG}).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-herRefs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("HER WORKFLOW -- STAGE 2: REFERENCES & ZPE PREP", 'bold'))
    print("-" * 60)

    output_root = args.directory
    sites_root = os.path.join(output_root, "sites")
    clean_slab_source = os.path.join(output_root, "clean_slab_source", "structure.fdf")
    if not os.path.isdir(sites_root) or not os.path.isfile(clean_slab_source):
        print(color_text(f"[ERROR] '{sites_root}' or '{clean_slab_source}' not found -- run "
                          "stb-her (Stage 1) first.", 'red'))
        sys.exit(1)

    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== HER STAGE 2 REPORT (REFERENCES & ZPE PREP) =====', 'magenta')}", f_out)

        print_section('[0] RUN METADATA', f_out)
        print_dual(f"Directory       : {output_root}", f_out)
        print_dual(f"ZPE mode        : {args.zpe_mode}", f_out)
        print_dual(f"Displacement    : {args.displacement} Ang", f_out)

        print_section('[1] WINNING SITE', f_out)
        winning_dir, winning_energy, all_results = find_winning_site(sites_root, args.file, f_out)
        for label, energy in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(sites_root, label) == winning_dir else ""
            energy_str = f"{energy:.6f} eV" if energy is not None else "(no energy)"
            print_dual(f"  {label:<28}{energy_str}{marker}", f_out)
        print_dual(f"Winning site    : {os.path.basename(winning_dir)} ({winning_energy:.6f} eV)", f_out)
        report_quality_diagnostics(os.path.basename(winning_dir),
                                    os.path.join(winning_dir, args.file), 0.05, f_out)

        winning_template = structure_io.read_fdf(os.path.join(winning_dir, "structure.fdf"))
        relaxed = read_relaxed_structure(os.path.join(winning_dir, args.file), winning_template)
        if relaxed is None:
            print_dual(color_text(
                f"[ERROR] Could not read relaxed coordinates from '{winning_dir}/{args.file}' -- "
                "did the relaxation finish?", 'red'), f_out)
            sys.exit(1)
        n_total = len(relaxed.atoms)
        h_index = n_total - 1  # H is always appended last by stb-her (AdsorbateSiteFinder.add_adsorbate)

        with open(winning_dir + "/calc.fdf") as f:
            site_calc_text = f.read()

        print_section('[2] REFERENCE FOLDERS', f_out)

        # 00_clean_slab: pristine geometry (Stage 1's own input, already
        # relaxed by the user BEFORE stb-her), single-point with the
        # winning site's exact numerical settings for consistency.
        clean_template = structure_io.read_fdf(clean_slab_source)
        clean_dir = os.path.join(output_root, "00_clean_slab")
        clean_calc = force_system_label(force_single_point(site_calc_text), "her_clean_slab")
        write_folder(clean_dir, clean_template, clean_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {clean_dir}", f_out)

        # 02_h2_molecule: gas-phase CHE reference, RELAXES (not single-point).
        # Derived from the winning site's own calc.fdf (same XC functional/
        # basis/mesh cutoff -- numerical consistency with the slab-side
        # calculations matters for an energy difference), just forced to
        # Gamma-only + spin-polarized.
        h2_dir = os.path.join(output_root, "02_h2_molecule")
        h2_structure = build_h2_structure(args.vacuum_box)
        h2_calc = force_system_label(
            force_spin_polarized(force_gamma_kgrid(force_relaxation(site_calc_text))),
            "her_h2_molecule")
        write_folder(h2_dir, h2_structure, h2_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {h2_dir} (Gamma-only, spin-polarized, relaxes)", f_out)

        # 03_slab_deformed: winning site's relaxed geometry minus H --
        # diagnostic only (E_deformed - E_clean), not part of the final
        # Delta-G_H* formula.
        deformed_dir = os.path.join(output_root, "03_slab_deformed")
        deformed_structure = remove_atom(relaxed, h_index)
        deformed_calc = force_system_label(force_single_point(site_calc_text), "her_slab_deformed")
        write_folder(deformed_dir, deformed_structure, deformed_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {deformed_dir} (diagnostic, not used in "
                    "Delta-G_H* itself)", f_out)

        # BSSE counterpoise triad, all single-point, all at the winning
        # site's relaxed geometry.
        ghost_dir = os.path.join(output_root, "04_slab_ghost")
        ghost_variant = make_ghost_variant(relaxed, h_index, n_total)  # ghost H
        ghost_calc = force_system_label(force_single_point(site_calc_text), "her_slab_ghost")
        write_folder(ghost_dir, ghost_variant, ghost_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {ghost_dir}", f_out)

        h_ghost_dir = os.path.join(output_root, "06_h_ghost_slab")
        h_ghost_variant = make_ghost_variant(relaxed, 0, h_index)  # ghost everything but H
        h_ghost_calc = force_system_label(force_single_point(site_calc_text), "her_h_ghost_slab")
        write_folder(h_ghost_dir, h_ghost_variant, h_ghost_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {h_ghost_dir}", f_out)

        h_iso_dir = os.path.join(output_root, "07_h_isolated")
        h_iso_structure = isolate_atom(relaxed, h_index)
        h_iso_calc = force_system_label(force_single_point(site_calc_text), "her_h_isolated")
        write_folder(h_iso_dir, h_iso_structure, h_iso_calc, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {h_iso_dir}", f_out)

        print_section('[3] ZPE PREPARATION', f_out)
        if args.zpe_mode == "standard":
            print_dual("Standard mode -- no folders generated. stb-herAnalysis will use the "
                        "fixed Norskov offset (Delta-ZPE - T*Delta-S ~= +0.24 eV).", f_out)
        elif args.zpe_mode == "local":
            zpe_dir = os.path.join(output_root, "05_zpe_calc")
            local_calc = force_system_label(force_single_point(site_calc_text), "her_zpe")
            write_local_zpe_folders(zpe_dir, relaxed, h_index, args.displacement, local_calc,
                                     args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} {zpe_dir}/disp_001..disp_006 (H only, "
                        "partial Hessian -- decoupled-oscillator approximation, ignores "
                        "adsorbate-substrate vibrational coupling)", f_out)
        else:  # full
            print_dual(color_text(
                "[NOTE] 'full' mode needs a full phonon calculation of BOTH the winning site AND "
                "the clean slab (see stb-herAnalysis --help for why the clean-slab subtraction is "
                "physically necessary) -- roughly double the cost of a single full phonon run.",
                'yellow'), f_out)

            site_fdf_path = os.path.join(output_root, "05_zpe_calc_site", "_reference.fdf")
            os.makedirs(os.path.dirname(site_fdf_path), exist_ok=True)
            structure_io.write_fdf(relaxed, site_fdf_path)
            from phonopy.interface.siesta import read_siesta
            site_unitcell = read_siesta(site_fdf_path)
            supercell_matrix = [[args.supercell[0], 0, 0], [0, args.supercell[1], 0], [0, 0, args.supercell[2]]]
            site_phonon, site_supercells = build_phonon_displacements(
                site_unitcell, supercell_matrix, args.displacement)
            site_zpe_calc = force_system_label(force_single_point(site_calc_text), "her_zpe_site")
            site_folders, site_yaml = write_displacement_folders(
                os.path.join(output_root, "05_zpe_calc_site"), site_phonon, site_supercells,
                "structure.fdf", winning_dir + "/calc.fdf", [])
            for d in site_folders:
                symbols = sorted({sym for sym, _ in relaxed.atoms})
                for sym in symbols:
                    copy_pseudo(args.pseudo_dir, sym, d)
                with open(os.path.join(d, "calc.fdf"), "w") as f:
                    f.write(site_zpe_calc)
            print_dual(f"  {color_text('[OK]', 'green')} {len(site_folders)} displacement folder(s) "
                        f"under 05_zpe_calc_site/", f_out)

            clean_fdf_path = os.path.join(output_root, "05_zpe_calc_clean", "_reference.fdf")
            os.makedirs(os.path.dirname(clean_fdf_path), exist_ok=True)
            structure_io.write_fdf(clean_template, clean_fdf_path)
            clean_unitcell = read_siesta(clean_fdf_path)
            clean_phonon, clean_supercells = build_phonon_displacements(
                clean_unitcell, supercell_matrix, args.displacement)
            clean_zpe_calc = force_system_label(force_single_point(site_calc_text), "her_zpe_clean")
            clean_folders, clean_yaml = write_displacement_folders(
                os.path.join(output_root, "05_zpe_calc_clean"), clean_phonon, clean_supercells,
                "structure.fdf", winning_dir + "/calc.fdf", [])
            for d in clean_folders:
                symbols = sorted({sym for sym, _ in clean_template.atoms})
                for sym in symbols:
                    copy_pseudo(args.pseudo_dir, sym, d)
                with open(os.path.join(d, "calc.fdf"), "w") as f:
                    f.write(clean_zpe_calc)
            print_dual(f"  {color_text('[OK]', 'green')} {len(clean_folders)} displacement folder(s) "
                        f"under 05_zpe_calc_clean/", f_out)

        print_section('[4] SUMMARY & NEXT STEPS', f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual("  1. Run SIESTA in every folder written above.", f_out)
        print_dual(f"  2. Once they're done, run: stb-herAnalysis --directory {output_root}", f_out)

        f_out.write("\nWinning site   : " + os.path.relpath(winning_dir, output_root) + "\n")
        f_out.write("ZPE mode       : " + args.zpe_mode + "\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Reference and ZPE folders ready for Stage 3 (stb-herAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
