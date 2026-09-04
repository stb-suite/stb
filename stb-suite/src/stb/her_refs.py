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
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.adsorption_sites import (
    CONFIG_EXTRA_FILE, FIXED_CELL_BLOCK, SPIN_POLARIZED_BLOCK, DIPOLE_CORRECTION_BLOCK,
    VDW_CORRECTION_BLOCK, SINGLE_POINT_BLOCK,
)
from stb.core.bsse import strip_config_extra_include
from stb.core.siesta_log import get_free_energy, get_outcell, check_scf_and_force, report_quality_diagnostics
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders

REPORT_FILE = "her_stage2.txt"
_DEFAULT_LOCAL_DISPLACEMENT_ANG = 0.015
_DEFAULT_VACUUM_BOX_ANG = 15.0
_H2_BOND_LENGTH_ANG = 0.741  # experimental equilibrium bond length -- CG relaxation refines it

_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)
_RELAXED_COORDS_RE = re.compile(r'outcoor:\s*Relaxed atomic coordinates\s*\(fractional\)', re.IGNORECASE)

# Local-mode displacement order (axis, sign) -- matches the row order
# written to disp_NNN/ folders and the sidecar metadata Stage 3 reads back.
_LOCAL_DISPLACEMENTS = [(0, 1.0), (0, -1.0), (1, 1.0), (1, -1.0), (2, 1.0), (2, -1.0)]

# Gamma-only k-grid + CG relaxation, forced into 02_h2_molecule/config_extra.fdf
# only -- H2's own bond length must reach ITS equilibrium (unlike every other
# reference folder here, which is forced single-point via SINGLE_POINT_BLOCK
# instead). HER-specific combination (not shared with stb-adsorb), so kept
# local rather than added to core/adsorption_sites.py's own block constants,
# same "don't pre-extract until a second consumer needs it" policy as the
# rest of core/.
_GAMMA_KGRID_BLOCK = (
    "# Auto-generated -- Gamma-only k-grid, always correct for an isolated\n"
    "# molecule in a vacuum box.\n"
    "kgrid.MonkhorstPack   [1  1  1]\n"
)
_H2_RELAXATION_BLOCK = (
    "# Auto-generated -- H2's own bond length must reach ITS equilibrium (unlike\n"
    "# every other reference folder here, which is forced single-point instead).\n"
    "# Uses MD.Steps, not the older MD.NumCGsteps spelling (deprecated, no longer\n"
    "# honored by current SIESTA).\n"
    "MD.TypeOfRun          CG\n"
    "MD.Steps              200\n"
)
# H2's ground state is an unambiguous closed-shell singlet (2 electrons filling
# the bonding sigma orbital, no radical/open-shell character at all) -- unlike
# the winning SITE (which can genuinely be open-shell, hence SPIN_POLARIZED_BLOCK
# there), forcing Spin polarized on the isolated H2 molecule is not just
# unnecessary but actively wrong: it explicitly OVERRIDES the correct
# non-polarized/restricted state instead of just "costing nothing" the way it
# does for a genuinely closed-shell slab (DIPOLE_CORRECTION_BLOCK/
# SPIN_POLARIZED_BLOCK's own docstrings) -- a spin-polarized SCF can converge
# to a spurious nonzero moment for H2 depending on the initial spin guess,
# silently biasing 0.5*E(H2) in Delta-G_H*. Forcing it explicitly (rather than
# just omitting Spin and relying on SIESTA's own non-polarized default) makes
# the choice self-documenting in config_extra.fdf and immune to whatever the
# user's own --calc template happens to set.
_SPIN_NONPOLARIZED_BLOCK = (
    "# Auto-generated -- forces a spin-unpolarized (restricted) SCF: H2's ground\n"
    "# state is an unambiguous closed-shell singlet, not an open question the\n"
    "# way the adsorbed site can be (see SPIN_POLARIZED_BLOCK elsewhere in this\n"
    "# file) -- a spin-polarized run here risks converging to a spurious nonzero\n"
    "# moment instead of the true singlet.\n"
    "Spin                non-polarized\n"
)

# config_extra.fdf content shared by every single-point derived folder here
# (00_clean_slab, 03_slab_deformed, the BSSE ghost triad, every local/full ZPE
# displacement folder): fixed cell (this is one independent, uncoupled sample
# point, not something SIESTA should be moving on its own) plus the same
# Slab.DipoleCorrection/Spin polarized/DFTD3 the winning site itself needed
# (her.py's write_site_folder, where all three are mandatory, not opt-in) --
# propagating the winning site's own numerical settings into every derived
# single-point folder is required for a physically meaningful energy
# difference (the "level-of-theory propagation" convention documented in
# CLAUDE.md); harmless where the true dipole/moment is already zero
# (DIPOLE_CORRECTION_BLOCK's own docstring -- same reasoning applies to
# SPIN_POLARIZED_BLOCK, converges to zero moment for a genuinely closed-shell
# folder like 00_clean_slab); DFTD3 must match too, since a dispersion
# correction applied on one side of an energy difference but not the other
# would bias it.
_SINGLE_POINT_CONFIG_EXTRA = (FIXED_CELL_BLOCK + DIPOLE_CORRECTION_BLOCK + SPIN_POLARIZED_BLOCK
                               + VDW_CORRECTION_BLOCK + SINGLE_POINT_BLOCK)
_H2_CONFIG_EXTRA = (FIXED_CELL_BLOCK + DIPOLE_CORRECTION_BLOCK + _GAMMA_KGRID_BLOCK
                     + _SPIN_NONPOLARIZED_BLOCK + VDW_CORRECTION_BLOCK + _H2_RELAXATION_BLOCK)


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
    ever appends). `symbol` here is already a Stage-1 fragment label
    ('<real>_slab'/'<real>_ads', see her.py's write_site_folder) rather
    than a bare element symbol, so the real Z is read straight out of
    `species_meta` (already declared for every label present) instead of
    constructing a pymatgen Element from the label text -- Element(symbol)
    would raise on a non-bare label like 'H_ads'.
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
    Reverts the atom's Stage-1 fragment label ('<real>_ads') back to the
    bare real element symbol via structure_io.real_element: with no slab
    atoms left in this single-fragment folder, there is nothing left to
    disambiguate from (same "bare label for a single-fragment folder"
    convention stb-adsorb's own write_reference_folder uses for its
    isolated-adsorbate reference).
    """
    symbol, pos = base_structure.atoms[index]
    real_symbol = structure_io.real_element(symbol, base_structure.species_meta)
    real_z = base_structure.species_meta[symbol]['Z']
    return structure_io.FdfStructure(
        lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
        species=[real_symbol], species_meta={real_symbol: {'id': '1', 'Z': real_z}},
        atoms=[(real_symbol, pos)], coord_format=base_structure.coord_format, raw_lines=[],
    )


def formula_summary(fdf_structure):
    """Returns a compact, human-readable formula string for a written
    folder's [2]/[3] report row, e.g. 'B9N9' or 'B9N9 +H(ghost)' -- real
    atoms grouped by element (via structure_io.real_element, so a Stage-1
    fragment label like 'B_slab' or a ghost label like 'H_ads_ghost' both
    collapse to their real element), ghost atoms (negative Z) called out
    separately since they contribute zero electrons/charge despite sharing
    the real pseudopotential.
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
    sidecar convention as stb-adsorb/stb-raman/stb-ir (menu 4.8/4.11/4.12):
    `config_extra_content` (this folder's own combination of forced
    directives -- single-point vs. relaxation, spin, k-grid, dipole
    correction, fixed cell) is written as-is to config_extra.fdf, and
    `%include config_extra.fdf` is prepended to the UNTOUCHED `calc_text`
    (structure_io.prepend_include) rather than editing directives into it
    in place. Handles both a Stage-1 fragment label
    ('<real>_slab'/'<real>_ads') and a ghost label stacked on top of one
    ('<real>_slab_ghost', from make_ghost_variant) transparently for
    pseudopotential copying: the real element behind ANY label is
    recovered via structure_io.real_element (Z-based, robust to any
    suffix or stack of suffixes -- naive string-slicing off '_ghost'
    alone, this function's previous approach, silently mis-resolved a
    plain fragment label like 'H_ads' to a nonexistent 'H_ads.psf' source
    pseudopotential).
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
    return fdf_structure


def write_local_zpe_folders(zpe_dir, relaxed_structure, h_index, displacement_ang, calc_text,
                             pp_path, config_extra_content):
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
        write_folder(disp_dir, disp_structure, calc_text, pp_path, config_extra_content)

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
        print_dual(f"SIESTA output   : {args.file} (scanned in every site_*/ and read back below)", f_out)
        print_dual(f"Pseudo dir      : {args.pseudo_dir or '(reuse the winning site itself has)'}", f_out)
        print_dual(f"ZPE mode        : {args.zpe_mode}", f_out)
        print_dual(f"Displacement    : {args.displacement} Ang", f_out)
        print_dual(f"Vacuum box      : {args.vacuum_box:.1f} Ang (H2 gas-phase reference)", f_out)
        if args.zpe_mode == "full":
            print_dual(f"Supercell       : {args.supercell[0]} {args.supercell[1]} {args.supercell[2]} "
                        "(--zpe-mode full's phonon calculation)", f_out)
        print_dual(f"Report          : {report_path}", f_out)

        print_section('[1] WINNING SITE', f_out)
        winning_dir, winning_energy, all_results = find_winning_site(sites_root, args.file, f_out)
        n_readable = sum(1 for _label, energy in all_results if energy is not None)
        print_dual(f"  {len(all_results)} site(s) scanned, {n_readable} with a readable FreeEng.",
                    f_out)
        for label, energy in all_results:
            marker = color_text(" <-- winner", 'green') if os.path.join(sites_root, label) == winning_dir else ""
            energy_str = f"{energy:.6f} eV" if energy is not None else "(no energy)"
            print_dual(f"  {label:<28}{energy_str}{marker}", f_out)
        readable = [e for _l, e in all_results if e is not None]
        if len(readable) > 1:
            spread = max(readable) - min(readable)
            print_dual(f"  Energy spread across readable sites: {spread:.4f} eV (max - min).", f_out)
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
            # The winning site's own calc.fdf (written by her.py's
            # write_site_folder) is itself '%include config_extra.fdf' +
            # the untouched user template -- strip that include before
            # using this text as the base for THIS stage's own derived
            # folders, each of which gets its OWN, different
            # config_extra.fdf (see write_folder). A no-op if the site
            # folder predates this convention (plain calc_text already).
            site_calc_text = strip_config_extra_include(f.read())

        print_section('[2] REFERENCE FOLDERS', f_out)
        print_dual("  Every folder below %includes its own config_extra.fdf (fixed cell + the "
                    "winning site's own mandatory Slab.DipoleCorrection/DFTD3, plus single-point "
                    "or the H2 molecule's own Gamma-kgrid/relaxation directives, as appropriate) "
                    "instead of editing your --calc template in place. Spin is the one exception: "
                    "polarized (the site's own setting) everywhere except 02_h2_molecule, forced "
                    "NON-polarized there instead -- H2 is an unambiguous closed-shell singlet.",
                    f_out)
        folder_rows = []  # (label, written_fdf_structure, run_type) -- table at the end

        # 00_clean_slab: pristine geometry (Stage 1's own input, already
        # relaxed by the user BEFORE stb-her), single-point with the
        # winning site's exact numerical settings for consistency.
        clean_template = structure_io.read_fdf(clean_slab_source)
        clean_dir = os.path.join(output_root, "00_clean_slab")
        clean_calc = force_system_label(site_calc_text, "her_clean_slab")
        written = write_folder(clean_dir, clean_template, clean_calc, args.pseudo_dir,
                               _SINGLE_POINT_CONFIG_EXTRA)
        folder_rows.append(("00_clean_slab", written, "single-point"))
        print_dual(f"  {color_text('[OK]', 'green')} {clean_dir}", f_out)

        # 02_h2_molecule: gas-phase CHE reference, RELAXES (not single-point).
        # Derived from the winning site's own calc.fdf (same XC functional/
        # basis/mesh cutoff -- numerical consistency with the slab-side
        # calculations matters for an energy difference), forced to Gamma-only
        # + spin-UNpolarized (H2's closed-shell singlet ground state, unlike the
        # site's own Spin polarized) via config_extra.fdf (_H2_CONFIG_EXTRA).
        h2_dir = os.path.join(output_root, "02_h2_molecule")
        h2_structure = build_h2_structure(args.vacuum_box)
        h2_calc = force_system_label(site_calc_text, "her_h2_molecule")
        written = write_folder(h2_dir, h2_structure, h2_calc, args.pseudo_dir, _H2_CONFIG_EXTRA)
        folder_rows.append(("02_h2_molecule", written, "CG relax (Gamma-only)"))
        print_dual(f"  {color_text('[OK]', 'green')} {h2_dir} (Gamma-only, spin-unpolarized "
                    "singlet, relaxes)", f_out)

        # 03_slab_deformed: winning site's relaxed geometry minus H --
        # diagnostic only (E_deformed - E_clean), not part of the final
        # Delta-G_H* formula.
        deformed_dir = os.path.join(output_root, "03_slab_deformed")
        deformed_structure = remove_atom(relaxed, h_index)
        deformed_calc = force_system_label(site_calc_text, "her_slab_deformed")
        written = write_folder(deformed_dir, deformed_structure, deformed_calc, args.pseudo_dir,
                               _SINGLE_POINT_CONFIG_EXTRA)
        folder_rows.append(("03_slab_deformed", written, "single-point (diagnostic)"))
        print_dual(f"  {color_text('[OK]', 'green')} {deformed_dir} (diagnostic, not used in "
                    "Delta-G_H* itself)", f_out)

        # BSSE counterpoise triad, all single-point, all at the winning
        # site's relaxed geometry.
        ghost_dir = os.path.join(output_root, "04_slab_ghost")
        ghost_variant = make_ghost_variant(relaxed, h_index, n_total)  # ghost H
        ghost_calc = force_system_label(site_calc_text, "her_slab_ghost")
        written = write_folder(ghost_dir, ghost_variant, ghost_calc, args.pseudo_dir,
                               _SINGLE_POINT_CONFIG_EXTRA)
        folder_rows.append(("04_slab_ghost", written, "single-point (BSSE)"))
        print_dual(f"  {color_text('[OK]', 'green')} {ghost_dir}", f_out)

        h_ghost_dir = os.path.join(output_root, "06_h_ghost_slab")
        h_ghost_variant = make_ghost_variant(relaxed, 0, h_index)  # ghost everything but H
        h_ghost_calc = force_system_label(site_calc_text, "her_h_ghost_slab")
        written = write_folder(h_ghost_dir, h_ghost_variant, h_ghost_calc, args.pseudo_dir,
                               _SINGLE_POINT_CONFIG_EXTRA)
        folder_rows.append(("06_h_ghost_slab", written, "single-point (BSSE)"))
        print_dual(f"  {color_text('[OK]', 'green')} {h_ghost_dir}", f_out)

        h_iso_dir = os.path.join(output_root, "07_h_isolated")
        h_iso_structure = isolate_atom(relaxed, h_index)
        h_iso_calc = force_system_label(site_calc_text, "her_h_isolated")
        written = write_folder(h_iso_dir, h_iso_structure, h_iso_calc, args.pseudo_dir,
                               _SINGLE_POINT_CONFIG_EXTRA)
        folder_rows.append(("07_h_isolated", written, "single-point (BSSE)"))
        print_dual(f"  {color_text('[OK]', 'green')} {h_iso_dir}", f_out)

        print_dual("", f_out)
        print_table(["Folder", "Atoms", "Formula", "Run type"],
                    [([label, str(len(fdf.atoms)), formula_summary(fdf), run_type], None)
                     for label, fdf, run_type in folder_rows], f_out)

        print_section('[3] ZPE PREPARATION', f_out)
        n_zpe_folders = 0
        if args.zpe_mode == "standard":
            print_dual("Standard mode -- no folders generated. stb-herAnalysis will use the "
                        "fixed Norskov offset (Delta-ZPE - T*Delta-S ~= +0.24 eV).", f_out)
        elif args.zpe_mode == "local":
            zpe_dir = os.path.join(output_root, "05_zpe_calc")
            local_calc = force_system_label(site_calc_text, "her_zpe")
            write_local_zpe_folders(zpe_dir, relaxed, h_index, args.displacement, local_calc,
                                     args.pseudo_dir, _SINGLE_POINT_CONFIG_EXTRA)
            print_dual(f"  {color_text('[OK]', 'green')} {zpe_dir}/disp_001..disp_006 (6 folders, "
                        f"{n_total} atom(s) each, +/-{args.displacement} Ang finite-difference "
                        "displacements of H only -- partial Hessian, decoupled-oscillator "
                        "approximation, ignores adsorbate-substrate vibrational coupling)", f_out)
            n_zpe_folders = 6
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
            n_supercell_atoms = n_total * args.supercell[0] * args.supercell[1] * args.supercell[2]
            print_dual(f"  Site supercell  : {args.supercell[0]}x{args.supercell[1]}x{args.supercell[2]} "
                        f"({n_supercell_atoms} atom(s)) -- {len(site_supercells)} raw candidate "
                        "displacement(s) (before symmetry reduction).", f_out)
            site_zpe_calc = structure_io.prepend_include(
                force_system_label(site_calc_text, "her_zpe_site"), CONFIG_EXTRA_FILE)
            site_folders, site_yaml = write_displacement_folders(
                os.path.join(output_root, "05_zpe_calc_site"), site_phonon, site_supercells,
                "structure.fdf", winning_dir + "/calc.fdf", [])
            for d in site_folders:
                # phonopy's own siesta writer (write_displacement_folders ->
                # write_siesta) rebuilds each disp-NNN/structure.fdf from
                # PhonopyAtoms, which tracks only atomic numbers -- it always
                # declares bare real-element species labels, never `relaxed`'s
                # Stage-1 fragment labels ('<real>_slab'/'<real>_ads'). Resolve
                # to the real element before copying, or a fragment-labeled
                # pseudopotential name (e.g. 'H_ads.psf', which doesn't exist)
                # would silently fail to copy.
                real_symbols = sorted({structure_io.real_element(sym, relaxed.species_meta)
                                        for sym, _ in relaxed.atoms})
                for sym in real_symbols:
                    copy_pseudo(args.pseudo_dir, sym, d)
                with open(os.path.join(d, CONFIG_EXTRA_FILE), "w") as f:
                    f.write(_SINGLE_POINT_CONFIG_EXTRA)
                with open(os.path.join(d, "calc.fdf"), "w") as f:
                    f.write(site_zpe_calc)
            print_dual(f"  {color_text('[OK]', 'green')} {len(site_folders)} displacement folder(s) "
                        f"under 05_zpe_calc_site/ (symmetry-reduced from the {len(site_supercells)} "
                        "raw candidates above)", f_out)

            clean_fdf_path = os.path.join(output_root, "05_zpe_calc_clean", "_reference.fdf")
            os.makedirs(os.path.dirname(clean_fdf_path), exist_ok=True)
            structure_io.write_fdf(clean_template, clean_fdf_path)
            clean_unitcell = read_siesta(clean_fdf_path)
            clean_phonon, clean_supercells = build_phonon_displacements(
                clean_unitcell, supercell_matrix, args.displacement)
            n_clean_atoms = len(clean_template.atoms)
            n_clean_supercell_atoms = n_clean_atoms * args.supercell[0] * args.supercell[1] * args.supercell[2]
            print_dual(f"  Clean supercell : {args.supercell[0]}x{args.supercell[1]}x{args.supercell[2]} "
                        f"({n_clean_supercell_atoms} atom(s)) -- {len(clean_supercells)} raw candidate "
                        "displacement(s) (before symmetry reduction).", f_out)
            clean_zpe_calc = structure_io.prepend_include(
                force_system_label(site_calc_text, "her_zpe_clean"), CONFIG_EXTRA_FILE)
            clean_folders, clean_yaml = write_displacement_folders(
                os.path.join(output_root, "05_zpe_calc_clean"), clean_phonon, clean_supercells,
                "structure.fdf", winning_dir + "/calc.fdf", [])
            for d in clean_folders:
                symbols = sorted({sym for sym, _ in clean_template.atoms})
                for sym in symbols:
                    copy_pseudo(args.pseudo_dir, sym, d)
                with open(os.path.join(d, CONFIG_EXTRA_FILE), "w") as f:
                    f.write(_SINGLE_POINT_CONFIG_EXTRA)
                with open(os.path.join(d, "calc.fdf"), "w") as f:
                    f.write(clean_zpe_calc)
            print_dual(f"  {color_text('[OK]', 'green')} {len(clean_folders)} displacement folder(s) "
                        f"under 05_zpe_calc_clean/ (symmetry-reduced from the {len(clean_supercells)} "
                        "raw candidates above)", f_out)
            n_zpe_folders = len(site_folders) + len(clean_folders)

        print_section('[4] SUMMARY & NEXT STEPS', f_out)
        print_dual(f"Winning site         : {os.path.basename(winning_dir)}", f_out)
        print_dual(f"Reference folders    : {len(folder_rows)} ({', '.join(label for label, _f, _r in folder_rows)})",
                    f_out)
        print_dual(f"ZPE folders          : {n_zpe_folders} ({args.zpe_mode} mode)", f_out)
        print_dual(f"Pseudo dir           : {args.pseudo_dir or '(reused per-folder from the winning site)'}",
                    f_out)
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
