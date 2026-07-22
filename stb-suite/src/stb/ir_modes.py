#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import glob
import json
import shutil
import argparse
from datetime import datetime
import numpy as np
import yaml
from ase import Atoms
from ase.io import write as ase_write
from phonopy.interface.siesta import write_siesta
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.calc_directives import force_single_point, force_born_charge_run
from stb.core.pseudopotentials import get_required_pseudos, resolve_pseudo_source
from stb.core import kspace
from stb.core.phonon_workflow import (
    detect_system_label, load_phonon_with_force_constants, get_gamma_modes, displace_along_mode,
    mode_eigendisplacement,
)
from stb.core.ir_symmetry import classify_ir_modes

# Same default as raman_modes.py -- see that module's own comment for the
# full reasoning (verified live with a relaxed H2O molecule via MACE-MP-0:
# rotational bands were ~-0.01 to 0.22 THz, real vibrational modes 44.5+
# THz -- 2.0 THz keeps a ~10x margin above observed rotational noise while
# staying well under realistic low-frequency vibrational modes).
_DEFAULT_ROTATIONAL_MODE_TOL_THZ = 2.0

REPORT_FILE = "ir_stage2.txt"

# Default Berry-phase k-point grid for the bulk (3D) Born-effective-charge
# path -- row-major 3x3, one row per lattice vector (diagonal = 1D
# line-integral k-points, off-diagonal = 2D surface-integral mesh), same
# example shape as SIESTA's own PolarizationGrids tutorial material.
_DEFAULT_POLARIZATION_GRID = [10, 4, 4, 4, 10, 4, 4, 4, 10]
_DEFAULT_FC_DISPL_BOHR = 0.02


def write_ir_folder(out_dir, atoms, structure_filename, calc_text, pseudos):
    """Writes one SIESTA input folder: structure + calc.fdf + copied
    pseudos -- used for both the non-bulk path's +/-delta dipole
    displacement folders and the bulk path's single equilibrium
    Born-effective-charge folder (same folder shape either way, just a
    different calc_text/atoms).
    """
    os.makedirs(out_dir, exist_ok=True)
    write_siesta(os.path.join(out_dir, os.path.basename(structure_filename)), atoms)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    for pseudo_path in pseudos:
        shutil.copy(pseudo_path, os.path.join(out_dir, os.path.basename(pseudo_path)))


def _phonopy_atoms_to_ase(patoms, internal_to_angstrom):
    """Same helper as raman_modes.py -- PhonopyAtoms (internal units) ->
    ase.Atoms in real Angstrom, ready for ase.io.write.
    """
    return Atoms(symbols=patoms.symbols,
                 positions=np.array(patoms.positions) * internal_to_angstrom,
                 cell=np.array(patoms.cell) * internal_to_angstrom,
                 pbc=True)


def write_mode_animation(phonon, band_index, amplitude_ang, internal_to_angstrom, out_path, n_frames=20):
    """Same as raman_modes.py's write_mode_animation -- a looping .axsf
    animation of one Gamma-point mode's eigendisplacement, independent of
    which IR path (dipole-derivative or Born-charge) is being used, since
    it's pure visualization of the phonon mode itself.
    """
    amplitudes = amplitude_ang * np.sin(2 * np.pi * np.arange(n_frames) / n_frames)
    frames = []
    for amp in amplitudes:
        displaced = displace_along_mode(phonon, band_index, float(amp), internal_to_angstrom, sign=1.0)
        frames.append(_phonopy_atoms_to_ase(displaced, internal_to_angstrom))
    ase_write(out_path, frames, format='xsf')


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 3: builds FORCE_SETS from Stage 1's phonon "
        "displacements, identifies the Gamma-point vibrational modes, and generates whichever "
        "follow-up SIESTA folder(s) the IR dipole-derivative calculation needs.", 'bold')}
The path is AUTO-SELECTED by structure dimensionality (core.kspace.detect_vacuum_axes on the
primitive cell), no flag needed:

  - Non-bulk (0D/1D/2D -- at least one vacuum-padded axis): SIESTA prints the total dipole
    moment automatically for a non-periodic direction, no special fdf block needed. Writes
    exactly 2 folders per selected mode (dipole_disp/mode_XX_plus/, mode_XX_minus/) -- a
    simple single-point SCF, no per-axis probing (unlike stb-ramanModes' Optical.Vector sweep):
    the dipole derivative dmu/dQ already comes out as a full 3-vector from one +/-delta pair.

  - Bulk (3D periodic, no vacuum axis): a naive dipole moment is gauge-ambiguous for a fully
    periodic system. Writes exactly ONE equilibrium folder (born_charge_disp/equilibrium/,
    the undisplaced structure) using SIESTA's own native Born-effective-charge automation
    (MD.TypeOfRun FC + BornCharge T + %block PolarizationGrids) -- independent of how many
    modes are selected. Also writes a small eigendisplacement JSON sidecar per selected mode
    (born_charge_disp/mode_XX_eigendisplacement.json) so stb-irAnalysis can combine it with the
    Born-charge tensors without reloading Phonopy itself.

Doesn't run SIESTA itself -- run each folder's calculation yourself, then use stb-irAnalysis.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory ir_study --calc calc.fdf\n"
               "  %(prog)s --directory ir_study --calc calc.fdf --modes 1 2 3\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="ir_study",
                        help="Root directory written by stb-ir (Stage 1) -- must contain "
                             "'phonon_disp/' (default: ir_study).")
    parser.add_argument("-l", "--label", type=str, default=None,
                        help="SystemLabel used in Stage 1's calc.fdf (default: auto-detected).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default=None, metavar="DIR",
                        help="Pseudopotentials source: a bundled bank or a folder path. Default: "
                             "reuse whatever Stage 1 already resolved and copied into its "
                             "phonon_disp/disp-*/ folders (the normal, self-contained case -- "
                             "no need to point at a source a second time). Only needed if that "
                             "folder is missing pseudopotentials for some reason (e.g. hand-built "
                             "phonon_disp/ instead of a real stb-ir run).")
    parser.add_argument("-c", "--calc", type=str, required=True,
                        help="calc.fdf template for the follow-up SIESTA calculation(s) -- can be "
                             "the same file used in Stage 1 or a different one, since these are "
                             "now single-point evaluations on the small unit cell, not the phonon "
                             "supercell.")
    parser.add_argument("--modes", type=int, nargs='+', default=None, metavar="N",
                        help="1-based mode index/indices to process (in ascending-frequency "
                             "order, acoustic modes already excluded -- see Stage 2's own "
                             "[1] PHONON MODES table for the numbering). Default: every "
                             "non-acoustic mode. On the bulk path this only controls which modes "
                             "appear in the report/MODE_TABLE and get an eigendisplacement "
                             "sidecar -- the single equilibrium folder is written regardless.")
    parser.add_argument("--freq-min", type=float, default=None,
                        help="Skip modes below this frequency (THz).")
    parser.add_argument("--freq-max", type=float, default=None,
                        help="Skip modes above this frequency (THz).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum gap (Ang) along an axis to consider it vacuum-padded, "
                             "same convention as stb-ir (default: 10.0). Used both for the 0D "
                             "extra-trivial-mode detection AND to auto-select the bulk vs. "
                             "non-bulk IR path.")
    parser.add_argument("--rotational-mode-tol", type=float, default=_DEFAULT_ROTATIONAL_MODE_TOL_THZ,
                        help="0D (molecule) only: a mode within the first 3 bands after the "
                             "translations is treated as free rotation (trivial, excluded) if "
                             "|frequency| is below this (THz) (default: "
                             f"{_DEFAULT_ROTATIONAL_MODE_TOL_THZ}). Lower this if your molecule "
                             "has genuine low-frequency skeletal/torsional vibrational modes you "
                             "want to make sure aren't mistaken for rotation.")
    parser.add_argument("--displacement", type=float, default=0.02,
                        help="Finite-difference displacement (Ang) along each mode's "
                             "eigendisplacement for the non-bulk dipole derivative, and the "
                             "default amplitude for --export-animations (default: 0.02). Has no "
                             "effect on the bulk path's folder (it's the equilibrium structure, "
                             "never displaced).")
    parser.add_argument("--use-symmetry", action="store_true",
                        help="Skip modes that group theory guarantees are IR-INACTIVE (their "
                             "irreducible representation at Gamma isn't contained in the vector "
                             "(x,y,z) representation -- their dipole derivative is exactly zero, "
                             "not just numerically small). Works for any of the 32 point groups. "
                             "Only applies to the default mode selection (every non-acoustic "
                             "mode) -- an explicit --modes list is never second-guessed. On the "
                             "bulk path this only reduces what's listed/sidecar'd, never the "
                             "folder count (always 1). Needs a PRIMITIVE cell (see stb-unitcell "
                             "--mode primitive); falls back to running every mode, with a "
                             "warning, if that's not available.")
    parser.add_argument("--skip-degenerate", action="store_true",
                        help="Non-bulk path only: for a degenerate group of modes (2 or more "
                             "bands sharing one irrep), write dipole-displacement folders for "
                             "only ONE representative band and skip the rest -- the IR intensity "
                             "|dmu/dQ|^2 is a rotational invariant, identical for every "
                             "degenerate partner regardless of phonopy's arbitrary eigenvector "
                             "choice within the subspace. stb-irAnalysis reuses the "
                             "representative's value automatically. No effect on the bulk path "
                             "(already exactly 1 folder regardless of mode count -- nothing to "
                             "save by skipping). Needs a PRIMITIVE cell, same as --use-symmetry.")
    parser.add_argument("--export-animations", action="store_true",
                        help="Write a looping animation of each selected mode's eigendisplacement "
                             "(a smooth 0/+A/0/-A/0 sweep) as a mode_XX_animation.axsf file -- "
                             "open in XCrySDen/VESTA to see which atoms move before interpreting "
                             "the spectrum. Same amplitude as --displacement.")
    parser.add_argument("--animation-frames", type=int, default=20,
                        help="Frames per mode animation, with --export-animations (default: 20).")
    parser.add_argument("--fc-displ", type=float, default=_DEFAULT_FC_DISPL_BOHR, metavar="BOHR",
                        help="Bulk path only: MD.FCDispl for SIESTA's own native Born-effective-"
                             f"charge finite-displacement automation, in Bohr (default: "
                             f"{_DEFAULT_FC_DISPL_BOHR}). Has no effect on the non-bulk path.")
    parser.add_argument("--polarization-grid", type=int, nargs=9,
                        default=_DEFAULT_POLARIZATION_GRID, metavar="N",
                        help="Bulk path only: the 9 integers (row-major 3x3, one row per lattice "
                             "vector) for SIESTA's %%block PolarizationGrids -- diagonal entries "
                             "are 1D line-integral k-point counts, off-diagonal entries are the "
                             "2D surface-integral mesh (default: "
                             f"{' '.join(str(v) for v in _DEFAULT_POLARIZATION_GRID)}). Manual, "
                             "not auto-generated from a density heuristic -- consult the SIESTA "
                             "manual for your structure. Has no effect on the non-bulk path.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-irModes {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("IR SPECTRUM WORKFLOW -- STAGE 2: MODES & DIPOLE/BORN-CHARGE DISPLACEMENTS", 'bold'))
    print("-" * 60)

    phonon_dir = os.path.join(args.directory, "phonon_disp")
    if not os.path.isdir(phonon_dir):
        print(color_text(f"[ERROR] '{phonon_dir}' not found -- run stb-ir (Stage 1) first.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
        sys.exit(1)

    yaml_file = os.path.join(phonon_dir, "phonopy_disp.yaml")
    if not os.path.exists(yaml_file):
        print(color_text(f"[ERROR] '{yaml_file}' not found -- did Stage 1 finish successfully?", 'red'))
        sys.exit(1)
    with open(yaml_file) as f:
        has_embedded_fc = "force_constants" in (yaml.safe_load(f) or {})

    system_label, label_source = None, None
    if not has_embedded_fc:
        if args.label is not None:
            system_label, label_source = args.label, "manual (-l/--label)"
        else:
            system_label = detect_system_label(phonon_dir) or "siesta"
            label_source = "auto-detected from calc.fdf"
            print(f"[INFO] Auto-detected SystemLabel '{system_label}' from calc.fdf "
                  "(pass -l/--label to override).")

    with open(args.calc) as f:
        calc_template = f.read()

    if args.pseudo_dir is not None:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    output_root = args.directory
    dipole_root = os.path.join(output_root, "dipole_disp")
    born_charge_root = os.path.join(output_root, "born_charge_disp")
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== IR STAGE 2 REPORT (MODES & DISPLACEMENTS) =====', 'magenta')}", f_out)

        print_section('[0] RUN METADATA', f_out)
        print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Directory         : {output_root}", f_out)
        print_dual(f"SystemLabel       : {'N/A (ML-computed force constants)' if has_embedded_fc else f'{system_label} ({label_source})'}", f_out)
        print_dual(f"Calc template     : {args.calc}", f_out)
        print_dual(f"Displacement      : {args.displacement} Ang (non-bulk path only)", f_out)
        print_dual(f"FC displacement   : {args.fc_displ} Bohr (bulk path only)", f_out)
        print_dual(f"Polarization grid : {' '.join(str(v) for v in args.polarization_grid)} "
                    "(bulk path only)", f_out)

        print_section('[1] PHONON MODES AT GAMMA', f_out)
        phonon, internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            phonon_dir, system_label, has_embedded_fc, f_out)
        try:
            unique_elements = list(set(phonon.primitive.symbols))
            lattice_ang = np.array(phonon.primitive.cell) * internal_to_angstrom
            frac_coords = np.array(phonon.primitive.scaled_positions)
            vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
            is_0d = all(vacuum_axes)
            is_bulk = not any(vacuum_axes)
            frequencies, mode_band_indices, n_extra_trivial = get_gamma_modes(
                phonon, exclude_acoustic=True,
                extra_trivial_tol_thz=args.rotational_mode_tol if is_0d else None)
            mode_symmetries, point_group, symmetry_error = classify_ir_modes(phonon)
        finally:
            os.chdir(original_dir)

        band_to_symmetry = {}
        for ms in mode_symmetries:
            for b in ms.band_indices:
                band_to_symmetry[b] = ms

        print_dual(f"Dimensionality    : {kspace.dimensionality_label(vacuum_axes)}", f_out)
        print_dual(f"Vacuum axes (abc) : {vacuum_axes[0]} {vacuum_axes[1]} {vacuum_axes[2]}", f_out)
        print_dual("Lattice vectors (Ang, for Stage 3's dipole-axis alignment check):", f_out)
        for name, vec in zip("abc", lattice_ang):
            print_dual(f"  {name}  {vec[0]:.6f}  {vec[1]:.6f}  {vec[2]:.6f}", f_out)
        if is_bulk:
            print_dual(color_text(
                "[Path] Bulk (3D periodic) -- one equilibrium Born-effective-charge SIESTA run "
                "(MD.TypeOfRun FC + BornCharge T + PolarizationGrids) covers every selected "
                "mode's IR intensity.", 'cyan'), f_out)
        else:
            print_dual(color_text(
                "[Path] Non-bulk (0D/1D/2D) -- a +/-delta dipole-moment displacement pair per "
                "selected mode (plain single-point SCF, SIESTA prints the total dipole "
                "automatically).", 'cyan'), f_out)
        if n_extra_trivial:
            print_dual(color_text(
                f"[NOTE] 0D structure detected -- excluded {n_extra_trivial} extra trivial "
                "mode(s) (free rotation) beyond the usual 3 translations. These aren't real "
                "vibrations and have no IR activity to speak of.", 'yellow'), f_out)

        n_imaginary = int((frequencies < 0).sum())
        print_dual(f"Non-acoustic Gamma modes : {len(frequencies)}", f_out)
        if n_imaginary:
            print_dual(color_text(
                f"[WARNING] {n_imaginary} mode(s) have imaginary (negative) frequency -- the "
                "structure/supercell may not be at a real energy minimum. Their IR intensity "
                "is not physically meaningful.", 'yellow'), f_out)
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            flag = color_text(" [IMAGINARY]", 'red') if freq < 0 else ""
            ms = band_to_symmetry.get(int(band_idx))
            sym_note = ""
            if args.use_symmetry and ms is not None:
                label_str = f" ({ms.label})" if ms.label else ""
                sym_note = (color_text(f"{label_str} [IR-active]", 'green') if ms.is_ir_active
                            else color_text(f"{label_str} [symmetry-forbidden]", 'yellow'))
            print_dual(f"  mode {k:3d} (band {int(band_idx):3d}) : {freq:10.4f} THz{flag} {sym_note}", f_out)

        if args.use_symmetry:
            print_section('[1b] SYMMETRY ANALYSIS', f_out)
            print_dual(f"Point group       : {point_group}", f_out)
            if symmetry_error:
                print_dual(color_text(
                    f"[WARNING] Symmetry classification unavailable ({symmetry_error}) -- "
                    "running every mode, same as without --use-symmetry. If your structure is "
                    "a conventional (non-primitive) cell, try reducing it first with "
                    "stb-unitcell --mode primitive.", 'yellow'), f_out)
            else:
                n_forbidden = sum(1 for band_idx in mode_band_indices
                                   if band_to_symmetry.get(int(band_idx)) is not None
                                   and not band_to_symmetry[int(band_idx)].is_ir_active)
                n_unknown = sum(1 for band_idx in mode_band_indices
                                 if band_to_symmetry.get(int(band_idx)) is None)
                print_dual(f"Symmetry-forbidden (IR-inactive) : {n_forbidden}/{len(mode_band_indices)}", f_out)
                if n_unknown:
                    print_dual(color_text(
                        f"[NOTE] {n_unknown} mode(s) could not be classified (label matching "
                        "failed) -- kept, never auto-skipped when uncertain.", 'yellow'), f_out)
                if args.modes is not None:
                    print_dual(color_text(
                        "--modes was given explicitly -- symmetry filtering is informational "
                        "only here, no mode is auto-skipped.", 'yellow'), f_out)

        print_section('[2] PSEUDOPOTENTIALS', f_out)
        print_dual(f"Elements needed   : {', '.join(sorted(unique_elements))}", f_out)
        if args.pseudo_dir is not None:
            pseudo_source = args.pseudo_dir
            print_dual(f"Source            : {pseudo_source} (-p/--pseudo-dir override)", f_out)
        else:
            disp_dirs = sorted(glob.glob(os.path.join(phonon_dir, "disp-*")))
            if not disp_dirs:
                print_dual(color_text(
                    f"[CRITICAL ERROR] No disp-*/ folders found in '{phonon_dir}' to reuse "
                    "pseudopotentials from, and no -p/--pseudo-dir override was given. Pass "
                    "-p/--pseudo-dir explicitly.", 'red'), f_out)
                sys.exit(1)
            pseudo_source = disp_dirs[0]
            print_dual(f"Source            : {pseudo_source} (reused from Stage 1's disp-* folders)", f_out)
        pseudos, missing = get_required_pseudos(unique_elements, pseudo_source)
        if missing:
            print_dual(color_text(
                f"[CRITICAL ERROR] Missing pseudopotential(s) for: {', '.join(sorted(missing))} "
                f"in '{pseudo_source}'.", 'red'), f_out)
            print_dual(color_text(
                "If this directory came from a real stb-ir (Stage 1) run this shouldn't "
                "happen -- Stage 1 requires every element's pseudopotential before it writes "
                "any folder. Re-run stb-ir, or pass -p/--pseudo-dir pointing at a folder "
                "that has all of the elements listed above.", 'yellow'), f_out)
            sys.exit(1)
        print_dual(f"Found all required : {', '.join(os.path.basename(p) for p in pseudos)}", f_out)

        apply_symmetry_skip = (args.use_symmetry and args.modes is None and not symmetry_error)
        selected = []
        n_actually_skipped = 0
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            if args.modes is not None and k not in args.modes:
                continue
            if args.freq_min is not None and freq < args.freq_min:
                continue
            if args.freq_max is not None and freq > args.freq_max:
                continue
            if apply_symmetry_skip:
                ms = band_to_symmetry.get(int(band_idx))
                if ms is not None and not ms.is_ir_active:
                    n_actually_skipped += 1
                    continue
            selected.append((k, freq, int(band_idx)))

        if apply_symmetry_skip and n_actually_skipped:
            print_dual(f"\n{color_text('[Symmetry]', 'cyan')} Skipped {n_actually_skipped} "
                        "mode(s) confirmed symmetry-forbidden from being IR-active (see "
                        "[1b] SYMMETRY ANALYSIS above).", f_out)

        if not selected:
            print_dual(color_text(
                "\n[ERROR] No modes selected after applying --modes/--freq-min/--freq-max"
                + ("/--use-symmetry" if apply_symmetry_skip else "")
                + " -- nothing to do.", 'red'), f_out)
            sys.exit(1)

        # --skip-degenerate (non-bulk path only -- see the flag's own
        # --help): within a degenerate group of selected modes, keep only
        # the first-encountered band as the group's representative. Same
        # rotational-invariant argument as stb-ramanModes: |dmu/dQ|^2 is
        # unchanged by any orthogonal transformation relating two
        # degenerate partners, so the representative's intensity applies
        # exactly to every other partner in its group.
        representative_of = {k: k for k, _, _ in selected}
        derived_groups = {}  # representative_k -> [derived_k, ...]
        skip_degenerate_applied = False
        if args.skip_degenerate:
            if is_bulk:
                print_dual(color_text(
                    "[NOTE] --skip-degenerate has no effect on the bulk path -- it already "
                    "writes exactly 1 folder regardless of mode count, so there is no folder "
                    "cost to save; every selected mode still gets its own reported IR "
                    "intensity.", 'yellow'), f_out)
            elif symmetry_error:
                print_dual(color_text(
                    f"[WARNING] --skip-degenerate requested but symmetry classification is "
                    f"unavailable ({symmetry_error}) -- running every mode individually, same "
                    "as without the flag.", 'yellow'), f_out)
            else:
                seen_groups = {}  # id(IRModeSymmetry) -> representative k
                for k, freq, band_idx in selected:
                    ms = band_to_symmetry.get(band_idx)
                    if ms is None or len(ms.band_indices) <= 1:
                        continue
                    group_key = id(ms)
                    if group_key not in seen_groups:
                        seen_groups[group_key] = k
                    else:
                        rep_k = seen_groups[group_key]
                        representative_of[k] = rep_k
                        derived_groups.setdefault(rep_k, []).append(k)
                        skip_degenerate_applied = True

        print_section('[3] IR DISPLACEMENT FOLDERS', f_out)
        print_dual(f"Selected modes    : {len(selected)}", f_out)
        if skip_degenerate_applied:
            n_derived_total = sum(len(v) for v in derived_groups.values())
            print_dual(f"Degenerate groups : {n_derived_total} mode(s) skipped -- their IR "
                        "intensity will be reused from a representative partner (see below).", f_out)
            for rep_k in sorted(derived_groups):
                rep_band = next(band_idx for k, _, band_idx in selected if k == rep_k)
                ms = band_to_symmetry.get(rep_band)
                label_str = f" ({ms.label})" if ms is not None and ms.label else ""
                derived_list = ", ".join(str(x) for x in sorted(derived_groups[rep_k]))
                print_dual(f"  mode {rep_k:3d}{label_str} -> also covers mode(s) {derived_list} "
                            "by symmetry (no folders written for them)", f_out)

        if args.export_animations:
            anim_root = born_charge_root if is_bulk else dipole_root
            os.makedirs(anim_root, exist_ok=True)
            print_dual(f"Mode animations   : {len(selected)} .axsf file(s), "
                        f"{args.animation_frames} frames each.", f_out)

        report_rows = []  # (label, mode_index, band_index, frequency, path, sign, dir, derived_from)
        structure_filename = "structure.fdf"

        if is_bulk:
            n_representative = sum(1 for k, _, _ in selected if representative_of[k] == k)
            print_dual(f"Folders to write  : 1 equilibrium Born-effective-charge SIESTA run "
                        f"(covers {n_representative} representative mode(s) out of "
                        f"{len(selected)} selected).", f_out)

            equilibrium_dir = os.path.join(born_charge_root, "equilibrium")
            calc_text = force_born_charge_run(calc_template, args.fc_displ, args.polarization_grid)
            write_ir_folder(equilibrium_dir, phonon.primitive, structure_filename, calc_text, pseudos)
            print_dual(f"  {color_text('[OK]', 'green')} {equilibrium_dir}", f_out)

            os.makedirs(born_charge_root, exist_ok=True)
            for k, freq, band_idx in selected:
                if args.export_animations:
                    animation_path = os.path.join(born_charge_root, f"mode_{k:02d}_animation.axsf")
                    write_mode_animation(phonon, band_idx, args.displacement, internal_to_angstrom,
                                          animation_path, n_frames=args.animation_frames)
                    print_dual(f"  {color_text('[OK]', 'green')} {animation_path}", f_out)

                eigendisp = mode_eigendisplacement(phonon, band_idx, internal_to_angstrom)
                sidecar_path = os.path.join(born_charge_root, f"mode_{k:02d}_eigendisplacement.json")
                with open(sidecar_path, "w") as f:
                    json.dump({"eigendisplacement": eigendisp.tolist()}, f)

                label = f"mode_{k:02d}_born_charge"
                report_rows.append((label, k, band_idx, freq, "BULK", "-", equilibrium_dir, None))

            n_real_folders = 1
        else:
            n_folders = sum(1 for k, _, _ in selected if representative_of[k] == k) * 2
            print_dual(f"Folders to write  : {n_folders} independent single-point SIESTA runs "
                        f"(2 signs x per-mode dipole displacement, across "
                        f"{n_folders // 2} representative mode(s) out of {len(selected)} "
                        "selected).", f_out)

            for k, freq, band_idx in selected:
                if args.export_animations:
                    animation_path = os.path.join(dipole_root, f"mode_{k:02d}_animation.axsf")
                    write_mode_animation(phonon, band_idx, args.displacement, internal_to_angstrom,
                                          animation_path, n_frames=args.animation_frames)
                    print_dual(f"  {color_text('[OK]', 'green')} {animation_path}", f_out)

                if representative_of[k] != k:
                    rep_k = representative_of[k]
                    label = f"mode_{k:02d}_derived"
                    report_rows.append((label, k, band_idx, freq, "NONBULK", "DERIVED", "-", rep_k))
                    continue

                for sign_name, sign_val in (("plus", 1.0), ("minus", -1.0)):
                    displaced = displace_along_mode(
                        phonon, band_idx, args.displacement, internal_to_angstrom, sign=sign_val)
                    label = f"mode_{k:02d}_{sign_name}"
                    mode_dir = os.path.join(dipole_root, label)
                    calc_text = force_single_point(calc_template)
                    write_ir_folder(mode_dir, displaced, structure_filename, calc_text, pseudos)
                    report_rows.append((label, k, band_idx, freq, "NONBULK", sign_name, mode_dir, None))
                    print_dual(f"  {color_text('[OK]', 'green')} {mode_dir}", f_out)

            n_real_folders = sum(1 for row in report_rows if row[5] != "DERIVED")

        print_section('[4] SUMMARY & NEXT STEPS', f_out)
        print_dual(f"{n_real_folders} folder(s) written under '{output_root}'.", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        if is_bulk:
            print_dual(f"  1. Run SIESTA in '{born_charge_root}/equilibrium/'.", f_out)
        else:
            print_dual(f"  1. Run SIESTA in every '{dipole_root}/mode_*/' folder.", f_out)
        print_dual(f"  2. Once they're done, run: stb-irAnalysis --directory {output_root}", f_out)

        f_out.write("\n# MODE_TABLE -- parsed by stb-irAnalysis, do not reorder the "
                     "first 7 columns\n")
        f_out.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
                     f"{'path':<8}{'sign':<8}{'dir':<24}{'derived_from'}\n")
        for label, k, band_idx, freq, path_kind, sign_name, mode_dir, derived_from in report_rows:
            derived_str = "-" if derived_from is None else str(derived_from)
            f_out.write(f"{label:<24}{k:<12}{band_idx:<12}{freq:<16.6f}{path_kind:<8}{sign_name:<8}"
                         f"{mode_dir}  {derived_str}\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Displacement folder(s) ready for Stage 3 (stb-irAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
