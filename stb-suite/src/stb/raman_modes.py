#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
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
from stb.core.cli import color_text, show_intro
from stb.core.calc_directives import force_single_point, build_optical_block
from stb.core.pseudopotentials import get_required_pseudos, resolve_pseudo_source
from stb.core import kspace
from stb.core.phonon_workflow import (
    detect_system_label, load_phonon_with_force_constants, get_gamma_modes, displace_along_mode,
)
from stb.core.raman_symmetry import classify_modes, tensor_form

# Default THz tolerance for the extra 0D trivial-mode (free rotation)
# search in get_gamma_modes. Deliberately conservative (LOW), not just
# "generously above the rotational noise seen in testing": the two
# failure directions aren't symmetric. Too low just means a rotational
# mode doesn't get caught and falls back to today's behavior (treated as
# an ordinary mode, same as before this fix -- harmless). Too high can
# silently swallow a genuine low-frequency vibrational mode (skeletal
# torsions in floppy/heavy molecules are a known case that can sit well
# under 5 THz / 167 cm^-1) before the user ever sees it in Stage 2's mode
# list -- unlike the group-theory --use-symmetry skip, this has no
# "[symmetry-forbidden]" audit trail to catch it. Verified live with a
# relaxed H2O molecule via MACE-MP-0: rotational bands were ~-0.01 to
# 0.22 THz, real vibrational modes 44.5+ THz -- 2.0 THz keeps a ~10x
# margin above the observed rotational noise while staying well under
# realistic low-frequency vibrational modes. Exposed as --rotational-
# mode-tol since it's a real judgment call for unusual molecules, not a
# pure numerical-noise constant.
_DEFAULT_ROTATIONAL_MODE_TOL_THZ = 2.0

REPORT_FILE = "raman_stage2.txt"

# Diagonal Raman tensor (Rxx, Ryy, Rzz) -- the default/cheap scope, 6
# Optical calculations per mode (2 signs x 3 axes).
_AXES_DIAGONAL = [("x", (1.0, 0.0, 0.0)), ("y", (0.0, 1.0, 0.0)), ("z", (0.0, 0.0, 1.0))]

# Off-diagonal probe directions for --full-tensor: the (i+j)/sqrt(2)
# face-diagonals. A scalar Optical.Vector measurement along direction n
# gives n^T R n; for n=(i+j)/sqrt(2) that's (Rii+Rjj)/2 + Rij, so combined
# with the already-measured diagonal components Rii/Rjj (from _AXES_DIAGONAL)
# stb-ramanAnalysis can recover Rij = R_(i+j)/sqrt(2) - (Rii+Rjj)/2. Doubles
# the per-mode folder count (12 instead of 6) -- opt-in via --full-tensor.
_INV_SQRT2 = 0.70710678118654752440
_AXES_OFFDIAG = [
    ("xy", (_INV_SQRT2, _INV_SQRT2, 0.0)),
    ("xz", (_INV_SQRT2, 0.0, _INV_SQRT2)),
    ("yz", (0.0, _INV_SQRT2, _INV_SQRT2)),
]
_SIGNS = [("plus", 1.0), ("minus", -1.0)]

# name -> vector lookup covering all 6 possible probe directions, for
# --reduce-components (core.raman_symmetry.tensor_form picks a probe_axis
# name from this same set; this is how the folder-writing loop turns that
# name back into the actual Optical.Vector).
_ALL_AXIS_VECTORS = dict(_AXES_DIAGONAL + _AXES_OFFDIAG)


def print_dual(text, file_handle=None):
    """Prints to stdout with color, writes to file without color. Same
    duplicated-per-tool helper as phonons_create.py/phonons_pos.py."""
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


def write_optical_folder(out_dir, displaced_atoms, structure_filename, calc_text, pseudos):
    os.makedirs(out_dir, exist_ok=True)
    write_siesta(os.path.join(out_dir, os.path.basename(structure_filename)), displaced_atoms)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    for pseudo_path in pseudos:
        shutil.copy(pseudo_path, os.path.join(out_dir, os.path.basename(pseudo_path)))


def _phonopy_atoms_to_ase(patoms, internal_to_angstrom):
    """PhonopyAtoms (internal units -- bohr for a SIESTA-sourced phonon
    object, already Angstrom for an ML-sourced one, same convention as
    everywhere else in this workflow) -> ase.Atoms in real Angstrom, ready
    for ase.io.write (xsf/axsf expects real physical units, no internal-
    unit ambiguity of its own).
    """
    return Atoms(symbols=patoms.symbols,
                 positions=np.array(patoms.positions) * internal_to_angstrom,
                 cell=np.array(patoms.cell) * internal_to_angstrom,
                 pbc=True)


def write_mode_animation(phonon, band_index, amplitude_ang, internal_to_angstrom, out_path, n_frames=20):
    """Writes a looping animation of one Gamma-point mode's eigendisplacement
    (a smooth 0 -> +A -> 0 -> -A -> 0 sweep, `n_frames` frames) as an
    animated XSF (.axsf) file -- readable directly by XCrySDen/VESTA to
    visually inspect which atoms move in this mode, before interpreting
    the spectrum. Reuses displace_along_mode (already used to build the
    real Optical-calculation folders) at each frame's signed amplitude --
    passing a NEGATIVE amplitude_ang directly (rather than using its
    `sign` parameter) works identically, since displace_along_mode only
    ever multiplies the two together internally.

    ase.io.write(..., format='xsf') accepts a list of Atoms directly and
    produces the ANIMSTEPS-header animated format on its own -- verified
    during planning (ase.io.formats.ioformats['xsf'].single is False) --
    so this doesn't hand-roll the AXSF format.
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
        "displacements, identifies the Gamma-point vibrational modes, and generates the +/-delta "
        "Optical-calculation displacement folders needed for the Raman tensor.", 'bold')}
For each selected mode, writes 6 folders by default (2 signs x 3
Optical.Vector axes: x, y, z -- a diagonal-only Raman tensor Rxx/Ryy/Rzz)
under <directory>/optical_disp/, or 12 with --full-tensor (adds the
(x+y)/sqrt(2), (x+z)/sqrt(2), (y+z)/sqrt(2) face-diagonal directions needed
to also recover Rxy/Rxz/Ryz). stb-ramanAnalysis auto-detects which scope
was used from the folders actually present. Doesn't run SIESTA itself --
run each folder's calculation yourself, then use stb-ramanAnalysis.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory raman_study --calc calc_optical.fdf\n"
               "  %(prog)s --directory raman_study --calc calc_optical.fdf --modes 1 2 3\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="raman_study",
                        help="Root directory written by stb-raman (Stage 1) -- must contain "
                             "'phonon_disp/' (default: raman_study).")
    parser.add_argument("-l", "--label", type=str, default=None,
                        help="SystemLabel used in Stage 1's calc.fdf (default: auto-detected).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default=None, metavar="DIR",
                        help="Pseudopotentials source: a bundled bank or a folder path. Default: "
                             "reuse whatever Stage 1 already resolved and copied into its "
                             "phonon_disp/disp-*/ folders (the normal, self-contained case -- "
                             "no need to point at a source a second time). Only needed if that "
                             "folder is missing pseudopotentials for some reason (e.g. hand-built "
                             "phonon_disp/ instead of a real stb-raman run).")
    parser.add_argument("-c", "--calc", type=str, required=True,
                        help="calc.fdf template for the Optical calculations -- can be the same "
                             "file used in Stage 1 or a different one (e.g. a denser k-grid, "
                             "since these are now single-point evaluations on the small unit "
                             "cell, not the phonon supercell).")
    parser.add_argument("--modes", type=int, nargs='+', default=None, metavar="N",
                        help="1-based mode index/indices to process (in ascending-frequency "
                             "order, acoustic modes already excluded -- see Stage 2's own "
                             "[1] PHONON MODES table for the numbering). Default: every "
                             "non-acoustic mode.")
    parser.add_argument("--freq-min", type=float, default=None,
                        help="Skip modes below this frequency (THz).")
    parser.add_argument("--freq-max", type=float, default=None,
                        help="Skip modes above this frequency (THz).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum gap (Ang) along an axis to consider it vacuum-padded, "
                             "same convention as stb-raman (default: 10.0). Only used to detect "
                             "a genuinely isolated (0D) structure, which has 3 extra trivial "
                             "(free-rotation) modes at Gamma beyond the usual 3 translations -- "
                             "irrelevant for 3D/2D/1D structures.")
    parser.add_argument("--rotational-mode-tol", type=float, default=_DEFAULT_ROTATIONAL_MODE_TOL_THZ,
                        help="0D (molecule) only: a mode within the first 3 bands after the "
                             "translations is treated as free rotation (trivial, excluded) if "
                             "|frequency| is below this (THz) (default: "
                             f"{_DEFAULT_ROTATIONAL_MODE_TOL_THZ}). Lower this if your molecule "
                             "has genuine low-frequency skeletal/torsional vibrational modes you "
                             "want to make sure aren't mistaken for rotation.")
    parser.add_argument("--displacement", type=float, default=0.02,
                        help="Finite-difference displacement (Ang) along each mode's "
                             "eigendisplacement for the Raman-tensor derivative (default: 0.02).")
    parser.add_argument("--full-tensor", action="store_true",
                        help="Compute the FULL symmetric Raman tensor (Rxx, Ryy, Rzz, Rxy, Rxz, "
                             "Ryz) instead of just the diagonal -- adds 3 more Optical.Vector "
                             "directions per sign (the face-diagonals), doubling the folder count "
                             "per mode (12 instead of 6). Off-diagonal components are recovered "
                             "by stb-ramanAnalysis from these mixed-direction measurements.")
    parser.add_argument("--use-symmetry", action="store_true",
                        help="Skip modes that group theory guarantees are Raman-INACTIVE (their "
                             "irreducible representation at Gamma isn't contained in the "
                             "quadratic-function representation x^2/y^2/z^2/xy/xz/yz -- their "
                             "Raman tensor is exactly zero, not just numerically small). Works "
                             "for any of the 32 point groups, not just centrosymmetric ones. "
                             "Only applies to the default mode selection (every non-acoustic "
                             "mode) -- an explicit --modes list is never second-guessed. Needs a "
                             "PRIMITIVE cell (see stb-unitcell --mode primitive); falls back to "
                             "running every mode, with a warning, if that's not available.")
    parser.add_argument("--reduce-components", action="store_true",
                        help="For a non-degenerate mode whose Raman tensor SHAPE is fixed by "
                             "symmetry up to one overall scale (a single-band irrep with "
                             "multiplicity exactly 1 in the quadratic-function representation), "
                             "probe only ONE Optical.Vector direction instead of all of them -- "
                             "stb-ramanAnalysis reconstructs the rest algebraically. Falls back "
                             "to full probing for degenerate modes or any mode whose tensor shape "
                             "isn't uniquely fixed (needs a PRIMITIVE cell, same as "
                             "--use-symmetry; degrades to no reduction, never a crash, if that's "
                             "not available).")
    parser.add_argument("--skip-degenerate", action="store_true",
                        help="For a degenerate group of modes (2 or more bands sharing one "
                             "irrep, e.g. an E or T representation), write Optical folders for "
                             "only ONE representative band and skip the rest entirely -- Raman "
                             "activity and the depolarization ratio are rotational invariants of "
                             "the tensor, so whatever orthogonal transformation relates two "
                             "degenerate partners (guaranteed to exist by group theory, even "
                             "though phonopy's returned eigenvectors don't tell you which one) "
                             "leaves both invariants unchanged -- the representative's computed "
                             "values apply exactly to every other partner in its group. "
                             "stb-ramanAnalysis reuses them automatically. Needs a PRIMITIVE "
                             "cell (same as --use-symmetry/--reduce-components; degrades to "
                             "running every mode, with a warning, if that's not available).")
    parser.add_argument("--optical-mesh", type=int, nargs=3, default=[10, 10, 10],
                        help="Optical.Mesh k-grid for the dielectric-function calculation "
                             "(default: 10 10 10).")
    parser.add_argument("--optical-nbands", type=int, default=None,
                        help="Optical.NumberOfBands (default: SIESTA's own default -- all "
                             "bands within the basis).")
    parser.add_argument("--optical-broaden", type=float, default=0.2,
                        help="Optical.Broaden in eV (default: 0.2).")
    parser.add_argument("--export-animations", action="store_true",
                        help="Write a looping animation of each selected mode's eigendisplacement "
                             "(a smooth 0/+A/0/-A/0 sweep) as optical_disp/mode_XX_animation.axsf "
                             "-- open in XCrySDen/VESTA to see which atoms move before interpreting "
                             "the spectrum. Same amplitude as --displacement.")
    parser.add_argument("--animation-frames", type=int, default=20,
                        help="Frames per mode animation, with --export-animations (default: 20).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-ramanModes {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("RAMAN SPECTRUM WORKFLOW -- STAGE 2: MODES & OPTICAL DISPLACEMENTS", 'bold'))
    print("-" * 60)

    phonon_dir = os.path.join(args.directory, "phonon_disp")
    if not os.path.isdir(phonon_dir):
        print(color_text(f"[ERROR] '{phonon_dir}' not found -- run stb-raman (Stage 1) first.", 'red'))
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
        optical_calc_template = f.read()

    if args.pseudo_dir is not None:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    axes = _AXES_DIAGONAL + _AXES_OFFDIAG if args.full_tensor else _AXES_DIAGONAL

    output_root = args.directory
    optical_root = os.path.join(output_root, "optical_disp")
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== RAMAN STAGE 2 REPORT (MODES & OPTICAL DISPLACEMENTS) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Directory         : {output_root}", f_out)
        print_dual(f"SystemLabel       : {'N/A (ML-computed force constants)' if has_embedded_fc else f'{system_label} ({label_source})'}", f_out)
        print_dual(f"Optical calc      : {args.calc}", f_out)
        print_dual(f"Displacement      : {args.displacement} Ang", f_out)
        print_dual(f"Optical mesh      : {args.optical_mesh[0]} x {args.optical_mesh[1]} x {args.optical_mesh[2]}", f_out)
        print_dual(f"Optical broaden   : {args.optical_broaden} eV", f_out)
        print_dual(f"Tensor scope      : {'full symmetric tensor' if args.full_tensor else 'diagonal only'}", f_out)

        print_dual(f"\n{color_text('[1] PHONON MODES AT GAMMA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        phonon, internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            phonon_dir, system_label, has_embedded_fc, f_out)
        try:
            unique_elements = list(set(phonon.primitive.symbols))
            lattice_ang = np.array(phonon.primitive.cell) * internal_to_angstrom
            frac_coords = np.array(phonon.primitive.scaled_positions)
            vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
            is_0d = all(vacuum_axes)
            frequencies, mode_band_indices, n_extra_trivial = get_gamma_modes(
                phonon, exclude_acoustic=True,
                extra_trivial_tol_thz=args.rotational_mode_tol if is_0d else None)
            mode_symmetries, point_group, symmetry_error = classify_modes(phonon)
        finally:
            os.chdir(original_dir)

        band_to_symmetry = {}
        for ms in mode_symmetries:
            for b in ms.band_indices:
                band_to_symmetry[b] = ms

        print_dual(f"Dimensionality    : {kspace.dimensionality_label(vacuum_axes)}", f_out)
        if n_extra_trivial:
            print_dual(color_text(
                f"[NOTE] 0D structure detected -- excluded {n_extra_trivial} extra trivial "
                "mode(s) (free rotation) beyond the usual 3 translations. These aren't real "
                "vibrations and have no Raman activity to speak of.", 'yellow'), f_out)

        n_imaginary = int((frequencies < 0).sum())
        print_dual(f"Non-acoustic Gamma modes : {len(frequencies)}", f_out)
        if n_imaginary:
            print_dual(color_text(
                f"[WARNING] {n_imaginary} mode(s) have imaginary (negative) frequency -- the "
                "structure/supercell may not be at a real energy minimum. Their Raman tensor "
                "is not physically meaningful.", 'yellow'), f_out)
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            flag = color_text(" [IMAGINARY]", 'red') if freq < 0 else ""
            ms = band_to_symmetry.get(int(band_idx))
            sym_note = ""
            if args.use_symmetry and ms is not None:
                label_str = f" ({ms.label})" if ms.label else ""
                sym_note = (color_text(f"{label_str} [Raman-active]", 'green') if ms.is_raman_active
                            else color_text(f"{label_str} [symmetry-forbidden]", 'yellow'))
            print_dual(f"  mode {k:3d} (band {int(band_idx):3d}) : {freq:10.4f} THz{flag} {sym_note}", f_out)

        if args.use_symmetry:
            print_dual(f"\n{color_text('[1b] SYMMETRY ANALYSIS', 'magenta')}", f_out)
            print_dual("-" * 60, f_out)
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
                                   and not band_to_symmetry[int(band_idx)].is_raman_active)
                n_unknown = sum(1 for band_idx in mode_band_indices
                                 if band_to_symmetry.get(int(band_idx)) is None)
                print_dual(f"Symmetry-forbidden (Raman-inactive) : {n_forbidden}/{len(mode_band_indices)}", f_out)
                if n_unknown:
                    print_dual(color_text(
                        f"[NOTE] {n_unknown} mode(s) could not be classified (label matching "
                        "failed) -- kept, never auto-skipped when uncertain.", 'yellow'), f_out)
                if args.modes is not None:
                    print_dual(color_text(
                        "--modes was given explicitly -- symmetry filtering is informational "
                        "only here, no mode is auto-skipped.", 'yellow'), f_out)

        print_dual(f"\n{color_text('[2] PSEUDOPOTENTIALS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
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
                "If this directory came from a real stb-raman (Stage 1) run this shouldn't "
                "happen -- Stage 1 requires every element's pseudopotential before it writes "
                "any folder. Re-run stb-raman, or pass -p/--pseudo-dir pointing at a folder "
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
                if ms is not None and not ms.is_raman_active:
                    n_actually_skipped += 1
                    continue
            selected.append((k, freq, int(band_idx)))

        if apply_symmetry_skip and n_actually_skipped:
            print_dual(f"\n{color_text('[Symmetry]', 'cyan')} Skipped {n_actually_skipped} "
                        "mode(s) confirmed symmetry-forbidden from being Raman-active (see "
                        "[1b] SYMMETRY ANALYSIS above).", f_out)

        if not selected:
            print_dual(color_text(
                "\n[ERROR] No modes selected after applying --modes/--freq-min/--freq-max"
                + ("/--use-symmetry" if apply_symmetry_skip else "")
                + " -- nothing to do.", 'red'), f_out)
            sys.exit(1)

        # --skip-degenerate: within a degenerate group of selected modes
        # (2+ bands sharing one irrep), keep only the first-encountered band
        # as the group's representative -- Raman activity and rho are
        # rotational invariants (trace + deviatoric norm of the tensor), so
        # they're identical for every partner in the group regardless of
        # which specific orthogonal transformation phonopy's arbitrary
        # eigenvector choice corresponds to (see --skip-degenerate's
        # --help). representative_of[k] == k for every non-derived mode.
        representative_of = {k: k for k, _, _ in selected}
        derived_groups = {}  # representative_k -> [derived_k, ...]
        skip_degenerate_applied = False
        if args.skip_degenerate:
            if symmetry_error:
                print_dual(color_text(
                    f"[WARNING] --skip-degenerate requested but symmetry classification is "
                    f"unavailable ({symmetry_error}) -- running every mode individually, same "
                    "as without the flag.", 'yellow'), f_out)
            else:
                seen_groups = {}  # id(ModeSymmetry) -> representative k
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

        # --reduce-components: for each selected mode whose Raman tensor
        # SHAPE is fixed by symmetry up to one overall scale (single-band,
        # multiplicity-1 in the quadratic representation -- see
        # core.raman_symmetry.tensor_form), probe only that one direction
        # instead of the full `axes` list; everything else falls back to
        # full probing (degenerate modes, multiplicity >= 2, or symmetry
        # unavailable -- tensor_form itself returns None in all of those,
        # never raises). Derived (--skip-degenerate) modes never get here --
        # they don't write folders at all.
        mode_axes = {}
        mode_reductions = {}
        for k, freq, band_idx in selected:
            if representative_of[k] != k:
                continue
            reduction = None
            if args.reduce_components:
                ms = band_to_symmetry.get(band_idx)
                is_degenerate = ms is not None and len(ms.band_indices) > 1
                if not is_degenerate:
                    reduction = tensor_form(phonon, band_idx, full_tensor=args.full_tensor)
            if reduction is not None:
                mode_axes[k] = [(reduction["probe_axis"], _ALL_AXIS_VECTORS[reduction["probe_axis"]])]
                mode_reductions[k] = reduction
            else:
                mode_axes[k] = axes

        print_dual(f"\n{color_text('[3] OPTICAL DISPLACEMENT FOLDERS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        n_folders = sum(len(mode_axes[k]) for k, _, _ in selected
                         if representative_of[k] == k) * len(_SIGNS)
        print_dual(f"Selected modes    : {len(selected)}", f_out)
        if skip_degenerate_applied:
            n_derived_total = sum(len(v) for v in derived_groups.values())
            print_dual(f"Degenerate groups : {n_derived_total} mode(s) skipped -- their Raman "
                        "tensor will be reused from a representative partner (see below).", f_out)
            for rep_k in sorted(derived_groups):
                rep_band = next(band_idx for k, _, band_idx in selected if k == rep_k)
                ms = band_to_symmetry.get(rep_band)
                label_str = f" ({ms.label})" if ms is not None and ms.label else ""
                derived_list = ", ".join(str(x) for x in sorted(derived_groups[rep_k]))
                print_dual(f"  mode {rep_k:3d}{label_str} -> also covers mode(s) {derived_list} "
                            "by symmetry (no folders written for them)", f_out)
        if mode_reductions:
            print_dual(f"Component reduction : {len(mode_reductions)}/{len(selected)} mode(s) "
                        "reduced to a single probe axis (tensor shape fixed by symmetry).", f_out)
            for k in sorted(mode_reductions):
                print_dual(f"  mode {k:3d} -> probe axis '{mode_reductions[k]['probe_axis']}' only "
                            "(other components reconstructed algebraically in Stage 3)", f_out)
        if mode_reductions or skip_degenerate_applied:
            n_representative = sum(1 for k, _, _ in selected if representative_of[k] == k)
            print_dual(f"Folders to write  : {n_folders} independent single-point SIESTA runs "
                        f"({len(_SIGNS)} signs x per-mode axis count, across {n_representative} "
                        f"representative mode(s) out of {len(selected)} selected).", f_out)
        else:
            print_dual(f"Folders to write  : {len(selected)} modes x {len(_SIGNS)} signs x "
                        f"{len(axes)} axes = {n_folders} independent single-point SIESTA runs.", f_out)

        if args.export_animations:
            os.makedirs(optical_root, exist_ok=True)
            print_dual(f"Mode animations   : {len(selected)} .axsf file(s), "
                        f"{args.animation_frames} frames each.", f_out)

        report_rows = []  # (label, mode_index, band_index, frequency, sign, axis, dir, derived_from)
        structure_filename = "structure.fdf"
        for k, freq, band_idx in selected:
            if args.export_animations:
                animation_path = os.path.join(optical_root, f"mode_{k:02d}_animation.axsf")
                write_mode_animation(phonon, band_idx, args.displacement, internal_to_angstrom,
                                      animation_path, n_frames=args.animation_frames)
                print_dual(f"  {color_text('[OK]', 'green')} {animation_path}", f_out)
            if representative_of[k] != k:
                rep_k = representative_of[k]
                label = f"mode_{k:02d}_derived"
                report_rows.append((label, k, band_idx, freq, "DERIVED", "-", "-", rep_k))
                continue
            for sign_name, sign_val in _SIGNS:
                displaced = displace_along_mode(
                    phonon, band_idx, args.displacement, internal_to_angstrom, sign=sign_val)
                for axis_name, axis_vec in mode_axes[k]:
                    label = f"mode_{k:02d}_{sign_name}_{axis_name}"
                    mode_dir = os.path.join(optical_root, label)
                    optical_block = build_optical_block(
                        args.optical_mesh, args.optical_broaden, axis_vec, args.optical_nbands)
                    calc_text = force_single_point(optical_calc_template) + "\n" + optical_block
                    write_optical_folder(mode_dir, displaced, structure_filename, calc_text, pseudos)
                    report_rows.append((label, k, band_idx, freq, sign_name, axis_name, mode_dir, None))
                    print_dual(f"  {color_text('[OK]', 'green')} {mode_dir}", f_out)
            if k in mode_reductions:
                t0 = mode_reductions[k]["T0"]
                sidecar_path = os.path.join(optical_root, f"mode_{k:02d}_tensor_form.json")
                with open(sidecar_path, "w") as f:
                    json.dump({"T0": t0.tolist(), "probe_axis": mode_reductions[k]["probe_axis"]}, f)

        n_real_folders = sum(1 for row in report_rows if row[4] != "DERIVED")
        print_dual(f"\n{color_text('[4] SUMMARY & NEXT STEPS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"{n_real_folders} folder(s) written under '{optical_root}'.", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual(f"  1. Run SIESTA in every '{optical_root}/mode_*/' folder.", f_out)
        print_dual(f"  2. Once they're done, run: stb-ramanAnalysis --directory {output_root}", f_out)

        f_out.write("\n# MODE_TABLE -- parsed by stb-ramanAnalysis, do not reorder the "
                     "first 6 columns\n")
        f_out.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
                     f"{'sign':<8}{'axis':<6}{'dir':<24}{'derived_from'}\n")
        for label, k, band_idx, freq, sign_name, axis_name, mode_dir, derived_from in report_rows:
            derived_str = "-" if derived_from is None else str(derived_from)
            f_out.write(f"{label:<24}{k:<12}{band_idx:<12}{freq:<16.6f}{sign_name:<8}{axis_name:<6}"
                         f"{mode_dir}  {derived_str}\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Optical displacement folders ready for Stage 3 (stb-ramanAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
