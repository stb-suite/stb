#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.4.1"  # [1b] SYMMETRY ANALYSIS's degenerate-groups listing no longer prints an empty
                    # "modes  (T1u)" line for a group made entirely of ACOUSTIC bands (not in
                    # band_to_k, e.g. NaCl's own T1u translational mode shares its irrep label with
                    # the real IR-active T1u optical mode) -- now built as (ms, mode_ids) pairs and
                    # filtered to drop any group with no displayed (non-acoustic) member.
                    # Also: new [LIMITATION] note after [1]'s per-mode frequency listing, for the
                    # BULK (3D) path only -- these Gamma frequencies have no non-analytic (LO-TO)
                    # correction applied, so for a genuinely POLAR bulk crystal they can land
                    # between the true TO and LO values rather than matching either. Root-caused
                    # (not fixed) on NaCl: T1u computed at 185.6 cm^-1 vs. ~164/~264 cm^-1
                    # (TO/LO, Raunio et al. 1969), landing almost exactly midway -- the textbook
                    # symptom of a finite-supercell frozen-phonon calculation missing the Gonze &
                    # Lee (1997) dipole-dipole subtraction, which needs eps_inf (not currently
                    # computed by any stb tool) alongside the Z* this stage already computes.
                    # Documented as a known limitation, not fixed -- a real fix needs a new
                    # eps_inf-producing SIESTA Optical-module run plus a dynamical-matrix-level
                    # change in core/phonon_workflow.py, planned for a future version.
                    # (previously 1.4.0: --symprec (default 0.01, pymatgen's own default) threaded through to
                    # Phonopy itself, same fix as stb-ramanModes (VERSION 1.2.0): Phonopy's own
                    # raw default (1e-5) is far tighter than any DFT relaxation's real numerical
                    # noise floor and can silently misdetect the true point group -- this is the
                    # actual tolerance driving [1b] SYMMETRY ANALYSIS/--use-symmetry, so a
                    # too-tight symprec here can artificially split a truly degenerate IR-active
                    # mode into separate near-identical frequencies and mislabel a symmetry-silent
                    # mode as IR-active. [1b] SYMMETRY ANALYSIS (space group, point group,
                    # symmetry op count, per-mode Mulliken label/activity) is now always printed,
                    # not just with --use-symmetry, along with an explicit [CHECK] warning to
                    # verify the reported group against the crystal's known symmetry and loosen
                    # --symprec if it looks wrong -- misdetection here is silent, never an error.
                    # Also: every candidate mode being symmetry-forbidden (e.g. a centrosymmetric
                    # crystal with no IR-active Gamma modes at all) now finishes as a [WARNING]
                    # with exit code 0, not a fatal [ERROR]/exit(1) -- it's a valid physical
                    # result, not a misconfiguration. Still a hard [ERROR]/exit(1) when
                    # --modes/--freq-min/--freq-max (not symmetry) are what emptied the selection.
                    # New HYBRID path for exactly-one-vacuum-axis (2D slab) structures. A naive
                    # real-space dipole moment (the old "non-bulk" path, still used as-is for
                    # 0D/1D) is only physically valid along a genuinely NON-periodic direction --
                    # confirmed live on monolayer h-BN: SIESTA itself printed an EXACT ZERO
                    # dipole change for the in-plane E' mode's +/-delta displacement (not a
                    # parsing bug), silently reporting zero IR intensity for a mode the
                    # literature puts as the DOMINANT IR peak (~1420 cm^-1). The two in-plane
                    # (periodic) axes of a slab need the same Born-effective-charge/Berry-phase
                    # machinery the "bulk" path already uses; the vacuum axis is still fine via
                    # the old dipole-difference method (verified: the out-of-plane A2'' mode's
                    # dmu_z/dQ came out correct, matching the literature TO/LO range). HYBRID
                    # writes BOTH a shared born_charge_disp/equilibrium/ folder (PolarizationGrids
                    # built vacuum-aware via build_slab_polarization_grid -- the vacuum axis's own
                    # row zeroed, since SIESTA skips the polarization calculation entirely for any
                    # zero grid-count direction) AND per-mode dipole_disp/mode_XX_plus|minus/
                    # folders; stb-irAnalysis combines one component from each. Empirically
                    # validated before writing this code: a manual BornCharge T run on the relaxed
                    # h-BN structure with this exact grid shape gave a sane, sign-correct,
                    # sum-rule-satisfying in-plane Z* (B ~+2.76/+2.78, N ~-2.76/-2.74 |e|, in-plane
                    # sum ~0.01-0.04 |e|, out-of-plane row correctly ~0). Deliberately scoped to 2D
                    # only (exactly 1 vacuum axis) -- 1D (2 vacuum axes) stays on the old non-bulk
                    # path for now, same limitation, no test case yet.
                    # Also: [1] now tags each mode with a Cartesian-polarization character label
                    # (mode_character_label, e.g. "98% z" or "51% x, 47% y", from its
                    # eigendisplacement -- purely descriptive, compare against the Vacuum
                    # axes/Lattice vectors already printed to interpret in-plane vs. out-of-plane
                    # for a given structure), and [1b] always lists degenerate mode groups by
                    # symmetry (band indices sharing one irrep -- Phonopy's own eigenvector basis
                    # choice within such a group is arbitrary, informational regardless of
                    # --use-symmetry/--skip-degenerate).)

import os
import sys
import glob
import json
import shutil
import argparse
from datetime import datetime
import numpy as np
import yaml
from phonopy.interface.siesta import write_siesta
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.calc_directives import force_single_point, force_born_charge_run
from stb.core.pseudopotentials import get_required_pseudos, resolve_pseudo_source
from stb.core import kspace
from stb.core.phonon_workflow import (
    detect_system_label, load_phonon_with_force_constants, get_gamma_modes, displace_along_mode,
    mode_eigendisplacement, build_mode_animation_frames, write_mode_animation,
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
# In-plane diagonal/cross-term Berry-phase k-point counts for the HYBRID
# (2D slab) path's PolarizationGrids -- same magnitude as
# _DEFAULT_POLARIZATION_GRID's own in-plane entries, empirically validated
# live on monolayer h-BN (see VERSION comment above) before this code was
# written.
_SLAB_IN_PLANE_GRID = 10
_SLAB_CROSS_TERM_GRID = 4


def build_slab_polarization_grid(vacuum_axes, in_plane=_SLAB_IN_PLANE_GRID,
                                  cross_term=_SLAB_CROSS_TERM_GRID):
    """Vacuum-aware %block PolarizationGrids for the HYBRID (2D slab) path:
    a row-major 3x3 grid (flattened to 9 ints, same shape
    force_born_charge_run already expects) where any entry touching the
    vacuum-padded lattice vector is forced to 1 (mirrors the slab SCF
    k-grid's own Nx x Ny x 1 convention -- a large real-space vacuum gives a
    tiny reciprocal extent), EXCEPT the vacuum axis's own diagonal entry,
    forced to exactly 0 -- SIESTA's documented behavior is to skip the
    polarization calculation entirely for any direction whose grid count is
    0, which is exactly what's wanted here: the out-of-plane component is
    physically ill-defined via Berry phase for a slab's non-periodic axis
    (the existing dipole-difference path already covers it correctly, see
    ir_analysis.py's HYBRID combination). The two genuinely periodic
    in-plane axes get `in_plane`/`cross_term`, matching
    _DEFAULT_POLARIZATION_GRID's own in-plane values.

    `vacuum_axes` must have exactly one True entry (a 2D slab) -- callers
    only invoke this for the HYBRID path, which is gated on exactly that.
    """
    grid = [[0, 0, 0] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            touches_vacuum = vacuum_axes[i] or vacuum_axes[j]
            if i == j:
                grid[i][j] = 0 if vacuum_axes[i] else in_plane
            else:
                grid[i][j] = 1 if touches_vacuum else cross_term
    return [v for row in grid for v in row]


def mode_character_label(eigendisp):
    """Short descriptive tag for a mode's dominant Cartesian polarization,
    from its (n_atoms, 3) eigendisplacement pattern -- e.g. "98% z" for a
    purely out-of-plane-polarized mode, "51% x, 47% y" for an in-plane one.
    Purely descriptive of the raw Cartesian weights (sum of squared
    displacement per axis, normalized to fractions of the total) -- makes
    no assumption about which Cartesian axis is this particular
    structure's own vacuum/surface-normal direction, so compare against
    the "Lattice vectors"/"Vacuum axes" already printed above to interpret
    it as in-plane vs. out-of-plane for THIS structure. Lists axes by
    descending weight until >=97% of the total is accounted for (never
    more than all 3), so a clearly one-axis-dominated mode gets a single
    short entry instead of three near-zero trailing ones.
    """
    weights = np.sum(np.asarray(eigendisp, dtype=float) ** 2, axis=0)
    total = float(np.sum(weights))
    if total <= 0:
        return ""
    fractions = weights / total
    order = np.argsort(fractions)[::-1]
    axis_names = "xyz"
    parts, cumulative = [], 0.0
    for idx in order:
        if cumulative >= 0.97 and parts:
            break
        parts.append(f"{fractions[idx] * 100:.0f}% {axis_names[idx]}")
        cumulative += fractions[idx]
    return ", ".join(parts)


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




def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 3: builds FORCE_SETS from Stage 1's phonon "
        "displacements, identifies the Gamma-point vibrational modes, and generates whichever "
        "follow-up SIESTA folder(s) the IR dipole-derivative calculation needs.", 'bold')}
The path is AUTO-SELECTED by structure dimensionality (core.kspace.detect_vacuum_axes on the
primitive cell), no flag needed:

  - Non-bulk (0D molecule, or 1D wire -- zero or two vacuum-padded axes): SIESTA prints the
    total dipole moment automatically for a non-periodic direction, no special fdf block
    needed. Writes exactly 2 folders per selected mode (dipole_disp/mode_XX_plus/,
    mode_XX_minus/) -- a simple single-point SCF, no per-axis probing (unlike
    stb-ramanModes' Optical.Vector sweep): the dipole derivative dmu/dQ already comes out as
    a full 3-vector from one +/-delta pair. NOTE: for a 1D wire this is only valid along the
    two genuinely non-periodic transverse directions -- the along-axis component is subject
    to the same periodic-direction caveat the HYBRID path below exists to fix for 2D; 1D isn't
    covered by that fix yet.

  - Bulk (3D periodic, zero vacuum axes): a naive dipole moment is gauge-ambiguous for a fully
    periodic system. Writes exactly ONE equilibrium folder (born_charge_disp/equilibrium/,
    the undisplaced structure) using SIESTA's own native Born-effective-charge automation
    (MD.TypeOfRun FC + BornCharge T + %block PolarizationGrids) -- independent of how many
    modes are selected. Also writes a small eigendisplacement JSON sidecar per selected mode
    (born_charge_disp/mode_XX_eigendisplacement.json) so stb-irAnalysis can combine it with the
    Born-charge tensors without reloading Phonopy itself.

  - Hybrid (2D slab, exactly one vacuum-padded axis): a naive dipole moment is only valid
    along the slab's non-periodic (vacuum) axis -- the two in-plane axes are genuinely
    periodic, same gauge-ambiguity problem the bulk path exists to solve. Writes BOTH the
    bulk path's single shared born_charge_disp/equilibrium/ folder (PolarizationGrids built
    vacuum-aware -- the vacuum axis's own row zeroed, so SIESTA skips its ill-defined
    out-of-plane polarization component entirely) AND the non-bulk path's per-mode
    dipole_disp/mode_XX_plus/mode_XX_minus/ folders. stb-irAnalysis takes the in-plane
    dmu/dQ components from the Born-charge run and the vacuum-axis component from the
    dipole-difference pair, combining them into one physically valid 3-vector per mode.

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
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry-detection tolerance (Ang) Phonopy uses internally when "
                             "loading the force constants -- drives both the reported point "
                             "group and, when --use-symmetry is on, the IR-active/inactive "
                             "classification and degenerate-mode grouping (default: 0.01, "
                             "pymatgen's own default -- matches the rest of the suite, and "
                             "deliberately NOT Phonopy's own raw default of 1e-5, which is far "
                             "too tight for a real DFT-relaxed structure and can misdetect the "
                             "true point group: a mode that should be exactly degenerate can "
                             "come out as two separate near-identical frequencies, and a "
                             "symmetry-silent mode can be misclassified as IR-active). Loosen "
                             "further (e.g. 0.02-0.05) for a structure relaxed with a looser "
                             "force tolerance.")
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
        print_dual(f"Displacement      : {args.displacement} Ang (non-bulk/hybrid path only)", f_out)
        print_dual(f"FC displacement   : {args.fc_displ} Bohr (bulk/hybrid path only)", f_out)
        print_dual(f"Polarization grid : {' '.join(str(v) for v in args.polarization_grid)} "
                    "(bulk path only -- hybrid path auto-builds its own vacuum-aware grid, "
                    "see [1])", f_out)
        print_dual(f"Symmetry tolerance: symprec={args.symprec:g} Ang", f_out)

        print_section('[1] PHONON MODES AT GAMMA', f_out)
        phonon, internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            phonon_dir, system_label, has_embedded_fc, f_out, symprec=args.symprec)
        try:
            unique_elements = list(set(phonon.primitive.symbols))
            lattice_ang = np.array(phonon.primitive.cell) * internal_to_angstrom
            frac_coords = np.array(phonon.primitive.scaled_positions)
            vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
            n_vacuum_axes = sum(vacuum_axes)
            is_0d = n_vacuum_axes == 3
            is_bulk = n_vacuum_axes == 0
            # Exactly one vacuum axis == a 2D slab -- the HYBRID path (see VERSION comment
            # above). Two vacuum axes (1D wire) deliberately stays on the plain non-bulk path
            # for now (out of scope, see the argparse description's own 1D caveat).
            is_hybrid_2d = n_vacuum_axes == 1
            frequencies, mode_band_indices, n_extra_trivial = get_gamma_modes(
                phonon, exclude_acoustic=True,
                extra_trivial_tol_thz=args.rotational_mode_tol if is_0d else None)
            mode_symmetries, point_group, symmetry_error = classify_ir_modes(phonon)
            space_group = phonon.symmetry.get_international_table() or "unknown"
            n_sym_ops = len(phonon.symmetry.symmetry_operations['rotations'])
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
        elif is_hybrid_2d:
            print_dual(color_text(
                "[Path] Hybrid (2D slab, one vacuum-padded axis) -- one shared equilibrium "
                "Born-effective-charge SIESTA run covers the two in-plane dmu/dQ components "
                "(PolarizationGrids built vacuum-aware, out-of-plane row zeroed), PLUS a "
                "+/-delta dipole-moment displacement pair per selected mode for the "
                "vacuum-axis component (a naive dipole is only valid along that non-periodic "
                "direction). See --help for why.", 'cyan'), f_out)
        else:
            print_dual(color_text(
                "[Path] Non-bulk (0D molecule / 1D wire) -- a +/-delta dipole-moment "
                "displacement pair per selected mode (plain single-point SCF, SIESTA prints "
                "the total dipole automatically).", 'cyan'), f_out)
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
        band_to_k = {int(band_idx): k for k, band_idx in enumerate(mode_band_indices, start=1)}
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            flag = color_text(" [IMAGINARY]", 'red') if freq < 0 else ""
            ms = band_to_symmetry.get(int(band_idx))
            sym_note = ""
            if ms is not None:
                label_str = f" ({ms.label})" if ms.label else ""
                sym_note = (color_text(f"{label_str} [IR-active]", 'green') if ms.is_ir_active
                            else color_text(f"{label_str} [symmetry-forbidden]", 'yellow'))
            eigendisp = mode_eigendisplacement(phonon, int(band_idx), internal_to_angstrom)
            character = mode_character_label(eigendisp)
            character_note = f"  [{character}]" if character else ""
            print_dual(f"  mode {k:3d} (band {int(band_idx):3d}) : {freq:10.4f} THz{flag} "
                       f"{sym_note}{character_note}", f_out)

        if is_bulk:
            print_dual(color_text(
                "[LIMITATION] These Gamma-point frequencies come from real-space "
                "finite-displacement force constants with NO non-analytic (long-range "
                "dipole-dipole, 'LO-TO') correction applied. For a non-polar bulk crystal "
                "(Born effective charges ~0, e.g. Si, diamond, graphene) this is exact. For a "
                "POLAR bulk crystal (ionic/partially-ionic bonding, nonzero Z*), a finite "
                "periodic supercell cannot fully decay the 1/r^3 dipole-dipole interaction "
                "between a displaced atom and its own periodic images, so the frequency "
                "reported here for an IR-active mode is contaminated and can land noticeably "
                "BETWEEN the true TO and LO frequencies rather than matching either -- verified "
                "on NaCl's T1u mode: 185.6 cm^-1 computed vs. ~164 cm^-1 (TO) / ~264 cm^-1 (LO) "
                "from inelastic-neutron-scattering literature (Raunio, Almqvist & Stedman, "
                "Phys. Rev. 178, 1496 (1969)), landing almost exactly midway. A correct "
                "treatment (Gonze & Lee, Phys. Rev. B 55, 10355 (1997); see also Togo & Tanaka, "
                "J. Phys. Soc. Jpn. 92, 012001 (2023) Sec. on non-analytic term correction) "
                "needs the high-frequency dielectric tensor (eps_inf, not currently computed by "
                "any stb tool) in addition to the Born charges already computed in this stage's "
                "[Path] Bulk run, and would subtract the spurious supercell-truncated "
                "dipole-dipole contribution from the force constants before diagonalizing. Not "
                "yet implemented in stb-ir -- planned for a future version. Until then, treat "
                "a bulk polar crystal's reported Gamma frequency as a rough (TO,LO) bracket "
                "midpoint, not a clean TO value, and prefer comparing to experiment only "
                "qualitatively (mode symmetry, activity, relative ordering) rather than "
                "quantitatively.", 'yellow'), f_out)

        print_section('[1b] SYMMETRY ANALYSIS', f_out)
        print_dual(f"Symmetry precision: symprec={args.symprec:g} Ang", f_out)
        print_dual(f"Space group       : {space_group}", f_out)
        print_dual(f"Point group       : {point_group}", f_out)
        print_dual(f"Symmetry ops      : {n_sym_ops}", f_out)
        # len(ms.band_indices) > 1 alone also matches a degenerate ACOUSTIC
        # group (e.g. NaCl's T1u translation at Gamma) -- exclude any group
        # with no band in band_to_k (i.e. nothing among the displayed
        # non-acoustic modes above), or it prints an empty "modes  (T1u)"
        # line for a group the user never sees listed anywhere else.
        degenerate_groups = [
            (ms, sorted(band_to_k[b] for b in ms.band_indices if b in band_to_k))
            for ms in mode_symmetries if len(ms.band_indices) > 1
        ]
        degenerate_groups = [(ms, ids) for ms, ids in degenerate_groups if ids]
        if degenerate_groups:
            print_dual(f"Degenerate groups : {len(degenerate_groups)} group(s) by symmetry "
                        "(same irrep -- Phonopy's own choice of basis within each group is "
                        "arbitrary, treat individual partners' directions accordingly):", f_out)
            for ms, mode_ids in degenerate_groups:
                label_str = f" ({ms.label})" if ms.label else ""
                print_dual(f"  modes {', '.join(str(m) for m in mode_ids)}{label_str}", f_out)
        if symmetry_error:
            print_dual(color_text(
                f"[WARNING] Symmetry classification unavailable ({symmetry_error}) -- "
                "every mode is kept/probed, same as without --use-symmetry. If your structure "
                "is a conventional (non-primitive) cell, try reducing it first with "
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
            if args.use_symmetry and args.modes is not None:
                print_dual(color_text(
                    "--modes was given explicitly -- symmetry filtering is informational "
                    "only here, no mode is auto-skipped.", 'yellow'), f_out)
            if not args.use_symmetry:
                print_dual(
                    "(--use-symmetry not given -- the labels/tags above are informational "
                    "only, every mode is still probed regardless of activity.)", f_out)
        print_dual(color_text(
            "\n[CHECK] Confirm the space/point group above is the one you actually expect for "
            "this crystal (e.g. from its known/published structure) BEFORE trusting the "
            "IR-active/forbidden labels or any near-degenerate mode split reported above. A "
            "too-tight --symprec silently detects a LOWER symmetry than the real one -- it "
            "will not raise an error, it will just misclassify some modes as active/forbidden "
            "and split a truly degenerate mode into two slightly different frequencies. If the "
            "reported group looks wrong (e.g. an orthorhombic label like 'mmm' for what should "
            "be a hexagonal/cubic crystal), rerun this stage with a looser --symprec (e.g. 0.02, "
            "0.05, 0.1) until the expected group is recovered -- and re-check that a much looser "
            "value doesn't overshoot into a HIGHER symmetry than the real (possibly slightly "
            "distorted) relaxed structure actually has.", 'yellow'), f_out)

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
        n_passed_filters = 0
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            if args.modes is not None and k not in args.modes:
                continue
            if args.freq_min is not None and freq < args.freq_min:
                continue
            if args.freq_max is not None and freq > args.freq_max:
                continue
            n_passed_filters += 1
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

        all_symmetry_forbidden = False
        if not selected:
            # Every candidate mode being symmetry-forbidden is a valid physical
            # result (e.g. a centrosymmetric crystal where every Gamma mode
            # happens to be IR-silent, per the rule of mutual exclusion) -- not
            # a misconfiguration, so it's a [WARNING] and the run finishes
            # cleanly with zero folders rather than exiting non-zero. Only
            # downgrade when --modes/--freq-min/--freq-max weren't ALSO
            # responsible (n_passed_filters > 0 means symmetry alone emptied
            # the selection).
            all_symmetry_forbidden = (apply_symmetry_skip and n_passed_filters > 0
                                       and n_actually_skipped == n_passed_filters)
            if all_symmetry_forbidden:
                print_dual(color_text(
                    "\n[WARNING] Every mode that passed --modes/--freq-min/--freq-max is "
                    "symmetry-forbidden from being IR-active -- nothing to compute. This is a "
                    "valid physical result (e.g. a centrosymmetric crystal with no IR-active "
                    "Gamma modes), not an error. No displacement folders written.", 'yellow'), f_out)
            else:
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

        if (is_bulk or is_hybrid_2d) and not selected:
            print_dual(
                "Folders to write  : 0 -- every candidate mode is symmetry-forbidden from "
                "being IR-active, so the equilibrium Born-effective-charge run would have "
                "nothing to report on.", f_out)
            n_real_folders = 0
        elif is_bulk:
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
                    frames = build_mode_animation_frames(
                        phonon, band_idx, args.displacement, internal_to_angstrom,
                        n_frames=args.animation_frames)
                    write_mode_animation(frames, animation_path)
                    print_dual(f"  {color_text('[OK]', 'green')} {animation_path}", f_out)

                eigendisp = mode_eigendisplacement(phonon, band_idx, internal_to_angstrom)
                sidecar_path = os.path.join(born_charge_root, f"mode_{k:02d}_eigendisplacement.json")
                with open(sidecar_path, "w") as f:
                    json.dump({"eigendisplacement": eigendisp.tolist()}, f)

                label = f"mode_{k:02d}_born_charge"
                report_rows.append((label, k, band_idx, freq, "BULK", "-", equilibrium_dir, None))

            n_real_folders = 1
        elif is_hybrid_2d:
            n_representative = sum(1 for k, _, _ in selected if representative_of[k] == k)
            print_dual(f"Folders to write  : 1 equilibrium Born-effective-charge SIESTA run "
                        f"(in-plane components) + {n_representative * 2} independent "
                        f"single-point SIESTA runs (2 signs x per-mode dipole displacement, "
                        f"vacuum-axis component), across {n_representative} representative "
                        f"mode(s) out of {len(selected)} selected).", f_out)

            equilibrium_dir = os.path.join(born_charge_root, "equilibrium")
            slab_grid = build_slab_polarization_grid(vacuum_axes)
            calc_text = force_born_charge_run(calc_template, args.fc_displ, slab_grid)
            write_ir_folder(equilibrium_dir, phonon.primitive, structure_filename, calc_text, pseudos)
            print_dual(f"  {color_text('[OK]', 'green')} {equilibrium_dir} "
                        f"(PolarizationGrids {slab_grid})", f_out)

            os.makedirs(born_charge_root, exist_ok=True)
            for k, freq, band_idx in selected:
                if args.export_animations:
                    animation_path = os.path.join(dipole_root, f"mode_{k:02d}_animation.axsf")
                    frames = build_mode_animation_frames(
                        phonon, band_idx, args.displacement, internal_to_angstrom,
                        n_frames=args.animation_frames)
                    write_mode_animation(frames, animation_path)
                    print_dual(f"  {color_text('[OK]', 'green')} {animation_path}", f_out)

                # The Born-charge/eigendisplacement combination is free (shared equilibrium
                # run, no extra SIESTA cost) -- write it for every selected mode regardless of
                # --skip-degenerate, same "no folder cost to save" reasoning the bulk path
                # already uses. Only the per-mode dipole PAIR (the real SIESTA cost) is
                # actually skippable for a derived mode.
                eigendisp = mode_eigendisplacement(phonon, band_idx, internal_to_angstrom)
                sidecar_path = os.path.join(born_charge_root, f"mode_{k:02d}_eigendisplacement.json")
                with open(sidecar_path, "w") as f:
                    json.dump({"eigendisplacement": eigendisp.tolist()}, f)
                report_rows.append((f"mode_{k:02d}_born_charge", k, band_idx, freq, "HYBRID",
                                     "-", equilibrium_dir, None))

                if representative_of[k] != k:
                    rep_k = representative_of[k]
                    report_rows.append((f"mode_{k:02d}_dipole_derived", k, band_idx, freq,
                                         "HYBRID", "DERIVED", "-", rep_k))
                    continue

                for sign_name, sign_val in (("plus", 1.0), ("minus", -1.0)):
                    displaced = displace_along_mode(
                        phonon, band_idx, args.displacement, internal_to_angstrom, sign=sign_val)
                    label = f"mode_{k:02d}_{sign_name}"
                    mode_dir = os.path.join(dipole_root, label)
                    dipole_calc_text = force_single_point(calc_template)
                    write_ir_folder(mode_dir, displaced, structure_filename, dipole_calc_text, pseudos)
                    report_rows.append((label, k, band_idx, freq, "HYBRID", sign_name, mode_dir, None))
                    print_dual(f"  {color_text('[OK]', 'green')} {mode_dir}", f_out)

            n_real_folders = 1 + sum(1 for row in report_rows
                                      if row[4] == "HYBRID" and row[5] in ("plus", "minus"))
        else:
            n_folders = sum(1 for k, _, _ in selected if representative_of[k] == k) * 2
            print_dual(f"Folders to write  : {n_folders} independent single-point SIESTA runs "
                        f"(2 signs x per-mode dipole displacement, across "
                        f"{n_folders // 2} representative mode(s) out of {len(selected)} "
                        "selected).", f_out)

            for k, freq, band_idx in selected:
                if args.export_animations:
                    animation_path = os.path.join(dipole_root, f"mode_{k:02d}_animation.axsf")
                    frames = build_mode_animation_frames(
                        phonon, band_idx, args.displacement, internal_to_angstrom,
                        n_frames=args.animation_frames)
                    write_mode_animation(frames, animation_path)
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
        if n_real_folders:
            print_dual(color_text("\nNext steps:", 'yellow'), f_out)
            if is_bulk:
                print_dual(f"  1. Run SIESTA in '{born_charge_root}/equilibrium/'.", f_out)
            elif is_hybrid_2d:
                print_dual(f"  1. Run SIESTA in '{born_charge_root}/equilibrium/' AND in every "
                            f"'{dipole_root}/mode_*/' folder.", f_out)
            else:
                print_dual(f"  1. Run SIESTA in every '{dipole_root}/mode_*/' folder.", f_out)
            print_dual(f"  2. Once they're done, run: stb-irAnalysis --directory {output_root}", f_out)
        else:
            print_dual(color_text(
                "\nNo SIESTA runs needed -- every candidate mode was symmetry-forbidden from "
                "being IR-active (see [1b] above). Nothing further to do for this structure's "
                "IR spectrum.", 'yellow'), f_out)

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
    if n_real_folders:
        print(color_text("Displacement folder(s) ready for Stage 3 (stb-irAnalysis).\n", 'bold'))
    else:
        print(color_text(
            "No displacement folders were needed -- every candidate mode is symmetry-forbidden "
            "from being IR-active.\n", 'bold'))


if __name__ == "__main__":
    main()
