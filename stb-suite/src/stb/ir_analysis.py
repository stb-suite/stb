#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.2.1"  # New [LIMITATION] note at the top of [2] IR-ACTIVE MODES SUMMARY, shown
                    # whenever any mode used the BULK path: those frequencies have no
                    # non-analytic (LO-TO) correction, so for a polar bulk crystal they can land
                    # between the true TO and LO values -- see stb-irModes' VERSION 1.4.1 comment
                    # for the full root-cause story and literature references. Documented as a
                    # known limitation only; a real fix is planned for a future version.
                    # (previously 1.2.0: New HYBRID path (2D slab, exactly one vacuum axis): combines the BULK
                    # path's Born-charge x eigendisplacement formula (in-plane/periodic axes)
                    # with the NONBULK path's dipole-difference formula (vacuum axis) into one
                    # physically valid dmu/dQ per mode -- see stb-irModes' own VERSION comment
                    # for the full root-cause story (a naive dipole moment silently came back
                    # exactly zero for h-BN's in-plane E' mode, a mode the literature puts as
                    # the dominant IR peak). The two contributions are additive by construction
                    # (each is ~0 exactly where the other is valid); axis_mask selects/verifies
                    # this when the vacuum axis's Cartesian alignment can be confirmed, and
                    # falls back to plain addition with a warning when it can't.
                    # Also: [3] SPECTRUM now always lists detected peaks (not just with
                    # --experimental), [2] reports combined intensity for near-degenerate mode
                    # groups (the basis-independent quantity, since individual dmu/dQ directions
                    # within a degenerate subspace are an arbitrary Phonopy eigenvector choice),
                    # the gnuplot .dat/.gplot pair is now opt-in via --save-gnuplot (written
                    # under <directory>/plot/, same convention stb-ramanAnalysis already uses)
                    # instead of always-on, --view opens an interactive matplotlib preview, and
                    # a new [3c] MODE VIBRATIONS section always reloads Stage 1's force constants
                    # to write an extended-XYZ animation per analyzed mode under
                    # <directory>/mode_animations/ (skipped gracefully, never a hard error, if
                    # phonon_disp/ is no longer present) -- --view-modes additionally opens each
                    # one interactively in ASE's viewer. build_mode_animation_frames/
                    # write_mode_animation/phonopy_atoms_to_ase moved to core/phonon_workflow.py
                    # (were duplicated in raman_modes.py/ir_modes.py) now that this module is a
                    # third consumer -- extract-on-second-use, per this suite's own policy.

import os
import re
import sys
import glob
import json
import argparse
from datetime import datetime
import numpy as np
import yaml
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.siesta_log import get_electric_dipole, check_scf_and_force
from stb.core.born_charges import read_born_charges
from stb.core.spectrum import (
    bose_einstein_weight, build_lorentzian_spectrum, read_experimental_spectrum,
    find_spectrum_peaks, match_peaks,
)
from stb.core.phonon_workflow import (
    load_phonon_with_force_constants, build_mode_animation_frames, write_mode_animation,
)
from stb.core.ase_view import view_structure_interactive

REPORT_FILE = "ir_stage3.txt"
# Same conventional spectroscopy unit Raman's Stage 3 uses -- Phonopy's own
# native unit (and Stage 2's MODE_TABLE) is THz throughout the rest of this
# workflow.
THZ_TO_CM1 = 33.35641

_LABEL_RE = re.compile(r'^\s*SystemLabel\s+(\S+)', re.IGNORECASE | re.MULTILINE)


def detect_ir_label(directory):
    """SystemLabel from the first available calc.fdf -- checks the
    non-bulk path's dipole_disp/mode_*/ folders first, then the bulk
    path's born_charge_disp/equilibrium/ folder, since a given
    stb-irModes run only ever wrote one or the other. Same regex-based
    approach as raman_analysis.detect_optical_label.
    """
    candidates = sorted(glob.glob(os.path.join(directory, "dipole_disp", "mode_*")))
    equilibrium = os.path.join(directory, "born_charge_disp", "equilibrium")
    if not candidates and os.path.isdir(equilibrium):
        candidates = [equilibrium]
    if not candidates:
        return None
    fdf_path = os.path.join(candidates[0], "calc.fdf")
    try:
        with open(fdf_path) as f:
            match = _LABEL_RE.search(f.read())
    except OSError:
        return None
    return match.group(1) if match else None


def read_mode_table(report_path):
    """Parses Stage 2's # MODE_TABLE (label, mode_index, band_index,
    frequency_thz, path, sign, dir, [derived_from]) into a list of dicts
    -- same "read the persisted # XXX_TABLE, don't infer from folder
    names" convention as raman_analysis.read_mode_table. `path` is
    "BULK", "NONBULK", or "HYBRID" (replaces Raman's per-axis "axis"
    column entirely -- IR has no probe-axis concept at all).
    """
    rows = []
    in_table = False
    with open(report_path) as f:
        for line in f:
            if line.startswith("# MODE_TABLE"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            rows.append({
                "label": parts[0],
                "mode_index": int(parts[1]),
                "band_index": int(parts[2]),
                "frequency_thz": float(parts[3]),
                "path": parts[4],
                "sign": parts[5],
                "dir": parts[6],
                "derived_from": (int(parts[7]) if len(parts) >= 8 and parts[7] != "-" else None),
            })
    return rows


def group_by_mode(rows):
    """{mode_index: {"frequency_thz": f, "band_index": b,
    "path": "BULK"/"NONBULK"/"HYBRID", "folders": {"plus": dir, "minus": dir}
    (non-bulk), {"-": dir} (bulk, a single shared equilibrium folder), or
    {"-": dir, "plus": dir, "minus": dir} (hybrid -- both at once),
    "derived_from": int or None}}. A DERIVED sign row (--skip-degenerate)
    contributes no folder of its own -- for HYBRID this means the mode's
    own "-" (Born-charge) entry is still present (always written, see
    stb-irModes) while "plus"/"minus" are missing (skipped, reused from the
    representative instead). `band_index` is the same 0-based Phonopy band
    for every row of a given mode_index -- carried through so a later stage
    (mode-animation export) can reload the phonon object and reconstruct
    this exact mode's eigendisplacement without re-deriving it from scratch.
    """
    modes = {}
    for row in rows:
        entry = modes.setdefault(row["mode_index"], {
            "frequency_thz": row["frequency_thz"], "band_index": row["band_index"],
            "path": row["path"], "folders": {}, "derived_from": None})
        if row["derived_from"] is not None:
            entry["derived_from"] = row["derived_from"]
        if row["sign"] != "DERIVED":
            entry["folders"][row["sign"]] = row["dir"]
    return modes


def dipole_axis_mask(vacuum_axes, lattice_ang, tol_deg=2.0):
    """For each vacuum-padded lattice vector, checks whether it's nearly
    parallel to a single Cartesian axis (true for the overwhelming
    majority of structures this suite builds -- stb-slab/stb-molecule/
    etc. are axis-aligned by construction). Returns a (3,) boolean array
    -- which CARTESIAN axes are non-periodic ("safe", the dipole
    component along them is the physically meaningful one SIESTA
    reports) -- or None if the alignment assumption fails for any
    vacuum-padded lattice vector, or the needed metadata wasn't
    recoverable at all. Callers must fall back to the raw, unmasked
    dipole vector (with a prominent warning) when this returns None,
    rather than masking with unverified assumptions.
    """
    if vacuum_axes is None or lattice_ang is None or len(lattice_ang) != 3:
        return None
    lattice_ang = np.asarray(lattice_ang, dtype=float)
    cartesian_nonperiodic = np.zeros(3, dtype=bool)
    for i, is_vac in enumerate(vacuum_axes):
        if not is_vac:
            continue
        vec = lattice_ang[i]
        norm = np.linalg.norm(vec)
        if norm < 1e-8:
            return None
        unit = np.abs(vec) / norm
        dominant = int(np.argmax(unit))
        if unit[dominant] < np.cos(np.radians(tol_deg)):
            return None
        cartesian_nonperiodic[dominant] = True
    return cartesian_nonperiodic


def ir_intensity(dmu_dq):
    """I_IR ~ |dmu/dQ|^2 -- the same formula for both the non-bulk
    (direct dipole finite-difference) and bulk (Born-charge x
    eigendisplacement projection) paths, the only place their differing
    derivations converge. Arbitrary units, same "relative spectrum, not
    an absolute cross-section" convention as Raman's activity.
    """
    return float(np.sum(np.asarray(dmu_dq, dtype=float) ** 2))


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 3 of 3: analyzes a stb-irModes displacement sweep and "
        "computes the IR spectrum.", 'bold')}
Auto-detects, per mode, which path Stage 2 used (from the MODE_TABLE's own
`path` column, never inferred from which files happen to exist):

  - NONBULK: reads the total dipole moment (SIESTA's own "Electric dipole
    (a.u.)" line) from each mode's +/-delta dipole_disp/mode_*/ folder,
    central-differences to get dmu/dQ directly as a 3-vector.

  - BULK: reads the Born effective charge tensors (SystemLabel.BC) from the
    single shared born_charge_disp/equilibrium/ folder, combines with each
    mode's own eigendisplacement (a JSON sidecar Stage 2 already wrote) via
    dmu_beta/dQ = sum_atom sum_tau Z*[atom,tau,beta] * e[atom,tau].

  - HYBRID (2D slab): both of the above at once -- the in-plane dmu/dQ
    components come from the same Born-charge x eigendisplacement formula as
    BULK (computed with a vacuum-aware PolarizationGrids that only covers
    the two periodic axes), the vacuum-axis component comes from the same
    dipole-difference formula as NONBULK. A naive dipole is only valid along
    a non-periodic direction, and Born-effective-charge/Berry-phase is only
    valid along a periodic one -- neither method alone covers all 3 axes of
    a slab correctly, so this combines one component from each.

Either way, IR intensity is |dmu/dQ|^2 -- same formula, only the derivation
differs. Builds a Lorentzian-summed spectrum and always lists its detected
peaks; the gnuplot .dat/.gplot pair itself is opt-in (--save-gnuplot, under
<directory>/plot/) and an on-screen matplotlib preview is available too
(--view). Also always reloads Stage 1's force constants to write an
extended-XYZ vibration animation per analyzed mode under
<directory>/mode_animations/ (skipped gracefully if phonon_disp/ is gone),
optionally opened one at a time in ASE's interactive viewer (--view-modes).
Near-degenerate mode groups get a combined, basis-independent intensity
alongside their individual (basis-dependent) values.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --directory ir_study\n"
               "  %(prog)s --directory ir_study --save-gnuplot --view\n"
               "  %(prog)s --directory ir_study --view-modes --animation-frames 30\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="ir_study",
                        help="Root directory written by stb-ir/stb-irModes (default: ir_study).")
    parser.add_argument("-l", "--label", type=str, default=None,
                        help="SystemLabel used in Stage 2's calc.fdf (default: auto-detected).")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each displacement folder (default: calc.out).")
    parser.add_argument("--linewidth", type=float, default=10.0,
                        help="Lorentzian linewidth (cm^-1) for the summed spectrum (default: 10.0).")
    parser.add_argument("--temperature", type=float, default=None, metavar="K",
                        help="Weight each mode's spectrum contribution by the Stokes "
                             "Bose-Einstein thermal population factor (n+1) at this temperature "
                             "(Kelvin). Default: no weighting (every mode contributes by its raw "
                             "intensity alone, as before this option existed).")
    parser.add_argument("-o", "--output", type=str, default="ir_spectrum",
                        help="Base filename (no extension) for the spectrum .dat/.gplot pair "
                             "(default: ir_spectrum).")
    parser.add_argument("--experimental", type=str, default=None, metavar="FILE",
                        help="Overlay a measured IR spectrum (2-column text file: frequency "
                             "cm^-1, intensity) on the simulated one, and report a nearest-"
                             "neighbor peak match (simulated vs. experimental peak positions, "
                             "one row per experimental peak). Plotted on a separate y-axis "
                             "(simulated intensity and real intensity aren't on the same scale).")
    parser.add_argument("--peak-prominence", type=float, default=0.01, metavar="FRACTION",
                        help="Minimum peak prominence, as a fraction of that spectrum's own max "
                             "intensity, for peak detection (both the standalone [3] SPECTRUM "
                             "peak list and --experimental's matching) (default: 0.01 -- filters "
                             "out sub-1%% noise bumps).")
    parser.add_argument("--degenerate-tol", type=float, default=1e-3, metavar="THZ",
                        help="Two modes within this frequency tolerance (THz) are reported as a "
                             "degenerate group with a combined intensity, in addition to their "
                             "individual values -- |dmu/dQ|^2 is not itself basis-independent "
                             "within a degenerate subspace (Phonopy's returned eigenvectors are "
                             "an arbitrary orthogonal basis of it), but the GROUP's summed "
                             "intensity is the physically meaningful, basis-independent quantity "
                             "(default: 0.001 THz, comfortably above float round-trip noise "
                             "through the MODE_TABLE's 6-decimal text format, comfortably below "
                             "any real frequency splitting between genuinely distinct modes).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <directory>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also save the spectrum (and, with --experimental, the overlaid "
                             "measured spectrum) as gnuplot .dat + .gplot scripts, under "
                             "'<directory>/plot/'. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="View the IR spectrum (and, with --experimental, the overlaid "
                             "measured spectrum on a second y-axis) interactively via matplotlib "
                             "now. Off by default. Needs a display. Independent of --save-gnuplot "
                             "-- matplotlib here is only an on-screen preview, never written to disk.")
    parser.add_argument("--view-modes", action="store_true",
                        help="Open each analyzed mode's eigendisplacement animation interactively "
                             "in ASE's 3D viewer, one at a time (close a window to advance to the "
                             "next) -- needs a display and Stage 1's 'phonon_disp/' folder to "
                             "still be present (reloads the same force constants Stage "
                             "2/this stage's own always-on XYZ export already use). Off by default.")
    parser.add_argument("--animation-frames", type=int, default=20,
                        help="Frames per mode animation, for both --view-modes and the always-on "
                             "extended-XYZ export below (default: 20).")
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry-detection tolerance (Ang), only used when reloading the "
                             "phonon force constants for mode animations (--view-modes and the "
                             "always-on XYZ export) -- same default as every other stage of this "
                             "workflow (0.01, pymatgen's own default). Does not affect the "
                             "IR-active/inactive classification (that was already decided in "
                             "Stage 2) or the spectrum itself.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-irAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("IR SPECTRUM WORKFLOW -- STAGE 3: ANALYSIS", 'bold'))
    print("-" * 60)

    stage2_report = os.path.join(args.directory, "ir_stage2.txt")
    if not os.path.isfile(stage2_report):
        print(color_text(f"[ERROR] '{stage2_report}' not found -- run stb-irModes (Stage 2) first.", 'red'))
        sys.exit(1)

    if args.experimental is not None and not os.path.isfile(args.experimental):
        print(color_text(f"[ERROR] --experimental file '{args.experimental}' not found.", 'red'))
        sys.exit(1)

    label = args.label or detect_ir_label(args.directory)
    if not label:
        print(color_text(
            f"[ERROR] Could not detect SystemLabel from '{args.directory}' -- pass -l/--label.", 'red'))
        sys.exit(1)

    # Recovered from Stage 2's own report text, same "read it back, don't
    # re-derive" convention raman_analysis.py already uses for its own
    # delta: the finite-difference displacement (non-bulk path), the raw
    # vacuum-axis booleans, and the lattice matrix (both needed for the
    # non-bulk path's dipole-axis alignment check -- see
    # dipole_axis_mask's own docstring).
    stage2_delta = None
    vacuum_axes = None
    lattice_ang = []
    in_lattice_block = False
    with open(stage2_report) as f:
        for line in f:
            if line.startswith("Displacement      :"):
                try:
                    stage2_delta = float(line.split(":")[1].split()[0])
                except (IndexError, ValueError):
                    pass
            elif line.startswith("Vacuum axes (abc) :"):
                try:
                    parts = line.split(":")[1].split()
                    vacuum_axes = [p == "True" for p in parts[:3]]
                except (IndexError, ValueError):
                    pass
            elif line.startswith("Lattice vectors (Ang"):
                in_lattice_block = True
            elif in_lattice_block:
                parts = line.split()
                if len(parts) == 4 and parts[0] in ("a", "b", "c"):
                    try:
                        lattice_ang.append([float(parts[1]), float(parts[2]), float(parts[3])])
                    except ValueError:
                        pass
                if len(lattice_ang) >= 3:
                    in_lattice_block = False
    if stage2_delta is None:
        print(color_text(
            f"[ERROR] Could not recover the Stage 2 displacement value from '{stage2_report}'.", 'red'))
        sys.exit(1)

    axis_mask = dipole_axis_mask(vacuum_axes, lattice_ang)

    rows = read_mode_table(stage2_report)
    if not rows:
        print(color_text(f"[ERROR] No # MODE_TABLE rows found in '{stage2_report}'.", 'red'))
        sys.exit(1)
    modes = group_by_mode(rows)

    born_charge_root = os.path.join(args.directory, "born_charge_disp")
    born_charges_cache = None  # (atom_labels, Z_star), read once, shared by every BULK mode

    report_path = os.path.join(args.directory, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    print_dual(f"{color_text('===== IR STAGE 3 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory         : {args.directory}", f_out)
    print_dual(f"SystemLabel       : {label}", f_out)
    print_dual(f"Modes in table    : {len(modes)}", f_out)
    print_dual(f"Displacement      : {stage2_delta} Ang (from Stage 2, non-bulk path only)", f_out)
    if axis_mask is None and vacuum_axes is not None and any(vacuum_axes):
        print_dual(color_text(
            "[WARNING] Could not verify that a vacuum-padded lattice vector is aligned to a "
            "single Cartesian axis -- the non-bulk dipole derivative below uses the FULL, "
            "unmasked 3-vector. Periodic-direction component(s) may be gauge-ambiguous.",
            'yellow'), f_out)

    print_section('[1] READING RUNS', f_out)

    mode_results = []  # (mode_index, freq_thz, dmu_dq, intensity)
    computed = {}  # mode_index -> (dmu_dq, intensity, scf_ok), for --skip-degenerate reuse
    any_bulk_path = False  # any mode using the BULK (3D polar-crystal) intensity path -- see [2]'s LO-TO caveat
    for mode_index in sorted(modes):
        entry = modes[mode_index]
        freq_thz = entry["frequency_thz"]
        print_dual(f"  Mode {mode_index} ({freq_thz:.4f} THz / {freq_thz * THZ_TO_CM1:.2f} cm^-1):", f_out)

        derived_from = entry.get("derived_from")
        if derived_from is not None:
            rep = computed.get(derived_from)
            if rep is None:
                print_dual(color_text(
                    f"    -> SKIP (representative mode {derived_from} has no usable dipole "
                    "derivative)", 'yellow'), f_out)
                continue
            dmu_dq, intensity, scf_ok = rep
            print_dual("    [Degenerate partner] Reusing mode "
                        f"{derived_from}'s dipole derivative by symmetry "
                        "(stb-irModes --skip-degenerate) -- |dmu/dQ|^2 is a rotational "
                        "invariant, identical for every partner in a degenerate group.", f_out)
            print_dual(f"    -> dmu/dQ=({dmu_dq[0]:.6f}, {dmu_dq[1]:.6f}, {dmu_dq[2]:.6f})  "
                        f"intensity~{intensity:.6f}  [derived]", f_out)
            mode_results.append((mode_index, freq_thz, dmu_dq, intensity))
            computed[mode_index] = rep
            continue

        if entry["path"] in ("BULK", "HYBRID"):
            if born_charges_cache is None:
                equilibrium_dir = entry["folders"].get("-")
                bc_path = os.path.join(equilibrium_dir, f"{label}.BC")
                result = read_born_charges(bc_path)
                if result is None:
                    print_dual(color_text(
                        f"    [SKIP] Could not read Born effective charges from "
                        f"'{bc_path}'.", 'yellow'), f_out)
                    born_charges_cache = (None, None)
                else:
                    atom_labels, Z_star = result
                    scf_ok_bulk, _max_force = check_scf_and_force(
                        os.path.join(equilibrium_dir, args.file))
                    if not scf_ok_bulk:
                        print_dual(color_text(
                            f"    [WARNING] Could not confirm SCF convergence for "
                            f"'{equilibrium_dir}' -- every mode's Born-charge-derived "
                            "intensity may be unreliable.", 'yellow'), f_out)
                    born_charges_cache = (Z_star, scf_ok_bulk)

            Z_star, scf_ok_bulk = born_charges_cache
            if Z_star is None:
                print_dual(color_text("    -> SKIP (Born effective charges unavailable)", 'yellow'), f_out)
                continue

            sidecar_path = os.path.join(born_charge_root, f"mode_{mode_index:02d}_eigendisplacement.json")
            try:
                with open(sidecar_path) as f:
                    eigendisp = np.array(json.load(f)["eigendisplacement"])
            except (OSError, ValueError, KeyError) as e:
                print_dual(color_text(
                    f"    [SKIP] Could not read {sidecar_path}: {e}", 'yellow'), f_out)
                continue

            if eigendisp.shape[0] != Z_star.shape[0]:
                print_dual(color_text(
                    f"    [SKIP] Atom count mismatch between eigendisplacement "
                    f"({eigendisp.shape[0]}) and Born charges ({Z_star.shape[0]}).",
                    'yellow'), f_out)
                continue

            dmu_dq_bec = np.einsum('atb,at->b', Z_star, eigendisp)

            if entry["path"] == "BULK":
                any_bulk_path = True
                dmu_dq = dmu_dq_bec
                scf_ok = scf_ok_bulk
                unconverged_note = "  [equilibrium run unconverged]" if not scf_ok else ""
            else:
                # HYBRID: dmu_dq_bec already carries the two in-plane (periodic)
                # components (the vacuum axis's own PolarizationGrids row was
                # zeroed in Stage 2, so SIESTA never computed a component there --
                # it comes back ~0 by construction). The vacuum-axis component
                # instead comes from the same dipole-difference formula NONBULK
                # uses below, which is ~0 along the periodic axes for the mirror-
                # image reason (a naive dipole is gauge-ambiguous there -- see
                # VERSION comment history). The two contributions are therefore
                # additive by construction; axis_mask (when available) is used
                # only to select/verify that, not as the sole mechanism, so this
                # degrades gracefully if the alignment check itself is uncertain.
                plus_dir = entry["folders"].get("plus")
                minus_dir = entry["folders"].get("minus")
                if plus_dir is None or minus_dir is None:
                    print_dual(color_text(
                        "    -> SKIP (missing +/-delta dipole folder(s) for the vacuum-axis "
                        "component)", 'yellow'), f_out)
                    continue

                scf_ok_dipole = True
                dipoles = {}
                for sign, folder in (("plus", plus_dir), ("minus", minus_dir)):
                    out_path = os.path.join(folder, args.file)
                    ok, _max_force = check_scf_and_force(out_path)
                    scf_ok_dipole = scf_ok_dipole and ok
                    if not ok:
                        print_dual(color_text(
                            f"    [WARNING] Could not confirm SCF convergence for {folder} -- "
                            "this dipole value may be unreliable.", 'yellow'), f_out)
                    dipole = get_electric_dipole(out_path)
                    if dipole is None:
                        print_dual(color_text(
                            f"    [SKIP] Could not read the electric dipole from '{out_path}'.",
                            'yellow'), f_out)
                        dipoles = None
                        break
                    dipoles[sign] = dipole
                if not dipoles:
                    continue

                dmu_dq_dipole_full = (dipoles["plus"] - dipoles["minus"]) / (2.0 * stage2_delta)

                if axis_mask is not None:
                    dmu_dq = np.where(axis_mask, dmu_dq_dipole_full, dmu_dq_bec)
                    bec_leak = dmu_dq_bec[axis_mask]
                    dipole_leak = dmu_dq_dipole_full[~axis_mask]
                    if (bec_leak.size and np.max(np.abs(bec_leak)) > 1e-3) or \
                       (dipole_leak.size and np.max(np.abs(dipole_leak)) > 1e-3):
                        print_dual(color_text(
                            "    [WARNING] Unexpected non-negligible signal in the Born-charge "
                            "component along the vacuum axis, or in the dipole-difference "
                            "component along an in-plane axis -- one of the two HYBRID methods "
                            "may not be behaving as expected for this structure; inspect "
                            "dmu/dQ below before trusting it.", 'yellow'), f_out)
                else:
                    dmu_dq = dmu_dq_bec + dmu_dq_dipole_full
                    print_dual(color_text(
                        "    [WARNING] Could not verify the vacuum axis's Cartesian alignment "
                        "-- combining the Born-charge and dipole-difference components by "
                        "plain addition (each should be ~0 where the other is valid) without "
                        "being able to check that assumption here.", 'yellow'), f_out)

                scf_ok = scf_ok_bulk and scf_ok_dipole
                unconverged_note = "  [some folders unconverged]" if not scf_ok else ""

            intensity = ir_intensity(dmu_dq)
            print_dual(f"    -> dmu/dQ=({dmu_dq[0]:.6f}, {dmu_dq[1]:.6f}, {dmu_dq[2]:.6f})  "
                        f"intensity~{intensity:.6f}"
                        + ("" if not unconverged_note else color_text(unconverged_note, 'yellow')), f_out)
            mode_results.append((mode_index, freq_thz, dmu_dq, intensity))
            computed[mode_index] = (dmu_dq, intensity, scf_ok)

        else:  # NONBULK
            plus_dir = entry["folders"].get("plus")
            minus_dir = entry["folders"].get("minus")
            if plus_dir is None or minus_dir is None:
                print_dual(color_text(
                    "    -> SKIP (missing +/-delta folder(s))", 'yellow'), f_out)
                continue

            scf_ok = True
            dipoles = {}
            for sign, folder in (("plus", plus_dir), ("minus", minus_dir)):
                out_path = os.path.join(folder, args.file)
                ok, _max_force = check_scf_and_force(out_path)
                scf_ok = scf_ok and ok
                if not ok:
                    print_dual(color_text(
                        f"    [WARNING] Could not confirm SCF convergence for {folder} -- "
                        "this dipole value may be unreliable.", 'yellow'), f_out)
                dipole = get_electric_dipole(out_path)
                if dipole is None:
                    print_dual(color_text(
                        f"    [SKIP] Could not read the electric dipole from '{out_path}'.",
                        'yellow'), f_out)
                    dipoles = None
                    break
                dipoles[sign] = dipole
            if not dipoles:
                continue

            dmu_dq_full = (dipoles["plus"] - dipoles["minus"]) / (2.0 * stage2_delta)
            if axis_mask is not None:
                dmu_dq = np.where(axis_mask, dmu_dq_full, 0.0)
                periodic_component = dmu_dq_full[~axis_mask]
                if periodic_component.size and np.max(np.abs(periodic_component)) > 1e-3:
                    print_dual(color_text(
                        "    [WARNING] Non-negligible dipole-derivative component along a "
                        "periodic direction (masked out) -- check the structure's cell "
                        "alignment.", 'yellow'), f_out)
            else:
                dmu_dq = dmu_dq_full

            intensity = ir_intensity(dmu_dq)
            print_dual(f"    -> dmu/dQ=({dmu_dq[0]:.6f}, {dmu_dq[1]:.6f}, {dmu_dq[2]:.6f})  "
                        f"intensity~{intensity:.6f}"
                        + ("" if scf_ok else color_text("  [some folders unconverged]", 'yellow')), f_out)
            mode_results.append((mode_index, freq_thz, dmu_dq, intensity))
            computed[mode_index] = (dmu_dq, intensity, scf_ok)

    if not mode_results:
        print_dual(color_text(
            "\n[ERROR] No mode had a usable dipole derivative -- nothing to report/plot. "
            "Make sure SIESTA finished in every displacement folder.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_section('[2] IR-ACTIVE MODES SUMMARY', f_out)
    if any_bulk_path:
        print_dual(color_text(
            "[LIMITATION] Frequencies below (for the BULK/3D-periodic path) have no "
            "non-analytic (LO-TO) correction applied -- for a POLAR bulk crystal (nonzero "
            "Born effective charges) they can land between the true TO and LO frequencies "
            "rather than matching either. See stb-irModes' Stage 2 report ([1] PHONON MODES "
            "AT GAMMA) for the full explanation and literature references; not yet "
            "implemented in stb-ir, planned for a future version.", 'yellow'), f_out)
    print_dual(f"  {'Mode':<6}{'THz':<10}{'cm^-1':<10}{'dmu/dQ_x':<12}{'dmu/dQ_y':<12}"
                f"{'dmu/dQ_z':<12}{'Intensity'}", f_out)
    for mode_index, freq_thz, dmu_dq, intensity in mode_results:
        print_dual(f"  {mode_index:<6}{freq_thz:<10.4f}{freq_thz * THZ_TO_CM1:<10.2f}"
                    f"{dmu_dq[0]:<12.6f}{dmu_dq[1]:<12.6f}{dmu_dq[2]:<12.6f}{intensity:.6f}", f_out)

    # Degenerate groups: modes within --degenerate-tol of each other come from
    # the same physical vibration -- Phonopy's own choice of eigenvector basis
    # within a degenerate subspace is arbitrary (any orthogonal rotation of it
    # is equally valid), so each individual partner's dmu/dQ DIRECTION (and
    # therefore its own reported intensity) is basis-dependent, not itself
    # physically meaningful in isolation. Only the GROUP's summed intensity is
    # basis-independent (sum of |projection|^2 over an orthonormal basis of a
    # subspace is invariant under any orthogonal change of basis within it) --
    # this is exactly the quantity that should be compared against a single
    # experimental peak when several modes land on top of each other.
    sorted_results = sorted(mode_results, key=lambda r: r[1])
    degenerate_groups = []
    for entry in sorted_results:
        if degenerate_groups and abs(entry[1] - degenerate_groups[-1][-1][1]) <= args.degenerate_tol:
            degenerate_groups[-1].append(entry)
        else:
            degenerate_groups.append([entry])
    degenerate_groups = [g for g in degenerate_groups if len(g) > 1]
    if degenerate_groups:
        print_dual(f"\n{color_text('[Degenerate groups]', 'cyan')} modes within "
                    f"{args.degenerate_tol:g} THz of each other -- individual dmu/dQ directions "
                    "are an arbitrary basis choice, the GROUP's summed intensity below is the "
                    "physically meaningful quantity to compare against a single experimental peak:",
                    f_out)
        for group in degenerate_groups:
            mode_ids = ", ".join(str(m[0]) for m in group)
            combined = sum(m[3] for m in group)
            avg_freq = sum(m[1] for m in group) / len(group)
            print_dual(f"  modes {mode_ids} ({avg_freq:.4f} THz / {avg_freq * THZ_TO_CM1:.2f} "
                        f"cm^-1) -> combined intensity {combined:.6f}", f_out)

    print_section('[3] SPECTRUM', f_out)
    freqs_cm1 = [freq_thz * THZ_TO_CM1 for _, freq_thz, _, _ in mode_results]
    intensities = [intensity for _, _, _, intensity in mode_results]
    if args.temperature is not None:
        # Same imaginary-mode-aware skip as Raman's Stage 3 -- Bose-
        # Einstein occupation is undefined for a non-positive frequency.
        n_imaginary_skipped = 0
        weights = []
        for _, freq_thz, _, _ in mode_results:
            if freq_thz <= 0:
                weights.append(1.0)
                n_imaginary_skipped += 1
            else:
                weights.append(bose_einstein_weight(freq_thz, args.temperature))
        intensities = [a * w for a, w in zip(intensities, weights)]
        print_dual(f"Thermal weighting : Bose-Einstein Stokes factor (n+1) at "
                    f"{args.temperature:.1f} K applied to every mode's contribution.", f_out)
        if n_imaginary_skipped:
            print_dual(color_text(
                f"[WARNING] {n_imaginary_skipped} imaginary-frequency mode(s) left unweighted "
                "(Bose-Einstein occupation is undefined for a non-real vibration).", 'yellow'), f_out)
    grid, intensity_grid = build_lorentzian_spectrum(freqs_cm1, intensities, args.linewidth)
    print_dual(f"Grid              : {grid[0]:.1f} - {grid[-1]:.1f} cm^-1, {len(grid)} points, "
                f"{args.linewidth:.2f} cm^-1 linewidth", f_out)

    sim_peaks_all = find_spectrum_peaks(grid, intensity_grid, args.peak_prominence)
    max_intensity = float(np.max(intensity_grid)) if len(intensity_grid) else 0.0
    print_dual(f"Peaks detected    : {len(sim_peaks_all)} (prominence >= "
                f"{args.peak_prominence * 100:.1f}% of max intensity)", f_out)
    if len(sim_peaks_all) and max_intensity > 0:
        peak_intensities = intensity_grid[np.searchsorted(grid, sim_peaks_all)]
        print_dual(f"  {'Position (cm^-1)':<20}{'Rel. intensity (%)'}", f_out)
        for peak_f, peak_i in zip(sim_peaks_all, peak_intensities):
            print_dual(f"  {peak_f:<20.2f}{100.0 * peak_i / max_intensity:.1f}", f_out)

    experimental, experimental_dat_path = None, None
    if args.experimental is not None:
        exp_freq, exp_intensity = read_experimental_spectrum(args.experimental)
        experimental = (exp_freq, exp_intensity)

    if args.save_gnuplot:
        plot_dir = os.path.join(args.directory, "plot")
        os.makedirs(plot_dir, exist_ok=True)
        dat_path = os.path.join(plot_dir, f"{args.output}.dat")
        gplot_path = os.path.join(plot_dir, f"{args.output}.gplot")
        if experimental is not None:
            experimental_dat_path = os.path.join(plot_dir, f"{args.output}_experimental.dat")
        write_ir_spectrum_plot(dat_path, gplot_path, grid, intensity_grid, args.temperature,
                                experimental, experimental_dat_path)
        print_dual(f"{color_text('[Saved]', 'cyan')} {dat_path}, {gplot_path} "
                    f"(cd {plot_dir} && gnuplot {os.path.basename(gplot_path)})", f_out)
    else:
        dat_path = gplot_path = None
        print_dual("Gnuplot data+script : not written (off by default -- pass --save-gnuplot to "
                    "write it, under '<directory>/plot/').", f_out)

    extra_files = ""
    if experimental is not None:
        print_section('[3b] EXPERIMENTAL COMPARISON', f_out)
        print_dual(f"Experimental file : {args.experimental} ({len(exp_freq)} point(s))", f_out)
        exp_peaks = find_spectrum_peaks(exp_freq, exp_intensity, args.peak_prominence)
        print_dual(f"Peaks found       : {len(sim_peaks_all)} simulated, {len(exp_peaks)} experimental", f_out)
        if len(exp_peaks) == 0:
            print_dual(color_text(
                "[WARNING] No experimental peaks found above the prominence threshold -- "
                "try lowering --peak-prominence.", 'yellow'), f_out)
        elif len(sim_peaks_all) == 0:
            print_dual(color_text(
                "[WARNING] No simulated peaks found -- nothing to match against.", 'yellow'), f_out)
        else:
            matches = match_peaks(sim_peaks_all, exp_peaks)
            print_dual(f"  {'Exp. peak (cm^-1)':<20}{'Sim. peak (cm^-1)':<20}{'Delta (cm^-1)'}", f_out)
            for exp_f, sim_f, delta in matches:
                print_dual(f"  {exp_f:<20.2f}{sim_f:<20.2f}{delta:+.2f}", f_out)
            mean_abs_delta = float(np.mean([abs(d) for _, _, d in matches]))
            print_dual(f"Mean |delta|      : {mean_abs_delta:.2f} cm^-1 (average across "
                        f"{len(matches)} matched pair(s))", f_out)
        if experimental_dat_path:
            print_dual(f"{color_text('[Saved]', 'cyan')} {experimental_dat_path}", f_out)
            extra_files = f", {experimental_dat_path}"

    print_section('[3c] MODE VIBRATIONS', f_out)
    phonon_dir = os.path.join(args.directory, "phonon_disp")
    animation_targets = []
    if not os.path.isdir(phonon_dir):
        print_dual(color_text(
            f"[NOTE] '{phonon_dir}' not found -- skipping mode-vibration export/--view-modes. "
            "This folder is only needed here to reload the force constants (to reconstruct each "
            "mode's eigendisplacement pattern); the spectrum above is unaffected.", 'yellow'), f_out)
    else:
        yaml_file = os.path.join(phonon_dir, "phonopy_disp.yaml")
        with open(yaml_file) as f:
            has_embedded_fc = "force_constants" in (yaml.safe_load(f) or {})
        phonon, internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            phonon_dir, label, has_embedded_fc, f_out, symprec=args.symprec)
        # load_phonon_with_force_constants leaves the process chdir'd into
        # phonon_dir -- restore it BEFORE writing anything below via a
        # relative path (args.directory-relative), same "chdir back
        # immediately, before any folder-writing" ordering stb-irModes
        # itself uses. Everything from here on (frame-building, file
        # writes) only touches the already-loaded `phonon` object in
        # memory, no further need to be inside phonon_dir.
        os.chdir(original_dir)
        anim_dir = os.path.join(args.directory, "mode_animations")
        os.makedirs(anim_dir, exist_ok=True)
        print_dual(f"Writing {len(mode_results)} mode animation(s) (extended XYZ, "
                    f"{args.animation_frames} frames each, {stage2_delta} Ang amplitude -- "
                    "same displacement Stage 2 used) under "
                    f"'{anim_dir}':", f_out)
        for mode_index, freq_thz, dmu_dq, intensity in mode_results:
            band_idx = modes[mode_index]["band_index"]
            frames = build_mode_animation_frames(
                phonon, band_idx, stage2_delta, internal_to_angstrom,
                n_frames=args.animation_frames)
            xyz_path = os.path.join(anim_dir, f"mode_{mode_index:02d}_animation.xyz")
            write_mode_animation(frames, xyz_path, fmt='extxyz')
            print_dual(f"  {color_text('[OK]', 'green')} {xyz_path}", f_out)
            animation_targets.append((mode_index, band_idx, freq_thz, frames))

    print_section('[4] SUMMARY & FILES', f_out)
    print_dual(f"Modes analyzed      : {len(mode_results)}/{len(modes)}", f_out)
    if report_path:
        print_dual(f"Report              : {report_path}", f_out)
    if dat_path:
        print_dual(f"Spectrum files      : {dat_path}, {gplot_path}{extra_files}", f_out)
    else:
        print_dual("Spectrum files      : none (pass --save-gnuplot to write the spectrum "
                    "data/script under '<directory>/plot/').", f_out)
    if animation_targets:
        print_dual(f"Mode animations     : {len(animation_targets)} file(s) under "
                    f"'{os.path.join(args.directory, 'mode_animations')}'.", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("IR analysis complete.\n", 'bold'))

    if args.view_modes:
        if not animation_targets:
            print(color_text(
                "\n--view-modes requested but no mode animations were built (see [3c] above).",
                'yellow'))
        else:
            print(color_text(f"\nOpening interactive animation viewer for {len(animation_targets)} "
                              "mode(s) -- close each ASE window to advance to the next.", 'cyan'))
            for mode_index, band_idx, freq_thz, frames in animation_targets:
                print(f"  mode {mode_index:3d} (band {band_idx:3d}, {freq_thz:10.4f} THz / "
                      f"{freq_thz * THZ_TO_CM1:9.3f} cm^-1) -- close window to continue...")
                view_structure_interactive(frames)

    if args.view:
        print(color_text("\n--view: opening the spectrum in matplotlib...", 'cyan'))
        view_ir_spectrum_plot(grid, intensity_grid, args.temperature, experimental)


def write_ir_spectrum_plot(dat_path, gplot_path, grid, intensity, temperature_k=None,
                            experimental=None, experimental_dat_path=None):
    """.dat/.gplot writer for the IR spectrum -- structurally similar to
    raman_analysis.write_spectrum_plot but genuinely simpler (no
    diagonal/full-tensor scope caveat to express, IR has no such
    distinction -- see ir_intensity's own docstring), so kept as its own
    small local function rather than force-generalizing Raman's version.
    """
    scope = "IR spectrum"
    if temperature_k is not None:
        scope += f", {temperature_k:.0f} K Bose weighting"
    with open(dat_path, "w") as f:
        f.write(f"# stb-irAnalysis -- {scope}\n")
        f.write("# columns: 1=frequency(cm^-1) 2=intensity(arb. units)\n")
        for freq, inten in zip(grid, intensity):
            f.write(f"{freq:12.4f} {inten:14.6f}\n")

    if experimental is not None:
        exp_freq, exp_intensity = experimental
        with open(experimental_dat_path, "w") as f:
            f.write("# stb-irAnalysis -- user-supplied experimental IR spectrum\n")
            f.write("# columns: 1=frequency(cm^-1) 2=intensity(arb. units, original scale)\n")
            for freq, inten in zip(exp_freq, exp_intensity):
                f.write(f"{freq:12.4f} {inten:14.6f}\n")

    pdf_name = os.path.splitext(os.path.basename(dat_path))[0] + ".pdf"
    with open(gplot_path, "w") as f:
        lines = [
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-irAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 8,5\n',
            f'set output "{pdf_name}"\n\n',
            f'set title "{scope}"\n',
            'set xlabel "Wavenumber (cm^{-1})"\n',
            'set ylabel "Simulated intensity (arb. units)"\n',
            'set grid\n',
        ]
        if experimental is not None:
            lines += [
                'set y2label "Experimental intensity (arb. units)"\n',
                'set y2tics\n',
                'set key top right\n',
                f'plot "{os.path.basename(dat_path)}" using 1:2 with lines lw 2 lc rgb "#2255cc" '
                'axes x1y1 title "Simulated", \\\n'
                f'     "{os.path.basename(experimental_dat_path)}" using 1:2 with lines lw 2 '
                'lc rgb "#cc5522" axes x1y2 title "Experimental"\n',
            ]
        else:
            lines += [
                'unset key\n',
                f'plot "{os.path.basename(dat_path)}" using 1:2 with lines lw 2 lc rgb "#2255cc"\n',
            ]
        f.writelines(lines)


def view_ir_spectrum_plot(grid, intensity, temperature_k=None, experimental=None):
    """Interactive matplotlib preview of the IR spectrum -- the on-screen
    counterpart to write_ir_spectrum_plot's saved gnuplot .dat/.gplot pair,
    same data, never written to disk (see CLAUDE.md: WORKFLOW_TOOLS writes
    gnuplot pairs, matplotlib here is preview-only, matching
    raman_analysis.view_spectrum_plot). Experimental overlay on a separate
    right-hand y-axis (ax.twinx()), same rationale as the gnuplot y2 axis
    in write_ir_spectrum_plot -- simulated intensity and real experimental
    intensity are on incomparable scales. matplotlib is imported lazily,
    only when --view was actually passed.
    """
    import matplotlib.pyplot as plt
    scope = "IR spectrum"
    if temperature_k is not None:
        scope += f", {temperature_k:.0f} K Bose weighting"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(grid, intensity, color='#2255cc', lw=2, label="Simulated")
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Simulated intensity (arb. units)")
    ax.set_title(scope)
    ax.grid(True, alpha=0.3)
    if experimental is not None:
        exp_freq, exp_intensity = experimental
        ax2 = ax.twinx()
        ax2.plot(exp_freq, exp_intensity, color='#cc5522', lw=2, label="Experimental")
        ax2.set_ylabel("Experimental intensity (arb. units)")
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
