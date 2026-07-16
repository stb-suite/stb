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
import shutil
import argparse
from datetime import datetime
import yaml
from phonopy.interface.siesta import write_siesta
from stb.core.cli import color_text, show_intro
from stb.core.calc_directives import force_single_point
from stb.core.pseudopotentials import get_required_pseudos, resolve_pseudo_source
from stb.core.phonon_workflow import (
    detect_system_label, load_phonon_with_force_constants, get_gamma_modes, displace_along_mode,
)

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


def print_dual(text, file_handle=None):
    """Prints to stdout with color, writes to file without color. Same
    duplicated-per-tool helper as phonons_create.py/phonons_pos.py."""
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


def build_optical_block(mesh, broaden_ev, axis_vec, nbands=None):
    """The %block Optical.Mesh/Optical.Vector + OpticalCalculation T fdf
    stanza for one Optical.Vector direction -- SIESTA's interband/RPA
    dielectric-function machinery, verified via SIESTA's own official
    tutorial to work for fully 3D periodic bulk crystals (not just
    vacuum-padded systems, unlike %block ExternalElectricField).
    """
    lines = [
        "OpticalCalculation T",
        f"Optical.Broaden        {broaden_ev} eV",
        "%block Optical.Mesh",
        f"  {mesh[0]}  {mesh[1]}  {mesh[2]}",
        "%endblock Optical.Mesh",
        "%block Optical.Vector",
        f"  {axis_vec[0]:.4f}  {axis_vec[1]:.4f}  {axis_vec[2]:.4f}",
        "%endblock Optical.Vector",
    ]
    if nbands:
        lines.append(f"Optical.NumberOfBands  {nbands}")
    return "\n".join(lines) + "\n"


def write_optical_folder(out_dir, displaced_atoms, structure_filename, calc_text, pseudos):
    os.makedirs(out_dir, exist_ok=True)
    write_siesta(os.path.join(out_dir, os.path.basename(structure_filename)), displaced_atoms)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    for pseudo_path in pseudos:
        shutil.copy(pseudo_path, os.path.join(out_dir, os.path.basename(pseudo_path)))


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
    parser.add_argument("--displacement", type=float, default=0.02,
                        help="Finite-difference displacement (Ang) along each mode's "
                             "eigendisplacement for the Raman-tensor derivative (default: 0.02).")
    parser.add_argument("--full-tensor", action="store_true",
                        help="Compute the FULL symmetric Raman tensor (Rxx, Ryy, Rzz, Rxy, Rxz, "
                             "Ryz) instead of just the diagonal -- adds 3 more Optical.Vector "
                             "directions per sign (the face-diagonals), doubling the folder count "
                             "per mode (12 instead of 6). Off-diagonal components are recovered "
                             "by stb-ramanAnalysis from these mixed-direction measurements.")
    parser.add_argument("--optical-mesh", type=int, nargs=3, default=[10, 10, 10],
                        help="Optical.Mesh k-grid for the dielectric-function calculation "
                             "(default: 10 10 10).")
    parser.add_argument("--optical-nbands", type=int, default=None,
                        help="Optical.NumberOfBands (default: SIESTA's own default -- all "
                             "bands within the basis).")
    parser.add_argument("--optical-broaden", type=float, default=0.2,
                        help="Optical.Broaden in eV (default: 0.2).")
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
            frequencies, mode_band_indices = get_gamma_modes(phonon, exclude_acoustic=True)
            unique_elements = list(set(phonon.primitive.symbols))
        finally:
            os.chdir(original_dir)

        n_imaginary = int((frequencies < 0).sum())
        print_dual(f"Non-acoustic Gamma modes : {len(frequencies)}", f_out)
        if n_imaginary:
            print_dual(color_text(
                f"[WARNING] {n_imaginary} mode(s) have imaginary (negative) frequency -- the "
                "structure/supercell may not be at a real energy minimum. Their Raman tensor "
                "is not physically meaningful.", 'yellow'), f_out)
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            flag = color_text(" [IMAGINARY]", 'red') if freq < 0 else ""
            print_dual(f"  mode {k:3d} (band {int(band_idx):3d}) : {freq:10.4f} THz{flag}", f_out)

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

        selected = []
        for k, (freq, band_idx) in enumerate(zip(frequencies, mode_band_indices), start=1):
            if args.modes is not None and k not in args.modes:
                continue
            if args.freq_min is not None and freq < args.freq_min:
                continue
            if args.freq_max is not None and freq > args.freq_max:
                continue
            selected.append((k, freq, int(band_idx)))

        if not selected:
            print_dual(color_text(
                "\n[ERROR] No modes selected after applying --modes/--freq-min/--freq-max "
                "-- nothing to do.", 'red'), f_out)
            sys.exit(1)

        print_dual(f"\n{color_text('[3] OPTICAL DISPLACEMENT FOLDERS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        n_folders = len(selected) * len(_SIGNS) * len(axes)
        print_dual(f"Selected modes    : {len(selected)}", f_out)
        print_dual(f"Folders to write  : {len(selected)} modes x {len(_SIGNS)} signs x "
                    f"{len(axes)} axes = {n_folders} independent single-point SIESTA runs.", f_out)

        report_rows = []  # (label, mode_index, band_index, frequency, sign, axis, dir)
        structure_filename = "structure.fdf"
        for k, freq, band_idx in selected:
            for sign_name, sign_val in _SIGNS:
                displaced = displace_along_mode(
                    phonon, band_idx, args.displacement, internal_to_angstrom, sign=sign_val)
                for axis_name, axis_vec in axes:
                    label = f"mode_{k:02d}_{sign_name}_{axis_name}"
                    mode_dir = os.path.join(optical_root, label)
                    optical_block = build_optical_block(
                        args.optical_mesh, args.optical_broaden, axis_vec, args.optical_nbands)
                    calc_text = force_single_point(optical_calc_template) + "\n" + optical_block
                    write_optical_folder(mode_dir, displaced, structure_filename, calc_text, pseudos)
                    report_rows.append((label, k, band_idx, freq, sign_name, axis_name, mode_dir))
                    print_dual(f"  {color_text('[OK]', 'green')} {mode_dir}", f_out)

        print_dual(f"\n{color_text('[4] SUMMARY & NEXT STEPS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"{len(report_rows)} folder(s) written under '{optical_root}'.", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual(f"  1. Run SIESTA in every '{optical_root}/mode_*/' folder.", f_out)
        print_dual(f"  2. Once they're done, run: stb-ramanAnalysis --directory {output_root}", f_out)

        f_out.write("\n# MODE_TABLE -- parsed by stb-ramanAnalysis, do not reorder the "
                     "first 6 columns\n")
        f_out.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
                     f"{'sign':<8}{'axis':<6}{'dir'}\n")
        for label, k, band_idx, freq, sign_name, axis_name, mode_dir in report_rows:
            f_out.write(f"{label:<24}{k:<12}{band_idx:<12}{freq:<16.6f}{sign_name:<8}{axis_name:<6}{mode_dir}\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Optical displacement folders ready for Stage 3 (stb-ramanAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
