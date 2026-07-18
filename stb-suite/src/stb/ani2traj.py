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
import argparse
from datetime import datetime
import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from ase import Atoms
from ase.io import write as ase_write

from stb.core import structure_io, siesta_log
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_ani2traj_report.txt"

# --out-format -> (ASE format name, default file extension). All three are
# ASE-native writers for a LIST of Atoms (a real multi-frame trajectory), and
# all three carry a full per-frame lattice -- unlike SIESTA's own .ANI file,
# which is a bare multi-frame XYZ with no cell information at all.
#   xsf: multi-frame "AXSF" (ANIMSTEPS + one PRIMVEC/PRIMCOORD block per
#        frame) -- read natively by both OVITO and VMD's xsf plugin.
#   pdb: multi-model PDB (CRYST1 + one MODEL/ENDMDL block per frame) -- VMD's
#        own default/preferred format; also read by OVITO.
#   xyz: extended XYZ (a Lattice="..." tag per frame) -- OVITO reads this
#        natively; VMD's plain xyz reader does NOT apply the Lattice tag, so
#        positions show but without periodic boundaries there. The ONLY one
#        of the three that can also carry the per-frame E_KS/Temp_ion/Time
#        enrichment below (ASE's extxyz writer serializes atoms.info scalars
#        into the header line; its xsf/pdb writers have no such slot).
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}

_MD_TIMESTEP_RE = re.compile(r"^\s*MD\.LengthTimeStep\s+([\d.eE+-]+)\s*(\w*)", re.IGNORECASE | re.MULTILINE)
_MD_INITIAL_STEP_RE = re.compile(r"^\s*MD\.InitialTimeStep\s+(\d+)", re.IGNORECASE | re.MULTILINE)
_TIME_UNIT_TO_FS = {"fs": 1.0, "ps": 1.0e3, "ns": 1.0e6, "": 1.0}


def read_static_lattice(label, f_out=None):
    """Returns the (3, 3) cell matrix (Angstrom) for this SIESTA run, read
    from <label>.XV (preferred -- SIESTA's own output geometry) or
    <label>.fdf (fallback -- the input geometry). Used as-is for every
    frame when a richer per-frame lattice (read_frame_lattices, from the
    .out log) isn't available. Returns None if neither is readable."""
    xv_file = f"{label}.XV"
    fdf_file = f"{label}.fdf"

    if os.path.exists(xv_file):
        try:
            geom = sisl.get_sile(xv_file).read_geometry()
            print_dual(f"[INFO] Lattice read from {xv_file}.", f_out)
            return np.array(geom.cell)
        except Exception as e:
            print_dual(color_text(
                f"[WARNING] Found '{xv_file}' but could not read it: {e}", 'yellow'), f_out)

    if os.path.exists(fdf_file):
        try:
            structure = structure_io.read_fdf(fdf_file)
            print_dual(f"[INFO] Lattice read from {fdf_file}.", f_out)
            return np.array(structure.lattice)
        except Exception as e:
            print_dual(color_text(
                f"[WARNING] Found '{fdf_file}' but could not read it: {e}", 'yellow'), f_out)

    return None


def read_frame_lattices(label, nframes, f_out=None):
    """Returns (cells, steps): `cells` is a list of `nframes` (3, 3) cell
    matrices, one per .ANI frame (BEFORE --stride is applied -- the caller
    strides both lists together so they stay aligned); `steps` is the raw
    per-MD-step data from siesta_log.get_md_trajectory (for the E_KS/
    Temp_ion enrichment below), or None if a per-frame lattice wasn't used.

    Prefers <label>.out's per-MD-step 'outcell:' blocks -- correct even for
    variable-cell (NPT/Parrinello-Rahman) runs, and simply constant frame to
    frame for a fixed-cell one, so there's no reason to prefer the static
    fallback when this is available. Falls back to a single static lattice
    (read_static_lattice, repeated for every frame) when '.out' is missing,
    has no MD steps at all (e.g. it's from an unrelated run), or doesn't
    have cell data for every one of the `nframes` frames.
    """
    out_file = f"{label}.out"
    if os.path.exists(out_file):
        steps = siesta_log.get_md_trajectory(out_file)
        cells = [s['cell'] for s in steps]
        if len(cells) >= nframes and all(c is not None for c in cells[:nframes]):
            print_dual(f"[INFO] Lattice read per-frame from {out_file} "
                       f"({len(cells)} MD step(s) found).", f_out)
            return cells[:nframes], steps[:nframes]
        elif steps:
            print_dual(color_text(
                f"[WARNING] '{out_file}' has MD step data but not a usable cell for "
                f"every one of the {nframes} frame(s) -- falling back to a single "
                f"static lattice.", 'yellow'), f_out)

    static_cell = read_static_lattice(label, f_out)
    if static_cell is None:
        return None, None
    return [static_cell] * nframes, None


def read_md_timestep_fs(label):
    """Returns (initial_step, dt_fs) parsed from <label>.fdf's
    MD.InitialTimeStep/MD.LengthTimeStep -- used to give each frame its real
    simulation time (fs) instead of a bare frame index, for correct-speed
    playback in OVITO. (1, None) if <label>.fdf is missing or doesn't have
    MD.LengthTimeStep (1 is SIESTA's own default MD.InitialTimeStep)."""
    fdf_file = f"{label}.fdf"
    if not os.path.exists(fdf_file):
        return 1, None
    try:
        with open(fdf_file) as f:
            text = f.read()
    except Exception:
        return 1, None

    initial = 1
    m = _MD_INITIAL_STEP_RE.search(text)
    if m:
        initial = int(m.group(1))

    m = _MD_TIMESTEP_RE.search(text)
    if not m:
        return initial, None
    value = float(m.group(1))
    unit = m.group(2).lower()
    return initial, value * _TIME_UNIT_TO_FS.get(unit, 1.0)


def unwrap_trajectory(positions, cells):
    """Returns `positions` (list of (natoms, 3) Cartesian arrays, one per
    frame) with each atom's periodic-boundary "jumps" removed, so a molecule
    that legitimately drifts across the box edge during the AIMD run reads
    as continuous motion instead of visually teleporting to the opposite
    side. For each frame after the first, the raw Cartesian displacement
    from the previous frame is converted to fractional coordinates (using
    that frame's own cell -- matters for a variable-cell run), wrapped into
    [-0.5, 0.5) (the minimum-image convention -- same assumption
    Atoms.get_all_distances(mic=True) makes elsewhere in this suite: no atom
    moves more than half a cell length in a single MD step, true for any
    physically reasonable timestep), converted back to Cartesian, and
    accumulated onto the previous (already-unwrapped) position.
    """
    unwrapped = [positions[0].copy()]
    for t in range(1, len(positions)):
        inv_cell = np.linalg.inv(cells[t])
        disp_frac = (positions[t] - positions[t - 1]) @ inv_cell
        disp_frac -= np.round(disp_frac)
        unwrapped.append(unwrapped[-1] + disp_frac @ cells[t])
    return unwrapped


def main():
    parser = argparse.ArgumentParser(
        description="Convert a SIESTA AIMD trajectory (<label>.ANI) into a multi-frame "
                     "trajectory format that also carries the lattice -- for opening in "
                     "OVITO, VMD, etc. .ANI itself is a bare multi-frame XYZ with no cell "
                     "information at all; the lattice is read per-frame from <label>.out "
                     "when available (correct for variable-cell/NPT runs too), or a single "
                     "static one from <label>.XV/<label>.fdf otherwise.",
        epilog="Example usage:\n"
               "  stb-ani2traj --label aimd --out-format xsf\n"
               "  stb-ani2traj -l siesta -of pdb -o traj.pdb --stride 5\n"
               "  stb-ani2traj -l siesta -of xyz --unwrap",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-l", "--label", required=True, help="SystemLabel (required).")
    parser.add_argument(
        "-of", "--out-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
        help="Output trajectory format (default: xsf):\n"
             "  xsf - multi-frame AXSF, read natively by OVITO and VMD's xsf plugin.\n"
             "  pdb - multi-model PDB (CRYST1), VMD's own default format, also read by OVITO.\n"
             "  xyz - extended XYZ (Lattice tag), OVITO-native; VMD's xyz reader ignores\n"
             "        the Lattice tag (positions show, but without periodic boundaries).\n"
             "        The only format that also carries per-frame E_KS/Temp_ion/Time."
    )
    parser.add_argument("-o", "--output", default=None,
                        help="Output filename. Default: <label>_traj.<ext>")
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep only every Nth frame (default: 1, every frame).")
    parser.add_argument(
        "--unwrap", action="store_true",
        help="Remove periodic-boundary 'jumps': each frame's positions are corrected by "
             "the minimum-image displacement relative to the previous frame, so atoms "
             "that legitimately cross the box edge read as continuous motion in OVITO/VMD "
             "instead of teleporting. Off by default (raw .ANI positions, as SIESTA wrote "
             "them, mod the cell)."
    )
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-ani2traj {VERSION}")

    args = parser.parse_args()

    ase_format, default_ext = OUTPUT_FORMATS[args.out_format]
    ani_file = f"{args.label}.ANI"
    output_file = args.output if args.output else f"{args.label}_traj{default_ext}"

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - AIMD Trajectory Converter",
            "Adds back the lattice SIESTA's own .ANI file doesn't carry",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("AIMD TRAJECTORY CONVERTER:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-ANI2TRAJ REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {ani_file}", f_out)
    print_dual(f"Output file       : {output_file} ({args.out_format})", f_out)
    print_dual(f"Stride            : {args.stride}", f_out)
    print_dual(f"Unwrap PBC        : {'yes' if args.unwrap else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING TRAJECTORY", f_out)

    if not os.path.exists(ani_file):
        fail(f"Input file '{ani_file}' not found.")
    if args.stride < 1:
        fail(f"--stride must be >= 1, got {args.stride}.")

    try:
        sile = sisl.get_sile(ani_file)
        all_frames = sile.read_geometry[:]()
    except Exception as e:
        fail(f"Could not read '{ani_file}': {e}")

    if not all_frames:
        fail(f"'{ani_file}' contains no frames.")

    all_cells, all_steps = read_frame_lattices(args.label, len(all_frames), f_out)
    if all_cells is None:
        fail(f"Could not find a readable '{args.label}.out', '{args.label}.XV' or "
             f"'{args.label}.fdf' -- a lattice is required ('{ani_file}' doesn't carry "
             f"one of its own).")

    # Apply --stride identically to the frames and their (possibly per-frame)
    # cells/MD-step data, so index i keeps meaning "the same MD step" in both.
    frames = all_frames[::args.stride]
    cells = all_cells[::args.stride]
    steps = all_steps[::args.stride] if all_steps is not None else None

    print_dual(f"[OK] Read {len(frames)} frame(s) from {ani_file} (stride {args.stride}).", f_out)
    print_dual(f"Atoms per frame   : {len(frames[0].atoms)}", f_out)
    if steps is not None and any(not np.allclose(cells[0], c) for c in cells[1:]):
        print_dual("[INFO] Cell varies frame to frame (variable-cell/NPT run).", f_out)
    print_dual("Lattice (Ang, first frame):", f_out)
    for row in cells[0]:
        print_dual(f"    [{row[0]:9.4f}, {row[1]:9.4f}, {row[2]:9.4f}]", f_out)

    positions = [g.xyz for g in frames]
    if args.unwrap:
        positions = unwrap_trajectory(positions, cells)
        print_dual("[OK] Unwrapped periodic-boundary jumps relative to the previous frame.", f_out)

    embed_properties = steps is not None and args.out_format == 'xyz'
    if steps is not None and args.out_format != 'xyz':
        print_dual(color_text(
            "[INFO] Per-frame energy/temperature (from the .out log) is only embedded "
            "for --out-format xyz -- skipping it for this format.", 'yellow'), f_out)

    initial_step, dt_fs = (None, None)
    if args.out_format == 'xyz':
        initial_step, dt_fs = read_md_timestep_fs(args.label)
        if dt_fs is not None:
            print_dual(f"[INFO] MD.LengthTimeStep found -- embedding real simulation time (fs) per frame.", f_out)

    ase_frames = []
    for i, g in enumerate(frames):
        atoms = Atoms(symbols=[a.symbol for a in g.atoms], positions=positions[i],
                       cell=cells[i], pbc=True)
        if embed_properties:
            step = steps[i]
            if step['E_KS'] is not None:
                atoms.info['E_KS'] = step['E_KS']
            if step['Temp_ion'] is not None:
                atoms.info['Temp_ion'] = step['Temp_ion']
        if dt_fs is not None:
            step_number = initial_step + i * args.stride
            atoms.info['Time'] = step_number * dt_fs
        ase_frames.append(atoms)

    print_section("[2] WRITING OUTPUT", f_out)
    try:
        ase_write(output_file, ase_frames, format=ase_format)
        print_dual(f"[OK] Saved {output_file} ({len(ase_frames)} frame(s), {args.out_format} format)", f_out)
    except Exception as e:
        fail(f"Failed to write output file: {e}")

    print_section("[3] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Output            : {output_file}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Trajectory converted -- now go make a nice movie of it.\n", 'bold'))


if __name__ == "__main__":
    main()
