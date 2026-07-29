"""AIMD trajectory helpers: per-frame lattice, real MD timestep, and PBC
unwrap. Extracted from ani2traj.py once stb-aimdAnalysis became a second
consumer needing the exact same per-frame cell/timestep/unwrap logic (same
"extract on second use" policy as structure_io.py/siesta_log.py/symmetry.py).
"""

import os
import re

import numpy as np

from stb.core.deps import require_sisl
from stb.core import structure_io, siesta_log
from stb.core.cli import color_text, print_dual

_MD_TIMESTEP_RE = re.compile(r"^\s*MD\.LengthTimeStep\s+([\d.eE+-]+)\s*(\w*)", re.IGNORECASE | re.MULTILINE)
_MD_INITIAL_STEP_RE = re.compile(r"^\s*MD\.InitialTimeStep\s+(\d+)", re.IGNORECASE | re.MULTILINE)
_TIME_UNIT_TO_FS = {"fs": 1.0, "ps": 1.0e3, "ns": 1.0e6, "": 1.0}


def read_static_lattice(label, f_out=None, fdf_path=None):
    """Returns the (3, 3) cell matrix (Angstrom) for this SIESTA run, read
    from <label>.XV (preferred -- SIESTA's own output geometry, always
    named after SystemLabel regardless of how the input file itself was
    named) or a .fdf (fallback -- the input geometry). Used as-is for
    every frame when a richer per-frame lattice (read_frame_lattices, from
    the .out log) isn't available. Returns None if neither is readable.

    `fdf_path`, if given, overrides the auto-derived <label>.fdf guess for
    the fallback -- unlike .XV/.ANI/.HSX/.WFSX (always literally named
    after SystemLabel by SIESTA itself), the INPUT .fdf's filename is
    chosen by the user and is very often NOT <label>.fdf (e.g. SystemLabel
    'siesta' with the real input file called calc.fdf) -- the same gap
    already closed elsewhere in this suite via --geometry-file (stb-sts/
    stb-coop/stb-ipr/stb-effmass/stb-spintexture).
    """
    sisl = require_sisl()
    xv_file = f"{label}.XV"
    fdf_file = fdf_path if fdf_path is not None else f"{label}.fdf"

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


def read_frame_lattices(label, nframes, f_out=None, fdf_path=None):
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

    `fdf_path`, if given, is passed straight through to read_static_lattice's
    own fallback -- see its docstring for why the INPUT .fdf's filename
    can't be reliably auto-derived from `label` the way .out/.XV/.ANI can.
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

    static_cell = read_static_lattice(label, f_out, fdf_path=fdf_path)
    if static_cell is None:
        return None, None
    return [static_cell] * nframes, None


def read_md_timestep_fs(label, fdf_path=None):
    """Returns (initial_step, dt_fs) parsed from a .fdf's
    MD.InitialTimeStep/MD.LengthTimeStep -- used to give each frame its real
    simulation time (fs) instead of a bare frame index, for correct-speed
    playback in OVITO. (1, None) if the .fdf is missing or doesn't have
    MD.LengthTimeStep (1 is SIESTA's own default MD.InitialTimeStep).

    `fdf_path`, if given, overrides the auto-derived <label>.fdf guess --
    see read_static_lattice's docstring for why this can't be reliably
    auto-derived from `label` alone in general.
    """
    fdf_file = fdf_path if fdf_path is not None else f"{label}.fdf"
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
