#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import io
import re
import sys
import glob
import argparse
import contextlib
from datetime import datetime
import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections
from stb.core.kspace import detect_vacuum_axes, dimensionality_label
from stb.core import siesta_log

# Purpose: write one of SIESTA's real-space grid files (<label>.RHO, .VT, ...) as a Gaussian cube file for VESTA,
# VMD, XCrySDen, Multiwfn or a Bader code -- and say, on the way, what the grid is: its quantity and units, the mesh
# and the cutoff behind it, the geometry it belongs to, its values (integral against what it must integrate to,
# extremes and where they are, isosurface levels), the spin channels, the dipole of a molecule or slab, and a
# read-back check of the cube that was written.
#
# Units. SIESTA writes densities in e/Bohr^3 and potentials in Ry, on a Bohr cell. sisl's reader scales them to
# e/Ang^3 and eV (`sile.grid_unit`). A cube file has Bohr lengths, so a density goes back to e/Bohr^3 (dividing by
# grid_unit, exactly 1/0.529177^3) -- before that correction every RHO cube paired Bohr coordinates with e/Ang^3
# values, 6.75x too large for any program that integrates it. A potential has no length in its unit; it is written in
# eV unless --potential-unit says Ry or Ha.

REPORT_FILE = "stb_cube_report.txt"
BOHR = 0.529177210903            # Ang
DEBYE_PER_E_ANG = 1.0 / 0.20819434
POTENTIAL_UNITS = {"eV": 1.0, "Ry": 13.605693122994, "Ha": 27.211386245988}

# kind 'density' -> written in e/Bohr^3; 'potential' -> --potential-unit. `expect`: what the grid must integrate to.
# `signed`: the values change sign by construction (no 'negative density' warning; +-isovalues).
GRID_TYPES = {
    "RHO":     dict(kind="density", flag="SaveRho", expect="electrons", signed=False,
                    what="valence electron density"),
    "DRHO":    dict(kind="density", flag="SaveDeltaRho", expect="zero", signed=True,
                    what="deformation density: valence density minus the sum of free-atom densities"),
    "BADER":   dict(kind="density", flag="SaveBaderCharge", expect="all-electron", signed=False,
                    what="valence density plus model core charges (input for a Bader analysis)"),
    "TOCH":    dict(kind="density", flag="SaveTotalCharge", expect="neutral", signed=True,
                    what="total charge: valence density plus the ionic (pseudo-core) charge, ions negative"),
    "IOCH":    dict(kind="density", flag="SaveIonicCharge", expect="ions", signed=True,
                    what="ionic (pseudo-core) charge, negative in SIESTA's sign convention"),
    "RHOXC":   dict(kind="density", flag="SaveRhoXC", expect=None, signed=False,
                    what="density seen by exchange-correlation (valence + partial core correction)"),
    "RHOINIT": dict(kind="density", flag="SaveInitialChargeDensity", expect="free-atoms", signed=False,
                    what="initial density: the sum of free-atom densities"),
    "LDOS":    dict(kind="density", flag="%block LocalDensityOfStates", expect=None, signed=False,
                    what="local density of states in the LocalDensityOfStates energy window"),
    "VT":      dict(kind="potential", flag="SaveTotalPotential", expect=None, signed=True,
                    what="total Kohn-Sham potential (Hartree + local pseudopotential + exchange-correlation)"),
    "VH":      dict(kind="potential", flag="SaveElectrostaticPotential", expect=None, signed=True,
                    what="electrostatic potential (Hartree + local part of the pseudopotentials)"),
    "VNA":     dict(kind="potential", flag="SaveNeutralAtomPotential", expect=None, signed=True,
                    what="neutral-atom potential (sum of the free atoms' electrostatic potentials)"),
}

SPIN_CHOICES = ["total", "up", "down", "diff"]
VACUUM_GAP = 4.0                 # Ang of empty space along an axis that make it a vacuum (non-periodic) axis
ISO_FRACTIONS = (0.5, 0.9, 0.99)
CUBE_FMT = "%13.5e"              # 6 significant digits, the usual cube precision
MAX_ATOM_ROWS = 12


# ------------------------------------------------------------------------------------------------- small helpers --
def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0


def fmt_num(x, digits=4):
    """A number with `digits` significant digits, fixed notation when that stays readable."""
    if x == 0 or not np.isfinite(x):
        return f"{x:.{digits}f}" if np.isfinite(x) else str(x)
    if 1e-3 <= abs(x) < 1e5:
        decimals = max(0, digits - 1 - int(np.floor(np.log10(abs(x)))))
        return f"{x:.{decimals}f}"
    return f"{x:.{digits - 1}e}"


_NEG_ZERO = re.compile(r"(?<![\w.])-(0\.0+)(?![\d])")


def nz(text):
    """'-0.0000' -> '0.0000' (a rounded tiny negative number reads as a sign that means nothing)."""
    return _NEG_ZERO.sub(r"\1", text)


def cell_parameters(cell):
    lengths = np.linalg.norm(cell, axis=1)

    def angle(i, j):
        c = np.dot(cell[i], cell[j]) / (lengths[i] * lengths[j])
        return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    return lengths, (angle(1, 2), angle(0, 2), angle(0, 1)), abs(np.linalg.det(cell))


def effective_mesh_cutoff(cell_ang, shape):
    """SIESTA's 'Mesh cutoff (used)' in Ry: min_i (n_i |b_i| / 2)^2, b_i the reciprocal vectors (2 pi included) in
    1/Bohr. Checked against InitMesh of a real Si run: 248.240 Ry from the grid, 248.240 Ry in the log."""
    b = 2 * np.pi * np.linalg.inv(cell_ang / BOHR).T
    per_axis = (np.asarray(shape) * np.linalg.norm(b, axis=1) / 2.0) ** 2
    return float(per_axis.min()), per_axis


def plane_spacing(cell_ang, shape):
    """Distance between neighbouring grid planes along each axis (Ang): the real resolution of the mesh."""
    recip = np.linalg.inv(cell_ang).T
    return 1.0 / np.linalg.norm(recip, axis=1) / np.asarray(shape)


def read_log_info(path):
    """What the SIESTA log says that bears on a grid. Every field None when not found."""
    info = dict(path=path, version=None, spin=None, electrons=None, mesh=None, cutoff_required=None,
                cutoff_used=None, dipole_debye=None, spin_moment=None, finished=False)
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if s.startswith("Version") and ":" in s and info["version"] is None:
                    info["version"] = s.split(":", 1)[1].strip()
                elif "redata: Spin configuration" in s and "=" in s:
                    info["spin"] = s.split("=", 1)[1].strip()
                elif s.startswith("Total number of electrons:"):
                    try:
                        info["electrons"] = float(s.split(":", 1)[1].split()[0])
                    except (IndexError, ValueError):
                        pass
                elif s.startswith("InitMesh: MESH =") or s.startswith("InitMesh: MESH="):
                    try:
                        dims = s.split("=", 1)[1].split("=")[0].replace("x", " ").split()
                        info["mesh"] = tuple(int(v) for v in dims[:3])
                    except (IndexError, ValueError):
                        pass
                elif "Mesh cutoff (required, used)" in s:
                    try:
                        vals = s.split("=", 1)[1].split()
                        info["cutoff_required"], info["cutoff_used"] = float(vals[0]), float(vals[1])
                    except (IndexError, ValueError):
                        pass
                elif "Electric dipole (Debye)" in s and "=" in s:
                    try:
                        info["dipole_debye"] = np.array([float(v) for v in s.split("=", 1)[1].split()[:3]])
                    except ValueError:
                        pass
                elif "End of run" in s or s.startswith("Job completed"):
                    info["finished"] = True
    except OSError:
        return None
    info["spin_moment"] = siesta_log.get_spin_moment(path)
    return info


def spin_kind_from_log(text):
    t = (text or "").lower()
    if not t:
        return None
    if "orbit" in t:
        return "spin-orbit"
    if "non" in t and "lin" in t:
        return "non-collinear"
    if "colinear" in t or "collinear" in t or "polarized" in t:
        return "polarized"
    if t in ("none", "unpolarized", "f", "false"):
        return "unpolarized"
    return None


def detect_spin_configuration(label):
    """Spin configuration from <label>.HSX (sisl reads the Hamiltonian; used only when no SIESTA log says it).
    Returns (kind or None, message or None). Verified on a real spin-polarized O2 run: Spin{polarized}."""
    hsx_file = f"{label}.HSX"
    if not os.path.exists(hsx_file):
        return None, None
    try:
        spin = sisl.get_sile(hsx_file).read_hamiltonian().spin
    except Exception as e:
        return None, f"found '{hsx_file}' but could not read its spin configuration: {e}"
    kind = ("spin-orbit" if spin.is_spinorbit else "non-collinear" if spin.is_noncolinear else
            "polarized" if spin.is_polarized else "unpolarized")
    return kind, None


def channel_weights(kind, nspin, channel):
    """sisl read_grid weights for one spin channel. A potential's 'total' is the AVERAGE of up and down (summing two
    potentials doubles them); a density's is the sum."""
    if nspin == 1 or channel == "total":
        if nspin == 1:
            return [1.0]
        w = 0.5 if kind == "potential" else 1.0
        return [w, w]
    return {"up": [1.0, 0.0], "down": [0.0, 1.0], "diff": [1.0, -1.0]}[channel]


def nearest_atom(point, geom, cell):
    """(index, distance in Ang) of the atom nearest to `point`, minimum image."""
    if geom is None or len(geom.xyz) == 0:
        return None, None
    inv = np.linalg.inv(cell)
    d = (geom.xyz - point) @ inv
    d -= np.round(d)
    dist = np.linalg.norm(d @ cell, axis=1)
    i = int(np.argmin(dist))
    return i, float(dist[i])


def iso_levels(values, dv, fractions, signed):
    """Isovalue enclosing each fraction of the integral of |values| (descending sort). For a signed grid the level is
    a +- pair on |values|."""
    v = np.abs(values.ravel()) if signed else values.ravel()
    v = np.sort(v[v > 0])[::-1]
    if v.size == 0:
        return []
    cum = np.cumsum(v) * dv
    total = cum[-1]
    out = []
    for frac in fractions:
        k = int(np.searchsorted(cum, frac * total))
        out.append((frac, float(v[min(k, v.size - 1)])))
    return out


def write_cube(path, data, cell_ang, numbers, positions_ang, comment1, comment2):
    """Gaussian cube: two comment lines (the second one Gaussian's loop-order line); natoms + origin; per axis n_i and the voxel vector (Bohr); one line per
    atom (Z, charge = Z, x y z in Bohr); values with x slowest and z fastest, every z-row starting a new line in
    lines of six (Gaussian's own layout, which every reader accepts)."""
    nx, ny, nz = data.shape
    with open(path, "w") as f:
        f.write(comment1.replace("\n", " ") + "\n")
        f.write(comment2.replace("\n", " ") + "\n")
        f.write(f"{len(numbers):5d} {0.0:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        for n, vec in zip((nx, ny, nz), cell_ang / np.array([[nx], [ny], [nz]]) / BOHR):
            f.write(f"{n:5d} {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
        for z, pos in zip(numbers, positions_ang / BOHR):
            f.write(f"{int(z):5d} {float(z):12.6f} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        full, rest = divmod(nz, 6)
        row_fmt = (" ".join([CUBE_FMT] * 6) + "\n") * full
        if rest:
            row_fmt += " ".join([CUBE_FMT] * rest) + "\n"
        rows = np.ascontiguousarray(data, dtype=float).reshape(nx * ny, nz)
        for row in rows:
            f.write(row_fmt % tuple(row))


def read_back_cube(path):
    """Reads a cube with ASE's reader (independent of the writer above)."""
    from ase.io.cube import read_cube_data
    data, atoms = read_cube_data(path)
    return data, atoms


def find_labels(grid_type):
    return sorted(f[:-(len(grid_type) + 1)] for f in glob.glob(f"*.{grid_type}"))


# --------------------------------------------------------------------------------------------------------- main --
def main():
    type_lines = "\n".join(f"  {k:<8}{v['what']} ({v['flag']})" for k, v in GRID_TYPES.items())
    parser = argparse.ArgumentParser(
        description="Convert a SIESTA real-space grid file (<label>.RHO, .VT, .VH, ...) to a Gaussian cube file\n"
                    "(VESTA, VMD, XCrySDen, Multiwfn, Bader codes), with a report on the grid: quantity and units,\n"
                    "mesh and effective cutoff, geometry, spin channels, values (integral vs what it must be,\n"
                    "extremes, isosurface levels), the dipole of a molecule or slab, and a read-back check of the\n"
                    "cube. Lengths are written in Bohr, densities in e/Bohr^3, potentials in --potential-unit.\n\n"
                    "Grid types (and the fdf flag that makes SIESTA write them):\n" + type_lines,
        epilog="Example usage:\n"
               "  stb-cube -t RHO                          # the only *.RHO in this folder\n"
               "  stb-cube -l si -t RHO                    # si.RHO -> si_RHO.cube\n"
               "  stb-cube -l run/h2o -t VH --potential-unit Ha -o h2o_vh.cube\n"
               "  stb-cube -l o2 -t RHO --spin diff        # magnetization density up - down\n"
               "  stb-cube -l si -t DRHO --dry-run         # report only, no cube",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-l", "--label",
                        help="SystemLabel of the run; may include a folder (run/si). Default: the only <label>.<type>\n"
                             "in the current folder.")
    parser.add_argument("-t", "--type", required=True, type=str.upper, choices=list(GRID_TYPES),
                        metavar="{" + ",".join(GRID_TYPES) + "}", help="Grid file to convert (see the list above).")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file. Default: <label>_<type>.cube (<label>_<type>_<spin>.cube for a spin\n"
                             "channel other than total), next to the input.")
    parser.add_argument(
        "--spin", choices=SPIN_CHOICES, default="total",
        help="Spin channel of a spin-polarized (collinear) grid: 'total' (default; up + down for a density, their\n"
             "average for a potential), 'up', 'down' or 'diff' (up - down: the magnetization density, or the\n"
             "exchange splitting of a potential). Ignored, with a warning, for an unpolarized grid and for a\n"
             "non-collinear/spin-orbit one. The spin configuration comes from the SIESTA log, else <label>.HSX.")
    parser.add_argument("--potential-unit", choices=list(POTENTIAL_UNITS), default="eV",
                        help="Unit of a potential (VT, VH, VNA) in the cube: eV (default), Ry (SIESTA's own) or Ha\n"
                             "(atomic units). Densities are always e/Bohr^3.")
    parser.add_argument("--dry-run", action="store_true", help="Report only; write no cube.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also write the report to {REPORT_FILE} (next to the input). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-cube {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Grid Converter",
            "SIESTA real-space grids to Gaussian cube",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    gtype = args.type
    spec = GRID_TYPES[gtype]
    warnings = []

    def warn(msg):
        warnings.append(msg)

    # ---- which file -------------------------------------------------------------------------------------------
    label = args.label
    if label is None:
        found = find_labels(gtype)
        if len(found) != 1:
            print(color_text(f"[FAIL] --label not given and {len(found)} *.{gtype} file(s) in this folder"
                             + (f" ({', '.join(f + '.' + gtype for f in found)})" if found else "")
                             + " -- name the run with -l.", 'red'))
            sys.exit(1)
        label = found[0]
    folder = os.path.dirname(label) or "."
    base = os.path.basename(label)
    input_file = f"{label}.{gtype}"
    xv_file = f"{label}.XV"
    default_output = (f"{label}_{gtype}.cube" if args.spin == "total" else f"{label}_{gtype}_{args.spin}.cube")
    output_file = args.output or default_output

    report_path = os.path.join(folder, REPORT_FILE) if args.save_report else None
    buf = io.StringIO() if report_path else None
    sec = Sections(buf)

    def out(text=""):
        print_dual(nz(text), buf)

    def fail(message):
        out(color_text(f"[FAIL] {message}", 'red'))
        if buf is not None:
            with open(report_path, "w") as fh:
                fh.write(buf.getvalue())
        sys.exit(1)

    print("\n" + color_text("GRID CONVERTER:", 'bold'))
    print("-" * 60)
    out(color_text("===== STB-CUBE REPORT =====", 'magenta'))

    # [0] ------------------------------------------------------------------------------------------------------
    sec("RUN METADATA")
    out(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out(f"stb / sisl        : {VERSION} / {sisl.__version__}")
    out(f"Input             : {input_file}" + ("" if args.label else "  (the only *.%s here)" % gtype))
    out(f"Grid type         : {gtype} -- {spec['what']}")
    out(f"Spin channel      : {args.spin} (requested)")
    if spec["kind"] == "potential":
        out(f"Potential unit    : {args.potential_unit} in the cube")
    out(f"Output            : {'none (--dry-run)' if args.dry_run else output_file}")

    # [1] ------------------------------------------------------------------------------------------------------
    sec("INPUT GRID")
    if not os.path.exists(input_file):
        here = sorted(glob.glob(os.path.join(folder, "*.*")))
        grids = sorted({os.path.splitext(f)[1][1:] for f in here if os.path.splitext(f)[1][1:] in GRID_TYPES})
        fail(f"Input file '{input_file}' not found."
             + (f" Grid files in {folder}: {', '.join(grids)}." if grids else "")
             + f" SIESTA writes it only with {spec['flag']}.")
    try:
        sile = sisl.get_sile(input_file)
        nspin, mesh = sile.read_grid_size()
        lattice = sile.read_lattice()
    except Exception as e:
        fail(f"Could not read grid: {e}")
    cell = np.array(lattice.cell, dtype=float)
    shape = tuple(int(m) for m in mesh)
    npts = int(np.prod(shape))
    st = os.stat(input_file)
    out(f"File              : {input_file} ({human_size(st.st_size)}, written "
        f"{datetime.fromtimestamp(st.st_mtime):%Y-%m-%d %H:%M})")
    out(f"Quantity          : {spec['what']}")
    out(f"Written by        : {spec['flag']} in the .fdf")
    native = "e/Bohr^3" if spec["kind"] == "density" else "Ry"
    out(f"Stored as         : {native} on a Bohr cell, single precision, {nspin} component(s)")
    out(f"Mesh              : {shape[0]} x {shape[1]} x {shape[2]} = {npts} points")
    spacing = plane_spacing(cell, shape)
    out(f"Plane spacing     : {spacing[0]:.4f} / {spacing[1]:.4f} / {spacing[2]:.4f} Ang along a / b / c")
    lengths, angles, volume = cell_parameters(cell)
    dv = volume / npts
    out(f"Voxel volume      : {dv:.5f} Ang^3 ({dv / BOHR**3:.5f} Bohr^3)")
    ecut, ecut_axes = effective_mesh_cutoff(cell, shape)
    out(f"Effective cutoff  : {ecut:.3f} Ry (min over a/b/c of {ecut_axes[0]:.1f} / {ecut_axes[1]:.1f} / "
        f"{ecut_axes[2]:.1f} Ry)")
    out(f"Cell lengths      : a = {lengths[0]:.4f}, b = {lengths[1]:.4f}, c = {lengths[2]:.4f} Ang")
    out(f"Cell angles       : alpha = {angles[0]:.2f}, beta = {angles[1]:.2f}, gamma = {angles[2]:.2f} deg")
    out(f"Cell volume       : {volume:.4f} Ang^3")
    out("Cell vectors (Ang):")
    for name, row in zip("abc", cell):
        out(f"   {name}  {row[0]:11.5f} {row[1]:11.5f} {row[2]:11.5f}")

    # [2] ------------------------------------------------------------------------------------------------------
    sec("SIESTA RUN")
    log_path = siesta_log.find_out_file(folder, base)
    log = read_log_info(log_path) if log_path else None
    if log is None:
        out("Log               : none found (<label>.out or the only *.out next to the grid) -- the checks against")
        out("                    SIESTA's own numbers (electrons, mesh, dipole, spin moment) are skipped")
    else:
        out(f"Log               : {log_path}" + (f" (SIESTA {log['version']})" if log["version"] else ""))
        out(f"Finished          : {'yes' if log['finished'] else 'NO -- the grid may be from an unfinished run'}")
        if not log["finished"]:
            warn(f"{log_path} did not finish: the grid is from the last step written, which may not be converged.")
        if log["electrons"] is not None:
            out(f"Electrons         : {log['electrons']:.6f} (valence)")
        if log["spin"]:
            out(f"Spin (log)        : {log['spin']}")
        if log["mesh"]:
            same = tuple(log["mesh"]) == shape
            out(f"Mesh (log)        : {log['mesh'][0]} x {log['mesh'][1]} x {log['mesh'][2]} -- "
                + ("same as the grid" if same else "DIFFERENT from the grid"))
            if not same:
                warn(f"the log's last mesh {log['mesh']} differs from the grid's {shape}: the grid and the log may "
                     "come from different runs.")
        if log["cutoff_used"] is not None:
            diff = abs(log["cutoff_used"] - ecut)
            out(f"Mesh cutoff (log) : required {log['cutoff_required']:.3f} Ry, used {log['cutoff_used']:.3f} Ry -- "
                f"from the grid {ecut:.3f} Ry" + (" (agrees)" if diff < 0.01 * ecut else " (differs)"))
            if diff >= 0.01 * ecut:
                warn(f"the effective cutoff of the grid ({ecut:.2f} Ry) differs from the log's last one "
                     f"({log['cutoff_used']:.2f} Ry); a variable-cell run changes it at every step -- or the files "
                     "come from different runs.")

    # [3] ------------------------------------------------------------------------------------------------------
    sec("SPIN")
    spin_kind, spin_source = None, None
    if log is not None and log["spin"]:
        spin_kind, spin_source = spin_kind_from_log(log["spin"]), os.path.basename(log_path)
    if spin_kind is None:
        hsx_kind, hsx_msg = detect_spin_configuration(label)
        if hsx_msg:
            warn(hsx_msg)
        if hsx_kind:
            spin_kind, spin_source = hsx_kind, f"{base}.HSX"
    expected_nspin = {"unpolarized": 1, "polarized": 2, "non-collinear": 4, "spin-orbit": 4}
    if spin_kind:
        out(f"Spin configuration (from {spin_source}): {spin_kind}")
        exp = expected_nspin[spin_kind]
        if exp != nspin:
            warn(f"{spin_source} reports {spin_kind} ({exp} component(s) expected) but '{input_file}' has {nspin} "
                 "-- they may not belong to the same calculation.")
    else:
        out("Spin configuration: not stated (no log line, no .HSX) -- taken from the grid's component count")
        spin_kind = {1: "unpolarized", 2: "polarized"}.get(nspin, "non-collinear")
    out(f"Components        : {nspin}")

    effective_spin = args.spin
    if nspin == 1:
        if args.spin != "total":
            warn(f"'{input_file}' is not spin-polarized (1 component) -- ignoring --spin {args.spin}, writing the "
                 "total.")
            effective_spin = "total"
    elif nspin != 2 or spin_kind in ("non-collinear", "spin-orbit"):
        if args.spin != "total":
            warn(f"--spin {args.spin} assumes a collinear polarized calculation; this grid is {spin_kind} "
                 f"({nspin} components: a charge/mx/my/mz spin texture, not an up/down split) -- writing the total.")
        effective_spin = "total"

    channels = {}
    try:
        if nspin == 2:
            for ch in SPIN_CHOICES:
                channels[ch] = sile.read_grid(index=channel_weights(spec["kind"], nspin, ch)).grid
        else:
            channels["total"] = sile.read_grid(index="total" if spec["kind"] == "density" else
                                               ([1.0] if nspin == 1 else [0.5, 0.5])).grid
    except Exception as e:
        fail(f"Could not read grid: {e}")

    if nspin == 2:
        rows = []
        moment = log["spin_moment"] if log else None
        nel = log["electrons"] if log else None
        for ch, label_ch in (("total", "total (up + down)" if spec["kind"] == "density" else "total (average)"),
                             ("up", "up"), ("down", "down"), ("diff", "up - down")):
            data = channels[ch]
            if spec["kind"] == "density":
                integ = data.sum() * dv
                ref = ""
                if gtype == "RHO" and nel is not None:
                    if ch == "total":
                        ref = f"{nel:.4f}"
                    elif moment is not None:
                        ref = {"up": f"{(nel + moment) / 2:.4f}", "down": f"{(nel - moment) / 2:.4f}",
                               "diff": f"{moment:.4f}"}[ch]
                rows.append(([label_ch, f"{integ:.4f}", ref, "<- written" if ch == effective_spin else ""], None))
            else:
                rows.append(([label_ch, f"{data.mean():.4f}", "<- written" if ch == effective_spin else ""], None))
        if spec["kind"] == "density":
            print_table(["channel", "integral (e)", "SIESTA", ""], rows, buf)
            m = channels["diff"]
            out(f"Absolute moment   : {np.abs(m).sum() * dv:.4f} e (integral of |up - down|; above the net moment "
                "where some region is polarized the other way)")
            if moment is not None:
                diff_int = channels["diff"].sum() * dv
                out(f"Spin moment (log) : {moment:.5f} -- grid {diff_int:.5f}"
                    + (" (agrees)" if abs(diff_int - moment) < 0.01 else " (differs)"))
                if abs(diff_int - moment) >= 0.01:
                    warn(f"the up - down integral ({diff_int:.4f}) differs from SIESTA's spin moment ({moment:.4f}).")
        else:
            print_table(["channel", "mean (eV)", ""], rows, buf)
            out("Note              : a potential's 'total' is the average of up and down, not their sum")
        if effective_spin == args.spin and effective_spin != "total":
            out(f"Converting        : the '{effective_spin}' channel")
    data_ang = channels[effective_spin]

    # [4] ------------------------------------------------------------------------------------------------------
    sec("GEOMETRY")
    geom = None
    zval = None
    if log_path:
        with contextlib.redirect_stdout(io.StringIO()):
            try:
                zval = siesta_log.get_zval_from_output(base, override_path=log_path)
            except Exception:
                zval = None
    if os.path.exists(xv_file):
        try:
            geom = sisl.get_sile(xv_file).read_geometry()
        except Exception as e:
            warn(f"could not read the geometry from {xv_file}: {e} -- the cube has no atoms.")
    else:
        warn(f"no {xv_file}: the cube has no atoms (the grid itself is complete).")
    vac = [False, False, False]
    if geom is None:
        out(f"Geometry          : none ({xv_file} missing or unreadable)")
    else:
        symbols = [a.symbol for a in geom.atoms]
        counts = {}
        for s_ in symbols:
            counts[s_] = counts.get(s_, 0) + 1
        formula = "".join(f"{k}{v if v > 1 else ''}" for k, v in counts.items())
        out(f"Source            : {xv_file} ({len(symbols)} atoms, {formula})")
        cell_diff = float(np.abs(np.asarray(geom.cell) - cell).max())
        if cell_diff < 1e-3:
            out(f"Cell vs grid      : identical (largest difference {cell_diff:.1e} Ang)")
        else:
            out(f"Cell vs grid      : DIFFERENT (largest difference {cell_diff:.4f} Ang)")
            warn(f"lattice mismatch between the {gtype} grid and '{xv_file}' (largest difference {cell_diff:.4f} "
                 "Ang) -- they may not belong to the same calculation; atom positions in the cube may be wrong.")
        ghosts = [i for i, a in enumerate(geom.atoms) if a.Z <= 0]
        frac = geom.xyz @ np.linalg.inv(cell)
        vac = detect_vacuum_axes(frac, cell, VACUUM_GAP)
        vac_names = [n for n, v in zip("abc", vac) if v]
        out(f"Dimensionality    : {dimensionality_label(vac)}"
            + (f" -- vacuum along {', '.join(vac_names)} (gap >= {VACUUM_GAP:g} Ang)" if vac_names else ""))
        rows = []
        for i, (sym, a, pos, fr) in enumerate(zip(symbols, geom.atoms, geom.xyz, frac)):
            if i >= MAX_ATOM_ROWS:
                break
            zv = f"{zval[sym]:g}" if zval and sym in zval else "-"
            rows.append(([str(i + 1), sym, str(a.Z), zv, f"{pos[0]:.4f}", f"{pos[1]:.4f}", f"{pos[2]:.4f}",
                          f"{fr[0]:.4f} {fr[1]:.4f} {fr[2]:.4f}"], None))
        print_table(["#", "el", "Z", "Z_val", "x (Ang)", "y (Ang)", "z (Ang)", "fractional"], rows, buf)
        if len(symbols) > MAX_ATOM_ROWS:
            out(f"   ... and {len(symbols) - MAX_ATOM_ROWS} more atoms")
        if ghosts:
            out(f"Ghost atoms       : {len(ghosts)} (Z <= 0) -- left out of the cube, which has no ghost species")

    # [5] ------------------------------------------------------------------------------------------------------
    sec("VALUES")
    if spec["kind"] == "density":
        to_cube = 1.0 / sile.grid_unit                  # e/Ang^3 -> e/Bohr^3
        cube_unit = "e/Bohr^3"
    else:
        to_cube = 1.0 / POTENTIAL_UNITS[args.potential_unit]
        cube_unit = args.potential_unit
    vals = data_ang * to_cube
    vmin, vmax = float(vals.min()), float(vals.max())
    imin = np.unravel_index(int(np.argmin(vals)), vals.shape)
    imax = np.unravel_index(int(np.argmax(vals)), vals.shape)

    def where(idx):
        point = (np.asarray(idx) / np.asarray(shape)) @ cell
        i, d = nearest_atom(point, geom, cell)
        loc = f"({point[0]:.3f}, {point[1]:.3f}, {point[2]:.3f}) Ang"
        return loc + (f", {d:.3f} Ang from {geom.atoms[i].symbol}{i + 1}" if i is not None else "")

    rows = [(["minimum", fmt_num(vmin, 5), where(imin)], None),
            (["maximum", fmt_num(vmax, 5), where(imax)], None),
            (["mean", fmt_num(float(vals.mean()), 5), "over the cell"], None),
            (["std. deviation", fmt_num(float(vals.std()), 5), ""], None)]
    if spec["kind"] == "potential":
        rows.append((["max - min", fmt_num(vmax - vmin, 5), ""], None))
    print_table(["", f"value ({cube_unit})", "where"], rows, buf)
    if spec["kind"] == "density":
        out(f"In e/Ang^3        : min {fmt_num(float(data_ang.min()), 5)}, max {fmt_num(float(data_ang.max()), 5)}, "
            f"mean {fmt_num(float(data_ang.mean()), 5)}")

    integral = float(data_ang.sum() * dv)
    if spec["kind"] == "density":
        pos_int = float(data_ang[data_ang > 0].sum() * dv)
        neg_int = float(data_ang[data_ang < 0].sum() * dv)
        neg_frac = float((data_ang < 0).mean() * 100)
        out(f"Integral          : {integral:.6f} e  (positive part {pos_int:.4f}, negative part {neg_int:.4f})")
        expected, why = None, None
        nel = log["electrons"] if log else None
        z_list = [a.Z for a in geom.atoms] if geom is not None else None
        zval_sum = (sum(zval.get(a.symbol, np.nan) for a in geom.atoms) if (geom is not None and zval) else None)
        if effective_spin == "total":
            if spec["expect"] == "electrons" and nel is not None:
                expected, why = nel, "valence electrons (log)"
            elif spec["expect"] == "free-atoms" and zval_sum is not None:
                expected, why = zval_sum, "sum of Z_val (free neutral atoms)"
            elif spec["expect"] == "zero":
                expected, why = 0.0, "zero by construction"
            elif spec["expect"] == "all-electron" and nel is not None and zval_sum is not None and z_list:
                expected = nel + sum(z_list) - zval_sum
                why = "valence electrons + core electrons (sum of Z - Z_val)"
            elif spec["expect"] == "neutral" and nel is not None and zval_sum is not None:
                expected, why = nel - zval_sum, "minus the net charge (0 for a neutral cell)"
            elif spec["expect"] == "ions" and zval_sum is not None:
                expected, why = -zval_sum, "minus the sum of Z_val"
        if expected is not None:
            dev = integral - expected
            # the ionic charges of TOCH/IOCH are Gaussians sampled on the mesh: ~1e-2 e is the mesh's own error
            tol = 0.05 if gtype in ("TOCH", "IOCH") else max(1e-3, 1e-3 * abs(expected))
            ok = abs(dev) < tol
            out(f"Expected          : {expected:.6f} e ({why}) -- difference {dev:+.2e} e"
                + (" (agrees)" if ok else ""))
            if gtype in ("TOCH", "IOCH"):
                out("                    the ionic charges are Gaussians sampled on the mesh: up to ~1e-2 e is the "
                    "mesh's own error")
            if not ok:
                warn(f"the {gtype} integral ({integral:.4f} e) differs from the expected {expected:.4f} e ({why}) "
                     f"by more than {tol:g} e -- a different run, an unconverged SCF, or a too coarse mesh.")
        elif spec["expect"] and effective_spin == "total":
            out("Expected          : unknown (needs the SIESTA log" + (" and Z_val" if spec["expect"] != "zero" else "")
                + ")")
        if gtype == "DRHO" and effective_spin == "total":
            out(f"Charge moved      : {pos_int:.4f} e gathered where Delta rho > 0, the same amount taken from "
                "where it is < 0")
        if not spec["signed"] and effective_spin in ("total", "up", "down"):
            if neg_frac > 0:
                out(f"Negative values   : {neg_frac:.2f} % of the points, {neg_int:.2e} e in all")
                if neg_int < -1e-3:
                    warn(f"{gtype} has negative values ({neg_int:.3e} e in {neg_frac:.1f} % of the points): a density "
                         "cannot be negative -- a basis/mesh artefact; isosurface programs ignore it.")
            else:
                out("Negative values   : none")
        levels = iso_levels(vals, dv / BOHR**3, ISO_FRACTIONS, spec["signed"] or effective_spin == "diff")
        if levels:
            signed = spec["signed"] or effective_spin == "diff"
            txt = ", ".join(f"{'+-' if signed else ''}{fmt_num(v, 3)} ({int(f * 100)} %)" for f, v in levels)
            out(f"Isosurface levels : {txt} {cube_unit}")
            out("                    each level encloses that fraction of the "
                + ("integral of |values|" if signed else "integral") + " -- a starting point in VESTA/VMD")
    else:
        pct = np.percentile(vals, [1, 50, 99])
        out(f"Percentiles       : 1 % {fmt_num(pct[0], 5)}, median {fmt_num(pct[1], 5)}, 99 % {fmt_num(pct[2], 5)} "
            f"{cube_unit}")
        out("Reference         : SIESTA's potentials are defined up to a constant; a periodic cell's average is "
            "not a vacuum level")

    # vacuum planes: what the grid looks like far from the atoms
    if geom is not None and any(vac):
        frac = geom.xyz @ np.linalg.inv(cell)
        rows, not_flat = [], []
        for ax in range(3):
            if not vac[ax]:
                continue
            f_ax = np.sort(np.mod(frac[:, ax], 1.0))
            gaps = np.diff(np.concatenate([f_ax, [f_ax[0] + 1.0]]))
            k = int(np.argmax(gaps))
            centre = (f_ax[k] + gaps[k] / 2.0) % 1.0
            plane = int(round(centre * shape[ax])) % shape[ax]
            other = tuple(i for i in range(3) if i != ax)
            profile = vals.mean(axis=other)
            gap_ang = gaps[k] * lengths[ax]
            half = max(1, int(round(min(1.0, gap_ang / 4) / spacing[ax])))
            window = [profile[(plane + j) % shape[ax]] for j in range(-half, half + 1)]
            slab = np.take(vals, plane, axis=ax)
            centre = 0.0 if centre > 0.9995 else centre
            variation = float(np.max(window) - np.min(window))
            if spec["kind"] == "potential" and variation > 0.01 * (POTENTIAL_UNITS["eV"] / POTENTIAL_UNITS[
                    args.potential_unit]):
                not_flat.append(("abc"[ax], variation))
            rows.append(([f"{'abc'[ax]}", f"{centre:.3f}", f"{gap_ang:.2f}", fmt_num(float(profile[plane]), 5),
                          fmt_num(variation, 3),
                          fmt_num(float(np.abs(slab).max()), 3)], None))
        out("Vacuum planes (the plane farthest from the atoms along each vacuum axis):")
        print_table(["axis", "frac", "gap (Ang)", f"plane mean ({cube_unit})", "variation +-1 Ang",
                     "max |value| in plane"], rows, buf)
        if spec["kind"] == "potential":
            out("                    the plane mean is the vacuum level only where the potential is flat "
                "(compare stb-workfunction)")
            for name, var in not_flat:
                out(f"Not flat along {name}  : the plane mean changes by {fmt_num(var, 3)} {cube_unit} within +-1 Ang -- "
                    "a dipole across a periodic box (a slab needs Slab.DipoleCorrection); no single vacuum level")
        else:
            worst = max(float(r[0][5]) for r in rows)
            if gtype in ("RHO", "BADER", "RHOXC", "RHOINIT") and worst > 1e-4:
                warn(f"the density in the vacuum plane reaches {worst:.1e} {cube_unit}: the box may be too small "
                     "for the orbitals' tails (or atoms are close to the boundary).")

    # dipole of a molecule/slab, along the non-periodic axes only
    if (geom is not None and any(vac) and spec["kind"] == "density" and effective_spin == "total"
            and gtype in ("RHO", "DRHO", "TOCH")):
        inv = np.linalg.inv(cell)
        centre_frac = (geom.xyz @ inv).mean(axis=0)
        idx = np.indices(shape).reshape(3, -1).T / np.asarray(shape)
        dfrac = idx - centre_frac
        dfrac -= np.round(dfrac)
        r = dfrac @ cell
        p = -(data_ang.reshape(-1)[:, None] * r).sum(axis=0) * dv       # electrons are negative
        point_ok = True
        if gtype == "RHO":
            if zval:
                dr_at = geom.xyz @ inv - centre_frac
                dr_at -= np.round(dr_at)
                charges = np.array([zval.get(a.symbol, np.nan) for a in geom.atoms])
                if np.isnan(charges).any():
                    point_ok = False
                else:
                    p = p + (charges[:, None] * (dr_at @ cell)).sum(axis=0)
            else:
                point_ok = False
        if point_ok:
            p_debye = p * DEBYE_PER_E_ANG
            rows = []
            for ax in range(3):
                if not vac[ax]:
                    continue
                normal = np.cross(cell[(ax + 1) % 3], cell[(ax + 2) % 3])
                normal /= np.linalg.norm(normal)
                if np.dot(normal, cell[ax]) < 0:
                    normal = -normal
                ours = float(np.dot(p_debye, normal))
                theirs = (float(np.dot(log["dipole_debye"], normal))
                          if (log and log["dipole_debye"] is not None) else None)
                cart = [c for c, e in zip("xyz", np.eye(3)) if abs(np.dot(normal, e)) > 0.9999]
                direction = cart[0] if cart else f"({normal[0]:.3f}, {normal[1]:.3f}, {normal[2]:.3f})"
                rows.append(([f"{'abc'[ax]}: along {direction}", nz(f"{ours:.5f}"),
                              nz(f"{theirs:.5f}") if theirs is not None else "-",
                              (f"{ours - theirs:+.1e}" if theirs is not None else "")], None))
                if theirs is not None and abs(ours - theirs) > max(0.01, 0.01 * abs(theirs)):
                    warn(f"the dipole from the grid along {'abc'[ax]} ({ours:.4f} D) differs from SIESTA's "
                         f"({theirs:.4f} D) -- density leaking across the box boundary, or a different run.")
            how = {"RHO": "density + valence point charges Z_val at the atoms",
                   "DRHO": "deformation density alone (free atoms have no dipole)",
                   "TOCH": "total charge (ions included)"}[gtype]
            out(f"Dipole            : from the grid -- {how}; positions about the atoms' centre, minimum image;")
            out("                    one component per vacuum axis, along the normal to the periodic plane")
            print_table(["component", "grid (Debye)", "SIESTA (Debye)", "difference"], rows, buf)
            if len(rows) == 3:
                out(f"|p|               : {np.linalg.norm(p_debye):.5f} D ({np.linalg.norm(p):.5f} e Ang)")
        else:
            out("Dipole            : skipped (no Z_val per species: needs the SIESTA log)")

    # [6] ------------------------------------------------------------------------------------------------------
    sec("OUTPUT")
    if args.dry_run:
        out("Dry run           : no cube written")
    else:
        numbers, positions = [], np.zeros((0, 3))
        if geom is not None:
            keep = [i for i, a in enumerate(geom.atoms) if a.Z > 0]
            numbers = [geom.atoms[i].Z for i in keep]
            positions = geom.xyz[keep]
        comment1 = (f"STB-CUBE {VERSION}: {base}.{gtype}" + (f" spin {effective_spin}" if nspin > 1 else "")
                    + f" -- {spec['what']} -- lengths in Bohr, values in {cube_unit}")
        comment2 = "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"      # Gaussian's wording; ASE parses it
        try:
            write_cube(output_file, vals, cell, numbers, positions, comment1, comment2)
        except Exception as e:
            fail(f"Failed to write output file: {e}")
        size = os.path.getsize(output_file)
        out(f"File              : {output_file} ({human_size(size)})")
        out(f"Units             : lengths Bohr, values {cube_unit}"
            + (" -- the integral over the cube's own voxels gives electrons" if spec["kind"] == "density" else ""))
        out(f"Header            : 2 comment lines, origin (0, 0, 0), {shape[0]} x {shape[1]} x {shape[2]} voxels, "
            f"{len(numbers)} atoms")
        out(f"Values            : {npts} numbers, 6 significant digits ({CUBE_FMT}), x slowest, z fastest")
        try:
            back, atoms_back = read_back_cube(output_file)
            problems = []
            if back.shape != vals.shape:
                problems.append(f"shape {back.shape} instead of {vals.shape}")
            scale = max(abs(vmin), abs(vmax), 1e-30)
            rel = float(np.abs(back - vals).max() / scale) if back.shape == vals.shape else np.inf
            if rel > 1e-5:
                problems.append(f"values differ by {rel:.1e} of the largest")
            cell_err = float(np.abs(np.asarray(atoms_back.cell) - cell).max())
            if cell_err > 1e-4:
                problems.append(f"cell differs by {cell_err:.1e} Ang")
            if len(atoms_back) != len(numbers):
                problems.append(f"{len(atoms_back)} atoms instead of {len(numbers)}")
            elif len(numbers):
                pos_err = float(np.abs(atoms_back.positions - positions).max())
                if pos_err > 1e-4:
                    problems.append(f"atoms moved by {pos_err:.1e} Ang")
                if list(atoms_back.numbers) != list(numbers):
                    problems.append("atomic numbers differ")
            out("Read-back check   : read again with ASE's cube reader -- "
                + ("ok" if not problems else "FAILED: " + "; ".join(problems)))
            out(f"                    values agree to {rel:.1e} of the largest; cell, {len(numbers)} atoms and "
                "atomic numbers as written")
            if spec["kind"] == "density":
                out(f"Integral (cube)   : {back.sum() * dv / BOHR**3:.6f} e  (values x voxel volume in Bohr^3)")
            if problems:
                warn("the cube read back does not match what was written: " + "; ".join(problems))
        except Exception as e:
            warn(f"could not read the cube back to check it: {e}")

    # [7] ------------------------------------------------------------------------------------------------------
    sec("SUMMARY & FILES")
    out(f"Converted         : {input_file} [{gtype}"
        + (f", {effective_spin}" if nspin > 1 else "") + f"] -> "
        + ("nothing (--dry-run)" if args.dry_run else output_file))
    if spec["kind"] == "density":
        out(f"Integral          : {integral:.4f} e")
    else:
        out(f"Range             : {fmt_num(vmin, 5)} to {fmt_num(vmax, 5)} {cube_unit}")
    if not args.dry_run:
        out(f"Open with         : VESTA {output_file}  |  vmd {output_file}  |  xcrysden --cube {output_file}")
    if report_path:
        out(f"Report            : {report_path}")

    sec("WARNINGS")
    if warnings:
        for i, w in enumerate(warnings, 1):
            out(color_text(f"  {i}. {w}", 'yellow'))
    else:
        out("  none")

    if buf is not None:
        with open(report_path, "w") as fh:
            fh.write(buf.getvalue())

    print("\n" + "-" * 60)
    print(color_text("Conversion complete. It's hip to be square (cube).\n", 'bold'))


if __name__ == "__main__":
    main()
