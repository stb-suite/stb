#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 2 of the Charge Density Difference workflow: Delta rho = rho_combined - sum_i rho_fragment_i.

Besides the map itself (2D slice / planar-averaged profile / 3D point cloud / Gaussian cube), the report gives
what is needed to trust and to interpret it:

- the consistency of the folders (same grid, same cell, same mesh cutoff and level of theory, every SCF
  converged, the electron count of each .RHO against SIESTA's own, charge conservation of Delta rho);
- the counterpoise-corrected interaction energy E(combined) - sum_i E(fragment_i with the others as ghosts),
  from the same runs;
- statistics of Delta rho (charge moved, extremes with their positions and nearest atoms, rms);
- the redistribution of charge between fragments and atoms, by a nearest-atom partition of the grid (a
  geometric partition, NOT a Bader analysis: it depends on the atoms' size only through the distance);
- a 4-panel dashboard (--png / --show), gnuplot files (--save-gnuplot) and the numbers as CSV/JSON (--save-data).

The cube file is written in e/bohr^3, the unit every cube reader assumes (sisl writes the grid unit-less: the
e/Ang^3 grid is converted here).
"""

from stb import __version__ as VERSION

import os
import re
import sys
import glob
import json
import argparse
from datetime import datetime

import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from stb.core import citations, structure_io, siesta_log, fdf_directives
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.rho_io import read_total, try_read_net_spin, find_one_rho, compute_colorbar_range
from stb.core.grid_export import (write_data_file, write_profile_data_file, write_gnuplot_script,
                                  check_planar_orthogonality, NON_ORTHOGONAL_WARN_DEG)

REPORT_FILE = "chargediff_analysis_report.txt"
BOHR_TO_ANG = 0.529177210903
BIB_FILE = "references.bib"
CHARGEDIFF_MANIFEST_FILE = "chargediff_manifest.json"
EV_TO_KCAL = 23.0605
EV_TO_KJ = 96.4853
ELECTRON_TOL = 1e-3          # e: conservation of the electron count
SETTINGS_KEYS = ["xcfunctional", "xcauthors", "paobasistype", "paobasissize", "paoenergyshift", "meshcutoff",
                 "spin", "electronictemperature", "dmtolerance", "dftd3"]
SETTINGS_LABELS = {"xcfunctional": "XC functional", "xcauthors": "XC authors", "paobasistype": "Basis type",
                   "paobasissize": "Basis size", "paoenergyshift": "PAO energy shift", "meshcutoff": "MeshCutoff",
                   "spin": "Spin", "electronictemperature": "Electronic temperature",
                   "dmtolerance": "DM tolerance", "dftd3": "DFT-D3", "kgrid": "SCF k-grid"}


def find_geometry_file(folder):
    """Returns a (path, sisl Geometry) for --cube's geometry, read from
    whichever finished '*.XV' is in `folder` (preferred, the actual
    converged geometry) or its 'structure.fdf' (written by
    stb-chargediffPrep) otherwise. Returns (None, None) if neither exists.
    """
    xv_matches = sorted(glob.glob(os.path.join(folder, "*.XV")))
    if xv_matches:
        path = xv_matches[0]
    else:
        path = os.path.join(folder, "structure.fdf")
        if not os.path.isfile(path):
            return None, None
    return path, sisl.get_sile(path).read_geometry()


# ---- reading the runs -------------------------------------------------------------------

def folder_settings(folder):
    """Settings that shape the energy surface as SIESTA reads them in `folder` (config_extra.fdf first),
    plus the k-grid and the set of pseudopotential hashes (labels differ between folders by design)."""
    import hashlib
    texts = []
    for name in ("config_extra.fdf", "calc.fdf"):
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            with open(path, errors="replace") as fh:
                texts.append(fh.read())
    d = fdf_directives.parse_directives(*texts)
    out = {k: " ".join(d[k].split()) for k in SETTINGS_KEYS if k in d}
    kg = fdf_directives.parse_kgrid(*texts)
    out["kgrid"] = "x".join(str(v) for v in kg) if kg else "Gamma"
    hashes = sorted({hashlib.md5(open(p, "rb").read()).hexdigest()[:8]
                     for p in glob.glob(os.path.join(folder, "*.psml")) + glob.glob(os.path.join(folder, "*.psf"))})
    return out, hashes


def read_run(folder):
    """Energy, SCF convergence, electron count, mesh and SIESTA version of the run in `folder` (None when absent)."""
    out_path = siesta_log.find_out_file(folder, None)
    info = {"out": out_path, "energy": None, "scf": None, "n_electrons": None, "mesh_used": None,
            "mesh": None, "version": None}
    if not out_path:
        return info
    info["energy"] = siesta_log.get_free_energy(out_path)
    info["scf"] = siesta_log.get_scf_convergence(out_path)[0]
    info["version"] = siesta_log.get_siesta_version(out_path)
    with open(out_path, errors="replace") as fh:
        for line in fh:
            m = re.match(r"\s*Total number of electrons:\s*([0-9.]+)", line)
            if m:
                info["n_electrons"] = float(m.group(1))
            m = re.match(r"\s*InitMesh: Mesh cutoff \(required, used\) =\s*([0-9.]+)\s+([0-9.]+)", line)
            if m:
                info["mesh"], info["mesh_used"] = float(m.group(1)), float(m.group(2))
    return info


def read_atoms(folders, fragment_folders):
    """[(label, element, fragment index or -1, Cartesian position in Ang)] of the combined system, in file order;
    the fragment of an atom is the folder where its label is not a ghost."""
    path = os.path.join(folders["combined"], "structure.fdf")
    st = structure_io.read_fdf(path)
    lat = np.asarray(st.lattice)
    pos = np.array([p for _, p in st.atoms], float)
    if st.coord_format != "cartesian":
        pos = pos @ lat
    frag_of = [-1] * len(st.atoms)
    own_label = [None] * len(st.atoms)
    inv = np.linalg.inv(lat)
    for k, name in enumerate(fragment_folders):
        fst = structure_io.read_fdf(os.path.join(folders[name], "structure.fdf"))
        fpos = np.array([p for _, p in fst.atoms], float)
        if fst.coord_format != "cartesian":
            fpos = fpos @ np.asarray(fst.lattice)
        # the folders list their atoms grouped by species label, so the order differs from the combined one:
        # match atoms by position (minimum image)
        for i in range(len(st.atoms)):
            d = (fpos - pos[i]) @ inv
            d -= np.round(d)
            j = int(np.argmin(np.linalg.norm(d @ lat, axis=1)))
            if not fst.atoms[j][0].endswith("_ghost"):
                frag_of[i] = k
                own_label[i] = fst.atoms[j][0]
    atoms = []
    for i, (label, _p) in enumerate(st.atoms):
        atoms.append((own_label[i] or label, structure_io.real_element(label, st.species_meta), frag_of[i], pos[i]))
    return atoms, lat


# ---- numerics ---------------------------------------------------------------------------

def nearest_atom_partition(diff, lattice, atom_pos):
    """Net quantity (e) of the grid assigned to each atom by nearest-atom partition (minimum image on the
    fractional coordinates, exact for orthogonal-ish cells). Returns None when the grid is too large."""
    nx, ny, nz = diff.shape
    n_atoms = len(atom_pos)
    if diff.size * n_atoms > 4e9:
        return None
    inv = np.linalg.inv(lattice)
    atom_frac = atom_pos @ inv
    dv = abs(np.linalg.det(lattice)) / diff.size
    totals = np.zeros(n_atoms)
    j, k = np.meshgrid(np.arange(ny) / ny, np.arange(nz) / nz, indexing="ij")
    chunk = max(1, int(2e6 // (ny * nz * max(n_atoms, 1)) or 1))
    for i0 in range(0, nx, chunk):
        i1 = min(nx, i0 + chunk)
        ii = (np.arange(i0, i1) / nx)[:, None, None]
        frac = np.stack(np.broadcast_arrays(ii, j[None], k[None]), axis=-1)       # (c, ny, nz, 3)
        best = np.full(frac.shape[:3], np.inf)
        owner = np.zeros(frac.shape[:3], dtype=np.int32)
        for a in range(n_atoms):
            d = frac - atom_frac[a]
            d -= np.round(d)
            dist2 = ((d @ lattice) ** 2).sum(axis=-1)
            better = dist2 < best
            best = np.where(better, dist2, best)
            owner = np.where(better, a, owner)
        totals += np.bincount(owner.ravel(), weights=diff[i0:i1].ravel(), minlength=n_atoms)
    return totals * dv


def nearest_atom(point_cart, lattice, atom_pos, atoms):
    inv = np.linalg.inv(lattice)
    d = (atom_pos @ inv) - (point_cart @ inv)
    d -= np.round(d)
    dist = np.linalg.norm(d @ lattice, axis=1)
    k = int(np.argmin(dist))
    return k, float(dist[k])


def voxel_to_cart(idx, shape, lattice, origin):
    return (np.array(idx, float) / np.array(shape)) @ lattice + np.asarray(origin)


def plane_basis(lattice, axis):
    """(e1, e2, normal): orthonormal in-plane basis of the plane spanned by the two lattice vectors other than `axis`."""
    ia, ib = [a for a in (0, 1, 2) if a != axis]
    a, b = lattice[ia], lattice[ib]
    e1 = a / np.linalg.norm(a)
    b_perp = b - np.dot(b, e1) * e1
    e2 = b_perp / np.linalg.norm(b_perp)
    n = np.cross(e1, e2)
    return e1, e2, n / np.linalg.norm(n)


# ---- dashboard --------------------------------------------------------------------------

def make_dashboard(ctx, png_path, show, fallback_png=None):
    """4-panel matplotlib figure: slice map with atoms, planar profile with the cumulative charge, charge change per
    atom (coloured by fragment) and a text summary. Saved to `png_path` and/or shown. Returns (produced, shown)."""
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(color_text("[WARNING] matplotlib is not installed -- skipping the dashboard.", 'yellow'))
        return False, False
    if show and (matplotlib.get_backend().lower() == "agg"
                 or (sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")))):
        print(color_text("[WARNING] --show needs a display (local desktop or `ssh -X`); saving the figure instead.", 'yellow'))
        show = False
        png_path = png_path or fallback_png
        plt.switch_backend("Agg")

    diff, lat, origin, axis = ctx["diff"], ctx["lattice"], ctx["origin"], ctx["axis"]
    n = diff.shape
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))

    # (0,0) slice
    a = ax[0, 0]
    idx = ctx["slice_idx"]
    sl = np.take(diff, idx, axis=axis)
    ia, ib = [t for t in (0, 1, 2) if t != axis]
    e1, e2, nrm = plane_basis(lat, axis)
    ua = (np.arange(n[ia]) / n[ia])[:, None] * np.ones((1, n[ib]))
    ub = np.ones((n[ia], 1)) * (np.arange(n[ib]) / n[ib])[None, :]
    P = ua[..., None] * lat[ia] + ub[..., None] * lat[ib]
    X, Y = P @ e1, P @ e2
    # a robust range (99.5th percentile of |Delta rho| in the cut): the largest values sit on top of the nuclei
    # and would wash out the bonding features
    nz = np.abs(sl[np.abs(sl) > 0])
    vmax = float(np.percentile(nz, 99.5)) if nz.size else 1.0
    if ctx["cb"] and ctx["cb"][1] < vmax:
        vmax = ctx["cb"][1]
    pcm = a.pcolormesh(X, Y, sl, cmap="RdBu", vmin=-vmax, vmax=vmax, shading="auto")
    if ctx["contour"]:
        a.contour(X, Y, sl, levels=np.linspace(-vmax, vmax, 13), colors="k", linewidths=0.4)
    pos_axis = ctx["slice_pos"]
    for (lab, el, frag, p) in ctx["atoms"]:
        d_n = np.dot(p - np.asarray(origin), nrm) - np.dot(np.asarray(origin) * 0 + (idx / n[axis]) * lat[axis], nrm)
        L = abs(np.dot(lat[axis], nrm))
        d_n -= L * np.round(d_n / L) if L > 0 else 0
        if abs(d_n) < 1.0:
            a.plot(np.dot(p, e1) % np.linalg.norm(lat[ia]) if False else np.dot(p, e1), np.dot(p, e2), "o",
                   ms=14 if el != "H" else 8, mfc=palette[frag % len(palette)] if frag >= 0 else "grey", mec="k", zorder=5)
    if ctx["atoms"]:
        pts = np.array([[np.dot(at[3], e1), np.dot(at[3], e2)] for at in ctx["atoms"]])
        lo, hi = pts.min(axis=0) - 3.0, pts.max(axis=0) + 3.0
        if hi[0] - lo[0] < X.max() - X.min() and hi[1] - lo[1] < Y.max() - Y.min():   # zoom on the atoms (isolated systems)
            a.set_xlim(lo[0], hi[0])
            a.set_ylim(lo[1], hi[1])
    a.set_aspect("equal")
    a.set_xlabel("in-plane coordinate 1 (Å)")
    a.set_ylabel("in-plane coordinate 2 (Å)")
    a.set_title(f"{ctx['quantity']}: plane normal to {'xyz'[axis]} at {pos_axis:.2f} Å (atoms within 1 Å)", fontsize=10)
    fig.colorbar(pcm, ax=a, label="e/Å³ (blue = gain, red = loss)")

    # (0,1) profile
    a = ax[0, 1]
    pax = ctx["prof_axis"]
    prof = diff.mean(axis=tuple(t for t in (0, 1, 2) if t != pax))
    length = float(np.linalg.norm(lat[pax]))
    x = (np.arange(n[pax]) / n[pax]) * length
    area = abs(np.linalg.det(lat)) / length
    planar = prof * area                                   # e/Å: integral of Delta rho over the plane
    cum = np.cumsum(planar) * (length / n[pax])
    a.plot(x, planar, color="#0072B2", lw=1.8, label="planar integral (e/Å)")
    a.plot(x, cum, color="#D55E00", lw=1.8, label="charge below x (e)")
    a.axhline(0, color="k", lw=0.7)
    for (lab, el, frag, p) in ctx["atoms"]:
        a.axvline(float(np.dot(p - np.asarray(origin), lat[pax] / length)) % length,
                  color=palette[frag % len(palette)] if frag >= 0 else "grey", lw=0.6, ls="--" if el != "H" else ":")
    a.set_xlabel(f"position along {'xyz'[pax]} (Å)")
    a.set_title("Planar integral and cumulative charge (dashed: heavy atoms)", fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (1,0) per-atom bars
    a = ax[1, 0]
    if ctx["per_atom"] is not None:
        names = [f"{el}{i + 1}" for i, (lab, el, fr, p) in enumerate(ctx["atoms"])]
        cols = [palette[fr % len(palette)] if fr >= 0 else "grey" for (_l, _e, fr, _p) in ctx["atoms"]]
        a.bar(range(len(names)), ctx["per_atom"], color=cols, edgecolor="k")
        a.set_xticks(range(len(names)))
        a.set_xticklabels(names, rotation=45, fontsize=8)
        a.axhline(0, color="k", lw=0.7)
        handles = [plt.Line2D([0], [0], marker="s", ls="", color=palette[k % len(palette)], label=f"{nm}: {q:+.4f} e")
                   for k, (nm, q) in enumerate(zip(ctx["fragment_names"], ctx["per_fragment"]))]
        a.legend(handles=handles, fontsize=8)
        a.set_ylabel("net Δ charge of the atom's region (e)")
        a.set_title("Charge change by atom (nearest-atom partition)")
        a.grid(alpha=0.3, axis="y")
    else:
        a.axis("off")
        a.text(0.5, 0.5, "per-atom partition skipped (grid too large)", ha="center")

    # (1,1) text
    a = ax[1, 1]
    a.axis("off")
    lines = [f"{ctx['title']}", "", f"{ctx['quantity']}:"]
    lines.append(f"  integral          : {ctx['total']:+.5f} e")
    lines.append(f"  charge moved      : {ctx['moved']:.4f} e  (= half of ∫|Δρ|)")
    lines.append(f"  max / min         : {ctx['max'][0]:+.4f} / {ctx['min'][0]:+.4f} e/Å³")
    lines.append(f"  rms               : {ctx['rms']:.5f} e/Å³")
    if ctx["e_int"] is not None:
        lines += ["", "Counterpoise interaction energy:",
                  f"  {ctx['e_int']:+.4f} eV = {ctx['e_int'] * EV_TO_KCAL:+.2f} kcal/mol"]
    if ctx["per_fragment"] is not None:
        lines += ["", "Net charge change by fragment:"]
        lines += [f"  {nm:<10s}: {q:+.4f} e" for nm, q in zip(ctx["fragment_names"], ctx["per_fragment"])]
    lines += ["", f"checks: {ctx['n_ok']} ok / {ctx['n_warn']} warn"]
    a.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=10)

    fig.suptitle(f"Charge density difference — {ctx['title']}", fontsize=14)
    fig.tight_layout()
    if png_path:
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return True, show


# ---- main -------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of the Charge Density Difference workflow: given just "
        "--dir (the same folder stb-chargediffPrep wrote into), auto-discovers the combined system "
        "and every fragment from chargediff_manifest.json, reads each one's own .RHO (once SIESTA "
        "has been run in every folder), and computes Delta rho = rho_combined - "
        "sum(rho_fragment_i) -- the standard way to visualize charge accumulation/depletion "
        "regions.", 'bold')}
Besides the map (2D slice / planar profile / 3D point cloud / Gaussian cube in e/bohr^3) the report checks the
consistency of the folders (grid, cell, mesh cutoff, level of theory, SCF, electron count), gives the
counterpoise-corrected interaction energy from the same runs, statistics of Delta rho and the charge change of
every fragment and atom (nearest-atom partition). --show / --png draw a 4-panel dashboard; --save-gnuplot writes
the .gplot; --save-data a CSV and a JSON. --spin computes the spin-density difference instead (only if the
combined system and every fragment are spin-polarized).""",
        epilog="Example usage:\n"
               "  %(prog)s --dir chargediff_run --png\n"
               "  %(prog)s --dir chargediff_run --profile --axis 2 --save-gnuplot\n"
               "  %(prog)s --dir chargediff_run --axis 2 --pos 7.1 --contour --show --save-data\n"
               "  %(prog)s --dir chargediff_run --3d --iso-min 0.01 --no-cube\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--dir", default="chargediff_run",
                         help="Directory stb-chargediffPrep wrote into -- must contain "
                              "'chargediff_manifest.json' plus the 'combined/' and fragment "
                              "subfolders it lists, each already run in SIESTA. Default: "
                              "chargediff_run (stb-chargediffPrep's own default --output-dir).")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write the data/gnuplot/cube/report files and "
                             "references.bib into (default: current directory). Created if "
                             "it doesn't exist.")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--3d", action="store_true", dest="full_3d", help="Export full 3D volume")
    mode_group.add_argument("--profile", action="store_true",
                        help="Write the planar average of Delta rho along --axis as a 1D profile "
                             "(average taken over the other two axes) instead of a 2D/3D map -- "
                             "useful for bilayers, slabs, and interfaces")
    parser.add_argument("-a", "--axis", type=int, default=None, choices=[0, 1, 2],
                        help="For a slice/profile: the axis normal to the cut plane (slice) or "
                             "the axis the profile varies along (profile) (0=X, 1=Y, 2=Z). Default: Z for the "
                             "slice; for the profile (and the dashboard's profile panel) the lattice axis along "
                             "which the fragments are farthest apart.")
    parser.add_argument("-p", "--pos", type=float,
                        help="Position (Angstrom) of the 2D cut along --axis (default: the plane where the planar "
                             "integral of |Delta rho| is largest; ignored in --profile mode)")
    parser.add_argument("--spin", action="store_true",
                        help="Compute the SPIN density difference instead of the charge density "
                             "difference -- only available if the combined system AND every "
                             "fragment are spin-polarized.")
    parser.add_argument("--iso-min", type=float, default=None, metavar="RHO",
                        help="For --3d: only export/preview points with |Delta rho| >= this "
                             "threshold (e/Ang^3) -- keeps the point cloud a manageable size "
                             "and visually meaningful instead of dominated by near-zero points")
    parser.add_argument("--no-cube", dest="cube", action="store_false", default=True,
                        help="Skip writing Delta rho's full 3D grid as a Gaussian .cube file (in e/bohr^3) -- "
                             "written by default (needs a '*.XV' or 'structure.fdf' in the "
                             "combined folder for the geometry), for real 3D isosurface "
                             "rendering in VESTA/VMD/Avogadro.")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Fix the colorbar/palette lower bound manually (e/Ang^3), e.g. to "
                             "compare several slices on the same scale.")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Fix the colorbar/palette upper bound manually (e/Ang^3).")
    parser.add_argument("--contour", action="store_true",
                        help="Overlay contour lines on a 2D slice map/dashboard (ignored for --3d/--profile)")
    parser.add_argument("--show", "--view", dest="show", action="store_true",
                        help="Open the 4-panel dashboard (slice with atoms, planar profile, charge change per atom, "
                             "summary) on screen with matplotlib. Needs a display. Off by default (--view is the "
                             "older name of this option).")
    parser.add_argument("--png", action="store_true",
                        help="Save the dashboard as <output-dir>/<name>_summary.png. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a .gplot script next to each .dat file this run "
                             "generates. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write <name>_atoms.csv (charge change per atom) and <name>_results.json (every number "
                             "of the report). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-chargediffAnalysis {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Charge Density Difference (Analysis)",
            "Delta rho = rho_combined - sum(rho_fragment_i)",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CHARGE DENSITY DIFFERENCE -- STAGE 2: ANALYSIS", 'bold'))
    print("-" * 60)

    manifest_path = os.path.join(args.dir, CHARGEDIFF_MANIFEST_FILE)
    if not os.path.isfile(manifest_path):
        print(color_text(
            f"[ERROR] '{manifest_path}' not found -- run stb-chargediffPrep with "
            f"--output-dir {args.dir!r} first (this file is where it records the combined/"
            "fragment folder names Stage 2 needs; it's never typed by hand).", 'red'))
        sys.exit(1)
    with open(manifest_path) as f:
        manifest = json.load(f)
    if "combined_folder" not in manifest or "fragment_folders" not in manifest:
        print(color_text(
            f"[ERROR] '{manifest_path}' is from an older stb-chargediffPrep (schema "
            f"{manifest.get('schema', 'unknown')!r}) that didn't yet write a 'combined/' folder "
            "or record 'combined_folder'/'fragment_folders' -- re-run stb-chargediffPrep with "
            f"--output-dir {args.dir!r} to regenerate it before using this Stage 2.", 'red'))
        sys.exit(1)
    combined_folder = manifest["combined_folder"]
    fragment_folders = manifest["fragment_folders"]

    folders = {"combined": os.path.join(args.dir, combined_folder)}
    for name in fragment_folders:
        folders[name] = os.path.join(args.dir, name)

    os.makedirs(args.output_dir, exist_ok=True)

    filenames = {}
    try:
        for name, folder in folders.items():
            filenames[name] = find_one_rho(folder)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    grids, charge_grids = {}, {}
    for name, filename in filenames.items():
        print(f"[INFO] Reading: {color_text(filename, 'cyan')}")
        try:
            sile = sisl.get_sile(filename)
            grids[name] = try_read_net_spin(sile, filename, required=True) if args.spin else read_total(sile)
            charge_grids[name] = grids[name] if not args.spin else read_total(sile)
        except Exception as e:
            print(color_text(f"[FATAL] {filename}: {e}", 'red'))
            sys.exit(1)

    combined_grid_obj = grids["combined"]
    lattice = np.asarray(combined_grid_obj.lattice.cell)
    origin = combined_grid_obj.origin
    combined_shape = combined_grid_obj.grid.shape

    for name in fragment_folders:
        if grids[name].grid.shape != combined_shape:
            print(color_text(
                f"[ERROR] Grid shape mismatch: '{filenames['combined']}' is "
                f"{combined_shape} but '{filenames[name]}' is {grids[name].grid.shape} -- every fragment "
                "must be on the same grid as the combined system to take a difference (see "
                "stb-chargediffPrep's --mesh-cutoff).", 'red'))
            sys.exit(1)
        if not np.allclose(grids[name].lattice.cell, lattice, atol=1e-3):
            print(color_text(
                f"[ERROR] Lattice mismatch: '{filenames[name]}' does not share the combined "
                "system's cell -- Delta rho is only meaningful point-by-point on an identical "
                "cell/grid.", 'red'))
            sys.exit(1)

    diff = combined_grid_obj.grid.copy()
    for name in fragment_folders:
        diff = diff - grids[name].grid

    quantity = "Delta spin density" if args.spin else "Delta charge density"
    prof_axis = 2
    mode = 'profile' if args.profile else ('3d' if args.full_3d else 'slice')
    stem = f"{os.path.basename(os.path.normpath(args.dir)) or 'chargediff'}_chargediff"
    base = os.path.join(args.output_dir, stem)

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warnings = []

    def warn(text):
        warnings.append(text)
        return color_text(f"[WARNING] {text}", 'yellow')

    dv = abs(np.linalg.det(lattice)) / diff.size
    print_dual(color_text("===== STB-CHARGEDIFFANALYSIS REPORT =====", 'magenta'), f_out)

    # ------------------------------------------------------------------ [0]
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Tool           : stb-chargediffAnalysis {VERSION}", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Prep dir       : {args.dir}  ({1 + len(fragment_folders)} folders from {CHARGEDIFF_MANIFEST_FILE})", f_out)
    print_dual(f"Quantity       : {quantity}", f_out)
    print_dual(f"Mode           : {mode}", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    if args.contour and mode != 'slice':
        print_dual(warn("--contour only applies to the 2D slice mode -- ignored."), f_out)
    if args.pos is not None and mode == 'profile':
        print_dual(warn("--pos has no effect in --profile mode -- the profile covers the whole axis."), f_out)

    # ------------------------------------------------------------------ [1]
    print_section("[1] SYSTEM", f_out)
    try:
        atoms, geo_lat = read_atoms(folders, fragment_folders)
    except Exception as e:
        atoms, geo_lat = [], lattice
        print_dual(warn(f"could not read the atoms from the prep folders ({e}): the per-atom analysis is skipped."), f_out)
    runs = {name: read_run(folder) for name, folder in folders.items()}
    if args.axis is not None:
        prof_axis = args.axis
    elif atoms and len(fragment_folders) > 1:
        cents = [np.mean([a[3] for a in atoms if a[2] == k], axis=0) for k in range(len(fragment_folders))
                 if any(a[2] == k for a in atoms)]
        units = [lattice[i] / np.linalg.norm(lattice[i]) for i in range(3)]
        sep = [max(abs(np.dot(c1 - c2, units[i])) for c1 in cents for c2 in cents) for i in range(3)]
        prof_axis = int(np.argmax(sep))
    lengths = [float(np.linalg.norm(v)) for v in lattice]
    spacing = [lengths[i] / combined_shape[i] for i in range(3)]
    print_dual(f"Cell (Ang)     : a={lengths[0]:.4f}  b={lengths[1]:.4f}  c={lengths[2]:.4f}; volume {abs(np.linalg.det(lattice)):.3f} Ang^3", f_out)
    print_dual(f"Real-space grid: {combined_shape[0]} x {combined_shape[1]} x {combined_shape[2]} = {diff.size:,} points, "
               f"spacing {spacing[0]:.4f} {spacing[1]:.4f} {spacing[2]:.4f} Ang, voxel {dv:.3e} Ang^3", f_out)
    versions = {r["version"] for r in runs.values() if r["version"]}
    if versions:
        print_dual(f"SIESTA version : {', '.join(sorted(versions))}", f_out)
    rows = []
    for name in folders:
        if name == "combined":
            n_at = len(atoms) if atoms else "-"
            kind = "all atoms real"
        else:
            k = fragment_folders.index(name)
            n_at = sum(1 for a in atoms if a[2] == k) if atoms else "-"
            kind = f"{n_at} real + {len(atoms) - n_at if atoms else '?'} ghost"
        q_int = float(charge_grids[name].grid.sum() * dv)
        rows.append(([name + "/", str(n_at) if name != "combined" else str(n_at), kind, f"{q_int:.4f}",
                      f"{runs[name]['n_electrons']:.4f}" if runs[name]['n_electrons'] is not None else "-",
                      f"{runs[name]['energy']:.6f}" if runs[name]['energy'] is not None else "-"], None))
    print_table(["Folder", "Atoms", "Basis", "∫rho (e)", "SIESTA electrons", "E_KS (eV)"], rows, f_out)

    # ------------------------------------------------------------------ [2]
    print_section("[2] CONSISTENCY CHECKS", f_out)
    checks = []

    def check(name, ok, detail):
        checks.append((name, bool(ok), detail))
    check("Same real-space grid in every folder", True, f"{combined_shape[0]}x{combined_shape[1]}x{combined_shape[2]}")
    check("Same cell in every folder", True, "within 1e-3 Ang")
    settings = {name: folder_settings(folder) for name, folder in folders.items()}
    ref_name = "combined"
    differing = []
    for key in SETTINGS_KEYS + ["kgrid"]:
        vals = {name: s[0].get(key, "(absent)") for name, s in settings.items()}
        if len(set(vals.values())) > 1:
            differing.append((key, vals))
    check("Same level of theory (functional, basis, smearing, k-grid, D3...)", not differing,
          "identical" if not differing else "differs: " + ", ".join(k for k, _ in differing))
    hashes = {name: tuple(s[1]) for name, s in settings.items()}
    check("Same pseudopotential files (by content) in every folder", len(set(hashes.values())) == 1,
          "identical" if len(set(hashes.values())) == 1 else "different sets")
    meshes = {name: r["mesh_used"] for name, r in runs.items() if r["mesh_used"]}
    check("Same mesh actually used by SIESTA", len(set(meshes.values())) <= 1,
          (f"{next(iter(meshes.values())):.1f} Ry" if meshes and len(set(meshes.values())) == 1 else
           ", ".join(f"{n}: {v:.1f}" for n, v in meshes.items()) or "not found in the logs"))
    not_conv = [n for n, r in runs.items() if r["scf"] is False]
    missing = [n for n, r in runs.items() if r["scf"] is None]
    check("SCF converged in every run", not not_conv and not missing,
          "yes" if not not_conv and not missing else f"not confirmed: {', '.join(not_conv + missing)}")
    if not args.spin:
        sum_frag = sum(float(charge_grids[n].grid.sum() * dv) for n in fragment_folders)
        q_comb = float(charge_grids["combined"].grid.sum() * dv)
        check(f"Electron count of the combined .RHO equals the sum over fragments (< {ELECTRON_TOL:g} e)",
              abs(q_comb - sum_frag) < ELECTRON_TOL, f"{q_comb:.4f} vs {sum_frag:.4f} e")
        for n in folders:
            ne = runs[n]["n_electrons"]
            q = float(charge_grids[n].grid.sum() * dv)
            if ne is not None and abs(q - ne) > 5e-3:
                check(f"Grid integral of {n}'s rho equals SIESTA's electron count", False, f"{q:.4f} vs {ne:.4f} e")
    total = float(diff.sum() * dv)
    if not args.spin:
        check(f"Net Delta rho vanishes (|integral| < {ELECTRON_TOL:g} e)", abs(total) < ELECTRON_TOL, f"{total:+.2e} e")
    rows = [([nm, "OK" if ok else "WARN", detail], None if ok else 'yellow') for nm, ok, detail in checks]
    print_table(["Check", "Status", "Detail"], rows, f_out)
    for key, vals in differing:
        print_dual(warn(f"{SETTINGS_LABELS.get(key, key)} differs between folders: "
                        + ", ".join(f"{n} '{v}'" for n, v in vals.items()) + " -- Delta rho is meaningless unless every "
                        "folder uses the same settings."), f_out)
    for nm, ok, detail in checks:
        if not ok and not nm.startswith("Same level"):
            print_dual(warn(f"{nm}: {detail}."), f_out)

    # ------------------------------------------------------------------ [3]
    print_section("[3] COUNTERPOISE INTERACTION ENERGY", f_out)
    e_int = None
    energies = {n: r["energy"] for n, r in runs.items()}
    if all(v is not None for v in energies.values()):
        e_int = energies["combined"] - sum(energies[n] for n in fragment_folders)
        print_dual("E_int = E(combined) - sum_i E(fragment_i with the other fragments as ghost atoms): every term has the "
                   "same basis, so the basis-set superposition error cancels.", f_out)
        print_dual(f"E(combined)    : {energies['combined']:.6f} eV", f_out)
        for n in fragment_folders:
            print_dual(f"E({n:<12s}): {energies[n]:.6f} eV", f_out)
        print_dual(f"E_int (CP)     : {color_text(f'{e_int:+.4f}', 'bold')} eV = {e_int * EV_TO_KCAL:+.2f} kcal/mol = "
                   f"{e_int * EV_TO_KJ:+.2f} kJ/mol", f_out)
        print_dual(color_text("[NOTE] At the geometry given: it includes no deformation energy of the fragments (the "
                              "fragments are not relaxed apart), so it is the interaction energy, not a binding energy.",
                              'cyan'), f_out)
        if e_int > 0:
            print_dual(warn("E_int is positive: the fragments repel at this geometry (or the level of theory differs)."), f_out)
    else:
        print_dual("Energies not found in the logs (no single *.out per folder): skipped.", f_out)

    # ------------------------------------------------------------------ [4]
    print_section("[4] " + quantity.upper() + ": STATISTICS", f_out)
    moved = float(diff[diff > 0].sum() * dv)
    lost = float(diff[diff < 0].sum() * dv)
    rms = float(np.sqrt(np.mean(diff ** 2)))
    imax = np.unravel_index(int(np.argmax(diff)), diff.shape)
    imin = np.unravel_index(int(np.argmin(diff)), diff.shape)
    stats = {}
    rows = []
    for tag, idx in (("maximum (gain)", imax), ("minimum (loss)", imin)):
        val = float(diff[idx])
        p = voxel_to_cart(idx, diff.shape, lattice, origin)
        near = ""
        if atoms:
            k, dist = nearest_atom(p - np.asarray(origin) if False else p, lattice, np.array([a[3] for a in atoms]), atoms)
            near = f"{atoms[k][1]}{k + 1} ({atoms[k][0]}), {dist:.2f} Ang"
        rows.append(([tag, f"{val:+.5f}", f"{p[0]:.2f} {p[1]:.2f} {p[2]:.2f}", near], None))
        stats[tag] = (val, p)
    print_dual(f"Integral       : {total:+.6f} e", f_out)
    print_dual(f"Charge moved   : {moved:.4f} e gained = {-lost:.4f} e lost (half of the integral of |Delta rho|)", f_out)
    print_dual(f"RMS            : {rms:.5f} e/Ang^3;  |Delta rho| > 1e-3 e/Ang^3 on {100 * np.mean(np.abs(diff) > 1e-3):.2f} % of the cell", f_out)
    print_table(["Extreme", "Value (e/Ang^3)", "Position (Ang)", "Nearest atom"], rows, f_out)
    print_dual(f"Isosurface hint: +/-{0.1 * max(abs(stats['maximum (gain)'][0]), abs(stats['minimum (loss)'][0])) * BOHR_TO_ANG ** 3:.4f} e/bohr^3 "
               "(10 % of the largest |Delta rho|, in the unit of the cube file) is a good first level in VESTA/VMD.", f_out)

    # ------------------------------------------------------------------ [5]
    print_section("[5] CHARGE CHANGE BY FRAGMENT AND ATOM (nearest-atom partition)", f_out)
    per_atom = per_fragment = None
    if atoms:
        apos = np.array([a[3] for a in atoms])
        per_atom = nearest_atom_partition(diff, lattice, apos)
    if per_atom is not None:
        rows = []
        for i, (lab, el, fr, p) in enumerate(atoms):
            rows.append(([str(i + 1), el, lab, fragment_folders[fr] if fr >= 0 else "-", f"{per_atom[i]:+.5f}"], None))
        print_table(["#", "Element", "Label", "Fragment", "Net Delta charge (e)"], rows, f_out)
        per_fragment = [float(sum(per_atom[i] for i, a in enumerate(atoms) if a[2] == k)) for k in range(len(fragment_folders))]
        for nm, q in zip(fragment_folders, per_fragment):
            print_dual(f"{nm:<14s}: {q:+.5f} e ({'gains' if q > 0 else 'loses'} electrons)", f_out)
        print_dual(color_text("[NOTE] The grid is split among the atoms by distance (a geometric, Voronoi-like partition, "
                              "not a Bader analysis): the numbers depend on that choice and are for comparing systems, not "
                              "absolute charges. Positive = the region gains electrons.", 'cyan'), f_out)
        if not args.spin and abs(sum(per_atom)) > 1e-3:
            print_dual(warn(f"the partition does not conserve charge ({sum(per_atom):+.4f} e)."), f_out)
    else:
        print_dual("Skipped (atoms unavailable or the grid is too large for the partition).", f_out)

    # ------------------------------------------------------------------ [6]
    print_section("[6] MAP AND FILES", f_out)
    axis = args.axis if args.axis is not None else (prof_axis if mode == 'profile' else 2)
    n_axis = diff.shape[axis]
    planar_abs = np.abs(diff).sum(axis=tuple(t for t in (0, 1, 2) if t != axis))
    if args.pos is not None:
        slice_idx = int(round((args.pos / np.linalg.norm(lattice[axis])) * n_axis)) % n_axis
        pos_note = "--pos"
    else:
        slice_idx = int(np.argmax(planar_abs))
        pos_note = "plane of the largest planar integral of |Delta rho|"
    slice_pos = slice_idx / n_axis * float(np.linalg.norm(lattice[axis]))
    plane_angle = check_planar_orthogonality(lattice, axis)
    if mode == 'slice':
        print_dual(f"Slice          : plane normal to {'xyz'[axis]} at {slice_pos:.3f} Ang (grid index {slice_idx}; {pos_note})", f_out)
        if abs(plane_angle - 90.0) > NON_ORTHOGONAL_WARN_DEG:
            print_dual(warn(f"the cut plane is skewed ({plane_angle:.1f} deg between its in-plane lattice vectors): gnuplot's "
                            "map may look distorted although the data coordinates are correct."), f_out)
    dat_path = f"{base}_{mode}.dat"
    if mode == 'profile':
        dmin, dmax = write_profile_data_file(diff, lattice, axis, dat_path)
    else:
        dmin, dmax = write_data_file(diff, lattice, origin, dat_path, mode, slice_idx, axis,
                                     iso_min=args.iso_min if mode == '3d' else None)
    cb = compute_colorbar_range(dmin, dmax, True, args.vmin, args.vmax) if mode in ('slice', '3d') else None
    files = [dat_path]
    print_dual(color_text(f"[OK] Data written to '{dat_path}'.", 'green'), f_out)
    if args.save_gnuplot:
        write_gnuplot_script(dat_path, dat_path, mode, quantity, True, axis, slice_pos, plane_angle, cb_range=cb,
                             contour=args.contour and mode == 'slice')
        gp = dat_path.rsplit('.', 1)[0] + ".gplot"
        files.append(gp)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gp}' (gnuplot {gp}).", 'green'), f_out)

    cube_path = None
    if args.cube:
        geo_file, geometry = find_geometry_file(folders["combined"])
        if geometry is None:
            print_dual(color_text(f"[ERROR] --cube requires a geometry file -- neither a '*.XV' nor 'structure.fdf' was "
                                  f"found in '{folders['combined']}'.", 'red'), f_out)
        elif not np.allclose(geometry.cell, lattice, atol=1e-2):
            print_dual(color_text(f"[ERROR] The cell in '{geo_file}' does not match the density grid's own cell -- "
                                  "refusing to write a cube file with mismatched geometry.", 'red'), f_out)
        else:
            cube_path = f"{base}.cube"
            # The grid is in e/Ang^3; sisl writes it unit-less next to a lattice in Bohr, and every cube reader
            # (VESTA, VMD, Avogadro) takes the values as e/bohr^3 -- convert, or isosurface levels and integrated
            # charges come out 6.75x too large.
            combined_grid_obj.grid = diff * BOHR_TO_ANG ** 3
            combined_grid_obj.set_geometry(geometry)
            combined_grid_obj.write(cube_path)
            combined_grid_obj.grid = diff
            files.append(cube_path)
            print_dual(color_text(f"[OK] Cube file saved to '{cube_path}' (values in e/bohr^3, the cube convention; open in "
                                  "VESTA/VMD/Avogadro for real 3D isosurfaces).", 'green'), f_out)

    n_ok = sum(1 for _, ok, _ in checks if ok)
    n_warn = len(checks) - n_ok
    title = f"{manifest.get('title') or os.path.basename(os.path.normpath(args.dir))}"
    if atoms:
        els = [a[1] for a in atoms]
        title = "".join(f"{e}{els.count(e) if els.count(e) > 1 else ''}" for e in sorted(set(els)))
    ctx = {"diff": diff, "lattice": lattice, "origin": origin, "axis": axis, "slice_idx": slice_idx, "slice_pos": slice_pos,
           "cb": cb, "prof_axis": prof_axis, "contour": args.contour, "atoms": atoms, "per_atom": per_atom, "per_fragment": per_fragment,
           "fragment_names": fragment_folders, "quantity": quantity, "title": title, "total": total, "moved": moved,
           "max": (stats["maximum (gain)"][0],), "min": (stats["minimum (loss)"][0],), "rms": rms, "e_int": e_int,
           "n_ok": n_ok, "n_warn": n_warn}
    if args.png or args.show:
        png_path = f"{base}_summary.png" if args.png else None
        produced, shown = make_dashboard(ctx, png_path, args.show, fallback_png=f"{base}_summary.png")
        if produced:
            if not png_path and not shown:        # --show without a display: the figure was saved instead
                png_path = f"{base}_summary.png"
            if png_path:
                files.append(png_path)
                print_dual(color_text(f"[OK] Dashboard saved to '{png_path}'.", 'green'), f_out)
            if shown:
                print_dual("[OK] Dashboard shown on screen.", f_out)
    if args.save_data:
        csv_path, json_path = f"{base}_atoms.csv", f"{base}_results.json"
        with open(csv_path, "w") as f:
            f.write("index,element,label,fragment,net_delta_charge_e\n")
            for i, (lab, el, fr, p) in enumerate(atoms):
                f.write(f"{i + 1},{el},{lab},{fragment_folders[fr] if fr >= 0 else ''},"
                        f"{'' if per_atom is None else format(per_atom[i], '.6f')}\n")
        results = {"tool": f"stb-chargediffAnalysis {VERSION}", "date": datetime.now().isoformat(timespec="seconds"),
                   "directory": args.dir, "quantity": quantity, "grid": list(combined_shape), "integral_e": total,
                   "charge_moved_e": moved, "rms_e_per_ang3": rms,
                   "maximum": {"value": stats["maximum (gain)"][0], "position_ang": stats["maximum (gain)"][1].tolist()},
                   "minimum": {"value": stats["minimum (loss)"][0], "position_ang": stats["minimum (loss)"][1].tolist()},
                   "energies_ev": energies, "e_int_cp_ev": e_int,
                   "per_fragment_e": (dict(zip(fragment_folders, per_fragment)) if per_fragment else None),
                   "per_atom_e": (per_atom.tolist() if per_atom is not None else None),
                   "checks": [{"name": nm, "ok": ok, "detail": dt} for nm, ok, dt in checks], "warnings": warnings}
        with open(json_path, "w") as f:
            json.dump(results, f, indent=1, default=float)
            f.write("\n")
        files += [csv_path, json_path]
        print_dual(color_text(f"[OK] Data files: '{csv_path}', '{json_path}'.", 'green'), f_out)

    # ------------------------------------------------------------------ [7]
    print_section("[7] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(f"[OK] Citations for the methods used in this run written to "
                          f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    # ------------------------------------------------------------------ [8]
    print_section("[8] SUMMARY & FILES", f_out)
    print_dual(f"Status         : {'OK' if not warnings else str(len(warnings)) + ' warning(s)'}; checks {n_ok} ok, {n_warn} warn", f_out)
    if e_int is not None:
        print_dual(f"E_int (CP)     : {e_int:+.4f} eV ({e_int * EV_TO_KCAL:+.2f} kcal/mol)", f_out)
    print_dual(f"Charge moved   : {moved:.4f} e", f_out)
    if per_fragment:
        print_dual("Fragments      : " + "; ".join(f"{nm} {q:+.4f} e" for nm, q in zip(fragment_folders, per_fragment)), f_out)
    for p in files + [os.path.join(args.output_dir, BIB_FILE)]:
        print_dual(f"File           : {p}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    # ------------------------------------------------------------------ [9]
    print_section("[9] WARNINGS", f_out)
    if warnings:
        for k, w in enumerate(warnings, 1):
            print_dual(color_text(f"{k:>2}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Charge density difference computed.\n", 'bold'))


if __name__ == "__main__":
    main()
