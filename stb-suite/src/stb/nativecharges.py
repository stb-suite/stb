#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Reads and compares the atomic charges SIESTA itself prints at the end of a run -- Mulliken (WriteMullikenPop),
Hirshfeld (Charge.Hirshfeld end) and Voronoi (Charge.Voronoi end) -- from the log: per atom side by side, per species,
the spread between methods, checks (the charges add up to the net charge, population + charge = Z_val, the Sz add up to
the spin moment) and, for a molecule or slab, the dipole of the point charges against SIESTA's own. No new partitioning
is computed: for the iterative Hirshfeld-I see stb-hirshfeldPrep/Ions/Analysis (5.17), for Bader stb-bader (2.6).
"""

from stb import __version__ as VERSION

import io
import os
import re
import sys
import csv
import glob
import argparse
import contextlib
from datetime import datetime

import numpy as np

from stb.core import siesta_log, structure_io
from stb.core.kspace import detect_vacuum_axes
from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections

REPORT_FILE = "stb_nativecharges_report.txt"
DEBYE_PER_E_ANG = 1.0 / 0.20819434
METHODS = [("Mulliken", "WriteMullikenPop 1", siesta_log.get_mulliken_charges),
           ("Hirshfeld", "Charge.Hirshfeld end", siesta_log.get_hirshfeld_charges),
           ("Voronoi", "Charge.Voronoi end", siesta_log.get_voronoi_charges)]


_NEG_ZERO = re.compile(r"(?<![\w.])-(0\.0+)(?![\d])")


def nz(text):
    """'-0.0000' -> '+0.0000' (a rounded tiny negative number reads as a sign that means nothing)."""
    return _NEG_ZERO.sub(r"+\1", text)


def table(headers, rows, buf):
    print_table(headers, [([nz(c) for c in cells], color) for cells, color in rows], buf)


def fdf_value(text, key):
    """First occurrence of an fdf label (SIESTA's first-occurrence-wins), case/./_/- insensitive."""
    norm = key.lower().replace(".", "").replace("_", "").replace("-", "")
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or line.startswith("%"):
            continue
        parts = line.split(None, 1)
        if parts[0].lower().replace(".", "").replace("_", "").replace("-", "") == norm:
            return parts[1].strip() if len(parts) > 1 else ""
    return None


def inspect_folder(path, explicit_label):
    """Everything about one folder's native charges; prints nothing."""
    label = structure_io.find_fdf_system_label(path, explicit_label)
    out_file = siesta_log.find_out_file(path, label)
    r = dict(path=path, label=label, out=out_file, methods={}, warnings=[], basis=None, net_charge=0.0, zval=None,
             status=None, xyz=None, cell=None, vac=None, dipole_siesta=None, spin_moment=None)
    try:
        from stb.status import inspect_folder as status_of
        r["status"] = status_of(path, explicit_label)
    except Exception:
        pass
    for fdf in sorted(glob.glob(os.path.join(path, "*.fdf"))):
        with open(fdf, errors="ignore") as fh:
            text = fh.read()
        if r["basis"] is None:
            r["basis"] = fdf_value(text, "PAO.BasisSize")
        nc = fdf_value(text, "NetCharge")
        if nc:
            try:
                r["net_charge"] = float(nc.split()[0])
            except ValueError:
                pass
    if not out_file:
        r["warnings"].append(f"{path}: no SIESTA log (<label>.out or the only *.out).")
        return r
    for name, flag, reader in METHODS:
        rows = reader(out_file)
        if rows:
            r["methods"][name] = rows
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            r["zval"] = siesta_log.get_zval_from_output(label or "siesta", override_path=out_file)
        except Exception:
            r["zval"] = None
    r["dipole_siesta"] = None
    with open(out_file, errors="ignore") as fh:
        for line in fh:
            if "Electric dipole (Debye)" in line and "=" in line:
                try:
                    r["dipole_siesta"] = np.array([float(v) for v in line.split("=", 1)[1].split()[:3]])
                except ValueError:
                    pass
    r["spin_moment"] = siesta_log.get_spin_moment(out_file)
    xv = os.path.join(path, f"{label}.XV") if label else None
    if xv and os.path.isfile(xv):
        try:
            import sisl
            g = sisl.get_sile(xv).read_geometry()
            r["xyz"], r["cell"] = np.array(g.xyz), np.array(g.cell)
            r["vac"] = detect_vacuum_axes(r["xyz"] @ np.linalg.inv(r["cell"]), r["cell"], 4.0)
        except Exception:
            pass

    w = r["warnings"]
    st = r["status"] or {}
    if st.get("state") and st["state"] != "FINISHED":
        w.append(f"{path}: the run is {st['state']} -- the charges are those of the last step printed.")
    missing = [f"{flag} ({name})" for name, flag, _ in METHODS if name not in r["methods"]]
    if not r["methods"]:
        w.append(f"{path}: no Mulliken, Hirshfeld or Voronoi block in the log -- add {', '.join(missing)} to the .fdf "
                 "and rerun.")
    for name, rows in r["methods"].items():
        total = sum(x["charge"] for x in rows)
        if abs(total + r["net_charge"]) > 0.01:
            w.append(f"{path}: the {name} charges add up to {total:+.4f} e, not {-r['net_charge']:+.4f} (the cell's net "
                     "charge with the sign of the atomic charges).")
    if "Mulliken" in r["methods"] and "Hirshfeld" in r["methods"]:
        flips = [f"{m['sym']}{m['id']}" for m, h in zip(r["methods"]["Mulliken"], r["methods"]["Hirshfeld"])
                 if m["charge"] * h["charge"] < 0 and abs(m["charge"]) > 0.05 and abs(h["charge"]) > 0.05]
        if flips:
            w.append(f"{path}: Mulliken and Hirshfeld disagree in SIGN on {', '.join(flips)} -- a basis-set artefact of "
                     "Mulliken (diffuse/polarization orbitals); prefer Hirshfeld/Voronoi (or Bader, Hirshfeld-I).")
    return r


def point_dipole(r, rows):
    """Dipole (Debye) of the point charges, along the normal of each vacuum axis: positions about the atoms' centre,
    minimum image. None for a bulk cell."""
    if r["xyz"] is None or not any(r["vac"] or []):
        return None
    inv = np.linalg.inv(r["cell"])
    frac = r["xyz"] @ inv
    centre = frac.mean(axis=0)
    d = frac - centre
    d -= np.round(d)
    pos = d @ r["cell"]
    q = np.array([x["charge"] for x in rows])
    return (q[:, None] * pos).sum(axis=0) * DEBYE_PER_E_ANG


def print_folder(r, out, sec, buf):
    st = r["status"] or {}
    info = st.get("info") or {}
    sec("RUN")
    out(f"Folder            : {r['path']}")
    out(f"Log               : {r['out'] or 'none'}" + (f" (SIESTA {info['version']})" if info.get("version") else ""))
    if st:
        out(f"State             : {st.get('state')}, {st.get('type')}")
    if st.get("formula"):
        out(f"Atoms             : {st['natoms']} ({st['formula']})")
    out(f"Basis             : {r['basis'] or 'DZP (SIESTA default)'}")
    if info.get("spin"):
        out(f"Spin              : {info['spin']}")
    if r["net_charge"]:
        out(f"Net charge        : {r['net_charge']:+g} e (NetCharge)")
    for name, flag, _ in METHODS:
        out(f"{name:<18}: " + ("found" if name in r["methods"] else f"not in the log (add {flag})"))
    if not r["methods"]:
        return

    names = list(r["methods"])
    rows0 = r["methods"][names[0]]
    spin = any(x["sz"] is not None for rows in r["methods"].values() for x in rows)
    sec("CHARGES PER ATOM")
    rows = []
    for i, base in enumerate(rows0):
        zv = r["zval"].get(base["sym"]) if r["zval"] else None
        qs = [r["methods"][n][i]["charge"] for n in names if i < len(r["methods"][n])]
        cells = [str(base["id"]), base["sym"], f"{zv:g}" if zv is not None else "-"]
        cells += [f"{q:+.4f}" for q in qs]
        cells += [f"{max(qs) - min(qs):.4f}" if len(qs) > 1 else "-"]
        if spin:
            cells += [f"{r['methods'][n][i]['sz']:+.4f}" if r["methods"][n][i]["sz"] is not None else "-" for n in names]
        rows.append((cells, None))
    headers = ["#", "el", "Z_val"] + [f"q {n}" for n in names] + ["spread"]
    if spin:
        headers += [f"Sz {n}" for n in names]
    table(headers, rows, buf)
    out("q in e (positive = electrons lost); spread = largest minus smallest q over the methods")

    sec("PER SPECIES")
    species = list(dict.fromkeys(x["sym"] for x in rows0))
    rows = []
    for sp in species:
        cells = [sp, str(sum(1 for x in rows0 if x["sym"] == sp))]
        for n in names:
            qs = np.array([x["charge"] for x in r["methods"][n] if x["sym"] == sp])
            cells.append(f"{qs.mean():+.4f}" + (f" +-{qs.std():.4f}" if len(qs) > 1 else ""))
        rows.append((cells, None))
    table(["species", "atoms"] + [f"mean q {n}" for n in names], rows, buf)
    out("+-: spread over the atoms of one species (equivalent atoms should agree)")

    sec("CHECKS")
    rows = []
    for n in names:
        rr = r["methods"][n]
        total = sum(x["charge"] for x in rr)
        dev = None
        if r["zval"]:
            devs = [abs(x["population"] + x["charge"] - r["zval"].get(x["sym"], np.nan)) for x in rr]
            dev = max(devs) if devs else None
        sz = [x["sz"] for x in rr if x["sz"] is not None]
        cells = [n, f"{total:+.6f}", f"{dev:.1e}" if dev is not None and np.isfinite(dev) else "-",
                 f"{sum(sz):.5f}" if sz else "-"]
        rows.append((cells, None))
    table(["method", "sum of q (e)", "max |pop + q - Z_val|", "sum of Sz"], rows, buf)
    expected = -r["net_charge"]
    out(f"Expected          : sum of q = {expected + 0.0:+g} e".replace("-0 e", "0 e") + (f"; sum of Sz = spin moment {r['spin_moment']:.5f}"
                                                            if r["spin_moment"] is not None and spin else ""))
    if len(names) > 1:
        spreads = [max(r["methods"][n][i]["charge"] for n in names) - min(r["methods"][n][i]["charge"] for n in names)
                   for i in range(len(rows0))]
        i = int(np.argmax(spreads))
        out(f"Method dependence : up to {spreads[i]:.3f} e ({rows0[i]['sym']}{rows0[i]['id']}) -- a partial charge is a "
            "convention, not an observable")

    dips = {n: point_dipole(r, r["methods"][n]) for n in names}
    if any(d is not None for d in dips.values()):
        sec("DIPOLE")
        normals = []
        for ax in range(3):
            if r["vac"][ax]:
                nrm = np.cross(r["cell"][(ax + 1) % 3], r["cell"][(ax + 2) % 3])
                nrm /= np.linalg.norm(nrm)
                if np.dot(nrm, r["cell"][ax]) < 0:
                    nrm = -nrm
                normals.append(("abc"[ax], nrm))
        ref = r["dipole_siesta"]
        tot_ref = np.linalg.norm(ref) if ref is not None else None
        rows = []
        for n in names:
            p = dips[n]
            cells = [n] + [f"{np.dot(p, nrm):+.4f}" for _, nrm in normals] + [f"{np.linalg.norm(p):.4f}"]
            cells.append(f"{100 * np.linalg.norm(p) / tot_ref:.0f} %" if tot_ref else "-")
            rows.append((cells, None))
        if ref is not None:
            rows.append((["SIESTA (density)"] + [f"{np.dot(ref, nrm):+.4f}" for _, nrm in normals]
                         + [f"{tot_ref:.4f}", "100 %"], None))
        table(["source"] + [f"along {a} (D)" for a, _ in normals] + ["|p| (D)", "of SIESTA's"], rows, buf)
        out("point charges at the nuclei miss the atoms' own dipoles (lone pairs, polarization): their dipole is not "
            "the molecule's")


def main():
    parser = argparse.ArgumentParser(
        description="Read and compare the atomic charges SIESTA prints in its log -- Mulliken (WriteMullikenPop 1),\n"
                    "Hirshfeld (Charge.Hirshfeld end) and Voronoi (Charge.Voronoi end): per atom side by side, per\n"
                    "species, the spread between methods, checks (sum = net charge, population + q = Z_val, sum of Sz =\n"
                    "spin moment) and, for a molecule or slab, the dipole of the point charges vs SIESTA's. Several\n"
                    "folders (a glob) give one table per species. No new partitioning is computed (Hirshfeld-I: 5.17,\n"
                    "Bader: 2.6).",
        epilog="Example usage:\n"
               "  stb-nativecharges                              # this folder\n"
               "  stb-nativecharges --path nacl                  # another folder\n"
               "  stb-nativecharges --path 'h2o_basis/*'         # a series: one table per species and method\n"
               "  stb-nativecharges --path 'h2o_basis/*' --csv charges.csv",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--path", default=".",
                        help="Directory to inspect, or a glob pattern (e.g. 'strain_*') matching several directories\n"
                             "(default: current directory).")
    parser.add_argument("-l", "--label", default=None,
                        help="SystemLabel override, applied to every directory. Auto-detected per directory from its\n"
                             ".fdf file(s) if not given.")
    parser.add_argument("--details", action="store_true",
                        help="With several folders, also print the full report of each one.")
    parser.add_argument("--csv", default=None, help="Also write every atom's charges (all folders, all methods) as CSV.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also write the report to {REPORT_FILE} (current folder). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nativecharges {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Native Charges",
            "Read and compare SIESTA's own Mulliken/Hirshfeld/Voronoi charges",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    is_glob = any(c in args.path for c in "*?[")
    matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p)) if is_glob else [args.path]
    if not matches:
        parser.error(f"No directory matched '{args.path}'.")
    if not is_glob and not os.path.isdir(args.path):
        parser.error(f"'{args.path}' is not a directory.")
    batch = len(matches) > 1

    buf = io.StringIO() if args.save_report else None
    sec = Sections(buf)

    def out(text=""):
        print_dual(text, buf)

    print("\n" + color_text("NATIVE CHARGES:", 'bold'))
    print("-" * 60)
    out(color_text("===== STB-NATIVECHARGES REPORT =====", 'magenta'))
    sec("RUN METADATA")
    out(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out(f"stb               : {VERSION}")
    out(f"Target            : {args.path}" + (f" ({len(matches)} folders)" if batch else ""))
    out("Methods           : Mulliken (WriteMullikenPop 1), Hirshfeld (Charge.Hirshfeld end), Voronoi "
        "(Charge.Voronoi end) -- as SIESTA printed them")

    results = [inspect_folder(d, args.label) for d in matches]

    if not batch:
        print_folder(results[0], out, sec, buf)
    else:
        sec("CHARGES PER SPECIES")
        names = [n for n, _, _ in METHODS]
        rows = []
        for r in results:
            if not r["methods"]:
                rows.append(([r["path"], r["basis"] or "-", "-", "-"] + ["-"] * len(names), "yellow"))
                continue
            first = next(iter(r["methods"].values()))
            for sp in dict.fromkeys(x["sym"] for x in first):
                cells = [r["path"], r["basis"] or "DZP", sp, str(sum(1 for x in first if x["sym"] == sp))]
                for n in names:
                    if n in r["methods"]:
                        qs = [x["charge"] for x in r["methods"][n] if x["sym"] == sp]
                        cells.append(f"{np.mean(qs):+.4f}")
                    else:
                        cells.append("-")
                rows.append((cells, None))
        table(["folder", "basis", "species", "atoms"] + [f"mean q {n}" for n in names], rows, buf)
        out("mean q in e per atom of the species; '-' = the method is not in that log")
        if args.details:
            for r in results:
                out("")
                out(color_text(f"===== {r['path']} =====", 'magenta'))
                print_folder(r, out, sec, buf)

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["folder", "basis", "atom", "species", "method", "charge", "population", "Sz"])
            for r in results:
                for n, rows in r["methods"].items():
                    for x in rows:
                        wr.writerow([r["path"], r["basis"] or "DZP", x["id"], x["sym"], n, x["charge"], x["population"],
                                     "" if x["sz"] is None else x["sz"]])

    sec("SUMMARY & FILES")
    for n, _, _ in METHODS:
        out(f"{n:<18}: in {sum(1 for r in results if n in r['methods'])} of {len(results)} folder(s)")
    if args.csv:
        out(f"CSV               : {args.csv}")
    if args.save_report:
        out(f"Report            : {REPORT_FILE}")

    sec("WARNINGS")
    allw = [w for r in results for w in r["warnings"]]
    if allw:
        for k, w in enumerate(allw, 1):
            out(color_text(f"  {k}. {w}", 'yellow'))
    else:
        out("  none")
    if buf is not None:
        with open(REPORT_FILE, "w") as fh:
            fh.write(buf.getvalue())

    print("\n" + "-" * 60)
    print(color_text("Native charges report complete.\n", 'bold'))


if __name__ == "__main__":
    main()
