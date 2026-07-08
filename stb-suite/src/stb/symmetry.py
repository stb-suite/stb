#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import argparse
from datetime import datetime
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from stb.core import structure_io
from stb.core.cli import color_text, show_intro

def compute_symmetry(structure, symprec=1e-3, angle_tolerance=5.0):
    """Pure computation -- no printing, no file I/O. Returns a results dict
    that format_report() turns into the console/file report.

    Uses pymatgen's SpacegroupAnalyzer (wrapping spglib) -- the same
    approach already shared by stb-unitcell/stb-fetch via
    core/symmetry.py::reduce_to_unitcell(), instead of calling spglib
    directly with a hand-rolled crystal-system classifier: SpacegroupAnalyzer
    provides the crystal system (and Pearson symbol, lattice type, per-orbit
    Wyckoff multiplicity) natively, so there's nothing left to duplicate.
    """
    sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
    ss = sga.get_symmetrized_structure()
    # site_symmetry_symbols is the per-atom point-group symmetry AT that
    # site (e.g. "m-3m" for a rock-salt ion sitting on a high-symmetry
    # position) -- same original atom order as ss.wyckoff_letters, verified
    # against ds.wyckoffs (which matches ss.wyckoff_letters exactly).
    site_symmetry = sga.get_symmetry_dataset().site_symmetry_symbols

    # ss.wyckoff_symbols is per ORBIT (e.g. "4a" for a 4-atom equivalence
    # class), not per atom -- ss.equivalent_indices maps each orbit to its
    # member atom indices (0-based, original atom order preserved). Rebuild
    # a per-atom "multiplicity+letter" label by looking up which orbit each
    # atom belongs to.
    orbit_of_atom = {}
    for orbit_idx, indices in enumerate(ss.equivalent_indices):
        for i in indices:
            orbit_of_atom[i] = orbit_idx

    sites = []
    for i, site in enumerate(structure):
        orbit_idx = orbit_of_atom[i]
        wyckoff = f"{len(ss.equivalent_indices[orbit_idx])}{ss.wyckoff_letters[i]}"
        sites.append({
            "atom_id": i + 1,
            "species": str(site.specie.symbol),
            "wyckoff": wyckoff,
            "site_symmetry": site_symmetry[i],
            "orbit": orbit_idx + 1,
            "frac_coords": site.frac_coords,
        })

    lattice = structure.lattice
    ops = sga.get_symmetry_operations()

    return {
        "space_group_symbol": sga.get_space_group_symbol(),
        "space_group_number": sga.get_space_group_number(),
        "hall_symbol": sga.get_hall(),
        "hall_number": sga.get_symmetry_dataset().hall_number,
        "point_group": sga.get_point_group_symbol(),
        "crystal_system": sga.get_crystal_system(),
        "lattice_type": sga.get_lattice_type(),
        "pearson_symbol": sga.get_pearson_symbol(),
        "symprec": symprec,
        "angle_tolerance": angle_tolerance,
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "volume": lattice.volume,
            "vectors": lattice.matrix,
        },
        "sites": sites,
        "n_distinct_sites": len(ss.equivalent_indices),
        "orbits": [
            {"wyckoff": ss.wyckoff_symbols[k], "species": str(structure[indices[0]].specie.symbol),
             "n_atoms": len(indices), "example_atom_id": indices[0] + 1,
             "site_symmetry": site_symmetry[indices[0]]}
            for k, indices in enumerate(ss.equivalent_indices)
        ],
        "symmetry_operations": [op.as_xyz_str() for op in ops],
    }

# Default tolerance sweep for --scan-symprec: log-spaced from tight (typical
# DFT-relaxation numerical noise) to loose enough to reveal near-symmetry a
# structure might have lost to that same noise.
DEFAULT_SYMPREC_SCAN = (1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1)

def scan_symprec(structure, symprec_values=DEFAULT_SYMPREC_SCAN, angle_tolerance=5.0):
    """Runs symmetry detection at each symprec in symprec_values (ascending)
    and returns a list of (symprec, space_group_symbol, space_group_number)
    -- lets a caller see whether a structure's *apparent* symmetry is being
    hidden by small numerical noise (e.g. from a DFT relaxation) that a
    tighter tolerance is picking up as a genuine distortion."""
    results = []
    for sp in symprec_values:
        sga = SpacegroupAnalyzer(structure, symprec=sp, angle_tolerance=angle_tolerance)
        results.append((sp, sga.get_space_group_symbol(), sga.get_space_group_number()))
    return results

def _describe_symprec_scan(scan_results):
    """One-line summary of the scan: the first tolerance (if any) at which
    the space group changes from the tightest-tolerance result. Only
    reports the first transition, not every subsequent one, if the space
    group changes more than once across the scanned range."""
    first_symbol, first_number = scan_results[0][1], scan_results[0][2]
    for symprec, symbol, number in scan_results[1:]:
        if number != first_number:
            return (f"Symmetry increases at symprec >= {symprec:g}: "
                    f"{first_symbol} (No. {first_number}) -> {symbol} (No. {number}).")
    return "No change in detected space group across the scanned tolerance range."

# --- Report formatting --------------------------------------------------
_WIDTH = 74

def _rule(char="-"):
    return char * _WIDTH

def format_report(results, source_file, fmt, show_operations=True, symprec_scan=None):
    lat = results["lattice"]
    lines = []
    lines.append(_rule("="))
    lines.append("CRYSTAL SYMMETRY REPORT - STB Suite".center(_WIDTH))
    lines.append(_rule("="))
    lines.append("")
    lines.append(f"Generated        : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Source file      : {source_file}  (format: {fmt})")
    lines.append(f"Symmetry precision: symprec={results['symprec']:g}, "
                 f"angle_tolerance={results['angle_tolerance']:g}°")

    lines.append("")
    lines.append(_rule())
    lines.append("SPACE GROUP")
    lines.append(_rule())
    lines.append(f"Space group      : {results['space_group_symbol']} (No. {results['space_group_number']})")
    lines.append(f"Hall symbol      : {results['hall_symbol']} (Hall No. {results['hall_number']})")
    lines.append(f"Point group      : {results['point_group']}")
    lines.append(f"Crystal system   : {results['crystal_system']}")
    lines.append(f"Lattice type     : {results['lattice_type']}")
    lines.append(f"Pearson symbol   : {results['pearson_symbol']}")
    lines.append(f"Symmetry operations: {len(results['symmetry_operations'])}")

    lines.append("")
    lines.append(_rule())
    lines.append("LATTICE")
    lines.append(_rule())
    lines.append(f"a = {lat['a']:.3f} Å   b = {lat['b']:.3f} Å   c = {lat['c']:.3f} Å")
    lines.append(f"alpha = {lat['alpha']:.2f}°   beta = {lat['beta']:.2f}°   gamma = {lat['gamma']:.2f}°")
    lines.append(f"Volume = {lat['volume']:.3f} Å³")
    lines.append("")
    lines.append("Lattice vectors (Å):")
    for i, vec in enumerate(lat["vectors"]):
        lines.append(f"  a_{i+1}: {vec[0]:12.6f}  {vec[1]:12.6f}  {vec[2]:12.6f}")

    if symprec_scan is not None:
        lines.append("")
        lines.append(_rule())
        lines.append("SYMPREC SENSITIVITY SCAN")
        lines.append(_rule())
        lines.append(f"{'symprec':<12}{'space group'}")
        for symprec, symbol, number in symprec_scan:
            lines.append(f"{symprec:<12g}{symbol} (No. {number})")
        lines.append("")
        lines.append(_describe_symprec_scan(symprec_scan))

    lines.append("")
    lines.append(_rule())
    lines.append(f"SYMMETRICALLY DISTINCT SITES: {results['n_distinct_sites']}")
    lines.append(_rule())
    lines.append(f"{'Wyckoff':<10}{'Species':<10}{'n atoms':<10}{'Site symmetry':<16}{'Example atom'}")
    for orbit in results["orbits"]:
        lines.append(f"{orbit['wyckoff']:<10}{orbit['species']:<10}{orbit['n_atoms']:<10}"
                     f"{orbit['site_symmetry']:<16}{orbit['example_atom_id']}")

    lines.append("")
    lines.append(_rule())
    lines.append("ATOMIC SITES (fractional coordinates)")
    lines.append(_rule())
    lines.append(f"{'Atom':>4}  {'Sp.':<3}  {'Wyckoff':<8}  {'Fractional coordinates':<32}{'Orbit'}")
    for site in results["sites"]:
        fc = site["frac_coords"]
        lines.append(f"{site['atom_id']:>4}  {site['species']:<3}  {site['wyckoff']:<8}  "
                     f"{fc[0]:10.6f}  {fc[1]:10.6f}  {fc[2]:10.6f}  {site['orbit']}")

    if show_operations:
        lines.append("")
        lines.append(_rule())
        lines.append(f"SYMMETRY OPERATIONS ({len(results['symmetry_operations'])}), in x,y,z notation")
        lines.append(_rule())
        for i, op_str in enumerate(results["symmetry_operations"]):
            lines.append(f"{i+1:>4}: {op_str}")

    lines.append(_rule("="))
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(
        description="Analyze crystal symmetry (space group, Wyckoff positions, symmetry "
                     "operations) from a SIESTA structure file.",
        epilog="Example usage:\n"
               "  stb-symmetry --file structure.fdf --format fdf\n"
               "  stb-symmetry --file siesta.STRUCT_OUT --format struct_out --symprec 1e-4",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["fdf", "struct_out"],
                        help="Input file format:\n"
                             "  fdf:        SIESTA structure input (%%block LatticeVectors etc.)\n"
                             "  struct_out: SIESTA post-relaxation output (.STRUCT_OUT)")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry-detection tolerance in Å (default: 1e-3).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry-detection angle tolerance in degrees (default: 5.0).")
    parser.add_argument("--no-operations", dest="show_operations", action="store_false",
                        help="Skip the full symmetry-operations list in the report -- keeps "
                             "just the count in the SPACE GROUP section. Useful for "
                             "high-symmetry structures where the full list is long.")
    parser.add_argument("--scan-symprec", action="store_true",
                        help="Also detect symmetry at a range of tolerances "
                             f"({', '.join(f'{v:g}' for v in DEFAULT_SYMPREC_SCAN)}) and report "
                             "where the space group changes -- reveals symmetry that's present "
                             "but hidden by small numerical noise (e.g. from a DFT relaxation) "
                             "at the --symprec you asked for.")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write symmetry.dat into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-symmetry {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.symprec <= 0:
        parser.error("--symprec must be positive.")
    if args.angle_tolerance <= 0:
        parser.error("--angle-tolerance must be positive.")

    os.makedirs(args.output_dir, exist_ok=True)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CRYSTAL SYMMETRY ANALYSIS:", 'bold'))
    print("-"*60)

    print("\n[INFO] Reading structure file...")
    try:
        structure = structure_io.read_siesta_structure(args.file, args.format)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    print("[INFO] Detecting symmetry...")
    try:
        results = compute_symmetry(structure, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
    except Exception as e:
        print(color_text(f"[ERROR] Symmetry detection failed: {e}", 'red'))
        sys.exit(1)

    symprec_scan = None
    if args.scan_symprec:
        print("[INFO] Scanning symprec sensitivity...")
        try:
            symprec_scan = scan_symprec(structure, angle_tolerance=args.angle_tolerance)
        except Exception as e:
            print(color_text(f"[WARNING] --scan-symprec failed: {e}", 'yellow'))

    report = format_report(results, args.file, args.format,
                           show_operations=args.show_operations, symprec_scan=symprec_scan)
    print("\n" + report)

    out_path = os.path.join(args.output_dir, "symmetry.dat")
    with open(out_path, "w") as f:
        f.write(report)

    print(f"[INFO] Job complete! Results saved to {out_path}")
    print("-"*60)

if __name__ == "__main__":
    main()
