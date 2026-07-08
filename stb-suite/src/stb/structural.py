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
import logging
import warnings
from datetime import datetime
import numpy as np
from scipy.spatial import ConvexHull, QhullError
from pymatgen.analysis.local_env import (
     JmolNN, MinimumDistanceNN, CrystalNN,
    BrunnerNNRelative, EconNN)
from ase.io import read as ase_read
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import color_text, show_intro

def warn_handler(message, category, filename, lineno, file=None, line=None):
    log_message = f"{category.__name__}: {message} (File: {filename}, Line: {lineno})"
    logging.warning(log_message)
    print("[WARNING] Warning detected! Check warnings.log for details.")

def read_structure(path, fmt):
    """Returns a pymatgen Structure for one of this suite's own SIESTA
    formats: .fdf (structure input, via the shared core/structure_io.py
    parser) or .STRUCT_OUT (post-relaxation output, via ASE -- there's no
    shared reader for that format elsewhere in the suite)."""
    if fmt == "fdf":
        fdf_structure = structure_io.read_fdf(path)
        return structure_io.to_pymatgen(fdf_structure)
    atoms = ase_read(path, format="struct_out")
    return AseAtomsAdaptor.get_structure(atoms)

def _safe_mean(values):
    values = [v for v in values if v is not None]
    return float(np.mean(values)) if values else None

def _bond_length_distortion(distances):
    """Baur's bond-length distortion index (%): mean absolute deviation of
    each bond length from this site's own mean bond length. Needs at least
    2 bonds; returns None otherwise (distortion is undefined for a single
    bond)."""
    if len(distances) < 2:
        return None
    d_avg = np.mean(distances)
    if d_avg == 0:
        return None
    return float(100.0 * np.mean([abs(d - d_avg) / d_avg for d in distances]))

def _bond_angle_variance(angles):
    """Bond-angle variance (deg²): variance of this site's L-center-L
    angles around their own mean. Two deliberate departures from Robinson
    et al.'s original (octahedron/tetrahedron) formula, both because the
    coordination number here can be any -- often fractional, "effective"
    -- value, not just 6 or 4:
      1. Uses the site's own observed mean angle as the reference instead
         of a theoretical ideal polyhedral angle (not well-defined for an
         arbitrary/fractional CN).
      2. Includes ALL pairwise ligand-center-ligand angles, not just the
         12 "cis" (~90°) angles Robinson's octahedral formula uses --
         there's no general way to separate "cis" from "trans" pairs
         outside an assumed octahedral topology. For a real 6-coordinate
         site this pulls in the ~180° "trans" angles too, inflating the
         variance well above textbook octahedral BAV values.
    Not directly comparable to literature BAV numbers as a result; still
    a valid, self-consistent distortion metric for comparing sites/species
    within one report. Needs at least 2 angles (3 neighbors); returns None
    otherwise."""
    if len(angles) < 2:
        return None
    theta_avg = np.mean(angles)
    n = len(angles)
    return float(sum((theta - theta_avg) ** 2 for theta in angles) / (n - 1))

def _polyhedron_volume(neighbor_coords):
    """Volume (Å³) of the convex hull of a site's neighbor positions --
    needs at least 4 non-coplanar neighbors to bound a 3D volume; returns
    None otherwise (too few neighbors, or they're coplanar/collinear)."""
    if len(neighbor_coords) < 4:
        return None
    try:
        return float(ConvexHull(np.array(neighbor_coords)).volume)
    except QhullError:
        return None

def compute_ecn(structure, mode, atoms_position=None):
    """Pure computation -- no printing, no file I/O. Returns a results dict
    that format_report() turns into the console/file report."""
    lattice = structure.lattice
    results = {
        "mode": mode,
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "volume": lattice.volume,
            "density": float(structure.density),
            "vectors": lattice.matrix,
        },
    }

    # Coordination-number methods. use_weights=True makes every one of
    # them -- EconNN included, its get_cn() respects the flag exactly
    # like the rest even though its formula (Hoppe's ECoN) is defined
    # as a continuous quantity: verified this atom-by-atom, e.g.
    # use_weights=False gives a plain integer count (6) where
    # use_weights=True gives the real ECoN value (5.98) -- return a
    # genuinely "effective" (continuous, neighbor-weighted) coordination
    # number. pymatgen's NearNeighbors.get_cn() defaults to
    # use_weights=False (a plain integer neighbor count) for all five.
    #
    # JmolNN's weights are on a different scale than the other four: its
    # reference distance is a fixed Jmol bonding-radius lookup table, not
    # this atom's own closest-neighbor distance, so its weight can exceed
    # 1.0 and its resulting CN is typically the outlier among the five.
    methods = {
        "JmolNN": JmolNN(),
        "MinDistNN": MinimumDistanceNN(),
        # weighted_cn=True must match the use_weights=True passed to
        # get_cn() below -- CrystalNN raises ValueError otherwise,
        # unlike the other methods (which don't require a matching
        # constructor flag).
        "CrystalNN": CrystalNN(weighted_cn=True),
        "BrunnerNN": BrunnerNNRelative(),
        "EconNN": EconNN()
    }
    ecn_results = {method: [] for method in methods}

    pos_atomics = [(i+1, str(site.specie.symbol), site.coords) for i, site in enumerate(structure)]
    results["atomic_positions"] = pos_atomics

    if mode == "mean":
        for i in range(len(structure)):
            for method_name, method in methods.items():
                try:
                    ecn_results[method_name].append(method.get_cn(structure, i, use_weights=True))
                except Exception as e:
                    ecn_results[method_name].append(None)
                    print(f"[WARNING] {method_name} failed for atom {i+1}: {e}")

        # Per-species averages -- a single structure-wide average mixes
        # chemically distinct sites (e.g. a cation's coordination with
        # an anion's) into one physically meaningless number. Reported
        # per species, plus an overall figure for reference.
        species_by_index = [pos_atomics[i][1] for i in range(len(structure))]
        species_list = sorted(set(species_by_index))

        cn_per_species = {}
        for sp in species_list:
            sp_idx = [i for i, s in enumerate(species_by_index) if s == sp]
            cn_per_species[sp] = {
                "n_atoms": len(sp_idx),
                "values": {
                    method_name: _safe_mean([values[i] for i in sp_idx])
                    for method_name, values in ecn_results.items()
                },
            }
        results["cn_per_species"] = cn_per_species
        results["cn_overall"] = {
            method_name: _safe_mean(values) for method_name, values in ecn_results.items()
        }

    elif mode == "list" and atoms_position:
        for i in atoms_position:
            for method_name, method in methods.items():
                try:
                    ecn_results[method_name].append(method.get_cn(structure, i-1, use_weights=True))
                except Exception as e:
                    ecn_results[method_name].append(None)
                    print(f"[WARNING] {method_name} failed for atom {i}: {e}")

        cn_per_atom = []
        for k, atom_index in enumerate(atoms_position):
            atom_id, species, position = pos_atomics[atom_index - 1]
            cn_per_atom.append({
                "atom_id": atom_id,
                "species": species,
                "position": position,
                "values": {method: values[k] for method, values in ecn_results.items()},
            })
        results["cn_per_atom"] = cn_per_atom

    # Bond distances/angles/polyhedron shape, all from one CrystalNN
    # neighbor list per atom -- reused for everything below instead of
    # adding yet another per-method comparison (the 5-method table above
    # is specifically about comparing coordination-number *algorithms*;
    # this geometric analysis just needs one consistent neighbor set, the
    # common choice in the literature).
    #
    # Distances/angles are pooled by species pair/triplet across the
    # structure (a single average lumping e.g. Sn-O with any Sn-Sn/O-O
    # found isn't meaningful on its own); a bond or angle can be counted
    # from more than one of its atoms' neighbor lists in "mean" mode (n=
    # makes that visible) -- this only skews the pooled average if the
    # neighbor relationship isn't symmetric between the atoms involved.
    #
    # Bond-length distortion (BLD), bond-angle variance (BAV), and
    # polyhedron volume are inherently per-atom (they characterize one
    # coordination polyhedron), then averaged per species like ECN.
    cnn = CrystalNN()
    distances_by_pair = {}
    all_distances = []
    angles_by_triplet = {}
    all_angles = []
    bld_by_atom = {}
    bav_by_atom = {}
    volume_by_atom = {}

    indices = range(len(structure)) if mode == "mean" else [i - 1 for i in atoms_position]

    for i in indices:
        try:
            neighbors = cnn.get_nn_info(structure, i)
            center_coords = structure[i].coords

            site_distances = []
            for neighbor in neighbors:
                dist = neighbor["site"].distance(structure[i])
                site_distances.append(dist)
                all_distances.append(dist)
                pair = tuple(sorted((pos_atomics[i][1], str(neighbor["site"].specie.symbol))))
                distances_by_pair.setdefault(pair, []).append(dist)
            bld = _bond_length_distortion(site_distances)
            if bld is not None:
                bld_by_atom[i] = bld
            volume = _polyhedron_volume([n["site"].coords for n in neighbors])
            if volume is not None:
                volume_by_atom[i] = volume

            site_angles = []
            for a in range(len(neighbors)):
                for b in range(a + 1, len(neighbors)):
                    v1 = neighbors[a]["site"].coords - center_coords
                    v2 = neighbors[b]["site"].coords - center_coords
                    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    angle_deg = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
                    site_angles.append(angle_deg)
                    all_angles.append(angle_deg)
                    ligand_pair = tuple(sorted((str(neighbors[a]["site"].specie.symbol),
                                                 str(neighbors[b]["site"].specie.symbol))))
                    triplet = (ligand_pair[0], pos_atomics[i][1], ligand_pair[1])
                    angles_by_triplet.setdefault(triplet, []).append(angle_deg)
            bav = _bond_angle_variance(site_angles)
            if bav is not None:
                bav_by_atom[i] = bav
        except Exception as e:
            print(f"[WARNING] Failed to compute distances/angles for atom {i+1}: {e}")

    results["bond_distances_by_pair"] = {
        pair: (float(np.mean(d)), len(d)) for pair, d in distances_by_pair.items()
    }
    results["bond_distances_overall"] = (
        (float(np.mean(all_distances)), len(all_distances)) if all_distances else None
    )
    results["bond_angles_by_triplet"] = {
        triplet: (float(np.mean(a)), len(a)) for triplet, a in angles_by_triplet.items()
    }
    results["bond_angles_overall"] = (
        (float(np.mean(all_angles)), len(all_angles)) if all_angles else None
    )

    if mode == "mean":
        distortion_per_species = {}
        for sp in species_list:
            sp_idx = [i for i, s in enumerate(species_by_index) if s == sp]
            distortion_per_species[sp] = {
                "n_atoms": len(sp_idx),
                "bld": _safe_mean([bld_by_atom.get(i) for i in sp_idx]),
                "bav": _safe_mean([bav_by_atom.get(i) for i in sp_idx]),
                "volume": _safe_mean([volume_by_atom.get(i) for i in sp_idx]),
            }
        results["distortion_per_species"] = distortion_per_species
    else:
        distortion_per_atom = []
        for atom_index in atoms_position:
            i = atom_index - 1
            distortion_per_atom.append({
                "atom_id": pos_atomics[i][0],
                "species": pos_atomics[i][1],
                "bld": bld_by_atom.get(i),
                "bav": bav_by_atom.get(i),
                "volume": volume_by_atom.get(i),
            })
        results["distortion_per_atom"] = distortion_per_atom

    return results

# --- Report formatting --------------------------------------------------
_WIDTH = 74

def _rule(char="-"):
    return char * _WIDTH

def _fmt(value, width=7, prec=3):
    return f"{value:{width}.{prec}f}" if value is not None else "N/A".rjust(width)

def format_report(results, source_file, fmt):
    lat = results["lattice"]
    lines = []
    lines.append(_rule("="))
    lines.append("STRUCTURAL PROPERTIES REPORT - STB Suite".center(_WIDTH))
    lines.append(_rule("="))
    lines.append("")
    lines.append(f"Generated        : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Source file      : {source_file}  (format: {fmt})")
    lines.append(f"Analysis mode    : {results['mode']}")

    lines.append("")
    lines.append(_rule())
    lines.append("LATTICE")
    lines.append(_rule())
    lines.append(f"a = {lat['a']:.3f} Å   b = {lat['b']:.3f} Å   c = {lat['c']:.3f} Å")
    lines.append(f"alpha = {lat['alpha']:.2f}°   beta = {lat['beta']:.2f}°   gamma = {lat['gamma']:.2f}°")
    lines.append(f"Volume = {lat['volume']:.3f} Å³   Density = {lat['density']:.3f} g/cm³")
    lines.append("")
    lines.append("Lattice vectors (Å):")
    for i, vec in enumerate(lat["vectors"]):
        lines.append(f"  a_{i+1}: {vec[0]:12.6f}  {vec[1]:12.6f}  {vec[2]:12.6f}")

    if results["mode"] == "mean":
        lines.append("")
        lines.append(_rule())
        lines.append("EFFECTIVE COORDINATION NUMBER (weighted), PER SPECIES")
        lines.append(_rule())
        for sp, data in results["cn_per_species"].items():
            lines.append(f"{sp} ({data['n_atoms']} atoms):")
            for method, avg in data["values"].items():
                lines.append(f"  {method:<15}: {_fmt(avg)}")
            lines.append("")
        lines.append("Overall average:")
        for method, avg in results["cn_overall"].items():
            lines.append(f"  {method:<15}: {_fmt(avg)}")

    elif results["mode"] == "list":
        lines.append("")
        lines.append(_rule())
        lines.append("EFFECTIVE COORDINATION NUMBER (weighted), SPECIFIED ATOMS")
        lines.append(_rule())
        for atom in results["cn_per_atom"]:
            pos = atom["position"]
            lines.append(f"Atom {atom['atom_id']} ({atom['species']}), position: "
                         f"{pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
            for method, val in atom["values"].items():
                lines.append(f"  {method:<15}: {_fmt(val)}")
            lines.append("")

    lines.append(_rule())
    lines.append("AVERAGE BOND DISTANCE, PER SPECIES PAIR")
    lines.append(_rule())
    if results["bond_distances_overall"] is not None:
        for pair in sorted(results["bond_distances_by_pair"]):
            avg, n = results["bond_distances_by_pair"][pair]
            label = f"{pair[0]}-{pair[1]}"
            lines.append(f"  {label:<10}: {avg:.4f} Å  (n={n})")
        overall_avg, overall_n = results["bond_distances_overall"]
        lines.append(f"  {'Overall':<10}: {overall_avg:.4f} Å  (n={overall_n})")
    else:
        lines.append("  No distances could be computed.")

    lines.append("")
    lines.append(_rule())
    lines.append("AVERAGE BOND ANGLE, PER LIGAND-CENTER-LIGAND TRIPLET")
    lines.append(_rule())
    if results["bond_angles_overall"] is not None:
        for triplet in sorted(results["bond_angles_by_triplet"]):
            avg, n = results["bond_angles_by_triplet"][triplet]
            label = f"{triplet[0]}-{triplet[1]}-{triplet[2]}"
            lines.append(f"  {label:<12}: {avg:7.3f}°  (n={n})")
        overall_avg, overall_n = results["bond_angles_overall"]
        lines.append(f"  {'Overall':<12}: {overall_avg:7.3f}°  (n={overall_n})")
    else:
        lines.append("  No angles could be computed (fewer than 2 neighbors per site).")

    lines.append("")
    lines.append(_rule())
    lines.append("COORDINATION POLYHEDRON DISTORTION")
    lines.append(_rule())
    lines.append("BLD: bond-length distortion (%, Baur). BAV: bond-angle variance (deg²),")
    lines.append("from ALL ligand-center-ligand angles (not just Robinson's 12 'cis' angles")
    lines.append("for an ideal octahedron -- there's no general cis/trans split for an")
    lines.append("arbitrary/fractional CN), around the site's own mean angle (not a")
    lines.append("theoretical ideal) -- not directly comparable to textbook BAV values, but")
    lines.append("self-consistent for comparing sites/species within this report. Volume:")
    lines.append("convex hull of the neighbor positions (Å³, needs >= 4 non-coplanar neighbors).")
    lines.append("")
    if results["mode"] == "mean":
        for sp, data in results["distortion_per_species"].items():
            lines.append(f"{sp} ({data['n_atoms']} atoms):")
            lines.append(f"  BLD   : {_fmt(data['bld'])} %")
            lines.append(f"  BAV   : {_fmt(data['bav'])} deg²")
            lines.append(f"  Volume: {_fmt(data['volume'])} Å³")
            lines.append("")
    else:
        for atom in results["distortion_per_atom"]:
            lines.append(f"Atom {atom['atom_id']} ({atom['species']}):")
            lines.append(f"  BLD   : {_fmt(atom['bld'])} %")
            lines.append(f"  BAV   : {_fmt(atom['bav'])} deg²")
            lines.append(f"  Volume: {_fmt(atom['volume'])} Å³")
            lines.append("")

    lines.append(_rule())
    lines.append("ATOMIC POSITIONS (Cartesian, Å)")
    lines.append(_rule())
    for atom_id, species, coords in results["atomic_positions"]:
        lines.append(f"{atom_id:>4}  {species:<3}  {coords[0]:12.6f}  {coords[1]:12.6f}  {coords[2]:12.6f}")

    lines.append(_rule("="))
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(
        description="Compute ECN and structural properties from a SIESTA structure file.",
        epilog="Example usage:\n"
               "  stb-structural --file structure.fdf --format fdf --mode mean\n"
               "  stb-structural --file siesta.STRUCT_OUT --format struct_out --mode list --list 1,4,5",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["fdf", "struct_out"],
                        help="Input file format:\n"
                             "  fdf:        SIESTA structure input (%%block LatticeVectors etc.)\n"
                             "  struct_out: SIESTA post-relaxation output (.STRUCT_OUT)")
    parser.add_argument("--mode", choices=["list", "mean"], required=True, help="Calculation mode: list or mean")
    parser.add_argument("--list", type=str, help="List of atom indices (comma-separated, 1-based). Example: 1,4,5,7 - Required for 'list' mode")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write structural_information.dat and warnings.log into "
                             "(default: current directory). Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-structural {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.mode == "list" and not args.list:
        parser.error("--list is required when --mode is 'list'")

    os.makedirs(args.output_dir, exist_ok=True)

    # Configure logger for warnings (done here, not at module level, so importing
    # stb.structural has no side effect of creating warnings.log on disk).
    # filemode='w' so stale warnings from a previous run in the same directory
    # don't linger forever (logging.basicConfig defaults to append mode).
    logging.basicConfig(filename=os.path.join(args.output_dir, "warnings.log"),
                         level=logging.WARNING, format="%(message)s", filemode='w')
    warnings.showwarning = warn_handler

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STRUCTURAL PROPERTIES:", 'bold'))
    print("-"*60)

    print("\n[INFO] Reading structure file...")
    atoms_position = list(map(int, args.list.strip('[]').split(','))) if args.list else None

    try:
        structure = read_structure(args.file, args.format)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    if atoms_position:
        invalid = [i for i in atoms_position if i < 1 or i > len(structure)]
        if invalid:
            print(color_text(f"[ERROR] Atom index/indices {invalid} out of range "
                              f"(structure has {len(structure)} atoms, 1-based).", 'red'))
            sys.exit(1)

    print("[INFO] Computing coordination numbers and bond distances...")
    results = compute_ecn(structure, args.mode, atoms_position)

    report = format_report(results, args.file, args.format)
    print("\n" + report)

    out_path = os.path.join(args.output_dir, "structural_information.dat")
    with open(out_path, "w") as f:
        f.write(report)

    print(f"[INFO] Job complete! Results saved to {out_path}")
    print("-"*60)

if __name__ == "__main__":
    main()
