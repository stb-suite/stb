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
import numpy as np
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

def compute_ecn(structure, mode, output_dir, atoms_position=None):
    out_path = os.path.join(output_dir, "structural_information.dat")
    with open(out_path, "w") as f:

        # Lattice parameters
        lattice = structure.lattice
        print("\n[INFO] Lattice parameters:")
        print(f"   a = {lattice.a:.3f} Å")
        print(f"   b = {lattice.b:.3f} Å")
        print(f"   c = {lattice.c:.3f} Å")
        print(f"   Alpha = {lattice.alpha:.2f}°")
        print(f"   Beta = {lattice.beta:.2f}°")
        print(f"   Gamma = {lattice.gamma:.2f}°")

        f.write("\n[INFO] Lattice parameters:\n")
        f.write(f"   a = {lattice.a:.3f} Å\n")
        f.write(f"   b = {lattice.b:.3f} Å\n")
        f.write(f"   c = {lattice.c:.3f} Å\n")
        f.write(f"   Alpha = {lattice.alpha:.2f}°\n")
        f.write(f"   Beta = {lattice.beta:.2f}°\n")
        f.write(f"   Gamma = {lattice.gamma:.2f}°\n")

        # Lattice vectors
        lattice_vectors = structure.lattice.matrix
        print("\n[INFO] Writing Lattice vectors.")
        f.write("\nLattice vectors:\n")
        for i, vector in enumerate(lattice_vectors):
            f.write(f"   a_{i+1}: {vector[0]}   {vector[1]}   {vector[2]}\n")

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

        # Atomic positions
        pos_atomics = [[i+1, str(site.specie.symbol), site.coords] for i, site in enumerate(structure)]

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

            print("\n[INFO] Calculating effective (weighted) CN per species.")
            f.write("\nEffective Coordination Number (weighted), per species:\n")
            for sp in species_list:
                sp_idx = [i for i, s in enumerate(species_by_index) if s == sp]
                print(f"   {sp} ({len(sp_idx)} atoms):")
                f.write(f" {sp} ({len(sp_idx)} atoms):\n")
                for method_name, values in ecn_results.items():
                    sp_values = [values[i] for i in sp_idx if values[i] is not None]
                    avg = np.mean(sp_values) if sp_values else float('nan')
                    print(f"      {method_name:15}: {avg:.3f}")
                    f.write(f"      {method_name:15}: {avg:.3f}\n")

            f.write("\nEffective Coordination Number (weighted), overall average:\n")
            for method_name, values in ecn_results.items():
                valid = [v for v in values if v is not None]
                avg = np.mean(valid) if valid else float('nan')
                f.write(f" {method_name:15}: {avg:.3f}\n")

        elif mode == "list" and atoms_position:
            for i in atoms_position:
                for method_name, method in methods.items():
                    try:
                        ecn_results[method_name].append(method.get_cn(structure, i-1, use_weights=True))
                    except Exception as e:
                        ecn_results[method_name].append(None)
                        print(f"[WARNING] {method_name} failed for atom {i}: {e}")

            print("\n[INFO] Calculating the effective (weighted) CN for specified atoms.")
            f.write("\nEffective Coordination Number (weighted) for specified atoms:\n")
            for i, atom_index in enumerate(atoms_position):
                f.write(f" Atom {pos_atomics[atom_index-1][0]}:\n")
                f.write(f"   Element: {pos_atomics[atom_index-1][1]}     Cartesian Position: {pos_atomics[atom_index-1][2]}\n")
                for method, values in ecn_results.items():
                    # values[i] is a raw float (~15 significant digits) or
                    # None (that method failed for this atom) -- format it
                    # instead of printing it raw, matching "mean" mode's
                    # .3f, and guard None since f"{None:.3f}" raises.
                    value_str = f"{values[i]:.3f}" if values[i] is not None else "N/A"
                    print(f"      {method:15}: {value_str}")
                    f.write(f"      {method:15}: {value_str}\n")

        # Bond distances via CrystalNN, broken down by species pair -- a
        # single average lumping e.g. Sn-O with any Sn-Sn/O-O neighbors
        # found isn't a physically meaningful number on its own. Note a
        # bond can be counted from both of its atoms' neighbor lists in
        # "mean" mode (n= makes that visible); this only skews the average
        # if the neighbor relationship isn't symmetric between the two.
        print("\n[INFO] Calculating average bond distance per species pair...")
        f.write("\nAverage bond distance, per species pair:\n")

        cnn = CrystalNN()
        distances_by_pair = {}
        all_distances = []

        indices = range(len(structure)) if mode == "mean" else [i - 1 for i in atoms_position]

        for i in indices:
            try:
                neighbors = cnn.get_nn_info(structure, i)
                for neighbor in neighbors:
                    dist = neighbor['site'].distance(structure[i])
                    all_distances.append(dist)
                    pair = tuple(sorted((pos_atomics[i][1], str(neighbor['site'].specie.symbol))))
                    distances_by_pair.setdefault(pair, []).append(dist)
            except Exception as e:
                print(f"[WARNING] Failed to compute distances for atom {i+1}: {e}")

        if all_distances:
            for pair in sorted(distances_by_pair):
                pair_distances = distances_by_pair[pair]
                avg = np.mean(pair_distances)
                label = f"{pair[0]}-{pair[1]}"
                print(f"   {label:10}: {avg:.4f} Å  (n={len(pair_distances)})")
                f.write(f"   {label:10}: {avg:.4f} Å  (n={len(pair_distances)})\n")
            overall = np.mean(all_distances)
            print(f"   {'Overall':10}: {overall:.4f} Å  (n={len(all_distances)})")
            f.write(f"   {'Overall':10}: {overall:.4f} Å  (n={len(all_distances)})\n")
        else:
            print("[WARNING] No distances could be computed.")
            f.write("   No distances could be computed.\n")

        # Atomic positions
        f.write("\nAtomic positions:\n")
        for atom_id, symbol, coords in pos_atomics:
            f.write(f"{atom_id}  {symbol} cartesian position: {coords}\n")

    return out_path

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

    out_path = compute_ecn(structure, args.mode, args.output_dir, atoms_position)

    print(f"\n[INFO] Job complete! Results saved to {out_path}")
    print("\n"+"-"*60)

if __name__ == "__main__":
    main()
