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

        # ECN Methods
        methods = {
            "JmolNN": JmolNN(),
            "MinDistNN": MinimumDistanceNN(),
            "CrystalNN": CrystalNN(),
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
                        ecn_results[method_name].append(method.get_cn(structure, i))
                    except Exception as e:
                        ecn_results[method_name].append(None)
                        print(f"[WARNING] {method_name} failed for atom {i+1}: {e}")

            ecn_avg = {method: np.nanmean([v for v in values if v is not None]) for method, values in ecn_results.items()}
            print("\n[INFO] Calculating the average ECN.")
            f.write("\nAverage ECN:\n")
            for method, value in ecn_avg.items():
                f.write(f"{method:15}: {value:.2f}\n")

        elif mode == "list" and atoms_position:
            for i in atoms_position:
                for method_name, method in methods.items():
                    try:
                        ecn_results[method_name].append(method.get_cn(structure, i-1))
                    except Exception as e:
                        ecn_results[method_name].append(None)
                        print(f"[WARNING] {method_name} failed for atom {i}: {e}")

            print("\n[INFO] Calculating the ECN for specified atoms.")
            f.write("\nECN for specified atoms:\n")
            for i, atom_index in enumerate(atoms_position):
                f.write(f" Atom {pos_atomics[atom_index-1][0]}:\n")
                f.write(f"   Element: {pos_atomics[atom_index-1][1]}     Cartesian Position: {pos_atomics[atom_index-1][2]}\n")
                for method, values in ecn_results.items():
                    print(f"      {method:15}: {values[i]}")
                    f.write(f"      {method:15}: {values[i]}\n")

        # Average bond distance calculation using CrystalNN
        print("\n[INFO] Calculating average bond distance...")
        f.write("\n[INFO] Average bond distance:\n")

        cnn = CrystalNN()
        distances = []

        indices = range(len(structure)) if mode == "mean" else [i - 1 for i in atoms_position]

        for i in indices:
            try:
                neighbors = cnn.get_nn_info(structure, i)
                for neighbor in neighbors:
                    dist = neighbor['site'].distance(structure[i])
                    distances.append(dist)
            except Exception as e:
                print(f"[WARNING] Failed to compute distances for atom {i+1}: {e}")

        if distances:
            avg_distance = np.mean(distances)
            print(f"   Average bond distance: {avg_distance:.4f} Å")
            f.write(f"   Average bond distance: {avg_distance:.4f} Å\n")
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
