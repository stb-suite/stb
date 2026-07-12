#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.1"    

from time import sleep
import argparse
import numpy as np
import sys
import os
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro
# Try to import ASE, required only for .cif files
try:
    import ase.io
except ImportError:
    pass


def parse_poscar(filename):
    """Reads a POSCAR file and returns (lattice, positions, is_cartesian).

    positions are as literally written (fractional under 'Direct', Cartesian
    Angstrom -- pre-scale -- under 'Cartesian'), one row per atom across all
    species, in file order.
    """
    with open(filename, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    scale = float(lines[1])
    vecs = []
    for i in range(2, 5):
        parts = lines[i].split()
        vec = [float(p) for p in parts]
        vecs.append(vec)
    lattice = np.array(vecs) * scale

    counts = [int(x) for x in lines[6].split()]
    n_atoms = sum(counts)
    coord_type = lines[7].strip().lower()
    is_cartesian = coord_type.startswith(('c', 'k'))  # Cartesian, or old VASP 'Kartesian'

    positions = []
    for i in range(8, 8 + n_atoms):
        parts = lines[i].split()
        positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
    positions = np.array(positions)

    return lattice, positions, is_cartesian

def parse_cif(filename):
    """Reads a CIF file using ASE and returns (lattice, fractional_positions, is_cartesian=False)."""
    try:
        atoms = ase.io.read(filename)
        lattice = atoms.get_cell()
        positions = atoms.get_scaled_positions()
        return lattice, positions, False
    except NameError:
        print("Error: The 'ase' library is required to read .cif files.")
        print("Please install it using: pip install ase")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CIF file '{filename}': {e}")
        sys.exit(1)

def parse_fhi(filename):
    """Reads a geometry.in (FHI-aims) file and returns (lattice, positions, is_cartesian).

    is_cartesian reflects whichever atom keyword ('atom_frac' vs 'atom') was
    actually used in the file; mixing both in the same file isn't supported.
    """
    vecs = []
    positions = []
    is_cartesian = False
    with open(filename, 'r') as f:
        for line in f:
            # Remove comments and whitespace
            cleaned_line = line.split('#', 1)[0].strip()
            if not cleaned_line:
                continue
            keyword = cleaned_line.split()[0]

            if keyword == 'lattice_vector':
                parts = cleaned_line.split()
                if len(parts) >= 4:
                    try:
                        # Extract x, y, z values (indices 1, 2, 3)
                        vec = [float(parts[1]), float(parts[2]), float(parts[3])]
                        vecs.append(vec)
                    except ValueError:
                        print(f"Error: 'lattice_vector' line malformed in {filename}: {line}")
                        sys.exit(1)
            elif keyword in ('atom_frac', 'atom'):
                parts = cleaned_line.split()
                if len(parts) >= 4:
                    try:
                        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
                        is_cartesian = (keyword == 'atom')
                    except ValueError:
                        print(f"Error: '{keyword}' line malformed in {filename}: {line}")
                        sys.exit(1)

    # Check if we found exactly 3 vectors
    if len(vecs) != 3:
        print(f"Error: Could not find 3 'lattice_vector' lines in file {filename}.")
        print(f"Found: {len(vecs)}")
        sys.exit(1)

    if not positions:
        print(f"Error: Could not find any 'atom_frac'/'atom' lines in file {filename}.")
        sys.exit(1)

    lattice = np.array(vecs)
    positions = np.array(positions)
    return lattice, positions, is_cartesian

def print_density_recommendation():
    """Prints a friendly k-point density recommendation table."""
    print("\n" + "="*65)
    print("              📐 K-Point Density Recommendation Guide              ")
    print("="*65)
    print("  Density (1/Å⁻¹)        Accuracy Level")
    print("  ---------------      --------------------------")
    print("  0.05 – 0.1           High precision")
    print("  0.10 – 0.30          Medium precision")
    print("  0.30 – 0.50          Low precision")
    print()
    print("  ⚠️  Tip: For most systems, a density between 0.2 and 0.3 is")
    print("     generally accurate enough while keeping cost reasonable.")
    print("="*65 + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="Compute the Monkhorst-Pack grid based on desired k-point density and a structure file."
    )
    parser.add_argument(
        "--density", "-d", type=float, required=True,
        help="Target k-point density (in 1/Å). Example: 0.03"
    )
    parser.add_argument(
        "--file", "-f", type=str, required=True,
        help="Path to the structure file."
    )
    
    parser.add_argument(
        "--type", "-t", type=str, required=True,
        # Changed 'geometry' to 'fhi'
        choices=['poscar', 'cif', 'fhi', 'fdf'],
        help="Type of the structure file. Currently supports: 'poscar', 'cif', 'fhi', 'fdf'."
    )
    parser.add_argument(
        "--vacuum-gap", type=float, default=10.0,
        help="Minimum empty span (in Angstrom) between atoms along an axis, wrapped "
             "periodically, to treat that axis as vacuum-padded (non-periodic) and force "
             "a single k-point there regardless of density. Default: 10.0"
    )


    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kgrid {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Suggested Monkhorst-Pack k-grid from structure :", 'bold'))
    print("-"*60)

    args = parser.parse_args()

    filename = args.file
    file_type = args.type.lower()
    
    try:
        # --- Updated Decision Logic ---
        if file_type == 'cif':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (using ASE)...")
            lattice, positions, is_cartesian = parse_cif(filename)

        elif file_type == 'poscar':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice, positions, is_cartesian = parse_poscar(filename)

        elif file_type == 'fhi': # Changed from 'geometry'
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice, positions, is_cartesian = parse_fhi(filename) # Renamed function call

        elif file_type == 'fdf':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            structure = structure_io.read_fdf(filename)
            lattice = structure.lattice
            positions = np.array([pos for _, pos in structure.atoms])
            is_cartesian = structure.coord_format == 'cartesian'
        # --- End of Update ---

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return

    # Detect vacuum-padded axes from the actual atomic layout, then compute divisions
    try:
        frac_coords = kspace.to_fractional(positions, lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice, args.vacuum_gap)
        divisions = kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], args.density, vacuum_axes)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Print recommendation table
    print_density_recommendation()

    # Print the suggested grid
    print(f"✅ Suggested Monkhorst-Pack grid: {divisions[0]} {divisions[1]} {divisions[2]}\n")

    # Analyze and print dimensionality
    kspace.analyze_dimensionality(vacuum_axes)

if __name__ == "__main__":
    main()
