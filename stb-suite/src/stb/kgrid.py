#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.8.0"    

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
    """Reads a POSCAR file and returns the lattice vectors as a 3x3 numpy array."""
    with open(filename, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    scale = float(lines[1])
    vecs = []
    for i in range(2, 5):
        parts = lines[i].split()
        vec = [float(p) for p in parts]
        vecs.append(vec)
    lattice = np.array(vecs) * scale
    return lattice

def parse_cif(filename):
    """Reads a CIF file using ASE and returns the lattice vectors."""
    try:
        atoms = ase.io.read(filename)
        lattice = atoms.get_cell()
        return lattice
    except NameError:
        print("Error: The 'ase' library is required to read .cif files.")
        print("Please install it using: pip install ase")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CIF file '{filename}': {e}")
        sys.exit(1)

def parse_fhi(filename):
    """Reads a geometry.in (FHI-aims) file and returns the lattice vectors as a 3x3 numpy array."""
    vecs = []
    with open(filename, 'r') as f:
        for line in f:
            # Remove comments and whitespace
            cleaned_line = line.split('#', 1)[0].strip()
            
            if cleaned_line.startswith('lattice_vector'):
                parts = cleaned_line.split()
                if len(parts) >= 4:
                    try:
                        # Extract x, y, z values (indices 1, 2, 3)
                        vec = [float(parts[1]), float(parts[2]), float(parts[3])]
                        vecs.append(vec)
                    except ValueError:
                        print(f"Error: 'lattice_vector' line malformed in {filename}: {line}")
                        sys.exit(1)

    # Check if we found exactly 3 vectors
    if len(vecs) != 3:
        print(f"Error: Could not find 3 'lattice_vector' lines in file {filename}.")
        print(f"Found: {len(vecs)}")
        sys.exit(1)
        
    lattice = np.array(vecs)
    return lattice

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

def analyze_dimensionality(divisions):
    """Analyzes the computed grid to suggest the system's dimensionality."""
    ones_count = divisions.count(1)
    
    print("--- Dimensionality Analysis ---")
    if ones_count == 3:
        # Grid is [1, 1, 1]
        print("System appears to be 0D (e.g., a molecule).")
        print("A 1x1x1 grid (Gamma point) is typically sufficient.")
    elif ones_count == 2:
        # Grid is [N, 1, 1] or [1, N, 1] or [1, 1, N]
        print("System appears to be 1D (e.g., a nanotube or polymer).")
        print("The '1's in the grid correspond to the vacuum-padded directions.")
    elif ones_count == 1:
        # Grid is [N, M, 1] or [N, 1, M] or [1, N, M]
        print("System appears to be 2D (e.g., a slab or surface).")
        print("The '1' in the grid corresponds to the vacuum-padded direction.")
    else:
        # Grid is [N, M, P]
        print("System appears to be 3D (bulk material).")
    print("---------------------------------\n")

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
            lattice = parse_cif(filename)
            
        elif file_type == 'poscar':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice = parse_poscar(filename)

        elif file_type == 'fhi': # Changed from 'geometry'
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice = parse_fhi(filename) # Renamed function call

        elif file_type == 'fdf':
            print(f"ℹ️  Reading file '{filename}' as type '{file_type}' (native method)...")
            lattice = structure_io.lattice_only(filename)
        # --- End of Update ---
            
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return

    # Calculate divisions
    try:
        divisions = kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], args.density)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Print recommendation table
    print_density_recommendation()
    
    # Print the suggested grid
    print(f"✅ Suggested Monkhorst-Pack grid: {divisions[0]} {divisions[1]} {divisions[2]}\n")
    
    # --- New feature ---
    # Analyze and print dimensionality
    analyze_dimensionality(divisions)
    # --- End of new feature ---

if __name__ == "__main__":
    main()
