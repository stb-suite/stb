#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import os
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro

def determine_strain_type(direction):
    """Determine strain type based on direction input."""
    if len(direction) == 1:  # x, y, z
        return 'uniaxial'
    elif len(direction) == 2:  # xy, xz, yx, yz, zx, zy
        return 'biaxial'
    else:
        raise ValueError("Invalid direction. Use x, y, z for uniaxial or combinations like xy, yz for biaxial.")

def normalize_direction(direction):
    """Normalize direction input (e.g., yx -> xy)."""
    if len(direction) == 2:
        return ''.join(sorted(direction.lower()))
    return direction.lower()

def apply_cartesian_strain(lattice_vectors, strain, direction):
    """
    Apply uniaxial or biaxial strain in Cartesian coordinates.
    
    Args:
        lattice_vectors: 3x3 numpy array of lattice vectors
        strain: strain value (positive for tension, negative for compression)
        direction: strain direction (x, y, z, xy, xz, yz, etc.)
        
    Returns:
        Strained lattice vectors
    """
    # Normalize and validate direction
    direction = normalize_direction(direction)
    valid_directions = {'x', 'y', 'z', 'xy', 'xz', 'yz'}
    if direction not in valid_directions:
        raise ValueError(f"Invalid direction '{direction}'. Use x, y, z, xy, xz, or yz.")

    # Create strain tensor
    strain_tensor = np.zeros((3, 3))
    
    if len(direction) == 1:  # Uniaxial
        if direction == 'x':
            strain_tensor[0, 0] = strain
        elif direction == 'y':
            strain_tensor[1, 1] = strain
        elif direction == 'z':
            strain_tensor[2, 2] = strain
            
    else:  # Biaxial
        if 'x' in direction:
            strain_tensor[0, 0] = strain
        if 'y' in direction:
            strain_tensor[1, 1] = strain
        if 'z' in direction:
            strain_tensor[2, 2] = strain
    
    # Apply strain transformation: new_vec = (I + ε) · vec
    identity = np.eye(3)
    transformation = identity + strain_tensor
    
    # Transform each lattice vector
    strained_vectors = np.dot(transformation, lattice_vectors.T).T
    
    return strained_vectors

def main():
    parser = argparse.ArgumentParser(
        description="Applies strain in Cartesian coordinates to a SIESTA FDF file. "
                   "Type (uniaxial/biaxial) is inferred from direction. "
                   "IMPORTANT: atomic coordinates must be in fractional."
    )
    parser.add_argument("--file", required=True, help="Input FDF file.")
    parser.add_argument("--stdir", required=True, 
                       help="Direction of strain: x, y, z for uniaxial; xy, xz, yz, etc. for biaxial.")
    parser.add_argument("--stmin", type=float, default=0,
                       help="Minimum strain percentage (default: 0). Can be negative for compression.")
    parser.add_argument("--stmax", type=float, default=25,
                       help="Maximum strain percentage (default: 25).")
    parser.add_argument("--step", type=float, default=1,
                       help="Strain step percentage (default: 1).")
    parser.add_argument("-o", "--output-dir", default="strain_runs",
                       help="Directory to write the 'strain_<dir>_<val>' run folders into "
                            "(default: strain_runs). Always used, even for a single direction, "
                            "so multiple stb-strain runs never collide in the same directory.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                       help="Vacuum gap threshold in Ang used to detect which lattice axes are "
                            "periodic vs. vacuum-padded (default: 10.0), same convention as "
                            "stb-kgrid/stb-kpath. Straining a vacuum-padded axis is physically "
                            "meaningless and is refused.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-strain {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    try:
        # Validate strain range
        if args.stmin > args.stmax:
            raise ValueError("Minimum strain cannot be greater than maximum strain.")
        if args.step == 0:
            raise ValueError("Step cannot be zero.")

        # Determine strain type from direction
        strain_type = determine_strain_type(args.stdir)
        norm_dir = normalize_direction(args.stdir)

        if args.intro == True:
            show_intro([
                "Siesta ToolBox Suite",
                "A comprehensive toolkit for SIESTA DFT simulations",
                f"Version {VERSION} | University of Brasilia - 2025",
                "Developed by Dr. Carlos M. O. Bastos"
            ])

        print("\n" + color_text("STRAIN:", 'bold'))
        print("-"*60)

        print(f"[INFO] Detected strain type: {strain_type} (from direction '{args.stdir}')")
        print(f"[INFO] Strain direction: {norm_dir}")
        print(f"[INFO] Strain range: {args.stmin}% to {args.stmax}% with step {args.step}%")
        print("[INFO] Read file")

        structure = structure_io.read_fdf(args.file)
        lattice_vectors = structure_io.raw_lattice_vectors(structure)

        # Detect which axes are periodic vs. vacuum-padded, same convention
        # as stb-kgrid/stb-kpath/stb-mlrelax -- straining a vacuum gap just
        # moves empty space, not the material.
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
        print(f"[INFO] Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}")

        # Only check letters that are actually valid axes here -- an invalid
        # direction (e.g. 'w') is still caught later by apply_cartesian_strain's
        # own validation, with its original, more specific error message.
        axis_index = {'x': 0, 'y': 1, 'z': 2}
        vacuum_requested = [c for c in dict.fromkeys(norm_dir)
                             if c in axis_index and vacuum_axes[axis_index[c]]]
        if vacuum_requested:
            periodic = [c for c, is_vac in zip('xyz', vacuum_axes) if not is_vac]
            periodic_str = ', '.join(periodic) if periodic else 'none (structure is fully isolated / 0D)'
            raise ValueError(
                f"Direction '{norm_dir}' includes vacuum-padded axis/axes "
                f"'{', '.join(vacuum_requested)}' (detected via --vacuum-gap {args.vacuum_gap} Ang); "
                "straining a vacuum gap doesn't correspond to any physical deformation. "
                f"Periodic axis/axes available for this structure: {periodic_str}."
            )

        print("[INFO] Generating FDF files with strain...")

        os.makedirs(args.output_dir, exist_ok=True)

        # Generate strained structures
        n_steps = int(round((args.stmax - args.stmin) / args.step)) + 1
        strain_values = list(np.linspace(args.stmin, args.stmin + (n_steps - 1) * args.step, n_steps) / 100)

        input_basename = os.path.basename(args.file)
        for strain in strain_values:
            # Handle negative strain (compression) in folder name
            strain_prefix = "m" if strain < 0 else ""
            folder = os.path.join(
                args.output_dir, f"strain_{norm_dir}_{strain_prefix}{abs((strain * 100)):.2f}")
            os.makedirs(folder, exist_ok=True)

            new_vectors = apply_cartesian_strain(lattice_vectors, strain, norm_dir)

            output_fdf = os.path.join(folder, input_basename)
            structure_io.rewrite_fdf_lattice(args.file, new_vectors, output_fdf)
            print(f"[OK] Generated: {output_fdf} (strain: {strain*100:.1f}%)")

        print("[INFO] Complete job!")
        print("\n"+"-"*60)
        print(color_text("This lattice has more tension than my last Zoom meeting.\n\n", 'bold'))
    except (ValueError, FileNotFoundError) as e:
        sys.exit(color_text(f"[ERROR] {e}", 'red'))

if __name__ == "__main__":
    main()
