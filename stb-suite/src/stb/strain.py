#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.10.0" # pretty axis-symmetry table (point group + operations + equivalence) replaces the one-line advisory

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
from stb.core import structure_io, kspace, symmetry
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

def print_axis_symmetry_table(requested_axis, groups, point_group, ops, vacuum_axes):
    """Pretty table of the 3 Cartesian axes' symmetry-equivalence groups,
    for a uniaxial strain request -- same visual convention as
    stb-elasticInputs' deformation-direction table (see elastic_inputs.py::
    _print_symmetry_table), reusing core/symmetry.py::operations_summary for
    the reduced point-group/operations line. Always shown (even if nothing
    is equivalent), same "give full information regardless" philosophy.

    Public (not underscore-prefixed): stb_suite.py's interactive wrapper
    imports this directly to build its own pre-run "generate the equivalent
    axis/axes too?" prompt, since stb-strain's CLI only ever strains one
    direction per invocation and the wrapper needs the same table before
    deciding whether to loop over extra directions.

    Returns the list of axis letters equivalent to `requested_axis` (periodic
    ones only -- vacuum-padded axes are excluded even if symmetry-equivalent,
    since they can't physically be strained).
    """
    axis_letters = ['x', 'y', 'z']
    requested_idx = axis_letters.index(requested_axis)
    group = next(g for g in groups if requested_idx in g)

    print()
    print("-" * 60)
    print(color_text(f"AXIS SYMMETRY (uniaxial direction '{requested_axis}')", 'cyan').center(60))
    print("-" * 60)
    print(f"  Detected symmetry : point group {point_group} -- {symmetry.operations_summary(ops)}")
    print("-" * 60)
    print(f"  {'Axis':<8}{'Status':<14}Equivalent to")
    print(f"  {'-' * 52}")
    equivalent = []
    for i, letter in enumerate(axis_letters):
        if i == requested_idx:
            status = color_text("REQUESTED".ljust(14), 'green')
            print(f"  {letter:<8}{status}--")
        elif vacuum_axes[i]:
            status = color_text("VACUUM".ljust(14), 'yellow')
            print(f"  {letter:<8}{status}(not periodic -- can't be strained)")
        elif i in group:
            status = color_text("EQUIVALENT".ljust(14), 'cyan')
            print(f"  {letter:<8}{status}{requested_axis}")
            equivalent.append(letter)
        else:
            print(f"  {letter:<8}{'INDEPENDENT':<14}--")
    print(f"  {'-' * 52}")
    if equivalent:
        letters = ', '.join(equivalent)
        verb = "is" if len(equivalent) == 1 else "are"
        print(f"  {letters} {verb} equivalent to '{requested_axis}' by symmetry -- straining "
              f"{'it' if len(equivalent) == 1 else 'them'} should give the same mechanical "
              "response; you may not need to compute both.")
    else:
        print(f"  No other periodic axis is equivalent to '{requested_axis}' for this point group.")
    print("-" * 60)
    return equivalent


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
    parser.add_argument("--symprec", type=float, default=1e-3,
                       help="Symmetry-detection tolerance (default: 1e-3, pymatgen's own default), "
                            "used only for the informational note about symmetry-equivalent "
                            "directions (uniaxial only) -- never blocks the run.")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                       help="Symmetry angle tolerance in degrees (default: 5.0, pymatgen's own "
                            "default), same use as --symprec.")
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

        # Advisory only (never blocks the run): if the requested uniaxial
        # direction is symmetry-equivalent to another periodic axis, straining
        # that axis instead/too would just repeat the same DFT calculation.
        # Biaxial directions are out of scope for this check -- see
        # core/symmetry.py::equivalent_cartesian_axes docstring.
        if strain_type == 'uniaxial' and norm_dir in ('x', 'y', 'z'):
            try:
                pmg_structure = structure_io.to_pymatgen(structure)
                groups, point_group = symmetry.equivalent_cartesian_axes(
                    pmg_structure, args.symprec, args.angle_tolerance)
                _, ops = symmetry.get_point_group_operations(
                    pmg_structure, args.symprec, args.angle_tolerance)
                print_axis_symmetry_table(norm_dir, groups, point_group, ops, vacuum_axes)
            except Exception:
                pass  # symmetry detection is advisory only, never blocks the tool

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
