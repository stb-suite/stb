#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

from time import sleep
import numpy as np
from ase.cell import Cell
from ase.dft.kpoints import parse_path_string
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro
import sys
import os
import argparse

def parse_fdf_to_structure(fdf_file="struct.fdf"):
    """Reads a SIESTA .fdf file and returns its FdfStructure."""
    try:
        return structure_io.read_fdf(fdf_file)
    except FileNotFoundError:
        print(color_text(f"Error: File '{fdf_file}' not found.", 'red'))
        return None
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        return None


def write_siesta_kpath_file(kpoints_dict, path_segments, num_points=50, output_filename="kpath_bs.fdf"):
    """
    Writes the k-path to a SIESTA-formatted FDF file for band structure.

    This function now writes ALL suggested path segments, even if disjointed.
    """

    # 1. Determine the full sequence of k-points in order
    if not path_segments:
        print(f"  {color_text('Warning: No k-path segments found. File will not be written.', 'yellow')}")
        return

    # --- START OF CORRECTION ---
    # Rebuilds the 'path_sequence' based on the correct data structure
    # (which is a list of paths, e.g: [['A', 'B', 'C'], ['D', 'E']])

    path_sequence = []

    for i, segment_list in enumerate(path_segments): # ex: segment_list = ['\Gamma', 'Y', 'H']
        if not segment_list:
            continue # Skip if segment is empty

        if i == 0:
            # For the first path, simply add all points
            path_sequence.extend(segment_list)
        else:
            # For subsequent paths
            last_point_in_sequence = path_sequence[-1]
            first_point_of_new_segment = segment_list[0]

            if last_point_in_sequence == first_point_of_new_segment:
                # The path is continuous (e.g., ends in Z, starts in Z)
                # Add all points, *except* the first one (which is duplicated)
                path_sequence.extend(segment_list[1:])
            else:
                # The path is disjoint (it's a "jump", e.g., ends in H_1, starts in M)
                print(f"  {color_text(f'Note: Disjointed path detected (jump from {last_point_in_sequence} to {first_point_of_new_segment}). Including full segment.', 'cyan')}")
                # Add ALL points from the new segment
                path_sequence.extend(segment_list)

    # --- END OF CORRECTION ---

    if not path_sequence:
        print(f"  {color_text('Error: Could not determine path sequence.', 'red')}")
        return

    # This line will now print the full path, including jumps
    path_str_display = '-'.join([r'\Gamma' if p == 'GAMMA' else p for p in path_sequence])
    print(f"  {color_text('SIESTA Path to be written:', 'cyan')} {path_str_display}")

    try:
        with open(output_filename, 'w') as f:
            f.write("### BANDS\n")
            f.write(" BandLinesScale  ReciprocalLatticeVectors\n\n")
            f.write("%block BandLines\n")

            # 2. Write the first point (with '1' point)
            first_label = path_sequence[0]
            first_coords = kpoints_dict[first_label]
            coord_str = " ".join([f"{c:14.10f}" for c in first_coords])
            # Fix 'GAMMA' to '\Gamma' label in file
            f_label = r'\Gamma' if first_label == 'GAMMA' else first_label
            f.write(f"1   {coord_str}   {f_label}\n")

            # 3. Write the rest of the points
            for label in path_sequence[1:]:
                coords = kpoints_dict[label]
                coord_str = " ".join([f"{c:14.10f}" for c in coords])
                # Fix 'GAMMA' to '\Gamma' label in file
                f_label = r'\Gamma' if label == 'GAMMA' else label
                # Use the number of points (ex: 50) for all following segments
                f.write(f"{num_points}   {coord_str}   {f_label}\n")

            f.write("%endblock BandLines\n")

        print(f"  {color_text('Success:', 'green')} SIESTA file '{color_text(output_filename, 'bold')}' has been created.")

    except KeyError as e:
        # This error happens if a label in the path (ex: '\Gamma') has a
        # different name in the dictionary (ex: 'GAMMA')

        # Tries to substitute '\Gamma' with 'GAMMA' if 'GAMMA' exists
        if str(e) == r"'\Gamma'" and 'GAMMA' in kpoints_dict:
            print("  " + color_text("Warning: Found '\\Gamma' label, attempting to use 'GAMMA' internally.", 'yellow'))
            # Recreate the sequence, substituting \Gamma with GAMMA
            path_sequence_fixed = ['GAMMA' if label == r'\Gamma' else label for label in path_sequence]
            # Try to write the file AGAIN
            write_siesta_kpath_file_fixed(kpoints_dict, path_sequence_fixed, num_points, output_filename)
        else:
             print(f"  {color_text(f'Error writing file: K-point label {e} found in path but not in k-points list.', 'red')}")
             print(f"  {color_text(f'Available labels: {list(kpoints_dict.keys())}', 'red')}")

    except Exception as e:
        print(f"  {color_text(f'Error writing {output_filename}: {e}', 'red')}")

# --- HELPER FUNCTION TO FIX LABEL NAMES ---
def write_siesta_kpath_file_fixed(kpoints_dict, path_sequence, num_points, output_filename):
    """
    Emergency "helper" function to write the file if the KeyError
    for '\Gamma' vs 'GAMMA' occurs.
    """
    print("  " + color_text("Retrying file write with 'GAMMA' label...", 'cyan'))
    try:
        with open(output_filename, 'w') as f:
            f.write("### BANDS\n")
            f.write(" BandLinesScale  ReciprocalLatticeVectors\n\n")
            f.write("%block BandLines\n")

            first_label = path_sequence[0]
            first_coords = kpoints_dict[first_label]
            coord_str = " ".join([f"{c:14.10f}" for c in first_coords])
            f_label = r'\Gamma' if first_label == 'GAMMA' else first_label
            f.write(f"1   {coord_str}   {f_label}\n")

            for label in path_sequence[1:]:
                coords = kpoints_dict[label]
                coord_str = " ".join([f"{c:14.10f}" for c in coords])
                f_label = r'\Gamma' if label == 'GAMMA' else label
                f.write(f"{num_points}   {coord_str}   {f_label}\n")

            f.write("%endblock BandLines\n")

        print(f"  {color_text('Success:', 'green')} SIESTA file '{color_text(output_filename, 'bold')}' has been created (with label fix).")
    except Exception as e:
        print(f"  {color_text(f'Error during retry: {e}', 'red')}")
# --- END OF SIESTA WRITE FUNCTIONS ---


def get_kpath_from_structure(fdf_structure, vacuum_gap=10.0, eps=0.0002):
    """
    Takes an FdfStructure, detects its dimensionality from vacuum-padded axes
    (see kspace.detect_vacuum_axes), and prints/writes its high-symmetry
    k-path using ASE's dimension-aware Bravais-lattice path finder
    (ase.cell.Cell.bandpath), then writes it to 'kpath_bs.fdf'.
    """
    try:
        lattice = fdf_structure.lattice
        positions = np.array([pos for _, pos in fdf_structure.atoms])
        is_cartesian = fdf_structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice, vacuum_gap)
        dimension = 3 - sum(vacuum_axes)

        pymatgen_structure = structure_io.to_pymatgen(fdf_structure)

        print(color_text("--- 1. Structure Analysis ---", 'bold'))
        print(f"  {color_text('Formula:', 'cyan')} {pymatgen_structure.composition.reduced_formula}")

        if dimension == 0:
            print(f"  {color_text('Dimensionality:', 'cyan')} {color_text('0D (isolated molecule)', 'bold')}")
            print("\n" + color_text(
                "No periodic direction detected (every axis is vacuum-padded) -- a k-path\n"
                "is not physically meaningful for an isolated molecule. Skipping file generation.",
                'yellow'))
            return

        pbc = tuple(not v for v in vacuum_axes)
        cell = Cell(lattice)
        bravais = cell.get_bravais_lattice(eps=eps, pbc=pbc)

        dim_label = {1: "1D", 2: "2D", 3: "3D"}[dimension]
        print(f"  {color_text('Dimensionality:', 'cyan')} {color_text(dim_label, 'bold')}")
        print(f"  {color_text('Bravais Lattice:', 'cyan')} {color_text(f'{bravais.longname} ({bravais.name})', 'bold')}")
        if dimension < 3:
            print(color_text(
                "  Note: this reflects only the periodic axes (vacuum-padded axes excluded);\n"
                "  it is not a 3D space group.", 'yellow'))
        print(f"  {color_text('Using vacuum gap threshold:', 'yellow')} {vacuum_gap} Ang")

        bp = cell.bandpath(pbc=pbc, npoints=0, eps=eps)

        # ASE labels the Gamma point 'G'; translate to 'GAMMA' so the writer's
        # existing '\Gamma' display/formatting logic keeps working unchanged.
        kpoints_dict = {('GAMMA' if label == 'G' else label): coords
                        for label, coords in bp.special_points.items()}
        path_segments = [[('GAMMA' if label == 'G' else label) for label in segment]
                         for segment in parse_path_string(bp.path)]

        print("\n" + color_text("--- 2. High-Symmetry K-Points (Fractional Coordinates) ---", 'bold'))
        for label, coords in kpoints_dict.items():
            coord_str = ", ".join([f"{c:8.5f}" for c in coords])
            display_label = r'\Gamma' if label == 'GAMMA' else label
            print(f"  {color_text(f'{display_label:<5}:', 'green')} ({coord_str})")

        print("\n" + color_text("--- 3. Suggested K-Path ---", 'bold'))
        path_segments_str_list = []
        for segment_list in path_segments:
            display_segment = [r'\Gamma' if label == 'GAMMA' else label for label in segment_list]
            path_segments_str_list.append("-".join(display_segment))

        print(f"  Path: {color_text(color_text(' | '.join(path_segments_str_list), 'bold'), 'green')}")

        # --- METHODOLOGY ---
        print("\n" + color_text("Methodology:", 'cyan'))
        print(f"  The path follows ASE's dimension-aware Bravais-lattice convention,")
        print(f"  which extends the Setyawan & Curtarolo (2010) 3D scheme to 1D/2D")
        print(f"  systems via their periodic-axes mask (pbc).")
        # --- END OF METHODOLOGY ---

        num_points_per_segment = 50

        # --- FILE GENERATION ---
        print("\n" + color_text("--- 4. Output File Generation ---", 'bold'))

        write_siesta_kpath_file(
            kpoints_dict,
            path_segments,
            num_points_per_segment,
            "kpath_bs.fdf" # Output filename
        )
        # --- END OF GENERATION ---

    except Exception as e:
        print("\n" + color_text(f"An error occurred while analyzing the structure: {e}", 'red'))
        print(color_text(f"It might be necessary to adjust the 'vacuum_gap' or 'eps' parameters (current: {vacuum_gap} Ang, {eps}).", 'yellow'))
        print(color_text("This often happens if the input structure is distorted or the FDF parse failed.", 'yellow'))


def main():

    # 1. Configure ArgParse
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Finds the high-symmetry k-path for FDF (SIESTA) files.", 'bold')}
Dimension-aware (0D/1D/2D/3D) via ASE's Bravais-lattice path finder.""",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Argument 1: Filename (now a required flag)
    parser.add_argument(
        "-f", "--file",
        dest="filename", # Stores the value in 'args.filename'
        type=str,
        required=True,  # Makes this flag mandatory
        help="The structure file name (e.g., struct.fdf)."
    )

    # Argument 2: Precision (optional)
    parser.add_argument(
        "-p", "--prec",
        dest="eps",
        type=float,
        default=0.0002,
        help="Tolerance (ase.cell.Cell's 'eps') for Bravais-lattice/symmetry detection.\n"
             "Default: 0.0002"
    )

    # Argument 3: Vacuum gap threshold (optional)
    parser.add_argument(
        "--vacuum-gap",
        type=float,
        default=10.0,
        help="Minimum empty span (in Angstrom) between atoms along an axis, wrapped "
             "periodically, to treat that axis as vacuum-padded (non-periodic) when "
             "detecting the system's dimensionality. Default: 10.0"
    )

    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kpath {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Suggested k-path from structure:", 'bold'))
    print("-"*60)


    # Read the arguments provided on the command line
    args = parser.parse_args()

    # 2. Process the arguments

    filename = args.filename
    eps = args.eps
    vacuum_gap = args.vacuum_gap

    # Check if the file exists
    if not os.path.exists(filename):
        print(color_text(f"Error: File '{filename}' not found.", 'red'))
        sys.exit(1)

    print(color_text(f"Reading file: {filename} (Type: fdf)", 'cyan') + "\n")

    # 3. Load the Structure object

    try:
        structure_obj = parse_fdf_to_structure(filename)
    except Exception as e:
        print(color_text(f"An error occurred while READING the file '{filename}':", 'red'))
        print(f"{e}")
        sys.exit(1)

    # 4. Run the analysis (if loading was successful)
    if structure_obj:
        get_kpath_from_structure(structure_obj, vacuum_gap=vacuum_gap, eps=eps)
    else:
        print(color_text("Error: Could not create the structure object from the file.", 'red'))

if __name__ == "__main__":
    main()
