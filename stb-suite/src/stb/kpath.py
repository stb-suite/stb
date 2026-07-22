#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import numpy as np
from ase.cell import Cell
from ase.dft.kpoints import parse_path_string
from stb.core import structure_io, kspace
from stb.core import symmetry as core_symmetry
from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section
import sys
import os
import argparse

REPORT_FILE = "stb_kpath_report.txt"
BIB_FILE = "references.bib"

# The actual band-path convention ASE's Cell.bandpath implements (extended to
# 1D/2D via its periodic-axes mask) -- kept local since this is currently the
# only consumer (extract-on-second-use, same policy as core/citations.py).
_BIB_SETYAWAN_CURTAROLO = ("SetyawanCurtarolo2010", """@article{SetyawanCurtarolo2010,
  author  = {Setyawan, Wahyu and Curtarolo, Stefano},
  title   = {High-throughput electronic band structure calculations: Challenges and tools},
  journal = {Computational Materials Science},
  year    = {2010},
  volume  = {49},
  number  = {2},
  pages   = {299--312},
  doi     = {10.1016/j.commatsci.2010.05.010}
}""")


def write_siesta_kpath_file(kpoints_dict, path_segments, num_points=50,
                             output_filename="kpath_bs.fdf", f_out=None):
    """
    Writes the k-path to a SIESTA-formatted FDF file for band structure.
    Writes ALL suggested path segments, even if disjointed. Returns True on
    success, False if the file could not be written (already reported).
    """
    if not path_segments:
        print_dual(color_text("Warning: No k-path segments found. File will not be written.", 'yellow'), f_out)
        return False

    # Rebuilds the 'path_sequence' based on the correct data structure
    # (which is a list of paths, e.g: [['A', 'B', 'C'], ['D', 'E']])
    path_sequence = []
    for i, segment_list in enumerate(path_segments):  # ex: segment_list = ['\Gamma', 'Y', 'H']
        if not segment_list:
            continue

        if i == 0:
            path_sequence.extend(segment_list)
        else:
            last_point_in_sequence = path_sequence[-1]
            first_point_of_new_segment = segment_list[0]

            if last_point_in_sequence == first_point_of_new_segment:
                # Continuous path (e.g., ends in Z, starts in Z) -- skip the duplicate.
                path_sequence.extend(segment_list[1:])
            else:
                # Disjoint path (a "jump", e.g., ends in H_1, starts in M).
                print_dual(color_text(
                    f"Note: Disjointed path detected (jump from {last_point_in_sequence} "
                    f"to {first_point_of_new_segment}). Including full segment.", 'cyan'), f_out)
                path_sequence.extend(segment_list)

    if not path_sequence:
        print_dual(color_text("[ERROR] Could not determine path sequence.", 'red'), f_out)
        return False

    path_str_display = '-'.join([r'\Gamma' if p == 'GAMMA' else p for p in path_sequence])
    print_dual(f"SIESTA path to be written : {path_str_display}", f_out)

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

        print_dual(color_text(
            f"[OK] SIESTA file '{output_filename}' has been created.", 'green'), f_out)
        return True

    except KeyError as e:
        # A label in the path (e.g. '\Gamma') has a different name in the
        # dictionary (e.g. 'GAMMA') -- retry once with that substitution.
        if str(e) == r"'\Gamma'" and 'GAMMA' in kpoints_dict:
            print_dual(color_text(
                "Warning: Found '\\Gamma' label, attempting to use 'GAMMA' internally.", 'yellow'), f_out)
            path_sequence_fixed = ['GAMMA' if label == r'\Gamma' else label for label in path_sequence]
            return _write_siesta_kpath_file_fixed(
                kpoints_dict, path_sequence_fixed, num_points, output_filename, f_out
            )
        print_dual(color_text(
            f"[ERROR] K-point label {e} found in path but not in k-points list.", 'red'), f_out)
        print_dual(color_text(f"Available labels: {list(kpoints_dict.keys())}", 'red'), f_out)
        return False

    except Exception as e:
        print_dual(color_text(f"[ERROR] Writing {output_filename}: {e}", 'red'), f_out)
        return False


def _write_siesta_kpath_file_fixed(kpoints_dict, path_sequence, num_points, output_filename, f_out=None):
    """Retries the write after substituting '\\Gamma' -> 'GAMMA' in the path sequence."""
    print_dual(color_text("Retrying file write with 'GAMMA' label...", 'cyan'), f_out)
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

        print_dual(color_text(
            f"[OK] SIESTA file '{output_filename}' has been created (with label fix).", 'green'), f_out)
        return True
    except Exception as e:
        print_dual(color_text(f"[ERROR] During retry: {e}", 'red'), f_out)
        return False


def run_kpath(filename, vacuum_gap, eps, symprec, angle_tolerance, output_filename, f_out=None):
    """
    Core logic: reads `filename`, analyzes its dimensionality and symmetry,
    computes its high-symmetry k-path (dimension-aware, via ASE's
    Cell.bandpath), and writes it to `output_filename`. Returns True on
    success, False on any handled failure (missing file, parse error, a 0D
    structure -- a k-path needs at least one periodic direction -- or an
    exception during detection).
    """
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file : {filename}", f_out)
    print_dual(f"Vacuum-gap threshold : {vacuum_gap} Ang", f_out)
    print_dual(f"Bravais-lattice tolerance (eps) : {eps}", f_out)
    print_dual(f"Symmetry tolerance : symprec={symprec}, angle={angle_tolerance} deg", f_out)

    if not os.path.exists(filename):
        print_dual(color_text(f"[ERROR] Structure file '{filename}' not found.", 'red'), f_out)
        return False

    try:
        fdf_structure = structure_io.read_fdf(filename)
    except Exception as e:
        print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
        return False

    try:
        print_section("[1] STRUCTURE ANALYSIS", f_out)
        lattice = fdf_structure.lattice
        positions = np.array([pos for _, pos in fdf_structure.atoms])
        is_cartesian = fdf_structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice, vacuum_gap)
        dimension = 3 - sum(vacuum_axes)

        pymatgen_structure = structure_io.to_pymatgen(fdf_structure)
        print_dual(f"Formula : {pymatgen_structure.composition.reduced_formula}", f_out)

        if dimension == 0:
            print_dual(color_text("Dimensionality : 0D (isolated molecule)", 'bold'), f_out)
            print_dual(color_text(
                "[ERROR] No periodic direction detected (every axis is vacuum-padded) -- "
                "a k-path is not physically meaningful for an isolated molecule (0D). "
                "No output file written.", 'red'), f_out)
            return False

        pbc = tuple(not v for v in vacuum_axes)
        cell = Cell(lattice)
        bravais = cell.get_bravais_lattice(eps=eps, pbc=pbc)

        dim_label = {1: "1D", 2: "2D", 3: "3D"}[dimension]
        print_dual(f"Dimensionality : {dim_label}", f_out)
        print_dual(f"Bravais lattice : {bravais.longname} ({bravais.name})", f_out)
        if dimension < 3:
            print_dual(color_text(
                "Note: the Bravais lattice above reflects only the periodic axes "
                "(vacuum-padded axes excluded); it is not a 3D space group.", 'yellow'), f_out)

        sg_label = core_symmetry.space_group_label(pymatgen_structure, symprec=symprec)
        crystal_sys, point_group = core_symmetry.crystal_system(
            pymatgen_structure, symprec=symprec, angle_tolerance=angle_tolerance
        )
        print_dual(f"Space group : {sg_label}", f_out)
        print_dual(f"Crystal system : {crystal_sys} (point group {point_group})", f_out)
        if dimension < 3:
            print_dual(color_text(
                "Note: this space group treats the vacuum-padded axis/axes as an "
                "ordinary periodic direction and may not reflect the true symmetry. "
                "Use stb-symmetry (code 3.5) for a dimension-aware layer-group/"
                "point-group analysis.", 'yellow'), f_out)

        bp = cell.bandpath(pbc=pbc, npoints=0, eps=eps)

        # ASE labels the Gamma point 'G'; translate to 'GAMMA' so the writer's
        # existing '\Gamma' display/formatting logic keeps working unchanged.
        kpoints_dict = {('GAMMA' if label == 'G' else label): coords
                        for label, coords in bp.special_points.items()}
        path_segments = [[('GAMMA' if label == 'G' else label) for label in segment]
                         for segment in parse_path_string(bp.path)]

        print_section("[2] HIGH-SYMMETRY K-POINTS", f_out)
        for label, coords in kpoints_dict.items():
            coord_str = ", ".join([f"{c:8.5f}" for c in coords])
            display_label = r'\Gamma' if label == 'GAMMA' else label
            print_dual(f"  {display_label:<5}: ({coord_str})", f_out)

        print_section("[3] SUGGESTED K-PATH", f_out)
        path_segments_str_list = []
        for segment_list in path_segments:
            display_segment = [r'\Gamma' if label == 'GAMMA' else label for label in segment_list]
            path_segments_str_list.append("-".join(display_segment))
        path_str = ' | '.join(path_segments_str_list)
        print_dual(f"Path : {path_str}", f_out)
        print_dual(
            "Methodology : ASE's dimension-aware Bravais-lattice convention, extending "
            "the Setyawan & Curtarolo (2010) 3D scheme to 1D/2D systems via their "
            "periodic-axes mask (pbc).", f_out)

        print_section("[4] WRITING OUTPUT FILE", f_out)
        if not write_siesta_kpath_file(kpoints_dict, path_segments, 50, output_filename, f_out):
            return False

        print_section("[5] REFERENCES", f_out)
        bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, _BIB_SETYAWAN_CURTAROLO]
        citations.write_bib_file(BIB_FILE, bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)

        print_section("[6] SUMMARY & FILES", f_out)
        print_dual("Status         : OK", f_out)
        print_dual(f"Path           : {path_str}", f_out)
        print_dual(f"Output file    : {output_filename}", f_out)
        print_dual(f"References     : {BIB_FILE}", f_out)
        return True

    except Exception as e:
        print_dual(color_text(f"[ERROR] An error occurred while analyzing the structure: {e}", 'red'), f_out)
        print_dual(color_text(
            f"It might be necessary to adjust 'vacuum_gap' or 'eps' (current: "
            f"{vacuum_gap} Ang, {eps}). This often happens if the input structure "
            "is distorted or the parse failed.", 'yellow'), f_out)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Finds the high-symmetry k-path for a SIESTA structure file (.fdf). "
                    "Dimension-aware (1D/2D/3D) via ASE's Bravais-lattice path finder.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-f", "--file",
        dest="filename",
        type=str,
        required=True,
        help="The structure file name (e.g., struct.fdf)."
    )
    parser.add_argument(
        "-p", "--prec",
        dest="eps",
        type=float,
        default=0.0002,
        help="Tolerance (ase.cell.Cell's 'eps') for Bravais-lattice/symmetry detection.\n"
             "Default: 0.0002"
    )
    parser.add_argument(
        "--vacuum-gap",
        type=float,
        default=10.0,
        help="Minimum empty span (in Angstrom) between atoms along an axis, wrapped "
             "periodically, to treat that axis as vacuum-padded (non-periodic) when "
             "detecting the system's dimensionality. Default: 10.0"
    )
    parser.add_argument(
        "--symprec",
        type=float,
        default=1e-3,
        help="Symmetry-detection tolerance for the crystallographic space-group "
             "analysis. Default: 0.001"
    )
    parser.add_argument(
        "--angle-tolerance",
        type=float,
        default=5.0,
        help="Symmetry angle tolerance (degrees) for the crystallographic space-group "
             "analysis. Default: 5.0"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="kpath_bs.fdf",
        help="Output .fdf file name (default: kpath_bs.fdf)."
    )
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kpath {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-KPATH REPORT =====", 'magenta'), f_out)

    ok = run_kpath(
        args.filename, args.vacuum_gap, args.eps, args.symprec, args.angle_tolerance,
        args.output, f_out
    )

    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
