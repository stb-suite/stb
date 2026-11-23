#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by the STB-SUITE team               #
#      stb-suite.github.io/stb                  #
#################################################

from stb import __version__ as VERSION
from stb._input_header import mark_input
from stb.core.cli import new_report_path

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


# ASE Bravais-lattice names compatible with a crystal system and centering letter of the
# space group (3D only).
_BRAVAIS_FOR_SYSTEM = {
    ("cubic", "P"): {"CUB"}, ("cubic", "F"): {"FCC"}, ("cubic", "I"): {"BCC"},
    ("tetragonal", "P"): {"TET"}, ("tetragonal", "I"): {"BCT"},
    ("orthorhombic", "P"): {"ORC"}, ("orthorhombic", "F"): {"ORCF"}, ("orthorhombic", "I"): {"ORCI"},
    ("orthorhombic", "C"): {"ORCC"}, ("orthorhombic", "A"): {"ORCC"}, ("orthorhombic", "B"): {"ORCC"},
    ("hexagonal", "P"): {"HEX"}, ("trigonal", "P"): {"HEX"}, ("trigonal", "R"): {"RHL"},
    ("monoclinic", "P"): {"MCL"}, ("monoclinic", "C"): {"MCLC"}, ("monoclinic", "A"): {"MCLC"},
    ("monoclinic", "B"): {"MCLC"}, ("monoclinic", "I"): {"MCLC"},
    ("triclinic", "P"): {"TRI"},
}
_SYSTEM_RANK = {"triclinic": 0, "monoclinic": 1, "orthorhombic": 2, "tetragonal": 3, "trigonal": 4,
                "hexagonal": 5, "cubic": 6}
_BRAVAIS_RANK = {"TRI": 0, "MCL": 1, "MCLC": 1, "ORC": 2, "ORCF": 2, "ORCI": 2, "ORCC": 2, "TET": 3,
                 "BCT": 3, "RHL": 4, "HEX": 5, "CUB": 6, "FCC": 6, "BCC": 6}
# eps values tried, in order, when the Bravais lattice found at the user's eps disagrees with
# the space group found at symprec.
_EPS_LADDER = (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05)


def expected_bravais(crystal_sys, sg_symbol):
    """The ASE Bravais-lattice names a 3D space group allows, or None when unknown."""
    if not crystal_sys or not sg_symbol:
        return None
    return _BRAVAIS_FOR_SYSTEM.get((crystal_sys.lower(), sg_symbol[0].upper()))


def bandlines_rows(path_segments):
    """[(label, n_points)] for SIESTA's %block BandLines: the first point gets 1, every point
    along a segment the line density, and the first point of a segment that does not start
    where the previous one ended (a jump, 'K|U') gets 1 -- SIESTA then starts a new line
    there instead of computing a K -> U segment that is no high-symmetry line. Before 2.1.0
    the jump point got the full count and the bands showed that segment as if it were one."""
    rows = []
    for segment in path_segments:
        if not segment:
            continue
        if not rows:
            rows.append((segment[0], 1))
        elif rows[-1][0] == segment[0]:
            pass                                   # continuous: the shared point is already there
        else:
            rows.append((segment[0], 1))           # jump: a new line starts here
        rows.extend((label, None) for label in segment[1:])
    return rows


def write_siesta_kpath_file(kpoints_dict, path_segments, num_points=50,
                             output_filename="kpath_bs.fdf", f_out=None):
    """
    Writes the k-path to a SIESTA-formatted FDF file for band structure.
    Writes ALL suggested path segments; a jump between segments starts a new
    line (1 point, see bandlines_rows). Returns True on success, False if the
    file could not be written (already reported).
    """
    if not path_segments:
        print_dual(color_text("Warning: No k-path segments found. File will not be written.", 'yellow'), f_out)
        return False

    rows = bandlines_rows(path_segments)
    if not rows:
        print_dual(color_text("[ERROR] Could not determine path sequence.", 'red'), f_out)
        return False
    for k in range(1, len(rows)):
        if rows[k][1] == 1:
            print_dual(color_text(
                f"Note: path jump {rows[k - 1][0]} | {rows[k][0]} -- written as a new line in "
                f"BandLines (1 point at {rows[k][0]}), not as a {rows[k - 1][0]}-{rows[k][0]} segment.",
                'cyan'), f_out)

    def show(label):
        return r'\Gamma' if label in ('GAMMA', r'\Gamma') else label

    display = show(rows[0][0])
    for label, n in rows[1:]:
        display += (" | " if n == 1 else "-") + show(label)
    print_dual(f"SIESTA path to be written : {display}", f_out)

    missing = [label for label, _n in rows if label not in kpoints_dict]
    if missing:
        print_dual(color_text(
            f"[ERROR] K-point label(s) {missing} found in path but not in k-points list.", 'red'), f_out)
        print_dual(color_text(f"Available labels: {list(kpoints_dict.keys())}", 'red'), f_out)
        return False
    try:
        mark_input(output_filename)
        with open(output_filename, 'w') as f:
            f.write("### BANDS\n")
            f.write(" BandLinesScale  ReciprocalLatticeVectors\n\n")
            f.write("%block BandLines\n")
            for label, n in rows:
                coord_str = " ".join([f"{c:14.10f}" for c in kpoints_dict[label]])
                f.write(f"{1 if n == 1 else num_points}   {coord_str}   {show(label)}\n")
            f.write("%endblock BandLines\n")
        print_dual(color_text(
            f"[OK] SIESTA file '{output_filename}' has been created.", 'green'), f_out)
        return True
    except Exception as e:
        print_dual(color_text(f"[ERROR] Writing {output_filename}: {e}", 'red'), f_out)
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
        sg_label = core_symmetry.space_group_label(pymatgen_structure, symprec=symprec)
        crystal_sys, point_group = core_symmetry.crystal_system(
            pymatgen_structure, symprec=symprec, angle_tolerance=angle_tolerance
        )
        eps_used = eps
        bravais_note = None
        allowed = expected_bravais(crystal_sys, sg_label.split()[0]) if dimension == 3 else None
        if allowed:
            # Only a PRIMITIVE cell's Bravais lattice must match the space group: a conventional
            # cell (cubic Si with 8 atoms, Fd-3m) legitimately has another one (CUB).
            try:
                n_prim = len(core_symmetry.SpacegroupAnalyzer(
                    pymatgen_structure, symprec=symprec, angle_tolerance=angle_tolerance).find_primitive())
                if n_prim != len(pymatgen_structure):
                    allowed = None
            except Exception:
                allowed = None
        # Act only when the lattice found is LESS symmetric than the space group needs (RHL for
        # Fm-3m): a lattice more symmetric than the crystal (cubic vectors, trigonal atoms) is legitimate.
        if allowed and _BRAVAIS_RANK.get(bravais.name, 99) >= _SYSTEM_RANK.get(crystal_sys.lower(), -1):
            allowed = None
        if allowed and bravais.name not in allowed:
            # A relaxed cubic cell with angles of 59.91-59.99 deg is RHL at eps 2e-4 while spglib
            # (symprec 0.01) still finds Fm-3m: the path then follows the rhombohedral lattice.
            # Use the smallest eps that agrees with the space group; the k-points stay
            # fractional coordinates of the user's own cell.
            found = None
            for trial in _EPS_LADDER:
                if trial <= eps:
                    continue
                try:
                    candidate = cell.get_bravais_lattice(eps=trial, pbc=pbc)
                except Exception:
                    continue
                if candidate.name in allowed:
                    found = (trial, candidate)
                    break
            angles = cell.cellpar()[3:]
            if found:
                bravais_note = (
                    f"The cell's Bravais lattice at eps = {eps} is {bravais.name}, but the space group "
                    f"{sg_label} (symprec {symprec}) is {crystal_sys} ({'/'.join(sorted(allowed))}): the "
                    f"cell is slightly distorted (angles {angles[0]:.4f}, {angles[1]:.4f}, {angles[2]:.4f} "
                    f"deg). The path uses eps = {found[0]}, which gives {found[1].name}. Symmetrize the "
                    "cell (or pass -p) to remove the distortion itself.")
                eps_used, bravais = found
            else:
                bravais_note = (
                    f"The cell's Bravais lattice ({bravais.name}) disagrees with the space group "
                    f"{sg_label} (symprec {symprec}, {crystal_sys}) for every eps up to "
                    f"{_EPS_LADDER[-1]}: the path follows {bravais.name}. Check the cell, or "
                    "symmetrize it first.")

        dim_label = {1: "1D", 2: "2D", 3: "3D"}[dimension]
        print_dual(f"Dimensionality : {dim_label}", f_out)
        print_dual(f"Bravais lattice : {bravais.longname} ({bravais.name})", f_out)
        if dimension < 3:
            print_dual(color_text(
                "Note: the Bravais lattice above reflects only the periodic axes "
                "(vacuum-padded axes excluded); it is not a 3D space group.", 'yellow'), f_out)
        if bravais_note:
            print_dual(color_text(f"[WARNING] {bravais_note}", 'yellow'), f_out)

        print_dual(f"Space group : {sg_label}", f_out)
        print_dual(f"Crystal system : {crystal_sys} (point group {point_group})", f_out)
        if dimension < 3:
            print_dual(color_text(
                "Note: this space group treats the vacuum-padded axis/axes as an "
                "ordinary periodic direction and may not reflect the true symmetry. "
                "Use stb-symmetry (code 2.5) for a dimension-aware layer-group/"
                "point-group analysis.", 'yellow'), f_out)

        bp = cell.bandpath(pbc=pbc, npoints=0, eps=eps_used)

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
        citations.cite(bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{citations.run_references_path()}' "
            f" (with the libraries the run loads; written when it ends).", 'green'), f_out)

        print_section("[6] SUMMARY & FILES", f_out)
        print_dual("Status         : OK", f_out)
        print_dual(f"Path           : {path_str}", f_out)
        print_dual(f"Output file    : {output_filename}", f_out)
        print_dual(f"References     : {citations.run_references_path()}", f_out)
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
        default=0.01,
        help="Symmetry-detection tolerance for the crystallographic space-group "
             "analysis. Default: 0.01 (pymatgen's own default). Loosen further "
             "(e.g. 0.02-0.05) for a structure relaxed with a looser force "
             "tolerance -- a tighter value than the relaxation's own residual "
             "numerical noise can misdetect the true space group."
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
                        help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-kpath {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | 2026",
            "Developed by the STB-SUITE team"
        ])

    report_path = new_report_path("kpath")
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
