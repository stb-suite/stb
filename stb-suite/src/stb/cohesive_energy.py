#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"  # stb-standard numbered report ([0]-[7], --save-report) and a
                    # single wrapping --output-dir (default cohesive_runs, matching
                    # stb-strain/stb-elasticInputs' own <property>_runs convention)
                    # instead of writing structure/atoms/atoms_bsse[_check] loose
                    # into the current directory by default.

import os
import sys
import re
import shutil
import argparse
import numpy as np
from stb.core import structure_io, kspace, symmetry
from stb.core.cli import color_text, print_dual, print_section, print_table, show_intro
from stb.core.pseudopotentials import BANKS, get_required_pseudos, resolve_pseudo_source, link_pseudo

REPORT_FILE = "cohesive_stage1.txt"

# Base Template
CALC_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
SystemLabel      siesta
SystemName       siesta

%include structure.fdf

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry
FilterCutoff          150  Ry

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
XC.Functional       GGA
XC.Authors          PBE

## ===================================================================
## VAN DER WAALS CORRECTION (GRIMME DFT-D3) -- same convention as
## stb-inputfile's own CALC_*_TEMPLATE blocks
## ===================================================================
DFTD3                   .false.

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
Spin                non-polarized

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## ===================================================================
## STARTING THE CALCULATION
## ===================================================================
DM.UseSaveDM            .true.
"""

def parse_structure_fdf(filename, vacuum_gap):
    """Parses .fdf and returns (FdfStructure, per-axis vacuum flags -- see
    kspace.detect_vacuum_axes), used to compute a dimensionality-aware k-grid
    and (for --bsse-correction) the real bulk geometry around each species.
    """
    try:
        structure = structure_io.read_fdf(filename)
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
        return structure, vacuum_axes
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{filename}' not found.", 'red'))
        sys.exit(1)
    except Exception as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

def generate_isolated_atom_fdf(symbol, z_num, out_path, vacuum):
    """Creates a structure.fdf for a single isolated atom in a large box.

    `vacuum` is the cubic box side, in Ang -- the atom sits at the
    fractional center (0.5, 0.5, 0.5), so it stays centered regardless of
    the box size substituted here.
    """
    content = f"""# Isolated {symbol} atom for cohesive energy
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   {z_num}   {symbol}
%endblock ChemicalSpeciesLabel

LatticeConstant 1.00 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 {vacuum:.6f}   0.000000   0.000000
  0.000000  {vacuum:.6f}   0.000000
  0.000000   0.000000  {vacuum:.6f}
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.500000000   0.500000000   0.500000000   1
%endblock AtomicCoordinatesAndAtomicSpecies
"""
    with open(out_path, 'w') as f:
        f.write(content)
    return

def find_ghost_neighbors(pmg_structure, anchor_index, cutoff):
    """Returns a list of (element_symbol, atomic_number, relative_cartesian_vector)
    for every periodic neighbor of `pmg_structure[anchor_index]` within
    `cutoff` Ang -- the real local coordination environment used to build a
    BSSE (counterpoise) ghost cluster around that site. Works unmodified for
    any dimensionality (0D/1D/2D/3D): a vacuum-padded axis simply has no
    periodic image within any physically sensible cutoff (a handful of Ang),
    since the vacuum gap is by construction much larger than that.
    """
    anchor = pmg_structure[anchor_index]
    neighbors = pmg_structure.get_neighbors(anchor, r=cutoff)
    return [(n.specie.symbol, n.specie.Z, n.coords - anchor.coords) for n in neighbors]


def build_ghost_cluster(anchor_symbol, anchor_Z, neighbors, vacuum_box):
    """Builds a synthetic FdfStructure for a BSSE (counterpoise) cluster: the
    real anchor atom at the box center, plus one ghost entry per unique
    neighboring element (SIESTA convention: Z negative, no valence charge,
    same basis as the real element -- see link_pseudo's dest_label), placed
    at the anchor's REAL local coordination geometry (from find_ghost_neighbors).
    `vacuum_box` is the cubic box side in Ang (anchor at the fractional
    center, same convention as generate_isolated_atom_fdf).
    """
    species = [anchor_symbol]
    species_meta = {anchor_symbol: {'id': '1', 'Z': anchor_Z}}
    atoms = [(anchor_symbol, np.array([0.5, 0.5, 0.5]))]

    next_id = 2
    for elem, z, rel_vec in neighbors:
        ghost_label = f"{elem}_ghost"
        if ghost_label not in species_meta:
            species.append(ghost_label)
            species_meta[ghost_label] = {'id': str(next_id), 'Z': -abs(z)}
            next_id += 1
        frac_pos = np.array([0.5, 0.5, 0.5]) + rel_vec / vacuum_box
        atoms.append((ghost_label, frac_pos))

    return structure_io.FdfStructure(
        lattice=np.eye(3) * vacuum_box,
        lattice_constant=1.0,
        species=species,
        species_meta=species_meta,
        atoms=atoms,
        coord_format="fractional",
        raw_lines=[],
    )


def _isolated_atom_calc(dispersion):
    """calc.fdf text shared by every isolated-atom-type calculation (plain
    isolated atom, BSSE ghost cluster, BSSE convergence-check cluster):
    always spin-polarized, Gamma-only k-grid, optionally with D3 dispersion.
    """
    calc = CALC_TEMPLATE.replace("Spin                non-polarized", "Spin                polarized")
    calc = re.sub(r'kgrid\.MonkhorstPack\s+\[.*?\]', 'kgrid.MonkhorstPack   [1  1  1]', calc)
    if dispersion:
        calc = calc.replace("DFTD3                   .false.", "DFTD3                   .true.")
    return calc


def generate_bsse_reference(pmg_structure, sym, anchor_index, z_num, cutoff, vacuum, out_dir,
                             pp_path, dispersion):
    """Builds and writes ONE BSSE (counterpoise) ghost-cluster reference
    (structure.fdf + calc.fdf + linked pseudopotentials, real anchor + real
    ghost neighbors) at `out_dir`, anchored at `pmg_structure[anchor_index]`.
    Shared by the main --bsse-cutoff reference and the optional --bsse-
    convergence-check reference (same construction, different cutoff/dir),
    and by every site of a multi-site species (see find_ghost_neighbors/
    build_ghost_cluster above -- this just wires them together with the
    calc.fdf template and pseudopotential linking).

    Returns the number of ghost neighbors found (0 if none -- the caller
    decides whether that's worth a warning, e.g. cutoff too small).

    Also writes a small 'bsse_cutoff.txt' sidecar with the numeric `cutoff`
    -- a flat 'atoms_bsse/'/'atoms_bsse_check/'-style folder name carries no
    cutoff information of its own (unlike the multi-scan 'atoms_bsse_check_
    <cutoff>/' naming), so stb-cohesiveAnalysis has no other way to recover
    it later for reporting/plotting (see core.cohesive_shared eventually, or
    cohesive_analysis.py::find_reference_cutoff for now).
    """
    neighbors = find_ghost_neighbors(pmg_structure, anchor_index, cutoff)
    cluster = build_ghost_cluster(sym, z_num, neighbors, vacuum)
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(cluster, os.path.join(out_dir, "structure.fdf"))

    with open(os.path.join(out_dir, "calc.fdf"), 'w') as f:
        f.write(_isolated_atom_calc(dispersion))
    with open(os.path.join(out_dir, "bsse_cutoff.txt"), 'w') as f:
        f.write(f"{cutoff}\n")

    link_pseudo(pp_path, sym, out_dir)
    for elem, _z, _rel in neighbors:
        link_pseudo(pp_path, elem, out_dir, dest_label=f"{elem}_ghost")
    return len(neighbors)


def resolve_bsse_sites(pmg_structure, sym, symprec, multi_site, fallback_index):
    """Decides which site(s) of species `sym` get their own BSSE ghost-
    cluster reference, and returns (sites_to_use, all_sites_found,
    space_group) where the first two are lists of (site_index,
    wyckoff_letter_or_None, multiplicity) and `space_group` is the
    "<international symbol> (No. <number>)" label (None if detection
    failed).

    Uses core/symmetry.py::find_inequivalent_sites (already shared with
    stb-defect/stb-hubbardu) to detect symmetrically distinct sites for this
    species. If detection fails or finds nothing (advisory only, never
    blocks BSSE generation), falls back to the single `fallback_index`
    (first occurrence) with unknown Wyckoff/multiplicity -- this tool's
    original, pre-multi-site behavior.

    When 2+ distinct sites are found: `multi_site=True` (default, see
    --bsse-multi-site) returns ALL of them, so the caller builds one ghost
    cluster per site (a real physical distinction -- e.g. octahedral vs.
    tetrahedral Fe -- would otherwise be averaged away into a single,
    possibly unrepresentative reference); `multi_site=False` returns only
    the first (by site index), preserving the flat single-reference layout
    while still telling the caller how many distinct sites exist so it can
    warn about the ones being skipped.
    """
    try:
        sites, space_group = symmetry.find_inequivalent_sites(pmg_structure, symprec, filter_species=sym)
    except Exception:
        sites = []
        space_group = None
    if not sites:
        sites = [(fallback_index, None, None)]

    if len(sites) > 1 and not multi_site:
        return [sites[0]], sites, space_group
    return sites, sites, space_group


def main():
    parser = argparse.ArgumentParser(
        description="Prepare folder structure for cohesive energy calculations.",
        epilog="Example usage:\n"
               "  stb_cohesive -s structure.fdf -k 0.15\n"
               "  stb_cohesive -s structure.fdf --spin -p /path/to/pseudos",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-s", "--structure", dest="structure", type=str, required=True,
                        help="Path to the initial structure.fdf file")

    parser.add_argument("-p", "--pp-path", dest="pp_path", type=str, default="",
                        help=f"Pseudopotentials source (optional): a bundled bank "
                             f"({', '.join(BANKS)}) or a folder path.")

    parser.add_argument("-k", "--k-density", dest="k_density", type=float, default=0.2,
                        help="K-point density in 1/Angstrom (default: 0.2)")

    parser.add_argument("--spin", dest="spin", action="store_true",
                        help="Set the full structure calculation to spin polarized")

    parser.add_argument("--vacuum", dest="vacuum", type=float, default=20.0,
                        help="Guaranteed buffer (Ang/2, see below) between the outermost atom "
                             "of each isolated-atom-type calculation and the edge of its cubic "
                             "simulation box (default: 20.0). For a plain isolated atom (no "
                             "neighbors) the box is exactly this value, same as before. For a "
                             "BSSE ghost cluster (see --bsse-cutoff) or its --bsse-convergence-"
                             "check counterpart, the box automatically GROWS with the cluster's "
                             "own radius (box side = 2*cutoff + --vacuum) so the buffer beyond "
                             "the outermost ghost atom stays exactly --vacuum/2 regardless of "
                             "how large --bsse-cutoff/--bsse-convergence-increment are -- "
                             "otherwise a bigger convergence-check cluster would eat into a "
                             "fixed-size box and risk interacting with its own periodic images, "
                             "contaminating the very convergence check it's meant to validate. "
                             "Increase for elements with unusually diffuse basis orbitals.")

    parser.add_argument("--vacuum-gap", dest="vacuum_gap", type=float, default=10.0,
                        help="Vacuum gap threshold in Ang used to detect which lattice axes "
                             "of the FULL structure are periodic vs. vacuum-padded (default: "
                             "10.0), same convention as stb-kgrid/stb-strain/stb-elasticInputs. "
                             "Vacuum-padded axes are forced to a single (Gamma-only) k-grid "
                             "division regardless of --k-density.")

    parser.add_argument("-O", "--output-dir", dest="output_dir", default="cohesive_runs",
                        help="Root directory (default: cohesive_runs, matching stb-strain's/"
                             "stb-elasticInputs' own <property>_runs convention) for the "
                             "'structure', 'atoms' and 'atoms_bsse' folders -- so "
                             "stb-cohesiveAnalysis can just 'cd <output-dir>' (or pass --dir) "
                             "and find everything with no extra flags. Set this to avoid "
                             "silently overwriting a previous run's setup when preparing more "
                             "than one structure from the same working directory.")

    parser.add_argument("--bsse-correction", dest="bsse_correction", action="store_true",
                        default=True,
                        help="Also generate a BSSE (Basis Set Superposition Error) -corrected "
                             "reference per species, in 'atoms_bsse/<symbol>/' (default: ON). "
                             "LCAO cohesive energies are otherwise systematically over-bound: "
                             "the isolated atom is computed with a smaller effective basis than "
                             "it has inside the solid (no neighboring atoms' orbitals "
                             "available). The correction embeds the same real atom in 'ghost' "
                             "copies (SIESTA convention: negative Z, no valence charge, same "
                             "basis as the real element) of its actual nearest neighbors from "
                             "the full structure (see --bsse-cutoff), giving it the same "
                             "effective basis it has in the solid without changing its physics. "
                             "Doubles the number of isolated-atom-type DFT calculations needed "
                             "(one plain + one BSSE-corrected per species).")
    parser.add_argument("--no-bsse-correction", dest="bsse_correction", action="store_false",
                        help="Skip the BSSE-corrected reference (halves the isolated-atom-type "
                             "calculation count, at the cost of a systematically over-bound "
                             "cohesive energy).")
    parser.add_argument("--bsse-cutoff", dest="bsse_cutoff", type=float, default=4.0,
                        help="Radius in Ang (default: 4.0) around each species' anchor atom "
                             "within which real neighbors become ghost atoms in its BSSE-"
                             "corrected reference (see --bsse-correction). How many "
                             "coordination shells this actually reaches depends on the "
                             "material's bond length -- for short-bonded covalent solids (e.g. "
                             "graphene, diamond, ~1.4-1.5 Ang bonds) it typically reaches "
                             "several shells; for longer-bonded ionic/metallic solids it may "
                             "only reach the 1st. Check the printed ghost-neighbor count per "
                             "species; the correction is converged once increasing --bsse-"
                             "cutoff no longer changes the reported BSSE-corrected energy "
                             "appreciably. Only used when --bsse-correction is on.")
    parser.add_argument("--symprec", dest="symprec", type=float, default=1e-3,
                        help="Symmetry-detection tolerance (default: 1e-3, pymatgen's own "
                             "default), used to find symmetrically distinct sites per species "
                             "for --bsse-multi-site.")
    parser.add_argument("--bsse-multi-site", dest="bsse_multi_site", action="store_true",
                        default=True,
                        help="When a species occupies 2+ symmetrically distinct sites (e.g. "
                             "octahedral vs. tetrahedral Fe), build a separate BSSE ghost-"
                             "cluster reference PER site (default: ON) -- "
                             "stb-cohesiveAnalysis then weights each site's BSSE-corrected "
                             "energy by its multiplicity, instead of using a single "
                             "(potentially unrepresentative) site for the whole species.")
    parser.add_argument("--no-bsse-multi-site", dest="bsse_multi_site", action="store_false",
                        help="Use only the first detected site as the BSSE reference for a "
                             "species with multiple distinct sites (prints a warning naming "
                             "the sites being skipped), instead of one cluster per site.")
    parser.add_argument("--bsse-convergence-check", dest="bsse_convergence_check",
                        action="store_true", default=False,
                        help="Also generate one or more extra BSSE references per species/site "
                             "at a larger cutoff (see --bsse-convergence-increment) -- "
                             "stb-cohesiveAnalysis then reports how much the BSSE-corrected "
                             "cohesive energy shifts as the cutoff grows, as a convergence "
                             "check (default: OFF -- on top of --bsse-correction's own "
                             "doubling, each extra cutoff adds one more isolated-atom-type "
                             "calculation per species/site). Requires --bsse-correction.")
    parser.add_argument("--bsse-convergence-increment", dest="bsse_convergence_increment",
                        type=float, default=[2.0], nargs='+',
                        help="One or more Ang values added to --bsse-cutoff (default: 2.0, a "
                             "single point). Pass several (e.g. --bsse-convergence-increment 2 "
                             "4 6) to scan multiple cutoffs in one run and see the correction's "
                             "full trend instead of a single before/after comparison -- "
                             "stb-cohesiveAnalysis will report and plot the whole curve. With "
                             "exactly one increment (the default), the extra reference is "
                             "written to 'atoms_bsse_check/' as before; with 2+, each is "
                             "written to its own 'atoms_bsse_check_<cutoff>/' instead.")

    parser.add_argument("--dispersion", dest="dispersion", action="store_true",
                        help="Enable Grimme DFT-D3 dispersion correction (SIESTA's own "
                             "DFTD3 directive, default: OFF) for every calculation (full "
                             "structure, isolated atoms, and BSSE references alike) -- relevant "
                             "for layered/molecular materials where van der Waals interactions "
                             "contribute non-negligibly to the cohesive energy. Safe to apply "
                             "uniformly: a genuinely isolated atom has no neighbor to form a "
                             "dispersion pair with, so D3 contributes ~0 there regardless.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by "
                             "default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-cohesive {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    if args.vacuum <= 0:
        parser.error("--vacuum must be positive.")
    if args.bsse_correction and args.bsse_cutoff <= 0:
        parser.error("--bsse-cutoff must be positive.")
    if args.bsse_convergence_check:
        if not args.bsse_correction:
            parser.error("--bsse-convergence-check requires --bsse-correction.")
        if any(inc <= 0 for inc in args.bsse_convergence_increment):
            parser.error("--bsse-convergence-increment must be positive.")

    # Computed up front (used by both [0] RUN METADATA and [5] BSSE
    # CORRECTION below): the sorted, de-duplicated set of extra cutoffs to
    # generate. Exactly one -> the classic single before/after check
    # ('atoms_bsse_check/', unchanged naming/behavior); two or more -> a full
    # convergence scan ('atoms_bsse_check_<cutoff>/' per point).
    check_cutoffs = (
        sorted({round(args.bsse_cutoff + inc, 6) for inc in args.bsse_convergence_increment})
        if args.bsse_convergence_check else []
    )
    multi_check = len(check_cutoffs) > 1

    if args.pp_path:
        try:
            args.pp_path = resolve_pseudo_source(args.pp_path)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "Cohesive Energy Workflow setup",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-COHESIVE STAGE 1 REPORT (COHESIVE ENERGY PREP) =====", 'magenta'), f_out)

    # --- [0] RUN METADATA ---
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file        : {args.structure}", f_out)
    print_dual(f"Pseudo source         : {args.pp_path or '(not given)'}", f_out)
    print_dual(f"K-point density       : {args.k_density} 1/Ang", f_out)
    print_dual(f"Spin (full structure) : {'polarized' if args.spin else 'non-polarized'}", f_out)
    print_dual(f"Dispersion (D3)       : {'ON' if args.dispersion else 'OFF'}", f_out)
    print_dual(f"Isolated-atom vacuum  : {args.vacuum} Ang (buffer beyond the outermost atom; "
               "BSSE cluster boxes auto-grow with --bsse-cutoff to keep this buffer constant)",
               f_out)
    print_dual(f"Vacuum-gap threshold  : {args.vacuum_gap} Ang", f_out)
    if args.bsse_correction:
        print_dual(f"BSSE correction       : ON (cutoff={args.bsse_cutoff} Ang, "
                   f"multi-site={'on' if args.bsse_multi_site else 'off'})", f_out)
        if args.bsse_convergence_check:
            if multi_check:
                cuts = ", ".join(f"{c:.1f}" for c in check_cutoffs)
                print_dual(f"BSSE convergence check: ON, scanning {len(check_cutoffs)} cutoffs "
                           f"({cuts} Ang)", f_out)
            else:
                inc = args.bsse_convergence_increment[0]
                print_dual(f"BSSE convergence check: ON (+{inc:g} Ang -> cutoff "
                           f"{check_cutoffs[0]:.1f} Ang)", f_out)
        else:
            print_dual("BSSE convergence check: OFF", f_out)
        print_dual(f"Symmetry tolerance    : symprec={args.symprec}", f_out)
    else:
        print_dual("BSSE correction       : OFF", f_out)
    print_dual(f"Output directory      : {args.output_dir}", f_out)
    print_dual(f"Save report           : {'yes' if args.save_report else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file           : {report_path}", f_out)

    # --- [1] INPUT STRUCTURE ---
    print_section("[1] INPUT STRUCTURE", f_out)
    structure, vacuum_axes = parse_structure_fdf(args.structure, args.vacuum_gap)
    lattice = structure.lattice
    species_meta = structure_io.species_dict(structure)
    counts = structure_io.atom_counts(structure)
    # Species declared in ChemicalSpeciesLabel but never placed in the
    # coordinates block need no isolated-atom (or BSSE) calculation at all --
    # they'd contribute exactly 0 to cohesive_analysis.py's sum regardless of
    # their own energy, so requiring one would only waste DFT time and (if
    # missing) incorrectly block the analysis (see core/structure_io.py::
    # atom_counts docstring).
    species = {sym: meta for sym, meta in species_meta.items() if counts.get(sym, 0) > 0}
    unused_species = [sym for sym in species_meta if counts.get(sym, 0) == 0]
    composition = ", ".join(f"{sym}{n}" for sym, n in counts.items())
    cell_volume = abs(np.linalg.det(lattice))
    print_table(["Quantity", "Value"], [
        (["Composition", composition], None),
        (["Total atoms", str(sum(counts.values()))], None),
        (["Species (need a calculation)", ", ".join(species.keys())], None),
        (["Coordinate format", structure.coord_format], None),
        (["Lattice constant", f"{structure.lattice_constant} Ang"], None),
        (["Cell volume", f"{cell_volume:.4f} Ang^3"], None),
    ], f_out)
    if unused_species:
        print_dual(color_text(
            f"[INFO] Declared but never placed (no isolated-atom calculation needed): "
            f"{', '.join(unused_species)}", 'cyan'), f_out)

    # --- [2] DIMENSIONALITY & K-GRID ---
    print_section("[2] DIMENSIONALITY & K-GRID", f_out)
    print_dual(f"Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)
    print_table(["Axis", "Status"], [
        ([letter, "vacuum" if is_vac else "periodic"], 'yellow' if is_vac else None)
        for letter, is_vac in zip('xyz', vacuum_axes)
    ], f_out)
    try:
        kgrid_divs = kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], args.k_density, vacuum_axes)
    except ValueError as e:
        print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"K-grid (full structure): {kgrid_divs[0]} {kgrid_divs[1]} {kgrid_divs[2]} "
               f"(density={args.k_density} 1/Ang)", f_out)
    print_dual("K-grid (isolated atoms/BSSE clusters): 1 1 1 (Gamma-only, always -- an "
               "isolated cluster in a large vacuum box has no dispersion to sample)", f_out)

    # --- [3] FULL STRUCTURE SETUP ---
    print_section("[3] FULL STRUCTURE SETUP", f_out)
    struct_dir = os.path.join(args.output_dir, "structure")
    os.makedirs(struct_dir, exist_ok=True)
    shutil.copy(args.structure, os.path.join(struct_dir, "structure.fdf"))

    struct_calc = CALC_TEMPLATE
    if args.spin:
        struct_calc = struct_calc.replace("Spin                non-polarized", "Spin                polarized")
    else:
        print_dual(color_text(
            "[NOTE] Spin polarization is OFF for the full structure. Every isolated-atom (and "
            "BSSE-corrected) reference is always spin-polarized regardless of this setting -- "
            "correct for closed-shell/non-magnetic bulk materials, but if this structure is "
            "expected to order magnetically, its non-polarized bulk energy will be compared "
            "against polarized atomic references, an inconsistent (not just incomplete) "
            "comparison. Pass --spin if in doubt.", 'yellow'), f_out)

    kgrid_str = f"kgrid.MonkhorstPack   [{kgrid_divs[0]}  {kgrid_divs[1]}  {kgrid_divs[2]}]"
    struct_calc = re.sub(r'kgrid\.MonkhorstPack\s+\[.*?\]', kgrid_str, struct_calc)
    if args.dispersion:
        struct_calc = struct_calc.replace("DFTD3                   .false.", "DFTD3                   .true.")

    with open(os.path.join(struct_dir, "calc.fdf"), 'w') as f:
        f.write(struct_calc)
    for sym in species.keys():
        link_pseudo(args.pp_path, sym, struct_dir)
    print_dual(f"Written: {struct_dir}/structure.fdf, {struct_dir}/calc.fdf", f_out)

    # --- [4] ISOLATED ATOMS ---
    print_section("[4] ISOLATED ATOMS", f_out)
    atoms_root = os.path.join(args.output_dir, "atoms")
    os.makedirs(atoms_root, exist_ok=True)

    atom_rows = []
    for sym, data in species.items():
        atom_dir = os.path.join(atoms_root, sym)
        os.makedirs(atom_dir, exist_ok=True)
        generate_isolated_atom_fdf(sym, data['Z'], os.path.join(atom_dir, "structure.fdf"), args.vacuum)
        with open(os.path.join(atom_dir, "calc.fdf"), 'w') as f:
            f.write(_isolated_atom_calc(args.dispersion))
        link_pseudo(args.pp_path, sym, atom_dir)
        atom_rows.append(([sym, str(data['Z']), f"{atoms_root}/{sym}/"], None))
    print_table(["Species", "Z", "Folder"], atom_rows, f_out)

    # --- [5] BSSE (COUNTERPOISE) CORRECTION ---
    print_section("[5] BSSE (COUNTERPOISE) CORRECTION", f_out)
    if not args.bsse_correction:
        print_dual(color_text(
            "[NOTE] --bsse-correction is OFF -- the resulting cohesive energy will NOT be "
            "corrected for Basis Set Superposition Error (a known LCAO bias that systematically "
            "over-binds). Run without --no-bsse-correction (the default) to get a corrected "
            "reference.", 'yellow'), f_out)
    else:
        # Box side for a ghost cluster at a given cutoff is 2*cutoff + vacuum
        # (see --vacuum's own help text) -- this guarantees a buffer of
        # exactly vacuum/2 beyond the outermost ghost atom REGARDLESS of how
        # large the cutoff is, so --bsse-convergence-check's larger cluster
        # (cutoff + --bsse-convergence-increment) gets a proportionally
        # bigger box instead of eating into a fixed-size one. A plain
        # isolated atom (cutoff=0) reduces to box=vacuum, unchanged from
        # before.
        buffer = args.vacuum / 2.0
        if buffer < 5.0:
            print_dual(color_text(
                f"[WARNING] --vacuum {args.vacuum} Ang gives only {buffer:.1f} Ang of buffer "
                "beyond the outermost atom of every isolated-atom-type calculation (plain atom, "
                "BSSE cluster, and BSSE convergence-check cluster alike -- the box now grows "
                "with --bsse-cutoff/--bsse-convergence-increment to keep this buffer constant, "
                "so increasing either no longer helps) -- consider a larger --vacuum to avoid "
                "interaction with periodic images.", 'yellow'), f_out)

        def check_root(cutoff):
            """'atoms_bsse_check/' for the classic single-point check (exactly
            one --bsse-convergence-increment, unchanged naming/behavior);
            'atoms_bsse_check_<cutoff>/' per point for a multi-cutoff scan."""
            name = "atoms_bsse_check" if not multi_check else f"atoms_bsse_check_{cutoff:.1f}"
            return os.path.join(args.output_dir, name)

        pmg_structure = structure_io.to_pymatgen(structure)
        atoms_bsse_root = os.path.join(args.output_dir, "atoms_bsse")

        any_empty = False
        bsse_rows = []
        space_groups = {}
        for sym, data in species.items():
            fallback_index = next(i for i, site in enumerate(pmg_structure) if site.specie.symbol == sym)
            sites_to_use, all_sites, space_group = resolve_bsse_sites(
                pmg_structure, sym, args.symprec, args.bsse_multi_site, fallback_index)
            space_groups[sym] = space_group
            used_indices = {idx for idx, _w, _m in sites_to_use}

            if len(all_sites) > 1:
                site_desc = ", ".join(f"{w or '?'} (x{m or '?'})" for _i, w, m in all_sites)
                print_dual(f"{color_text('[INFO]', 'cyan')} {len(all_sites)} symmetrically distinct "
                           f"site(s) detected for {sym}: {site_desc}", f_out)
                if not args.bsse_multi_site:
                    skipped = ", ".join(f"{w or '?'} (x{m or '?'})" for _i, w, m in all_sites[1:])
                    print_dual(color_text(
                        f"[WARNING] --no-bsse-multi-site: using only the first site as the BSSE "
                        f"reference for {sym} -- skipping {skipped}, which may not be "
                        "representative if their local environment differs. Pass "
                        "--bsse-multi-site (the default) for a reference per site.", 'yellow'), f_out)

            single_site = len(sites_to_use) == 1
            for idx, wyckoff, mult in sites_to_use:
                if single_site:
                    site_label = "--"
                    bsse_rel = f"atoms_bsse/{sym}/"
                    bsse_dir = os.path.join(atoms_bsse_root, sym)
                    site_suffix = sym
                else:
                    site_label = f"{wyckoff} (x{mult})"
                    bsse_rel = f"atoms_bsse/{sym}/site_{wyckoff}_x{mult}/"
                    bsse_dir = os.path.join(atoms_bsse_root, sym, f"site_{wyckoff}_x{mult}")
                    site_suffix = os.path.join(sym, f"site_{wyckoff}_x{mult}")

                n_ghosts = generate_bsse_reference(
                    pmg_structure, sym, idx, data['Z'], args.bsse_cutoff,
                    2 * args.bsse_cutoff + args.vacuum,
                    bsse_dir, args.pp_path, args.dispersion)
                if n_ghosts == 0:
                    any_empty = True
                bsse_rows.append((
                    [sym, site_label, str(n_ghosts), f"{args.bsse_cutoff} Ang", "yes", bsse_rel],
                    'yellow' if n_ghosts == 0 else None))

                if args.bsse_convergence_check:
                    for cutoff in check_cutoffs:
                        generate_bsse_reference(
                            pmg_structure, sym, idx, data['Z'], cutoff,
                            2 * cutoff + args.vacuum,
                            os.path.join(check_root(cutoff), site_suffix),
                            args.pp_path, args.dispersion)

            for idx, wyckoff, mult in all_sites:
                if idx not in used_indices:
                    site_label = f"{wyckoff} (x{mult})" if wyckoff is not None else "--"
                    bsse_rows.append(([sym, site_label, "--", "--", "no (skipped)", "--"], None))

        print_table(["Species", "Site", "Ghost nb.", "Cutoff", "Used", "Folder"], bsse_rows, f_out)

        if any_empty:
            print_dual(color_text(
                "[WARNING] At least one species/site found ZERO ghost neighbors within "
                f"--bsse-cutoff {args.bsse_cutoff} Ang -- its BSSE-corrected reference will be "
                "IDENTICAL to the uncorrected one (no ghost basis added). Increase "
                "--bsse-cutoff if it does have real neighbors further out.", 'yellow'), f_out)
        print_dual(color_text(
            "How many coordination shells 'Ghost nb.' actually reaches depends on the "
            "material's bond length (see --bsse-cutoff --help) -- the correction is converged "
            "once a larger --bsse-cutoff no longer changes the BSSE-corrected energy "
            "appreciably; consider checking that before trusting the result "
            + ("(see 'atoms_bsse_check/', generated below)." if args.bsse_convergence_check
               else "(--bsse-convergence-check can automate this)."), 'cyan'), f_out)
        if args.bsse_convergence_check and not multi_check:
            c = check_cutoffs[0]
            print_dual(color_text(
                f"Also generated 'atoms_bsse_check/' at cutoff {c:.1f} Ang "
                "(--bsse-cutoff + --bsse-convergence-increment) -- stb-cohesiveAnalysis will "
                "report the shift between the two as a convergence diagnostic. Its box is "
                f"{2 * c + args.vacuum:.1f} Ang wide (vs. "
                f"{2 * args.bsse_cutoff + args.vacuum:.1f} Ang for the --bsse-cutoff reference) "
                f"-- both keep the same {args.vacuum / 2.0:.1f} Ang buffer beyond their "
                "respective outermost ghost atom, so this shift reflects the correction's real "
                "cutoff dependence, not a box-size artifact.", 'cyan'), f_out)
        elif args.bsse_convergence_check and multi_check:
            print_dual(color_text(
                f"Also generated {len(check_cutoffs)} BSSE convergence-check references -- a "
                "full cutoff scan from --bsse-convergence-increment's multiple values -- so "
                "stb-cohesiveAnalysis can report and plot the cohesive energy vs. cutoff trend "
                f"instead of a single before/after point. Every box keeps the same "
                f"{args.vacuum / 2.0:.1f} Ang buffer beyond its own outermost ghost atom (box "
                "side = 2*cutoff + --vacuum), regardless of cutoff, so the trend reflects the "
                "correction's real cutoff dependence, not a box-size artifact.", 'cyan'), f_out)
            print_table(["Cutoff", "Box side", "Folder"], [
                ([f"{c:.1f} Ang", f"{2 * c + args.vacuum:.1f} Ang", f"atoms_bsse_check_{c:.1f}/"], None)
                for c in check_cutoffs
            ], f_out)
        print_dual(color_text(
            "Space group(s) detected (per species, same structure -- should agree): "
            + ", ".join(f"{sym}={sg or '(detection failed)'}" for sym, sg in space_groups.items()),
            'cyan'), f_out)
        print_dual(color_text(
            "stb-cohesiveAnalysis will auto-detect 'atoms_bsse/' and report both the "
            "uncorrected and BSSE-corrected cohesive energy.", 'cyan'), f_out)

    # --- [6] PSEUDOPOTENTIALS ---
    print_section("[6] PSEUDOPOTENTIALS", f_out)
    if args.pp_path:
        found, missing = get_required_pseudos(species.keys(), args.pp_path)
        print_dual(f"Source            : {args.pp_path}", f_out)
        print_table(["Species", "Status"], [
            ([sp, "MISSING" if sp in missing else "found"], 'yellow' if sp in missing else None)
            for sp in species.keys()
        ], f_out)
        if missing:
            print_dual(color_text(
                f"[WARNING] Missing pseudopotential(s) for: {', '.join(missing)} -- these will "
                "need to be added manually to every generated folder.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[OK] All required pseudopotentials found -- linked into every generated "
                "folder.", 'green'), f_out)
    else:
        print_dual("Not given (pass -p/--pp-path -- a bundled bank or a folder path -- to "
                   "link the required pseudopotential for every species into every generated "
                   "folder). Pseudopotentials will need to be added manually.", f_out)

    # --- [7] SUMMARY & NEXT STEPS ---
    print_section("[7] SUMMARY & NEXT STEPS", f_out)
    folder_list = [f"{struct_dir}/", f"{atoms_root}/ ({len(species)} species)"]
    if args.bsse_correction:
        folder_list.append(f"{os.path.join(args.output_dir, 'atoms_bsse')}/")
        if args.bsse_convergence_check:
            if multi_check:
                for c in check_cutoffs:
                    folder_list.append(f"{os.path.join(args.output_dir, f'atoms_bsse_check_{c:.1f}')}/")
            else:
                folder_list.append(f"{os.path.join(args.output_dir, 'atoms_bsse_check')}/")
    print_dual(f"Written under '{args.output_dir}':", f_out)
    for entry in folder_list:
        print_dual(f"  {entry}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every generated folder (structure/, atoms/<species>/"
               + (", atoms_bsse/<species>/..." if args.bsse_correction else "") + ").", f_out)
    print_dual(f"  2. Once they're done: cd {args.output_dir} && stb-cohesiveAnalysis -o calc.out "
               f"(or, from elsewhere: stb-cohesiveAnalysis -o calc.out --dir {args.output_dir}).", f_out)

    if f_out:
        f_out.close()

if __name__ == "__main__":
    main()
