#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.11.0" # D3 dispersion, multi-site BSSE (default on), and --bsse-convergence-check

import os
import sys
import re
import shutil
import argparse
import numpy as np
from stb.core import structure_io, kspace, symmetry
from stb.core.cli import color_text, show_intro
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source

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

def link_pseudo(pp_path, symbol, target_dir, dest_label=None):
    """Symlinks the pseudopotential (.psml preferred, .psf fallback) if the
    directory is provided -- same priority order as
    inputfile.py::copy_pseudopotentials, the other consumer of a pp_path
    folder in this suite.

    `dest_label` (default: `symbol`) is the species LABEL the destination
    filename must match -- SIESTA resolves a species' pseudopotential file by
    its declared ChemicalSpeciesLabel label, not its Z. Used for BSSE ghost
    species (e.g. real pseudopotential 'C.psf' symlinked as 'C_ghost.psf' so
    the ghost label 'C_ghost' -- Z=-6, no valence charge, same basis as real
    carbon -- resolves to the same real pseudopotential file.
    """
    if not pp_path:
        return
    dest_label = dest_label or symbol
    pp_path_abs = os.path.abspath(pp_path)
    for ext in ("psml", "psf"):
        src = os.path.join(pp_path_abs, f"{symbol}.{ext}")
        if os.path.exists(src):
            dst = os.path.join(target_dir, f"{dest_label}.{ext}")
            try:
                os.symlink(src, dst)
            except FileExistsError:
                pass
            return
    print(f"[WARNING] Pseudopotential '{symbol}.psml' or '{symbol}.psf' not found in {pp_path}")
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
    """
    neighbors = find_ghost_neighbors(pmg_structure, anchor_index, cutoff)
    cluster = build_ghost_cluster(sym, z_num, neighbors, vacuum)
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(cluster, os.path.join(out_dir, "structure.fdf"))

    with open(os.path.join(out_dir, "calc.fdf"), 'w') as f:
        f.write(_isolated_atom_calc(dispersion))

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


def write_inequivalent_sites_report(path, struct_file, symprec, bsse_cutoff, check_cutoff,
                                     convergence_check, multi_site, space_groups, rows):
    """Writes the full non-equivalent-site / BSSE-reference table (species,
    Wyckoff letter, multiplicity, fractional coordinates, whether it was
    actually used as a BSSE reference, and its ghost-neighbor count) plus
    the symmetry context (space group per species, symprec, cutoffs) to a
    plain-text file at `path`. Complements the same table already printed
    to stdout during setup -- this is the persistent copy so it survives
    past the terminal, for citing in a paper/report or re-checking later
    which Wyckoff site a given atoms_bsse/<sym>/site_*/ folder corresponds
    to.
    """
    lines = []
    lines.append("=" * 80)
    lines.append("NON-EQUIVALENT SITES / BSSE REFERENCE REPORT")
    lines.append("=" * 80)
    lines.append(f"Structure file                             : {struct_file}")
    lines.append(f"Symmetry tolerance (--symprec)              : {symprec}")
    lines.append(f"BSSE ghost-neighbor cutoff (--bsse-cutoff)  : {bsse_cutoff} Ang")
    if convergence_check:
        lines.append(f"BSSE convergence check (--bsse-convergence-check): ON "
                      f"(check cutoff {check_cutoff} Ang)")
    else:
        lines.append("BSSE convergence check                      : OFF")
    lines.append(f"Multi-site BSSE (--bsse-multi-site)         : {'ON' if multi_site else 'OFF'}")
    lines.append("")
    lines.append("Space group(s) detected (per species, same structure -- should agree):")
    for sym, sg in space_groups.items():
        lines.append(f"  {sym}: {sg or '(detection failed/unavailable)'}")
    lines.append("")

    header = f"{'Species':<10}{'Wyckoff':<10}{'Mult.':<8}{'Fractional coordinates':<32}{'Used':<7}{'Ghost nb.'}"
    sep = "-" * 80
    lines.append(sep)
    lines.append(header)
    lines.append(sep)
    for sym, wyckoff, mult, frac_coords, used, n_ghosts, _bsse_rel in rows:
        wy = wyckoff or "?"
        m = str(mult) if mult is not None else "?"
        coord_str = f"{frac_coords[0]:.6f}  {frac_coords[1]:.6f}  {frac_coords[2]:.6f}"
        used_str = "yes" if used else "no"
        ghost_str = str(n_ghosts) if n_ghosts is not None else "--"
        lines.append(f"{sym:<10}{wy:<10}{m:<8}{coord_str:<32}{used_str:<7}{ghost_str}")
    lines.append(sep)
    lines.append("")

    used_rows = [r for r in rows if r[4]]
    if used_rows:
        lines.append("BSSE reference directories:")
        for sym, wyckoff, mult, _frac, _used, _n, bsse_rel in used_rows:
            site_desc = sym if wyckoff is None else f"{sym}, site {wyckoff} (x{mult})"
            entry = f"  {site_desc}: {bsse_rel}"
            if convergence_check:
                entry += f"  (check: {bsse_rel.replace('atoms_bsse/', 'atoms_bsse_check/', 1)})"
            lines.append(entry)
        lines.append("")

    with open(path, 'w') as f:
        f.write("\n".join(lines) + "\n")


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
                        help="Cubic box side (Ang) for each isolated-atom calculation "
                             "(default: 20.0). Increase for elements with unusually "
                             "diffuse basis orbitals.")

    parser.add_argument("--vacuum-gap", dest="vacuum_gap", type=float, default=10.0,
                        help="Vacuum gap threshold in Ang used to detect which lattice axes "
                             "of the FULL structure are periodic vs. vacuum-padded (default: "
                             "10.0), same convention as stb-kgrid/stb-strain/stb-elasticInputs. "
                             "Vacuum-padded axes are forced to a single (Gamma-only) k-grid "
                             "division regardless of --k-density.")

    parser.add_argument("-O", "--output-dir", dest="output_dir", default=".",
                        help="Root directory (default: current directory) for the 'structure', "
                             "'atoms' and 'atoms_bsse' folders -- set this to avoid silently "
                             "overwriting a previous run's setup when preparing more than one "
                             "structure from the same working directory.")

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
                        help="Also generate a second BSSE reference per species/site at a "
                             "larger cutoff (see --bsse-convergence-increment), in "
                             "'atoms_bsse_check/' -- stb-cohesiveAnalysis then reports how "
                             "much the BSSE-corrected cohesive energy shifts between the two "
                             "cutoffs, as a convergence check (default: OFF -- on top of "
                             "--bsse-correction's own doubling, this triples the isolated-"
                             "atom-type calculation count). Requires --bsse-correction.")
    parser.add_argument("--bsse-convergence-increment", dest="bsse_convergence_increment",
                        type=float, default=2.0,
                        help="Ang added to --bsse-cutoff for the --bsse-convergence-check "
                             "reference (default: 2.0).")

    parser.add_argument("--dispersion", dest="dispersion", action="store_true",
                        help="Enable Grimme DFT-D3 dispersion correction (SIESTA's own "
                             "DFTD3 directive, default: OFF) for every calculation (full "
                             "structure, isolated atoms, and BSSE references alike) -- relevant "
                             "for layered/molecular materials where van der Waals interactions "
                             "contribute non-negligibly to the cohesive energy. Safe to apply "
                             "uniformly: a genuinely isolated atom has no neighbor to form a "
                             "dispersion pair with, so D3 contributes ~0 there regardless.")

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
        if args.bsse_convergence_increment <= 0:
            parser.error("--bsse-convergence-increment must be positive.")

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

    print("\n" + color_text("COHESIVE ENERGY:", 'bold'))
    print("-"*60)

    # Extract structure data
    print("\n[INFO] Read structure file ...")
    structure, vacuum_axes = parse_structure_fdf(args.structure, args.vacuum_gap)
    print(f"{color_text('[INFO]', 'green')} Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}")
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
    print(f"[INFO] Detected species: {', '.join(species.keys())}")
    if unused_species:
        print(f"{color_text('[INFO]', 'cyan')} Declared but never placed (no isolated-atom "
              f"calculation needed): {', '.join(unused_species)}")

    # Calculate K-grid for the full structure
    print("[INFO] Calculate Monkhorst-Pack grid ...")
    try:
        kgrid_divs = kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], args.k_density, vacuum_axes)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)
    print(f"[INFO] Calculated K-grid for full structure: {kgrid_divs[0]} {kgrid_divs[1]} {kgrid_divs[2]} (density = {args.k_density})")

    # 1. Setup full structure directory
    print("\n[INFO] Setting up full structure directory ...")
    struct_dir = os.path.join(args.output_dir, "structure")
    os.makedirs(struct_dir, exist_ok=True)

    shutil.copy(args.structure, os.path.join(struct_dir, "structure.fdf"))

    struct_calc = CALC_TEMPLATE
    if args.spin:
        struct_calc = struct_calc.replace("Spin                non-polarized", "Spin                polarized")
        print("[INFO] Spin polarization ENABLED for the full structure.")
    else:
        print("[INFO] Spin polarization DISABLED for the full structure.")
        print(color_text(
            "[NOTE] Every isolated-atom (and BSSE-corrected) reference is always spin-"
            "polarized regardless of this setting -- correct for closed-shell/non-magnetic "
            "bulk materials, but if this structure is expected to order magnetically, its "
            "non-polarized bulk energy will be compared against polarized atomic references, "
            "an inconsistent (not just incomplete) comparison. Pass --spin if in doubt.",
            'yellow'))

    # Apply calculated K-grid
    kgrid_str = f"kgrid.MonkhorstPack   [{kgrid_divs[0]}  {kgrid_divs[1]}  {kgrid_divs[2]}]"
    struct_calc = re.sub(r'kgrid\.MonkhorstPack\s+\[.*?\]', kgrid_str, struct_calc)

    if args.dispersion:
        struct_calc = struct_calc.replace("DFTD3                   .false.", "DFTD3                   .true.")
        print(f"{color_text('[INFO]', 'cyan')} DFT-D3 dispersion correction ENABLED.")

    with open(os.path.join(struct_dir, "calc.fdf"), 'w') as f:
        f.write(struct_calc)

    for sym in species.keys():
        link_pseudo(args.pp_path, sym, struct_dir)

    # 2. Setup isolated atoms directories
    print("\n[INFO] Setting up isolated atoms directories ...")
    atoms_root = os.path.join(args.output_dir, "atoms")
    os.makedirs(atoms_root, exist_ok=True)

    for sym, data in species.items():
        print(f"[INFO] Setting up isolated {sym} ...")
        atom_dir = os.path.join(atoms_root, sym)
        os.makedirs(atom_dir, exist_ok=True)

        # Isolated atom structure
        generate_isolated_atom_fdf(sym, data['Z'], os.path.join(atom_dir, "structure.fdf"), args.vacuum)

        # Calc file for isolated atom (Always polarized + Gamma point only)
        with open(os.path.join(atom_dir, "calc.fdf"), 'w') as f:
            f.write(_isolated_atom_calc(args.dispersion))

        link_pseudo(args.pp_path, sym, atom_dir)

    # 3. BSSE (counterpoise) -corrected reference: one ghost cluster per
    # species, or per symmetrically distinct site if the species occupies
    # more than one (see --bsse-multi-site).
    if args.bsse_correction:
        print("\n" + "-"*60)
        print(color_text("BSSE (COUNTERPOISE) CORRECTION", 'cyan').center(60))
        print("-"*60)
        half_box = args.vacuum / 2.0
        if half_box - args.bsse_cutoff < 5.0:
            print(color_text(
                f"[WARNING] --vacuum {args.vacuum} Ang leaves only "
                f"{half_box - args.bsse_cutoff:.1f} Ang of buffer beyond the ghost cluster's "
                f"outermost neighbor (--bsse-cutoff {args.bsse_cutoff} Ang) -- consider a "
                "larger --vacuum or a smaller --bsse-cutoff to avoid the cluster interacting "
                "with its own periodic images.", 'yellow'))

        check_cutoff = args.bsse_cutoff + args.bsse_convergence_increment
        if args.bsse_convergence_check and half_box - check_cutoff < 5.0:
            print(color_text(
                f"[WARNING] --vacuum {args.vacuum} Ang leaves only "
                f"{half_box - check_cutoff:.1f} Ang of buffer for the --bsse-convergence-check "
                f"cluster (cutoff {check_cutoff} Ang) -- consider a larger --vacuum.", 'yellow'))

        pmg_structure = structure_io.to_pymatgen(structure)
        atoms_bsse_root = os.path.join(args.output_dir, "atoms_bsse")
        atoms_bsse_check_root = os.path.join(args.output_dir, "atoms_bsse_check")

        print(f"  {'Species':<10}{'Site':<14}{'Ghost neighbors':<18}Cutoff")
        print(f"  {'-'*55}")
        any_empty = False
        report_rows = []
        space_groups = {}
        for sym, data in species.items():
            fallback_index = next(i for i, site in enumerate(pmg_structure) if site.specie.symbol == sym)
            sites_to_use, all_sites, space_group = resolve_bsse_sites(
                pmg_structure, sym, args.symprec, args.bsse_multi_site, fallback_index)
            space_groups[sym] = space_group
            used_indices = {idx for idx, _w, _m in sites_to_use}

            if len(all_sites) > 1:
                site_desc = ", ".join(f"{w or '?'} (x{m or '?'})" for _i, w, m in all_sites)
                print(f"{color_text('[INFO]', 'cyan')} {len(all_sites)} symmetrically distinct "
                      f"site(s) detected for {sym}: {site_desc}")
                if not args.bsse_multi_site:
                    skipped = ", ".join(f"{w or '?'} (x{m or '?'})" for _i, w, m in all_sites[1:])
                    print(color_text(
                        f"[WARNING] --no-bsse-multi-site: using only the first site as the BSSE "
                        f"reference for {sym} -- skipping {skipped}, which may not be "
                        "representative if their local environment differs. Pass "
                        "--bsse-multi-site (the default) for a reference per site.", 'yellow'))

            single_site = len(sites_to_use) == 1
            for idx, wyckoff, mult in sites_to_use:
                if single_site:
                    site_label = "--"
                    bsse_rel = f"atoms_bsse/{sym}/"
                    bsse_dir = os.path.join(atoms_bsse_root, sym)
                    check_dir = os.path.join(atoms_bsse_check_root, sym)
                else:
                    site_label = f"{wyckoff} (x{mult})"
                    bsse_rel = f"atoms_bsse/{sym}/site_{wyckoff}_x{mult}/"
                    bsse_dir = os.path.join(atoms_bsse_root, sym, f"site_{wyckoff}_x{mult}")
                    check_dir = os.path.join(atoms_bsse_check_root, sym, f"site_{wyckoff}_x{mult}")

                n_ghosts = generate_bsse_reference(
                    pmg_structure, sym, idx, data['Z'], args.bsse_cutoff, args.vacuum,
                    bsse_dir, args.pp_path, args.dispersion)
                print(f"  {sym:<10}{site_label:<14}{n_ghosts:<18}{args.bsse_cutoff} Ang")
                if n_ghosts == 0:
                    any_empty = True
                report_rows.append((sym, wyckoff, mult, pmg_structure[idx].frac_coords, True, n_ghosts, bsse_rel))

                if args.bsse_convergence_check:
                    generate_bsse_reference(
                        pmg_structure, sym, idx, data['Z'], check_cutoff, args.vacuum,
                        check_dir, args.pp_path, args.dispersion)

            for idx, wyckoff, mult in all_sites:
                if idx not in used_indices:
                    report_rows.append((sym, wyckoff, mult, pmg_structure[idx].frac_coords, False, None, None))
        print(f"  {'-'*55}")
        if any_empty:
            print(color_text(
                "  [WARNING] At least one species/site found ZERO ghost neighbors within "
                f"--bsse-cutoff {args.bsse_cutoff} Ang -- its BSSE-corrected reference will be "
                "IDENTICAL to the uncorrected one (no ghost basis added). Increase "
                "--bsse-cutoff if it does have real neighbors further out.", 'yellow'))
        print(color_text(
            "  How many coordination shells 'Ghost neighbors' actually reaches depends on the "
            "material's bond length (see --bsse-cutoff --help) -- the correction is converged "
            "once a larger --bsse-cutoff no longer changes the BSSE-corrected energy "
            "appreciably; consider checking that before trusting the result "
            + ("(see 'atoms_bsse_check/', generated below)." if args.bsse_convergence_check
               else "(--bsse-convergence-check can automate this).") , 'cyan'))
        if args.bsse_convergence_check:
            print(color_text(
                f"  Also generated 'atoms_bsse_check/' at cutoff {check_cutoff} Ang "
                "(--bsse-cutoff + --bsse-convergence-increment) -- stb-cohesiveAnalysis will "
                "report the shift between the two as a convergence diagnostic.", 'cyan'))
        print(color_text(
            "  stb-cohesiveAnalysis will auto-detect 'atoms_bsse/' and report both the "
            "uncorrected and BSSE-corrected cohesive energy.", 'cyan'))

        report_path = os.path.join(atoms_bsse_root, "inequivalent_sites.txt")
        write_inequivalent_sites_report(
            report_path, args.structure, args.symprec, args.bsse_cutoff, check_cutoff,
            args.bsse_convergence_check, args.bsse_multi_site, space_groups, report_rows)
        print(color_text(f"  Non-equivalent site / symmetry report written to: {report_path}", 'cyan'))
        print("-"*60)
    else:
        print(color_text(
            "\n[NOTE] --bsse-correction is OFF -- the resulting cohesive energy will NOT be "
            "corrected for Basis Set Superposition Error (a known LCAO bias that "
            "systematically over-binds). Run without --no-bsse-correction (the default) to "
            "get a corrected reference.", 'yellow'))

    print("\n[INFO] Complete job!")
    print("\n"+"-"*60)
    print(color_text("Setup complete! Folders 'structure' and 'atoms' are ready.\n\n", 'bold'))

if __name__ == "__main__":
    main()
