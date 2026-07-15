#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import shutil
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.pseudopotentials import resolve_pseudo_source, link_pseudo
from stb.core.deps import require_mace

REPORT_FILE = "adsorption_sites.txt"
_KGRID_RE = re.compile(r'kgrid\.MonkhorstPack\s+\[.*?\]', re.IGNORECASE)
_SPIN_RE = re.compile(r'(Spin\s+)(\S+)', re.IGNORECASE)


def resolve_adsorbate(name):
    """Returns a pymatgen Molecule for `name`: a molecule from ASE's bundled
    G2 database (same technique as stb-molecule.py, which doesn't expose an
    importable function -- ase.build.molecule is called directly here), or a
    single chemical element (one atom at the origin). Tries G2 first since
    G2 names are compound formulas ('H2O', 'CH4', ...) that never collide
    with a bare element symbol.
    """
    from ase.collections import g2
    from ase.build import molecule as ase_build_molecule

    if name in g2.names:
        atoms = ase_build_molecule(name)
        return AseAtomsAdaptor.get_molecule(atoms)

    try:
        Element(name)
    except ValueError:
        raise ValueError(
            f"'{name}' is not a recognized element symbol or ASE G2 molecule name. "
            "Pass --list to see the available G2 molecules."
        )
    return Molecule([name], [[0.0, 0.0, 0.0]])


def reorient_vacuum_to_c(structure, vacuum_axes):
    """Permutes the lattice vectors (and the matching fractional-coordinate
    components) so the single vacuum-padded axis becomes c -- a pure
    relabeling, not a literal 3D rotation: AdsorbateSiteFinder only cares
    about which lattice VECTOR is "out of plane" (mi_vec = cross(a, b)), not
    any particular Cartesian alignment, same convention stb-slab already
    writes (SlabGenerator's reorient_lattice=True, slab.py:21-29). Restores
    a right-handed lattice (positive determinant) if the permutation itself
    would flip handedness, by swapping the two in-plane axes -- the physical
    crystal is identical either way (same "relabeling, not a different
    crystal" spirit as core/symmetry.py's reduce_to_unitcell docstring).

    `vacuum_axes` must have exactly one True entry (caller's responsibility
    -- see resolve_slab_orientation). Returns (new_structure, new_vacuum_axes).
    """
    vacuum_axis = vacuum_axes.index(True)
    order = [i for i in range(3) if i != vacuum_axis] + [vacuum_axis]
    new_lattice = structure.lattice[order]
    if np.linalg.det(new_lattice) < 0:
        order[0], order[1] = order[1], order[0]
        new_lattice = structure.lattice[order]

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    new_frac_coords = frac_coords[:, order]
    new_atoms = [(sym, new_frac_coords[i]) for i, (sym, _pos) in enumerate(structure.atoms)]

    new_structure = structure_io.FdfStructure(
        lattice=new_lattice,
        lattice_constant=structure.lattice_constant,
        species=structure.species,
        species_meta=structure.species_meta,
        atoms=new_atoms,
        coord_format="fractional",
        raw_lines=[],
    )
    return new_structure, [False, False, True]


def resolve_slab_orientation(structure, vacuum_gap):
    """Validates the input has a single well-defined vacuum axis, and
    reorients it to c (see reorient_vacuum_to_c) if it isn't already there.
    Exits with a clear error only for the cases that can't be auto-fixed:
    no vacuum axis at all (a bulk 3D structure, not a slab) or vacuum on 2-3
    axes (a wire or an isolated molecule, no single well-defined surface to
    adsorb onto).
    """
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
    if sum(vacuum_axes) != 1:
        detected = ', '.join(axis for axis, is_vac in zip('abc', vacuum_axes) if is_vac) or 'none'
        print(color_text(
            f"[ERROR] stb-adsorb needs a slab/2D material with vacuum along exactly one "
            f"axis; detected vacuum axis/axes: {detected}. A bulk 3D structure (no vacuum) "
            "or a wire/molecule (vacuum on 2-3 axes) doesn't have a single well-defined "
            "surface to adsorb onto.", 'red'))
        sys.exit(1)

    if not vacuum_axes[2]:
        old_axis = 'abc'[vacuum_axes.index(True)]
        structure, vacuum_axes = reorient_vacuum_to_c(structure, vacuum_axes)
        print(color_text(
            f"[INFO] Input structure's vacuum axis was '{old_axis}', not c -- relabeled the "
            "lattice vectors so vacuum is on c (the convention stb-adsorb/stb-slab/"
            "AdsorbateSiteFinder assume). This is a pure relabeling, not a rotation: the "
            "physical structure is unchanged, but every written structure.fdf below uses "
            "this relabeled a/b/c order, not the input file's original one.", 'yellow'))

    return structure


def force_gamma_kgrid(calc_text):
    """Substitutes the kgrid.MonkhorstPack tag for Gamma-only [1 1 1] --
    always correct for the isolated-adsorbate reference (an isolated
    molecule/atom in a large vacuum box has no real periodicity to sample),
    regardless of what density the caller's calc.fdf template used for the
    slab. Same substitution style as convergence.py's substitute_numeric_tag/
    substitute_kgrid_tag (regex on the tag, not a fixed string match, since
    this template is user-supplied and its exact formatting isn't known).
    """
    new_text, count = _KGRID_RE.subn('kgrid.MonkhorstPack   [1  1  1]', calc_text)
    if count == 0:
        raise ValueError("Could not find a 'kgrid.MonkhorstPack' tag in the calc.fdf template.")
    return new_text


def force_spin_polarized(calc_text):
    """Forces Spin polarized for the isolated-adsorbate reference -- same
    reasoning as cohesive_energy.py's isolated-atom references: a genuinely
    isolated atom (or a molecule containing one, e.g. NO) commonly has a net
    magnetic moment that a spin-restricted calculation gets wrong. Costs
    nothing for a closed-shell molecule (converges to zero moment). Appends
    the tag if the template doesn't have one at all, rather than erroring --
    SIESTA itself defaults to non-polarized when the tag is absent, so an
    absent tag is a normal, valid template, not a malformed one.
    """
    new_text, count = _SPIN_RE.subn(r'\g<1>polarized', calc_text)
    if count == 0:
        return calc_text + "\nSpin                polarized\n"
    return new_text


def isolated_adsorbate_structure(molecule_pmg, vacuum_box):
    """Centers `molecule_pmg` (pymatgen Molecule) in a cubic vacuum box and
    returns it as a periodic FdfStructure -- same convention as
    cohesive_energy.py::generate_isolated_atom_fdf, generalized from one atom
    to any number of atoms. Built as a real pymatgen Structure first (not
    handed to structure_io.from_pymatgen directly as a Molecule) because a
    Molecule's sites have no .frac_coords -- only a genuinely periodic
    Structure with the new box lattice can compute them.
    """
    lattice = Lattice(np.eye(3) * vacuum_box)
    center = molecule_pmg.center_of_mass
    cart_coords = molecule_pmg.cart_coords - center + vacuum_box / 2.0
    species = [site.specie for site in molecule_pmg]
    pmg_structure = Structure(lattice, species, cart_coords, coords_are_cartesian=True)
    return structure_io.from_pymatgen(pmg_structure, coord_format="fractional")


def write_reference_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path):
    """Writes structure.fdf + calc.fdf + linked pseudopotentials for one
    reference/candidate folder (clean_slab/, adsorbate/, or a sites/site_*/).
    Returns the FdfStructure written, so callers that need the same geometry
    again (e.g. to derive BSSE ghost variants) don't have to rebuild it.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta, coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = {site.specie.symbol for site in pmg_structure}
    for sym in sorted(symbols):
        link_pseudo(pp_path, sym, out_dir)
    return fdf_structure


def make_ghost_variant(base_structure, ghost_start, ghost_end):
    """Returns a copy of `base_structure` (an FdfStructure, atoms in a KNOWN
    order: slab atoms first, adsorbate atom(s) appended after -- guaranteed
    by AdsorbateSiteFinder.add_adsorbate/adsorb_both_surfaces, which only
    ever append) with the atoms in [ghost_start, ghost_end) turned into
    ghost species: '<symbol>_ghost' label, negative Z, no valence charge,
    same basis (from the same real pseudopotential file, see link_pseudo's
    dest_label) as the real element -- SIESTA's standard ghost-atom
    convention, already used by cohesive_energy.py's BSSE ghost clusters
    there for "one atom's real local neighbors". Here it's applied to a
    whole fragment (the slab, or the adsorbate) instead of a local
    neighbor shell -- the standard Boys-Bernardi counterpoise scheme for a
    2-fragment (slab + adsorbate) interaction, exact by construction (no
    cutoff to truncate the correction, unlike cohesive_energy.py's --bsse-
    cutoff, since both fragments here are already complete/finite).
    """
    species_meta = dict(base_structure.species_meta)
    new_atoms = []
    for i, (symbol, pos) in enumerate(base_structure.atoms):
        if ghost_start <= i < ghost_end:
            label = f"{symbol}_ghost"
            if label not in species_meta:
                real_z = Element(symbol).Z
                used_ids = {str(info['id']) for info in species_meta.values()}
                next_id = 1
                while str(next_id) in used_ids:
                    next_id += 1
                species_meta[label] = {'id': str(next_id), 'Z': -abs(real_z)}
        else:
            label = symbol
        new_atoms.append((label, pos))

    species = list(dict.fromkeys(sym for sym, _ in new_atoms))
    return structure_io.FdfStructure(
        lattice=base_structure.lattice,
        lattice_constant=base_structure.lattice_constant,
        species=species,
        species_meta=species_meta,
        atoms=new_atoms,
        coord_format=base_structure.coord_format,
        raw_lines=[],
    )


def write_site_plot(pmg_structure, out_path):
    """Saves a top-view PNG of the slab with every ontop/bridge/hollow
    adsorption site marked (pymatgen.analysis.adsorption.plot_slab, already
    confirmed to work on a plain Structure -- no Slab wrapper needed here,
    unlike adsorb_both_surfaces -- checked by reading its source: the
    adsorption-site markers come from AdsorbateSiteFinder.find_adsorption_
    sites() with pymatgen's own defaults, not this run's --height/--symprec/
    --site-type, so treat it as a quick sanity check of the site layout, not
    a literal preview of what will be written below).
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_slab(pmg_structure, ax, adsorption_sites=True)
    ax.set_title("Candidate adsorption sites (ontop/bridge/hollow)")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_bsse_folders(site_dir, site_fdf, n_substrate, calc_text, pp_path):
    """Writes the two ghost-fragment references a Boys-Bernardi counterpoise
    (BSSE) correction of this site's adsorption energy needs: 'bsse_slab/'
    (real slab + ghost adsorbate, same geometry as the site) and
    'bsse_adsorbate/' (ghost slab + real adsorbate, same geometry) --
    stb-adsorbAnalysis then computes E_ads(BSSE) = E_site - E_bsse_slab -
    E_bsse_adsorbate. Same calc.fdf as the site itself (same cell/k-grid,
    no regeneration needed).
    """
    n_total = len(site_fdf.atoms)
    slab_variant = make_ghost_variant(site_fdf, n_substrate, n_total)  # ghost the adsorbate part
    ads_variant = make_ghost_variant(site_fdf, 0, n_substrate)         # ghost the slab part

    for sub_dir, variant in [("bsse_slab", slab_variant), ("bsse_adsorbate", ads_variant)]:
        out_dir = os.path.join(site_dir, sub_dir)
        os.makedirs(out_dir, exist_ok=True)
        structure_io.write_fdf(variant, os.path.join(out_dir, "structure.fdf"))
        with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
            f.write(calc_text)
        present_labels = sorted({symbol for symbol, _ in variant.atoms})
        for label in present_labels:
            real_symbol = label[:-len("_ghost")] if label.endswith("_ghost") else label
            link_pseudo(pp_path, real_symbol, out_dir, dest_label=label)


def min_adsorbate_slab_distance(pmg_structure, n_substrate):
    """Minimum periodic distance between any slab atom (index < n_substrate)
    and any adsorbate atom (index >= n_substrate) in a site structure --
    catches a --height/--position that places the adsorbate unphysically
    close to (or literally inside) the slab. Deliberately narrower than
    structure_io.min_pairwise_distance, which minimizes over ALL pairs
    (including slab-slab and adsorbate-adsorbate ones, not useful here --
    a legitimately short slab-slab bond shouldn't trigger this warning).
    """
    n_total = len(pmg_structure)
    if n_substrate == 0 or n_substrate >= n_total:
        return None
    dm = pmg_structure.distance_matrix
    return float(dm[:n_substrate, n_substrate:].min())


def molecule_extent(molecule_pmg):
    """Largest pairwise distance between any two atoms in an (isolated,
    non-periodic) pymatgen Molecule -- used to sanity-check --vacuum-box
    against the adsorbate's own size (a molecule whose extent approaches
    the box size would self-interact with its periodic images once boxed
    up for the isolated-adsorbate reference calculation).
    """
    if len(molecule_pmg) < 2:
        return 0.0
    coords = molecule_pmg.cart_coords
    diffs = coords[:, None, :] - coords[None, :, :]
    return float(np.linalg.norm(diffs, axis=-1).max())


def build_sweep_values(min_v, max_v, step):
    """Evenly spaced values from min_v to max_v (inclusive) in steps of
    step -- same construction as convergence.py's build_values, duplicated
    (not imported) since it's a 3-line generic numeric helper, not logic
    worth a cross-import between two otherwise-unrelated workflow tools.
    """
    n_steps = int(round((max_v - min_v) / step)) + 1
    return list(np.linspace(min_v, min_v + (n_steps - 1) * step, n_steps))


def parse_adsorbates(spec):
    """Splits a comma-separated --adsorbate value into [(name, Molecule),
    ...], each resolved via resolve_adsorbate. Exits with a clear error
    (naming the offending entry) if any of them doesn't resolve.
    """
    names = [s.strip() for s in spec.split(",") if s.strip()]
    resolved = []
    for name in names:
        try:
            resolved.append((name, resolve_adsorbate(name)))
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)
    return resolved


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Prepares SIESTA folders for an adsorption-energy study: a clean slab/2D "
        "material, an isolated adsorbate, and one folder per candidate adsorption site.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate O --site-type ontop\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate H2O --all-sites --site-type all\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate H --all-sites --ml-rank --top-k 3\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate O --no-bsse-correction\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate O,N,C --site-type ontop\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate O --height-sweep 1.5 3.0 0.5\n"
    )

    parser.add_argument("-s", "--structure", type=str, help="Input slab/2D structure.fdf (vacuum along c).")
    parser.add_argument("-c", "--calc", type=str,
                         help="calc.fdf template already configured for the slab (kgrid, basis, XC, "
                              "%%include structure.fdf, etc.) -- copied as-is into 'clean_slab/' and "
                              "every 'sites/site_*/' folder (same lattice, so the same k-grid applies "
                              "to all of them); a Gamma-only, spin-polarized derivative goes into "
                              "'adsorbate/'.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (optional): a bundled bank or a folder path.")

    parser.add_argument("--adsorbate", type=str, help="Element symbol (single atom) or ASE G2 "
                         "molecule name (e.g. H2O, CO2, CH4). Comma-separated for more than one "
                         "(e.g. 'O,N,C') -- the same site search runs once per adsorbate, useful "
                         "for comparing binding affinity across species at the same site(s).")
    parser.add_argument("--list", action="store_true", help="List the available G2 molecule names and exit.")

    parser.add_argument("--site-type", choices=["ontop", "bridge", "hollow", "all"], default="ontop",
                         help="Which family of adsorption sites to search (default: ontop).")
    parser.add_argument("--site-index", type=int, default=0,
                         help="Which symmetrically distinct site of --site-type to use, 0-based "
                              "(default: 0). Ignored with --all-sites or --position.")
    parser.add_argument("--all-sites", action="store_true",
                         help="Write one folder per symmetrically distinct site of --site-type, "
                              "instead of just --site-index.")
    parser.add_argument("--position", type=float, nargs=2, metavar=("X", "Y"), default=None,
                         help="Manual Cartesian (X, Y) position in Ang, overriding automatic "
                              "site-finding entirely (--site-type/--site-index/--all-sites are "
                              "ignored when this is given).")
    parser.add_argument("--height", type=float, default=2.0,
                         help="Adsorption distance in Ang above the surface (default: 2.0). Ignored "
                              "if --height-sweep is given.")
    parser.add_argument("--height-sweep", type=float, nargs=3, default=None,
                         metavar=("MIN", "MAX", "STEP"),
                         help="Generate one site per height in this range (inclusive), instead of "
                              "a single fixed --height -- an approach-curve scan. "
                              "stb-adsorbAnalysis will report E_ads vs. height for it. Same min/"
                              "max/step convention as stb-convergence. Incompatible with "
                              "--both-sides.")
    parser.add_argument("--both-sides", action="store_true",
                         help="Adsorb on both exposed faces (for a free-standing 2D material with "
                              "vacuum on both sides of the c-axis). Requires a concrete --site-type "
                              "(not 'all'); incompatible with --position, --ml-rank and "
                              "--height-sweep.")
    parser.add_argument("--symprec", type=float, default=0.01,
                         help="Symmetry-reduction tolerance for site-finding (default: 0.01, "
                              "pymatgen's own default).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang, used to validate the input's "
                              "vacuum axis is c (default: 10.0), same convention as stb-kgrid/stb-slab.")
    parser.add_argument("--vacuum-box", type=float, default=20.0,
                         help="Cubic box side (Ang) for the isolated-adsorbate reference calculation "
                              "(default: 20.0).")

    parser.add_argument("--ml-prerelax", action="store_true",
                         help="Relax the isolated adsorbate's geometry (positions only) with "
                              "MACE-MP-0 before writing its reference folder -- same idea as "
                              "stb-phononsCreate --ml-prerelax. Mainly useful for adsorbates whose "
                              "starting geometry isn't already at equilibrium (e.g. a hand-built "
                              "one rather than an ASE G2 entry, which already is). Needs the "
                              "optional 'ml' extra.")

    parser.add_argument("--ml-rank", action="store_true",
                         help="Pre-screen candidate sites with a MACE-MP-0 relax (substrate fixed, "
                              "only the adsorbate moves) before writing SIESTA folders, ranking by "
                              "relaxed energy -- a fast pre-screen for which site(s) to prioritize "
                              "for real DFT, not a replacement for stb-adsorbAnalysis. Needs the "
                              "optional 'ml' extra. Only valid with --all-sites.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--ml-fmax", type=float, default=0.05)
    parser.add_argument("--top-k", type=int, default=None,
                         help="With --ml-rank: only write SIESTA folders for the N best-ranked "
                              "sites, instead of all of them.")

    parser.add_argument("--bsse-correction", dest="bsse_correction", action="store_true", default=True,
                         help="Also generate a BSSE (Basis Set Superposition Error) -corrected "
                              "reference per site, in 'sites/site_*/bsse_slab/' and "
                              "'sites/site_*/bsse_adsorbate/' (default: ON). LCAO adsorption "
                              "energies are otherwise systematically over-bound: the combined "
                              "slab+adsorbate calculation benefits from a larger effective basis "
                              "(each fragment 'borrows' the other's orbitals) than either "
                              "reference computed alone. The correction re-evaluates each "
                              "fragment at the SAME geometry as the site, with the other "
                              "fragment's atoms present as 'ghosts' (SIESTA convention: negative "
                              "Z, no valence charge, same basis as the real element) -- the "
                              "standard Boys-Bernardi counterpoise scheme, exact here (no cutoff "
                              "to truncate, unlike cohesive_energy.py's --bsse-cutoff, since both "
                              "fragments are already complete). Doubles the number of DFT "
                              "calculations needed per site.")
    parser.add_argument("--no-bsse-correction", dest="bsse_correction", action="store_false",
                         help="Skip the BSSE-corrected reference (halves the per-site calculation "
                              "count, at the cost of a systematically over-bound adsorption energy).")

    parser.add_argument("-O", "--output-dir", type=str, default=".",
                         help="Root directory (default: current directory) for 'clean_slab', "
                              "'adsorbate' and 'sites'.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-adsorb {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.list:
        from ase.collections import g2
        print(color_text("Available ASE G2 molecule names:", 'bold'))
        for name in sorted(g2.names):
            print(f"  {name}")
        sys.exit(0)

    if not args.structure or not args.calc or not args.adsorbate:
        parser.error("-s/--structure, -c/--calc and --adsorbate are required (unless --list).")
    if args.both_sides and args.site_type == "all":
        parser.error("--both-sides requires a concrete --site-type (not 'all').")
    if args.both_sides and args.position is not None:
        parser.error("--both-sides and --position are mutually exclusive.")
    if args.both_sides and args.ml_rank:
        parser.error("--both-sides and --ml-rank are not supported together yet.")
    if args.both_sides and args.height_sweep is not None:
        parser.error("--both-sides and --height-sweep are not supported together yet.")
    if args.top_k is not None and not args.ml_rank:
        parser.error("--top-k is only valid with --ml-rank.")
    if args.height_sweep is not None:
        hmin, hmax, hstep = args.height_sweep
        if hstep <= 0:
            parser.error("--height-sweep STEP must be > 0.")
        if hmax <= hmin:
            parser.error("--height-sweep MAX must be greater than MIN.")

    print("\n" + color_text("Prepare an adsorption-energy study:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
        sys.exit(1)

    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    structure = structure_io.read_fdf(args.structure)
    structure = resolve_slab_orientation(structure, args.vacuum_gap)
    pmg_structure = structure_io.to_pymatgen(structure)
    slab_species_meta = structure_io.species_dict(structure)
    n_substrate = len(pmg_structure)

    adsorbates = parse_adsorbates(args.adsorbate)
    multi_adsorbate = len(adsorbates) > 1
    heights = build_sweep_values(*args.height_sweep) if args.height_sweep is not None else [args.height]
    multi_height = len(heights) > 1

    with open(args.calc) as f:
        calc_text = f.read()

    finder = AdsorbateSiteFinder(pmg_structure)
    site_types = ["ontop", "bridge", "hollow"] if args.site_type == "all" else [args.site_type]

    output_root = args.output_dir
    clean_slab_dir = os.path.join(output_root, "clean_slab")
    sites_root = os.path.join(output_root, "sites")
    os.makedirs(sites_root, exist_ok=True)
    report_path = os.path.join(sites_root, REPORT_FILE)

    with open(report_path, 'w') as f_out:
        # --- [0] Run metadata: every input that shapes the run, persisted up
        # front so the report is self-describing even if a later stage exits
        # early (e.g. no valid sites found, --both-sides on an asymmetric slab). ---
        print_dual(f"\n{color_text('[0] RUN METADATA', 'bold')}", f_out)
        print_dual(f"  Structure       : {args.structure}", f_out)
        print_dual(f"  Calc template   : {args.calc}", f_out)
        print_dual(f"  Pseudo dir      : {args.pseudo_dir or '(none)'}", f_out)
        print_dual(f"  Output dir      : {output_root}", f_out)
        for name, mol in adsorbates:
            print_dual(f"  {color_text('Adsorbate:', 'cyan')} {name} ({len(mol)} atom(s))", f_out)
        print_dual(f"  Site type       : {args.site_type}", f_out)
        if multi_height:
            print_dual(f"  {color_text('Height sweep:', 'cyan')} "
                        f"{', '.join(f'{h:.2f}' for h in heights)} Ang", f_out)
        else:
            print_dual(f"  Height          : {heights[0]:.2f} Ang", f_out)
        print_dual(f"  Both sides      : {'yes' if args.both_sides else 'no'}", f_out)
        print_dual(f"  BSSE correction : {'ON' if args.bsse_correction else 'OFF'}", f_out)
        print_dual(f"  ML pre-relax    : {'yes' if args.ml_prerelax else 'no'}", f_out)
        print_dual(f"  ML rank         : {'yes' if args.ml_rank else 'no'}"
                    + (f" (top {args.top_k})" if args.ml_rank and args.top_k else ""), f_out)
        print_dual(f"  Vacuum gap/box  : {args.vacuum_gap:.1f} / {args.vacuum_box:.1f} Ang", f_out)
        print_dual(f"  symprec         : {args.symprec}", f_out)

        plot_path = os.path.join(sites_root, "adsorption_sites.png")
        write_site_plot(pmg_structure, plot_path)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {plot_path} (sanity check of the candidate "
                    "site layout)", f_out)

        # --- Reference folders: clean slab (once) + one isolated-adsorbate
        # reference per requested adsorbate ---
        print_dual(f"\n{color_text('[1] REFERENCE FOLDERS', 'bold')}", f_out)
        write_reference_folder(clean_slab_dir, pmg_structure, calc_text, slab_species_meta, args.pseudo_dir)
        print_dual(f"  {color_text('[OK]', 'green')} {clean_slab_dir}", f_out)

        adsorbate_calc_text = force_spin_polarized(force_gamma_kgrid(calc_text))
        adsorbate_dirs = {}  # name -> dir, used below when writing sites
        for i, (name, mol) in enumerate(adsorbates):
            if args.ml_prerelax:
                require_mace()
                from stb.core import mace_relax
                print_dual(f"  {color_text('ML pre-relax:', 'cyan')} relaxing '{name}' with MACE-MP-0 "
                            "(positions only) ...", f_out)
                calc_mace_ads = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
                ads_atoms = AseAtomsAdaptor.get_atoms(mol)
                ads_atoms.pbc = False
                converged, steps = mace_relax.relax(ads_atoms, calc_mace_ads, fmax=args.ml_fmax, max_steps=200)
                mol = AseAtomsAdaptor.get_molecule(ads_atoms)
                adsorbates[i] = (name, mol)
                print_dual(f"  {'Converged' if converged else 'Hit step cap, not fully converged'} after "
                            f"{steps} step(s).", f_out)

            extent = molecule_extent(mol)
            if extent > 0.5 * args.vacuum_box:
                print_dual(color_text(
                    f"  [WARNING] '{name}' spans {extent:.2f} Ang, more than half of --vacuum-box "
                    f"({args.vacuum_box:.1f} Ang) -- it may self-interact with its own periodic "
                    "images in the isolated-adsorbate reference cell. Increase --vacuum-box.", 'yellow'), f_out)
            if len(mol) == 1:
                print_dual(color_text(
                    f"  [NOTE] '{name}' is a single atom: its isolated reference is forced spin-"
                    "polarized below (many atoms have a non-zero ground-state spin, e.g. O, N), but "
                    "the COMBINED slab+adsorbate calc.fdf is used as given, unmodified. Adsorption "
                    "often (not always) quenches the adsorbate's spin -- verify whether Spin "
                    "polarized is also needed for 'clean_slab/'/'sites/site_*/' before trusting "
                    "E_ads for an open-shell adsorbate.", 'yellow'), f_out)

            ads_dir = os.path.join(output_root, f"adsorbate_{name}") if multi_adsorbate \
                else os.path.join(output_root, "adsorbate")
            isolated_structure = isolated_adsorbate_structure(mol, args.vacuum_box)
            pmg_isolated = structure_io.to_pymatgen(isolated_structure)
            write_reference_folder(ads_dir, pmg_isolated, adsorbate_calc_text, {}, args.pseudo_dir)
            adsorbate_dirs[name] = ads_dir
            print_dual(f"  {color_text('[OK]', 'green')} {ads_dir} (Gamma-only, spin-polarized, "
                        f"{args.vacuum_box:.1f} Ang box)", f_out)

        # --- Candidate adsorption sites: adsorbate x site x height ---
        print_dual(f"\n{color_text('[2] ADSORPTION SITES', 'bold')}", f_out)
        site_records = []  # (label, pmg_structure_with_adsorbate, adsorbate_name, height)

        if args.both_sides:
            # AdsorbateSiteFinder.adsorb_both_surfaces internally calls
            # Slab.get_symmetric_site (to mirror each site onto the opposite
            # face), a method that only exists on pymatgen's Slab class, not on
            # a plain Structure -- confirmed by hitting the AttributeError live
            # and reading pymatgen's source. The other code paths (single-side
            # find_adsorption_sites/add_adsorbate) work fine on a plain
            # Structure, so this wrapping is only needed here. get_symmetric_site
            # itself only needs SpacegroupAnalyzer(self) to work, not any of the
            # Slab-specific provenance metadata, so miller_index/shift/
            # scale_factor below are harmless placeholders -- c is already
            # validated to be the vacuum/surface-normal axis. --height-sweep is
            # rejected earlier (parser.error) for this mode, so heights == [h].
            from pymatgen.core.surface import Slab
            slab_obj = Slab(
                lattice=pmg_structure.lattice,
                species=[site.specie for site in pmg_structure],
                coords=pmg_structure.frac_coords,
                miller_index=(0, 0, 1),
                oriented_unit_cell=pmg_structure,
                shift=0.0,
                scale_factor=np.eye(3),
                reorient_lattice=False,
            )
            (height,) = heights
            for name, mol in adsorbates:
                try:
                    both_structs = AdsorbateSiteFinder(slab_obj).adsorb_both_surfaces(
                        mol, repeat=(1, 1, 1),
                        find_args={"distance": height, "symm_reduce": args.symprec, "positions": site_types})
                except RuntimeError as e:
                    # Slab.get_symmetric_site (see the wrapping note above) raises
                    # this when it can't find a symmetry operation mapping the
                    # top face onto the bottom -- physically, --both-sides needs
                    # a genuinely symmetric slab (a Laue-symmetric, free-standing
                    # 2D material or a centrosymmetric slab), not an arbitrary
                    # one-sided/asymmetric-termination slab.
                    print_dual(color_text(
                        f"[ERROR] Could not place '{name}' on both faces: {e}. --both-sides needs a "
                        "slab with a symmetry operation mapping top to bottom (a free-standing 2D "
                        "material, or a centrosymmetric slab) -- for an asymmetric/one-sided slab, "
                        "adsorb on one face at a time instead.", 'red'), f_out)
                    sys.exit(1)
                if not both_structs:
                    print_dual(color_text(f"[ERROR] No {site_types[0]} sites found.", 'red'), f_out)
                    sys.exit(1)
                ads_suffix = f"_{name}" if multi_adsorbate else ""
                for i, s in enumerate(both_structs, start=1):
                    site_records.append(
                        (f"site_{i}_{site_types[0]}_bothsides{ads_suffix}", s, name, height))
        else:
            # Candidate site XY coordinates depend only on the slab + height, not
            # on the adsorbate -- found once per height and shared across every
            # requested adsorbate below.
            if args.position is not None:
                # Offset along the surface normal (finder.mvec = cross(a, b),
                # normalized -- same vector find_adsorption_sites itself uses
                # internally for every auto-found site's height), not a naive
                # Cartesian-z shift: those only coincide when the slab's in-
                # plane lattice vectors happen to already lie in the xy-plane.
                # For a structure whose in-plane vectors have genuine
                # out-of-plane Cartesian components (e.g. fetched from a
                # database without being pre-aligned), a plain z-offset silently
                # places the adsorbate off-normal -- reproduced live with a
                # deliberately tilted test lattice: mvec came out as
                # [-0.12, -0.16, 0.98], an 11.5-degree error from a naive
                # Cartesian z-hat assumption.
                x, y = args.position
                z_ref = float(np.max(pmg_structure.cart_coords[:, 2]))
                mvec = np.asarray(finder.mvec, dtype=float)
                base_point = np.array([x, y, z_ref])
                selected_by_height = {h: [("manual", base_point + h * mvec)] for h in heights}
                if abs(mvec[2]) < 0.999:
                    angle_deg = np.degrees(np.arccos(np.clip(abs(mvec[2]), -1.0, 1.0)))
                    print_dual(color_text(
                        f"  [INFO] Surface normal is not parallel to Cartesian z (tilted "
                        f"{angle_deg:.2f} deg) -- --position's height is measured along the true "
                        f"normal {tuple(round(float(v), 4) for v in mvec)}, not a plain z-offset.", 'yellow'), f_out)
            else:
                candidates_by_height = {}
                for h in heights:
                    found = finder.find_adsorption_sites(distance=h, symm_reduce=args.symprec,
                                                          positions=site_types)
                    candidates_by_height[h] = [(st, coord) for st in site_types for coord in found[st]]

                n_candidates = len(candidates_by_height[heights[0]])
                if n_candidates == 0:
                    print_dual(color_text(f"[ERROR] No {'/'.join(site_types)} sites found.", 'red'), f_out)
                    sys.exit(1)
                if args.all_sites:
                    selected_indices = list(range(n_candidates))
                else:
                    if args.site_index >= n_candidates:
                        print_dual(color_text(
                            f"[ERROR] --site-index {args.site_index} out of range: only "
                            f"{n_candidates} '{args.site_type}' site(s) found.", 'red'), f_out)
                        sys.exit(1)
                    selected_indices = [args.site_index]
                selected_by_height = {
                    h: [candidates_by_height[h][i] for i in selected_indices] for h in heights
                }

            n_selected = len(selected_by_height[heights[0]])
            variants = []  # (slot, height, site_type, coord)
            for slot in range(n_selected):
                for h in heights:
                    st, coord = selected_by_height[h][slot]
                    variants.append((slot + 1, h, st, coord))

            for name, mol in adsorbates:
                ads_suffix = f"_{name}" if multi_adsorbate else ""

                if args.ml_rank:
                    require_mace()
                    from stb.core import mace_relax
                    from ase.constraints import FixAtoms
                    print_dual(f"  {color_text('ML pre-screen:', 'cyan')} relaxing each candidate for "
                                f"'{name}' with MACE-MP-0 (substrate fixed) ...", f_out)
                    calc_mace = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
                    scored = []
                    for slot, h, st, coord in variants:
                        ads_struct = finder.add_adsorbate(mol, coord)
                        ase_atoms = AseAtomsAdaptor.get_atoms(ads_struct)
                        ase_atoms.set_constraint(FixAtoms(indices=list(range(n_substrate))))
                        mace_relax.relax(ase_atoms, calc_mace, fmax=args.ml_fmax, max_steps=200)
                        energy = ase_atoms.get_potential_energy()
                        relaxed_struct = AseAtomsAdaptor.get_structure(ase_atoms)
                        scored.append((slot, h, st, energy, relaxed_struct))
                    scored.sort(key=lambda r: r[3])
                    e_min = scored[0][3]
                    print_dual(f"\n  {'Rank':<5}{'Site':<6}{'Height':<9}{'Type':<9}{'Energy (eV)':<14}{'dE (eV)':<10}", f_out)
                    for rank, (slot, h, st, energy, _s) in enumerate(scored, start=1):
                        print_dual(f"  {rank:<5}{slot:<6}{h:<9.2f}{st:<9}{energy:<14.4f}{energy - e_min:<10.4f}", f_out)
                    print_dual(color_text(
                        "  Note: a relative comparison from a fast ML potential, not an absolute DFT "
                        "adsorption energy -- use it to prioritize which site(s) to relax with SIESTA "
                        "(stb-adsorbAnalysis), not as a final answer.", 'yellow'), f_out)
                    if args.top_k is not None:
                        scored = scored[:args.top_k]
                    for slot, h, st, _energy, relaxed_struct in scored:
                        label = f"site_{slot}_{st}{ads_suffix}" + (f"_h{h:.2f}" if multi_height else "")
                        site_records.append((label, relaxed_struct, name, h))
                else:
                    for slot, h, st, coord in variants:
                        ads_struct = finder.add_adsorbate(mol, coord)
                        label = f"site_{slot}_{st}{ads_suffix}" + (f"_h{h:.2f}" if multi_height else "")
                        site_records.append((label, ads_struct, name, h))

        report_rows = []  # (label, adsorbate_name, height, dir)
        for label, ads_struct, ads_name, height in site_records:
            min_dist = min_adsorbate_slab_distance(ads_struct, n_substrate)
            if min_dist is not None and min_dist < 0.7:
                print_dual(color_text(
                    f"  [WARNING] {label}: closest slab-adsorbate distance is only {min_dist:.3f} Ang "
                    "-- likely overlapping atoms (height too small, or a bad --position). Check "
                    "before running SIESTA.", 'yellow'), f_out)

            site_dir = os.path.join(sites_root, label)
            site_fdf = write_reference_folder(site_dir, ads_struct, calc_text, slab_species_meta, args.pseudo_dir)
            if args.bsse_correction:
                write_bsse_folders(site_dir, site_fdf, n_substrate, calc_text, args.pseudo_dir)
                print_dual(f"  {color_text('[OK]', 'green')} {site_dir} (+ bsse_slab/, bsse_adsorbate/)", f_out)
            else:
                print_dual(f"  {color_text('[OK]', 'green')} {site_dir}", f_out)
            report_rows.append((label, ads_name, height, site_dir))

        if args.bsse_correction:
            print_dual(color_text(
                "  BSSE (counterpoise) correction: ON -- each site got 2 extra ghost-fragment "
                "references (bsse_slab/, bsse_adsorbate/), tripling its calculation count. "
                "stb-adsorbAnalysis will auto-detect them and report a corrected adsorption "
                "energy alongside the uncorrected one.", 'cyan'), f_out)

        # --- [3] Summary + machine-parseable site table (written directly to
        # f_out, not print_dual, so its fixed-width columns stay exact for
        # stb-adsorbAnalysis::read_site_table -- no ANSI codes ever reach it). ---
        print_dual(f"\n{color_text('[3] SUMMARY', 'bold')}", f_out)
        print_dual(f"  {len(site_records)} site folder(s) written under '{sites_root}'.", f_out)
        print_dual("  Run SIESTA in 'clean_slab/', every 'adsorbate*/' and every 'sites/site_*/' "
                    "folder, then use stb-adsorbAnalysis.", f_out)

        f_out.write("\n# SITE_TABLE -- parsed by stb-adsorbAnalysis, do not reorder the columns\n")
        f_out.write(f"# {'label':<30}{'adsorbate':<14}{'height':<10}{'dir'}\n")
        for label, ads_name, height, site_dir in report_rows:
            f_out.write(f"{label:<32}{ads_name:<14}{height:<10.4f}{site_dir}\n")

    print(f"\n{color_text('Success:', 'green')} {len(site_records)} site folder(s) written under "
          f"'{sites_root}'.")
    print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
