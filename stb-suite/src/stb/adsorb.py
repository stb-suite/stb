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
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from stb.core import structure_io, kspace
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                           capture_library_noise, print_progress_line, finish_progress_line)
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.deps import require_mace
from stb.core.ase_view import view_structure_interactive
from stb.core.adsorption_sites import (
    FIXED_CELL_BLOCK, CONFIG_EXTRA_FILE, resolve_adsorbate, resolve_slab_orientation,
    center_slab_in_vacuum, write_reference_folder, write_site_plot, min_adsorbate_image_distance,
    cluster_candidate_coords, _MIN_LATERAL_IMAGE_SEPARATION_ANG,
    generate_systematic_orientations, deduplicate_orientations,
)
from stb.neb import wrap_into_cell

# REPORT_FILE is always written -- it carries the machine-readable
# "# SITE_TABLE" block stb-adsorbAnalysis::read_site_table depends on, so
# it is never gated behind --save-report (same split as stb-hubbardu's
# always-on run_manifest.json vs. its own --save-report-gated narrative).
# PREP_REPORT_FILE is the optional, human narrative counterpart ([0]-[6]),
# written only with --save-report, matching every other tool's convention.
REPORT_FILE = "adsorption_sites.txt"
PREP_REPORT_FILE = "adsorption_prep_report.txt"
_KGRID_RE = re.compile(r'kgrid\.MonkhorstPack\s+\[.*?\]', re.IGNORECASE)
_SPIN_RE = re.compile(r'(Spin\s+)(\S+)', re.IGNORECASE)


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


def resolve_ml_device(device, mace_relax_module):
    """Returns (device, warning) for a --ml-device request: `device`
    unchanged with warning=None when valid, or ('cpu', <message>) when
    'cuda' was requested but no usable GPU is available. Unlike
    core/mace_relax.py's own resolve_device() (which raises -- the
    deliberate default everywhere else in the suite, so a bad --device cuda
    request never silently gives a much slower run with no indication
    anything is wrong), stb-adsorb's --ml-rank/--ml-prerelax step is only
    ever a fast pre-screening/pre-relax stage ahead of the real SIESTA jobs
    that follow, so degrading to CPU here (with a visible warning in the
    report) is a more convenient default than aborting the whole run.
    Returns the warning text instead of printing it directly, since both
    call sites resolve this from inside a capture_library_noise() block,
    whose stdout redirect would otherwise silently swallow a print_dual
    call made in there.
    """
    if device != "cuda":
        return device, None
    available, detail = mace_relax_module.gpu_available()
    if available:
        return device, None
    return "cpu", f"--ml-device cuda requested, but {detail} -- falling back to CPU."


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


def plot_ml_ranking(scored, ads_name, out_path, show=False):
    """Bar chart of the --ml-rank pre-screen's relative energies (dE = energy
    - lowest energy found) per candidate, one bar per (site, height, type,
    orientation). Single-series version of mladsorb.py::plot_site_ranking
    (no foundation-vs-fine-tuned comparison here, since --ml-rank always
    uses one model). `scored` is a list of (slot, height, site_type, energy,
    relaxed_struct, relaxed_dist, orientation) tuples, already sorted by
    energy (ascending) -- `orientation` is 0 when orientation sampling
    wasn't used (see generate_systematic_orientations), else the 1-based
    rank among that site's unique orientations.

    `show=True` (--view-plots) additionally blocks on plt.show() before the
    figure is closed -- see write_site_plot's own docstring for why this
    isn't just called --view.
    """
    e_min = scored[0][3]
    labels = [f"{st}\n#{slot}@{h:.1f}A" + (f"-o{orient}" if orient else "")
              for slot, h, st, _e, _s, _d, orient in scored]
    dE = [energy - e_min for _slot, _h, _st, energy, _s, _d, _orient in scored]
    fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(labels)), 5))
    colors = ['tab:green' if d == 0.0 else 'tab:blue' for d in dE]
    ax.bar(range(len(labels)), dE, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("dE relative to best site (eV)")
    ax.set_title(f"ML pre-screen ranking -- adsorbate {ads_name}")
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


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
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate O,N,C --site-type ontop\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate O --height-sweep 1.5 3.0 0.5\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                         help="Input slab/2D structure.fdf (vacuum along c). Default: structure.fdf.")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                         help="calc.fdf template already configured for the slab (kgrid, basis, XC, "
                              "%%include structure.fdf, etc.) -- copied as-is into 'clean_slab/' and "
                              "every 'sites/site_*/' folder (same lattice, so the same k-grid applies "
                              "to all of them); a Gamma-only, spin-polarized derivative goes into "
                              "'adsorbate/'. A config_extra.fdf sidecar is %%included on top, forcing "
                              "MD.VariableCell false (see [1]/[4] in the report). Default: calc.fdf.")
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
    parser.add_argument("--force-spin", dest="force_spin", action="store_true", default=True,
                         help="Force Spin polarized (via config_extra.fdf, overriding --calc's own "
                              "Spin tag if any) in every 'sites/site_*/' folder too, not just the "
                              "isolated-adsorbate reference (default: ON) -- a single adsorbate "
                              "atom (or a molecule containing one) bonded to a slab commonly leaves "
                              "the combined system with a net magnetic moment (most simply, an odd "
                              "total valence-electron count, which a spin-restricted calculation "
                              "cannot represent at all). Costs nothing for a genuinely closed-shell "
                              "site (converges to zero moment). Does not affect 'clean_slab/' (no "
                              "adsorbate present).")
    parser.add_argument("--no-force-spin", dest="force_spin", action="store_false",
                         help="Leave --calc's own Spin tag untouched in 'sites/site_*/' (the "
                              "isolated-adsorbate reference is still always forced polarized).")
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

    parser.add_argument("--n-orientations-polar", type=int, default=1,
                         help="With --ml-rank and a multi-atom --adsorbate: systematically sample "
                              "this many initial molecule orientations PER SITE before relaxing -- "
                              "which part of the molecule ends up facing the surface (e.g. O-down "
                              "vs. H-down for H2O) is otherwise fixed by the molecule's default G2 "
                              "geometry alone. Default 1 = current behavior, no sampling. Combined "
                              "with --n-orientations-azimuthal below (total orientations per site = "
                              "polar x azimuthal); a modest 4x2=8 is a reasonable starting point -- "
                              "cost scales linearly with this product (one full MACE relax per "
                              "orientation per site). Works with a single --site-index too, not just "
                              "--all-sites -- every sampled orientation of the SAME site is then what "
                              "gets ranked. Every relaxed orientation (not just the kept/unique ones) "
                              "is saved to sites/orientation_trajectory.xyz as a multi-frame "
                              "trajectory (viewable in VMD/OVITO).")
    parser.add_argument("--n-orientations-azimuthal", type=int, default=1,
                         help="With --ml-rank: this many evenly-spaced in-plane rotations sampled "
                              "per polar direction (see --n-orientations-polar). Default 1.")
    parser.add_argument("--orientation-top-k", type=int, default=None,
                         help="With --ml-rank and orientation sampling: keep only the N best "
                              "unique (post-deduplication) orientations per site, instead of all "
                              "surviving ones. Unset (default) keeps every unique orientation "
                              "found.")
    parser.add_argument("--orientation-rmsd-tol", type=float, default=0.3,
                         help="With --ml-rank and orientation sampling: two relaxed orientations "
                              "at the same site are treated as duplicates (only the better-energy "
                              "one kept) when their adsorbate-atom RMSD is below this (Ang) AND "
                              "their energy differs by less than 0.01 eV. Default 0.3.")

    parser.add_argument("-O", "--output-dir", type=str, default="adsorption_run",
                         help="Root directory (default: adsorption_run) for 'clean_slab', "
                              "'adsorbate' and 'sites' -- everything this study needs lives under "
                              "one folder, same convention as stb-hubbardu's hubbardu_runs.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the full report narrative to {PREP_REPORT_FILE} (under "
                              "--output-dir). Off by default. 'sites/adsorption_sites.txt' (the "
                              "SITE_TABLE stb-adsorbAnalysis reads) is always written regardless.")
    parser.add_argument("--view", action="store_true",
                         help="Open every generated structure (clean_slab, isolated adsorbate(s), "
                              "and every written site) in ASE's interactive 3D viewer as a multi-"
                              "frame browser, after everything else has finished. Needs a display.")
    parser.add_argument("--view-plots", action="store_true",
                         help="Also show the matplotlib plot(s) this run generates -- the "
                              "candidate-site layout, and the ML pre-screen ranking chart if "
                              "--ml-rank is used -- on screen (one blocking window per plot), "
                              "instead of only saving them as PNGs. Off by default. A separate "
                              "flag from --view since that one already means the ASE structure "
                              "viewer here.")
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

    if not args.adsorbate:
        parser.error("--adsorbate is required (unless --list).")
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
    orientation_flags_used = (
        args.n_orientations_polar != 1 or args.n_orientations_azimuthal != 1
        or args.orientation_top_k is not None
    )
    if orientation_flags_used and not args.ml_rank:
        parser.error("--n-orientations-polar/--n-orientations-azimuthal/--orientation-top-k are "
                      "only valid with --ml-rank.")
    if args.n_orientations_polar < 1 or args.n_orientations_azimuthal < 1:
        parser.error("--n-orientations-polar/--n-orientations-azimuthal must be >= 1.")
    if args.ml_rank and not args.all_sites and not orientation_flags_used:
        parser.error("--ml-rank is only valid with --all-sites (screening a single --site-index "
                      "site has nothing to rank against), UNLESS orientation sampling "
                      "(--n-orientations-polar/--n-orientations-azimuthal) is also given -- a "
                      "single site then has plenty to rank: its own sampled orientations.")
    if args.height_sweep is not None:
        hmin, hmax, hstep = args.height_sweep
        if hstep <= 0:
            parser.error("--height-sweep STEP must be > 0.")
        if hmax <= hmin:
            parser.error("--height-sweep MAX must be greater than MIN.")

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
    try:
        structure, relabel_note = resolve_slab_orientation(structure, args.vacuum_gap)
    except ValueError as e:
        print(color_text(f"[ERROR] stb-adsorb {e}", 'red'))
        sys.exit(1)
    if relabel_note:
        print(color_text(f"[INFO] {relabel_note}", 'yellow'))
    pmg_structure = structure_io.to_pymatgen(structure)
    pmg_structure, slab_recentered, recenter_shift_ang = center_slab_in_vacuum(pmg_structure)
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
    site_table_path = os.path.join(sites_root, REPORT_FILE)
    report_path = os.path.join(output_root, PREP_REPORT_FILE) if args.save_report else None
    f_out = open(report_path, 'w') if report_path else None

    library_warnings = []  # collected via capture_library_noise, reported in [6]
    view_frames = []  # (label, ase.Atoms) for --view, built up as folders are prepared

    print_dual(color_text("===== STB-ADSORB REPORT =====", 'magenta'), f_out)

    # --- [0] RUN METADATA: every input that shapes the run. ---
    print_section("[0] RUN METADATA", f_out)
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
    print_dual(f"  ML pre-relax    : {'yes' if args.ml_prerelax else 'no'}", f_out)
    print_dual(f"  ML rank         : {'yes' if args.ml_rank else 'no'}"
                + (f" (top {args.top_k})" if args.ml_rank and args.top_k else ""), f_out)
    print_dual(f"  Vacuum gap/box  : {args.vacuum_gap:.1f} / {args.vacuum_box:.1f} Ang", f_out)
    print_dual(f"  symprec         : {args.symprec}", f_out)
    print_dual(f"  Force spin      : "
                + ("yes -> Spin polarized (isolated ref + site folders)" if args.force_spin
                   else "isolated ref only (--calc as-is for site folders)"), f_out)
    print_dual("  Dipole correction: yes -> Slab.DipoleCorrection T (clean_slab/ + site folders, "
                "not the isolated ref)", f_out)
    print_dual("  van der Waals   : yes -> DFTD3 T (clean_slab/, isolated ref, and site folders, "
                "always on)", f_out)
    print_dual(f"  Save report     : {'yes -> ' + report_path if args.save_report else 'no'}", f_out)

    # --- [1] REFERENCE FOLDERS: clean slab (once) + one isolated-adsorbate
    # reference per requested adsorbate. ---
    print_section("[1] REFERENCE FOLDERS", f_out)
    if slab_recentered:
        print_dual(color_text(
            f"  [INFO] Slab was near a cell boundary along c (uneven vacuum split) -- shifted by "
            f"{recenter_shift_ang:+.3f} Ang to center it in the vacuum gap. This avoids a real "
            "risk: an atom starting near frac z=0 can drift slightly negative during SIESTA "
            "relaxation and wrap around to frac z~1 (the top of the cell) instead, which reads "
            "as an enormous, spurious displacement to anything comparing this relaxed geometry "
            "against another snapshot of the same atom (e.g. [2b]'s own bond-length-change check "
            "below, or stb-adsorbBsse's atom-count/species match).", 'yellow'), f_out)
    else:
        print_dual("  Slab already centered in the vacuum gap (no boundary-wrap risk).", f_out)
    write_reference_folder(clean_slab_dir, pmg_structure, calc_text, slab_species_meta, args.pseudo_dir,
                            force_dipole=True, force_vdw=True)
    print_dual(f"  {color_text('[OK]', 'green')} {clean_slab_dir}", f_out)
    view_frames.append(("clean_slab", AseAtomsAdaptor.get_atoms(pmg_structure)))

    adsorbate_calc_text = force_spin_polarized(force_gamma_kgrid(calc_text))
    adsorbate_dirs = {}  # name -> dir, used below when writing sites
    for i, (name, mol) in enumerate(adsorbates):
        if args.ml_prerelax:
            print_dual(f"  {color_text('ML pre-relax:', 'cyan')} relaxing '{name}' with MACE-MP-0 "
                        "(positions only) ...", f_out)
            # require_mace()/the mace_relax import both trigger mace/torch/e3nn's
            # own import-time prints and warnings the FIRST time this runs in the
            # process -- wrapped together with get_calculator() so none of that
            # leaks out here instead of landing in [6].
            with capture_library_noise(library_warnings, "MACE (--ml-prerelax calculator)"):
                require_mace()
                from stb.core import mace_relax
                ml_device, gpu_warning = resolve_ml_device(args.ml_device, mace_relax)
                calc_mace_ads = mace_relax.get_calculator(model=args.ml_model, device=ml_device)
            if gpu_warning:
                print_dual(color_text(f"  [WARNING] {gpu_warning}", 'yellow'), f_out)
            ads_atoms = AseAtomsAdaptor.get_atoms(mol)
            ads_atoms.pbc = False
            with capture_library_noise(library_warnings, "MACE (--ml-prerelax relax)"):
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

        lateral_dist = min_adsorbate_image_distance(mol.cart_coords, pmg_structure.lattice.matrix)
        if lateral_dist < _MIN_LATERAL_IMAGE_SEPARATION_ANG:
            print_dual(color_text(
                f"  [WARNING] '{name}' on this slab: only {lateral_dist:.2f} Ang from its own "
                f"periodic image in the ab-plane (recommended >= {_MIN_LATERAL_IMAGE_SEPARATION_ANG:.0f} "
                "Ang) -- the slab's in-plane supercell may be too small, letting the adsorbate "
                "spuriously interact with its own lateral copies instead of sitting in true "
                "isolation. Use a bigger in-plane supercell (stb-supercell) if so.", 'yellow'), f_out)
        if args.force_spin:
            print_dual(color_text(
                f"  [NOTE] '{name}': its isolated reference is forced spin-polarized below "
                "(many atoms/molecules have a non-zero ground-state spin, e.g. O, N), and so is "
                "every 'sites/site_*/' folder written in [4] (pass --no-force-spin to use "
                "--calc's own Spin setting there instead). 'clean_slab/' still uses --calc's own "
                "Spin setting unmodified (no adsorbate present, no clear universal reason to "
                "force it).", 'yellow'), f_out)
        else:
            print_dual(color_text(
                f"  [NOTE] '{name}': its isolated reference is forced spin-polarized below "
                "(many atoms/molecules have a non-zero ground-state spin, e.g. O, N), but "
                "--no-force-spin means the COMBINED slab+adsorbate calc.fdf ('sites/site_*/') is "
                "used as given, unmodified. Adsorption often (not always) quenches the "
                "adsorbate's spin -- verify whether Spin polarized is also needed there before "
                "trusting E_ads for an open-shell adsorbate.", 'yellow'), f_out)

        ads_dir = os.path.join(output_root, f"adsorbate_{name}") if multi_adsorbate \
            else os.path.join(output_root, "adsorbate")
        isolated_structure = isolated_adsorbate_structure(mol, args.vacuum_box)
        pmg_isolated = structure_io.to_pymatgen(isolated_structure)
        write_reference_folder(ads_dir, pmg_isolated, adsorbate_calc_text, {}, args.pseudo_dir,
                                force_vdw=True)
        adsorbate_dirs[name] = ads_dir
        print_dual(f"  {color_text('[OK]', 'green')} {ads_dir} (Gamma-only, spin-polarized, "
                    f"{args.vacuum_box:.1f} Ang box)", f_out)
        view_frames.append((f"adsorbate({name})", AseAtomsAdaptor.get_atoms(pmg_isolated)))

    # --- [2] ADSORPTION SITES: FINDING & COUNT -- discovers/counts candidate
    # sites (and, for --position, the manual override point); writes nothing
    # yet. Separated from materialization/writing (now [4]) so the ML pre-
    # screen (now [3]) can report live progress in between the two, instead
    # of being buried mid-way through a single combined section. ---
    print_section("[2] ADSORPTION SITES: FINDING & COUNT", f_out)

    with capture_library_noise(library_warnings, "pymatgen plot_slab"):
        plot_path = os.path.join(sites_root, "adsorption_sites.png")
        write_site_plot(pmg_structure, plot_path, show=args.view_plots)
    print_dual(f"  {color_text('[Saved]', 'cyan')} {plot_path} (sanity check of the candidate "
                "site layout, pymatgen's own defaults -- not necessarily this run's --height/"
                "--symprec/--site-type)", f_out)

    variants = None  # (slot, height, site_type, coord) -- position/auto-find modes

    if args.both_sides:
        # AdsorbateSiteFinder.adsorb_both_surfaces needs a pymatgen Slab (see
        # [4] below for why) -- but the raw-vs-reduced candidate count on the
        # top face alone is a plain find_adsorption_sites call, same as the
        # single-side path, so it's probed here for the same "how many
        # configurations, how many are symmetry-equivalent" reporting.
        (height,) = heights
        with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
            probe_raw = finder.find_adsorption_sites(distance=height, symm_reduce=0, positions=site_types)
            probe_reduced = finder.find_adsorption_sites(distance=height, symm_reduce=args.symprec,
                                                          positions=site_types)
        n_raw, n_reduced = len(probe_raw[site_types[0]]), len(probe_reduced[site_types[0]])
        print_table(["Site type", "Raw candidates", "After symmetry reduction", "Reduced away"],
                    [([site_types[0], str(n_raw), str(n_reduced), str(n_raw - n_reduced)], None)], f_out)
        print_dual("  (top face only -- --both-sides mirrors each selected site onto the bottom "
                    "face too, doubling the final count; the exact mirrored count is confirmed "
                    "in [4] once the symmetry mapping is actually applied.)", f_out)
        n_selected_estimate = n_reduced
    elif args.position is not None:
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
        print_dual("  Manual override (--position): no site-finding performed, exactly 1 "
                    "candidate position per height.", f_out)
        x, y = args.position
        z_ref = float(np.max(pmg_structure.cart_coords[:, 2]))
        mvec = np.asarray(finder.mvec, dtype=float)
        base_point = np.array([x, y, z_ref])
        variants = [(1, h, "manual", base_point + h * mvec) for h in heights]
        if abs(mvec[2]) < 0.999:
            angle_deg = np.degrees(np.arccos(np.clip(abs(mvec[2]), -1.0, 1.0)))
            print_dual(color_text(
                f"  [INFO] Surface normal is not parallel to Cartesian z (tilted "
                f"{angle_deg:.2f} deg) -- --position's height is measured along the true "
                f"normal {tuple(round(float(v), 4) for v in mvec)}, not a plain z-offset.", 'yellow'), f_out)
    else:
        # Candidate site XY coordinates depend only on the slab + height, not
        # on the adsorbate -- found once per height and shared across every
        # requested adsorbate below. A raw (symm_reduce=0, which pymatgen
        # treats as "skip reduction entirely") probe plus the real
        # (symm_reduce=args.symprec) call, one pair per --site-type, gives an
        # explicit "how many candidates, how many are symmetry-equivalent"
        # count -- pymatgen's own find_adsorption_sites already collapses
        # symmetry-equivalent duplicates internally before returning, so
        # without this second probe call there would be no visibility at all
        # into how many were merged away.
        rows, raw_counts, reduced_counts = [], {}, {}
        with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
            for st in site_types:
                raw = finder.find_adsorption_sites(distance=heights[0], symm_reduce=0, positions=[st])
                reduced = finder.find_adsorption_sites(distance=heights[0], symm_reduce=args.symprec,
                                                        positions=[st])
                raw_counts[st], reduced_counts[st] = len(raw[st]), len(reduced[st])
        for st in site_types:
            rows.append(([st, str(raw_counts[st]), str(reduced_counts[st]),
                          str(raw_counts[st] - reduced_counts[st])], None))
        total_raw, total_reduced = sum(raw_counts.values()), sum(reduced_counts.values())
        if len(site_types) > 1:
            rows.append((["TOTAL", str(total_raw), str(total_reduced),
                          str(total_raw - total_reduced)], 'cyan'))
        print_table(["Site type", "Raw candidates", "After symm. reduction", "Reduced away (equivalent)"],
                    rows, f_out)

        candidates_by_height = {}
        with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
            for h in heights:
                found = finder.find_adsorption_sites(distance=h, symm_reduce=args.symprec,
                                                      positions=site_types)
                # Re-wrap every candidate (across all site types) into the same periodic
                # image, anchored on the first one found -- pymatgen's own symm_reduce/
                # put_coord_inside wrap each candidate independently with no shared
                # reference point, so two symmetrically-equivalent candidates can land in
                # different periodic images and appear far apart despite being the same
                # physical site (see cluster_candidate_coords's own docstring).
                anchor = next((c for st in site_types for c in found[st]), None)
                if anchor is not None:
                    for st in site_types:
                        found[st] = cluster_candidate_coords(found[st], pmg_structure.lattice.matrix,
                                                              reference=anchor)
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
        print_dual(f"  Selected: {n_selected} of {n_candidates} symmetrically distinct "
                    f"'{args.site_type}' site(s)"
                    + (" (--all-sites)" if args.all_sites else f" (--site-index {args.site_index})"),
                    f_out)

    if variants is not None:
        # n_slots is always 1 for --position (every variant shares slot=1,
        # one per height) and equals the number of selected sites otherwise
        # -- so this one expression covers both modes without branching.
        n_slots = max((v[0] for v in variants), default=0)
        n_site_folders = n_slots * len(heights) * len(adsorbates)
        print_dual(f"  Configuration count : {n_slots} site(s) x {len(heights)} height(s) x "
                    f"{len(adsorbates)} adsorbate(s) = {n_site_folders} site folder(s), plus "
                    f"clean_slab/ and {len(adsorbates)} adsorbate reference(s)", f_out)

    # --- [3] ML PRE-SCREENING: only meaningful with --ml-rank (fixed section
    # number regardless, so the report layout never shifts). Runs BEFORE any
    # SIESTA folder is written, since it decides (with --top-k) which
    # candidates are worth writing at all. ---
    print_section("[3] ML PRE-SCREENING", f_out)
    ranked_by_adsorbate = {}  # name -> list of (slot, h, st, energy, relaxed_struct, relaxed_dist,
                              # orientation), post-top_k (orientation: 0 = sampling not used, else
                              # 1-based rank among that site's unique orientations)
    if not args.ml_rank:
        print_dual("  Not requested (pass --ml-rank to pre-screen every candidate with MACE-MP-0 "
                    "before committing to real SIESTA folders).", f_out)
    else:
        from ase.constraints import FixAtoms
        import ase.io as ase_io
        # require_mace()/the mace_relax import both trigger mace/torch/e3nn's
        # own import-time prints and warnings the FIRST time this runs in the
        # process -- wrapped together with get_calculator() so none of that
        # leaks out here instead of landing in [6].
        with capture_library_noise(library_warnings, "MACE (--ml-rank calculator)"):
            require_mace()
            from stb.core import mace_relax
            ml_device, gpu_warning = resolve_ml_device(args.ml_device, mace_relax)
            calc_mace = mace_relax.get_calculator(model=args.ml_model, device=ml_device)
        if gpu_warning:
            print_dual(color_text(f"  [WARNING] {gpu_warning}", 'yellow'), f_out)
        orientation_sampling = args.n_orientations_polar > 1 or args.n_orientations_azimuthal > 1
        for name, mol in adsorbates:
            ads_suffix = f"_{name}" if multi_adsorbate else ""
            orient_note = (f", {args.n_orientations_polar}x{args.n_orientations_azimuthal} "
                            "orientation(s)/site" if orientation_sampling else "")
            print_dual(f"  {color_text('Relaxing candidates for', 'cyan')} '{name}' with MACE-MP-0 "
                        f"(substrate fixed, {len(variants)} candidate(s){orient_note}) ...", f_out)
            scored = []
            orientation_frames = []  # every relaxed orientation, every site (before dedup/top-k) --
                                      # only populated when orientation_sampling, written to a .xyz
                                      # trajectory below for visual inspection (e.g. in VMD).
            for i, (slot, h, st, coord) in enumerate(variants, start=1):
                orientations = generate_systematic_orientations(
                    mol, args.n_orientations_polar, args.n_orientations_azimuthal)
                orient_results = []  # (energy, ase_atoms, converged, steps)
                for orient_j, orient_mol in enumerate(orientations, start=1):
                    ads_struct = finder.add_adsorbate(orient_mol, coord)
                    ase_atoms = AseAtomsAdaptor.get_atoms(ads_struct)
                    ase_atoms.set_constraint(FixAtoms(indices=list(range(n_substrate))))
                    with capture_library_noise(library_warnings, "MACE (--ml-rank relax)"):
                        converged, steps = mace_relax.relax(ase_atoms, calc_mace, fmax=args.ml_fmax,
                                                              max_steps=200)
                    energy = ase_atoms.get_potential_energy()
                    orient_results.append((energy, ase_atoms, converged, steps))
                    if orientation_sampling:
                        # Succinct, self-overwriting counter (or periodic full lines when not
                        # a tty, e.g. output redirected to a log) -- the per-orientation MACE
                        # relax below has no other progress feedback, and a full print_dual
                        # line per attempt (as opposed to per KEPT orientation, printed after
                        # dedup below) would flood the report for a large polar x azimuthal grid.
                        print_progress_line(
                            f"    site {i}/{len(variants)} ({st}, h={h:.2f}): orientation "
                            f"{orient_j}/{len(orientations)} relaxed (E = {energy:.4f} eV, "
                            f"{'converged' if converged else 'hit step cap'})",
                            orient_j, len(orientations))
                        ase_atoms.info.update({
                            "site": slot, "site_type": st, "height": round(h, 4),
                            "orientation": orient_j, "energy_eV": round(energy, 6),
                            "converged": bool(converged),
                        })
                        # Detach the calculator before saving -- its own attached results
                        # (e.g. 'energy') collide with ase's extxyz writer trying to fold
                        # them into atoms.info too, on top of the fields set just above.
                        ase_atoms.calc = None
                        orientation_frames.append(ase_atoms)
                if orientation_sampling:
                    finish_progress_line()
                orient_results.sort(key=lambda r: r[0])

                if orientation_sampling:
                    kept = deduplicate_orientations(
                        [(e, a) for e, a, _c, _s in orient_results], n_substrate,
                        rmsd_tol=args.orientation_rmsd_tol)
                    if args.orientation_top_k is not None:
                        kept = kept[:args.orientation_top_k]
                else:
                    kept = [0]

                for rank, idx in enumerate(kept, start=1):
                    energy, ase_atoms, converged, steps = orient_results[idx]
                    relaxed_struct = AseAtomsAdaptor.get_structure(ase_atoms)
                    relaxed_dist = min_adsorbate_slab_distance(relaxed_struct, n_substrate)
                    orient_id = rank if orientation_sampling else 0
                    scored.append((slot, h, st, energy, relaxed_struct, relaxed_dist, orient_id))
                    orient_str = f", orientation {rank}/{len(kept)}" if orientation_sampling else ""
                    print_dual(f"    [{i}/{len(variants)}{orient_str}] site {slot} ({st}, h={h:.2f}): "
                                f"E = {energy:.4f} eV ({'converged' if converged else 'hit step cap'}, "
                                f"{steps} step(s)), relaxed distance = {relaxed_dist:.3f} Ang", f_out)
                if orientation_sampling:
                    print_dual(f"    site {slot} ({st}, h={h:.2f}): {len(orientations)} "
                                f"orientation(s) sampled -> {len(kept)} unique kept"
                                + (f" (--orientation-top-k {args.orientation_top_k})"
                                   if args.orientation_top_k is not None else ""), f_out)
            scored.sort(key=lambda r: r[3])
            e_min = scored[0][3]
            rows = [([str(rank), str(slot), f"{h:.2f}", f"{dist:.3f}", st,
                      (str(orient) if orient else "--"), f"{energy:.4f}", f"{energy - e_min:.4f}"],
                     'green' if rank == 1 else None)
                    for rank, (slot, h, st, energy, _s, dist, orient) in enumerate(scored, start=1)]
            print_table(["Rank", "Site", "Init. h (Ang)", "Relaxed dist (Ang)", "Type",
                         "Orientation", "Energy (eV)", "dE (eV)"], rows, f_out)
            print_dual(color_text(
                "  Note: a relative comparison from a fast ML potential, not an absolute DFT "
                "adsorption energy -- use it to prioritize which site(s) to relax with SIESTA "
                "(stb-adsorbAnalysis), not as a final answer. 'Init. h' is the pre-relaxation "
                "guess used to build each candidate; 'Relaxed dist' is the actual closest "
                "adsorbate-slab distance AFTER this MACE-MP-0 relax (substrate fixed, adsorbate "
                "free) -- this relaxed geometry, not the initial guess, is what gets written to "
                "the SIESTA folders in [4] below. 'Orientation' is '--' when orientation sampling "
                "wasn't requested (--n-orientations-polar/--n-orientations-azimuthal both 1, the "
                "default); otherwise it's the 1-based rank (by energy) among that site's unique, "
                "post-deduplication orientations -- each becomes its own SIESTA folder in [4].",
                'yellow'), f_out)

            plot_path = os.path.join(sites_root, f"ml_rank_ranking{ads_suffix}.png")
            with capture_library_noise(library_warnings, "matplotlib plot_ml_ranking"):
                plot_ml_ranking(scored, name, plot_path, show=args.view_plots)
            print_dual(f"  {color_text('[Saved]', 'cyan')} {plot_path}", f_out)

            if orientation_sampling and orientation_frames:
                traj_path = os.path.join(sites_root, f"orientation_trajectory{ads_suffix}.xyz")
                with capture_library_noise(library_warnings, "ase.io.write (orientation trajectory)"):
                    ase_io.write(traj_path, orientation_frames, format="extxyz")
                print_dual(f"  {color_text('[Saved]', 'cyan')} {traj_path} ({len(orientation_frames)} "
                            "relaxed orientation(s) across every sampled site for this adsorbate, "
                            "including duplicates/discarded ones -- NOT just the unique/top-k subset "
                            "written to SIESTA folders in [4] below. Open as a multi-frame trajectory "
                            "in VMD/OVITO; each frame's info records site/site_type/height/"
                            "orientation/energy/converged.)", f_out)

            if args.top_k is not None:
                scored = scored[:args.top_k]
                print_dual(f"  --top-k {args.top_k}: keeping only the {len(scored)} best-ranked "
                            "site(s) for the SIESTA folders written below.", f_out)
            ranked_by_adsorbate[name] = scored

    # --- [4] WRITING SITE FOLDERS: materializes every selected candidate
    # (from [3]'s ranking if --ml-rank, otherwise fresh from [2]'s variants)
    # and writes it. BSSE ghost-fragment references are NOT generated here
    # any more -- see the [5] note below for why. ---
    print_section("[4] WRITING SITE FOLDERS", f_out)
    if args.ml_rank:
        print_dual("  --ml-rank was used: every site folder below is written from the MACE-MP-0 "
                    "-relaxed geometry reported in [3] ('Relaxed dist'), not the pre-relaxation "
                    "'Init. h' guess -- SIESTA will start from that already-relaxed adsorbate "
                    "position.", f_out)
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
                with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder (both-sides)"):
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
            print_dual(f"  {len(both_structs)} configuration(s) after mirroring onto both faces "
                        f"(top-face probe in [2] found {n_selected_estimate} symmetrically distinct "
                        "site(s)).", f_out)
            ads_suffix = f"_{name}" if multi_adsorbate else ""
            for i, s in enumerate(both_structs, start=1):
                site_records.append(
                    (f"site_{i}_{site_types[0]}_bothsides{ads_suffix}", s, name, height))
    else:
        for name, mol in adsorbates:
            ads_suffix = f"_{name}" if multi_adsorbate else ""
            if args.ml_rank:
                for slot, h, st, _energy, relaxed_struct, _dist, orient in ranked_by_adsorbate[name]:
                    label = f"site_{slot}_{st}{ads_suffix}" + (f"_h{h:.2f}" if multi_height else "") \
                        + (f"_orient{orient}" if orient else "")
                    site_records.append((label, relaxed_struct, name, h))
            else:
                for slot, h, st, coord in variants:
                    ads_struct = finder.add_adsorbate(mol, coord)
                    label = f"site_{slot}_{st}{ads_suffix}" + (f"_h{h:.2f}" if multi_height else "")
                    site_records.append((label, ads_struct, name, h))

    report_rows = []  # (label, adsorbate_name, height, dir)
    for label, ads_struct, ads_name, height in site_records:
        # AdsorbateSiteFinder.add_adsorbate/adsorb_both_surfaces (and the
        # MACE-relaxed structure from --ml-rank, built from the same
        # unwrapped Cartesian placement) place the adsorbate atom(s) via a
        # Cartesian coordinate with no to_unit_cell=True -- confirmed
        # against pymatgen's own source -- so the resulting fractional
        # coordinate can legitimately land outside [0, 1) (e.g. an
        # adsorbate offset pushing z past the cell's own c length). Wrapped
        # here, once, right before it's used for anything (the overlap
        # check below and write_reference_folder), covering all 3 possible
        # origins of ads_struct (plain add_adsorbate, --ml-rank relaxed,
        # --both-sides) since they all funnel through this one loop.
        # min_adsorbate_slab_distance's own distance_matrix is already
        # periodic-aware, so wrapping first doesn't change that result --
        # same wrap_into_cell() already reused (not duplicated) by
        # mlneb.py/mladsorb.py/mldiffusion.py.
        ads_struct = wrap_into_cell(ads_struct)
        min_dist = min_adsorbate_slab_distance(ads_struct, n_substrate)
        if min_dist is not None and min_dist < 0.7:
            print_dual(color_text(
                f"  [WARNING] {label}: closest slab-adsorbate distance is only {min_dist:.3f} Ang "
                "-- likely overlapping atoms (height too small, or a bad --position). Check "
                "before running SIESTA.", 'yellow'), f_out)

        site_dir = os.path.join(sites_root, label)
        write_reference_folder(site_dir, ads_struct, calc_text, slab_species_meta, args.pseudo_dir,
                                force_spin=args.force_spin, force_dipole=True, force_vdw=True)
        print_dual(f"  {color_text('[OK]', 'green')} {site_dir}", f_out)
        report_rows.append((label, ads_name, height, site_dir))
        view_frames.append((label, AseAtomsAdaptor.get_atoms(ads_struct)))

    # --- [5] SUMMARY & NEXT STEPS. ---
    print_section("[5] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"  {len(site_records)} site folder(s) written under '{sites_root}'.", f_out)
    print_dual(f"  Site table      : {site_table_path} (machine-readable, read by stb-adsorbAnalysis)", f_out)
    if report_path:
        print_dual(f"  Report          : {report_path}", f_out)
    print_dual("  Run SIESTA in 'clean_slab/', every 'adsorbate*/' and every 'sites/site_*/' "
                "folder. For a BSSE (counterpoise) -corrected adsorption energy, run "
                "stb-adsorbBsse next -- it needs each site's RELAXED geometry, so it can only "
                "run once these SIESTA jobs are finished, not now. Then use stb-adsorbAnalysis "
                "(with or without BSSE results).", f_out)

    # --- [6] LIBRARY WARNINGS: always last. ---
    print_section("[6] LIBRARY WARNINGS", f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (pymatgen, MACE/torch) during this run -- "
            "collected here instead of interleaved with the report above; harmless in almost "
            "every case, but worth a look if a section above looks suspicious.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    # Machine-readable site table -- always written, independent of
    # --save-report, since stb-adsorbAnalysis::read_site_table depends on it
    # (it skips everything before the "# SITE_TABLE" marker, so the exact
    # format below is unchanged from before this report restructuring).
    with open(site_table_path, 'w') as f_table:
        f_table.write("# Machine-readable site table written by stb-adsorb -- see "
                       f"{PREP_REPORT_FILE if report_path else '--save-report'} (or the console "
                       "log) for the full human-readable report.\n")
        f_table.write("\n# SITE_TABLE -- parsed by stb-adsorbAnalysis, do not reorder the columns\n")
        f_table.write(f"# {'label':<30}{'adsorbate':<14}{'height':<10}{'dir'}\n")
        for label, ads_name, height, site_dir in report_rows:
            f_table.write(f"{label:<32}{ads_name:<14}{height:<10.4f}{site_dir}\n")

    print(f"\n{color_text('Success:', 'green')} {len(site_records)} site folder(s) written under "
          f"'{sites_root}'.")
    print(f"Site table: {site_table_path}")
    if report_path:
        print(f"Full report: {report_path}")

    if args.view:
        print(f"\n{color_text('--view:', 'cyan')} opening {len(view_frames)} frame(s) in ASE's "
              "interactive viewer (use the frame slider/menu to step through them):")
        for i, (label, _atoms) in enumerate(view_frames):
            print(f"  {i} = {label}")
        view_structure_interactive([atoms for _label, atoms in view_frames])


if __name__ == "__main__":
    main()
