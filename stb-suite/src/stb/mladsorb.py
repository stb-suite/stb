#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Standalone MACE-driven adsorption-site screening -- the tenth tool in
the ML Simulations category. stb-adsorb already has a `--ml-rank` mode
(enumerate candidate ontop/bridge/hollow sites, relax each with MACE-MP-0
substrate-fixed, rank by energy) used to prioritize which site(s) to
actually send to SIESTA, but it's bolted onto a DFT-oriented workflow: it
always needs `--calc`/pseudopotentials even though the ML path never
touches either, only reports a RELATIVE energy (dE from the best site, not
an absolute adsorption energy), and its `--ml-model` only accepts
`small`/`medium`/`large` (no `--custom-model`, the same gap stb-mlneb
found and fixed for stb-neb). This tool is the pure-MACE extraction of
that same capability, additionally computing a genuine ML adsorption
energy per site (E_ads = E_site - E_slab - E_isolated_adsorbate, all three
MACE energies) rather than just a relative ranking number.

Deliberately reuses (not duplicates) adsorb.py's own pure-geometry helpers
-- resolve_slab_orientation, parse_adsorbates (itself built on
resolve_adsorbate), isolated_adsorbate_structure, molecule_extent,
write_site_plot, min_adsorbate_slab_distance -- since finding candidate
sites/building the isolated-adsorbate reference is exactly the same
problem whether the result goes to SIESTA or straight into a MACE
relaxation; only write_reference_folder/write_bsse_folders (SIESTA input/
pseudopotential-specific) aren't reused, since this tool writes no DFT
input at all.

Four features on top of the original single-model site ranking:
- Foundation-model comparison (--custom-model + overlay/table against the
  raw foundation model), same pattern as stb-mlphonons/stb-mlelastic.
- --n-orientations: for a multi-atom adsorbate, cheaply screens several
  random rotations at each site (single-point energy, no relax) before
  the expensive full relax, keeping the best-looking one -- the molecule's
  own z-axis is always re-aligned to the surface normal by
  AdsorbateSiteFinder.add_adsorbate regardless of input orientation (see
  its own docstring), so this samples the physically meaningful remaining
  degree of freedom (rotation about that normal / out-of-plane tilt).
- --height-sweep now also plots a proper binding-energy-vs-height curve
  per site (plot_height_curves), not just a flat ranking of every
  (site, height) pair together.
- --diffusion-barrier: estimates the adsorbate's own surface-diffusion
  barrier between its 2 most stable, distinct sites via a real climbing-
  image NEB, reusing stb-mlneb's run_single_neb directly (not
  duplicated) -- distinct from stb-mldiffusion (vacancy migration in the
  bulk), this is migration of an ADSORBATE across a surface.
"""

VERSION = "1.1.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase.constraints import FixAtoms
from pymatgen.core import Molecule
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.adsorb import (
    resolve_slab_orientation, parse_adsorbates, isolated_adsorbate_structure,
    molecule_extent, write_site_plot, min_adsorbate_slab_distance, build_sweep_values,
)
from stb.neb import wrap_into_cell
from stb.mlneb import run_single_neb

require_mace()

REPORT_FILE = "stb_mladsorb_report.txt"
_SERIES_COLORS = ['tab:blue', 'tab:orange']


def pick_best_orientation(finder, mol, coord, calc, n_orientations, rng):
    """Tries `n_orientations` random 3D rotations of `mol` about its own
    center of mass, does a cheap SINGLE-POINT energy evaluation of each at
    `coord` (no relax -- just screening which orientation looks best before
    the expensive full relax), and returns the best-scoring rotated
    Molecule. AdsorbateSiteFinder.add_adsorbate itself always re-centers/
    re-aligns the molecule's own z-axis to the surface normal (see its own
    docstring) regardless of the input orientation, so this samples the
    physically meaningful remaining degree of freedom: rotation about that
    normal, and any out-of-plane tilt for a non-linear molecule.
    """
    from scipy.spatial.transform import Rotation

    best_energy, best_mol = None, mol
    for _ in range(n_orientations):
        rot = Rotation.random(rng=rng)
        center = mol.center_of_mass
        rotated_coords = rot.apply(mol.cart_coords - center) + center
        rotated_mol = Molecule(mol.species, rotated_coords)
        trial_struct = finder.add_adsorbate(rotated_mol, coord)
        trial_atoms = AseAtomsAdaptor.get_atoms(trial_struct)
        trial_atoms.calc = calc
        energy = trial_atoms.get_potential_energy()
        if best_energy is None or energy < best_energy:
            best_energy, best_mol = energy, rotated_mol
    return best_mol


def screen_sites(pmg_structure, n_substrate, adsorbates, heights, site_types, args, calc,
                  rng, f_out, tag=""):
    """Runs the full site-screening pipeline (reference energies + every
    candidate site relaxed, substrate fixed) with ONE calculator. Returns
    (scored, e_slab, e_isolated): scored is a list of (slot, height, name,
    site_type, e_ads, ase_atoms) sorted ascending by e_ads (most stable
    first).
    """
    prefix = f"[{tag}] " if tag else ""
    slab_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    slab_atoms.calc = calc
    e_slab = slab_atoms.get_potential_energy()
    print_dual(f"{prefix}Clean slab energy : {e_slab:.6f} eV (single-point, as given)", f_out)

    e_isolated = {}
    for name, mol in adsorbates:
        extent = molecule_extent(mol)
        if extent > 0.5 * args.vacuum_box:
            print_dual(color_text(
                f"{prefix}[WARNING] '{name}' spans {extent:.2f} Ang, more than half of "
                f"--vacuum-box ({args.vacuum_box:.1f} Ang) -- may self-interact with its own "
                "periodic images. Increase --vacuum-box.", 'yellow'), f_out)
        isolated_structure = isolated_adsorbate_structure(mol, args.vacuum_box)
        pmg_isolated = structure_io.to_pymatgen(isolated_structure)
        iso_atoms = AseAtomsAdaptor.get_atoms(pmg_isolated)
        iso_atoms.calc = calc
        converged, steps = mace_relax.relax(iso_atoms, calc, fmax=args.fmax, max_steps=args.max_steps)
        e_isolated[name] = iso_atoms.get_potential_energy()
        print_dual(f"{prefix}Isolated '{name}' energy: {e_isolated[name]:.6f} eV "
                   f"({'converged' if converged else 'hit step cap'}, {steps} steps)", f_out)

    finder = AdsorbateSiteFinder(pmg_structure)
    candidates_by_height = {}
    for h in heights:
        found = finder.find_adsorption_sites(distance=h, symm_reduce=args.symprec, positions=site_types)
        candidates_by_height[h] = [(st, coord) for st in site_types for coord in found[st]]
    n_candidates = len(candidates_by_height[heights[0]])
    if n_candidates == 0:
        return None, e_slab, e_isolated

    scored = []
    for name, mol in adsorbates:
        for h in heights:
            for slot, (st, coord) in enumerate(candidates_by_height[h], start=1):
                mol_variant = mol
                if args.n_orientations > 1 and len(mol) > 1:
                    mol_variant = pick_best_orientation(finder, mol, coord, calc,
                                                        args.n_orientations, rng)
                ads_struct = finder.add_adsorbate(mol_variant, coord)
                ase_atoms = AseAtomsAdaptor.get_atoms(ads_struct)
                min_dist = min_adsorbate_slab_distance(ads_struct, n_substrate)
                if min_dist is not None and min_dist < 0.7:
                    print_dual(color_text(
                        f"{prefix}[WARNING] site {slot} ({st}, h={h:.2f}, {name}): closest "
                        f"slab-adsorbate distance is only {min_dist:.3f} Ang before "
                        "relaxation.", 'yellow'), f_out)
                ase_atoms.calc = calc
                ase_atoms.set_constraint(FixAtoms(indices=list(range(n_substrate))))
                mace_relax.relax(ase_atoms, calc, fmax=args.fmax, max_steps=args.max_steps)
                energy = ase_atoms.get_potential_energy()
                e_ads = energy - e_slab - e_isolated[name]
                scored.append((slot, h, name, st, e_ads, ase_atoms))

    scored.sort(key=lambda r: r[4])
    return scored, e_slab, e_isolated


def run_diffusion_barrier(scored, adsorbate_name, args, calc, f_out):
    """After ranking, estimates the adsorbate's own surface-diffusion
    barrier between its 2 most stable, DISTINCT sites (different `slot`,
    so a genuine hop rather than the same site at two different heights)
    via a real climbing-image NEB -- reusing stb-mlneb's run_single_neb
    directly. freeze_substrate=True is always used here: the two site
    endpoints already share bit-identical substrate positions (both site
    relaxations only ever moved the adsorbate, substrate held fixed by
    FixAtoms in screen_sites), so run_single_neb's own compute_frozen_
    indices auto-detects and freezes every substrate atom with zero extra
    bookkeeping needed. Distinct from stb-mldiffusion (vacancy migration
    in the bulk) -- this is migration of an ADSORBATE across a surface.

    Returns (barrier, dE, slot_a, slot_b) or None if fewer than 2 distinct
    sites are available for this adsorbate.
    """
    candidates = [r for r in scored if r[2] == adsorbate_name]
    chosen, seen_slots = [], []
    for r in candidates:
        if r[0] not in seen_slots:
            chosen.append(r)
            seen_slots.append(r[0])
        if len(chosen) == 2:
            break
    if len(chosen) < 2:
        print_dual(color_text(
            f"[WARNING] Not enough distinct sites for '{adsorbate_name}' to compute a "
            "diffusion barrier (need at least 2). Skipping.", 'yellow'), f_out)
        return None

    (slot_a, h_a, _, st_a, _e_a, atoms_a), (slot_b, h_b, _, st_b, _e_b, atoms_b) = chosen
    print_dual(f"Diffusion path    : site {slot_a} ({st_a}) -> site {slot_b} ({st_b}), "
               f"both for '{adsorbate_name}'", f_out)

    initial_pmg = wrap_into_cell(AseAtomsAdaptor.get_structure(atoms_a))
    final_pmg = wrap_into_cell(AseAtomsAdaptor.get_structure(atoms_b))

    neb_args = argparse.Namespace(
        n_images=args.diffusion_n_images, autosort_tol=0.5, idpp=True,
        idpp_fmax=0.1, idpp_steps=100, pre_relax_endpoints=False,
        freeze_substrate=True, freeze_threshold=1e-3, k=0.1,
        fmax=args.fmax, max_steps=args.max_steps,
    )
    barrier, dE, _ase_images, _pmg_images, converged = run_single_neb(
        initial_pmg, final_pmg, calc, neb_args, f_out, tag="diffusion")
    print_dual(f"Diffusion barrier : {color_text(f'{barrier:.4f}', 'bold')} eV "
               f"({'converged' if converged else 'hit step cap'})", f_out)
    return barrier, dE, slot_a, slot_b


def plot_site_ranking(series, out_path):
    """series: list of (label, scored, color). One bar group per series,
    same site/height/adsorbate/type combination side by side when
    comparing 2 models -- lets a direct visual check of whether fine-
    tuning changed which site ranks best.
    """
    label0, scored0, _ = series[0]
    labels = [f"{s[2]}\n#{s[0]}@{s[1]:.1f}A" for s in scored0]
    n_series = len(series)
    width = 0.8 / n_series
    fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(labels)), 5))
    for i, (label, scored, color) in enumerate(series):
        e_ads = [s[4] for s in scored]
        offsets = np.arange(len(labels)) + (i - (n_series - 1) / 2.0) * width
        ax.bar(offsets, e_ads, width=width, color=color, label=label)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("ML adsorption energy (eV)")
    ax.set_title("Adsorption-site ranking" if n_series == 1
                 else "Adsorption-site ranking: fine-tuned vs. foundation")
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
    if n_series > 1:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_height_curves(scored, out_path):
    """When multiple heights were screened (--height-sweep), plots E_ads
    vs. height per (site, adsorbate) -- the actual binding-energy curve
    shape, instead of just flatly ranking every (site, height) pair
    together. One line per distinct (slot, site_type, adsorbate).
    """
    groups = {}
    for slot, h, name, st, e_ads, _atoms in scored:
        groups.setdefault((slot, st, name), []).append((h, e_ads))

    fig, ax = plt.subplots(figsize=(8, 5))
    for (slot, st, name), points in sorted(groups.items()):
        points.sort()
        hs = [p[0] for p in points]
        es = [p[1] for p in points]
        ax.plot(hs, es, 'o-', label=f"site {slot} ({st}, {name})")
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
    ax.set_xlabel("Height above surface (Ang)")
    ax.set_ylabel("ML adsorption energy (eV)")
    ax.set_title("Binding-energy curve per site")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Standalone adsorption-site screening driven entirely by a MACE "
                     "potential -- the MACE-MP-0 foundation model, or a custom model "
                     "fine-tuned on your own SIESTA data via stb-mlffAnalysis "
                     "(--custom-model). Enumerates candidate ontop/bridge/hollow sites "
                     "(pymatgen's AdsorbateSiteFinder), relaxes each (substrate fixed) with "
                     "MACE, and reports a genuine ML adsorption energy per site (relative to "
                     "the clean slab + isolated, relaxed adsorbate) -- no SIESTA and no "
                     "calc.fdf/pseudopotentials needed. Also supports a foundation-model "
                     "comparison, random-orientation sampling for molecular adsorbates, a "
                     "binding-energy-vs-height curve (--height-sweep), and a surface-"
                     "diffusion-barrier estimate between the top sites (--diffusion-barrier). "
                     "For the real DFT-level pipeline, see stb-adsorb (and its own --ml-rank "
                     "preview mode, which this tool's machinery was extracted from).",
        epilog="Example usage:\n"
               "  stb-mladsorb --slab slab.fdf --adsorbate H\n"
               "  stb-mladsorb --slab slab.fdf --adsorbate O,H --site-type hollow --top-k 3\n"
               "  stb-mladsorb --slab slab.fdf --adsorbate H2O --n-orientations 8 "
               "--diffusion-barrier\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--slab", required=True, help="Slab/2D-material structure (.fdf), vacuum along one axis.")
    parser.add_argument("--adsorbate", required=True,
                        help="Adsorbate(s): an ASE G2 molecule name (e.g. 'H2O') or element "
                             "symbol, comma-separated for multiple (e.g. 'H,O').")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="If --custom-model is given, also screen every site with the raw "
                             "(non-fine-tuned) foundation model (--model) and overlay both on "
                             "the ranking plot plus a side-by-side table. On by default; pass "
                             "this to skip it.")
    parser.add_argument("--site-type", choices=["all", "ontop", "bridge", "hollow"], default="all",
                        help="Which candidate site type(s) to screen (default: all).")
    parser.add_argument("--height", type=float, default=2.0,
                        help="Adsorbate height above the surface, Ang (default: 2.0).")
    parser.add_argument("--height-sweep", type=float, nargs=3, default=None, metavar=("MIN", "MAX", "STEP"),
                        help="Sweep multiple heights instead of a single --height -- also "
                             "produces a binding-energy-vs-height curve per site.")
    parser.add_argument("--n-orientations", type=int, default=1,
                        help="For a multi-atom adsorbate, cheaply screen this many random "
                             "orientations per site (single-point energy, no relax) and keep "
                             "the best-looking one before the full relax (default: 1, i.e. "
                             "no sampling -- use the adsorbate's default G2 orientation). No "
                             "effect on single-atom adsorbates (no orientation to sample).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for --n-orientations reproducibility (default: unseeded).")
    parser.add_argument("--diffusion-barrier", action="store_true",
                        help="Also estimate the adsorbate's own surface-diffusion barrier "
                             "between its 2 most stable, distinct sites via a real climbing-"
                             "image NEB. Off by default.")
    parser.add_argument("--diffusion-adsorbate", default=None, metavar="NAME",
                        help="Which adsorbate's top sites to use for --diffusion-barrier when "
                             "multiple were given (default: the first one in --adsorbate).")
    parser.add_argument("--diffusion-n-images", type=int, default=5,
                        help="Total NEB images (including both endpoints) for --diffusion-barrier (default: 5).")
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry tolerance for site-finding/reduction (default: 0.01).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect the slab's vacuum axis (default: 10.0).")
    parser.add_argument("--vacuum-box", type=float, default=20.0,
                        help="Cubic box size (Ang) for the isolated-adsorbate reference (default: 20.0).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for all relaxations, eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max optimizer steps per relaxation (default: 200).")
    parser.add_argument("--top-k", type=int, default=None,
                        help="Only write out the N best-ranked sites as bare .fdf (default: all).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mladsorb_out",
                        help="Output directory for all files (default: mladsorb_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw ranking data (.dat). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mladsorb {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Adsorption Screening",
            "Adsorption-site ranking driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML ADSORPTION-SITE SCREENING:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLADSORB REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Slab file         : {args.slab}", f_out)
    print_dual(f"Adsorbate(s)      : {args.adsorbate}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"Site type         : {args.site_type}", f_out)
    compare = bool(args.custom_model) and not args.skip_foundation_comparison
    if compare:
        print_dual(f"Comparison        : also evaluating with foundation model ({args.model})", f_out)
    if args.n_orientations > 1:
        print_dual(f"Orientation sampling: {args.n_orientations} random orientation(s)/site", f_out)
    if args.diffusion_barrier:
        print_dual("Diffusion barrier : yes", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURE", f_out)

    if not os.path.isfile(args.slab):
        fail(f"Slab file '{args.slab}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")

    try:
        structure = structure_io.read_fdf(args.slab)
    except Exception as e:
        fail(f"Could not read '{args.slab}': {e}")

    try:
        structure, relabel_note = resolve_slab_orientation(structure, args.vacuum_gap)
    except ValueError as e:
        fail(str(e))
    if relabel_note:
        print_dual(color_text(f"[INFO] {relabel_note}", 'yellow'), f_out)
    pmg_structure = structure_io.to_pymatgen(structure)
    slab_species_meta = structure_io.species_dict(structure)
    n_substrate = len(pmg_structure)
    print_dual(f"[OK] Read {n_substrate} substrate atom(s)", f_out)

    # parse_adsorbates already prints its own clear error and exits(1) if an
    # adsorbate name doesn't resolve -- not wrapped in fail() here to avoid
    # a redundant second error message on top of that specific one.
    adsorbates = parse_adsorbates(args.adsorbate)
    diffusion_adsorbate = args.diffusion_adsorbate or adsorbates[0][0]
    if args.diffusion_barrier and diffusion_adsorbate not in [n for n, _ in adsorbates]:
        fail(f"--diffusion-adsorbate '{diffusion_adsorbate}' is not among --adsorbate {args.adsorbate}.")

    heights = build_sweep_values(*args.height_sweep) if args.height_sweep is not None else [args.height]
    site_types = ["ontop", "bridge", "hollow"] if args.site_type == "all" else [args.site_type]
    rng = np.random.default_rng(args.seed)

    plot_path = os.path.join(args.output_dir, "adsorption_sites.png")
    write_site_plot(pmg_structure, plot_path)
    print_dual(f"[OK] Saved {plot_path} (candidate site layout)", f_out)

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)

    print_section("[2] REFERENCE ENERGIES & SITE SCREENING", f_out)
    scored, e_slab, e_isolated = screen_sites(
        pmg_structure, n_substrate, adsorbates, heights, site_types, args, calc, rng, f_out,
        tag="fine-tuned" if args.custom_model else args.model)
    if scored is None:
        fail(f"No {'/'.join(site_types)} sites found.")
    print_dual(f"[OK] {len(scored)} (site, height, adsorbate) combination(s) screened", f_out)

    scored_found = None
    if compare:
        print_dual(f"\n[INFO] Also screening every site with the raw foundation model "
                   f"({args.model}) for comparison...", f_out)
        found_calc = mace_relax.get_calculator(model=args.model, device=args.device)
        scored_found, _e_slab_f, _e_iso_f = screen_sites(
            pmg_structure, n_substrate, adsorbates, heights, site_types, args, found_calc,
            rng, f_out, tag=f"foundation ({args.model})")

    print_section("[3] RANKING", f_out)
    print_dual(f"{'Rank':<6}{'Site':<6}{'Height':<9}{'Adsorbate':<12}{'Type':<9}{'E_ads (eV)':<12}", f_out)
    for rank, (slot, h, name, st, e_ads, _atoms) in enumerate(scored, start=1):
        print_dual(f"{rank:<6}{slot:<6}{h:<9.2f}{name:<12}{st:<9}{e_ads:<12.4f}", f_out)
    print_dual(color_text(
        "[NOTE] An ML-level estimate from a fast surrogate potential, not a DFT-level "
        "adsorption energy -- use it to prioritize which site(s) to relax with SIESTA "
        "(stb-adsorb/stb-adsorbAnalysis), not as a final answer.", 'yellow'), f_out)

    if scored_found is not None:
        print_section("[3b] FOUNDATION MODEL COMPARISON", f_out)
        best_ft = scored[0]
        best_found_slot = scored_found[0][0]
        same_best = best_ft[0] == best_found_slot
        print_dual(f"Best site (fine-tuned) : site {best_ft[0]} ({best_ft[3]}), "
                   f"E_ads = {best_ft[4]:.4f} eV", f_out)
        print_dual(f"Best site (foundation) : site {scored_found[0][0]} ({scored_found[0][3]}), "
                   f"E_ads = {scored_found[0][4]:.4f} eV", f_out)
        print_dual(f"Same best site         : {'yes' if same_best else 'no'}", f_out)
        if not same_best:
            print_dual(color_text(
                "[INFO] Fine-tuning changed which site ranks best -- worth double-checking "
                "against real DFT before trusting either ranking.", 'cyan'), f_out)

    if len(heights) > 1:
        height_plot = os.path.join(args.output_dir, "height_curve.png")
        plot_height_curves(scored, height_plot)
        print_dual(f"[OK] Saved {height_plot}", f_out)

    if args.diffusion_barrier:
        print_section("[3c] SURFACE DIFFUSION BARRIER", f_out)
        run_diffusion_barrier(scored, diffusion_adsorbate, args, calc, f_out)

    print_section("[4] SUMMARY & FILES", f_out)
    series = [(("Fine-tuned" if args.custom_model else f"MACE-MP-0 ({args.model})"), scored,
               _SERIES_COLORS[0])]
    if scored_found is not None:
        series.append((f"Foundation ({args.model})", scored_found, _SERIES_COLORS[1]))
    ranking_plot = os.path.join(args.output_dir, "site_ranking.png")
    plot_site_ranking(series, ranking_plot)
    print_dual(f"[OK] Saved {ranking_plot}", f_out)
    if args.save_data:
        data_path = os.path.join(args.output_dir, "site_ranking.dat")
        with open(data_path, "w") as f:
            f.write("# rank site height adsorbate type E_ads(eV)\n")
            for rank, (slot, h, name, st, e_ads, _atoms) in enumerate(scored, start=1):
                f.write(f"{rank} {slot} {h:.4f} {name} {st} {e_ads:.6f}\n")
        print_dual(f"[OK] Saved {data_path}", f_out)

    write_list = scored[:args.top_k] if args.top_k is not None else scored
    for rank, (slot, h, name, st, e_ads, ase_atoms) in enumerate(write_list, start=1):
        pmg_relaxed = AseAtomsAdaptor.get_structure(ase_atoms)
        fdf_structure = structure_io.from_pymatgen(pmg_relaxed, species_meta=slab_species_meta,
                                                    coord_format="fractional")
        site_path = os.path.join(args.output_dir, f"rank{rank:02d}_site{slot}_{st}_{name}.fdf")
        structure_io.write_fdf(fdf_structure, site_path)
    print_dual(f"[OK] Saved {len(write_list)} ranked site structure(s) to {args.output_dir}/", f_out)

    print_dual("Status            : OK", f_out)
    print_dual(f"Best site         : rank 1, site {scored[0][0]} ({scored[0][3]}), "
               f"E_ads = {scored[0][4]:.4f} eV", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML adsorption-site screening complete -- a fast heuristic, not a "
                      "substitute for a DFT-level adsorption energy.\n", 'bold'))


if __name__ == "__main__":
    main()
