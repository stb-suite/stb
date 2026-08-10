#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.calc_directives import force_single_point
from stb.core.heterostructure import find_zsl_match, build_stacked_structure
from stb.core.deps import require_mace

REPORT_FILE = "stackingfault_setup.txt"


def read_layer(path):
    """Reads one monolayer .fdf into a pymatgen Structure, exiting with a
    clear error on any parse failure -- same "fail loud, don't propagate a
    cryptic traceback" convention as adsorb.py/neb.py's own file reading.
    """
    try:
        return structure_io.to_pymatgen(structure_io.read_fdf(path))
    except Exception as e:
        print(color_text(f"[ERROR] Could not read '{path}': {e}", 'red'))
        sys.exit(1)


def min_interlayer_distance(pmg_structure, n_layer1):
    """Minimum periodic distance between any layer-1 atom (index <
    n_layer1) and any layer-2 atom (index >= n_layer1) -- same style as
    adsorb.py's min_adsorbate_slab_distance, generalized from
    adsorbate-vs-slab to layer-vs-layer. Purely informational here (unlike
    adsorb.py's overlap WARNING): a close contact at some grid point is
    often exactly the physics a stacking-fault sweep is meant to sample
    (an eclipsed, high-energy registry), not a mistake to flag.
    """
    n_total = len(pmg_structure)
    if n_layer1 == 0 or n_layer1 >= n_total:
        return None
    dm = pmg_structure.distance_matrix
    return float(dm[:n_layer1, n_layer1:].min())


def optimize_gap(layer1_pmg, layer2_pmg, t_mat1, t_mat2, shift_x, shift_y, calc,
                  nominal_gap, window=1.0, n_scan=8, target_vacuum=None, strain_mode='top'):
    """Finds the MACE-MP-0-optimal RIGID interlayer gap for one grid point
    (only the gap distance is varied -- both layers stay internally rigid,
    matching the "rigid-shift" protocol's own spirit, just no longer with
    a single gap value forced on every registry). A coarse scan over
    `n_scan` points in [nominal_gap - window, nominal_gap + window] first,
    then a bounded local polish around the scan's minimum -- same global-
    search-then-polish strategy as neb_analysis.py's fit_spline_barrier,
    for the same reason: a single bounded scalar search can miss the true
    minimum if it starts in a purely repulsive region far from the actual
    well (e.g. --gap set too small for this specific registry).

    Returns (optimal_gap, energy_at_optimal_gap) -- the energy is returned
    too so a caller also wanting an --ml-preview value doesn't need a
    second MACE evaluation at the same (now-known) geometry.
    """
    from scipy.optimize import minimize_scalar

    def energy_at_gap(gap):
        hetero, _n1, _strain = build_stacked_structure(
            layer1_pmg, layer2_pmg, t_mat1, t_mat2, shift_x, shift_y, gap,
            target_vacuum=target_vacuum, strain_mode=strain_mode)
        atoms = AseAtomsAdaptor.get_atoms(hetero)
        atoms.calc = calc
        return atoms.get_potential_energy()

    gap_min = max(nominal_gap - window, 0.5)
    gaps = list(np.linspace(gap_min, nominal_gap + window, n_scan))
    energies = [energy_at_gap(g) for g in gaps]
    i_min = int(np.argmin(energies))
    lo = gaps[max(i_min - 1, 0)]
    hi = gaps[min(i_min + 1, len(gaps) - 1)]
    result = minimize_scalar(energy_at_gap, bounds=(lo, hi), method='bounded')
    return float(result.x), float(result.fun)


def write_ml_preview_plot(ml_grid, grid_n, out_path):
    """2D heatmap preview of the MACE-MP-0 single-point energy at every
    grid point (E - E_min, eV) -- a fast sanity check of the gamma-
    surface's shape before committing to the full N**2 DFT single-points.
    Plain matplotlib imshow, NOT the gnuplot pm3d convention
    stb-stackingfaultAnalysis uses for the real DFT-level map (that one
    reuses core/grid_export.py, shared with density.py's charge-density
    plots) -- self-contained here instead, same "prep tool writes its own
    quick PNG directly" pattern as neb.py's write_ml_preview_plot /
    adsorb.py's write_site_plot, since this is the only consumer.
    """
    energies = np.full((grid_n, grid_n), np.nan)
    for i, j, energy in ml_grid:
        energies[i, j] = energy
    e_rel = energies - np.nanmin(energies)

    fig, ax = plt.subplots()
    im = ax.imshow(e_rel.T, origin='lower', extent=(0, 1, 0, 1), cmap='inferno', aspect='auto')
    ax.set_xlabel("shift_x (fractional)")
    ax.set_ylabel("shift_y (fractional)")
    ax.set_title("ML (MACE-MP-0) stacking-fault energy preview")
    fig.colorbar(im, ax=ax, label="E - E_min (eV)")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_grid_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path):
    """Writes structure.fdf + calc.fdf + copied pseudopotentials for one
    shift_II_JJ/ grid-point folder -- same shape as adsorb.py's
    write_reference_folder / neb.py's write_image_folder. `species_meta`
    only needs to cover layer1's declared species -- any symbol unique to
    layer2 gets a fresh id automatically via structure_io.from_pymatgen's
    ensure_species_id, same as adsorb.py passing the slab's species_meta
    alone when the adsorbate introduces a new element.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta,
                                                coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = {site.specie.symbol for site in pmg_structure}
    for sym in sorted(symbols):
        copy_pseudo(pp_path, sym, out_dir)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Prepares SIESTA folders for a 2D stacking-fault (generalized "
        "stacking fault energy / gamma-surface) study: rigidly slides one layer of a bilayer "
        "across a 2D grid of lateral offsets and writes one single-point SIESTA folder per grid "
        "point.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -n 7\n"
               "  %(prog)s -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf -n 9 -g 3.4\n"
               "\n"
               "Pass the SAME file for -l1/-l2 for the canonical use case (a material sliding "
               "against itself, e.g. graphite ABA vs. ABC stacking) -- already the precedented "
               "'identical layers' case for stb-2Dstacking's own ZSL match. Different files "
               "study an interlayer/heterostructure sliding energy landscape instead.\n"
    )

    parser.add_argument("-l1", "--layer1", type=str, required=True,
                         help="Bottom monolayer .fdf (stays fixed).")
    parser.add_argument("-l2", "--layer2", type=str, required=True,
                         help="Top monolayer .fdf (rigidly slid across the grid). Pass the same "
                              "file as --layer1 for a material sliding against itself.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template (kgrid, basis, XC, %%include structure.fdf, "
                              "etc.) -- copied into every shift_II_JJ/ folder with MD.TypeOfRun/"
                              "MD.Steps forced to a single-point evaluation.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (optional): a bundled bank or a folder path.")

    parser.add_argument("-a", "--max_area", type=float, default=150.0,
                         help="Maximum ZSL commensurate-supercell area, Ang^2 (default: 150.0, "
                              "same as stb-2Dstacking).")
    parser.add_argument("-s", "--max_strain", type=float, default=0.05,
                         help="Maximum allowed ZSL match strain fraction (default: 0.05).")
    parser.add_argument("-id", "--match_id", type=int, default=0,
                         help="Which ZSL match to use, 0-based, best (lowest strain) first "
                              "(default: 0).")
    parser.add_argument("-g", "--gap", type=float, default=3.2,
                         help="Interlayer gap in Ang, FIXED across the whole grid (default: 3.2, "
                              "same as stb-2Dstacking) -- a rigid-shift protocol, standard in the "
                              "literature; per-point interlayer-distance relaxation is not done "
                              "here (some high-energy registries may end up with a non-relaxed "
                              "distance -- see the [INFO] closest-contact line in the report).")
    parser.add_argument("--vacuum", type=float, default=None,
                         help="Target vacuum space in Ang. Inherits --layer1's by default.")
    parser.add_argument("-sm", "--strain_mode", choices=["top", "bottom", "sym"], default="top",
                         help="Strain distribution mode for the ZSL match (default: top).")
    parser.add_argument("-t", "--twist", type=float, default=0.0,
                         help="Twist angle of layer 2 in degrees (default: 0.0), FIXED for the "
                              "whole grid -- a different twist is a different physical system, "
                              "not swept here.")

    parser.add_argument("-n", "--grid-n", type=int, default=7,
                         help="Grid resolution: an N x N grid of lateral shifts covering [0, 1) "
                              "in each direction (default: 7).")

    parser.add_argument("--ml-prerelax-layers", action="store_true",
                         help="Relax each monolayer's positions (cell fixed) with MACE-MP-0 "
                              "before stacking -- a cheap safety net when you're not fully "
                              "certain --layer1/--layer2 are already at equilibrium. Needs the "
                              "optional 'ml' extra.")
    parser.add_argument("--ml-relax-gap", action="store_true",
                         help="Find the MACE-MP-0-optimal RIGID interlayer gap at EACH grid "
                              "point instead of using the single fixed --gap everywhere -- a "
                              "real physics improvement over the rigid-shift default (some "
                              "registries are more repulsive/attractive than others at a given "
                              "gap). Needs the optional 'ml' extra; slower (a small gap scan per "
                              "grid point, not just one MACE evaluation).")
    parser.add_argument("--ml-relax-gap-window", type=float, default=1.0,
                         help="With --ml-relax-gap: scan +/- this many Ang around --gap "
                              "(default: 1.0).")
    parser.add_argument("--ml-relax-gap-steps", type=int, default=8,
                         help="With --ml-relax-gap: number of coarse scan points before the "
                              "local polish (default: 8).")
    parser.add_argument("--ml-preview", action="store_true",
                         help="Evaluate the whole grid's single-point energy on MACE-MP-0 (fast, "
                              "no SIESTA) and write stackingfault_ml_preview.png -- a quick sanity "
                              "check of the gamma-surface's shape before committing to the full "
                              "N x N DFT single-points. Written alongside the normal SIESTA "
                              "folders, not instead of them. Needs the optional 'ml' extra.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--ml-fmax", type=float, default=0.05,
                         help="Force convergence for --ml-prerelax-layers (default: 0.05).")

    parser.add_argument("-O", "--output-dir", type=str, default=".",
                         help="Root directory (default: current directory) for every "
                              "shift_II_JJ/.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-stackingfault {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.grid_n < 2:
        parser.error("-n/--grid-n must be >= 2 (need at least 2 points per axis for a sweep).")

    print("\n" + color_text("Prepare a 2D stacking-fault study:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.layer1):
        print(color_text(f"[ERROR] Layer 1 file '{args.layer1}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.layer2):
        print(color_text(f"[ERROR] Layer 2 file '{args.layer2}' not found.", 'red'))
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

    layer1_structure = structure_io.read_fdf(args.layer1)
    species_meta = structure_io.species_dict(layer1_structure)
    layer1_pmg = read_layer(args.layer1)
    layer2_pmg = read_layer(args.layer2)

    with open(args.calc) as f:
        calc_text = f.read()
    single_point_calc_text = force_single_point(calc_text)

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, 'w') as f_out:
        print_dual(f"\n{color_text('[0] RUN METADATA', 'bold')}", f_out)
        print_dual(f"  Layer 1         : {args.layer1}", f_out)
        print_dual(f"  Layer 2         : {args.layer2}", f_out)
        print_dual(f"  Calc template   : {args.calc}", f_out)
        print_dual(f"  Pseudo dir      : {args.pseudo_dir or '(none)'}", f_out)
        print_dual(f"  Output dir      : {output_root}", f_out)
        print_dual(f"  Grid            : {args.grid_n} x {args.grid_n}", f_out)
        print_dual(f"  Gap (fixed)     : {args.gap} Ang"
                    + (" (overridden per-point by --ml-relax-gap below)" if args.ml_relax_gap else ""),
                    f_out)
        print_dual(f"  Twist (fixed)   : {args.twist} deg", f_out)
        print_dual(f"  max_area/max_strain/match_id: {args.max_area} / {args.max_strain} / {args.match_id}", f_out)
        print_dual(f"  ML pre-relax layers : {'yes' if args.ml_prerelax_layers else 'no'}", f_out)
        print_dual(f"  ML relax gap        : {'yes' if args.ml_relax_gap else 'no'}", f_out)
        print_dual(f"  ML preview          : {'yes' if args.ml_preview else 'no'}", f_out)

        calc_mace = None
        if args.ml_prerelax_layers or args.ml_relax_gap or args.ml_preview:
            require_mace()
            from stb.core import mace_relax
            calc_mace = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)

        if args.ml_prerelax_layers:
            print_dual(f"\n{color_text('ML pre-relax:', 'cyan')} relaxing each monolayer "
                        "(positions only) with MACE-MP-0 ...", f_out)
            for label, pmg in (("layer 1", layer1_pmg), ("layer 2", layer2_pmg)):
                ase_atoms = AseAtomsAdaptor.get_atoms(pmg)
                converged, steps = mace_relax.relax(ase_atoms, calc_mace, fmax=args.ml_fmax, max_steps=200)
                relaxed_pmg = AseAtomsAdaptor.get_structure(ase_atoms)
                if label == "layer 1":
                    layer1_pmg = relaxed_pmg
                else:
                    layer2_pmg = relaxed_pmg
                print_dual(f"  {'Converged' if converged else 'Hit step cap, not fully converged'} "
                            f"({label}) after {steps} step(s).", f_out)

        print_dual(f"\n{color_text('[1] ZSL MATCH', 'bold')}", f_out)
        t_mat1, t_mat2, best_match_data = find_zsl_match(
            layer1_pmg, layer2_pmg, max_area=args.max_area, max_strain=args.max_strain,
            match_id=args.match_id, interactive=False, twist_angle=args.twist)
        print_dual(f"  Selected match ID {args.match_id}: area {best_match_data['area']:.2f} Ang^2, "
                    f"strain {best_match_data['strain']:.2f}%, angular strain "
                    f"{best_match_data['angle_strain']:.2f} deg.", f_out)
        n_layer1_supercell = len(layer1_pmg) * round(np.linalg.det(t_mat1))
        n_layer2_supercell = len(layer2_pmg) * round(np.linalg.det(t_mat2))
        n_total = n_layer1_supercell + n_layer2_supercell
        n_grid_points = args.grid_n * args.grid_n
        print_dual(f"  Each grid point: {n_total} atoms ({n_layer1_supercell} layer 1 + "
                    f"{n_layer2_supercell} layer 2). {n_grid_points} grid point(s) total -- "
                    f"{n_grid_points} independent single-point SIESTA runs at {n_total} atoms each.",
                    f_out)

        shifts = list(np.linspace(0.0, 1.0, args.grid_n, endpoint=False))

        print_dual(f"\n{color_text('[2] GRID FOLDERS', 'bold')}", f_out)
        if args.ml_relax_gap:
            print_dual(f"  {color_text('ML relax gap:', 'cyan')} scanning +/- {args.ml_relax_gap_window} "
                        f"Ang around {args.gap} Ang at every grid point with MACE-MP-0 ...", f_out)
        report_rows = []  # (label, i, j, shift_x, shift_y, gap_used, dir)
        closest_contact = None
        closest_label = None
        ml_grid = []  # (i, j, energy) -- only populated if --ml-preview
        gaps_used = []
        for i, shift_x in enumerate(shifts):
            for j, shift_y in enumerate(shifts):
                label = f"shift_{i:02d}_{j:02d}"

                ml_energy = None
                if args.ml_relax_gap:
                    gap_used, ml_energy = optimize_gap(
                        layer1_pmg, layer2_pmg, t_mat1, t_mat2, shift_x, shift_y, calc_mace,
                        args.gap, window=args.ml_relax_gap_window, n_scan=args.ml_relax_gap_steps,
                        target_vacuum=args.vacuum, strain_mode=args.strain_mode)
                else:
                    gap_used = args.gap
                gaps_used.append(gap_used)

                hetero, n_layer1_atoms, _max_strain_val = build_stacked_structure(
                    layer1_pmg, layer2_pmg, t_mat1, t_mat2, shift_x, shift_y, gap_used,
                    target_vacuum=args.vacuum, strain_mode=args.strain_mode)

                if args.ml_preview:
                    if ml_energy is None:
                        atoms = AseAtomsAdaptor.get_atoms(hetero)
                        atoms.calc = calc_mace
                        ml_energy = atoms.get_potential_energy()
                    ml_grid.append((i, j, ml_energy))

                contact = min_interlayer_distance(hetero, n_layer1_atoms)
                if contact is not None and (closest_contact is None or contact < closest_contact):
                    closest_contact = contact
                    closest_label = label

                grid_dir = os.path.join(output_root, label)
                write_grid_folder(grid_dir, hetero, single_point_calc_text, species_meta,
                                   args.pseudo_dir)
                gap_note = f", gap: {gap_used:.3f} Ang (ML)" if args.ml_relax_gap else ""
                print_dual(f"  {color_text('[OK]', 'green')} {grid_dir} "
                            f"(shift: {shift_x:.4f}, {shift_y:.4f}{gap_note})", f_out)
                report_rows.append((label, i, j, shift_x, shift_y, gap_used, grid_dir))

        print_dual(f"\n{color_text('[3] SUMMARY', 'bold')}", f_out)
        print_dual(f"  {len(report_rows)} grid folder(s) written under '{output_root}'.", f_out)
        if closest_contact is not None:
            print_dual(f"  [INFO] Closest interlayer contact anywhere in the grid: "
                        f"{closest_contact:.3f} Ang at {closest_label} -- expected to be small at "
                        "eclipsed/high-energy registries, not necessarily a problem.", f_out)
        if args.ml_relax_gap:
            print_dual(f"  [INFO] ML-optimized gap ranged {min(gaps_used):.3f} - "
                        f"{max(gaps_used):.3f} Ang across the grid (nominal --gap was {args.gap}).",
                        f_out)
        if args.ml_preview and ml_grid:
            preview_path = os.path.join(output_root, "stackingfault_ml_preview.png")
            write_ml_preview_plot(ml_grid, args.grid_n, preview_path)
            ml_min = min(ml_grid, key=lambda r: r[2])
            ml_max = max(ml_grid, key=lambda r: r[2])
            print_dual(f"  {color_text('[Saved]', 'cyan')} {preview_path}", f_out)
            print_dual(f"  ML preview: predicted equilibrium at shift_{ml_min[0]:02d}_{ml_min[1]:02d}, "
                        f"predicted corrugation {ml_max[2] - ml_min[2]:.4f} eV (MACE-MP-0 -- NOT a "
                        "DFT-level number, use stb-stackingfaultAnalysis for the real one).", f_out)
        print_dual("  Run SIESTA (single-point: MD.Steps forced to 0) in every "
                    "shift_II_JJ/ folder, then use stb-stackingfaultAnalysis.", f_out)

        f_out.write("\n# GRID_TABLE -- parsed by stb-stackingfaultAnalysis, do not reorder the "
                     "first 5 columns\n")
        f_out.write(f"# {'label':<16}{'i':<5}{'j':<5}{'shift_x':<12}{'shift_y':<12}{'gap':<10}{'dir'}\n")
        for label, i, j, shift_x, shift_y, gap_used, grid_dir in report_rows:
            f_out.write(f"{label:<18}{i:<5}{j:<5}{shift_x:<12.6f}{shift_y:<12.6f}{gap_used:<10.4f}{grid_dir}\n")

    print(f"\n{color_text('Success:', 'green')} {len(report_rows)} grid folder(s) written under "
          f"'{output_root}'.")
    print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
