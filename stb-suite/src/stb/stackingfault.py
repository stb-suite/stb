#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io, adsorption_sites
from stb.core.cli import color_text, show_intro, print_dual, print_section, capture_library_noise
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.heterostructure import find_zsl_match, build_stacked_structure
from stb.core.deps import require_mace

# RUN_SUBDIR mirrors stb-neb's neb_run/ (and stb-adsorb's sites/) convention:
# every artifact of one run -- grid folders, the manifest, the optional
# narrative report, the optional ML preview -- lives together under one
# self-contained subfolder of --output-dir, not scattered directly into it.
RUN_SUBDIR = "sf_run"
# POSITIONS_SUBDIR nests every shift_II_JJ/ grid folder one level deeper,
# under <run_root>/positions/ -- a sibling of stb-stackingfaultBsse's own
# <run_root>/bsse/ tree, so the run root only ever has 2 "kind of content"
# subfolders (positions/, bsse/) plus its own manifest/report/preview,
# instead of dozens of shift_II_JJ/ folders mixed in at the same level as
# those run-level files -- same "one folder per kind of artifact" tidiness
# stb-adsorb already gets from 'sites/' + 'bsse/' both living under its own
# --output-dir root.
POSITIONS_SUBDIR = "positions"
# REPORT_FILE is the optional human narrative ([0]-[3]), written only with
# --save-report. MANIFEST_FILE is the always-on machine-readable handoff to
# stb-stackingfaultAnalysis (grid shape + per-point shift/gap), independent
# of whether a narrative report was requested -- same split as stb-adsorb's
# always-on site table vs. its own --save-report-gated narrative, and
# stb-hubbardu's always-on run_manifest.json.
REPORT_FILE = "stackingfault_setup.txt"
MANIFEST_FILE = "sf_manifest.json"


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


def effective_gap(pmg_structure, n_layer1):
    """Cartesian z-gap (Ang) actually realized between layer 1's topmost
    atom and layer 2's bottommost atom -- the inverse of the quantity
    build_stacked_structure uses to place the nominal --gap in the first
    place (z_max_l1, z_min_l2 in core/heterostructure.py). Used to report
    the REAL gap reached after a --mode 2 (MACE) z-relaxation, since the
    single scalar --gap no longer describes the geometry once individual
    atoms have moved independently.
    """
    n_total = len(pmg_structure)
    if n_layer1 == 0 or n_layer1 >= n_total:
        return None
    z_max_l1 = max(pmg_structure[i].coords[2] for i in range(n_layer1))
    z_min_l2 = min(pmg_structure[i].coords[2] for i in range(n_layer1, n_total))
    return float(z_min_l2 - z_max_l1)


def relax_gap_ml(hetero, calc, fmax=0.05, max_steps=200):
    """--mode 2: relaxes ONLY the z-coordinate of every atom (x, y frozen
    via ase.constraints.FixedLine(direction=[0, 0, 1])) with MACE-MP-0 --
    same physical definition as --mode 1's real-SIESTA restricted
    relaxation (build_z_relax_block's Geometry.Constraints block), just a
    cheaper engine, used as a pre-step before writing a single-point
    SIESTA folder at the resulting geometry. core.mace_relax.relax()
    never touches atoms.constraints, so setting the constraint before
    calling it is enough -- no changes needed there.

    Returns (relaxed_pmg_structure, energy) -- the energy is returned too
    so --ml-preview doesn't need a second MACE evaluation at the same
    (now-relaxed) geometry.
    """
    from ase.constraints import FixedLine
    from stb.core import mace_relax

    atoms = AseAtomsAdaptor.get_atoms(hetero)
    atoms.set_constraint(FixedLine(list(range(len(atoms))), direction=[0, 0, 1]))
    mace_relax.relax(atoms, calc, fmax=fmax, max_steps=max_steps)
    energy = float(atoms.get_potential_energy())
    relaxed_pmg = AseAtomsAdaptor.get_structure(atoms)
    return relaxed_pmg, energy


def build_z_relax_block(n_atoms, md_steps):
    """--mode 1's config_extra.fdf: a real, restricted SIESTA CG
    relaxation where every atom's x, y is frozen (both layers, including
    layer 1 -- the nominal reference -- as a belt-and-suspenders guard
    against numerical drift) and only z is free. Two SEPARATE
    'position ... vx vy vz' lines (one for x, one for y), NOT one combined
    line with two 1.0's -- confirmed via the SIESTA mailing list
    ("Geometry constraints" thread) and a manual citation, both agreeing
    that a single combined line does not reliably constrain both axes.
    NOT verified against a real siesta run in this environment (no SIESTA
    binary available here) -- check the .out/.FA of ONE grid point before
    trusting a full --mode 1 production sweep (x, y displacement should be
    ~0, z should have moved).
    """
    return (
        "# Auto-generated -- restricted relaxation: only z is free (x, y frozen for\n"
        "# EVERY atom) -- see build_z_relax_block's docstring for the source/caveat\n"
        "# on this Geometry.Constraints syntax.\n"
        "MD.TypeOfRun          CG\n"
        f"MD.Steps              {md_steps}\n"
        "MD.VariableCell        false\n"
        "%block Geometry.Constraints\n"
        f"  position from 1 to {n_atoms} 1.0 0.0 0.0\n"
        f"  position from 1 to {n_atoms} 0.0 1.0 0.0\n"
        "%endblock Geometry.Constraints\n"
    )


def build_scan_points(scan, grid_nx, grid_ny, scan_n):
    """Builds the list of (label, i, j, shift_x, shift_y) points to sample
    for every --scan mode. 'surface' is today's full Cartesian-product grid
    (Nx x Ny, unchanged). The three 1D scans (x/y/xy) are special cases of
    the SAME (i, j) / 'shift_II_JJ' indexing scheme -- 'x' holds j fixed at
    0 (shift_y always 0.0), 'y' holds i fixed at 0 (shift_x always 0.0),
    and 'xy' walks the diagonal i == j (shift_x == shift_y) -- so
    stb-stackingfaultAnalysis's existing 'shift_(\\d+)_(\\d+)$' folder
    parser and manifest row schema need no changes at all to support any
    of the three; only which (i, j) pairs get generated changes.
    """
    if scan == "surface":
        shifts_x = list(np.linspace(0.0, 1.0, grid_nx, endpoint=False))
        shifts_y = list(np.linspace(0.0, 1.0, grid_ny, endpoint=False))
        return [(f"shift_{i:02d}_{j:02d}", i, j, shift_x, shift_y)
                for i, shift_x in enumerate(shifts_x)
                for j, shift_y in enumerate(shifts_y)]
    shifts = list(np.linspace(0.0, 1.0, scan_n, endpoint=False))
    if scan == "x":
        return [(f"shift_{i:02d}_00", i, 0, shift, 0.0) for i, shift in enumerate(shifts)]
    if scan == "y":
        return [(f"shift_00_{j:02d}", 0, j, 0.0, shift) for j, shift in enumerate(shifts)]
    # scan == "xy": diagonal, shift_x == shift_y
    return [(f"shift_{i:02d}_{i:02d}", i, i, shift, shift) for i, shift in enumerate(shifts)]


def write_ml_preview_line_plot(ml_grid, scan, out_path):
    """1D line-plot preview of the MACE-MP-0 single-point energy along a
    --scan x/y/xy path (E - E_min, eV vs. fractional distance along the
    path) -- the --scan {x,y,xy} analog of write_ml_preview_plot's 2D
    heatmap for --scan surface. Fractional (not Cartesian) x-axis, same
    convention as write_ml_preview_plot's own fractional shift_x/shift_y
    axes -- Stage 1's preview always stays fractional; only Stage 2's real
    DFT-level plot (stb-stackingfaultAnalysis) converts to Cartesian Ang.
    ml_grid entries are (i, j, energy); by build_scan_points's construction
    the swept index (i for x/xy, j for y) already runs 0..n-1 in path order.
    """
    running_idx = 0 if scan in ("x", "xy") else 1
    points = sorted(ml_grid, key=lambda r: r[running_idx])
    n = len(points)
    frac = [r[running_idx] / n for r in points]
    e_rel = np.array([r[2] for r in points])
    e_rel = e_rel - e_rel.min()

    fig, ax = plt.subplots()
    ax.plot(frac, e_rel, marker='o', color='#cc5522')
    xlabel = "shift_x = shift_y (fractional)" if scan == "xy" else f"shift_{scan} (fractional)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("E - E_min (eV)")
    ax.set_title(f"ML (MACE-MP-0) stacking-fault energy preview ({scan} scan)")
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_ml_preview_plot(ml_grid, grid_nx, grid_ny, out_path):
    """2D heatmap preview of the MACE-MP-0 single-point energy at every
    grid point (E - E_min, eV) -- a fast sanity check of the gamma-
    surface's shape before committing to the full Nx*Ny DFT single-points.
    Plain matplotlib imshow, NOT the gnuplot pm3d convention
    stb-stackingfaultAnalysis uses for the real DFT-level map (that one
    reuses core/grid_export.py, shared with density.py's charge-density
    plots) -- self-contained here instead, same "prep tool writes its own
    quick PNG directly" pattern as neb.py's write_ml_preview_plot /
    adsorb.py's write_site_plot, since this is the only consumer. Supports
    an asymmetric grid_nx != grid_ny (independent resolution per shift
    axis) -- the (i, j) grid indices from the main loop map directly onto
    this array's shape regardless.
    """
    energies = np.full((grid_nx, grid_ny), np.nan)
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


def write_grid_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path,
                       d3=False, mode=3, relax_z_steps=200):
    """Writes structure.fdf + config_extra.fdf + calc.fdf + copied
    pseudopotentials for one shift_II_JJ/ grid-point folder -- same shape
    as adsorb.py's write_reference_folder / neb.py's write_image_folder.

    config_extra.fdf's core content depends on `mode`: modes 2 and 3 force
    a plain single-point evaluation (core.adsorption_sites.
    SINGLE_POINT_BLOCK -- each grid point is one independent, uncoupled
    sample; under mode 2 `pmg_structure` has already been z-relaxed by
    MACE before this call, so SIESTA only needs to confirm/evaluate it
    once) while mode 1 instead writes a real restricted SIESTA relaxation
    (build_z_relax_block -- SIESTA itself relaxes z here, so MD.Steps is
    NOT forced to 0). In every mode, `d3=True` additionally appends the
    Grimme DFT-D3 dispersion correction (VDW_CORRECTION_BLOCK) --
    physically relevant here since the interlayer binding being scanned
    across the grid is a van-der-Waals-dominated interaction that plain
    GGA misses. Goes through config_extra.fdf rather than editing --calc's
    own text in place (structure_io.prepend_include prepends '%include
    config_extra.fdf' at the very top of the written calc.fdf) -- same
    override mechanism, and the same reason neb.py's write_image_folder
    moved off of core.calc_directives.force_single_point, as every other
    config_extra.fdf consumer in the suite. `species_meta` only needs to
    cover layer1's declared species -- any symbol unique to layer2 gets a
    fresh id automatically via structure_io.from_pymatgen's
    ensure_species_id, same as adsorb.py passing the slab's species_meta
    alone when the adsorbate introduces a new element.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta,
                                                coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, adsorption_sites.CONFIG_EXTRA_FILE), "w") as f:
        if mode == 1:
            f.write(build_z_relax_block(len(pmg_structure), relax_z_steps))
        else:
            f.write(adsorption_sites.SINGLE_POINT_BLOCK)
        if d3:
            f.write(adsorption_sites.VDW_CORRECTION_BLOCK)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(structure_io.prepend_include(calc_text, adsorption_sites.CONFIG_EXTRA_FILE))
    symbols = {site.specie.symbol for site in pmg_structure}
    for sym in sorted(symbols):
        copy_pseudo(pp_path, sym, out_dir)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Prepares SIESTA folders for a 2D stacking-fault (generalized "
        "stacking fault energy / gamma-surface) study: rigidly slides one layer of a bilayer "
        "across a 2D grid of lateral offsets and writes one SIESTA folder per grid point -- a "
        "plain single-point by default, or a real/ML-assisted z-only relaxation, see --mode.",
        'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 7 -ny 7\n"
               "  %(prog)s -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf -nx 9 -ny 5 -g 3.4\n"
               "  %(prog)s -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 15\n"
               "\n"
               "Pass the SAME file for -l1/-l2 for the canonical use case (a material sliding "
               "against itself, e.g. graphite ABA vs. ABC stacking) -- already the precedented "
               "'identical layers' case for stb-2Dstacking's own ZSL match. Different files "
               "study an interlayer/heterostructure sliding energy landscape instead.\n"
               "\n"
               "-nx/-ny may differ (an asymmetric grid) -- useful for an anisotropic 2D lattice "
               "where one in-plane direction needs finer sampling than the other (--scan surface "
               "only).\n"
               "\n"
               "--scan x/y/xy trade the full 2D gamma-surface for a cheaper 1D line through it "
               "(along shift_x, shift_y, or the shift_x=shift_y diagonal) -- useful for a quick "
               "high-symmetry-line profile before committing to the full surface.\n"
    )

    parser.add_argument("-l1", "--layer1", type=str, required=True,
                         help="Bottom monolayer .fdf (stays fixed).")
    parser.add_argument("-l2", "--layer2", type=str, required=True,
                         help="Top monolayer .fdf (rigidly slid across the grid). Pass the same "
                              "file as --layer1 for a material sliding against itself.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template (kgrid, basis, XC, %%include structure.fdf, "
                              "etc.) -- copied into every shift_II_JJ/ folder; single-point SCF "
                              "enforcement (and the optional --d3 dispersion correction) is "
                              "added via a config_extra.fdf %%include, not embedded inline.")
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
                         help="Interlayer gap in Ang (default: 3.2, same as stb-2Dstacking). "
                              "Under --mode 3 (default) this is the FIXED value used across the "
                              "whole grid -- a rigid-shift protocol, standard in the literature "
                              "(see the [INFO] closest-contact line in the report for how close "
                              "high-energy registries got at this fixed distance). Under --mode "
                              "1/2 this is only the STARTING guess -- SIESTA (mode 1) or "
                              "MACE-MP-0 (mode 2) relaxes it per grid point from here.")
    parser.add_argument("--vacuum", type=float, default=None,
                         help="Target vacuum space in Ang. Inherits --layer1's by default.")
    parser.add_argument("-sm", "--strain_mode", choices=["top", "bottom", "sym"], default="top",
                         help="Strain distribution mode for the ZSL match (default: top).")
    parser.add_argument("-t", "--twist", type=float, default=0.0,
                         help="Twist angle of layer 2 in degrees (default: 0.0), FIXED for the "
                              "whole grid -- a different twist is a different physical system, "
                              "not swept here.")
    parser.add_argument("--d3", dest="d3", action="store_true", default=True,
                         help="Force the Grimme DFT-D3 dispersion correction (via "
                              "config_extra.fdf) in every generated calc.fdf -- physically "
                              "relevant here since the interlayer binding swept across the grid "
                              "is a van-der-Waals-dominated interaction that plain GGA misses "
                              "(same correction stb-adsorb/stb-neb force by default elsewhere in "
                              "the suite). ON by default. Applies in every --mode.")
    parser.add_argument("--no-d3", dest="d3", action="store_false",
                         help="Disable the DFT-D3 dispersion correction (not recommended for "
                              "this workflow -- see --d3).")

    parser.add_argument("--mode", type=int, choices=[1, 2, 3], default=1,
                         help="Interlayer-gap strategy, per grid point (default: 1):\n"
                              "  1 = SIESTA relaxes z for real: writes one restricted CG "
                              "relaxation per grid point (x,y frozen for every atom, only z "
                              "free -- see build_z_relax_block/--relax-z-steps). Most accurate "
                              "(same level of theory as the reported energy), most expensive "
                              "(a real relaxation, not one SCF, per point).\n"
                              "  2 = MACE-MP-0 relaxes z first (same x,y-frozen/z-free "
                              "definition as mode 1, cheap ML engine), THEN writes a plain "
                              "single-point SIESTA folder at the resulting geometry -- middle "
                              "ground between cost and accuracy. Needs the optional 'ml' extra.\n"
                              "  3 = fixed --gap, plain single-point SIESTA, no relaxation of "
                              "any kind (today's behavior, the cheapest of the three).")
    parser.add_argument("--relax-z-steps", type=int, default=200,
                         help="With --mode 1: MD.Steps for the restricted SIESTA z-relaxation "
                              "(default: 200). Ignored in modes 2/3.")

    parser.add_argument("--scan", choices=["surface", "x", "y", "xy"], default="surface",
                         help="Sweep shape (default: surface):\n"
                              "  surface = full 2D grid (-nx x -ny points), the gamma-surface -- "
                              "output is the pm3d map/heatmap already documented above.\n"
                              "  x       = 1D sweep along shift_x only (shift_y fixed at 0) -- "
                              "-n points, output is a line plot (E vs. distance along x).\n"
                              "  y       = 1D sweep along shift_y only (shift_x fixed at 0) -- "
                              "-n points, output is a line plot (E vs. distance along y).\n"
                              "  xy      = 1D sweep along the diagonal shift_x == shift_y -- "
                              "-n points, output is a line plot (E vs. distance along the "
                              "diagonal).\n"
                              "-nx/-ny only apply to 'surface'; -n/--scan-points only applies to "
                              "'x'/'y'/'xy'.")
    parser.add_argument("-n", "--scan-points", type=int, default=15,
                         help="With --scan x/y/xy: number of points along the 1D sweep (default: "
                              "15). Ignored under --scan surface (see -nx/-ny instead).")

    parser.add_argument("-nx", "--grid-nx", type=int, default=7,
                         help="With --scan surface: grid resolution along shift_x, this many "
                              "lateral shifts covering [0, 1) (default: 7). Independent of "
                              "--grid-ny -- an anisotropic lattice (a != b, e.g. a rectangular "
                              "or low-symmetry 2D cell) often warrants a different sampling "
                              "density per axis rather than forcing a square grid. Ignored under "
                              "--scan x/y/xy (see -n/--scan-points instead).")
    parser.add_argument("-ny", "--grid-ny", type=int, default=7,
                         help="With --scan surface: grid resolution along shift_y (default: 7). "
                              "See --grid-nx. Ignored under --scan x/y/xy.")

    parser.add_argument("--ml-prerelax-layers", action="store_true",
                         help="Relax each monolayer's positions (cell fixed) with MACE-MP-0 "
                              "before stacking -- a cheap safety net when you're not fully "
                              "certain --layer1/--layer2 are already at equilibrium. Needs the "
                              "optional 'ml' extra.")
    parser.add_argument("--ml-preview", action="store_true",
                         help="Evaluate the whole grid's single-point energy on MACE-MP-0 (fast, "
                              "no SIESTA) and write stackingfault_ml_preview.png -- a quick sanity "
                              "check of the gamma-surface's shape before committing to the full "
                              "Nx x Ny DFT single-points. Written alongside the normal SIESTA "
                              "folders, not instead of them. Needs the optional 'ml' extra.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small",
                         help="MACE-MP-0 model size (default: small). Ignored if "
                              "--ml-custom-model is given.")
    parser.add_argument("--ml-custom-model", default=None, metavar="PATH",
                         help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                              "potential -- e.g. one fine-tuned on your own SIESTA data via "
                              "stb-mlffAnalysis. Overrides --ml-model.")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--ml-fmax", type=float, default=0.05,
                         help="Force convergence for --ml-prerelax-layers and --mode 2's "
                              "z-only MACE relaxation (default: 0.05).")

    parser.add_argument("-O", "--output-dir", type=str, default=".",
                         help=f"Root directory (default: current directory) -- everything lives "
                              f"under its '{RUN_SUBDIR}/' subfolder (see --save-report): every "
                              f"'shift_II_JJ/' grid folder under '{RUN_SUBDIR}/{POSITIONS_SUBDIR}/' "
                              f"(a sibling of stb-stackingfaultBsse's own '{RUN_SUBDIR}/bsse/'), "
                              f"the manifest, and the optional report directly under "
                              f"'{RUN_SUBDIR}/' -- same self-contained-run-folder convention as "
                              "stb-adsorb's 'sites/'+'bsse/' and stb-neb's 'neb_run/'.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the full run narrative to <output-dir>/{RUN_SUBDIR}/"
                              f"{REPORT_FILE}. Off by default -- the machine-readable "
                              f"{MANIFEST_FILE} that stb-stackingfaultAnalysis actually reads is "
                              "always written regardless of this flag.")
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

    if args.scan == "surface":
        if args.grid_nx < 2:
            parser.error("-nx/--grid-nx must be >= 2 (need at least 2 points for a sweep).")
        if args.grid_ny < 2:
            parser.error("-ny/--grid-ny must be >= 2 (need at least 2 points for a sweep).")
    elif args.scan_points < 2:
        parser.error("-n/--scan-points must be >= 2 (need at least 2 points for a sweep).")
    if args.ml_custom_model and not os.path.isfile(args.ml_custom_model):
        parser.error(f"--ml-custom-model file not found: {args.ml_custom_model}")

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
    calc_text, include_status, old_include_name = structure_io.fix_structure_include(calc_text)

    output_root = args.output_dir
    run_root = os.path.join(output_root, RUN_SUBDIR)
    positions_root = os.path.join(run_root, POSITIONS_SUBDIR)
    os.makedirs(positions_root, exist_ok=True)
    report_path = os.path.join(run_root, REPORT_FILE) if args.save_report else None

    f_out = open(report_path, 'w') if report_path else None

    print_dual(f"{color_text('===== STACKING-FAULT (GAMMA-SURFACE) SETUP REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"  Layer 1         : {args.layer1}", f_out)
    print_dual(f"  Layer 2         : {args.layer2}", f_out)
    print_dual(f"  Calc template   : {args.calc}", f_out)
    if include_status == "fixed":
        print_dual(color_text(
            f"  [NOTE] --calc's structure include was '%include {old_include_name}' -- "
            "auto-corrected to '%include structure.fdf' in every generated shift_II_JJ/calc.fdf "
            "(every folder's structure is always written under that exact name).", 'yellow'), f_out)
    elif include_status == "not_found":
        print_dual(color_text(
            "  [NOTE] Could not find a '%include ...struct....fdf'-style line in --calc to "
            "auto-correct -- make sure it '%include structure.fdf' itself so SIESTA picks up "
            "the generated geometry.", 'yellow'), f_out)
    print_dual(f"  Pseudo dir      : {args.pseudo_dir or '(none)'}", f_out)
    print_dual(f"  Output dir      : {output_root} (run folder: {run_root})", f_out)
    if args.scan == "surface":
        print_dual(f"  Grid            : {args.grid_nx} x {args.grid_ny} (surface)"
                    + (" (asymmetric)" if args.grid_nx != args.grid_ny else ""), f_out)
    else:
        scan_desc = {"x": "shift_x only", "y": "shift_y only", "xy": "diagonal shift_x=shift_y"}
        print_dual(f"  Grid            : {args.scan_points} points, 1D ({scan_desc[args.scan]})",
                    f_out)
    mode_labels = {
        1: "SIESTA relaxes z for real (restricted CG, x/y frozen)",
        2: "MACE-MP-0 relaxes z, then SIESTA single-point",
        3: "fixed gap, plain SIESTA single-point",
    }
    print_dual(f"  Mode            : {args.mode} ({mode_labels[args.mode]})", f_out)
    if args.mode == 1:
        print_dual(color_text(
            "  [NOTE] --mode 1's Geometry.Constraints syntax (freezing x,y per atom) was NOT "
            "verified against a real SIESTA run in this environment -- check the .out/.FA of "
            "ONE grid point (x,y displacement should be ~0, z should have moved) before "
            "trusting a full production sweep.", 'yellow'), f_out)
    gap_role = "starting guess, relaxed per-point" if args.mode in (1, 2) else "fixed across the grid"
    print_dual(f"  Gap             : {args.gap} Ang ({gap_role})", f_out)
    print_dual(f"  Twist (fixed)   : {args.twist} deg", f_out)
    print_dual(f"  D3 dispersion   : {'yes (config_extra.fdf)' if args.d3 else 'no'}", f_out)
    print_dual(f"  max_area/max_strain/match_id: {args.max_area} / {args.max_strain} / {args.match_id}", f_out)
    print_dual(f"  ML pre-relax layers : {'yes' if args.ml_prerelax_layers else 'no'}", f_out)
    print_dual(f"  ML preview          : {'yes' if args.ml_preview else 'no'}", f_out)
    print_dual(f"  Save report         : {'yes -> ' + report_path if args.save_report else 'no'}", f_out)

    library_warnings = []  # collected via capture_library_noise, reported in the last section

    calc_mace = None
    uses_mace = args.ml_prerelax_layers or args.mode == 2 or args.ml_preview
    ml_model_arg = args.ml_custom_model if args.ml_custom_model else args.ml_model
    if uses_mace:
        print_dual(f"  ML model        : "
                    + (f"custom ({args.ml_custom_model})" if args.ml_custom_model
                       else f"MACE-MP-0 ({args.ml_model})")
                    + (" + D3(BJ) dispersion (matches --d3)" if args.d3 else " (no D3 dispersion)"),
                    f_out)
        with capture_library_noise(library_warnings, "MACE calculator setup"):
            require_mace()
            from stb.core import mace_relax
            calc_mace = mace_relax.get_calculator(model=ml_model_arg, device=args.ml_device,
                                                   dispersion=args.d3)

    if args.ml_prerelax_layers:
        print_dual(f"\n{color_text('ML pre-relax:', 'cyan')} relaxing each monolayer "
                    "(positions only) with MACE-MP-0 ...", f_out)
        for label, pmg in (("layer 1", layer1_pmg), ("layer 2", layer2_pmg)):
            ase_atoms = AseAtomsAdaptor.get_atoms(pmg)
            with capture_library_noise(library_warnings, f"MACE pre-relax ({label})"):
                converged, steps = mace_relax.relax(ase_atoms, calc_mace, fmax=args.ml_fmax, max_steps=200)
            relaxed_pmg = AseAtomsAdaptor.get_structure(ase_atoms)
            if label == "layer 1":
                layer1_pmg = relaxed_pmg
            else:
                layer2_pmg = relaxed_pmg
            print_dual(f"  {'Converged' if converged else 'Hit step cap, not fully converged'} "
                        f"({label}) after {steps} step(s).", f_out)

    print_section('[1] ZSL MATCH', f_out)
    # NOT wrapped in capture_library_noise: find_zsl_match (core/heterostructure.py)
    # prints its own user-facing [INFO] status lines via plain print(), which
    # capture_library_noise would otherwise swallow into [4] LIBRARY WARNINGS
    # instead of leaving them visible here where they belong.
    t_mat1, t_mat2, best_match_data = find_zsl_match(
        layer1_pmg, layer2_pmg, max_area=args.max_area, max_strain=args.max_strain,
        match_id=args.match_id, interactive=False, twist_angle=args.twist)
    print_dual(f"  Selected match ID {args.match_id}: area {best_match_data['area']:.2f} Ang^2, "
                f"strain {best_match_data['strain']:.2f}%, angular strain "
                f"{best_match_data['angle_strain']:.2f} deg.", f_out)
    n_layer1_supercell = len(layer1_pmg) * round(np.linalg.det(t_mat1))
    n_layer2_supercell = len(layer2_pmg) * round(np.linalg.det(t_mat2))
    n_total = n_layer1_supercell + n_layer2_supercell
    scan_points = build_scan_points(args.scan, args.grid_nx, args.grid_ny, args.scan_points)
    n_grid_points = len(scan_points)
    print_dual(f"  Each grid point: {n_total} atoms ({n_layer1_supercell} layer 1 + "
                f"{n_layer2_supercell} layer 2). {n_grid_points} grid point(s) total -- "
                f"{n_grid_points} independent single-point SIESTA runs at {n_total} atoms each.",
                f_out)

    print_section('[2] GRID FOLDERS', f_out)
    if args.mode == 2:
        print_dual(f"  {color_text('Mode 2:', 'cyan')} relaxing z (x,y frozen) with MACE-MP-0 "
                    f"from {args.gap} Ang at every grid point ...", f_out)
    elif args.mode == 1:
        print_dual(f"  {color_text('Mode 1:', 'cyan')} writing a restricted SIESTA z-relaxation "
                    f"(x,y frozen, up to {args.relax_z_steps} CG steps) at every grid point, "
                    f"starting from {args.gap} Ang ...", f_out)
    manifest_rows = []  # (label, i, j, shift_x, shift_y, gap_used)
    closest_contact = None
    closest_label = None
    ml_grid = []  # (i, j, energy) -- only populated if --ml-preview
    gaps_used = []
    for label, i, j, shift_x, shift_y in scan_points:
        hetero, n_layer1_atoms, _max_strain_val = build_stacked_structure(
            layer1_pmg, layer2_pmg, t_mat1, t_mat2, shift_x, shift_y, args.gap,
            target_vacuum=args.vacuum, strain_mode=args.strain_mode)

        ml_energy = None
        if args.mode == 2:
            with capture_library_noise(library_warnings, "MACE z-relax (mode 2)"):
                hetero, ml_energy = relax_gap_ml(hetero, calc_mace, fmax=args.ml_fmax)
            gap_used = effective_gap(hetero, n_layer1_atoms)
        else:
            gap_used = args.gap
        gaps_used.append(gap_used)

        if args.ml_preview:
            if ml_energy is None:
                with capture_library_noise(library_warnings, "MACE preview evaluation"):
                    atoms = AseAtomsAdaptor.get_atoms(hetero)
                    atoms.calc = calc_mace
                    ml_energy = atoms.get_potential_energy()
            ml_grid.append((i, j, ml_energy))

        contact = min_interlayer_distance(hetero, n_layer1_atoms)
        if contact is not None and (closest_contact is None or contact < closest_contact):
            closest_contact = contact
            closest_label = label

        grid_dir = os.path.join(positions_root, label)
        write_grid_folder(grid_dir, hetero, calc_text, species_meta, args.pseudo_dir,
                           d3=args.d3, mode=args.mode, relax_z_steps=args.relax_z_steps)
        gap_note = f", gap: {gap_used:.3f} Ang (ML-relaxed)" if args.mode == 2 else ""
        print_dual(f"  {color_text('[OK]', 'green')} {grid_dir} "
                    f"(shift: {shift_x:.4f}, {shift_y:.4f}{gap_note})", f_out)
        manifest_rows.append((label, i, j, shift_x, shift_y, gap_used))

    print_section('[3] SUMMARY & NEXT STEPS', f_out)
    print_dual(f"  {len(manifest_rows)} grid folder(s) written under '{run_root}'.", f_out)
    if closest_contact is not None:
        print_dual(f"  [INFO] Closest interlayer contact anywhere in the grid: "
                    f"{closest_contact:.3f} Ang at {closest_label} -- expected to be small at "
                    "eclipsed/high-energy registries, not necessarily a problem.", f_out)
    if args.mode == 2:
        print_dual(f"  [INFO] MACE-relaxed gap ranged {min(gaps_used):.3f} - "
                    f"{max(gaps_used):.3f} Ang across the grid (starting guess was {args.gap} "
                    "Ang).", f_out)
    if args.ml_preview and ml_grid:
        preview_path = os.path.join(run_root, "stackingfault_ml_preview.png")
        with capture_library_noise(library_warnings, "matplotlib preview plot"):
            if args.scan == "surface":
                write_ml_preview_plot(ml_grid, args.grid_nx, args.grid_ny, preview_path)
            else:
                write_ml_preview_line_plot(ml_grid, args.scan, preview_path)
        ml_min = min(ml_grid, key=lambda r: r[2])
        ml_max = max(ml_grid, key=lambda r: r[2])
        print_dual(f"  {color_text('[Saved]', 'cyan')} {preview_path}", f_out)
        print_dual(f"  ML preview: predicted equilibrium at shift_{ml_min[0]:02d}_{ml_min[1]:02d}, "
                    f"predicted corrugation {ml_max[2] - ml_min[2]:.4f} eV (MACE-MP-0 -- NOT a "
                    "DFT-level number, use stb-stackingfaultAnalysis for the real one).", f_out)

    manifest = {
        "version": 3,
        "grid_nx": args.grid_nx,
        "grid_ny": args.grid_ny,
        "scan": args.scan,
        "scan_points": args.scan_points,
        "gap_nominal": args.gap,
        "twist_deg": args.twist,
        "mode": args.mode,
        "n_layer1_atoms": int(n_layer1_supercell),
        "rows": [
            {"label": label, "i": i, "j": j, "shift_x": shift_x, "shift_y": shift_y, "gap": gap_used}
            for label, i, j, shift_x, shift_y, gap_used in manifest_rows
        ],
    }
    manifest_path = os.path.join(run_root, MANIFEST_FILE)
    with open(manifest_path, 'w') as f_manifest:
        json.dump(manifest, f_manifest, indent=2)
    print_dual(f"  {color_text('[Saved]', 'cyan')} Manifest      -> {manifest_path} "
                "(machine-readable, read by stb-stackingfaultAnalysis)", f_out)
    if report_path:
        print_dual(f"  {color_text('[Saved]', 'cyan')} Report        -> {report_path}", f_out)
    run_instruction = ("restricted relaxation: up to " + str(args.relax_z_steps) + " CG steps, "
                        "x/y frozen, via config_extra.fdf" if args.mode == 1 else
                        "single-point: MD.Steps forced to 0 via config_extra.fdf")
    print_dual(f"  Run SIESTA ({run_instruction}) in every '{positions_root}/shift_II_JJ/' "
                f"folder, then run: stb-stackingfaultAnalysis --dir {run_root}", f_out)

    print_section('[4] LIBRARY WARNINGS', f_out)
    if library_warnings:
        print_dual(color_text(
            "  Messages emitted by external libraries (pymatgen, MACE/torch, matplotlib) during "
            "this run -- collected here instead of interleaved with the report above; harmless "
            "in almost every case, but worth a look if a section above looks suspicious.",
            'cyan'), f_out)
        for entry in library_warnings:
            print_dual(f"  {entry}", f_out)
    else:
        print_dual("  No library warnings.", f_out)

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} {len(manifest_rows)} grid folder(s) written under "
          f"'{run_root}'.")
    print(f"Manifest: {manifest_path}")
    if report_path:
        print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
