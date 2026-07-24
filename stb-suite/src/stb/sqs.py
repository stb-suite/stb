#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import sys
import os
import math
import time
import argparse
from datetime import datetime
from fractions import Fraction

import numpy as np
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table, run_with_spinner
from stb.core.deps import require_icet, require_mace

require_icet()
from icet import ClusterSpace
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.advanced_transformations import SQSTransformation

REPORT_FILE = "stb_sqs_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-fetch/stb-kgrid/stb-mlrelax/
# stb-supercell (core/kspace.py's other callers).
VACUUM_GAP_ANG = 10.0


def parse_composition(spec):
    """Parses 'Fe:0.5,Ni:0.5' into {'Fe': 0.5, 'Ni': 0.5}, validated to sum to 1.0."""
    composition = {}
    for entry in spec.split(","):
        try:
            symbol, fraction = entry.split(":")
            composition[symbol.strip()] = float(fraction)
        except ValueError:
            raise ValueError(f"--composition entries must be 'Symbol:fraction', got '{entry}' in '{spec}'.")

    total = sum(composition.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"--composition fractions must sum to 1.0, got {total} ({spec}).")
    return composition


def parse_cluster_cutoffs(spec):
    """Parses '2:3,3:2' into {2: 3, 3: 2} (cluster size -> neighbor shell)."""
    cutoffs = {}
    for entry in spec.split(","):
        try:
            size, shell = entry.split(":")
            cutoffs[int(size)] = int(shell)
        except ValueError:
            raise ValueError(f"--cluster-cutoffs entries must be 'size:shell', got '{entry}' in '{spec}'.")
    return cutoffs


def minimal_scaling(pmg, composition, cluster_size_and_shell):
    """Smallest integer --scaling that gives every species in `composition` an exact
    atom count in the primitive cell icet auto-detects for `pmg` (which may be smaller
    than the literal input structure, e.g. a 4-atom fcc cell reduces to a 1-atom
    primitive). Using a --scaling that isn't a multiple of this value makes the target
    concentration unreachable at that supercell size: icet's enumeration search then
    finds zero candidate structures and pymatgen crashes with an unrelated
    AttributeError instead of a clear message, and its monte_carlo search silently
    skips that size. `pmg` must already have the disordered composition applied
    (i.e. after the `pmg.replace(idx, composition)` loop).

    Returns a dict with the minimal scaling plus the icet internals behind it
    (n_primitive, mult_factor), so callers can explain the number to the user.
    """
    clusters = SQSTransformation._sqs_cluster_estimate(pmg, cluster_size_and_shell)
    cutoffs_list = []
    for i in range(2, max(clusters.keys()) + 1):
        clusters.setdefault(i, 0.0)
        cutoffs_list.append(clusters[i])

    ordered = pmg.copy()
    dummy_symbol = next(iter(ordered.composition))
    # Per-site replace (not replace_species) -- replace_species remaps
    # existing species symbols but leaves a site's occupancy dict merged
    # from whatever partial species it already had; for a 3+-species
    # --composition on the disordered sublattice, that leaves the dummy
    # structure still fractionally multi-occupied on those sites (only
    # ever collapsed to a single species by chance for exactly 2 species),
    # and AseAtomsAdaptor.get_atoms() then raises "ASE Atoms only supports
    # ordered structures". A direct per-site replace always fully
    # overwrites each site with 100% dummy_symbol, regardless of how many
    # species it was disordered with.
    for i in range(len(ordered)):
        ordered.replace(i, dummy_symbol)
    ordered_atoms = AseAtomsAdaptor.get_atoms(ordered)

    chemical_symbols = [list(site.species.as_dict()) for site in pmg]
    cluster_space = ClusterSpace(structure=ordered_atoms, cutoffs=cutoffs_list, chemical_symbols=chemical_symbols)

    n_primitive = len(cluster_space.primitive_structure)
    sublattices = cluster_space.get_sublattices(cluster_space.primitive_structure)
    target_symbols = set(composition)
    matching = [sl for sl in sublattices if set(sl.chemical_symbols) == target_symbols]
    mult_factor = Fraction(len(matching[0].indices), n_primitive) if matching else Fraction(1, 1)

    denominators = [
        (Fraction(fraction).limit_denominator(1000) * mult_factor).denominator
        for fraction in composition.values()
    ]
    return {
        "minimal": math.lcm(*denominators),
        "n_primitive": n_primitive,
        "mult_factor": mult_factor,
    }


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as
    stacking2D.py/supercell.py/slab.py/nanotube.py/defect.py's own _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def _describe_structure(pmg_structure, vacuum_axes, f_out):
    """Formula/atoms/dimensionality/cell parameters -- same field set as
    supercell.py's own [1] INPUT STRUCTURE section."""
    print_dual(f"Formula        : {pmg_structure.composition.reduced_formula}", f_out)
    print_dual(f"Atoms          : {len(pmg_structure)}", f_out)
    print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    a, b, c, alpha, beta, gamma = pmg_structure.lattice.parameters
    print_dual(f"Cell a,b,c     : {a:.4f}, {b:.4f}, {c:.4f} Ang", f_out)
    print_dual(f"Cell angles    : {alpha:.2f}, {beta:.2f}, {gamma:.2f} deg", f_out)
    print_dual(f"Cell volume    : {pmg_structure.lattice.volume:.4f} Ang^3", f_out)
    if sum(vacuum_axes) == 0:
        print_dual(f"Density        : {len(pmg_structure) / pmg_structure.lattice.volume:.4f} atoms/Ang^3", f_out)


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a
    space-group label -- same shape as supercell.py/slab.py/nanotube.py/
    defect.py's own _validate_structure(), wrapped in try/except by the
    caller (a validation failure is reported, never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group    : {sg_label}", f_out)
    return sg_label


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Generates a Special Quasirandom Structure (SQS) for a substitutional alloy.", 'bold')}
Picks one sublattice (--sublattice) and disorders it with a target
composition, then searches for the atomic arrangement that best mimics an
ideal random alloy (via icet's Monte Carlo or enumeration SQS search).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5\n"
               "  %(prog)s -f ni.fdf --sublattice Ni --composition Ni:0.75,Fe:0.25 --scaling 8 --method enumeration\n"
               "  %(prog)s -f ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --ml-relax --save-report --view\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf).")
    parser.add_argument("--sublattice", type=str, required=True,
                        help="Existing species whose sites become the disordered alloy sublattice.")
    parser.add_argument("--composition", type=str, required=True,
                        help="Target composition on that sublattice, e.g. 'Fe:0.5,Ni:0.5' (must sum to 1.0). "
                             "Practically limited to 2 species: pymatgen's own SQSTransformation/IcetSQS "
                             "wrapper (pymatgen/io/icet.py) builds its internal reference structure the same "
                             "way this tool's minimal_scaling() used to (a replace_species() call that leaves "
                             "a site's occupancy dict merged rather than collapsed for 3+ species), and raises "
                             "'ASE Atoms only supports ordered structures' during the search itself for a 3+ "
                             "-species composition -- verified upstream in the installed pymatgen, not "
                             "something this tool can work around without bypassing SQSTransformation "
                             "entirely.")
    parser.add_argument("--scaling", type=int, default=None,
                        help="icet supercell-size control (positive integer). Multiplies icet's own "
                             "internally-detected primitive unit, not necessarily the literal input atom "
                             "count -- read the actual resulting atom count from the output. Default: "
                             "auto-detect the smallest scaling that gives every species in --composition an "
                             "exact atom count (a --scaling that doesn't divide evenly is rejected, since icet "
                             "would otherwise fail to find any candidate structure at that size).")
    parser.add_argument("--method", choices=["monte_carlo", "enumeration"], default="monte_carlo",
                        help="SQS search method (default: monte_carlo).")
    parser.add_argument("--instances", type=int, default=1,
                        help="Number of parallel search instances/processes (default: 1, serial). For "
                             "--method enumeration, icet's own enumeration of candidate structures is single "
                             "-threaded regardless of this value -- only the (comparatively cheap) scoring of "
                             "already-generated candidates is split across instances, so raising it rarely "
                             "helps. For monte_carlo, each instance runs an independent full search and the "
                             "best is kept, which can help on a machine with real spare cores -- but every "
                             "instance re-pays icet/pymatgen's startup cost (imports, primitive-cell detection "
                             "via spglib), so on a loaded or virtualized machine more instances can end up "
                             "*slower* in wall time despite using more CPU overall. Benchmark before raising "
                             "this; don't assume it helps.")
    parser.add_argument("--temperature", type=float, default=1.0,
                        help="Monte Carlo starting temperature (default: 1.0).")
    parser.add_argument("--mc-steps", type=int, default=None,
                        help="Caps the number of Monte Carlo steps (--method monte_carlo only). "
                             "Default: icet's own default.")
    parser.add_argument("--cluster-cutoffs", type=str, default=None,
                        help="Cluster size:shell pairs, e.g. '2:3,3:2' (default: icet's own {2:3, 3:2, 4:1}).")
    parser.add_argument("-sp", "--symprec", type=float, default=0.01,
                        help="Symmetry tolerance in Angstroms for the before/after "
                             "symmetry analysis (default: 0.01).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax the SQS structure found with a MACE potential "
                             "(needs the optional 'ml' extra: pip install stb_suite[ml]) "
                             "before writing it out -- positions only by default. "
                             "Off by default.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the cell -- any vacuum-padded "
                             "axis always stays exactly fixed. Only valid together with "
                             "--ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax, "
                             "instead of a MACE-MP-0 foundation size.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the input "
                             "structure and the final SQS structure (page through frames in "
                             "ase-gui) after writing the output file. Needs a display. "
                             "Off by default.")

    parser.add_argument("-o", "--output", type=str, default="sqs.fdf",
                        help="Output .fdf file name (default: sqs.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-sqs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.scaling is not None and args.scaling < 1:
        parser.error("--scaling must be >= 1.")
    if args.ml_relax_cell and not args.ml_relax:
        parser.error("--ml-relax-cell is only valid together with --ml-relax.")
    if (args.custom_model or args.model != "small") and not args.ml_relax:
        parser.error("--model/--custom-model are only valid together with --ml-relax.")
    if args.ml_relax:
        require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-SQS REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Sublattice     : {args.sublattice}", f_out)
    print_dual(f"Composition    : {args.composition}", f_out)
    print_dual(f"Scaling        : {args.scaling if args.scaling is not None else 'auto-detect'}", f_out)
    print_dual(f"Method         : {args.method}", f_out)
    print_dual(f"Instances      : {args.instances}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"ML pre-relax   : {model_desc if args.ml_relax else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell     : {'yes (vacuum axes fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    if not os.path.exists(args.filename):
        _fail(f"File '{args.filename}' not found.", f_out)

    try:
        composition = parse_composition(args.composition)
    except ValueError as e:
        _fail(str(e), f_out)

    cluster_size_and_shell = None
    if args.cluster_cutoffs:
        try:
            cluster_size_and_shell = parse_cluster_cutoffs(args.cluster_cutoffs)
        except ValueError as e:
            _fail(str(e), f_out)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        _fail(str(e), f_out)

    pmg_structure = structure_io.to_pymatgen(structure)

    print_section("[1] INPUT STRUCTURE", f_out)
    frac_coords_before = [site.frac_coords for site in pmg_structure]
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_coords_before, pmg_structure.lattice.matrix, VACUUM_GAP_ANG)
    _describe_structure(pmg_structure, vacuum_axes_before, f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    try:
        _validate_structure(pmg_structure, vacuum_axes_before, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[3] DISORDERED SUBLATTICE", f_out)
    site_indices = [i for i, site in enumerate(pmg_structure) if site.specie.symbol == args.sublattice]
    if not site_indices:
        _fail(f"species '{args.sublattice}' not found in the structure.", f_out)

    print_dual(f"Species            : {args.sublattice}", f_out)
    print_dual(f"Sites              : {len(site_indices)} (of {len(pmg_structure)} atoms in the input cell)", f_out)
    print_dual(f"Target composition : {composition}", f_out)

    pmg = pmg_structure.copy()
    for idx in site_indices:
        pmg.replace(idx, composition)

    print_section("[4] SUPERCELL SIZING (icet)", f_out)
    scaling_info = minimal_scaling(pmg, composition, cluster_size_and_shell)
    min_scaling = scaling_info["minimal"]
    print_dual(f"icet primitive cell    : {scaling_info['n_primitive']} atom(s)", f_out)
    print_dual(f"Sublattice mult.       : {scaling_info['mult_factor']} "
               "(sublattice sites / primitive atoms)", f_out)
    print_dual(f"Minimal valid scaling  : {min_scaling} "
               "(smallest multiplier giving every species an exact atom count)", f_out)

    if args.scaling is None:
        scaling = min_scaling
        print_dual(f"Scaling used           : {scaling} (auto-detected)", f_out)
    elif args.scaling % min_scaling != 0:
        _fail(
            f"--scaling {args.scaling} is incompatible with composition '{args.composition}': "
            f"icet's auto-detected primitive cell needs a multiple of {min_scaling} "
            f"(e.g. {min_scaling}, {2 * min_scaling}, {3 * min_scaling}, ...) for every species to get an "
            f"exact atom count. Rerun with --scaling {min_scaling} or omit --scaling to auto-detect it.",
            f_out)
    else:
        scaling = args.scaling
        print_dual(f"Scaling used           : {scaling} (user-specified)", f_out)

    print_section("[5] SQS SEARCH", f_out)
    print_dual(f"Method             : {args.method}", f_out)
    instances = max(1, args.instances)
    print_dual(f"Instances requested: {instances}", f_out)
    if args.method == "monte_carlo":
        print_dual(f"MC start temperature: {args.temperature}", f_out)
        print_dual(f"MC steps cap        : {args.mc_steps if args.mc_steps is not None else 'icet default'}", f_out)
    print_dual(f"Cluster cutoffs    : "
               f"{cluster_size_and_shell if cluster_size_and_shell else 'icet default (2:3, 3:2, 4:1)'}", f_out)

    icet_sqs_kwargs = {}
    if args.mc_steps is not None:
        icet_sqs_kwargs["n_steps"] = args.mc_steps

    def search(n_instances):
        transformation = SQSTransformation(
            scaling=scaling,
            sqs_method=f"icet-{args.method}",
            instances=n_instances,
            temperature=args.temperature,
            cluster_size_and_shell=cluster_size_and_shell,
            icet_sqs_kwargs=icet_sqs_kwargs,
        )
        return run_with_spinner(transformation.apply_transformation, pmg, return_ranked_list=1, label="Searching")

    used_instances = instances
    instances_fallback_note = None

    print_dual(color_text("Searching for SQS...", 'yellow'), f_out)
    try:
        ranked = search(instances)
    except ValueError as e:
        _fail(str(e), f_out)
    except AttributeError as e:
        # Known icet/pymatgen edge case: when --instances outnumbers the enumerated
        # candidate structures, the leftover processes get an empty chunk and
        # IcetSQS.run() unconditionally tries to convert their (None) "best" structure,
        # crashing with this unrelated-looking AttributeError instead of a clear message.
        if instances <= 1 or "get_chemical_symbols" not in str(e):
            raise
        instances_fallback_note = (
            f"{instances} parallel instances outnumbered the enumerated candidate "
            "structures for this small a cell; retried with 1 instance."
        )
        print_dual(color_text(f"[WARNING] {instances_fallback_note}", 'yellow'), f_out)
        used_instances = 1
        try:
            ranked = search(1)
        except ValueError as e2:
            _fail(str(e2), f_out)

    sqs_structure = ranked[0]["structure"]
    objective_function = ranked[0]["objective_function"]

    print_dual(f"Instances used     : {used_instances}", f_out)
    print_dual(f"Output formula     : {sqs_structure.composition.reduced_formula}", f_out)
    print_dual(f"Output atoms       : {len(sqs_structure)}", f_out)
    print_dual(f"Objective function : {objective_function:.6f}", f_out)

    ml_relax_info = None
    if args.ml_relax:
        print_section("[6] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc}", f_out)
        print_dual(f"Cell relaxation : "
                   f"{'in cell (vacuum axes fixed)' if args.ml_relax_cell else 'positions only'}", f_out)

        frac_coords_sqs = [site.frac_coords for site in sqs_structure]
        vacuum_axes_sqs = kspace.detect_vacuum_axes(frac_coords_sqs, sqs_structure.lattice.matrix, VACUUM_GAP_ANG)

        model_arg = args.custom_model if args.custom_model else args.model
        calc = mace_relax.get_calculator(model_arg)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        atoms = AseAtomsAdaptor.get_atoms(sqs_structure)
        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        f0 = float(np.abs(atoms.get_forces()).max())
        a0, b0, c0, _, _, _ = atoms.cell.cellpar()
        vol0 = atoms.get_volume()

        cell_mask = mace_relax.build_cell_mask(vacuum_axes_sqs) if args.ml_relax_cell else None
        t0 = time.time()
        converged, steps_used = mace_relax.relax(atoms, calc, cell_mask=cell_mask, fmax=0.05, max_steps=200)
        wall_time = time.time() - t0

        e1 = atoms.get_potential_energy()
        f1 = float(np.abs(atoms.get_forces()).max())
        a1, b1, c1, _, _, _ = atoms.cell.cellpar()

        print_dual(f"Steps used : {steps_used} "
                   f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
        print_dual(f"Wall time  : {wall_time:.1f} s", f_out)

        n_atoms = len(atoms)
        rows = [
            (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
              f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
            (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
        ]
        if args.ml_relax_cell:
            rows.append((["Lattice a,b,c (Ang)", f"{a0:.4f}, {b0:.4f}, {c0:.4f}",
                          f"{a1:.4f}, {b1:.4f}, {c1:.4f}",
                          f"max {100 * max(abs(a1 - a0) / a0, abs(b1 - b0) / b0, abs(c1 - c0) / c0):+.2f}%"], None))
            vol1 = atoms.get_volume()
            rows.append((["Volume (Ang^3)", f"{vol0:.4f}", f"{vol1:.4f}",
                          f"{100 * (vol1 - vol0) / vol0:+.2f}%"], None))
        print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

        sqs_structure = AseAtomsAdaptor.get_structure(atoms)
        ml_relax_info = (converged, steps_used, e0, e1)

    final_atoms = AseAtomsAdaptor.get_atoms(sqs_structure)

    print_section("[7] STRUCTURE VALIDATION (post-transform)", f_out)
    frac_coords_final = [site.frac_coords for site in sqs_structure]
    vacuum_axes_final = kspace.detect_vacuum_axes(frac_coords_final, sqs_structure.lattice.matrix, VACUUM_GAP_ANG)
    try:
        _validate_structure(sqs_structure, vacuum_axes_final, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[8] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    print_dual("Before = ordered input structure. After = final (disordered) SQS structure.", f_out)
    before_info = core_symmetry.symmetry_summary(pmg_structure, args.symprec, VACUUM_GAP_ANG)
    after_info = core_symmetry.symmetry_summary(sqs_structure, args.symprec, VACUUM_GAP_ANG)
    if "Error" in before_info or "Error" in after_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
        print_dual(f"  Before: {before_info.get('Error', 'OK')}", f_out)
        print_dual(f"  After : {after_info.get('Error', 'OK')}", f_out)
    else:
        properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
        rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                for prop in properties]
        print_table(["Property", "Before", "After"], rows, f_out)

    print_section("[9] WRITING OUTPUT FILE", f_out)
    header_comment = [
        "SQS structure built by stb-sqs from an input structure.",
        f"Input file: {args.filename}",
        f"Disordered sublattice: {args.sublattice} ({len(site_indices)} site(s)) -> {composition}.",
        f"Scaling: {scaling}x icet primitive cell ({scaling_info['n_primitive']} atom(s)), "
        f"{len(pmg_structure)} -> {len(sqs_structure)} atoms.",
        f"SQS search: icet-{args.method}, {used_instances} instance(s), "
        f"objective function = {objective_function:.6f}.",
    ]
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )

    species_meta = dict(structure.species_meta)
    for symbol in composition:
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(sqs_structure, species_meta=species_meta)
    structure_io.write_fdf(new_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[10] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.ICET]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[11] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # both frames (input structure vs. final SQS structure) so the user can
    # page through the actual comparison in ase-gui.
    if args.view:
        input_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        view_structure_interactive([input_atoms, final_atoms])


if __name__ == "__main__":
    main()
