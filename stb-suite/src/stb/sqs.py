#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.5"

import sys
import os
import math
import argparse
from fractions import Fraction
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro, run_with_spinner
from stb.core.deps import require_icet

require_icet()
from icet import ClusterSpace
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.advanced_transformations import SQSTransformation


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
    original_composition = ordered.composition.as_dict()
    dummy_symbol = next(iter(ordered.composition))
    ordered.replace_species({sym: dummy_symbol for sym in original_composition if sym != dummy_symbol})
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


def generate_txt_report(args, input_formula, input_atoms, site_indices, composition,
                         scaling, scaling_info, instances, used_instances, instances_fallback_note,
                         cluster_size_and_shell, sqs_structure, objective_function, output_atoms):
    def row(label, value):
        return f"{label:<22}: {value}"

    lines = []
    lines.append("========================================")
    lines.append("      SQS GENERATION REPORT")
    lines.append("========================================")
    lines.append(row("Input file", args.filename))
    lines.append(row("Input formula", input_formula))
    lines.append(row("Input atoms", input_atoms) + "\n")

    lines.append("== Disordered Sublattice ==")
    lines.append(row("Species", args.sublattice))
    lines.append(row("Sites", len(site_indices)))
    lines.append(row("Target composition", composition) + "\n")

    lines.append("== Supercell Size ==")
    lines.append(row("icet primitive cell", f"{scaling_info['n_primitive']} atom(s)"))
    lines.append(row("Sublattice mult.", f"{scaling_info['mult_factor']} (sublattice sites / primitive atoms)"))
    lines.append(row("Minimal valid scaling",
                      f"{scaling_info['minimal']} (smallest multiplier giving every species an exact atom count)"))
    lines.append(row("Scaling used",
                      f"{scaling} ({'auto-detected' if args.scaling is None else 'user-specified'})") + "\n")

    lines.append("== Search Settings ==")
    lines.append(row("Method", args.method))
    lines.append(row("Instances requested", instances))
    lines.append(row("Instances used", used_instances))
    if instances_fallback_note:
        lines.append(row("Note", instances_fallback_note))
    if args.method == "monte_carlo":
        lines.append(row("MC start temperature", args.temperature))
        lines.append(row("MC steps cap", args.mc_steps if args.mc_steps is not None else "icet default"))
    lines.append(row("Cluster cutoffs",
                      cluster_size_and_shell if cluster_size_and_shell else "icet default (2:3, 3:2, 4:1)") + "\n")

    lines.append("== Result ==")
    lines.append(row("Output formula", sqs_structure.composition.reduced_formula))
    lines.append(row("Output atoms", output_atoms))
    lines.append(row("Objective function", f"{objective_function:.6f}"))
    lines.append(row("Output structure", args.output))
    return "\n".join(lines) + "\n"


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
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf).")
    parser.add_argument("--sublattice", type=str, required=True,
                        help="Existing species whose sites become the disordered alloy sublattice.")
    parser.add_argument("--composition", type=str, required=True,
                        help="Target composition on that sublattice, e.g. 'Fe:0.5,Ni:0.5' (must sum to 1.0).")
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
    parser.add_argument("-o", "--output", type=str, default="sqs.fdf",
                        help="Output .fdf file name (default: sqs.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-sqs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.scaling is not None and args.scaling < 1:
        parser.error("--scaling must be >= 1.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Generate a Special Quasirandom Structure (SQS):", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    try:
        composition = parse_composition(args.composition)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    cluster_size_and_shell = None
    if args.cluster_cutoffs:
        try:
            cluster_size_and_shell = parse_cluster_cutoffs(args.cluster_cutoffs)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    pmg = structure_io.to_pymatgen(structure)
    input_formula = pmg.composition.reduced_formula
    input_atoms = len(pmg)
    print(f"  {color_text('Input formula:', 'cyan')} {input_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {input_atoms}")

    site_indices = [i for i, site in enumerate(pmg) if site.specie.symbol == args.sublattice]
    if not site_indices:
        print(color_text(f"Error: species '{args.sublattice}' not found in the structure.", 'red'))
        sys.exit(1)

    print(f"  {color_text('Sublattice:', 'cyan')} {args.sublattice} ({len(site_indices)} site(s))")
    print(f"  {color_text('Target composition:', 'cyan')} {composition}")

    for idx in site_indices:
        pmg.replace(idx, composition)

    scaling_info = minimal_scaling(pmg, composition, cluster_size_and_shell)
    min_scaling = scaling_info["minimal"]
    if args.scaling is None:
        scaling = min_scaling
        print(f"  {color_text('Scaling:', 'cyan')} {scaling} (auto-detected)  "
              f"{color_text('Method:', 'cyan')} {args.method}")
    elif args.scaling % min_scaling != 0:
        print(color_text(
            f"Error: --scaling {args.scaling} is incompatible with composition '{args.composition}': "
            f"icet's auto-detected primitive cell needs a multiple of {min_scaling} "
            f"(e.g. {min_scaling}, {2 * min_scaling}, {3 * min_scaling}, ...) for every species to get an "
            f"exact atom count. Rerun with --scaling {min_scaling} or omit --scaling to auto-detect it.",
            'red'))
        sys.exit(1)
    else:
        scaling = args.scaling
        print(f"  {color_text('Scaling:', 'cyan')} {scaling}  {color_text('Method:', 'cyan')} {args.method}")

    instances = max(1, args.instances)
    print(f"  {color_text('Instances:', 'cyan')} {instances}")

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

    print(f"\n  {color_text('Searching for SQS...', 'yellow')}")
    try:
        ranked = search(instances)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)
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
        print(color_text(f"  Note: {instances_fallback_note}", 'yellow'))
        used_instances = 1
        try:
            ranked = search(1)
        except ValueError as e2:
            print(color_text(f"Error: {e2}", 'red'))
            sys.exit(1)

    sqs_structure = ranked[0]["structure"]
    objective_function = ranked[0]["objective_function"]

    species_meta = dict(structure.species_meta)
    for symbol in composition:
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(sqs_structure, species_meta=species_meta)

    output_atoms = len(new_structure.atoms)
    print(f"\n  {color_text('Output formula:', 'cyan')} {sqs_structure.composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {output_atoms}")
    print(f"  {color_text('Objective function:', 'cyan')} {objective_function:.6f}")

    structure_io.write_fdf(new_structure, args.output)

    report_file = "sqs_report.txt"
    report_text = generate_txt_report(
        args, input_formula, input_atoms, site_indices, composition,
        scaling, scaling_info, instances, used_instances, instances_fallback_note,
        cluster_size_and_shell, sqs_structure, objective_function, output_atoms,
    )
    with open(report_file, "w") as f:
        f.write(report_text)

    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")
    print(f"{color_text('[Saved]', 'cyan')} Report    -> {report_file}")


if __name__ == "__main__":
    main()
