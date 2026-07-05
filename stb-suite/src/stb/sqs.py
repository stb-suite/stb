#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import sys
import os
import argparse
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.deps import require_icet

require_icet()
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


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Generates a Special Quasirandom Structure (SQS) for a substitutional alloy.", 'bold')}
Picks one sublattice (--sublattice) and disorders it with a target
composition, then searches for the atomic arrangement that best mimics an
ideal random alloy (via icet's Monte Carlo or enumeration SQS search).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4\n"
               "  %(prog)s -f ni.fdf --sublattice Ni --composition Ni:0.75,Fe:0.25 --scaling 8 --method enumeration\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf).")
    parser.add_argument("--sublattice", type=str, required=True,
                        help="Existing species whose sites become the disordered alloy sublattice.")
    parser.add_argument("--composition", type=str, required=True,
                        help="Target composition on that sublattice, e.g. 'Fe:0.5,Ni:0.5' (must sum to 1.0).")
    parser.add_argument("--scaling", type=int, required=True,
                        help="icet supercell-size control (positive integer). Multiplies icet's own "
                             "internally-detected primitive unit, not necessarily the literal input atom "
                             "count -- read the actual resulting atom count from the output.")
    parser.add_argument("--method", choices=["monte_carlo", "enumeration"], default="monte_carlo",
                        help="SQS search method (default: monte_carlo).")
    parser.add_argument("--instances", type=int, default=1,
                        help="Number of parallel search instances (default: 1).")
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

    if args.scaling < 1:
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
    print(f"  {color_text('Input formula:', 'cyan')} {pmg.composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {len(pmg)}")

    site_indices = [i for i, site in enumerate(pmg) if site.specie.symbol == args.sublattice]
    if not site_indices:
        print(color_text(f"Error: species '{args.sublattice}' not found in the structure.", 'red'))
        sys.exit(1)

    print(f"  {color_text('Sublattice:', 'cyan')} {args.sublattice} ({len(site_indices)} site(s))")
    print(f"  {color_text('Target composition:', 'cyan')} {composition}")
    print(f"  {color_text('Scaling:', 'cyan')} {args.scaling}  {color_text('Method:', 'cyan')} {args.method}")

    for idx in site_indices:
        pmg.replace(idx, composition)

    icet_sqs_kwargs = {}
    if args.mc_steps is not None:
        icet_sqs_kwargs["n_steps"] = args.mc_steps

    transformation = SQSTransformation(
        scaling=args.scaling,
        sqs_method=f"icet-{args.method}",
        instances=max(1, args.instances),
        temperature=args.temperature,
        cluster_size_and_shell=cluster_size_and_shell,
        icet_sqs_kwargs=icet_sqs_kwargs,
    )

    print(f"\n  {color_text('Searching for SQS...', 'yellow')}")
    try:
        ranked = transformation.apply_transformation(pmg, return_ranked_list=1)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    sqs_structure = ranked[0]["structure"]
    objective_function = ranked[0]["objective_function"]

    species_meta = dict(structure.species_meta)
    for symbol in composition:
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(sqs_structure, species_meta=species_meta)

    print(f"\n  {color_text('Output formula:', 'cyan')} {sqs_structure.composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {len(new_structure.atoms)}")
    print(f"  {color_text('Objective function:', 'cyan')} {objective_function:.6f}")

    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
