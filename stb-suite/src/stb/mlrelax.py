#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.15.0"

import sys
import os
import argparse
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.deps import require_mace

require_mace()

OPTIMIZERS = ["FIRE", "BFGS", "LBFGS"]


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Fast ML pre-relaxation with the MACE-MP-0 foundation potential.", 'bold')}
Cleans up obviously-wrong geometry (bad lattice guesses, freshly-added
defect/passivant atoms at the wrong bond length) cheaply, so a follow-up
SIESTA relaxation starts much closer to its own equilibrium and needs far
fewer ionic steps. This is a pre-relaxation heuristic, NOT a substitute for
a real DFT relaxation -- universal MLIPs are trained mostly on
Materials-Project-like inorganic crystals and are less reliable far from
that distribution.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f defect.fdf\n"
               "  %(prog)s -f bulk.fdf --relax-cell -o bulk_relaxed.fdf\n"
               "  %(prog)s -f slab2d.fdf --relax-cell -o slab2d_relaxed.fdf\n"
               "\n"
               "Note: an initial guess that's very far from equilibrium (e.g. a lattice "
               "constant off by 15-20%% or more) can occasionally lead a universal MLIP to "
               "a spurious low-energy state instead of the real minimum -- sanity-check a "
               "large --relax-cell change before trusting it.\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf).")
    parser.add_argument("--relax-cell", action="store_true",
                        help="Also relax the lattice, not just atomic positions. "
                             "Automatically adapts to the structure's periodicity: for a "
                             "slab/wire (vacuum-padded axes, e.g. from stb-slab), only the "
                             "genuinely periodic direction(s) relax and the vacuum "
                             "thickness is held exactly fixed -- relaxing stress along a "
                             "vacuum direction would be physically meaningless. A fully "
                             "isolated structure (vacuum in all 3 directions, e.g. a "
                             "molecule) has nothing to relax; --relax-cell is ignored.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes for --relax-cell "
                             "(default: 10.0, matches stb-kgrid).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size: speed/accuracy tradeoff (default: small). "
                             "Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned on your own SIESTA data via "
                             "stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold, eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Maximum optimizer steps (default: 200).")
    parser.add_argument("--optimizer", choices=list(OPTIMIZERS), default="FIRE",
                        help="ASE optimizer to use (default: FIRE).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output", type=str, default="relaxed.fdf",
                        help="Output .fdf file name (default: relaxed.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlrelax {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Fast ML pre-relaxation (MACE-MP-0) -- not a substitute for DFT:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    pmg = structure_io.to_pymatgen(structure)
    print(f"  {color_text('Input formula:', 'cyan')} {pmg.composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {len(pmg)}")

    cell_mask = None
    if args.relax_cell:
        frac_coords = [site.frac_coords for site in pmg]
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, pmg.lattice.matrix, args.vacuum_gap)
        axis_names = ("a", "b", "c")
        if all(vacuum_axes):
            print(color_text(
                "  Note: every axis is vacuum-padded (a fully isolated structure, e.g. a "
                "molecule) -- there is no periodic direction to relax. Ignoring --relax-cell, "
                "doing positions-only.", 'yellow'))
        else:
            cell_mask = mace_relax.build_cell_mask(vacuum_axes)
            fixed = [axis_names[i] for i, v in enumerate(vacuum_axes) if v]
            relaxed = [axis_names[i] for i, v in enumerate(vacuum_axes) if not v]
            if fixed:
                print(f"  {color_text('Vacuum axis/axes held fixed:', 'cyan')} {', '.join(fixed)} "
                      f"(gap >= {args.vacuum_gap} Ang detected)")
            print(f"  {color_text('Relaxing lattice along:', 'cyan')} {', '.join(relaxed)}")

    atoms = AseAtomsAdaptor.get_atoms(pmg)

    print(f"  {color_text('Mode:', 'cyan')} {'positions + cell' if cell_mask is not None else 'positions only'}")
    model_arg = args.custom_model if args.custom_model else args.model
    if args.custom_model and not os.path.isfile(args.custom_model):
        sys.exit(f"{color_text('[ERROR]', 'red')} --custom-model file not found: {args.custom_model}")
    print(f"  {color_text('Model:', 'cyan')} "
          f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}")

    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    atoms.calc = calc

    e0 = atoms.get_potential_energy()
    f0 = np.abs(atoms.get_forces()).max()
    print(f"\n  {color_text('Before:', 'cyan')} E = {e0:.4f} eV, max|F| = {f0:.4f} eV/Ang")

    positions_before = atoms.get_positions().copy()
    cell_before = atoms.get_cell().copy()

    print(f"\n  {color_text('Relaxing...', 'yellow')} (fmax={args.fmax}, max {args.max_steps} steps)")
    converged, steps_used = mace_relax.relax(
        atoms, calc, cell_mask=cell_mask, optimizer=args.optimizer,
        fmax=args.fmax, max_steps=args.max_steps)

    e1 = atoms.get_potential_energy()
    f1 = np.abs(atoms.get_forces()).max()
    print(f"\n  {color_text('After:', 'cyan')} E = {e1:.4f} eV, max|F| = {f1:.4f} eV/Ang")
    print(f"  {color_text('Steps used:', 'cyan')} {steps_used} "
          f"({color_text('converged', 'green') if converged else color_text('hit step cap, NOT converged', 'yellow')})")

    if cell_mask is not None:
        abc_before = cell_before.cellpar()[:3]
        abc_after = atoms.get_cell().cellpar()[:3]
        print(f"  {color_text('Lattice a,b,c before:', 'cyan')} "
              f"{abc_before[0]:.4f}, {abc_before[1]:.4f}, {abc_before[2]:.4f} Ang")
        print(f"  {color_text('Lattice a,b,c after: ', 'cyan')} "
              f"{abc_after[0]:.4f}, {abc_after[1]:.4f}, {abc_after[2]:.4f} Ang")
        rel_change = np.abs(abc_after - abc_before) / abc_before
        if rel_change.max() > 0.15:
            print(color_text(
                "  Warning: a relaxed lattice parameter changed by more than 15% -- "
                "universal MLIPs can occasionally settle into a spurious state this far "
                "from the initial guess. Sanity-check this result before trusting it.",
                'yellow'))
    else:
        max_disp = np.linalg.norm(atoms.get_positions() - positions_before, axis=1).max()
        print(f"  {color_text('Max atomic displacement:', 'cyan')} {max_disp:.4f} Ang")

    new_pmg = AseAtomsAdaptor.get_structure(atoms)
    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
