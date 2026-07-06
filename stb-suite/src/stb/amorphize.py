#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.16.0"

import sys
import os
import time
import argparse
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core.cli import COLORS, color_text, show_intro
from stb.core.deps import require_mace

require_mace()


def bond_angle_stats(atoms, cutoff=None):
    """Mean/std (degrees) of nearest-neighbor bond angles -- a cheap,
    concrete diagnostic for whether a structure looks amorphous: a large
    std increase relative to the crystalline input signals loss of
    long-range order (crystals give a std near 0; a melt-quenched network
    gives a broad spread while the mean often stays near the crystal's
    short-range bond angle).
    """
    from ase.geometry import get_distances

    positions = atoms.get_positions()
    cell = atoms.get_cell()
    if cutoff is None:
        _, all_dists = get_distances(positions, cell=cell, pbc=True)
        np.fill_diagonal(all_dists, np.inf)
        cutoff = 1.25 * all_dists.min()

    vectors, lengths = get_distances(positions, cell=cell, pbc=True)
    n = len(atoms)
    angles = []
    for i in range(n):
        neighbors = [j for j in range(n) if j != i and lengths[i, j] < cutoff]
        for a in range(len(neighbors)):
            for b in range(a + 1, len(neighbors)):
                j, k = neighbors[a], neighbors[b]
                v1, v2 = vectors[i, j], vectors[i, k]
                cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                angles.append(np.degrees(np.arccos(np.clip(cos_t, -1, 1))))
    angles = np.array(angles) if angles else np.array([0.0])
    return angles.mean(), angles.std()


def print_progress(stage, step, total, temperature, start_time):
    elapsed = time.monotonic() - start_time
    line = f"  {stage}: step {step}/{total}, T={temperature:6.0f} K, {elapsed:6.1f}s elapsed"
    if sys.stderr.isatty():
        sys.stderr.write(f"\r{line}  ")
        sys.stderr.flush()
    elif step == total or step % max(1, total // 10) == 0:
        sys.stderr.write(line + "\n")
        sys.stderr.flush()


def finish_progress_line():
    if sys.stderr.isatty():
        sys.stderr.write("\r" + " " * 90 + "\r")
        sys.stderr.flush()


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Melt-quench amorphous structure generator using the MACE-MP-0 potential.", 'bold')}
Heats a crystalline structure above its melting point long enough to erase
crystalline memory, then cools it back down -- a fast heuristic starting
guess for an amorphous/glassy structure, meant to give SIESTA a much better
starting point than an ad-hoc random-displacement structure. NOT a
substitute for a slower, production-quality quench (more steps, slower
cooling) or DFT verification. Bulk (3D periodic) structures only --
melting a vacuum-padded slab/wire/molecule is physically meaningless.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f bulk_supercell.fdf\n"
               "  %(prog)s -f bulk_supercell.fdf --melt-temp 4000 --quench-steps 2000\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input crystalline structure file (.fdf), "
                             "typically a supercell (e.g. from stb-supercell).")
    parser.add_argument("--melt-temp", type=float, default=3000.0,
                        help="Melt temperature in K (default: 3000.0).")
    parser.add_argument("--melt-steps", type=int, default=500,
                        help="MD steps held at --melt-temp (default: 500).")
    parser.add_argument("--quench-temp", type=float, default=300.0,
                        help="Final target temperature in K (default: 300.0).")
    parser.add_argument("--quench-steps", type=int, default=1000,
                        help="MD steps for a linear cooling ramp from --melt-temp "
                             "to --quench-temp (default: 1000).")
    parser.add_argument("--timestep", type=float, default=1.0,
                        help="MD timestep in fs (default: 1.0).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size: speed/accuracy tradeoff (default: small).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Seed for the initial Maxwell-Boltzmann velocities "
                             "(default: not reproducible).")
    parser.add_argument("--no-final-relax", dest="final_relax", action="store_false",
                        help="Skip the final static (position+cell) relax. On by "
                             "default: the MD quench alone leaves real residual "
                             "thermal energy, not a clean local minimum.")
    parser.add_argument("-o", "--output", type=str, default="amorphous.fdf",
                        help="Output .fdf file name (default: amorphous.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-amorphize {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Melt-quench amorphous structure generator (MACE-MP-0):", 'bold'))
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

    frac_coords = [site.frac_coords for site in pmg]
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, pmg.lattice.matrix, 10.0)
    if any(vacuum_axes):
        print(color_text(
            "Error: a vacuum-padded axis was detected (e.g. a slab/wire/molecule). "
            "stb-amorphize only supports bulk (3D periodic) structures -- melting/"
            "NPT-relaxing a vacuum-padded structure is physically meaningless.", 'red'))
        sys.exit(1)

    mean0, std0 = bond_angle_stats(AseAtomsAdaptor.get_atoms(pmg))
    print(f"  {color_text('Input bond-angle mean/std:', 'cyan')} {mean0:.2f} / {std0:.2f} deg")

    print(f"  {color_text('Protocol:', 'cyan')} melt {args.melt_temp:.0f} K for {args.melt_steps} steps, "
          f"then ramp to {args.quench_temp:.0f} K over {args.quench_steps} steps "
          f"(timestep {args.timestep} fs)")

    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen

    atoms = AseAtomsAdaptor.get_atoms(pmg)
    atoms.calc = mace_relax.get_calculator(model=args.model, device=args.device, dtype="float32")

    MaxwellBoltzmannDistribution(atoms, temperature_K=min(args.melt_temp, 300.0), rng=(
        np.random.default_rng(args.seed) if args.seed is not None else None))

    dyn = NPTBerendsen(atoms, timestep=args.timestep * 1.0, temperature_K=args.melt_temp,
                       pressure_au=0.0, taut=50.0, taup=200.0, compressibility_au=4.57e-5)

    print(f"\n  {color_text('Melting...', 'yellow')}")
    start = time.monotonic()
    for step in range(1, args.melt_steps + 1):
        dyn.run(1)
        print_progress("Melt", step, args.melt_steps, atoms.get_temperature(), start)
    finish_progress_line()

    print(f"\n  {color_text('Quenching...', 'yellow')}")
    start = time.monotonic()
    for step in range(1, args.quench_steps + 1):
        frac = step / args.quench_steps
        target_t = args.melt_temp + (args.quench_temp - args.melt_temp) * frac
        dyn.set_temperature(temperature_K=target_t)
        dyn.run(1)
        print_progress("Quench", step, args.quench_steps, atoms.get_temperature(), start)
    finish_progress_line()

    print(f"\n  {color_text('After MD:', 'cyan')} T = {atoms.get_temperature():.0f} K")

    mean1, std1 = bond_angle_stats(atoms)
    print(f"  {color_text('Post-quench bond-angle mean/std:', 'cyan')} {mean1:.2f} / {std1:.2f} deg "
          f"(input was {mean0:.2f} / {std0:.2f} deg)")
    if std1 > std0:
        print(color_text(
            "  Bond-angle spread increased relative to the crystalline input -- "
            "consistent with a loss of long-range order.", 'green'))
    else:
        print(color_text(
            "  Bond-angle spread did NOT increase -- the structure may not have "
            "amorphized. Consider a higher --melt-temp or more --melt-steps.", 'yellow'))

    if args.final_relax:
        print(f"\n  {color_text('Final static relax...', 'yellow')} (float64, position+cell)")
        relax_calc = mace_relax.get_calculator(model=args.model, device=args.device, dtype="float64")
        converged, steps_used = mace_relax.relax(atoms, relax_calc, cell_mask=None)
        print(f"  {color_text('Steps used:', 'cyan')} {steps_used} "
              f"({color_text('converged', 'green') if converged else color_text('hit step cap', 'yellow')})")

    new_pmg = AseAtomsAdaptor.get_structure(atoms)
    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")
    print(color_text(
        "Note: a fast melt-quench heuristic, not a substitute for a slower "
        "production-quality quench or DFT verification.", 'yellow'))


if __name__ == "__main__":
    main()
