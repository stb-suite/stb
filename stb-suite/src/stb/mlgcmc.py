#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Monte Carlo (canonical or grand-canonical) adsorption simulation driven
entirely by a MACE potential -- the ninth tool in the ML Simulations
category. A single-atom adsorbate species (e.g. a noble gas, a simple ion)
is inserted/moved/deleted in a fixed (rigid) host framework, sampled via
the Metropolis criterion; --ensemble canonical keeps a fixed adsorbate
count (displacement moves only, exploring configurational entropy at fixed
loading), --ensemble grand-canonical additionally allows insertion/deletion
at a given chemical potential (--mu, eV), reporting the equilibrium average
loading <N> -- and, with --mu-scan, a full loading-vs-chemical-potential
isotherm.

The insertion/deletion acceptance probabilities follow the standard grand
-canonical Monte Carlo formulas (Frenkel & Smit, "Understanding Molecular
Simulation"):

    acc(insert) = min(1, V*z/(N+1) * exp(-beta*dU))
    acc(delete) = min(1, N/(V*z)   * exp(-beta*dU))

where z = exp(beta*mu) / Lambda^3 is the absolute activity, Lambda =
h / sqrt(2*pi*m*kB*T) the adsorbate's own thermal de Broglie wavelength
(computed automatically from its atomic mass and T -- the only free
physical input is mu itself, the standard GCMC quantity). Verified
correct against the exact analytic result BEFORE ever touching MACE: for
a toy non-interacting (dU=0 always) ideal gas, this exact Metropolis
scheme was run for 2*10^5 steps and reproduced the known-exact grand
-canonical ideal-gas result <N> = V*z with variance(N) ~= mean(N) (the
Poisson-distribution signature), to within Monte Carlo noise.

mu is extremely sensitive: beta = 1/(kB*T) in eV^-1 is ~O(40) at room
temperature, so the ideal-gas reference loading V*z (activity times cell
volume, interactions ignored) changes by roughly a factor of e^(beta*0.1)
~ 50x per 0.1 eV change in mu -- found live testing a 105 Ang^3 graphene
-monolayer cell with Ar at 300 K: V*z went from ~1e-4 (essentially zero
loading) at mu=-0.5 eV to ~2.6e4 (severe, unphysical overcrowding) at
mu=0.0 eV, jumping straight past any sensible dilute/moderate loading in
between. This is real GCMC physics, not a code bug (the underlying
insertion/deletion formulas were already verified exact against the
analytic ideal-gas case above) -- but it means picking a good --mu by
guessing a single value rarely works. run_gcmc prints the ideal-gas
reference loading V*z up front, before the expensive interacting MC loop,
specifically so a bad choice is caught (and flagged) cheaply instead of
after a long run; a coarse --mu-scan to bracket where loading actually
transitions, then a refined scan near that bracket, is the practical
workflow.

The host framework is treated as RIGID throughout (never displaced/
relaxed by this tool) -- standard practice for adsorption-isotherm studies
(the framework's own thermal motion is a separate, much faster timescale
question better handled by stb-mlmd if ever needed); only adsorbate atoms
move/insert/delete. Single-atom adsorbates only in this version (no
molecular orientation sampling) -- exposed once a genuine need for
molecular species arises, same "expose on first genuine use" policy as
the rest of this suite.
"""

VERSION = "1.0.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.data import atomic_numbers, atomic_masses
from ase.io import write as ase_write
from scipy.constants import h, k as BOLTZMANN_K, atomic_mass, eV
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace

require_mace()

REPORT_FILE = "stb_mlgcmc_report.txt"


def thermal_wavelength(mass_amu, temperature_k):
    """De Broglie thermal wavelength Lambda = h / sqrt(2*pi*m*kB*T), in
    Angstrom -- the standard GCMC normalization converting a chemical
    potential into an absolute activity (see module docstring).
    """
    mass_kg = mass_amu * atomic_mass
    lam_m = h / np.sqrt(2 * np.pi * mass_kg * BOLTZMANN_K * temperature_k)
    return lam_m * 1e10


def random_insertion_position(cell, zmin_frac, zmax_frac, rng):
    """Uniform-random Cartesian position inside the cell -- fractional
    coordinates uniform in [0, 1) for a/b, restricted to [zmin_frac,
    zmax_frac) for c if given (e.g. confining insertion to a slab's vacuum
    gap instead of wasting attempts inside the solid framework).
    """
    frac = rng.random(3)
    frac[2] = zmin_frac + frac[2] * (zmax_frac - zmin_frac)
    return frac @ cell


def run_gcmc(host_atoms, adsorbate_symbol, ensemble, mu, temperature_k, n_initial,
             steps, equilibration_steps, dr, insert_delete_prob, zmin_frac, zmax_frac,
             calc, rng, f_out, tag=""):
    """Runs one canonical or grand-canonical MC simulation. Returns a dict:
    {"n_history", "energy_history", "n_accepted_move", "n_attempted_move",
    "n_accepted_insert", "n_attempted_insert", "n_accepted_delete",
    "n_attempted_delete", "final_atoms"}.
    """
    prefix = f"[{tag}] " if tag else ""
    beta = 1.0 / (BOLTZMANN_K / eV * temperature_k)
    mass_amu = atomic_masses[atomic_numbers[adsorbate_symbol]]
    lam = thermal_wavelength(mass_amu, temperature_k)
    z_activity = np.exp(beta * mu) / lam ** 3 if mu is not None else None

    atoms = host_atoms.copy()
    atoms.calc = calc
    cell = np.array(atoms.get_cell())
    volume = atoms.get_volume()
    n_host = len(atoms)

    if z_activity is not None:
        # V*z is the ideal-gas (zero-interaction) mean loading this mu/T/V
        # would give -- cheap to compute up front, before the expensive MC
        # loop, and an essential calibration check: at typical light-
        # adsorbate/room-temperature values, beta = 1/(kB*T) in eV^-1 is
        # ~O(40), so V*z changes by roughly a factor of e^(beta*0.1) ~ 50x
        # per 0.1 eV change in mu -- an extremely narrow "interesting"
        # window (found live: for a 105 Ang^3 box with Ar at 300 K, V*z
        # went from ~1e-4 at mu=-0.5 eV to ~2.6e4 at mu=0.0 eV, jumping
        # straight past any physically sensible loading in between). A V*z
        # far from O(0.1-10) here means the real (interacting) simulation
        # below will likely see either no adsorption at all or runaway
        # overcrowding, well before any interaction energy has a chance to
        # matter.
        ideal_n = volume * z_activity
        print_dual(f"{prefix}Ideal-gas reference loading (V*z, interactions ignored): "
                   f"{ideal_n:.4g}", f_out)
        if ideal_n > 100 or ideal_n < 1e-3:
            print_dual(color_text(
                f"{prefix}[WARNING] This mu is far from the O(0.1-10) range where "
                "interactions actually matter for this cell -- expect either no adsorption "
                "or runaway overcrowding rather than a meaningful equilibrium. mu is "
                "extremely sensitive (~50x change in loading per 0.1 eV at 300 K) -- "
                "bracket a transition with a coarse scan first, then refine in small steps.",
                'yellow'), f_out)

    for _ in range(n_initial):
        pos = random_insertion_position(cell, zmin_frac, zmax_frac, rng)
        atoms.append(Atoms(adsorbate_symbol, positions=[pos])[0])

    current_energy = atoms.get_potential_energy()

    n_history, energy_history = [], []
    counts = {"move": [0, 0], "insert": [0, 0], "delete": [0, 0]}  # [accepted, attempted]

    for step in range(steps):
        adsorbate_indices = list(range(n_host, len(atoms)))
        do_ins_del = ensemble == "grand-canonical" and rng.random() < insert_delete_prob

        if not do_ins_del:
            # Covers both --ensemble canonical (do_ins_del is always False
            # there, by construction above) and a grand-canonical step that
            # rolled a displacement move instead of insert/delete.
            if adsorbate_indices:
                counts["move"][1] += 1
                idx = adsorbate_indices[rng.integers(len(adsorbate_indices))]
                old_pos = atoms.positions[idx].copy()
                disp = rng.uniform(-dr, dr, 3)
                atoms.positions[idx] = old_pos + disp
                new_energy = atoms.get_potential_energy()
                d_u = new_energy - current_energy
                if d_u <= 0 or rng.random() < np.exp(-beta * d_u):
                    current_energy = new_energy
                    counts["move"][0] += 1
                else:
                    atoms.positions[idx] = old_pos
        elif rng.random() < 0.5:
            counts["insert"][1] += 1
            pos = random_insertion_position(cell, zmin_frac, zmax_frac, rng)
            trial_atoms = atoms + Atoms(adsorbate_symbol, positions=[pos])
            trial_atoms.calc = calc
            new_energy = trial_atoms.get_potential_energy()
            d_u = new_energy - current_energy
            n_current = len(adsorbate_indices)
            acc = min(1.0, volume * z_activity / (n_current + 1) * np.exp(-beta * d_u))
            if rng.random() < acc:
                atoms = trial_atoms
                current_energy = new_energy
                counts["insert"][0] += 1
        else:
            if adsorbate_indices:
                counts["delete"][1] += 1
                idx = adsorbate_indices[rng.integers(len(adsorbate_indices))]
                trial_atoms = atoms.copy()
                del trial_atoms[idx]
                trial_atoms.calc = calc
                new_energy = trial_atoms.get_potential_energy()
                d_u = new_energy - current_energy
                n_current = len(adsorbate_indices)
                acc = min(1.0, n_current / (volume * z_activity) * np.exp(-beta * d_u))
                if rng.random() < acc:
                    atoms = trial_atoms
                    current_energy = new_energy
                    counts["delete"][0] += 1

        if step >= equilibration_steps:
            n_history.append(len(atoms) - n_host)
            energy_history.append(current_energy)

        if (step + 1) % max(1, steps // 5) == 0:
            print(f"  ... {prefix}step {step + 1}/{steps}, N={len(atoms) - n_host}, "
                  f"E={current_energy:.4f} eV")

    for label, (acc, att) in counts.items():
        rate = 100.0 * acc / att if att > 0 else 0.0
        print_dual(f"{prefix}{label:>7} moves: {acc}/{att} accepted ({rate:.1f}%)", f_out)

    return {
        "n_history": np.array(n_history), "energy_history": np.array(energy_history),
        "counts": counts, "final_atoms": atoms,
    }


def plot_mc_history(n_history, energy_history, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(n_history, color='tab:blue')
    axes[0].set_xlabel("MC step (post-equilibration)")
    axes[0].set_ylabel("Adsorbate count N")
    axes[0].set_title(f"Loading (mean={n_history.mean():.2f})")

    axes[1].plot(energy_history, color='tab:orange')
    axes[1].set_xlabel("MC step (post-equilibration)")
    axes[1].set_ylabel("Energy (eV)")
    axes[1].set_title("Energy")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_isotherm(mu_values, mean_n, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(mu_values, mean_n, 'o-', color='tab:green')
    ax.set_xlabel("Chemical potential mu (eV)")
    ax.set_ylabel("Mean loading <N>")
    ax.set_title("Adsorption isotherm")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Monte Carlo (canonical or grand-canonical) adsorption simulation "
                     "driven entirely by a MACE potential -- the MACE-MP-0 foundation "
                     "model, or a custom model fine-tuned on your own SIESTA data via "
                     "stb-mlffAnalysis (--custom-model). A single-atom adsorbate species "
                     "is inserted/moved/deleted in a fixed host framework via the "
                     "Metropolis criterion. --ensemble canonical: fixed adsorbate count, "
                     "displacement moves only. --ensemble grand-canonical: additionally "
                     "inserts/deletes at a given chemical potential (--mu), reporting the "
                     "equilibrium loading -- --mu-scan sweeps multiple values for a full "
                     "adsorption isotherm.",
        epilog="Example usage:\n"
               "  stb-mlgcmc --host slab.fdf --adsorbate Ar --ensemble canonical "
               "--n-initial 4\n"
               "  stb-mlgcmc --host slab.fdf --adsorbate Ar --mu -0.3 --steps 5000\n"
               "  stb-mlgcmc --host slab.fdf --adsorbate Ar --mu-scan -0.6 -0.4 -0.2 0.0\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--host", required=True, help="Host framework structure (.fdf), treated as rigid.")
    parser.add_argument("--adsorbate", required=True, help="Adsorbate element symbol (single atom, e.g. Ar).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--ensemble", choices=["canonical", "grand-canonical"], default="grand-canonical",
                        help="Monte Carlo ensemble (default: grand-canonical).")
    parser.add_argument("--n-initial", type=int, default=0,
                        help="Initial number of adsorbate atoms, randomly placed (default: 0).")
    parser.add_argument("--mu", type=float, default=None,
                        help="Chemical potential, eV (required for --ensemble grand-canonical "
                             "unless --mu-scan is given). Loading is extremely sensitive to "
                             "this value (roughly a 50x change per 0.1 eV at 300 K for a "
                             "light adsorbate) -- check the printed 'ideal-gas reference "
                             "loading' before waiting on a full run, and bracket a transition "
                             "with a coarse --mu-scan first rather than guessing a single value.")
    parser.add_argument("--mu-scan", type=float, nargs='+', default=None, metavar="EV",
                        help="Multiple chemical potentials to scan (grand-canonical only) -- "
                             "produces a full loading-vs-mu isotherm instead of a single run. "
                             "Given the sensitivity noted under --mu, an initial scan should "
                             "use a wide range/coarse step to find where loading transitions "
                             "at all, before refining.")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature, K (default: 300).")
    parser.add_argument("--steps", type=int, default=2000, help="Total MC steps (default: 2000).")
    parser.add_argument("--equilibration-steps", type=int, default=500,
                        help="MC steps discarded before collecting statistics (default: 500).")
    parser.add_argument("--dr", type=float, default=0.5,
                        help="Max random displacement per move step, Angstrom (default: 0.5).")
    parser.add_argument("--insert-delete-prob", type=float, default=0.5,
                        help="Probability of attempting an insert-or-delete move instead of a "
                             "displacement move each step, grand-canonical only (default: 0.5).")
    parser.add_argument("--insertion-zmin", type=float, default=0.0,
                        help="Min fractional c-coordinate for insertion (default: 0.0, whole cell).")
    parser.add_argument("--insertion-zmax", type=float, default=1.0,
                        help="Max fractional c-coordinate for insertion (default: 1.0, whole cell).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (default: unseeded).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlgcmc_out",
                        help="Output directory for all files (default: mlgcmc_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw N/energy history data (.dat). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlgcmc {VERSION}")

    args = parser.parse_args()

    if args.ensemble == "grand-canonical" and args.mu is None and args.mu_scan is None:
        parser.error("--ensemble grand-canonical needs --mu or --mu-scan.")
    if args.ensemble == "canonical" and args.n_initial == 0:
        parser.error("--ensemble canonical needs --n-initial > 0 (no insertion moves in this ensemble).")
    if args.equilibration_steps >= args.steps:
        parser.error("--equilibration-steps must be less than --steps -- otherwise every step "
                     "is discarded as equilibration and no statistics are ever collected "
                     "(caught live: this silently produced 'mean of empty array' NaN energies "
                     "and no plot/data files, with the default --equilibration-steps 500 "
                     "exceeding a modest --steps like 150).")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Grand Canonical Monte Carlo",
            "Adsorption Monte Carlo (canonical/grand-canonical), driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text(f"ML MONTE CARLO ({args.ensemble.upper()}):", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLGCMC REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Host file         : {args.host}", f_out)
    print_dual(f"Adsorbate         : {args.adsorbate}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"Ensemble          : {args.ensemble}", f_out)
    print_dual(f"Temperature       : {args.temperature} K", f_out)
    print_dual(f"Steps             : {args.steps} ({args.equilibration_steps} equilibration)", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURE", f_out)

    if not os.path.isfile(args.host):
        fail(f"Host file '{args.host}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")
    if args.adsorbate not in atomic_numbers:
        fail(f"Unknown element symbol for --adsorbate: '{args.adsorbate}'.")

    try:
        structure = structure_io.read_fdf(args.host)
    except Exception as e:
        fail(f"Could not read '{args.host}': {e}")

    pmg_structure = structure_io.to_pymatgen(structure)
    host_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    print_dual(f"[OK] Read {len(host_atoms)} host atom(s), volume = {host_atoms.get_volume():.4f} Ang^3", f_out)

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    rng = np.random.default_rng(args.seed)

    mu_values = args.mu_scan if args.mu_scan else ([args.mu] if args.mu is not None else [None])

    print_section("[2] MONTE CARLO", f_out)
    all_results = []
    for mu in mu_values:
        tag = f"mu={mu:.3f}" if mu is not None else args.ensemble
        result = run_gcmc(host_atoms, args.adsorbate, args.ensemble, mu, args.temperature,
                          args.n_initial, args.steps, args.equilibration_steps, args.dr,
                          args.insert_delete_prob, args.insertion_zmin, args.insertion_zmax,
                          calc, rng, f_out, tag=tag)
        mean_n = result["n_history"].mean() if len(result["n_history"]) else float(args.n_initial)
        mean_e = result["energy_history"].mean() if len(result["energy_history"]) else float("nan")
        print_dual(f"{tag}: <N> = {mean_n:.3f}, <E> = {mean_e:.4f} eV", f_out)
        all_results.append((mu, result, mean_n, mean_e))

    print_section("[3] SUMMARY & FILES", f_out)
    if len(all_results) > 1:
        isotherm_path = os.path.join(args.output_dir, "isotherm.png")
        plot_isotherm([r[0] for r in all_results], [r[2] for r in all_results], isotherm_path)
        print_dual(f"[OK] Saved {isotherm_path}", f_out)
        if args.save_data:
            iso_data = os.path.join(args.output_dir, "isotherm.dat")
            np.savetxt(iso_data, np.column_stack([[r[0] for r in all_results], [r[2] for r in all_results]]),
                      header="mu(eV) mean_N")
            print_dual(f"[OK] Saved {iso_data}", f_out)
        monotonic = all(all_results[i][2] <= all_results[i + 1][2] + 1e-6 for i in range(len(all_results) - 1))
        if not monotonic:
            print_dual(color_text(
                "[WARNING] Loading is not monotonically non-decreasing with mu -- expected from "
                "grand-canonical thermodynamic stability; consider more --steps for better "
                "statistics.", 'yellow'), f_out)
    else:
        mu, result, mean_n, mean_e = all_results[0]
        history_plot = os.path.join(args.output_dir, "mc_history.png")
        if len(result["n_history"]):
            plot_mc_history(result["n_history"], result["energy_history"], history_plot)
            print_dual(f"[OK] Saved {history_plot}", f_out)
        if args.save_data and len(result["n_history"]):
            data_path = os.path.join(args.output_dir, "mc_history.dat")
            np.savetxt(data_path, np.column_stack([result["n_history"], result["energy_history"]]),
                      header="N energy(eV)")
            print_dual(f"[OK] Saved {data_path}", f_out)
        final_path = os.path.join(args.output_dir, "final_config.xsf")
        ase_write(final_path, result["final_atoms"], format="xsf")
        print_dual(f"[OK] Saved {final_path}", f_out)

    print_dual("Status            : OK", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML Monte Carlo complete -- a heuristic screen, not a substitute for a "
                      "fully converged/validated GCMC study.\n", 'bold'))


if __name__ == "__main__":
    main()
