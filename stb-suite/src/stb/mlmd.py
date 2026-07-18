#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace

require_mace()

REPORT_FILE = "stb_mlmd_report.txt"

# Same 3 formats/extensions stb-ani2traj already writes a SIESTA AIMD
# trajectory into -- reused here so an MLMD trajectory opens in the same
# viewers (OVITO/VMD) the same way, no separate format story for the suite.
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}


def build_dynamics(atoms, ensemble, timestep_fs, temperature_k, friction_fs, pressure_gpa, seed,
                   taut_fs=100.0, taup_fs=500.0, compressibility_au=4.57e-5):
    """Returns a configured ASE MD driver, already given Maxwell-Boltzmann
    initial velocities (net momentum removed) at `temperature_k`. Timestep/
    friction/pressure are given in physically intuitive units (fs, fs^-1,
    GPa) and converted to ASE's internal unit system here -- NOT done by
    simply passing the raw numbers through (a real, verified bug in
    stb-amorphize's own MD call: passing timestep in fs directly as ASE
    time units silently runs at timestep/units.fs =~ 10.18x the intended
    step size; fixed there alongside adding this tool, see amorphize.py's
    own comment).

    `compressibility_au` (eV/Ang^3, NPT only) defaults to water's value (the
    same generic placeholder stb-amorphize uses) -- NOT material-specific;
    --compressibility lets a caller override it with something closer to the
    real material's 1/bulk-modulus for a better-calibrated barostat response.
    """
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

    rng = np.random.default_rng(seed) if seed is not None else None
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_k, rng=rng)
    Stationary(atoms)

    dt = timestep_fs * units.fs
    if ensemble == "nve":
        from ase.md.verlet import VelocityVerlet
        return VelocityVerlet(atoms, timestep=dt)
    if ensemble == "nvt":
        from ase.md.langevin import Langevin
        return Langevin(atoms, timestep=dt, temperature_K=temperature_k,
                         friction=friction_fs / units.fs, rng=rng)
    from ase.md.nptberendsen import NPTBerendsen
    return NPTBerendsen(atoms, timestep=dt, temperature_K=temperature_k,
                        pressure_au=pressure_gpa * units.GPa,
                        taut=taut_fs * units.fs, taup=taup_fs * units.fs,
                        compressibility_au=compressibility_au)


def main():
    parser = argparse.ArgumentParser(
        description="Molecular dynamics driven by a MACE potential -- the MACE-MP-0 "
                     "foundation model, or a custom model fine-tuned on your own SIESTA "
                     "data via stb-mlffAnalysis (--custom-model). NVE (microcanonical), NVT "
                     "(Langevin thermostat), or NPT (Berendsen thermostat+barostat, bulk-only) "
                     "ensembles. A fast heuristic exploration tool, not a substitute for a "
                     "real ab initio MD run.",
        epilog="Example usage:\n"
               "  stb-mlmd --file structure.fdf --ensemble nvt --temperature 300 --n-steps 5000\n"
               "  stb-mlmd --file structure.fdf --ensemble nvt --equilibration-steps 200 --n-steps 5000\n"
               "  stb-mlmd --file structure.fdf --custom-model mlff_model.model --ensemble npt \\\n"
               "                           --temperature 500 --pressure 0.0 --n-steps 10000",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--ensemble", choices=["nve", "nvt", "npt"], default="nvt",
                        help="MD ensemble (default: nvt). npt requires a bulk (no vacuum-padded "
                             "axis) structure -- cell dynamics on a slab/wire/molecule would be "
                             "physically meaningless.")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="Target temperature, K (default: 300).")
    parser.add_argument("--pressure", type=float, default=0.0,
                        help="Target pressure, GPa (npt only, default: 0.0 -- ambient).")
    parser.add_argument("--timestep", type=float, default=1.0, help="MD timestep, fs (default: 1.0).")
    parser.add_argument("--n-steps", type=int, default=1000,
                        help="Number of PRODUCTION MD steps, after --equilibration-steps (default: 1000).")
    parser.add_argument("--equilibration-steps", type=int, default=0,
                        help="Steps to run and discard BEFORE production (default: 0). Verified live "
                             "that a freshly Maxwell-Boltzmann-seeded structure relaxes for ~100-200 "
                             "steps before settling -- letting that transient into the saved "
                             "trajectory/statistics biases them towards the initial condition.")
    parser.add_argument("--friction", type=float, default=0.01,
                        help="Langevin friction coefficient, fs^-1 (nvt only, default: 0.01).")
    parser.add_argument("--taut", type=float, default=100.0,
                        help="Thermostat coupling time, fs (nvt/npt, default: 100).")
    parser.add_argument("--taup", type=float, default=500.0,
                        help="Barostat coupling time, fs (npt only, default: 500).")
    parser.add_argument("--compressibility", type=float, default=4.57e-5,
                        help="Compressibility, eV/Ang^3 (npt only, default: 4.57e-5 -- water's "
                             "value, a generic placeholder; override with something closer to "
                             "1/bulk-modulus for the actual material for a better-calibrated "
                             "barostat response).")
    parser.add_argument("--stride", type=int, default=10,
                        help="Save one trajectory frame every Nth MD step (default: 10).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0, "
                             "matches stb-kgrid/stb-mlrelax).")
    parser.add_argument("-of", "--out-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Trajectory format (default: xsf). Same 3 choices as stb-ani2traj.")
    parser.add_argument("-o", "--output", default=None,
                        help="Trajectory filename. Default: <input stem>_md_traj.<ext>")
    parser.add_argument("--output-structure", default="md_final.fdf",
                        help="Final-frame structure output (default: md_final.fdf).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for initial velocities.")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw energy/temperature-vs-time data (.dat) behind "
                             "the diagnostics plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlmd {VERSION}")

    args = parser.parse_args()

    ase_format, default_ext = OUTPUT_FORMATS[args.out_format]
    stem = os.path.splitext(os.path.basename(args.file))[0]
    traj_path = args.output if args.output else f"{stem}_md_traj{default_ext}"

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Molecular Dynamics",
            "MD driven by MACE-MP-0 or your own fine-tuned potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML MOLECULAR DYNAMICS:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLMD REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Ensemble          : {args.ensemble.upper()}", f_out)
    print_dual(f"Temperature       : {args.temperature:.1f} K", f_out)
    if args.ensemble == "npt":
        print_dual(f"Pressure          : {args.pressure:.3f} GPa", f_out)
    print_dual(f"Timestep / steps  : {args.timestep} fs / {args.n_steps} "
               f"(+ {args.equilibration_steps} equilibration)", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURE", f_out)

    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")

    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    print_dual(f"[OK] Read {len(structure.atoms)} atom(s), "
               f"dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)

    if args.ensemble == "npt" and any(vacuum_axes):
        fail("--ensemble npt requires a bulk (fully periodic) structure -- a vacuum-padded "
             "axis was detected. Cell dynamics on a slab/wire/molecule would be physically "
             "meaningless. Use --ensemble nve/nvt instead.")

    pmg_structure = structure_io.to_pymatgen(structure)
    atoms = AseAtomsAdaptor.get_atoms(pmg_structure)

    model_arg = args.custom_model if args.custom_model else args.model
    atoms.calc = mace_relax.get_calculator(model=model_arg, device=args.device, dtype="float32")

    dyn = build_dynamics(atoms, args.ensemble, args.timestep, args.temperature,
                         args.friction, args.pressure, args.seed,
                         taut_fs=args.taut, taup_fs=args.taup, compressibility_au=args.compressibility)

    if args.equilibration_steps > 0:
        print_section("[2] EQUILIBRATION (discarded)", f_out)
        print_dual(f"Running {args.equilibration_steps} step(s) before production starts -- "
                   "a freshly Maxwell-Boltzmann-seeded structure needs to relax before its "
                   "trajectory/statistics are representative (verified live: ~100-200 steps "
                   "of E_pot relaxation from the initial condition in a 64-atom test case).", f_out)
        dyn.run(args.equilibration_steps)
        print_dual(f"[OK] Equilibration done. E = {atoms.get_potential_energy():.4f} eV, "
                   f"T = {atoms.get_temperature():.1f} K", f_out)

    print_section("[3] RUNNING MD (production)", f_out)

    e0 = atoms.get_potential_energy()
    print_dual(f"Initial state     : E = {e0:.4f} eV, T = {atoms.get_temperature():.1f} K", f_out)

    frames = []
    times, epots, ekins, temps = [], [], [], []
    volumes = [] if args.ensemble == "npt" else None
    n_frames_expected = args.n_steps // args.stride
    progress_every = max(n_frames_expected // 10, 1)
    for step in range(1, args.n_steps + 1):
        dyn.run(1)
        if step % args.stride == 0:
            epot = atoms.get_potential_energy()
            ekin = atoms.get_kinetic_energy()
            temp = atoms.get_temperature()

            times.append(step * args.timestep)
            epots.append(epot)
            ekins.append(ekin)
            temps.append(temp)
            if volumes is not None:
                volumes.append(atoms.get_volume())

            frame = atoms.copy()
            if args.out_format == "xyz":
                frame.info['Epot'] = epot
                frame.info['Temp'] = temp
                frame.info['Time'] = step * args.timestep
            frames.append(frame)
            n_saved = len(frames)
            if n_saved % progress_every == 0 or step == args.n_steps:
                print_dual(f"  step {step:>7}/{args.n_steps}  "
                          f"T = {temp:7.1f} K  E = {epot:.4f} eV", f_out)

    times = np.array(times)
    epots = np.array(epots)
    ekins = np.array(ekins)
    temps = np.array(temps)
    etots = epots + ekins

    ef = atoms.get_potential_energy()
    print_dual(f"Final state       : E = {ef:.4f} eV, T = {atoms.get_temperature():.1f} K", f_out)
    print_dual(f"[OK] Collected {len(frames)} frame(s) (stride {args.stride}).", f_out)

    print_section("[4] ENERGY / TEMPERATURE DIAGNOSTICS", f_out)
    if args.ensemble == "nve":
        drift = (etots[-1] - etots[0]) / abs(etots[0]) if etots[0] != 0 else float('nan')
        print_dual(f"Total energy      : E_tot(0) = {etots[0]:.6f} eV, E_tot(end) = {etots[-1]:.6f} eV", f_out)
        print_dual(f"Relative drift    : {drift * 100:.4f}% (E_pot+E_kin -- the real NVE-conserved "
                   f"quantity, not E_pot alone, which fluctuates by design as it trades with E_kin)", f_out)
        if abs(drift) > 0.01:
            print_dual(color_text(
                "[WARNING] Total energy drifted by more than 1% -- consider a smaller --timestep.",
                'yellow'), f_out)
    else:
        print_dual(f"Temperature       : mean = {temps.mean():.1f} K, std = {temps.std():.1f} K "
                   f"(target {args.temperature:.1f} K)", f_out)
        if volumes is not None:
            vol_arr = np.array(volumes)
            print_dual(f"Cell volume       : mean = {vol_arr.mean():.3f} Ang^3, "
                       f"std = {vol_arr.std():.3f} Ang^3", f_out)

    diag_plot_path = f"{stem}_md_diagnostics.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(times, etots)
    axes[0].set_xlabel("Time (fs)")
    axes[0].set_ylabel("E_total = E_pot + E_kin (eV)")
    axes[0].set_title("Total energy")
    axes[1].plot(times, temps)
    # NVE has no thermostat, so args.temperature is only the initial-velocity
    # seed, not something the dynamics maintains -- label it accordingly
    # instead of "target" (which would only be true for nvt/npt).
    ref_label = "initial (NVE has no target)" if args.ensemble == "nve" else "target"
    axes[1].axhline(args.temperature, color='gray', linestyle='--', linewidth=0.8, label=ref_label)
    axes[1].set_xlabel("Time (fs)")
    axes[1].set_ylabel("Temperature (K)")
    axes[1].set_title("Temperature")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(diag_plot_path, dpi=150)
    plt.close(fig)
    print_dual(f"[OK] Saved diagnostics plot: {diag_plot_path}", f_out)

    diag_data_path = None
    if args.save_data:
        diag_data_path = f"{stem}_md_diagnostics.dat"
        np.savetxt(diag_data_path, np.column_stack([times, epots, ekins, etots, temps]),
                  header="time(fs) Epot(eV) Ekin(eV) Etotal(eV) Temp(K)")
        print_dual(f"[OK] Saved {diag_data_path}", f_out)

    print_section("[5] WRITING OUTPUT", f_out)
    from ase.io import write as ase_write
    try:
        ase_write(traj_path, frames, format=ase_format)
        print_dual(f"[OK] Saved trajectory: {traj_path} ({len(frames)} frame(s), {args.out_format} format)", f_out)
    except Exception as e:
        fail(f"Failed to write trajectory: {e}")

    new_pmg = AseAtomsAdaptor.get_structure(atoms)
    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output_structure)
    print_dual(f"[OK] Saved final-frame structure: {args.output_structure}", f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Trajectory        : {traj_path}", f_out)
    print_dual(f"Final structure   : {args.output_structure}", f_out)
    print_dual(f"Diagnostics plot  : {diag_plot_path}", f_out)
    if diag_data_path:
        print_dual(f"Diagnostics data  : {diag_data_path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("MD run complete -- a fast heuristic, not a substitute for ab initio MD.\n", 'bold'))


if __name__ == "__main__":
    main()
