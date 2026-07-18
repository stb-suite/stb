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

from stb.core import structure_io, kspace, mace_relax, mace_phonons
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.phonons_pos import build_band_path, band_path_to_phonopy_format, band_tick_positions, pretty_label

require_mace()

REPORT_FILE = "stb_mlphonons_report.txt"


def plot_bands(bs, labels, path_connections, out_path):
    """Matplotlib band-structure plot (same visual convention as the rest
    of the suite's newer tools -- aimd_analysis.py/mlff_analysis.py/
    mlmd.py -- a PNG, not phonons_pos.py's own gnuplot .dat+.gplot pair).
    Reuses bs.frequencies/.distances (already computed by
    phonon.run_band_structure(), no recomputation) and
    band_tick_positions() for the high-symmetry tick labels.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    offset = 0.0
    for seg_dist, seg_freq in zip(bs.distances, bs.frequencies):
        x = seg_dist
        for band in seg_freq.T:
            ax.plot(x, band, color='tab:blue', linewidth=0.8)
        offset = x[-1]
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)

    ticks = band_tick_positions(bs.distances, labels, path_connections)
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([pretty_label(t[1]) for t in ticks])
    for t in ticks:
        ax.axvline(t[0], color='gray', linewidth=0.5, alpha=0.5)
    ax.set_xlim(bs.distances[0][0], bs.distances[-1][-1])
    ax.set_ylabel("Frequency (THz)")
    ax.set_title("Phonon band structure (MACE)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Standalone phonon calculation driven entirely by a MACE potential -- the "
                     "MACE-MP-0 foundation model, or a custom model fine-tuned on your own "
                     "SIESTA data via stb-mlffAnalysis (--custom-model). One shot: displacements, "
                     "force constants, band structure, DOS, and thermal properties, no SIESTA "
                     "and no separate analysis tool needed. A fast heuristic preview, not a "
                     "substitute for the real SIESTA-based stb-phononsCreate/stb-phononsPos "
                     "workflow -- ML potentials are markedly less reliable for phonons "
                     "(sign of low-frequency/imaginary modes especially) than DFT.",
        epilog="Example usage:\n"
               "  stb-mlphonons --file structure.fdf --dim 2 2 2\n"
               "  stb-mlphonons --file structure.fdf --custom-model mlff_model.model --dim 3 3 3",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--dim", type=int, nargs=3, default=[2, 2, 2], metavar=("NA", "NB", "NC"),
                        help="Supercell dimensions for the finite-displacement method (default: 2 2 2).")
    parser.add_argument("--distance", type=float, default=0.01,
                        help="Displacement distance, Angstrom (default: 0.01).")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the MACE pre-relax (positions only) before generating "
                             "displacements. On by default -- a phonon calculation assumes "
                             "~zero net force at the reference geometry; skipping this on an "
                             "unrelaxed structure is a common cause of spurious imaginary modes.")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the pre-relax, eV/Ang (default: 0.05).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0).")
    parser.add_argument("--band-points", type=int, default=101,
                        help="q-points per band-structure segment (default: 101).")
    parser.add_argument("--mesh", type=int, nargs=3, default=[20, 20, 20], metavar=("NA", "NB", "NC"),
                        help="q-point mesh for DOS/thermal properties (default: 20 20 20).")
    parser.add_argument("--tmin", type=float, default=0.0, help="Min temperature, K (default: 0).")
    parser.add_argument("--tmax", type=float, default=1000.0, help="Max temperature, K (default: 1000).")
    parser.add_argument("--tstep", type=float, default=10.0, help="Temperature step, K (default: 10).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlphonons_out",
                        help="Output directory for all plots/data (default: mlphonons_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw numeric data (.dat) behind each plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlphonons {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Phonons",
            "Bands, DOS, and thermal properties, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML PHONONS:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLPHONONS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
    print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
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

    pmg_structure = structure_io.to_pymatgen(structure)
    atoms = AseAtomsAdaptor.get_atoms(pmg_structure)

    lattice_ang = np.array(atoms.get_cell())
    frac_coords = atoms.get_scaled_positions()
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
    print_dual(f"[OK] Read {len(atoms)} atom(s), "
               f"dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)

    vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                              if is_vac and args.dim['abc'.index(axis)] > 1]
    if vacuum_dims_requested:
        print_dual(color_text(
            f"[WARNING] --dim requests more than 1 repetition along vacuum-padded "
            f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
            "Replicating a supercell across vacuum only multiplies the MACE cost "
            "without adding real periodicity -- consider --dim 1 on that axis.", 'yellow'), f_out)

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    atoms.calc = calc
    f0 = np.abs(atoms.get_forces()).max()
    print_dual(f"Residual force on input structure: {f0:.4f} eV/Ang", f_out)
    if not args.relax:
        print_dual(color_text(
            "[WARNING] --no-relax: generating displacements from the structure as given. "
            "If it isn't already at a MACE-relaxed minimum, expect spurious imaginary "
            "modes downstream.", 'yellow'), f_out)

    print_section("[2] DISPLACEMENTS & FORCE CONSTANTS", f_out)

    phonon, _, relax_info = mace_phonons.generate_ml_displacements(
        atoms, args.dim, args.distance, calc, relax=args.relax, fmax=args.fmax)
    if relax_info is not None:
        _, f1, converged, steps_used = relax_info
        print_dual(f"After relax: max|F| = {f1:.4f} eV/Ang "
                   f"({'converged' if converged else 'hit step cap, not fully converged'}, "
                   f"{steps_used} steps)", f_out)

    n_disp = len(phonon.supercells_with_displacements)
    print_dual(f"Displacements needed: {n_disp}", f_out)
    report_every = max(1, n_disp // 10)

    def _progress(i, n):
        if i % report_every == 0 or i == n:
            print(f"  ... {i}/{n}")

    mace_phonons.compute_force_constants(phonon, calc, progress_callback=_progress)
    yaml_path = os.path.join(args.output_dir, "phonopy_disp.yaml")
    phonon.save(yaml_path, settings={'force_constants': True})
    print_dual(f"[OK] Force constants computed. Saved {yaml_path}", f_out)

    print_section("[3] BAND STRUCTURE", f_out)
    bands_plot = None
    path = build_band_path(phonon.primitive, 1.0, args.vacuum_gap)
    if path is None:
        print_dual(color_text(
            "Skipped: every axis is vacuum-padded (0D/isolated system) -- a q-path "
            "isn't physically meaningful here.", 'yellow'), f_out)
    else:
        kpoints_dict, path_segments, bravais_name = path
        print_dual(f"Bravais lattice   : {bravais_name}", f_out)
        path_str = " | ".join("-".join(pretty_label(l) for l in seg) for seg in path_segments)
        print_dual(f"Path              : {path_str}", f_out)

        bands, band_labels, path_connections = band_path_to_phonopy_format(
            kpoints_dict, path_segments, args.band_points)
        phonon.run_band_structure(bands, with_group_velocities=False,
                                  path_connections=path_connections,
                                  labels=band_labels, is_legacy_plot=False)
        bs = phonon.band_structure

        band_min_freq = np.concatenate(bs.frequencies).min()
        print_dual(f"Minimum band-path frequency: {band_min_freq:.4f} THz", f_out)
        if band_min_freq < 0:
            print_dual(color_text(
                "[WARNING] Negative (imaginary) frequency found along the band path -- "
                "sign of an unstable/unconverged structure, or a limitation of the ML "
                "potential for this system. Try --no-relax off (default) with a tighter "
                "--fmax, or verify with the real SIESTA-based phonon workflow.", 'red'), f_out)
        else:
            print_dual(color_text("No imaginary modes found along the band path.", 'green'), f_out)

        bands_plot = os.path.join(args.output_dir, "bands.png")
        plot_bands(bs, band_labels, path_connections, bands_plot)
        print_dual(f"[OK] Saved {bands_plot}", f_out)
        if args.save_data:
            bands_data = os.path.join(args.output_dir, "bands.dat")
            with open(bands_data, "w") as f:
                f.write("# distance(1/Ang) freq_band1 freq_band2 ... (THz, blank line = segment break)\n")
                for seg_dist, seg_freq in zip(bs.distances, bs.frequencies):
                    for d, row in zip(seg_dist, seg_freq):
                        f.write(f"{d:.6f} " + " ".join(f"{f:.6f}" for f in row) + "\n")
                    f.write("\n")
            print_dual(f"[OK] Saved {bands_data}", f_out)

    print_section("[4] DENSITY OF STATES", f_out)
    phonon.run_mesh(args.mesh)
    phonon.run_total_dos()
    dos_dict = phonon.get_total_dos_dict()
    freqs, dos = dos_dict['frequency_points'], dos_dict['total_dos']

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(freqs, dos)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("DOS (states/THz/cell)")
    ax.set_title("Phonon density of states (MACE)")
    fig.tight_layout()
    dos_plot = os.path.join(args.output_dir, "dos.png")
    fig.savefig(dos_plot, dpi=150)
    plt.close(fig)
    print_dual(f"[OK] Saved {dos_plot}", f_out)
    if args.save_data:
        dos_data = os.path.join(args.output_dir, "dos.dat")
        np.savetxt(dos_data, np.column_stack([freqs, dos]), header="freq(THz) DOS")
        print_dual(f"[OK] Saved {dos_data}", f_out)

    print_section("[5] THERMAL PROPERTIES", f_out)
    phonon.run_thermal_properties(t_min=args.tmin, t_max=args.tmax, t_step=args.tstep)
    tp = phonon.get_thermal_properties_dict()
    temps = tp['temperatures']
    free_energy = tp['free_energy']
    entropy = tp['entropy']
    heat_capacity = tp['heat_capacity']

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].plot(temps, free_energy)
    axes[0].set_xlabel("Temperature (K)")
    axes[0].set_ylabel("Free energy (kJ/mol)")
    axes[0].set_title("Free energy")
    axes[1].plot(temps, entropy)
    axes[1].set_xlabel("Temperature (K)")
    axes[1].set_ylabel("Entropy (J/K/mol)")
    axes[1].set_title("Entropy")
    axes[2].plot(temps, heat_capacity)
    axes[2].set_xlabel("Temperature (K)")
    axes[2].set_ylabel("Heat capacity (J/K/mol)")
    axes[2].set_title("Heat capacity (Cv)")
    fig.tight_layout()
    thermal_plot = os.path.join(args.output_dir, "thermal.png")
    fig.savefig(thermal_plot, dpi=150)
    plt.close(fig)
    print_dual(f"[OK] Saved {thermal_plot}", f_out)
    if args.save_data:
        thermal_data = os.path.join(args.output_dir, "thermal.dat")
        np.savetxt(thermal_data, np.column_stack([temps, free_energy, entropy, heat_capacity]),
                  header="T(K) FreeEnergy(kJ/mol) Entropy(J/K/mol) HeatCapacity(J/K/mol)")
        print_dual(f"[OK] Saved {thermal_data}", f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Force constants   : {yaml_path}", f_out)
    if bands_plot:
        print_dual(f"Band structure    : {bands_plot}", f_out)
    print_dual(f"DOS               : {dos_plot}", f_out)
    print_dual(f"Thermal properties: {thermal_plot}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML phonons complete -- a fast heuristic, not a substitute for DFT phonons.\n", 'bold'))


if __name__ == "__main__":
    main()
