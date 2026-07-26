#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.1.0"

import sys
import os
import time
import argparse
from datetime import datetime

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_mace

REPORT_FILE = "stb_amorphize_report.txt"
BIB_FILE = "references.bib"

# Same 3-format trajectory choice as stb-ani2traj/stb-mlmd, for OVITO/VMD
# viewer consistency: xsf (multi-frame AXSF, native to OVITO/VMD's xsf
# plugin), pdb (multi-model, VMD's own default), xyz (extended XYZ, OVITO
# -native but ignored as a lattice source by VMD's plain xyz reader).
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}

# Same default vacuum-gap threshold as stb-unitcell/stb-supercell/stb-fetch/
# stb-kgrid/stb-mlrelax (core/kspace.py's other callers) -- used only to
# VALIDATE that the input has no vacuum-padded axis (melting/NPT-relaxing a
# slab/wire/molecule is physically meaningless): this tool is bulk-only, so
# unlike stb-unitcell/stb-supercell there is no vacuum-aware cell mask to
# build -- the final relax always uses cell_mask=[True]*6 (see main()).
VACUUM_GAP_ANG = 10.0


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as stb-unitcell/
    stb-supercell/stb-passivate/stb-crystalbuilder/stb-molecule's own
    _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


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


def _describe_structure(pmg_structure, vacuum_axes, f_out):
    """Formula/atoms/dimensionality/cell parameters -- same field set as
    stb-supercell/stb-crystalbuilder's own [1] section."""
    print_dual(f"Formula        : {pmg_structure.composition.reduced_formula}", f_out)
    print_dual(f"Atoms          : {len(pmg_structure)}", f_out)
    print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    a, b, c, alpha, beta, gamma = pmg_structure.lattice.parameters
    print_dual(f"Cell a,b,c     : {a:.4f}, {b:.4f}, {c:.4f} Ang", f_out)
    print_dual(f"Cell angles    : {alpha:.2f}, {beta:.2f}, {gamma:.2f} deg", f_out)
    print_dual(f"Cell volume    : {pmg_structure.lattice.volume:.4f} Ang^3", f_out)
    print_dual(f"Density        : {len(pmg_structure) / pmg_structure.lattice.volume:.4f} atoms/Ang^3", f_out)


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a
    space-group label -- same shape as stb-supercell/stb-crystalbuilder's
    own _validate_structure(), wrapped in try/except by the caller (a
    validation failure is reported, never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group    : {sg_label}", f_out)
    return sg_label


def write_md_diagnostics(stem, steps, phases, times_fs, temps, epots, ekins, etots,
                          volumes, melt_temp, quench_temp, transition_step):
    """Writes <stem>_md_diagnostics.dat + the companion .gplot, the same
    .dat+.gplot convention used throughout the suite (e.g.
    strain_analysis.py::write_strain_gplot) rather than mlmd.py's matplotlib
    PNG -- gnuplot was explicitly requested for this tool. Two gnuplot
    'index' blocks separate the melt phase (index 0) from the quench phase
    (index 1) in the same file, so either phase can be plotted on its own
    or both together (index 0:1). Returns (dat_path, gplot_path).
    """
    dat_path = f"{stem}_md_diagnostics.dat"
    gplot_path = f"{stem}_md_diagnostics.gplot"

    with open(dat_path, "w") as f:
        f.write("# stb-amorphize melt-quench MD diagnostics\n")
        f.write("# step  time_fs  temp_K  epot_eV  ekin_eV  etot_eV  volume_Ang3\n")
        for phase_name in ("melt", "quench"):
            if phase_name != "melt":
                f.write("\n\n")
            f.write(f"# index {'0' if phase_name == 'melt' else '1'}: {phase_name}\n")
            for i in range(len(steps)):
                if phases[i] != phase_name:
                    continue
                f.write(f"{steps[i]:8d}  {times_fs[i]:12.4f}  {temps[i]:10.3f}  "
                        f"{epots[i]:14.6f}  {ekins[i]:14.6f}  {etots[i]:14.6f}  "
                        f"{volumes[i]:12.4f}\n")

    lines = [
        '# --- STB Plot Configuration ---\n',
        '# Generated by stb-amorphize\n',
        'set terminal pdfcairo enhanced color font "Arial,14" size 8,8\n',
        f'set output "{gplot_path.rsplit(".", 1)[0]}.pdf"\n\n',
        'set multiplot layout 2,1\n\n',
        'set grid\n',
        'set key top right\n',
        f'set arrow from {transition_step}, graph 0 to {transition_step}, graph 1 nohead '
        'lc rgb "#888888" dt 2 front\n\n',
        'set title "Melt-Quench MD: Total Energy"\n',
        'set xlabel "Step"\n',
        'set ylabel "E_{total} (eV)"\n',
        f'plot "{dat_path}" index 0:1 using 1:6 with lines lc rgb "#2255cc" lw 1.5 '
        'title "E_{total}"\n\n',
        'set title "Melt-Quench MD: Temperature"\n',
        'set xlabel "Step"\n',
        'set ylabel "Temperature (K)"\n',
        f'plot "{dat_path}" index 0:1 using 1:3 with lines lc rgb "#cc2222" lw 1.5 '
        'title "T (instantaneous)", \\\n'
        f'     {melt_temp} with lines lc rgb "#888888" dt 2 title "melt target ({melt_temp:.0f} K)", \\\n'
        f'     {quench_temp} with lines lc rgb "#22aa55" dt 2 title "quench target ({quench_temp:.0f} K)"\n\n',
        'unset multiplot\n',
    ]
    with open(gplot_path, "w") as f:
        f.writelines(lines)

    return dat_path, gplot_path


def write_md_trajectory(frames, stem, traj_format):
    """Writes a multi-frame trajectory (xsf/pdb/xyz) of the melt-quench MD
    for viewing in OVITO/VMD -- same OUTPUT_FORMATS/ase.io.write convention
    as stb-ani2traj/stb-mlmd. Returns the path written."""
    from ase.io import write as ase_write

    ase_format, default_ext = OUTPUT_FORMATS[traj_format]
    traj_path = f"{stem}_md_traj{default_ext}"
    ase_write(traj_path, frames, format=ase_format)
    return traj_path


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Melt-quench amorphous structure generator using the MACE-MP-0 potential.", 'bold')}
Heats a crystalline structure above its melting point long enough to erase
crystalline memory, then cools it back down -- a fast heuristic starting
guess for an amorphous/glassy structure, meant to give SIESTA a much better
starting point than an ad-hoc random-displacement structure. NOT a
substitute for a slower, production-quality quench (more steps, slower
cooling) or DFT verification.

Bulk (3D periodic) structures only -- melting a vacuum-padded slab/wire/
molecule is physically meaningless. The input is rejected outright if any
vacuum-padded axis is detected.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f bulk_supercell.fdf\n"
               "  %(prog)s -f bulk_supercell.fdf --melt-temp 4000 --quench-steps 2000\n"
               "  %(prog)s -f bulk_supercell.fdf --save-report --view\n"
               "  %(prog)s -f bulk_supercell.fdf --save-data --save-traj --traj-format xsf\n"
               "\n"
               "Notes:\n"
               "  - Bulk (3D periodic) structures only -- see the description above.\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input crystalline structure file (.fdf), "
                             "typically a supercell (e.g. from stb-supercell). Must be "
                             "bulk (3D periodic, no vacuum-padded axis).")
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
    parser.add_argument("--taut", type=float, default=50.0,
                        help="Berendsen temperature-coupling time constant, fs "
                             "(default: 50.0).")
    parser.add_argument("--taup", type=float, default=200.0,
                        help="Berendsen pressure-coupling time constant, fs "
                             "(default: 200.0).")
    parser.add_argument("--compressibility", type=float, default=4.57e-5,
                        help="Compressibility in eV/Ang^3 for the Berendsen barostat "
                             "(default: 4.57e-5, water's value -- a generic placeholder; "
                             "override with your material's real compressibility if known).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size: speed/accuracy tradeoff (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file, instead of a "
                             "MACE-MP-0 foundation size.")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Seed for the initial Maxwell-Boltzmann velocities "
                             "(default: not reproducible).")
    parser.add_argument("--no-final-relax", dest="final_relax", action="store_false",
                        help="Skip the final static (position+cell) relax. On by "
                             "default: the MD quench alone leaves real residual "
                             "thermal energy, not a clean local minimum.")

    parser.add_argument("--stride", type=int, default=10,
                        help="Sample one diagnostics point / trajectory frame every "
                             "Nth MD step, across both the melt and quench stages "
                             "(default: 10). Only used if --save-data and/or "
                             "--save-traj is given.")
    parser.add_argument("--save-data", action="store_true",
                        help="Write <stem>_md_diagnostics.dat + a companion .gplot "
                             "(step, temperature, E_pot/E_kin/E_total, volume) for the "
                             "whole melt-quench process, plottable with gnuplot. Off "
                             "by default.")
    parser.add_argument("--save-traj", action="store_true",
                        help="Write a multi-frame trajectory of the melt-quench MD "
                             "(see --traj-format) for viewing in OVITO/VMD. Off by "
                             "default.")
    parser.add_argument("--traj-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Trajectory format for --save-traj (default: xsf). Same "
                             "3 choices as stb-ani2traj/stb-mlmd.")

    parser.add_argument("-sp", "--symprec", type=float, default=0.01,
                        help="Symmetry tolerance in Angstroms for the before/after "
                             "symmetry analysis (default: 0.01).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the input "
                             "crystalline structure and the final amorphized structure "
                             "(page through frames in ase-gui) after writing the output "
                             "file. Needs a display. Off by default.")

    parser.add_argument("-o", "--output", type=str, default="amorphous.fdf",
                        help="Output .fdf file name (default: amorphous.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-amorphize {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

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

    print_dual(color_text("===== STB-AMORPHIZE REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"
    model_arg = args.custom_model if args.custom_model else args.model

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Melt protocol  : {args.melt_temp:.0f} K for {args.melt_steps} steps", f_out)
    print_dual(f"Quench protocol: ramp to {args.quench_temp:.0f} K over {args.quench_steps} steps", f_out)
    print_dual(f"Timestep       : {args.timestep} fs", f_out)
    print_dual(f"taut/taup      : {args.taut} fs / {args.taup} fs", f_out)
    print_dual(f"Compressibility: {args.compressibility} eV/Ang^3", f_out)
    print_dual(f"Model          : {model_desc} ({args.device})", f_out)
    print_dual(f"Final relax    : {'yes (position+cell)' if args.final_relax else 'no'}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    if args.save_data or args.save_traj:
        print_dual(f"Sampling stride: every {args.stride} step(s) "
                   f"(save-data={args.save_data}, save-traj={args.save_traj}"
                   f"{', traj-format=' + args.traj_format if args.save_traj else ''})", f_out)

    if not os.path.exists(args.filename):
        _fail(f"File '{args.filename}' not found.", f_out)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        _fail(str(e), f_out)

    pmg_structure = structure_io.to_pymatgen(structure)

    print_section("[1] INPUT STRUCTURE", f_out)
    frac_coords_before = [site.frac_coords for site in pmg_structure]
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_coords_before, pmg_structure.lattice.matrix, VACUUM_GAP_ANG)
    _describe_structure(pmg_structure, vacuum_axes_before, f_out)
    if any(vacuum_axes_before):
        _fail(
            "a vacuum-padded axis was detected (e.g. a slab/wire/molecule). "
            "stb-amorphize only supports bulk (3D periodic) structures -- melting/"
            "NPT-relaxing a vacuum-padded structure is physically meaningless.", f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    try:
        _validate_structure(pmg_structure, vacuum_axes_before, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    mean0, std0 = bond_angle_stats(AseAtomsAdaptor.get_atoms(pmg_structure))
    print_dual(f"Bond-angle mean/std : {mean0:.2f} / {std0:.2f} deg", f_out)

    print_section("[3] MELT-QUENCH MD (MACE)", f_out)
    for line in mace_relax.describe_model(model_arg, mace_relax.get_calculator(model_arg, args.device, "float32")):
        print_dual(line, f_out)

    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.nptberendsen import NPTBerendsen

    atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    atoms.calc = mace_relax.get_calculator(model_arg, args.device, "float32")

    MaxwellBoltzmannDistribution(atoms, temperature_K=min(args.melt_temp, 300.0), rng=(
        np.random.default_rng(args.seed) if args.seed is not None else None))

    # taut/taup/timestep are all in fs (--help) and must be converted to ASE's
    # internal time units (* units.fs) before being passed to NPTBerendsen --
    # the same real bug already found and fixed for `timestep` alone in an
    # earlier pass (units.fs =~ 0.0982, so a raw, unconverted number here
    # silently runs at ~10.18x the requested time constant). taut=50.0/
    # taup=200.0 (fs) still had this exact bug until this pass: fixed here,
    # and both are now also exposed as --taut/--taup overrides (matching
    # stb-mlmd's own --taut/--taup/--compressibility, since these values are
    # a generic, non-material-specific placeholder either way -- see
    # --compressibility's own help text).
    dyn = NPTBerendsen(atoms, timestep=args.timestep * units.fs, temperature_K=args.melt_temp,
                       pressure_au=0.0, taut=args.taut * units.fs, taup=args.taup * units.fs,
                       compressibility_au=args.compressibility)

    # Per-step diagnostics/trajectory collection: a plain manual loop (same
    # pattern as stb-mlmd, no dyn.attach() observer), sampled every --stride
    # steps across both stages so the .dat/trajectory carry one continuous
    # step/time axis through the melt->quench transition.
    collect = args.save_data or args.save_traj
    global_step = 0
    transition_step = None
    d_steps, d_phases, d_times_fs, d_temps = [], [], [], []
    d_epots, d_ekins, d_etots, d_volumes = [], [], [], []
    traj_frames = []

    def _record(phase):
        d_steps.append(global_step)
        d_phases.append(phase)
        d_times_fs.append(global_step * args.timestep)
        d_temps.append(atoms.get_temperature())
        epot = atoms.get_potential_energy()
        ekin = atoms.get_kinetic_energy()
        d_epots.append(epot)
        d_ekins.append(ekin)
        d_etots.append(epot + ekin)
        d_volumes.append(atoms.get_volume())
        if args.save_traj:
            traj_frames.append(atoms.copy())

    print_dual(color_text("Melting...", 'yellow'), f_out)
    start = time.monotonic()
    for step in range(1, args.melt_steps + 1):
        dyn.run(1)
        global_step += 1
        print_progress("Melt", step, args.melt_steps, atoms.get_temperature(), start)
        if collect and global_step % args.stride == 0:
            _record("melt")
    finish_progress_line()

    transition_step = global_step

    print_dual(color_text("Quenching...", 'yellow'), f_out)
    start = time.monotonic()
    for step in range(1, args.quench_steps + 1):
        frac = step / args.quench_steps
        target_t = args.melt_temp + (args.quench_temp - args.melt_temp) * frac
        dyn.set_temperature(temperature_K=target_t)
        dyn.run(1)
        global_step += 1
        print_progress("Quench", step, args.quench_steps, atoms.get_temperature(), start)
        if collect and global_step % args.stride == 0:
            _record("quench")
    finish_progress_line()

    print_dual(f"After MD: T = {atoms.get_temperature():.0f} K", f_out)

    if args.save_data:
        stem = os.path.splitext(os.path.basename(args.output))[0]
        dat_path, gplot_path = write_md_diagnostics(
            stem, d_steps, d_phases, d_times_fs, d_temps, d_epots, d_ekins, d_etots,
            d_volumes, args.melt_temp, args.quench_temp, transition_step)
        print_dual(f"[OK] MD diagnostics ({len(d_steps)} samples, every {args.stride} step(s)) "
                   f"written to '{dat_path}' / '{gplot_path}' (run gnuplot on the latter).", f_out)

    if args.save_traj:
        stem = os.path.splitext(os.path.basename(args.output))[0]
        traj_path = write_md_trajectory(traj_frames, stem, args.traj_format)
        print_dual(f"[OK] MD trajectory ({len(traj_frames)} frame(s), {args.traj_format} format) "
                   f"written to '{traj_path}' -- open in OVITO/VMD.", f_out)

    mean1, std1 = bond_angle_stats(atoms)
    print_dual(f"Bond-angle mean/std : {mean1:.2f} / {std1:.2f} deg (input was {mean0:.2f} / {std0:.2f} deg)", f_out)
    if std1 > std0:
        print_dual(color_text(
            "[OK] Bond-angle spread increased relative to the crystalline input -- "
            "consistent with a loss of long-range order.", 'green'), f_out)
    else:
        print_dual(color_text(
            "[WARNING] Bond-angle spread did NOT increase -- the structure may not have "
            "amorphized. Consider a higher --melt-temp or more --melt-steps.", 'yellow'), f_out)

    if args.final_relax:
        print_section("[4] FINAL STATIC RELAX (MACE)", f_out)
        print_dual("Relaxing positions and cell (float64) to a clean local minimum...", f_out)
        relax_calc = mace_relax.get_calculator(model_arg, args.device, "float64")
        e0 = atoms.get_potential_energy()
        f0 = float(np.abs(atoms.get_forces()).max())
        t0 = time.time()
        # cell_mask=[True]*6, not None: the input was already validated bulk
        # (no vacuum axis) above, so every cell direction is free to relax --
        # a real, verified bug until this pass: cell_mask=None disables cell
        # relaxation entirely in core.mace_relax.relax() (positions only),
        # contradicting this section's own "position+cell" label and the
        # tool's documented intent.
        converged, steps_used = mace_relax.relax(atoms, relax_calc, cell_mask=[True] * 6, fmax=0.05, max_steps=200)
        wall_time = time.time() - t0
        e1 = atoms.get_potential_energy()
        f1 = float(np.abs(atoms.get_forces()).max())
        print_dual(f"Steps used : {steps_used} "
                   f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
        print_dual(f"Wall time  : {wall_time:.1f} s", f_out)
        n_atoms = len(atoms)
        rows = [
            (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
              f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
            (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
        ]
        print_table(["Quantity", "Before", "After", "Change"], rows, f_out)
        final_relax_info = (converged, steps_used, e0, e1)
    else:
        final_relax_info = None

    final_pmg = AseAtomsAdaptor.get_structure(atoms)

    print_section("[5] STRUCTURE VALIDATION (post-transform)", f_out)
    frac_coords_final = final_pmg.frac_coords
    vacuum_axes_final = kspace.detect_vacuum_axes(frac_coords_final, final_pmg.lattice.matrix, VACUUM_GAP_ANG)
    try:
        _validate_structure(final_pmg, vacuum_axes_final, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[6] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    print_dual("Before = crystalline input. After = melt-quenched structure -- UNLIKE "
               "stb-unitcell's own before/after table, these are EXPECTED to differ: a "
               "genuine amorphization loses the crystal's long-range symmetry (often down "
               "to P1, no symmetry at all), which is itself a sign the melt-quench worked.", f_out)
    # No "Layer Group" row here (unlike stb-unitcell/stb-slab) -- this tool is
    # bulk-only (vacuum rejected above), so it would always read
    # "N/A (not 2D-periodic)", never actual information -- same reasoning as
    # stb-crystalbuilder's own before/after table.
    before_info = core_symmetry.symmetry_summary(pmg_structure, args.symprec, VACUUM_GAP_ANG)
    after_info = core_symmetry.symmetry_summary(final_pmg, args.symprec, VACUUM_GAP_ANG)
    if "Error" in before_info or "Error" in after_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
        print_dual(f"  Before: {before_info.get('Error', 'OK')}", f_out)
        print_dual(f"  After : {after_info.get('Error', 'OK')}", f_out)
    else:
        properties = ["Crystal System", "Space Group", "Point Group", "Hall Symbol"]
        rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                for prop in properties]
        print_table(["Property", "Before", "After"], rows, f_out)
        if after_info.get("Space Group", "").startswith("P1"):
            print_dual(color_text(
                "[OK] Final space group is P1 (no symmetry) -- consistent with a genuine "
                "amorphous structure.", 'green'), f_out)

    print_section("[7] WRITING OUTPUT FILE", f_out)
    header_comment = [
        "Amorphous structure generated by stb-amorphize (melt-quench, MACE-MP-0).",
        f"Input file: {args.filename}.",
        f"Protocol: melt {args.melt_temp:.0f} K for {args.melt_steps} steps, "
        f"quench to {args.quench_temp:.0f} K over {args.quench_steps} steps "
        f"(timestep {args.timestep} fs, taut={args.taut} fs, taup={args.taup} fs).",
        f"Bond-angle mean/std: {mean0:.2f}/{std0:.2f} deg (input) -> {mean1:.2f}/{std1:.2f} deg (after MD).",
    ]
    if final_relax_info is not None:
        converged, steps_used, e0, e1 = final_relax_info
        header_comment.append(
            f"Final static relax with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )
    header_comment.append(
        "NOTE: a fast melt-quench heuristic, not a substitute for a slower "
        "production-quality quench or DFT verification."
    )
    new_structure = structure_io.from_pymatgen(final_pmg, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[8] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.MACE]
    if not args.custom_model:
        bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[9] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)
    print_dual(color_text(
        "Note: a fast melt-quench heuristic, not a substitute for a slower "
        "production-quality quench or DFT verification.", 'yellow'), f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # both frames (crystalline input vs. final amorphized structure) so the
    # user can page through the actual comparison in ase-gui.
    if args.view:
        input_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        final_atoms = AseAtomsAdaptor.get_atoms(final_pmg)
        view_structure_interactive([input_atoms, final_atoms])


if __name__ == "__main__":
    main()
