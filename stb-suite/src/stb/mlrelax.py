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
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
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

require_mace()

OPTIMIZERS = ["FIRE", "BFGS", "LBFGS"]

REPORT_FILE = "stb_mlrelax_report.txt"
BIB_FILE = "references.bib"


def _structure_validation(pmg_structure, vacuum_axes, f_out):
    """The shared malformation checklist (core.structure_checks, also used by
    stb-inputfile/stb-fetch) plus a space-group label. Shared shape between
    the pre- and post-relaxation validation sections below -- wrapped in
    try/except by the caller, same convention as stb-inputfile/stb-fetch (a
    validation failure is reported, never fatal).
    """
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)

    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group    : {sg_label}", f_out)
    return sg_label


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
               "  %(prog)s -f defect.fdf --save-report --save-data --view\n"
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
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw per-step relaxation convergence trace "
                             "(<stem>_relax_convergence.dat). Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the structure "
                             "before and after relaxation (page through frames in ase-gui). "
                             "Needs a display. Off by default.")
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

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLRELAX REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")
    print_dual(f"Model          : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"Device         : {args.device}", f_out)
    print_dual("Precision      : float64 (MACE's own guidance for geometry optimization, "
               "not MD)", f_out)
    print_dual(f"Optimizer      : {args.optimizer}", f_out)
    print_dual(f"fmax           : {args.fmax} eV/Ang", f_out)
    print_dual(f"Max steps      : {args.max_steps}", f_out)
    print_dual(f"Relax cell     : {'yes' if args.relax_cell else 'no (positions only)'}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)

    if not os.path.exists(args.filename):
        fail(f"File '{args.filename}' not found.")

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        fail(str(e))

    pmg = structure_io.to_pymatgen(structure)

    print_section("[1] INPUT STRUCTURE", f_out)
    frac_coords = [site.frac_coords for site in pmg]
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, pmg.lattice.matrix, args.vacuum_gap)
    print_dual(f"Formula        : {pmg.composition.reduced_formula}", f_out)
    print_dual(f"Atoms          : {len(pmg)}", f_out)
    print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    a, b, c, alpha, beta, gamma = pmg.lattice.parameters
    print_dual(f"Cell a,b,c     : {a:.4f}, {b:.4f}, {c:.4f} Ang", f_out)
    print_dual(f"Cell angles    : {alpha:.2f}, {beta:.2f}, {gamma:.2f} deg", f_out)
    print_dual(f"Cell volume    : {pmg.lattice.volume:.4f} Ang^3", f_out)
    if sum(vacuum_axes) == 0:
        print_dual(f"Density        : {len(pmg) / pmg.lattice.volume:.4f} atoms/Ang^3", f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-relaxation)", f_out)
    try:
        sg_before = _structure_validation(pmg, vacuum_axes, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)
        sg_before = None

    print_section("[3] MACE MODEL & SIMULATION SETUP", f_out)
    cell_mask = None
    if args.relax_cell:
        axis_names = ("a", "b", "c")
        if all(vacuum_axes):
            print_dual(color_text(
                "Note: every axis is vacuum-padded (a fully isolated structure, e.g. a "
                "molecule) -- there is no periodic direction to relax. Ignoring --relax-cell, "
                "doing positions-only.", 'yellow'), f_out)
        else:
            cell_mask = mace_relax.build_cell_mask(vacuum_axes)
            fixed = [axis_names[i] for i, v in enumerate(vacuum_axes) if v]
            relaxed = [axis_names[i] for i, v in enumerate(vacuum_axes) if not v]
            if fixed:
                print_dual(f"Vacuum axis/axes held fixed: {', '.join(fixed)} "
                           f"(gap >= {args.vacuum_gap} Ang detected)", f_out)
            print_dual(f"Relaxing lattice along: {', '.join(relaxed)}", f_out)

    atoms = AseAtomsAdaptor.get_atoms(pmg)
    atoms_before = atoms.copy()

    print_dual(f"Mode           : {'positions + cell' if cell_mask is not None else 'positions only'}", f_out)
    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    for line in mace_relax.describe_model(model_arg, calc):
        print_dual(line, f_out)
    atoms.calc = calc

    e0 = atoms.get_potential_energy()
    f0 = np.abs(atoms.get_forces()).max()

    print_section("[4] RELAXATION", f_out)
    print_dual(f"Before         : E = {e0:.4f} eV, max|F| = {f0:.4f} eV/Ang", f_out)
    print_dual(f"Relaxing... (fmax={args.fmax}, max {args.max_steps} steps)", f_out)

    history = []
    t0 = time.time()
    converged, steps_used = mace_relax.relax(
        atoms, calc, cell_mask=cell_mask, optimizer=args.optimizer,
        fmax=args.fmax, max_steps=args.max_steps, step_history=history)
    wall_time = time.time() - t0

    e1 = atoms.get_potential_energy()
    f1 = np.abs(atoms.get_forces()).max()
    print_dual(f"After          : E = {e1:.4f} eV, max|F| = {f1:.4f} eV/Ang", f_out)
    print_dual(f"Steps used     : {steps_used} "
               f"({color_text('converged', 'green') if converged else color_text('hit step cap, NOT converged', 'yellow')})", f_out)
    print_dual(f"Wall time      : {wall_time:.1f} s", f_out)

    stem = os.path.splitext(os.path.basename(args.output))[0]
    if len(history) >= 2:
        steps_arr = np.array([h[0] for h in history])
        energies_arr = np.array([h[1] for h in history])
        fmax_arr = np.array([h[2] for h in history])

        conv_plot_path = f"{stem}_relax_convergence.png"
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].plot(steps_arr, energies_arr, '-', color='tab:blue')
        axes[0].set_xlabel("Optimizer step")
        axes[0].set_ylabel("Energy (eV)")
        axes[0].set_title("Energy convergence")
        axes[1].plot(steps_arr, fmax_arr, '-', color='tab:red')
        axes[1].axhline(args.fmax, color='gray', linestyle='--', label=f"fmax target ({args.fmax})")
        axes[1].set_yscale('log')
        axes[1].set_xlabel("Optimizer step")
        axes[1].set_ylabel("max|F| (eV/Ang)")
        axes[1].set_title("Force convergence")
        axes[1].legend()
        fig.tight_layout()
        fig.savefig(conv_plot_path, dpi=150)
        plt.close(fig)
        print_dual(f"[OK] Saved convergence plot: {conv_plot_path}", f_out)

        if args.save_data:
            conv_data_path = f"{stem}_relax_convergence.dat"
            with open(conv_data_path, "w") as f_dat:
                f_dat.write("# step energy(eV) max_force(eV/Ang)\n")
                for step, energy, fmax_step in history:
                    f_dat.write(f"{step} {energy:.8f} {fmax_step:.6f}\n")
            print_dual(f"[OK] Saved convergence data: {conv_data_path}", f_out)
    else:
        conv_plot_path = None
        print_dual("Note: fewer than 2 optimizer steps taken -- no convergence plot saved "
                   "(structure was already essentially at this potential's minimum).", f_out)

    print_section("[5] BEFORE / AFTER COMPARISON (detailed)", f_out)
    n_atoms = len(atoms)
    rows = [
        (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
          f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
        (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
    ]

    lattice_warn = False
    if cell_mask is not None:
        cellpar_before = atoms_before.cell.cellpar()
        cellpar_after = atoms.cell.cellpar()
        rel_change = np.abs(cellpar_after[:3] - cellpar_before[:3]) / cellpar_before[:3]
        lattice_warn = rel_change.max() > 0.15
        rows.append((["Lattice a,b,c (Ang)",
                      f"{cellpar_before[0]:.4f}, {cellpar_before[1]:.4f}, {cellpar_before[2]:.4f}",
                      f"{cellpar_after[0]:.4f}, {cellpar_after[1]:.4f}, {cellpar_after[2]:.4f}",
                      f"max {100 * rel_change.max():+.2f}%"], None))
        rows.append((["Lattice angles (deg)",
                      f"{cellpar_before[3]:.2f}, {cellpar_before[4]:.2f}, {cellpar_before[5]:.2f}",
                      f"{cellpar_after[3]:.2f}, {cellpar_after[4]:.2f}, {cellpar_after[5]:.2f}",
                      ""], None))
        vol_before, vol_after = atoms_before.get_volume(), atoms.get_volume()
        rows.append((["Volume (Ang^3)", f"{vol_before:.4f}", f"{vol_after:.4f}",
                      f"{100 * (vol_after - vol_before) / vol_before:+.2f}%"], None))
        if sum(vacuum_axes) == 0:
            rows.append((["Density (atoms/Ang^3)", f"{n_atoms / vol_before:.4f}",
                          f"{n_atoms / vol_after:.4f}",
                          f"{n_atoms / vol_after - n_atoms / vol_before:+.4f}"], None))

    disp = np.linalg.norm(atoms.get_positions() - atoms_before.get_positions(), axis=1)
    max_disp_idx = int(np.argmax(disp))
    rows.append((["Displacement (Ang)", "--", "--",
                  f"mean={disp.mean():.4f} RMS={np.sqrt((disp ** 2).mean()):.4f} "
                  f"max={disp.max():.4f} (atom #{max_disp_idx} = {atoms[max_disp_idx].symbol})"], None))

    pmg_before = AseAtomsAdaptor.get_structure(atoms_before)
    pmg_after = AseAtomsAdaptor.get_structure(atoms)

    try:
        dist_before = structure_io.min_pairwise_distance(pmg_before)
        dist_after = structure_io.min_pairwise_distance(pmg_after)
        if dist_before is not None and dist_after is not None:
            rows.append((["Min. distance (Ang)", f"{dist_before:.4f}", f"{dist_after:.4f}",
                          f"{dist_after - dist_before:+.4f}"], None))
    except Exception:
        pass

    sg_after = None
    try:
        sg_after = core_symmetry.space_group_label(pmg_after)
        rows.append((["Space group", sg_before, sg_after,
                      "unchanged" if sg_before == sg_after else "CHANGED"], None))
    except Exception:
        pass

    # Best-effort only: a custom fine-tuned model may lack a stress head
    # (mace_run_train can be run without stress targets), so a failure here
    # is reported and skipped, never fatal.
    max_stress_after = None
    try:
        from ase import units
        atoms_before.calc = calc
        max_stress_before = np.abs(atoms_before.get_stress(voigt=False)).max() / units.GPa
        max_stress_after = np.abs(atoms.get_stress(voigt=False)).max() / units.GPa
        rows.append((["Max |stress| (GPa)", f"{max_stress_before:.4f}", f"{max_stress_after:.4f}",
                      f"{max_stress_after - max_stress_before:+.4f}"], None))
    except Exception as e:
        rows.append((["Max |stress| (GPa)", "n/a", "n/a", f"not available ({e})"], None))

    print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

    if lattice_warn:
        print_dual(color_text(
            "Warning: a relaxed lattice parameter changed by more than 15% -- "
            "universal MLIPs can occasionally settle into a spurious state this far "
            "from the initial guess. Sanity-check this result before trusting it.",
            'yellow'), f_out)
    if cell_mask is not None:
        print_dual(
            "Note: Displacement above also reflects lattice strain (--relax-cell was "
            "used), not only ionic motion.", f_out)

    print_section("[6] STRUCTURE VALIDATION (post-relaxation)", f_out)
    try:
        _structure_validation(pmg_after, vacuum_axes, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[7] WRITING OUTPUT FILE", f_out)
    new_pmg = AseAtomsAdaptor.get_structure(atoms)
    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=structure.species_meta)
    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"
    stress_desc = f"{max_stress_after:.4f} GPa" if max_stress_after is not None else "not available"
    header_comment = [
        f"Structure relaxed by stb-mlrelax using {model_desc}.",
        f"Optimizer: {args.optimizer}, {steps_used} step(s), "
        f"{'converged' if converged else 'NOT converged (hit step cap)'}.",
        f"Convergence: max stress = {stress_desc}.",
        f"Energy: {e1:.6f} eV (delta {e1 - e0:+.6f} eV from the input structure).",
    ]
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
    if conv_plot_path:
        print_dual(f"Convergence plot: {conv_plot_path}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # both frames (before/after) so the user can page through the actual
    # comparison in ase-gui, not just the final structure.
    if args.view:
        view_structure_interactive([atoms_before, atoms])


if __name__ == "__main__":
    main()
