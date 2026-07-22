#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Vacancy migration-barrier screening via MACE-driven NEB -- the eighth
tool in the ML Simulations category. Where stb-mlneb needs two already
-built endpoint structures, this tool builds them itself: given a PERFECT
bulk reference structure and the index of one atom to remove (--vacancy
-index), it auto-detects every neighbor within the first coordination
shell (--shell-tolerance around the nearest-neighbor distance), builds the
(vacancy-at-site, vacancy-at-neighbor-site) endpoint pair for each
candidate hop, and runs a real climbing-image NEB (via stb-mlneb's own
run_single_neb -- import-don't-duplicate, not a second copy of the NEB
pipeline) on every one of them, reporting a migration barrier per
candidate. Distinct from stb-mlmd/stb-aimdAnalysis's indirect, MSD-based
diffusion coefficient: this gives the actual migration barrier directly,
useful in the rare-event/low-temperature regime where an MD trajectory
would need an infeasibly long run to see even one real hop.

Scope: single-vacancy migration in a bulk (3D periodic) monoatomic-shell
sense only -- the neighbor search doesn't distinguish chemical species, so
for a multi-species material every neighbor within the shell (regardless
of species) is treated as a candidate hop target; a hop onto a site of a
different species is still built and run (composition is preserved, since
only the SAME vacancy moves), but is not a meaningful migration path for
most real materials and should be interpreted with that in mind.
"""

VERSION = "1.0.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.neb import wrap_into_cell, write_ml_preview_plot, write_path_trajectory
from stb.mlneb import run_single_neb

require_mace()

REPORT_FILE = "stb_mldiffusion_report.txt"


def find_neighbor_shell(atoms, vacancy_index, shell_tolerance):
    """Returns the sorted list of atom indices (excluding vacancy_index
    itself) within `shell_tolerance` fraction of the nearest-neighbor
    distance FROM vacancy_index specifically (not the structure's global
    minimum distance -- matters for a multi-site/defective structure where
    the vacancy's own local environment may differ slightly from the bulk
    average), plus that reference distance itself. E.g. shell_tolerance=
    0.15 keeps every neighbor within 15% of the closest one -- catches the
    full first coordination shell (all physically equivalent nearest
    neighbors) without also pulling in the next (2nd-shell) atoms, for a
    typical crystal's shell-distance separation.
    """
    distances = atoms.get_distances(vacancy_index, range(len(atoms)), mic=True)
    distances[vacancy_index] = np.inf
    d_min = distances.min()
    cutoff = d_min * (1.0 + shell_tolerance)
    neighbors = [i for i in range(len(atoms)) if distances[i] <= cutoff]
    return neighbors, d_min


def build_vacancy_hop_endpoints(atoms, vacancy_index, neighbor_index):
    """Builds the (initial, final) ASE Atoms pair for one candidate vacancy
    hop: `initial` has the vacancy at `vacancy_index` (that atom removed,
    every other atom -- including the one at `neighbor_index` -- exactly
    where it started); `final` keeps the exact same atom-index
    correspondence (critical for NEB), except the atom that was at
    `neighbor_index` is moved to `vacancy_index`'s old position, leaving
    the vacancy at `neighbor_index` instead -- the same construction
    verified live for stb-mlneb's own vacancy-migration test fixture.
    """
    frac = atoms.get_scaled_positions()
    numbers = atoms.numbers
    cell = atoms.get_cell()

    remaining = [i for i in range(len(atoms)) if i != vacancy_index]
    initial_frac = frac[remaining].copy()
    final_frac = initial_frac.copy()
    list_pos = remaining.index(neighbor_index)
    final_frac[list_pos] = frac[vacancy_index]

    initial_atoms = Atoms(numbers=numbers[remaining], scaled_positions=initial_frac, cell=cell, pbc=True)
    final_atoms = Atoms(numbers=numbers[remaining], scaled_positions=final_frac, cell=cell, pbc=True)
    return initial_atoms, final_atoms


def plot_barriers(labels, barriers, distances, out_path):
    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(labels)), 5))
    bars = ax.bar(labels, barriers, color='tab:blue')
    best_idx = int(np.argmin(barriers))
    bars[best_idx].set_color('tab:orange')
    ax.set_xlabel("Candidate hop (neighbor atom index)")
    ax.set_ylabel("Migration barrier (eV)")
    ax.set_title("Vacancy migration barriers (orange = lowest)")
    ax.tick_params(axis='x', rotation=30)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Vacancy migration-barrier screening via MACE-driven NEB -- given a "
                     "PERFECT bulk reference structure and the index of one atom to remove "
                     "(--vacancy-index), auto-detects every neighbor in the first "
                     "coordination shell, builds the vacancy-hop endpoint pair for each, and "
                     "runs a real climbing-image NEB (reusing stb-mlneb's own pipeline) on "
                     "every candidate, reporting a migration barrier per hop. No SIESTA and "
                     "no separate endpoint structures to build by hand.",
        epilog="Example usage:\n"
               "  stb-mldiffusion --file bulk.fdf --vacancy-index 0\n"
               "  stb-mldiffusion --file bulk.fdf --vacancy-index 0 --custom-model "
               "mlff_model.model\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Perfect bulk reference structure (.fdf).")
    parser.add_argument("--vacancy-index", type=int, required=True,
                        help="0-based index of the atom to remove, creating the vacancy.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--shell-tolerance", type=float, default=0.15,
                        help="Fraction above the nearest-neighbor distance still counted as "
                             "the first coordination shell (default: 0.15).")
    parser.add_argument("--n-images", type=int, default=5,
                        help="Total images including both endpoints, per candidate hop (default: 5).")
    parser.add_argument("--autosort-tol", type=float, default=0.5,
                        help="Atom-correspondence tolerance, Ang (default: 0.5).")
    parser.add_argument("--no-idpp", dest="idpp", action="store_false",
                        help="Skip IDPP refinement of interior images. On by default.")
    parser.add_argument("--idpp-fmax", type=float, default=0.1, help="IDPP force convergence (default: 0.1).")
    parser.add_argument("--idpp-steps", type=int, default=100, help="Max IDPP steps (default: 100).")
    parser.add_argument("--no-pre-relax-endpoints", dest="pre_relax_endpoints", action="store_false",
                        help="Skip the MACE relax (positions only) of each hop's endpoints. On by default.")
    parser.add_argument("--freeze-substrate", action="store_true",
                        help="Freeze atoms that barely move between each hop's endpoints "
                             "during the NEB relax. Off by default.")
    parser.add_argument("--freeze-threshold", type=float, default=0.3,
                        help="Displacement threshold for --freeze-substrate, Ang (default: 0.3).")
    parser.add_argument("--k", type=float, default=0.1, help="NEB spring constant, eV/Ang^2 (default: 0.1).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for relax/NEB, eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200, help="Max steps per stage (default: 200).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mldiffusion_out",
                        help="Output directory for all files (default: mldiffusion_out).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mldiffusion {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Vacancy Diffusion",
            "Vacancy migration-barrier screening via MACE-driven NEB",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML VACANCY MIGRATION-BARRIER SCREENING:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLDIFFUSION REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Vacancy index     : {args.vacancy_index}", f_out)
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
    print_dual(f"[OK] Read {len(atoms)} atom(s)", f_out)

    if not (0 <= args.vacancy_index < len(atoms)):
        fail(f"--vacancy-index {args.vacancy_index} out of range (structure has {len(atoms)} atoms, 0-based).")

    print_section("[2] NEIGHBOR SHELL DETECTION", f_out)
    neighbors, d_min = find_neighbor_shell(atoms, args.vacancy_index, args.shell_tolerance)
    print_dual(f"Nearest-neighbor distance from site {args.vacancy_index}: {d_min:.4f} Ang", f_out)
    print_dual(f"Candidate hop target(s) (within {args.shell_tolerance*100:.0f}% of that): "
               f"{neighbors}", f_out)
    if not neighbors:
        fail("No neighbor found within the shell tolerance -- check --vacancy-index or "
             "widen --shell-tolerance.")

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)

    print_section("[3] MIGRATION-BARRIER SCAN", f_out)
    results = []
    for neighbor_index in neighbors:
        tag = f"hop->{neighbor_index}"
        print_dual(f"\n--- Candidate hop: site {args.vacancy_index} -> site {neighbor_index} ---", f_out)
        initial_atoms, final_atoms = build_vacancy_hop_endpoints(atoms, args.vacancy_index, neighbor_index)
        initial_pmg = wrap_into_cell(AseAtomsAdaptor.get_structure(initial_atoms))
        final_pmg = wrap_into_cell(AseAtomsAdaptor.get_structure(final_atoms))

        barrier, dE, ase_images, pmg_images, converged = run_single_neb(
            initial_pmg, final_pmg, calc, args, f_out, tag=tag)
        results.append({
            "neighbor_index": neighbor_index, "barrier": barrier, "dE": dE,
            "converged": converged, "ase_images": ase_images, "pmg_images": pmg_images,
        })

    print_section("[4] SUMMARY & FILES", f_out)
    best = min(results, key=lambda r: r["barrier"])
    print_dual(f"{'Hop target':>12}{'Barrier (eV)':>16}{'Reaction E (eV)':>18}{'Converged':>12}", f_out)
    for r in results:
        marker = color_text(" <- lowest", 'green') if r is best else ""
        print_dual(f"{r['neighbor_index']:>12}{r['barrier']:16.4f}{r['dE']:18.4f}"
                   f"{str(r['converged']):>12}{marker}", f_out)

    labels = [str(r["neighbor_index"]) for r in results]
    barriers = [r["barrier"] for r in results]
    distances = [d_min] * len(results)
    plot_path = os.path.join(args.output_dir, "migration_barriers.png")
    plot_barriers(labels, barriers, distances, plot_path)
    print_dual(f"[OK] Saved {plot_path}", f_out)

    best_preview = os.path.join(args.output_dir, "lowest_barrier_preview.png")
    write_ml_preview_plot(best["ase_images"], best_preview)
    print_dual(f"[OK] Saved {best_preview} (lowest-barrier hop: -> site {best['neighbor_index']})", f_out)

    best_traj = os.path.join(args.output_dir, "lowest_barrier_path.xyz")
    write_path_trajectory(best["pmg_images"], best_traj)
    print_dual(f"[OK] Saved {best_traj}", f_out)

    print_dual("Status            : OK", f_out)
    print_dual(f"Lowest barrier    : {best['barrier']:.4f} eV (site {args.vacancy_index} -> "
               f"site {best['neighbor_index']})", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML vacancy migration-barrier screening complete -- a fast heuristic, "
                      "not a substitute for a DFT-level barrier.\n", 'bold'))


if __name__ == "__main__":
    main()
