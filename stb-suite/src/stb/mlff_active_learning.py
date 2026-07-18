#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import glob
import argparse
from datetime import datetime

import numpy as np
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io
from stb.core.pseudopotentials import get_required_pseudos
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.mlff import build_rattled_configs, write_config_folder

REPORT_FILE = "stb_mlffActiveLearning_report.txt"

_CONFIG_DIR_RE = re.compile(r"^mlff_config_(\d+)$")


def next_config_index(output_dir):
    """Returns 1 + the highest NNN already used by an existing
    mlff_config_NNN/ folder in `output_dir` (0 if none) -- so a round of
    active-learning-suggested configs continues the numbering from a prior
    stb-mlff/stb-mlffActiveLearning run instead of overwriting it.
    """
    highest = 0
    for d in glob.glob(os.path.join(output_dir, "mlff_config_*")):
        m = _CONFIG_DIR_RE.match(os.path.basename(d.rstrip(os.sep)))
        if m:
            highest = max(highest, int(m.group(1)))
    return highest + 1


def score_candidates(atoms_ref, coord_format, positions_list, cells, calc):
    """Evaluates `calc` (the fine-tuned MACE model) on each candidate
    (position array + optional cell) and returns a list of (score,
    max|F|, positions) sorted by score descending.

    Score is simply the maximum predicted atomic force magnitude after a
    single-point evaluation at the (unrelaxed, rattled) geometry -- NOT a
    rigorous uncertainty estimate (that would need a committee of several
    independently-trained models and disagreement between them). This is a
    much cheaper proxy: a configuration where the potential itself predicts
    an unusually large local force is a reasonable, if imperfect, signal
    that this region of configuration space is poorly represented in the
    current training set -- large real DFT forces at a genuinely
    rattled/strained geometry are expected, but an especially large
    predicted force (relative to the rest of the candidate batch) is still
    informative about where the model is least confident, at a fraction of
    the cost of training multiple models.
    """
    scored = []
    for positions, cell in zip(positions_list, cells):
        if coord_format == "fractional":
            cart = positions @ cell
        else:
            cart = positions
        atoms = Atoms(symbols=atoms_ref.get_chemical_symbols(), positions=cart, cell=cell, pbc=True)
        atoms.calc = calc
        try:
            forces = atoms.get_forces()
        except Exception:
            continue
        max_f = float(np.abs(forces).max())
        scored.append((max_f, positions))
    scored.sort(key=lambda x: x[0], reverse=True)
    return scored


def main():
    parser = argparse.ArgumentParser(
        description="Stage 3 (Active Learning) of the Custom ML Force Field workflow: uses "
                     "an already fine-tuned model (stb-mlffAnalysis) to screen a fresh batch "
                     "of candidate configurations and writes the ones with the largest "
                     "predicted atomic force as new mlff_config_NNN/ folders, ready for the "
                     "user to label with real SIESTA and feed back into another round of "
                     "stb-mlffAnalysis. NOT rigorous uncertainty quantification (that needs a "
                     "committee of independently-trained models) -- a cheap, single-model "
                     "heuristic for where the current model is likely weakest.",
        epilog="Example usage:\n"
               "  stb-mlffActiveLearning --model mlff_model_compiled.model --file calc.fdf \\\n"
               "                          --n-candidates 50 --stdev 0.15 --top-k 5",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--model", required=True, help="Fine-tuned MACE model file (from stb-mlffAnalysis).")
    parser.add_argument("--file", required=True,
                        help="Reference SIESTA .fdf (same one stb-mlff used, or any structure "
                             "of the same system) -- basis/pseudos/SCF settings are reused "
                             "verbatim in every candidate written out.")
    parser.add_argument("--pseudo-dir", default=".",
                        help="Directory to search for <element>.psf/.psml pseudopotentials "
                             "(default: current directory).")
    parser.add_argument("--n-candidates", type=int, default=50,
                        help="Number of candidate configurations to screen (default: 50).")
    parser.add_argument("--stdev", type=float, nargs="+", default=[0.15],
                        help="Rattling standard deviation(s) in Angstrom for the candidate "
                             "batch (default: 0.15 -- typically larger than the amplitude(s) "
                             "used for the original training set, to probe configuration "
                             "space not yet covered).")
    parser.add_argument("--top-k", type=int, default=5,
                        help="Number of highest-scoring candidates to write out (default: 5).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for rattling (default: 42).")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"],
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default=".",
                        help="Directory under which new mlff_config_NNN/ folders are created "
                             "(default: current directory -- continues the numbering of any "
                             "already present there).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlffActiveLearning {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Active Learning)",
            "Suggest new SIESTA configurations to label, using the fine-tuned model itself",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CUSTOM ML FORCE FIELD -- ACTIVE LEARNING:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLFFACTIVELEARNING REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Model             : {args.model}", f_out)
    print_dual(f"Reference file    : {args.file}", f_out)
    print_dual(f"Candidates/top-k  : {args.n_candidates} / {args.top_k} (stdev {args.stdev})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING REFERENCE STRUCTURE & MODEL", f_out)

    if not os.path.isfile(args.file):
        fail(f"Reference file '{args.file}' not found.")
    if not os.path.isfile(args.model):
        fail(f"Model file '{args.model}' not found.")

    try:
        ref_structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")

    unique_elements = sorted(ref_structure.species)
    print_dual(f"[OK] Read {len(ref_structure.atoms)} atom(s), species: {', '.join(unique_elements)}", f_out)

    pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)
    if missing:
        fail(f"Missing pseudopotential(s) for: {', '.join(missing)} (searched '{args.pseudo_dir}').")

    pmg_structure = structure_io.to_pymatgen(ref_structure)
    ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    ref_basename = os.path.basename(args.file)

    try:
        from mace.calculators import MACECalculator
    except ImportError:
        fail("mace-torch not found. This tool needs the optional 'ml' extra: pip install stb_suite[ml]")
    calc = MACECalculator(model_paths=[args.model], device=args.device)
    print_dual(f"[OK] Loaded model {args.model}", f_out)

    print_section("[2] SCREENING CANDIDATES", f_out)

    candidate_positions = build_rattled_configs(
        ase_atoms, ref_structure.coord_format, args.n_candidates, args.stdev, args.seed)
    cells = [np.array(ase_atoms.cell)] * len(candidate_positions)

    scored = score_candidates(ase_atoms, ref_structure.coord_format, candidate_positions, cells, calc)
    if not scored:
        fail("None of the candidate configurations could be evaluated by the model.")

    print_dual(f"[OK] Screened {len(scored)}/{args.n_candidates} candidate(s).", f_out)
    print_dual(f"Predicted max|F| range: {scored[-1][0]:.3f} - {scored[0][0]:.3f} eV/Ang", f_out)

    top_k = scored[:args.top_k]

    print_section("[3] WRITING SUGGESTED CONFIGURATIONS", f_out)
    os.makedirs(args.output_dir, exist_ok=True)
    start_index = next_config_index(args.output_dir)
    folders = []
    for offset, (score, positions) in enumerate(top_k):
        index = start_index + offset
        folder = write_config_folder(
            args.output_dir, index, args.file, ref_basename, positions, None, pseudos_to_copy)
        folders.append((folder, score))
        print_dual(f"  - {folder}  (predicted max|F| = {score:.3f} eV/Ang)", f_out)

    print_section("[4] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Suggested configs : {len(folders)}", f_out)
    print_dual("Next steps        : run SIESTA in each folder above, then re-run "
               "stb-mlffAnalysis with a --path matching all mlff_config_* folders "
               "(old + new) to refine the model.", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Candidates suggested -- label them, then retrain.\n", 'bold'))


if __name__ == "__main__":
    main()
