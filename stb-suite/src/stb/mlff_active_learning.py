#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import glob
import math
import shutil
import argparse
import tempfile
from datetime import datetime

import numpy as np
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, mlff_workflow
from stb.core.fdf_directives import parse_kgrid, parse_directives
from stb.core.pseudopotentials import get_required_pseudos
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.mlff import build_rattled_configs, write_config_folder
from stb.mlmd import build_dynamics

REPORT_FILE = "stb_mlffActiveLearning_report.txt"

_CONFIG_DIR_RE = re.compile(r"^mlff_config_(\d+)$")

# Atomic-environment fingerprint used for the novelty/diversity filter: the
# sorted distances from an atom to its FP_NEIGHBOURS nearest neighbours
# (periodic images included) within FP_CUTOFF Ang, padded with the cutoff.
# Model-independent (MACECalculator.get_descriptors does not work on the
# compiled model stb-mlffAnalysis also writes) and comparable across cell
# sizes, which a supercell candidate against an 8-atom training set needs.
FP_NEIGHBOURS = 16
FP_CUTOFF = 6.0
DEFAULT_MIN_NOVELTY = 0.5


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


# --------------------------------------------------------------------------
# Scoring: committee disagreement (2+ models) or max predicted force (1 model)
# --------------------------------------------------------------------------

def evaluate_candidate(atoms, calcs):
    """Single-point forces of every model in `calcs` on `atoms`. Returns
    {'max_f': max per-atom |F| of the committee-mean force (eV/Ang),
     'sigma': max over atoms of the committee force spread (eV/Ang) -- the
              RMS over models of |F_model - F_mean| per atom -- or None for a
              single model,
     'score': sigma with a committee, max_f with one model}
    or None if any model fails (e.g. a NaN force).

    The committee spread is the actual uncertainty signal: models trained
    on the same data with different seeds agree where the data constrain
    them and drift apart where they do not. One model's own max force only
    says the geometry is far from equilibrium, which selects the most
    distorted candidates rather than the ones the model gets wrong.
    """
    forces = []
    for calc in calcs:
        a = atoms.copy()
        a.calc = calc
        try:
            f = a.get_forces()
        except Exception:
            return None
        if not np.all(np.isfinite(f)):
            return None
        forces.append(f)
    forces = np.array(forces)
    mean = forces.mean(axis=0)
    max_f = float(np.linalg.norm(mean, axis=1).max())
    sigma = None
    if len(calcs) > 1:
        spread = np.sqrt(np.mean(np.sum((forces - mean) ** 2, axis=2), axis=0))
        sigma = float(spread.max())
    return {"max_f": max_f, "sigma": sigma, "score": sigma if sigma is not None else max_f}


def rattle_candidates(base_atoms, n_candidates, stdevs, seed):
    """Independent Gaussian-displaced copies of `base_atoms` (the reference
    cell, or its --supercell), stdevs cycled round-robin."""
    positions = build_rattled_configs(base_atoms, "cartesian", n_candidates, stdevs, seed)
    out = []
    for i, pos in enumerate(positions):
        a = base_atoms.copy()
        a.positions = pos
        out.append((a, f"rattle {stdevs[i % len(stdevs)]:g} Ang"))
    return out


def md_candidates(base_atoms, calc, n_steps, temperatures_k, stride, seed, friction_fs=0.01):
    """Frames from one short NVT trajectory per temperature, driven by
    `calc` (the first model) through stb-mlmd's own build_dynamics -- same
    verified fs/GPa unit conversions, not a re-implementation. Every
    `stride`-th frame is kept as a candidate; scoring is done afterwards by
    every committee member.

    Temperatures are deliberately meant to be higher than a production
    temperature, purely to cover more configuration space per step (a
    common MLIP training-set trick, not a claim about the material's real
    operating conditions). Several temperatures widen the sampled amplitudes.
    A poorly-fitted model can make this MD numerically unstable -- a step
    that raises (e.g. a NaN force) ends that trajectory, keeping the frames
    already collected, rather than aborting the whole scan.
    """
    out = []
    for t_index, temperature in enumerate(temperatures_k):
        atoms = base_atoms.copy()
        atoms.calc = calc
        dyn = build_dynamics(atoms, "nvt", timestep_fs=1.0, temperature_k=temperature,
                             friction_fs=friction_fs, pressure_gpa=0.0, seed=seed + t_index)
        for step in range(1, n_steps + 1):
            try:
                dyn.run(1)
            except Exception:
                break
            if step % stride == 0:
                if not np.all(np.isfinite(atoms.get_positions())):
                    break
                frame = atoms.copy()
                frame.calc = None
                out.append((frame, f"md {temperature:g} K step {step}"))
    return out


# --------------------------------------------------------------------------
# Novelty / diversity filter
# --------------------------------------------------------------------------

def environment_fingerprints(atoms, n_neighbours=FP_NEIGHBOURS, cutoff=FP_CUTOFF):
    """(n_atoms, n_neighbours) array: sorted distances (Ang) from each atom to
    its nearest neighbours within `cutoff`, padded with `cutoff`."""
    from ase.neighborlist import neighbor_list
    i_idx, dists = neighbor_list("id", atoms, cutoff)
    fp = np.full((len(atoms), n_neighbours), cutoff)
    order = np.lexsort((dists, i_idx))
    i_idx, dists = i_idx[order], dists[order]
    starts = np.searchsorted(i_idx, np.arange(len(atoms)))
    ends = np.searchsorted(i_idx, np.arange(len(atoms)), side="right")
    for a, (s, e) in enumerate(zip(starts, ends)):
        d = dists[s:min(e, s + n_neighbours)]
        fp[a, :len(d)] = d
    return fp


class EnvironmentSet:
    """Atomic-environment fingerprints of a set of configurations, kept per
    centre species (an environment is only compared to environments of the
    same element) and tagged with the configuration they came from."""

    def __init__(self):
        self.fps = {}      # symbol -> list of (n, d) arrays
        self.cfg = {}      # symbol -> list of config-id arrays
        self.n_configs = 0
        self._trees = {}

    def add(self, atoms):
        fp = environment_fingerprints(atoms)
        symbols = np.array(atoms.get_chemical_symbols())
        for sym in set(symbols):
            mask = symbols == sym
            self.fps.setdefault(sym, []).append(fp[mask])
            self.cfg.setdefault(sym, []).append(np.full(mask.sum(), self.n_configs))
            self._trees.pop(sym, None)
        self.n_configs += 1

    def _tree(self, sym):
        from scipy.spatial import cKDTree
        if sym not in self._trees:
            self._trees[sym] = (cKDTree(np.vstack(self.fps[sym])), np.concatenate(self.cfg[sym]))
        return self._trees[sym]

    def typical_spacing(self, k=64):
        """Median distance from an environment to the nearest one belonging
        to a DIFFERENT configuration -- the natural unit of 'how far apart
        two independent samples of this set are'. None with < 2 configs."""
        if self.n_configs < 2:
            return None
        nearest = []
        for sym in self.fps:
            tree, cfg = self._tree(sym)
            data = tree.data
            kk = min(k, len(data))
            dist, idx = tree.query(data, k=kk)
            dist, idx = np.atleast_2d(dist), np.atleast_2d(idx)
            for row in range(len(data)):
                other = cfg[idx[row]] != cfg[row]
                if other.any():
                    nearest.append(dist[row][other][0])
        positive = [d for d in nearest if d > 0]
        return float(np.median(positive)) if positive else None

    def novelty(self, atoms, scale):
        """Max over the atoms of `atoms` of the distance from its environment
        to the nearest one in this set (same element), in units of `scale`.
        inf when an element of `atoms` is absent from the set."""
        if self.n_configs == 0:
            return math.inf
        fp = environment_fingerprints(atoms)
        symbols = np.array(atoms.get_chemical_symbols())
        worst = 0.0
        for sym in set(symbols):
            if sym not in self.fps:
                return math.inf
            tree, _ = self._tree(sym)
            dist, _ = tree.query(fp[symbols == sym], k=1)
            worst = max(worst, float(np.max(dist)))
        return worst / scale


def read_training_configs(dirs, ref_basename):
    """ASE Atoms of every folder in `dirs` (the configs already labelled, or
    written and waiting to be), read from their own .fdf.
    Returns (atoms_list, n_unreadable)."""
    atoms_list, n_bad = [], 0
    for d in dirs:
        if not os.path.isdir(d):
            continue
        fdf = os.path.join(d, ref_basename)
        if not os.path.isfile(fdf):
            others = [f for f in sorted(glob.glob(os.path.join(d, "*.fdf")))
                      if os.path.basename(f) != "config_extra.fdf"]
            fdf = others[0] if others else None
        try:
            atoms_list.append(AseAtomsAdaptor.get_atoms(
                structure_io.to_pymatgen(structure_io.read_fdf(fdf))))
        except Exception:
            n_bad += 1
    return atoms_list, n_bad


def select_candidates(candidates, evaluations, top_k, reference, scale, min_novelty, max_force):
    """Greedy selection by decreasing score: a candidate is skipped when its
    predicted max|F| exceeds `max_force` (if set), or when its novelty --
    relative to the training set AND to the candidates already picked -- is
    below `min_novelty` (0 disables). Every pick is added to `reference`, so
    the next one must also differ from it. Returns (picked, n_force, n_redundant),
    picked = list of (candidate index, evaluation, novelty)."""
    order = sorted(range(len(candidates)), key=lambda i: evaluations[i]["score"], reverse=True)
    picked, n_force, n_redundant = [], 0, 0
    for i in order:
        if len(picked) >= top_k:
            break
        ev = evaluations[i]
        if max_force and ev["max_f"] > max_force:
            n_force += 1
            continue
        atoms = candidates[i][0]
        nov = reference.novelty(atoms, scale) if scale else math.inf
        if min_novelty > 0 and nov < min_novelty:
            n_redundant += 1
            continue
        picked.append((i, ev, nov))
        reference.add(atoms)
    return picked, n_force, n_redundant


def supercell_kgrid_block(calc_text, dims):
    """(config_extra block, note) for a supercell of the reference: the
    template's own Monkhorst-Pack divisions divided by `dims` (rounded up),
    so the supercell samples the Brillouin zone like the reference cell
    does and its energies stay on the same footing as the training set."""
    kgrid = parse_kgrid(calc_text)
    if kgrid is None:
        if "kgridcutoff" in parse_directives(calc_text):
            return "", ("template sets kgrid.Cutoff -- SIESTA sizes the supercell grid itself "
                        "(same real-space cutoff), no override needed.")
        return "", ("[WARNING] template has no k-grid (Gamma only): the supercell then samples "
                    "the Brillouin zone more densely than the reference cell, so its energies "
                    "are not on exactly the same footing as the training set.")
    new = [max(1, math.ceil(k / n)) for k, n in zip(kgrid, dims)]
    exact = all(k % n == 0 for k, n in zip(kgrid, dims))
    block = ("# Supercell k-grid (stb-mlffActiveLearning --supercell): the reference's\n"
             f"# {kgrid[0]}x{kgrid[1]}x{kgrid[2]} divided by {dims[0]}x{dims[1]}x{dims[2]}.\n"
             f"kgrid.MonkhorstPack [{new[0]} {new[1]} {new[2]}]\n\n")
    note = (f"k-grid {kgrid[0]}x{kgrid[1]}x{kgrid[2]} -> {new[0]}x{new[1]}x{new[2]} "
            + ("(equivalent sampling)" if exact else
               "([WARNING] not divisible -- rounded up, slightly denser than the reference)"))
    return block, note


def main():
    parser = argparse.ArgumentParser(
        description="Stage 4 (Active Learning) of the Custom ML Force Field workflow: uses "
                     "already fine-tuned model(s) (stb-mlffAnalysis) to screen a fresh batch "
                     "of candidate configurations (independent rattled snapshots, or frames "
                     "from short NVT trajectories driven by the model itself via "
                     "--sampling-method md, reusing stb-mlmd's own dynamics), optionally in a "
                     "--supercell of the reference, and writes the most informative ones as new "
                     "mlff_config_NNN/ folders, ready for the user to label with real SIESTA and "
                     "feed back into another round of stb-mlffAnalysis. With 2+ --model (a "
                     "committee, e.g. the same data trained with different --seed) the score is "
                     "the committee's force disagreement -- an actual uncertainty signal; with "
                     "one model it falls back to the largest predicted force, a cheap heuristic "
                     "that favours the most distorted geometries. Picks are filtered for "
                     "novelty: one whose atomic environments are all close to the training set "
                     "or to an earlier pick is skipped.",
        epilog="Example usage:\n"
               "  stb-mlffActiveLearning                    # latest round's committee from mlff_workflow.json\n"
               "  stb-mlffActiveLearning --model m_s1.model m_s2.model m_s3.model --file calc.fdf \\\n"
               "                          --n-candidates 50 --stdev 0.15 --top-k 5\n"
               "  stb-mlffActiveLearning --model m_s1.model m_s2.model m_s3.model --file calc.fdf \\\n"
               "                          --sampling-method md --md-temperature 300 600 900 \\\n"
               "                          --supercell 2 2 2",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--workflow", default=mlff_workflow.WORKFLOW_FILE, metavar="JSON",
                        help=f"Workflow state file (default: {mlff_workflow.WORKFLOW_FILE}): gives the "
                             "committee, the reference .fdf, the training set for the novelty "
                             "filter and the output directory when those options are omitted; the "
                             "picked folders are recorded in it.")
    parser.add_argument("--round", type=int, default=None, metavar="K",
                        help="Take the committee of round K from --workflow (default: the latest).")
    parser.add_argument("--model", default=None, nargs="+", metavar="PATH",
                        help="Fine-tuned MACE model file(s), overriding --workflow's committee. Two "
                             "or more form a committee: the score becomes their force disagreement "
                             "(stb-mlffAnalysis --committee N trains one). With one model the score "
                             "is its largest predicted force.")
    parser.add_argument("--file", default=None,
                        help="Reference SIESTA .fdf (default: the one in --workflow) -- basis/"
                             "pseudos/SCF settings are reused verbatim in every candidate written out.")
    parser.add_argument("--pseudo-dir", default=None,
                        help="Directory (or bundled bank) to search for <element>.psf/.psml "
                             "pseudopotentials (default: the one stb-mlff recorded in --workflow, "
                             "else the current directory).")
    parser.add_argument("--supercell", type=int, nargs=3, default=[1, 1, 1], metavar=("NA", "NB", "NC"),
                        help="Sample candidates in an NA x NB x NC supercell of --file (default: "
                             "1 1 1, the reference cell itself). A small training cell only "
                             "constrains the interactions its own periodicity allows -- e.g. an "
                             "8-atom diamond cell sees the Gamma and X phonons but not L -- and a "
                             "supercell adds the rest. Each folder gets the reference k-grid "
                             "divided by the supercell size in its config_extra.fdf. SIESTA cost "
                             "grows steeply with the atom count.")
    parser.add_argument("--sampling-method", choices=["rattle", "md"], default="rattle",
                        help="How to generate candidate configurations (default: rattle):\n"
                             "  rattle - independent, static Gaussian-displaced snapshots "
                             "(cheap: 1 evaluation per candidate and model).\n"
                             "  md     - short NVT trajectories with the (first) model itself "
                             "(stb-mlmd's own dynamics, reused directly); more physically "
                             "representative of the configuration space the system actually "
                             "explores, at the cost of running real MD steps.")
    parser.add_argument("--n-candidates", type=int, default=50,
                        help="Number of candidate configurations to screen (--sampling-method "
                             "rattle only, default: 50).")
    parser.add_argument("--stdev", type=float, nargs="+", default=[0.15],
                        help="Rattling standard deviation(s) in Angstrom for the candidate "
                             "batch (--sampling-method rattle only, default: 0.15 -- typically "
                             "larger than the amplitude(s) used for the original training set, "
                             "to probe configuration space not yet covered).")
    parser.add_argument("--md-steps", type=int, default=500,
                        help="NVT steps per temperature (--sampling-method md only, default: 500).")
    parser.add_argument("--md-temperature", type=float, nargs="+", default=[800.0],
                        help="NVT temperature(s), K -- one trajectory each (--sampling-method md "
                             "only, default: 800 -- deliberately higher than a typical production "
                             "temperature, purely to cover more configuration space per step; "
                             "not a claim about the material's real operating conditions).")
    parser.add_argument("--stride", type=int, default=10,
                        help="Keep one frame every Nth MD step as a candidate (--sampling-method md "
                             "only, default: 10).")
    parser.add_argument("--top-k", type=int, default=5,
                        help="Number of candidates to write out (default: 5).")
    parser.add_argument("--min-novelty", type=float, default=DEFAULT_MIN_NOVELTY,
                        help="Skip a candidate whose atomic environments are all closer than this "
                             "to the training set or to an earlier pick (default: "
                             f"{DEFAULT_MIN_NOVELTY:g}; 0 disables). Measured on a model-independent "
                             f"fingerprint (sorted distances to the {FP_NEIGHBOURS} nearest neighbours) "
                             "in units of the training set's own typical spacing -- 1 means as far "
                             "as two independent training configs are from each other.")
    parser.add_argument("--train-path", default=None, metavar="GLOB",
                        help="Folders the novelty filter compares against (default: the training "
                             "configs in --workflow, else '<output-dir>/mlff_config_*' -- labelled "
                             "or not). Quote the pattern.")
    parser.add_argument("--max-force", type=float, default=0.0, metavar="EV_ANG",
                        help="Skip candidates whose predicted max|F| exceeds this, eV/Ang -- "
                             "geometries too distorted to matter for the intended use (default: 0, "
                             "no limit).")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for rattling / initial MD velocities (default: 42).")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"],
                        help="Device to run the model(s) on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Directory under which new mlff_config_NNN/ folders are created "
                             "(default: --workflow's directory, else the current one -- continues "
                             "the numbering of any already present there).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlffActiveLearning {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Active Learning)",
            "Suggest new SIESTA configurations to label, using the fine-tuned model(s)",
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

    try:
        state = mlff_workflow.load(args.workflow)
    except ValueError as e:
        fail(str(e))
    wf_base = mlff_workflow.base_dir(args.workflow)
    round_entry = None
    if args.model is None or args.file is None:
        if state is None:
            fail(f"'{args.workflow}' not found -- pass --model and --file, or run the earlier "
                 "stages (stb-mlff, stb-mlffAnalysis) to create it.")
    if args.model is None:
        round_entry = mlff_workflow.get_round(state, args.round)
        if round_entry is None:
            fail(f"No training round{' ' + str(args.round) if args.round else ''} in '{args.workflow}' "
                 "-- run stb-mlffAnalysis (4.17.2) first, or pass --model.")
        args.model = mlff_workflow.committee_models(round_entry, wf_base)
        if not args.model:
            fail(f"Round {round_entry['round']} in '{args.workflow}' has no model file.")
    elif state is not None and args.round is not None:
        round_entry = mlff_workflow.get_round(state, args.round)
    if args.file is None:
        ref = state.get("reference", {}).get("file")
        if not ref:
            fail(f"'{args.workflow}' records no reference .fdf -- pass --file.")
        args.file = mlff_workflow.resolve(ref, wf_base)
    if args.output_dir is None:
        args.output_dir = os.path.relpath(wf_base) if state is not None else "."
    if args.pseudo_dir is None:
        recorded = (state or {}).get("reference", {}).get("pseudo_dir")
        if recorded:
            candidate = mlff_workflow.resolve(recorded, wf_base)
            args.pseudo_dir = candidate if os.path.isdir(candidate) else recorded
        else:
            args.pseudo_dir = "."
    proposer_round = round_entry["round"] if round_entry else (
        mlff_workflow.get_round(state)["round"] if state and state["rounds"] else None)

    committee = len(args.model) > 1
    dims = list(args.supercell)
    is_supercell = dims != [1, 1, 1]
    if args.train_path:
        train_dirs = sorted(d for d in glob.glob(args.train_path) if os.path.isdir(d))
        train_pattern = args.train_path
    elif state is not None and mlff_workflow.folders(state, wf_base, "train"):
        train_dirs = [d for d in mlff_workflow.folders(state, wf_base, "train") if os.path.isdir(d)]
        train_pattern = f"{args.workflow} (training configs)"
    else:
        train_pattern = os.path.join(args.output_dir, "mlff_config_*")
        train_dirs = sorted(d for d in glob.glob(train_pattern) if os.path.isdir(d))

    print_dual(color_text("===== STB-MLFFACTIVELEARNING REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if round_entry is not None:
        print_dual(f"Committee from    : {args.workflow}, round {round_entry['round']}", f_out)
    print_dual(f"Model(s)          : {args.model[0]}", f_out)
    for m in args.model[1:]:
        print_dual(f"                    {m}", f_out)
    print_dual("Score             : " + (
        f"committee force disagreement ({len(args.model)} models; max over atoms of the "
        "per-atom spread)" if committee else
        "largest predicted force (single model -- pass 2+ --model for a committee)"), f_out)
    print_dual(f"Reference file    : {args.file}", f_out)
    print_dual(f"Supercell         : {dims[0]} x {dims[1]} x {dims[2]}", f_out)
    print_dual(f"Sampling method   : {args.sampling_method}", f_out)
    if args.sampling_method == "rattle":
        print_dual(f"Candidates/top-k  : {args.n_candidates} / {args.top_k} (stdev "
                   f"{' '.join(f'{s:g}' for s in args.stdev)})", f_out)
    else:
        print_dual(f"MD steps/temp     : {args.md_steps} per temperature / "
                   f"{', '.join(f'{t:g}' for t in args.md_temperature)} K "
                   f"(stride {args.stride}) / top-k {args.top_k}", f_out)
    print_dual(f"Novelty filter    : " + (f"min {args.min_novelty:g} vs '{train_pattern}' + earlier picks"
                                          if args.min_novelty > 0 else "off"), f_out)
    print_dual(f"Max force         : " + (f"{args.max_force:g} eV/Ang" if args.max_force else "no limit"), f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING REFERENCE STRUCTURE & MODEL(S)", f_out)
    if not os.path.isfile(args.file):
        fail(f"Reference file '{args.file}' not found.")
    for m in args.model:
        if not os.path.isfile(m):
            fail(f"Model file '{m}' not found.")
    for name, value in (("--top-k", args.top_k), ("--n-candidates", args.n_candidates),
                        ("--md-steps", args.md_steps), ("--stride", args.stride)):
        if value < 1:
            fail(f"{name} must be >= 1 (got {value}).")
    if min(dims) < 1:
        fail(f"--supercell dimensions must be >= 1 (got {dims}).")
    try:
        ref_structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")
    unique_elements = sorted(ref_structure.species)
    print_dual(f"[OK] Read {len(ref_structure.atoms)} atom(s), species: {', '.join(unique_elements)}", f_out)

    pseudos_to_copy, missing = get_required_pseudos(unique_elements, args.pseudo_dir)
    if missing:
        fail(f"Missing pseudopotential(s) for: {', '.join(missing)} (searched '{args.pseudo_dir}').")

    ase_atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(ref_structure))
    ref_basename = os.path.basename(args.file)
    template = args.file
    extra_block = ""
    tmp_dir = None
    if is_supercell:
        tmp_dir = tempfile.mkdtemp(prefix="stb_mlffal_")
        template = os.path.join(tmp_dir, ref_basename)
        try:
            n_super = structure_io.replicate_fdf_cell(args.file, dims, template)
        except Exception as e:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            fail(f"Could not build the {dims} supercell of '{args.file}': {e}")
        with open(args.file) as f:
            extra_block, k_note = supercell_kgrid_block(f.read(), dims)
        print_dual(f"[OK] Supercell {dims[0]}x{dims[1]}x{dims[2]}: {n_super} atoms "
                   f"(reference {len(ase_atoms)})", f_out)
        print_dual(f"     {k_note}", f_out)
    base_atoms = ase_atoms.repeat(dims)

    try:
        from mace.calculators import MACECalculator
    except ImportError:
        fail("mace-torch not found. This tool needs the optional 'ml' extra: pip install stb_suite[ml]")
    calcs = [MACECalculator(model_paths=[m], device=args.device) for m in args.model]
    print_dual(f"[OK] Loaded {len(calcs)} model(s)", f_out)

    reference = EnvironmentSet()
    if args.min_novelty > 0:
        train_atoms, n_bad = read_training_configs(train_dirs, ref_basename)
        for a in train_atoms:
            reference.add(a)
        print_dual(f"[OK] Novelty reference: {len(train_atoms)} config(s) from '{train_pattern}'"
                   + (f" ({n_bad} unreadable, skipped)" if n_bad else ""), f_out)

    print_section("[2] SCREENING CANDIDATES", f_out)
    if args.sampling_method == "rattle":
        candidates = rattle_candidates(base_atoms, args.n_candidates, args.stdev, args.seed)
    else:
        candidates = md_candidates(base_atoms, calcs[0], args.md_steps, args.md_temperature,
                                   args.stride, args.seed)
    n_attempted = (args.n_candidates if args.sampling_method == "rattle"
                   else len(args.md_temperature) * (args.md_steps // args.stride))
    evaluated = []
    for cand in candidates:
        ev = evaluate_candidate(cand[0], calcs)
        if ev is not None:
            evaluated.append((cand, ev))
    if not evaluated:
        if tmp_dir:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        fail("None of the candidate configurations could be evaluated by the model(s).")
    candidates = [c for c, _ in evaluated]
    evaluations = [e for _, e in evaluated]
    print_dual(f"[OK] Screened {len(candidates)}/{n_attempted} candidate(s).", f_out)
    max_fs = [e["max_f"] for e in evaluations]
    print_dual(f"Predicted max|F| range: {min(max_fs):.3f} - {max(max_fs):.3f} eV/Ang", f_out)
    if committee:
        sig = [e["sigma"] for e in evaluations]
        print_dual(f"Committee spread range: {min(sig):.3f} - {max(sig):.3f} eV/Ang", f_out)

    print_section("[3] SELECTING & WRITING CONFIGURATIONS", f_out)
    scale = None
    if args.min_novelty > 0:
        scale = reference.typical_spacing()
        scale_source = "training set"
        if scale is None:
            pool = EnvironmentSet()
            for a, _ in candidates:
                pool.add(a)
            scale = pool.typical_spacing()
            scale_source = "candidate pool (training set has < 2 configs)"
        if scale:
            print_dual(f"Novelty unit      : {scale:.4f} Ang (typical spacing of the {scale_source})", f_out)
        else:
            print_dual(color_text("[WARNING] Could not size the novelty unit -- filter disabled.", 'yellow'), f_out)
    picked, n_force, n_redundant = select_candidates(
        candidates, evaluations, args.top_k, reference, scale, args.min_novelty, args.max_force)
    if n_force:
        print_dual(f"Skipped (max|F| > {args.max_force:g} eV/Ang) : {n_force}", f_out)
    if n_redundant:
        print_dual(f"Skipped (novelty < {args.min_novelty:g})      : {n_redundant}", f_out)
    if len(picked) < args.top_k:
        print_dual(color_text(
            f"[WARNING] Only {len(picked)} of the requested {args.top_k} candidate(s) passed the "
            "filters -- the training set already covers this sampling; try larger --stdev, higher "
            "--md-temperature, a --supercell, or a lower --min-novelty.", 'yellow'), f_out)

    os.makedirs(args.output_dir, exist_ok=True)
    start_index = next_config_index(args.output_dir)
    folders = []
    picked_entries = []
    for offset, (i, ev, nov) in enumerate(picked):
        atoms, origin = candidates[i]
        positions = (atoms.get_scaled_positions(wrap=False) if ref_structure.coord_format == "fractional"
                     else atoms.get_positions())
        folder = write_config_folder(
            args.output_dir, start_index + offset, template, ref_basename, positions, None,
            pseudos_to_copy, extra_block=extra_block)
        folders.append(folder)
        nov_txt = "n/a" if math.isinf(nov) else f"{nov:.2f}"
        picked_entries.append({"folder": folder, "set": "train",
                               "kind": f"al-{args.sampling_method}",
                               "origin": origin + (f", supercell {dims[0]}x{dims[1]}x{dims[2]}" if is_supercell else ""),
                               "stage": "active-learning", "round": proposer_round,
                               "n_atoms": len(atoms), "score": ev["score"], "max_f": ev["max_f"],
                               "novelty": None if math.isinf(nov) else nov})
        print_dual(f"  - {folder}  (score {ev['score']:.3f} eV/Ang | max|F| {ev['max_f']:.3f} eV/Ang"
                   f" | novelty {nov_txt} | {origin})", f_out)
    if tmp_dir:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    if state is not None:
        mlff_workflow.add_configs(state, picked_entries, wf_base)
        target = mlff_workflow.get_round(state, proposer_round) if proposer_round else None
        if target is not None:
            target.setdefault("active_learning", []).append({
                "time": datetime.now().isoformat(timespec="seconds"),
                "sampling": args.sampling_method, "supercell": dims,
                "committee_size": len(args.model),
                "picked": [mlff_workflow.rel(f, wf_base) for f in folders]})
        mlff_workflow.record(state, "active-learning", "stb-mlffActiveLearning", sys.argv[1:])
        mlff_workflow.save(state, args.workflow)

    print_section("[4] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Suggested configs : {len(folders)}", f_out)
    if state is not None:
        print_dual(f"Workflow file     : {args.workflow} ({len(folders)} config(s) added to the training set)", f_out)
    if not committee:
        print_dual(color_text(
            "[NOTE] Single model: the score is its largest predicted force, which favours the "
            "most distorted geometries rather than the ones it gets wrong. For an uncertainty "
            "score, train 2-3 more models on the same configs with stb-mlffAnalysis --seed N "
            "--work-dir <dir> and pass them all to --model.", 'yellow'), f_out)
    print_dual("Next steps        : run SIESTA in each folder above, then train the next "
               "round with stb-mlffAnalysis (4.17.2) -- it takes the configs from "
               "mlff_workflow.json -- and score it with stb-mlffValidate (4.17.3).", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Candidates suggested -- label them, then retrain.\n", 'bold'))


if __name__ == "__main__":
    main()
