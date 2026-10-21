"""The `mlff_workflow.json` state file shared by the Custom ML Force Field stages.

One file per project, inside the run folder (`mlff_runs/` by default), read and
updated by each stage, so a stage never has to be told again what an earlier
one did. Layout of the run folder:

    configs_round_01/mlff_config_001..   Stage 1's training configs (train round 1)
    round_01/                            round 1's committee (Stage 2)
    configs_round_02/mlff_config_011..   Stage 4's picks from round 1 (train round 2)
    round_02/ ...
    mlff_test/test_001..                 Stage 1's held-out test set
    md_check/rNN_<T>K/frame_000..        Stage 3 --md-check frames

The mlff_config_NNN numbers stay unique project-wide; an older flat layout
(mlff_config_* directly in the run folder) is still read (find_config_dirs).

- stb-mlff (4.17.1, Prep)            -- creates it: reference .fdf, every config
                                        folder written (training and test sets),
                                        with its origin.
- stb-mlffAnalysis (4.17.2, Train)   -- adds a round: the configs it trained on,
                                        the committee (one entry per seed: model
                                        path, validation metrics), settings.
- stb-mlffValidate (4.17.3, Validate)-- stores each round's test-set metrics; with
                                        --md-check also writes MD-frame folders
                                        (set "md_check").
- stb-mlffActiveLearning (4.17.4)    -- reads a round's committee, its validation
                                        (calibration, verdict) and the reference;
                                        adds the picked folders (set "train") and
                                        an "active_learning" entry to the round.

Shared thresholds: UNPHYSICAL_FMAX / UNPHYSICAL_DE below -- one definition of an
unphysical configuration for every stage.

Paths are stored relative to the file's own directory, so a project folder can
be moved or copied to a cluster as a whole. Writes go through a temporary file
and os.replace, so an interrupted stage never leaves a truncated file.

Schema (version 1)::

    {"format": "stb-mlff-workflow", "version": 1, "created": ..., "updated": ...,
     "reference": {"file", "calc", "n_atoms", "species", "pseudo_dir"},
                  -- "file" is the structure .fdf (or the one self-contained
                  .fdf); "calc" the calc.fdf template that %includes it,
                  absent for a self-contained file.
     "configs":   [{"folder", "set": "train"|"test"|"md_check", "kind", "origin",
                    "stage", "round", "n_atoms",
                    "geometry": {volume, dV_pct, u_rms, u_max, d_min, novelty}  (Stage 1),
                    "ml_screen": {model, f_max, dE_meV_atom}                    (Stage 1),
                    "score", "max_f", "novelty", "d_min"                        (Stage 4 picks),
                    "md": {run, time_fs, temperature, t_inst, round, model}     (md_check)}],
                  -- kinds: rattle, strain (Stage 1); al-rattle, al-md, al-strain
                  (Stage 4); md_check (Stage 3 frames, also after --md-check-to-train).
     "rounds":    [{"round", "dir", "created", "configs": [folders], "n_usable",
                    "committee": [{"seed", "dir", "model", "n_train", "n_valid",
                                   "valid": {metric: value}, ...}],
                    "settings": {...}, "init_from",
                    "training_errors": {out_of_sample, foundation, gain_pct,
                                        committee, cv_folds},          (Stage 2)
                    "model_checks": {stability: [...], speed: {...}},  (Stage 2)
                    "validation": {overall, groups, members, calibration
                                   {error_over_spread, ...}, phonons, per_config,
                                   elements, eos, md_check, foundation, verdict,
                                   learning_curve, compare},           (Stage 3)
                    "active_learning": [{sampling, settings, score, selection,
                                         pool, skipped, picked, another_round,
                                         auto_sampling, ...}]}],        (Stage 4)
     "history":   [{"time", "stage", "tool", "argv"}]}
"""

import os
import json
import shutil
from datetime import datetime

# One definition of an 'unphysical' configuration for every stage: Prep flags
# it from the foundation-model screen, Train warns about SIESTA labels beyond
# it, Active learning skips candidates above it by default (--max-force) and
# --auto-sampling keeps its median below it. Far above anything MD visits, such
# configs dominate the force loss.
UNPHYSICAL_FMAX = 15.0       # eV/Ang, largest atomic force
UNPHYSICAL_DE = 1000.0       # meV/atom above the reference

WORKFLOW_FILE = "mlff_workflow.json"
RUN_DIR = "mlff_runs"   # stb-mlff's default --output-dir, like the other workflows' *_runs
FORMAT = "stb-mlff-workflow"
VERSION = 1


def _now():
    return datetime.now().isoformat(timespec="seconds")


def default_path(explicit=None):
    """The --workflow file a stage uses: `explicit` when given, else
    mlff_runs/mlff_workflow.json -- or ./mlff_workflow.json when only that one
    exists (a stage run from inside the run folder, or an older layout)."""
    if explicit:
        return explicit
    inside = os.path.join(RUN_DIR, WORKFLOW_FILE)
    if not os.path.isfile(inside) and os.path.isfile(WORKFLOW_FILE):
        return WORKFLOW_FILE
    return inside


def new_state(reference_file=None, n_atoms=None, species=None, base_dir="."):
    state = {"format": FORMAT, "version": VERSION, "created": _now(), "updated": _now(),
             "reference": {}, "configs": [], "rounds": [], "history": []}
    if reference_file:
        state["reference"] = {"file": rel(reference_file, base_dir), "n_atoms": n_atoms,
                              "species": list(species or [])}
    return state


def load(path):
    """The state dict, or None when `path` does not exist. Raises ValueError
    for a file that is not a workflow state (wrong format / unreadable)."""
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            state = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise ValueError(f"'{path}' is not a readable JSON file: {e}")
    if not isinstance(state, dict) or state.get("format") != FORMAT:
        raise ValueError(f"'{path}' is not an {FORMAT} file.")
    if state.get("version", 0) > VERSION:
        raise ValueError(f"'{path}' was written by a newer stb (version {state['version']}).")
    for key, default in (("reference", {}), ("configs", []), ("rounds", []), ("history", [])):
        state.setdefault(key, default)
    return state


def save(state, path):
    state["updated"] = _now()
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(state, f, indent=1)
        f.write("\n")
    os.replace(tmp, path)


def backup(path):
    """Copies an existing file to <path>.bak (returns that path, or None)."""
    if not os.path.isfile(path):
        return None
    bak = path + ".bak"
    shutil.copy2(path, bak)
    return bak


def base_dir(path):
    return os.path.dirname(os.path.abspath(path))


def rel(path, base):
    return os.path.relpath(os.path.abspath(path), os.path.abspath(base))


def resolve(relative, base):
    """A stored path made usable from the current directory (relative to it)."""
    return os.path.relpath(os.path.normpath(os.path.join(base, relative)))


def record(state, stage, tool, argv):
    state["history"].append({"time": _now(), "stage": stage, "tool": tool, "argv": list(argv)})


def add_configs(state, entries, base):
    """Adds (or replaces, by folder) config entries. Each entry is a dict with
    at least 'folder' (any path; stored relative to `base`) and 'set'."""
    by_folder = {c["folder"]: i for i, c in enumerate(state["configs"])}
    for entry in entries:
        e = dict(entry)
        e["folder"] = rel(e["folder"], base)
        if e["folder"] in by_folder:
            state["configs"][by_folder[e["folder"]]] = e
        else:
            by_folder[e["folder"]] = len(state["configs"])
            state["configs"].append(e)


def folders(state, base, which="train"):
    """Absolute paths of the config folders of one set ('train' or 'test')."""
    return [resolve(c["folder"], base) for c in state["configs"] if c.get("set") == which]


def config_entry(state, folder, base):
    key = rel(folder, base)
    return next((c for c in state["configs"] if c["folder"] == key), None)


CONFIGS_PREFIX = "configs_round_"   # <run>/configs_round_NN/mlff_config_MMM: the configs that train round NN
CONFIG_PREFIX = "mlff_config"


def configs_dir(base, round_no):
    """Folder for the configs that will train round `round_no` (Stage 1: round 1;
    Stage 4: the round after the latest), next to that round's round_NN/."""
    return os.path.join(base, f"{CONFIGS_PREFIX}{round_no:02d}")


def config_globs(base):
    """Glob patterns of every training-config folder under `base`: the per-round
    layout and the older flat one (mlff_config_* directly in the run folder)."""
    return [os.path.join(base, f"{CONFIGS_PREFIX}*", f"{CONFIG_PREFIX}_*"),
            os.path.join(base, f"{CONFIG_PREFIX}_*")]


def find_config_dirs(base):
    """Every training-config folder under `base`, either layout, in numbering order."""
    import glob
    import re
    found = {d for pat in config_globs(base) for d in glob.glob(pat) if os.path.isdir(d)}

    def key(d):
        m = re.search(r"_(\d+)$", os.path.basename(os.path.normpath(d)))
        return (int(m.group(1)) if m else 10 ** 9, d)
    return sorted(found, key=key)


def next_config_index(base):
    """1 + the highest mlff_config_NNN number under `base` (either layout), so
    a new batch never reuses a name -- the number stays unique project-wide."""
    import re
    highest = 0
    for d in find_config_dirs(base):
        m = re.match(rf"^{CONFIG_PREFIX}_(\d+)$", os.path.basename(os.path.normpath(d)))
        if m:
            highest = max(highest, int(m.group(1)))
    return highest + 1


def siesta_finished(folder):
    """True when `folder` holds a finished SIESTA single point (.FA and .XV written)."""
    import glob
    return bool(glob.glob(os.path.join(folder, "*.FA"))) and bool(glob.glob(os.path.join(folder, "*.XV")))


def pending_by_dir(folders):
    """[(parent folder, n folders, n without a finished SIESTA run)] for `folders`,
    grouped by the folder they live in (e.g. configs_round_02/, mlff_test/) -- what
    the reports print as 'run SIESTA in ...'."""
    groups = {}
    for f in folders:
        parent = os.path.dirname(os.path.normpath(f)) or "."
        g = groups.setdefault(parent, [0, 0])
        g[0] += 1
        if not siesta_finished(f):
            g[1] += 1
    return [(p, n, k) for p, (n, k) in sorted(groups.items())]


def siesta_loop(parent, folders, main_fdf="calc.fdf", label="siesta"):
    """Shell loop running SIESTA in `folders` (all inside `parent`)."""
    names = sorted({os.path.basename(os.path.normpath(f)).rsplit("_", 1)[0] for f in folders})
    pat = " ".join(os.path.join(os.path.relpath(parent), f"{n}_*") for n in names) or os.path.relpath(parent)
    return f'for d in {pat}; do (cd "$d" && siesta < {main_fdf} > {label}.out); done'


def next_round_number(state):
    return 1 + max((r["round"] for r in state["rounds"]), default=0)


def get_round(state, number=None):
    """Round `number`, or the latest one when None; None if there is none."""
    if not state["rounds"]:
        return None
    if number is None:
        return max(state["rounds"], key=lambda r: r["round"])
    return next((r for r in state["rounds"] if r["round"] == number), None)


def add_round(state, entry):
    state["rounds"] = [r for r in state["rounds"] if r["round"] != entry["round"]]
    state["rounds"].append(entry)
    state["rounds"].sort(key=lambda r: r["round"])


def committee_models(round_entry, base):
    return [resolve(m["model"], base) for m in round_entry.get("committee", []) if m.get("model")]
