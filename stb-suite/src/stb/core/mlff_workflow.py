"""The `mlff_workflow.json` state file shared by the Custom ML Force Field stages.

One file per project (next to the `mlff_config_*` folders), read and updated by
each stage, so a stage never has to be told again what an earlier one did:

- stb-mlff (4.17.1, Prep)            -- creates it: reference .fdf, every config
                                        folder written (training and test sets),
                                        with its origin.
- stb-mlffAnalysis (4.17.2, Train)   -- adds a round: the configs it trained on,
                                        the committee (one entry per seed: model
                                        path, validation metrics), settings.
- stb-mlffValidate (4.17.3, Validate)-- stores each round's test-set metrics.
- stb-mlffActiveLearning (4.17.4)    -- reads the latest round's committee and the
                                        reference; adds the picked folders.

Paths are stored relative to the file's own directory, so a project folder can
be moved or copied to a cluster as a whole. Writes go through a temporary file
and os.replace, so an interrupted stage never leaves a truncated file.

Schema (version 1)::

    {"format": "stb-mlff-workflow", "version": 1, "created": ..., "updated": ...,
     "reference": {"file", "n_atoms", "species"},
     "configs":   [{"folder", "set": "train"|"test", "kind", "origin",
                    "stage", "round", "n_atoms"}],
     "rounds":    [{"round", "dir", "created", "configs": [folders],
                    "committee": [{"seed", "dir", "model", "n_train", "n_valid",
                                   "valid": {metric: value}}],
                    "settings": {...}, "validation": {...},
                    "active_learning": [{...}]}],
     "history":   [{"time", "stage", "tool", "argv"}]}
"""

import os
import json
import shutil
from datetime import datetime

WORKFLOW_FILE = "mlff_workflow.json"
FORMAT = "stb-mlff-workflow"
VERSION = 1


def _now():
    return datetime.now().isoformat(timespec="seconds")


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
