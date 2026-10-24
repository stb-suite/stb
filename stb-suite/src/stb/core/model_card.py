"""The model card: the JSON file written next to every model stb-mlffAnalysis trains.

A `.model` file holds weights and nothing about where they come from. The card
holds the rest -- what SIESTA was told (level of theory), which elements,
volumes, distances and forces the training data covered, how the model was
trained and how well it scored, its energy scale, whether the labels already
include DFT-D3 -- so every tool that runs the model can check that it is being
used inside the range it was fitted for, and say so when it is not.

The contract of the suite: **a fine-tuned model is always handed over together
with its card** (`--custom-model model.model --custom-model-card card.json`).
`setup_from_args()` is the one function a tool calls for it; it loads and
verifies the card (the model file's SHA-256 must be one of the card's) and
registers the pair, and `mace_relax.get_calculator()` -- the single place every
tool loads a model -- refuses a custom model file that has no registered card.
It also uses the card: a model whose labels include DFT-D3 is never given a
second dispersion correction, and the calculator it returns warns (once per
message) when a structure has an element the model never saw, a cell or a bond
length outside the trained range, or forces far beyond the training forces.

Card layout (version 1; every section may be extended, readers must ignore
unknown keys)::

    {"format": "stb-mlff-model-card", "version": 1, "created", "stb_version",
     "model":   {"name", "file", "files": [{"file", "kind", "sha256", "size_bytes"}],
                 "seed", "member", "architecture": {r_max, n_params, num_interactions,
                 atomic_numbers, ...}, "foundation": {...}, "init_from"},
     "committee": {"size", "round", "members": [{"seed", "card", "model"}]},
     "elements": ["Si", ...],
     "energy": {"scale", "e0_eV": {"Si": ...}, "e0_source", "note"},
     "level_of_theory": {"code", "siesta_version", "settings": {...}, "dftd3": {...},
                         "pseudopotentials": {...}, "reference_files": {...}},
     "reference_structure": {formula, n_atoms, cell, dimensionality, ...},
     "training_data": {n_configs, n_train, n_valid, kinds, elements_count, labels, ...,
                       "configs": [{folder, kind, set}]},
     "domain":  {the ranges used by check_structure()},
     "recipe":  {epochs, batch_size, lr, e0s, loss, stress_weight, valid_fraction, ...},
     "metrics": {train, valid, out_of_sample, foundation, gain_pct, committee, history},
     "checks":  {stability, speed},
     "validation": {filled by stb-mlffValidate: test-set scores, E(V), phonons, verdict},
     "usage":   {recommended_dtype, notes: [...], how_to_use},
     "provenance": {workflow, round, argv, software}}
"""

import os
import sys
import json
import hashlib
from datetime import datetime

import numpy as np

FORMAT = "stb-mlff-model-card"
VERSION = 1
CARD_FILE = "mlff_model_card.json"

# Tolerances of check_structure(): how far outside the trained range a structure
# may be before the card warns. A cell is "outside" when its area/volume per atom
# differs from the trained extremes by more than MEASURE_TOL (relative); a bond is
# "closer than trained" below D_MIN_FACTOR x the shortest trained distance; forces
# are "beyond training" above FORCE_FACTOR x the largest training force.
MEASURE_TOL = 0.02
D_MIN_FACTOR = 0.95
FORCE_FACTOR = 1.5
VACUUM_GAP = 5.0       # Ang: an empty gap this wide marks a non-periodic axis

# model files (resolved path) -> card, filled by register()
_REGISTRY = {}
_CARD_PATHS = {}
_WARNED = set()


class ModelCardError(ValueError):
    """The card is missing, unreadable, not a card, or does not belong to the model."""


# ---- files ----------------------------------------------------------------

def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _now():
    return datetime.now().isoformat(timespec="seconds")


def write_card(card, path):
    """Writes `card` to `path` (temporary file + os.replace, like the workflow file)."""
    card["updated"] = _now()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(_jsonable(card), f, indent=1)
        f.write("\n")
    os.replace(tmp, path)


def _jsonable(obj):
    """numpy scalars/arrays -> plain Python, non-finite floats -> None."""
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return _jsonable(obj.item())
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def load_card(path):
    """The card dict at `path`; raises ModelCardError if it is not a readable model card."""
    if not os.path.isfile(path):
        raise ModelCardError(f"model card '{path}' not found.")
    try:
        with open(path) as f:
            card = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise ModelCardError(f"'{path}' is not a readable JSON file: {e}")
    if not isinstance(card, dict) or card.get("format") != FORMAT:
        raise ModelCardError(f"'{path}' is not a {FORMAT} file (the card stb-mlffAnalysis writes next to "
                             f"each model as {CARD_FILE}).")
    if card.get("version", 0) > VERSION:
        raise ModelCardError(f"'{path}' was written by a newer stb (card version {card['version']}).")
    return card


def verify_model(card, model_path):
    """Raises ModelCardError unless `model_path` is one of the model files the card describes
    (by SHA-256): a card of another model, or of the same model retrained, is refused."""
    known = {f.get("sha256"): f.get("file") for f in (card.get("model") or {}).get("files", [])}
    digest = sha256_file(model_path)
    if digest not in known:
        names = ", ".join(sorted(n for n in known.values() if n)) or "none"
        raise ModelCardError(
            f"the model file '{model_path}' is not the one this card describes (SHA-256 mismatch; the "
            f"card is for: {names}). Pass the card that was written next to this model.")
    return known[digest]


def species(card):
    return list(card.get("elements") or [])


# ---- geometry measures shared by training and checking ----------------------

def periodic_measure(atoms, vacuum_gap=VACUUM_GAP):
    """(name, value, vacuum_flags) of the size of the cell per atom along its PERIODIC axes only:
    volume/atom (3D), area/atom (2D slab or layer), length/atom (1D), or (None, None, flags) for a
    molecule -- a vacuum-padded cell would otherwise make the number meaningless. Same function
    for the training data and for the structure being checked, so they are comparable."""
    from stb.core.kspace import detect_vacuum_axes
    cell = np.asarray(atoms.cell)
    try:
        vac = detect_vacuum_axes(atoms.get_scaled_positions(wrap=True), cell, vacuum_gap)
    except Exception:
        vac = [False, False, False]
    per = [i for i in range(3) if not vac[i]]
    n = len(atoms)
    if len(per) == 3:
        return "volume per atom (Ang^3)", abs(float(np.linalg.det(cell))) / n, vac
    if len(per) == 2:
        return "area per atom (Ang^2)", float(np.linalg.norm(np.cross(cell[per[0]], cell[per[1]]))) / n, vac
    if len(per) == 1:
        return "length per atom (Ang)", float(np.linalg.norm(cell[per[0]])) / n, vac
    return None, None, vac


def shortest_distance(atoms, cutoff=6.0):
    """Shortest interatomic distance (Ang) with periodic images, None for a single isolated atom."""
    from ase.neighborlist import neighbor_list
    try:
        d = neighbor_list("d", atoms, cutoff)
    except Exception:
        return None
    return float(d.min()) if len(d) else None


def _range(values):
    v = [float(x) for x in values if x is not None and np.isfinite(x)]
    return {"min": min(v), "max": max(v)} if v else None


def training_domain(atoms_list):
    """The ranges the training data covered, from ASE Atoms that carry REF_energy / REF_forces
    (and optionally REF_stress): what check_structure() compares a new structure with."""
    symbols = sorted({s for a in atoms_list for s in a.get_chemical_symbols()},
                     key=lambda s: _z(s))
    names, measures, dmins, vacs = set(), [], [], []
    for a in atoms_list:
        name, val, vac = periodic_measure(a)
        if name:
            names.add(name)
            measures.append(val)
        vacs.append(sum(vac))
        dmins.append(shortest_distance(a))
    e_atom = [a.info["REF_energy"] / len(a) for a in atoms_list]
    fmax = [float(np.linalg.norm(a.arrays["REF_forces"], axis=1).max()) for a in atoms_list]
    press = []
    for a in atoms_list:
        s_ = a.info.get("REF_stress")
        if s_ is not None:
            s_ = np.asarray(s_, float).reshape(3, 3)
            press.append(-float(np.trace(s_)) / 3 * 160.21766208)
    vac_count = max(set(vacs), key=vacs.count) if vacs else 0
    dom = {"elements": symbols,
           "measure_name": sorted(names)[0] if len(names) == 1 else None,
           "measure": _range(measures) if len(names) == 1 else None,
           "n_vacuum_axes": vac_count,
           "dimensionality": {0: "3D", 1: "2D", 2: "1D", 3: "0D"}.get(vac_count),
           "d_min_ang": _range(dmins),
           "n_atoms": {"min": min(len(a) for a in atoms_list), "max": max(len(a) for a in atoms_list)},
           "energy_per_atom_eV": _range(e_atom),
           "force_max_eV_per_ang": _range(fmax),
           "pressure_GPa": _range(press)}
    if len(names) > 1:
        dom["warning"] = "training configs have different dimensionality; the cell-size range is not given"
    return dom


def _z(symbol):
    from ase.data import atomic_numbers
    return atomic_numbers.get(symbol, 999)


# ---- using a card ---------------------------------------------------------

def check_structure(card, atoms):
    """List of messages (empty = inside the trained range) about running the model on `atoms`."""
    msgs = []
    dom = card.get("domain") or {}
    trained = set(dom.get("elements") or card.get("elements") or [])
    present = sorted(set(atoms.get_chemical_symbols()), key=_z)
    missing = [s for s in present if trained and s not in trained]
    if missing:
        msgs.append(f"element(s) {', '.join(missing)} are NOT in the training data (trained on "
                    f"{', '.join(sorted(trained, key=_z))}): the model has no SIESTA-fitted parameters for them.")
    name, val, vac = periodic_measure(atoms)
    dname, rng = dom.get("measure_name"), dom.get("measure")
    nvac = sum(vac)
    if dom.get("n_vacuum_axes") is not None and nvac != dom["n_vacuum_axes"]:
        msgs.append(f"the structure is {_dim(nvac)} (vacuum along {nvac} axis/axes) but the model was trained on "
                    f"{dom.get('dimensionality')} structures.")
    elif name and dname == name and rng:
        lo, hi = rng["min"] * (1 - MEASURE_TOL), rng["max"] * (1 + MEASURE_TOL)
        if not lo <= val <= hi:
            msgs.append(f"{name} is {val:.3f}, outside the trained range {rng['min']:.3f} .. {rng['max']:.3f} "
                        f"(+/- {100 * MEASURE_TOL:g} %): the model is extrapolating in cell size.")
    d = shortest_distance(atoms)
    dd = dom.get("d_min_ang")
    if d is not None and dd and d < D_MIN_FACTOR * dd["min"]:
        msgs.append(f"shortest interatomic distance {d:.3f} Ang is below the shortest trained one "
                    f"({dd['min']:.3f} Ang): bonds are closer than anything in the training data.")
    return msgs


def _dim(n_vacuum):
    return {0: "3D", 1: "2D", 2: "1D", 3: "0D"}.get(n_vacuum, "?")


_COLLECTOR = None


def collect_warnings(collector):
    """Send the card's run-time warnings (structure outside the trained range, forces beyond the
    training forces) to `collector`, a list, instead of stderr -- for a tool that lists every
    warning in one final section. None restores printing."""
    global _COLLECTOR
    _COLLECTOR = collector


def _warn_once(message, tag="model card"):
    if message in _WARNED:
        return
    _WARNED.add(message)
    if _COLLECTOR is not None:
        _COLLECTOR.append(f"{tag}: {message}")
        return
    print(f"[WARNING] {tag}: {message}", file=sys.stderr)


def summary_lines(card, card_path=None):
    """A short human summary printed by every tool that loads the model."""
    m, dom, rec = card.get("model") or {}, card.get("domain") or {}, card.get("recipe") or {}
    lt, met = card.get("level_of_theory") or {}, card.get("metrics") or {}
    lines = [f"Model card        : {card_path or '(loaded)'}"
             f" -- {m.get('name', '?')}, seed {m.get('seed', '?')}, round {(card.get('committee') or {}).get('round', '?')}"
             f" ({(card.get('committee') or {}).get('size', '?')}-model committee)"]
    lines.append(f"Trained on        : {', '.join(card.get('elements') or ['?'])}, "
                 f"{(card.get('training_data') or {}).get('n_configs', '?')} SIESTA configs, "
                 f"{dom.get('dimensionality', '?')}"
                 + (f", {dom['measure_name']} {dom['measure']['min']:.2f}..{dom['measure']['max']:.2f}"
                    if dom.get("measure") and dom.get("measure_name") else ""))
    s = lt.get("settings") or {}
    theory = ", ".join(f"{k} {s[k]}" for k in ("XC.functional", "XC.authors", "PAO.BasisSize", "Mesh.Cutoff") if k in s)
    lines.append(f"Level of theory   : SIESTA {lt.get('siesta_version') or ''} {theory}".rstrip()
                 + (", DFT-D3 included in the labels" if (lt.get("dftd3") or {}).get("on") else ""))
    oos = met.get("out_of_sample_committee") or met.get("valid") or {}
    val = card.get("validation") or {}
    if oos.get("F_rmse_meV") is not None:
        lines.append(f"Training accuracy : out-of-sample force RMSE {oos['F_rmse_meV']:.1f} meV/Ang"
                     + (f"; test set {val['overall']['F_rmse_meV']:.1f} meV/Ang ({val.get('verdict') or '-'})"
                        if (val.get("overall") or {}).get("F_rmse_meV") is not None else
                        " (not validated on a test set yet: stb-mlffValidate)"))
    lines.append(f"Energy scale      : SIESTA absolute total energies (eV) -- not comparable with the "
                 "foundation model's absolute energies; differences are what is physical")
    return lines


def register(model_path, card_path, echo=True):
    """Loads and verifies `card_path` against `model_path`, remembers the pair for
    get_calculator() and (echo) prints its summary. Returns the card."""
    if not os.path.isfile(model_path):
        raise ModelCardError(f"model file '{model_path}' not found.")
    card = load_card(card_path)
    verify_model(card, model_path)
    _REGISTRY[os.path.realpath(model_path)] = card
    _CARD_PATHS[os.path.realpath(model_path)] = card_path
    if echo:
        for line in summary_lines(card, card_path):
            print(line)
    return card


def registered_card(model_path):
    return _REGISTRY.get(os.path.realpath(model_path))


_ANNOUNCED = set()


def announce(model_path):
    """Prints the card's summary the first time its model is loaded in this process."""
    key = os.path.realpath(model_path)
    if key in _ANNOUNCED or key not in _REGISTRY:
        return
    _ANNOUNCED.add(key)
    for line in summary_lines(_REGISTRY[key], _CARD_PATHS.get(key)):
        print(line)


def mark_announced(model_path):
    """Skip announce()'s short summary for this model -- for a tool that prints the card in full itself."""
    _ANNOUNCED.add(os.path.realpath(model_path))


def card_flag(model_flag):
    """'--custom-model' -> '--custom-model-card' (also '--ml-custom-model', '--mace-custom-model')."""
    return model_flag + "-card"


def add_card_argument(parser, model_flag="--custom-model"):
    """Adds the card option that must accompany `model_flag`."""
    parser.add_argument(card_flag(model_flag), default=None, metavar="JSON",
                        help=f"Model card of the {model_flag} file: the {CARD_FILE} that stb-mlffAnalysis "
                             "wrote next to it (training data, level of theory, accuracy, energy scale). "
                             f"Required with {model_flag} -- the card is checked against the model file and "
                             "used to warn when a structure is outside what the model was trained on.")


def setup_from_args(model, card, model_flag="--custom-model", error=None):
    """The one call a tool makes after parsing: validates that a custom model comes with its card,
    loads/verifies/registers the pair and prints the summary. `model`/`card` are the parsed values;
    `error` is the tool's own failure function (parser.error by default behaviour: raises SystemExit).
    Returns the card (None when no custom model was given)."""
    def _fail(msg):
        if error:
            error(msg)
        raise SystemExit(f"error: {msg}")
    flag = card_flag(model_flag)
    if not model:
        if card:
            _fail(f"{flag} was given without {model_flag}.")
        return None
    if not card:
        _fail(f"{model_flag} needs {flag} <JSON>: a fine-tuned model is always used together with the model "
              f"card ({CARD_FILE}) that stb-mlffAnalysis wrote next to it.")
    try:
        return register(model, card, echo=False)
    except ModelCardError as e:
        _fail(f"{model_flag}/{flag}: {e}")


def checked_mace_calculator(model, card, device, dtype):
    """A MACECalculator that, on each new composition/cell it meets, compares the structure with the
    card's trained range (check_structure) and warns once per message, and warns once if the forces it
    predicts are far above the largest training force."""
    from mace.calculators import MACECalculator

    class _Checked(MACECalculator):
        _seen = None

        def calculate(self, atoms=None, properties=None, system_changes=None):
            if system_changes is None:
                from ase.calculators.calculator import all_changes
                system_changes = all_changes
            known = set(card.get("elements") or [])
            absent = sorted({s_ for s_ in atoms.get_chemical_symbols() if s_ not in known}, key=_z) if known else []
            if absent:
                raise ValueError(
                    f"this model was trained only on {', '.join(sorted(known, key=_z))} and cannot evaluate "
                    f"{', '.join(absent)} (see its model card): train a model with those elements in the "
                    "SIESTA data (4.17) or use the MACE-MP-0 foundation model.")
            super().calculate(atoms, properties, system_changes)
            try:
                key = (tuple(atoms.get_chemical_symbols()), tuple(np.round(np.asarray(atoms.cell).ravel(), 2)))
                if self._seen != key:
                    self._seen = key
                    for msg in check_structure(card, atoms):
                        _warn_once(msg)
                fm = ((card.get("domain") or {}).get("force_max_eV_per_ang") or {}).get("max")
                f = self.results.get("forces")
                if fm and f is not None:
                    top = float(np.abs(f).max())
                    if top > FORCE_FACTOR * fm and top > 1.0:
                        _warn_once(f"predicted forces reach {top:.1f} eV/Ang, well above the largest training "
                                   f"force ({fm:.1f} eV/Ang): the geometry is outside what the model was fitted "
                                   "on -- treat the result as unreliable.")
            except Exception:
                pass

    return _Checked(model_paths=[model], device=device, default_dtype=dtype)
