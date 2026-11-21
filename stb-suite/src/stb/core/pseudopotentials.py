"""Bundled pseudopotential banks shipped with the package (see
stb-suite/src/stb/pseudopotentials/), plus resolution of a user-supplied
--pseudo-dir/--pp-path/-p value into an actual directory: either one of the
bundled bank names below, or a literal filesystem path. Shared by every tool
that needs a pseudopotential copied into a run folder (hubbardu.py,
cohesive_energy.py, inputfile.py, phonons_create.py, adsorb.py, ...).

The PseudoDojo bank holds several tables, one folder each under
pseudopotentials/dojo/ (`<relativity>_<xc>_<version>_<accuracy>`, see the
README.md there). A user picks one with a single string, `dojo` followed by
hyphen-separated words in any order -- relativity (sr, fr/soc, sr3plus),
functional (pbe, pbesol, pw/lda) and accuracy (standard, stringent) --, each
missing word taking its default (sr, pbe, standard): `dojo` is
`dojo-sr-pbe-standard` (the table this bank held before it was split),
`dojo-fr-pbesol` the full-relativistic PBEsol standard table. One string keeps
every tool's existing -p option working with no new flag.

`check_pseudo_level` / `warn_pseudo_level` compare the pseudopotentials copied
into a run folder with that folder's own .fdf (functional, spin-orbit). They
only WARN -- a mismatch never stops a tool.

Pseudopotentials are always COPIED (shutil.copy2), never symlinked -- a
symlink into a bundled bank or an external --pseudo-dir points at an
absolute path outside the run folder, which silently breaks the moment that
folder is archived/rsynced/scp'd elsewhere (e.g. to an HPC cluster) without
also bringing the link target along; a real copy travels with the folder.
"""

import os
import re
import shutil
import atexit
import importlib.resources

from stb.core.cli import color_text

PSEUDO_DOJO_CITATION = (
    "van Setten, Giantomassi, Bousquet, Verstraete, Hamann, Gonze, Rignanese, "
    "\"The PseudoDojo: Training and grading a 85 element optimized norm-conserving "
    "pseudopotential table\", Comput. Phys. Commun. 226, 39-54 (2018), "
    "doi:10.1016/j.cpc.2018.01.012")

# Bank families. 'dojo' stands for its default table (dojo-sr-pbe-standard);
# the other tables are reached through DOJO_TABLES below.
BANKS = {
    "dojo": {
        "description": "PseudoDojo (PSML format): scalar/full-relativistic, PBE/PBEsol/LDA, standard/stringent",
        "citation": PSEUDO_DOJO_CITATION,
        "url": "http://www.pseudo-dojo.org/",
    },
    "virtual_vault": {
        "description": "SIESTA Pseudopotentials Virtual Vault (PSF format)",
        "citation": "NNIN/C Pseudopotential Virtual Vault, Cornell NanoScale Science and "
                    "Technology Facility",
        "url": "https://nninc.cnf.cornell.edu/",
    },
}

# (relativity, xc, accuracy) -> (folder under pseudopotentials/dojo/, PseudoDojo table)
DOJO_TABLES = {
    ("sr", "pbe", "standard"): ("sr_pbe_0.4_standard", "nc-sr-04 pbe standard"),
    ("sr", "pbe", "stringent"): ("sr_pbe_0.4_stringent", "nc-sr-04 pbe stringent"),
    ("sr", "pbesol", "standard"): ("sr_pbesol_0.4_standard", "nc-sr-04 pbesol standard"),
    ("sr", "pbesol", "stringent"): ("sr_pbesol_0.4_stringent", "nc-sr-04 pbesol stringent"),
    ("sr", "pw", "standard"): ("sr_pw_0.4_standard", "nc-sr-04 pw standard"),
    ("sr", "pw", "stringent"): ("sr_pw_0.4_stringent", "nc-sr-04 pw stringent"),
    ("fr", "pbe", "standard"): ("fr_pbe_0.4_standard", "nc-fr-04 pbe standard"),
    ("fr", "pbe", "stringent"): ("fr_pbe_0.4_stringent", "nc-fr-04 pbe stringent"),
    ("fr", "pbesol", "standard"): ("fr_pbesol_0.4_standard", "nc-fr-04 pbesol standard"),
    ("fr", "pbesol", "stringent"): ("fr_pbesol_0.4_stringent", "nc-fr-04 pbesol stringent"),
    ("sr3plus", "pbe", "standard"): ("sr3plus_pbe_0.4_standard", "nc-sr-04-3plus pbe standard"),
}
DOJO_DEFAULT = ("sr", "pbe", "standard")

# Every bundled dojo folder holds this file: the PseudoDojo table and archive it came from, and
# per file its md5, valence, relativity, lj / scalar-relativistic sets and ekb = 0 count
# (written by test/pseudo_bank/make_manifest.py; read by bank_manifest / stb-pseudo --list-banks).
BANK_MANIFEST = "manifest.json"

RELATIVITY_WORDS = {"sr": "sr", "scalar": "sr", "fr": "fr", "soc": "fr", "full": "fr",
                    "sr3plus": "sr3plus", "3plus": "sr3plus"}
XC_WORDS = {"pbe": "pbe", "pbesol": "pbesol", "pw": "pw", "lda": "pw", "pw92": "pw"}
ACCURACY_WORDS = {"standard": "standard", "stringent": "stringent"}

RELATIVITY_TEXT = {"sr": "scalar-relativistic", "fr": "full-relativistic (spin-orbit)",
                   "sr3plus": "scalar-relativistic, lanthanides with 4f in the core (Ln3+)"}
XC_TEXT = {"pbe": "PBE", "pbesol": "PBEsol", "pw": "LDA (PW92)"}

# One line for every --help that takes a bank name.
BANK_HELP = "dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent], virtual_vault"


def _dojo_name(key):
    return "dojo-" + "-".join(key)


def short_bank_name(value):
    """Shortest name of a bank: the dojo words that differ from the defaults
    (dojo-sr-pbe-standard -> dojo, dojo-fr-pbesol-standard -> dojo-fr-pbesol)."""
    key = parse_dojo_name(value)
    if key is None:
        return value
    words = [w for w, d in zip(key, DOJO_DEFAULT) if w != d]
    if key[0] == "sr3plus":
        words = ["sr3plus"]
    return "-".join(["dojo"] + words)


def parse_dojo_name(value):
    """(relativity, xc, accuracy) of a `dojo...` bank name, or None when `value`
    is not one. Raises ValueError for a dojo name with an unknown word, a word
    given twice, or a combination PseudoDojo does not publish."""
    if not isinstance(value, str):
        return None
    words = [w for w in re.split(r"[-_:\s]+", value.strip().lower()) if w]
    if not words or words[0] != "dojo":
        return None
    picked = {}
    for word in words[1:]:
        if re.fullmatch(r"v?0\.\d", word):       # a table version, e.g. 0.4 -- implied by the rest
            continue
        for field, table in (("rel", RELATIVITY_WORDS), ("xc", XC_WORDS), ("acc", ACCURACY_WORDS)):
            if word in table:
                if field in picked and picked[field] != table[word]:
                    raise ValueError(f"'{value}': '{word}' conflicts with '{picked[field]}'.")
                picked[field] = table[word]
                break
        else:
            raise ValueError(f"'{value}': unknown word '{word}'. Bank names: {BANK_HELP}.")
    key = (picked.get("rel", DOJO_DEFAULT[0]), picked.get("xc", DOJO_DEFAULT[1]),
           picked.get("acc", DOJO_DEFAULT[2]))
    if key not in DOJO_TABLES:
        if key[0] == "fr" and key[1] == "pw":
            why = "PseudoDojo has no full-relativistic LDA table."
        elif key[0] == "sr3plus":
            why = "the 4f-in-core lanthanide table (sr3plus) exists only for PBE, standard."
        else:
            why = "PseudoDojo has no such table."
        raise ValueError(f"'{value}' ({_dojo_name(key)}): {why} Available: "
                         + ", ".join(_dojo_name(k) for k in DOJO_TABLES))
    return key


def canonical_bank(value):
    """Canonical bank name of `value` (`dojo-sr-pbe-standard`, `virtual_vault`), None
    when `value` is not a bank name at all; ValueError for an impossible dojo name."""
    if value == "virtual_vault":
        return value
    key = parse_dojo_name(value)
    return _dojo_name(key) if key else None


def is_bank(value):
    """True when `value` names a bundled bank (any valid dojo-... spelling included)."""
    try:
        return canonical_bank(value) is not None
    except ValueError:
        return False


def bank_names():
    """Every canonical bank name: the dojo tables, then virtual_vault."""
    return [_dojo_name(k) for k in DOJO_TABLES] + ["virtual_vault"]


def bank_info(name):
    """{'name', 'description', 'citation', 'url', 'table', 'folder'} of a bank name
    (any accepted spelling). Raises ValueError for an unknown or impossible name."""
    canon = canonical_bank(name)
    if canon is None:
        raise ValueError(f"'{name}' is not a recognized pseudopotential bank ({BANK_HELP}).")
    if canon == "virtual_vault":
        return dict(BANKS["virtual_vault"], name=canon, table=None, folder="virtual_vault")
    key = parse_dojo_name(canon)
    folder, table = DOJO_TABLES[key]
    accuracy = key[2]
    description = (f"PseudoDojo {table} (PSML, {RELATIVITY_TEXT[key[0]]}, {XC_TEXT[key[1]]}, {accuracy})")
    return dict(BANKS["dojo"], name=canon, description=description, table=table,
                folder=os.path.join("dojo", folder))


def bank_path(name):
    """Absolute path to a bundled pseudopotential bank folder (any accepted spelling)."""
    info = bank_info(name)
    return str(importlib.resources.files("stb").joinpath("pseudopotentials", *info["folder"].split(os.sep)))


def resolve_pseudo_source(value):
    """`value` is either a bundled bank name (see BANK_HELP) or a
    filesystem path ('~' expanded). Returns the resolved, existing
    directory path. Raises ValueError with a clear message if it's neither,
    or if it is a dojo name for a table PseudoDojo does not publish.

    Prints the bank's origin/citation when `value` is a bundled bank (never
    for a plain filesystem path) -- every tool that accepts
    --pseudo-dir/--pp-path/-p calls this to resolve it, so the notice lives
    here once instead of being copied into each caller.
    """
    canon = canonical_bank(value)          # ValueError for an impossible dojo-... name
    if canon is not None:
        info = bank_info(canon)
        print(color_text(
            f"[INFO] Using bundled pseudopotential bank '{value}'"
            f"{'' if value == canon else f' ({canon})'}: {info['description']}. "
            f"Source: {info['citation']} -- {info['url']}", 'cyan'))
        return bank_path(canon)
    path = os.path.expanduser(value)
    if not os.path.isdir(path):
        raise ValueError(
            f"'{value}' is not a recognized pseudopotential bank "
            f"({BANK_HELP}) nor an existing directory."
        )
    return path


def detect_pseudo_bank(pp_path):
    """If `pp_path` resolves to one of the bundled bank folders, returns its
    canonical bank name (so the matching citation can be added); None for a
    blank path or a user-supplied folder -- there's no reliable citation for
    pseudopotentials this suite didn't curate itself.
    """
    if not pp_path:
        return None
    normalized = os.path.normpath(pp_path)
    base = os.path.basename(normalized)
    parent = os.path.basename(os.path.dirname(normalized))
    grandparent = os.path.basename(os.path.dirname(os.path.dirname(normalized)))
    if parent == "pseudopotentials" and base == "virtual_vault":
        return "virtual_vault"
    if parent == "dojo" and grandparent == "pseudopotentials":
        for key, (folder, _table) in DOJO_TABLES.items():
            if folder == base:
                return _dojo_name(key)
    return None


def bank_manifest(name):
    """The manifest.json of a bundled dojo bank (any accepted spelling) as a dict, or None
    (virtual_vault, or a missing/unreadable file)."""
    import json
    path = os.path.join(bank_path(name), BANK_MANIFEST)
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def verify_bank(name):
    """Problems (list of str) of a bundled dojo bank against its manifest: missing manifest, files
    missing or not listed, md5 mismatches, KB projectors with ekb = 0. [] when everything matches."""
    import hashlib
    manifest = bank_manifest(name)
    if manifest is None:
        return [] if canonical_bank(name) == "virtual_vault" else [f"{BANK_MANIFEST} missing"]
    folder = bank_path(name)
    present = sorted(n for n in os.listdir(folder) if n.endswith(".psml"))
    listed = manifest.get("files", {})
    problems = [f"{n}: not in {BANK_MANIFEST}" for n in present if n not in listed]
    problems += [f"{n}: listed but missing" for n in listed if n not in present]
    for n in present:
        if n not in listed:
            continue
        with open(os.path.join(folder, n), "rb") as f:
            if hashlib.md5(f.read()).hexdigest() != listed[n].get("md5"):
                problems.append(f"{n}: md5 differs from {BANK_MANIFEST}")
        if listed[n].get("ekb_zero"):
            problems.append(f"{n}: {listed[n]['ekb_zero']} KB projector(s) with ekb = 0")
    return problems


def list_bank_elements(name):
    """Element symbols available in bundled bank `name` (any accepted spelling),
    sorted -- derived from the .psf/.psml filenames actually present, so it
    can never drift from get_required_pseudos()'s own lookup convention.
    Raises ValueError if `name` isn't a known bank.
    """
    path = bank_path(name)
    elements = set()
    for fname in os.listdir(path):
        stem, ext = os.path.splitext(fname)
        if ext in (".psf", ".psml"):
            elements.add(stem)
    return sorted(elements)



def ekb_zero_message(path, count):
    """The [ERROR] text for a pseudopotential file with KB projectors of ekb = 0."""
    return (f"{os.path.basename(path)} has {count} KB projector(s) with ekb = 0 (the PseudoDojo nc-sr-05 "
            "regenerated elements: Ba, Bi, I, Pb, Po, Rb, Rn, S, Te, Tl, Xe). With SIESTA 5.4.2 such a file gives "
            "wrong equilibrium structures with any number of MPI ranks (Pb fcc volume -53 %) and aborts in "
            "parallel with 'propor: ERROR: IMAX = 0'. Use the nc-sr-04 file (the bundled -p dojo / "
            "dojo-stringent).")


def refuse_unusable_pseudo(path):
    """Stops the tool ([ERROR], exit 1) when the pseudopotential file `path` cannot be used:
    KB projectors with ekb = 0. Called before a file is copied or listed for copying, so a run
    folder never receives one; a no-op for any other file."""
    if not path or not path.endswith(".psml"):
        return
    info = read_pseudo_header(path)
    if info and info.get("n_ekb_zero"):
        print(color_text(f"[ERROR] {ekb_zero_message(path, info['n_ekb_zero'])} "
                         f"(file: {path})", 'red'))
        raise SystemExit(1)


def scalar_part_missing_message(names):
    """The [ERROR] text for full-relativistic .psml files with no scalar-relativistic part used
    without spin-orbit."""
    return (f"{', '.join(names)} {'has' if len(names) == 1 else 'have'} only lj projectors (no scalar-relativistic "
            "part, like 68 of the 71 PseudoDojo nc-fr-04 PBEsol files): without 'Spin spin-orbit' SIESTA 5.4.2 "
            "stops with 'generation of scalar-relativistic part is not available yet'. Use -p dojo-pbesol for "
            "a calculation without spin-orbit, or ask for spin-orbit.")


def files_without_scalar_part(paths):
    """Base names of the .psml files in `paths` that have lj projectors and no scalar-relativistic
    set: SIESTA can use them only with spin-orbit."""
    out = []
    for path in paths:
        if path and path.endswith(".psml"):
            info = read_pseudo_header(path)
            if info and info.get("lj") and info.get("sr_part") is False:
                out.append(os.path.basename(path))
    return out


def get_required_pseudos(symbols, pseudo_dir):
    """Checks whether a pseudopotential file (.psf preferred, .psml
    fallback) exists in `pseudo_dir` for every symbol in `symbols`.
    Returns (found_pseudos, missing_elements) -- found_pseudos is the list
    of resolved file paths (for copying, e.g. into disp-NNN/ folders),
    missing_elements is any symbol with neither extension present. Moved
    here from phonons_create.py once stb-raman (Stage 1, same "copy the
    exact pseudos a structure needs into every generated folder" need)
    became a second consumer, same extract-on-second-use policy as
    structure_io.py/core/symmetry.py.
    """
    # sorted(), not a plain set iteration: a bare `set` has no guaranteed
    # order (varies with Python's per-process string-hash randomization),
    # which made the "Found all required"/"Missing pseudopotentials" report
    # lines every caller prints from these lists non-deterministic run to
    # run -- undesirable for a persisted report.
    unique_elements = sorted(set(symbols))
    found_pseudos = []
    missing_elements = []

    for element in unique_elements:
        psf_path = os.path.join(pseudo_dir, f"{element}.psf")
        psml_path = os.path.join(pseudo_dir, f"{element}.psml")

        if os.path.exists(psf_path):
            found_pseudos.append(psf_path)
        elif os.path.exists(psml_path):
            refuse_unusable_pseudo(psml_path)
            found_pseudos.append(psml_path)
        else:
            missing_elements.append(element)

    return found_pseudos, missing_elements


def copy_pseudo(pp_path, symbol, target_dir, dest_label=None):
    """Copies (shutil.copy2, not a symlink -- see the module docstring) the
    pseudopotential (.psml preferred, .psf fallback) for `symbol`, found
    under `pp_path` (already resolved via resolve_pseudo_source), into
    `target_dir` -- no-op if `pp_path` is falsy. A no-op too if the
    destination file already exists (e.g. a rerun over the same folder),
    same idempotency this had back when it symlinked.

    `dest_label` (default: `symbol`) is the species LABEL the destination
    filename must match -- SIESTA resolves a species' pseudopotential file by
    its declared ChemicalSpeciesLabel label, not its Z. Used by
    cohesive_energy.py for BSSE ghost species (e.g. real pseudopotential
    'C.psf' copied as 'C_ghost.psf' so the ghost label 'C_ghost' -- Z=-6, no
    valence charge, same basis as real carbon -- resolves to the same real
    pseudopotential file); every other caller leaves it at the default.
    Moved here once adsorb.py became a second consumer of the exact same
    symbol -> per-folder-file logic cohesive_energy.py already had (named
    link_pseudo at the time, when it still symlinked).
    """
    if not pp_path:
        return
    dest_label = dest_label or symbol
    pp_path_abs = os.path.abspath(pp_path)
    for ext in ("psml", "psf"):
        src = os.path.join(pp_path_abs, f"{symbol}.{ext}")
        if os.path.exists(src):
            refuse_unusable_pseudo(src)
            dst = os.path.join(target_dir, f"{dest_label}.{ext}")
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
            note_pseudo_folder(target_dir)
            return
    print(color_text(f"[WARNING] Pseudopotential '{symbol}.psml' or '{symbol}.psf' not found in {pp_path}", 'yellow'))
    return


def pseudo_source_name(pp_path, label, species_meta):
    """The name (no extension) under `pp_path` whose .psml/.psf species `label` uses: the
    label itself without a '_ghost' suffix when such a file exists (an stb-adsorb site
    folder or a bsse/ folder holds 'Ag_slab.psml', 'O_ads.psml'), else its real element
    ('Ag' -- a bank or a folder of bare-element files). Returns the element when neither
    exists, so copy_pseudo's own warning names the file a user would add."""
    base = label[:-len("_ghost")] if label.endswith("_ghost") else label
    from stb.core import structure_io
    element = structure_io.real_element(label, species_meta)
    if pp_path:
        for name in dict.fromkeys((base, element)):
            if any(os.path.exists(os.path.join(pp_path, f"{name}.{ext}")) for ext in ("psml", "psf")):
                return name
    return element


_PSML_ZVAL_RE = re.compile(r'z-pseudo\s*=\s*"([0-9.eE+-]+)"')


def pseudo_valence(label, pseudo_dir):
    """(valence electrons, source note) of species `label`: the `z-pseudo` of its
    .psml when it can be found under `pseudo_dir`, else the textbook valence
    (which differs for semicore pseudopotentials), else (None, ...). Shared by
    stb-optical and stb-eosInputs for their reports."""
    from stb.core import siesta_log
    if pseudo_dir:
        path = os.path.join(pseudo_dir, f"{label}.psml")
        if os.path.isfile(path):
            try:
                with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                    m = _PSML_ZVAL_RE.search(f.read(20000))
                if m:
                    return float(m.group(1)), f"{label}.psml"
            except OSError:
                pass
    zval = siesta_log.FALLBACK_VALENCE.get(label)
    if zval is not None:
        return zval, "textbook valence (no pseudopotential read; semicore pseudos differ)"
    return None, "unknown"



# ---------------------------------------------------------------------------
# Level of theory: the pseudopotentials of a run folder vs its .fdf.
# Only warnings: SIESTA runs a PBE calculation with LDA pseudopotentials without
# complaint, and that is the trap. Spin-orbit with scalar files is the one case
# SIESTA itself stops on ('The PSML file does not contain lj projs').
# ---------------------------------------------------------------------------

# libxc ids (exchange, correlation) -> functional
_LIBXC_PAIRS = {
    frozenset({1, 9}): "LDA", frozenset({1, 12}): "LDA", frozenset({1, 8}): "LDA",
    frozenset({101, 130}): "PBE", frozenset({116, 133}): "PBEsol",
    frozenset({102, 130}): "revPBE", frozenset({117, 130}): "RPBE",
    frozenset({118, 130}): "WC", frozenset({106, 131}): "BLYP",
}
# oncvpsp iexc codes (positive: built-in; negative: -XXXCCC libxc pair)
_ONCV_CODES = {1: "LDA", 2: "LDA", 3: "LDA", 4: "PBE"}
# ATOM's PSF header code (2nd word of the first line)
_PSF_XC = {"ca": "LDA", "pz": "LDA", "pw": "LDA", "pb": "PBE", "ps": "PBEsol",
           "rp": "revPBE", "rv": "RPBE", "wc": "WC", "am": "AM05", "bl": "BLYP",
           "wp": "PW91"}
# SIESTA XC.Authors -> functional (XC.Functional LDA / GGA)
_FDF_AUTHORS = {"ca": "LDA", "pz": "LDA", "pw92": "LDA", "pbe": "PBE", "pbesol": "PBEsol",
                "revpbe": "revPBE", "rpbe": "RPBE", "wc": "WC", "am05": "AM05", "blyp": "BLYP",
                "pw91": "PW91", "pbejsjrlo": "PBEJsJrLO", "pbegcgxlo": "PBEGcGxLO",
                "pbegcgxhego": "PBEGcGxHEG"}

_PSEUDO_CACHE = {}

# A KB projector with ekb exactly zero (PseudoDojo nc-sr-05's 11 regenerated elements).
_EKB_ZERO_RE = re.compile(r'\bekb="[-+]?(?:0+(?:\.0*)?|\.0+)(?:[eEdD][-+]?\d+)?"')


def read_pseudo_header(path):
    """Header facts of a .psml or .psf file:
    'xc': functional or None; 'relativity': 'full'|'scalar'|'non-relativistic'|None;
    'lj': True/False/None (spin-orbit projectors present);
    'sr_part': True/False/None (a .psml's scalar-relativistic projector set: without it, SIESTA
    5.4.2 cannot run the file without spin-orbit; a .psf: None);
    'n_ekb_zero': number of KB projectors with ekb = 0 (.psml; 0 for a .psf);
    'z_valence': z-pseudo (float) or None;
    'valence': [(n, l, occupation), ...] of a .psml's <valence-configuration>, else []."""
    try:
        st = os.stat(path)
    except OSError:
        return None
    key = (os.path.basename(path), st.st_size, int(st.st_mtime))
    if key in _PSEUDO_CACHE:
        return _PSEUDO_CACHE[key]
    info = {"xc": None, "relativity": None, "lj": None, "sr_part": None, "n_ekb_zero": 0,
            "z_valence": None, "valence": []}
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()
    except OSError:
        return None
    if path.endswith(".psml"):
        head = text[:30000]
        ids = {int(i) for i in re.findall(r'<functional\b[^>]*\bid="(\d+)"', head)}
        if ids:
            info["xc"] = _LIBXC_PAIRS.get(frozenset(ids), "libxc " + "+".join(str(i) for i in sorted(ids)))
        else:
            m = re.search(r'oncvpsp-xc-code="(-?\d+)"', head)
            if m:
                code = int(m.group(1))
                if code < 0:
                    pair = frozenset({-code // 1000, -code % 1000})
                    info["xc"] = _LIBXC_PAIRS.get(pair, f"libxc {-code // 1000}+{-code % 1000}")
                else:
                    info["xc"] = _ONCV_CODES.get(code)
        m = re.search(r'\brelativity="([a-z-]+)"', head)
        if m:
            info["relativity"] = {"dirac": "full", "scalar": "scalar", "no": "non-relativistic"}.get(
                m.group(1), m.group(1))
        info["lj"] = 'set="lj"' in text
        info["sr_part"] = 'set="scalar_relativistic"' in text
        info["n_ekb_zero"] = len(_EKB_ZERO_RE.findall(text))
        m = _PSML_ZVAL_RE.search(head)
        if m:
            info["z_valence"] = float(m.group(1))
        m = re.search(r"<valence-configuration\b.*?</valence-configuration>", head, re.S)
        if m:
            info["valence"] = [(int(n), l, float(occ)) for n, l, occ in re.findall(
                r'<shell\s+n="(\d+)"\s+l="([spdf])"\s+occupation="([0-9.eE+-]+)"', m.group(0))]
    else:
        first = text.splitlines()[0].split() if text.strip() else []
        if len(first) >= 3:
            info["xc"] = _PSF_XC.get(first[1].lower(), first[1])
            info["relativity"] = {"rel": "full", "nrl": "scalar or non-relativistic",
                                  "isp": "scalar or non-relativistic"}.get(first[2].lower(), first[2])
            info["lj"] = first[2].lower() == "rel"
    _PSEUDO_CACHE[key] = info
    return info


def _expand_fdf(path, depth=0, seen=None):
    """The text of an .fdf with every `%include` replaced by the included file, in place
    (SIESTA's reading order); missing includes are skipped."""
    seen = seen if seen is not None else set()
    real = os.path.realpath(path)
    if depth > 8 or real in seen:
        return ""
    seen.add(real)
    try:
        lines = open(path, errors="replace").read().splitlines()
    except OSError:
        return ""
    out = []
    for raw in lines:
        m = re.match(r"\s*%include\s+(\S+)", raw, re.I)
        if m:
            inc = os.path.join(os.path.dirname(path), m.group(1))
            out.append(_expand_fdf(inc, depth + 1, seen))
        else:
            out.append(raw)
    return "\n".join(out)


def _root_fdfs(folder):
    """The .fdf files of `folder` that no other .fdf there includes (calc.fdf first)."""
    fdfs = sorted(f for f in os.listdir(folder) if f.endswith(".fdf"))
    included = set()
    for f in fdfs:
        try:
            for m in re.finditer(r"^\s*%include\s+(\S+)", open(os.path.join(folder, f), errors="replace").read(),
                                 re.I | re.M):
                included.add(os.path.normpath(m.group(1)))
        except OSError:
            pass
    roots = [f for f in fdfs if f not in included]
    roots.sort(key=lambda f: (f != "calc.fdf", f))
    return roots


def _fdf_level(text):
    """(functional or None when not comparable, how the input sets it, spin-orbit on?, species labels)."""
    from stb.core.fdf_directives import parse_directives, normalize_label
    d = parse_directives(text)
    species = []
    m = re.search(r"^[ \t]*%block\s+chemical[._-]?species[._-]?label\s*\n(.*?)^[ \t]*%endblock", text, re.I | re.S | re.M)
    if m:
        for line in m.group(1).splitlines():
            parts = line.split("#", 1)[0].split()
            if len(parts) >= 3:
                species.append(parts[2])
    structure_only = {"numberofatoms", "numberofspecies", "latticeconstant", "atomiccoordinatesformat",
                      "atomiccoordinatesoffset", "systemname", "systemlabel"}
    if not (set(d) - structure_only):
        return None, None, False, list(dict.fromkeys(species))   # a structure file, not a calculation input
    func = (d.get("xcfunctional") or "").split()
    auth = (d.get("xcauthors") or "").split()
    func = func[0] if func else ""
    auth = auth[0] if auth else ""
    if re.search(r"^[ \t]*%block\s+xc[._-]?mix", text, re.I | re.M) or func.upper() in ("VDW", "LIBXC"):
        xc, how = None, None
    elif not func and not auth:
        xc, how = "LDA", "no XC.Functional/XC.Authors: SIESTA's default, LDA"
    else:
        xc = _FDF_AUTHORS.get(normalize_label(auth)) if auth else ("LDA" if func.upper() == "LDA" else None)
        how = f"XC.Functional {func or '(default)'} / XC.Authors {auth or '(default)'}"
    spin = normalize_label((d.get("spin") or "").split()[0] if d.get("spin") else "")
    soc = spin.startswith("spinorbit") or spin in ("so", "soc", "sooffsite", "soonsite")
    old = (d.get("spinorbit") or "").split()
    if old and old[0].lower().strip(".") in ("true", "t", "yes"):
        soc = True
    return xc, how, soc, list(dict.fromkeys(species))


# Elements whose perturbative polarization orbital SIESTA fails to generate with the bundled
# pseudopotentials ("Generation of polarization orbital failed ... Rc no bigger than 1.32 Bohr":
# the Li 2s shell, rc ~11 Bohr, with the PseudoDojo Li; seen in LiF and in Li on graphene).
# Others (Na, K, ...) are not verified.
NONPERTURBATIVE_POLARIZATION_Z = {3: "Li"}


def polarization_scheme_block(labels):
    """%block PAO.PolarizationScheme giving `labels` the non-perturbative scheme."""
    return ("%block PAO.PolarizationScheme\n" + "".join(f"  {lb} non-perturbative\n" for lb in labels)
            + "%endblock PAO.PolarizationScheme\n")


def missing_nonperturbative(text):
    """Species labels (Z in NONPERTURBATIVE_POLARIZATION_Z, ghosts included) of an fdf text that
    get a polarization orbital (PAO.BasisSize with P, DZP by default; no explicit PAO.Basis entry)
    by the perturbative scheme, which fails for them. [] for a structure-only file."""
    from stb.core.fdf_directives import parse_directives, normalize_label
    d = parse_directives(text)
    m = re.search(r"^[ \t]*%block\s+chemical[._-]?species[._-]?label\s*\n(.*?)^[ \t]*%endblock", text, re.I | re.S | re.M)
    if not m:
        return []
    labels = []
    for line in m.group(1).splitlines():
        parts = line.split("#", 1)[0].split()
        if len(parts) >= 3 and parts[1].lstrip("+-").isdigit() and abs(int(parts[1])) in NONPERTURBATIVE_POLARIZATION_Z:
            labels.append(parts[2])
    if not labels:
        return []
    size = (d.get("paobasissize") or "DZP").split()[0].upper()
    if "P" not in size:
        return []
    flag = (d.get("paopolarizationnonperturbative") or "").split()
    if flag and flag[0].lower().strip(".") in ("true", "t", "yes"):
        return []
    blocks = {}
    for blk in re.finditer(r"^[ \t]*%block\s+(\S+)\s*\n(.*?)^[ \t]*%endblock", text, re.I | re.S | re.M):
        blocks.setdefault(normalize_label(blk.group(1)), blk.group(2))   # SIESTA keeps the FIRST block
    scheme = blocks.get("paopolarizationscheme", "")
    explicit = blocks.get("paobasis", "")
    missing = []
    for lb in labels:
        if re.search(rf"^\s*{re.escape(lb)}\s+non[-_]?perturbative", scheme, re.I | re.M):
            continue
        if re.search(rf"^\s*{re.escape(lb)}\s", explicit, re.M):
            continue
        missing.append(lb)
    return missing


def add_nonperturbative_block(calc_text, species_z):
    """(calc_text with a PAO.PolarizationScheme block for every label of `species_z`
    ({label: Z}, ghosts with negative Z) that needs it, labels added). Inserted right after the
    PAO.EnergyShift line (end of the text when there is none); unchanged when nothing is missing."""
    species_block = ("\n%block ChemicalSpeciesLabel\n"
                     + "".join(f" {i} {z} {lb}\n" for i, (lb, z) in enumerate(species_z.items(), 1))
                     + "%endblock ChemicalSpeciesLabel\n")
    labels = missing_nonperturbative(calc_text + species_block)
    if not labels:
        return calc_text, []
    block = ("# Li: SIESTA's perturbative polarization orbital fails to generate (2s, rc ~11 Bohr);\n"
             "# the non-perturbative scheme builds it from a real atomic calculation.\n"
             + polarization_scheme_block(labels))
    lines = calc_text.splitlines(keepends=True)
    at = next((k + 1 for k, ln in enumerate(lines) if ln.strip().lower().startswith("pao.energyshift")),
              len(lines))
    if at == len(lines) and lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    lines[at:at] = [block]
    return "".join(lines), labels


def _suggest_table(xc, soc):
    """' (-p dojo-...)' for a functional PseudoDojo has a table of, else ''."""
    xc_word = {"PBE": "", "PBEsol": "-pbesol", "LDA": "-pw"}.get(xc)
    if xc_word is None or (soc and xc == "LDA"):
        return ""
    return f" (-p dojo{'-fr' if soc else ''}{xc_word})"


def check_pseudo_level(folder):
    """Problems of one run folder, as (level, message) pairs, nothing printed: its pseudopotentials
    against each other and against its .fdf -- functional, spin-orbit with files that have no lj
    projectors, polarization orbitals SIESTA cannot build. level is 'error' when SIESTA cannot run
    or gives wrong results (a file with ekb = 0 projectors), 'warning' otherwise. Empty list when
    the folder has no .fdf or no pseudopotential."""
    try:
        names = os.listdir(folder)
    except OSError:
        return []
    pseudo_files = {os.path.splitext(n)[0]: n for n in sorted(names, key=lambda n: n.endswith(".psf"))
                    if n.endswith((".psml", ".psf"))}
    if not pseudo_files:
        return []
    warnings = []
    roots = _root_fdfs(folder)
    if roots:
        for fname in sorted(pseudo_files.values()):
            info = read_pseudo_header(os.path.join(folder, fname))
            if info and info.get("n_ekb_zero"):
                warnings.append(("error", ekb_zero_message(fname, info["n_ekb_zero"])))
    for root in roots:
        text = _expand_fdf(os.path.join(folder, root))
        xc, how, soc, species = _fdf_level(text)
        nonpert = missing_nonperturbative(text) if how is not None or xc is not None else []
        if nonpert:
            warnings.append(("warning",
                f"{root} gives {', '.join(nonpert)} a polarization orbital by SIESTA's default perturbative "
                "scheme, which fails for Li with the PseudoDojo pseudopotential ('Generation of polarization "
                "orbital failed ... Rc no bigger than 1.32 Bohr'). Add "
                f"'%block PAO.PolarizationScheme' with '{nonpert[0]} non-perturbative' (one line per Li label, "
                "ghosts too) or 'PAO.Polarization.NonPerturbative true'."))
        labels = [s for s in species if s in pseudo_files] or sorted(pseudo_files)
        headers = {}
        for label in labels:
            info = read_pseudo_header(os.path.join(folder, pseudo_files[label]))
            if info:
                headers[pseudo_files[label]] = info
        if not headers:
            continue
        by_xc = {}
        for fname, info in headers.items():
            by_xc.setdefault(info["xc"] or "unknown", []).append(fname)
        known = {k: v for k, v in by_xc.items() if k != "unknown"}
        if len(known) > 1:
            warnings.append(("warning", "the pseudopotentials were generated with different functionals: "
                            + "; ".join(f"{k}: {', '.join(v)}" for k, v in sorted(known.items()))))
        if xc is not None:
            wrong = {k: v for k, v in known.items() if k != xc}
            if wrong:
                warnings.append(("warning",
                    f"{root} uses {xc} ({how}), but "
                    + "; ".join(f"{', '.join(v)} {'was' if len(v) == 1 else 'were'} generated with {k}"
                                for k, v in sorted(wrong.items()))
                    + ". SIESTA runs anyway, with an inconsistent level of theory: use pseudopotentials of the "
                      f"calculation's functional{_suggest_table(xc, soc)} or change XC.Functional/XC.Authors."))
        no_scalar = [f for f, info in headers.items() if info.get("lj") and info.get("sr_part") is False]
        if no_scalar and not soc and (how is not None or xc is not None):
            warnings.append(("error", f"{root} does not ask for spin-orbit, but "
                             + scalar_part_missing_message(no_scalar)))
        rel = {}
        for fname, info in headers.items():
            rel.setdefault("full" if info["lj"] else "scalar", []).append(fname)
        if soc and "scalar" in rel:
            warnings.append(("warning",
                f"{root} asks for spin-orbit, but {', '.join(rel['scalar'])} "
                f"{'has' if len(rel['scalar']) == 1 else 'have'} no lj (full-relativistic) projectors: SIESTA stops "
                "with 'The PSML file does not contain lj projs'. Use a full-relativistic table (-p dojo-fr, "
                "dojo-fr-pbesol)."))
        elif len(rel) > 1:
            warnings.append(("warning",
                f"full-relativistic ({', '.join(rel['full'])}) and scalar-relativistic "
                f"({', '.join(rel['scalar'])}) pseudopotentials are mixed: do not compare absolute energies with "
                "runs made with the other kind."))
    return warnings


_NOTED_FOLDERS = []
_CHECKED_FOLDERS = set()


def note_pseudo_folder(folder):
    """Records a run folder that received pseudopotentials; when the process ends, every
    recorded folder not checked already is compared with its .fdf (warn_pseudo_level), so
    each tool that writes run folders gets the level-of-theory warning without its own call.
    copy_pseudo records its target itself; a tool that copies the files another way calls this."""
    path = os.path.abspath(folder)
    if not _NOTED_FOLDERS:
        atexit.register(_check_noted_folders)
    if path not in _NOTED_FOLDERS:
        _NOTED_FOLDERS.append(path)


def _check_noted_folders():
    pending = [f for f in _NOTED_FOLDERS if f not in _CHECKED_FOLDERS]
    if pending:
        try:
            warn_pseudo_level(pending, recursive=False)
        except Exception:
            pass


def warn_pseudo_level(paths, recursive=True, max_folders=5000):
    """Runs check_pseudo_level on every folder under `paths` (a folder or a list; walked
    recursively unless `recursive=False`, hidden folders skipped) that holds an .fdf and a pseudopotential, prints
    each distinct problem once, naming where it was found -- errors first as a red [ERROR], then
    warnings as a yellow [WARNING] -- and returns the distinct messages. Never raises: a check that cannot run is skipped silently."""
    if isinstance(paths, str):
        paths = [paths]
    found = {}
    count = 0
    for top in paths:
        if not top or not os.path.isdir(top):
            continue
        for root, dirs, files in os.walk(top):
            dirs[:] = sorted(d for d in dirs if not d.startswith(".")) if recursive else []
            if not any(f.endswith((".psml", ".psf")) for f in files) or not any(f.endswith(".fdf") for f in files):
                continue
            count += 1
            if count > max_folders:
                break
            _CHECKED_FOLDERS.add(os.path.abspath(root))
            try:
                messages = check_pseudo_level(root)
            except Exception:       # a report check must never break a tool
                continue
            for item in messages:
                found.setdefault(item, []).append(root)
    for (level, msg), where in sorted(found.items(), key=lambda kv: kv[0][0] != "error"):
        place = os.path.relpath(where[0])
        more = f" and {len(where) - 1} more folder(s)" if len(where) > 1 else ""
        if level == "error":
            print(color_text(f"[ERROR] Pseudopotentials ({place}{more}): {msg}", 'red'))
        else:
            print(color_text(f"[WARNING] Pseudopotentials vs input ({place}{more}): {msg}", 'yellow'))
    return [msg for _level, msg in found]
