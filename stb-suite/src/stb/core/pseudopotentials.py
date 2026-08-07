"""Bundled pseudopotential banks shipped with the package (see
stb-suite/src/stb/pseudopotentials/), plus resolution of a user-supplied
--pseudo-dir/--pp-path/-p value into an actual directory: either one of the
bundled bank names below, or a literal filesystem path. Shared by every tool
that needs a pseudopotential copied into a run folder (hubbardu.py,
cohesive_energy.py, inputfile.py, phonons_create.py, adsorb.py).

Pseudopotentials are always COPIED (shutil.copy2), never symlinked -- a
symlink into a bundled bank or an external --pseudo-dir points at an
absolute path outside the run folder, which silently breaks the moment that
folder is archived/rsynced/scp'd elsewhere (e.g. to an HPC cluster) without
also bringing the link target along; a real copy travels with the folder.
"""

import os
import shutil
import importlib.resources

from stb.core.cli import color_text

# Keys double as both the folder name under stb/pseudopotentials/ and the
# value users type/select to pick a bank -- keep them in sync with the
# directories on disk. Each entry also carries the citation/origin printed
# by resolve_pseudo_source() whenever the bank is actually used, so credit
# goes to the people who generated/curated these pseudopotentials, not just
# to this suite for bundling them.
BANKS = {
    "dojo": {
        "description": "PseudoDojo v0.5 (PBE, standard, PSML format)",
        "citation": "van Setten, Giantomassi, Bousquet, Verstraete, Hamann, Gonze, Rignanese, "
                    "\"The PseudoDojo: Training and grading a 85 element optimized norm-conserving "
                    "pseudopotential table\", Comput. Phys. Commun. 226, 39-54 (2018), "
                    "doi:10.1016/j.cpc.2018.01.012",
        "url": "http://www.pseudo-dojo.org/",
    },
    "virtual_vault": {
        "description": "SIESTA Pseudopotentials Virtual Vault (PSF format)",
        "citation": "NNIN/C Pseudopotential Virtual Vault, Cornell NanoScale Science and "
                    "Technology Facility",
        "url": "https://nninc.cnf.cornell.edu/",
    },
}


def bank_path(name):
    """Absolute path to a bundled pseudopotential bank folder."""
    return str(importlib.resources.files("stb").joinpath("pseudopotentials", name))


def resolve_pseudo_source(value):
    """`value` is either a bundled bank name (a key in BANKS) or a
    filesystem path ('~' expanded). Returns the resolved, existing
    directory path. Raises ValueError with a clear message if it's neither.

    Prints the bank's origin/citation when `value` is one of the bundled
    BANKS (never for a plain filesystem path) -- every tool that accepts
    --pseudo-dir/--pp-path/-p (hubbardu.py, cohesive_energy.py, inputfile.py,
    phonons_create.py) calls this to resolve it, so putting the notice here
    once covers all of them instead of needing it copied into each caller.
    """
    if value in BANKS:
        bank = BANKS[value]
        print(color_text(
            f"[INFO] Using bundled pseudopotential bank '{value}': {bank['description']}. "
            f"Source: {bank['citation']} -- {bank['url']}", 'cyan'))
        return bank_path(value)
    path = os.path.expanduser(value)
    if not os.path.isdir(path):
        raise ValueError(
            f"'{value}' is not a recognized pseudopotential bank "
            f"({', '.join(BANKS)}) nor an existing directory."
        )
    return path


def detect_pseudo_bank(pp_path):
    """If `pp_path` resolves to one of the bundled BANKS folders, returns its
    bank name (so the matching citation can be added); None for a blank path
    or a user-supplied folder -- there's no reliable citation for
    pseudopotentials this suite didn't curate itself. Moved here from
    inputfile.py once stb-pseudo became a second consumer of the exact same
    "which bundled bank is this?" check, same extract-on-second-use policy
    as the rest of core/.
    """
    if not pp_path:
        return None
    normalized = os.path.normpath(pp_path)
    base = os.path.basename(normalized)
    parent = os.path.basename(os.path.dirname(normalized))
    if parent == "pseudopotentials" and base in BANKS:
        return base
    return None


def list_bank_elements(name):
    """Element symbols available in bundled bank `name` (a BANKS key),
    sorted -- derived from the .psf/.psml filenames actually present, so it
    can never drift from get_required_pseudos()'s own lookup convention.
    Raises ValueError via bank_path()'s caller contract if `name` isn't a
    known bank (checked explicitly here since bank_path() itself doesn't
    validate its argument).
    """
    if name not in BANKS:
        raise ValueError(f"'{name}' is not a recognized pseudopotential bank ({', '.join(BANKS)}).")
    path = bank_path(name)
    elements = set()
    for fname in os.listdir(path):
        stem, ext = os.path.splitext(fname)
        if ext in (".psf", ".psml"):
            elements.add(stem)
    return sorted(elements)


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
            dst = os.path.join(target_dir, f"{dest_label}.{ext}")
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
            return
    print(color_text(f"[WARNING] Pseudopotential '{symbol}.psml' or '{symbol}.psf' not found in {pp_path}", 'yellow'))
    return
