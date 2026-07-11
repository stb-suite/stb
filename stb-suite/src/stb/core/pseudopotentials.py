"""Bundled pseudopotential banks shipped with the package (see
stb-suite/src/stb/pseudopotentials/), plus resolution of a user-supplied
--pseudo-dir/--pp-path/-p value into an actual directory: either one of the
bundled bank names below, or a literal filesystem path. Shared by every tool
that needs pseudopotentials copied/symlinked into a run folder
(hubbardu.py, cohesive_energy.py, inputfile.py, phonons_create.py).
"""

import os
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
