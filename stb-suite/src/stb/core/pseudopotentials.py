"""Bundled pseudopotential banks shipped with the package (see
stb-suite/src/stb/pseudopotentials/), plus resolution of a user-supplied
--pseudo-dir/--pp-path/-p value into an actual directory: either one of the
bundled bank names below, or a literal filesystem path. Shared by every tool
that needs pseudopotentials copied/symlinked into a run folder
(hubbardu.py, cohesive_energy.py, inputfile.py, phonons_create.py).
"""

import os
import importlib.resources

# Keys double as both the folder name under stb/pseudopotentials/ and the
# value users type/select to pick a bank -- keep them in sync with the
# directories on disk.
BANKS = {
    "dojo": "PseudoDojo v0.5 (PBE, standard, PSML format)",
    "virtual_vault": "SIESTA Pseudopotentials Virtual Vault (PSF format)",
}


def bank_path(name):
    """Absolute path to a bundled pseudopotential bank folder."""
    return str(importlib.resources.files("stb").joinpath("pseudopotentials", name))


def resolve_pseudo_source(value):
    """`value` is either a bundled bank name (a key in BANKS) or a
    filesystem path ('~' expanded). Returns the resolved, existing
    directory path. Raises ValueError with a clear message if it's neither.
    """
    if value in BANKS:
        return bank_path(value)
    path = os.path.expanduser(value)
    if not os.path.isdir(path):
        raise ValueError(
            f"'{value}' is not a recognized pseudopotential bank "
            f"({', '.join(BANKS)}) nor an existing directory."
        )
    return path
