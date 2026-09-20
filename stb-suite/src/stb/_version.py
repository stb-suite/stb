"""The suite's version: one number for every tool, `1.9.<number of commits>`.

`1.9` is the development series (the suite is still under development); the
last part is the number of commits, so it grows with every commit and can be
recomputed from any clone. Every tool reads it through `stb.__version__` (its
`--version` and its banner); no tool carries a version of its own.

Where the number comes from, in order:

1. `STB_VERSION`, if set: an explicit override, for builds that have no git
   history (for example a packaging recipe).
2. The number of commits, when this file sits in a git checkout of the suite:
   an editable/development install therefore always shows the current number,
   with no reinstall.
3. The installed package's metadata: a regular install has no git checkout, so
   it reports the number that was computed when it was built.
4. `1.9.0`, if none of these is available (a source archive without git
   history, built without `STB_VERSION`).

At build time `pyproject.toml` reads `__version__` (case 2 or 1), which is what
puts the number in the package metadata used by case 3.
"""

import os
import subprocess
from importlib import metadata
from pathlib import Path

SERIES = "1.9"

_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[3] if len(_HERE.parents) > 3 else None


def _commit_count():
    """Commits reachable from HEAD if this file is in a checkout of the suite's
    git repository (`<repo>/stb-suite/src/stb/_version.py`), else None."""
    if _REPO is None or _HERE != _REPO / "stb-suite" / "src" / "stb" / "_version.py" \
            or not (_REPO / ".git").exists():
        return None
    try:
        run = subprocess.run(["git", "rev-list", "--count", "HEAD"], cwd=_REPO, capture_output=True,
                             text=True, timeout=10, check=True)
        return int(run.stdout.strip())
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


def get_version():
    override = os.environ.get("STB_VERSION")
    if override:
        return override
    count = _commit_count()
    if count is not None:
        return f"{SERIES}.{count}"
    try:
        return metadata.version("stb_suite")
    except metadata.PackageNotFoundError:
        return f"{SERIES}.0"


__version__ = get_version()
