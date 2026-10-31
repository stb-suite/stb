"""The suite's version: one number for every tool.

Every tool reads it through `stb.__version__` (its `--version` and its banner); no
tool carries a version of its own.

Two kinds of branch, set by the constants below:

- a release branch (`main`, `DEV = False`): the version is `RELEASE` exactly, for
  example `2.0.0`; a bug-fix release raises the last number (`2.0.1`);
- the development branch (`dev`, `DEV = True`): `RELEASE` is the NEXT version and
  the version is `RELEASE.dev<N>` (PEP 440), N being the number of commits since
  the tag `BASE_TAG` of the last release, so it grows with every commit and can be
  recomputed from any clone.

Where the number comes from, in order:

1. `STB_VERSION`, if set: an explicit override (for example a packaging recipe).
2. On a release branch, `RELEASE`.
3. On the development branch, in a git checkout of the suite: `RELEASE.dev<N>`.
4. The installed package's metadata (an installed copy has no checkout, so it
   reports the number computed when it was built).
5. `RELEASE.dev0`, if none of these is available.

At build time `pyproject.toml` reads `__version__`, which is what puts the number
in the package metadata used by case 4.
"""

import os
import subprocess
from importlib import metadata
from pathlib import Path

RELEASE = "2.1.0"
DEV = True
BASE_TAG = "v2.0.0"

_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[3] if len(_HERE.parents) > 3 else None


def _commits_since_tag():
    """Commits since BASE_TAG if this file is in a git checkout of the suite
    (`<repo>/stb-suite/src/stb/_version.py`) that has the tag, else None."""
    if _REPO is None or _HERE != _REPO / "stb-suite" / "src" / "stb" / "_version.py" \
            or not (_REPO / ".git").exists():
        return None
    try:
        run = subprocess.run(["git", "rev-list", "--count", f"{BASE_TAG}..HEAD"], cwd=_REPO,
                             capture_output=True, text=True, timeout=10, check=True)
        return int(run.stdout.strip())
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


def get_version():
    override = os.environ.get("STB_VERSION")
    if override:
        return override
    if not DEV:
        return RELEASE
    count = _commits_since_tag()
    if count is not None:
        return f"{RELEASE}.dev{count}"
    try:
        return metadata.version("stb_suite")
    except metadata.PackageNotFoundError:
        return f"{RELEASE}.dev0"


__version__ = get_version()
