#!/usr/bin/env python3
"""Writes docs/reference_help.json: the `--help` text of every stb-* command.

The documentation site's Reference pages show these texts. They are captured
here, once, and committed -- instead of at site-build time -- because printing
`--help` means importing each tool, which needs the whole scientific stack and,
for the ML tools, PyTorch and MACE (they refuse to import without them). The CI
that builds the site installs none of that.

Run it whenever a command's options or help text change:

    pip install -e "stb-suite[ml]"        # the package and the ml extra
    python docs/dump_help.py              # rewrites docs/reference_help.json
    python docs/dump_help.py --check      # only compares; exit 1 if out of date

The output is made deterministic (fixed terminal width and hash seed, no colour,
no warnings), and the banner some ML libraries print to stdout before argparse's
own "usage:" line is dropped. The text is otherwise exactly what the tool prints:
many tools emit a whole paragraph as one long line, which the Reference pages wrap
with CSS (docs/assets/reference.css) rather than by rewriting it here.
"""

import json
import os
import re
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "stb-suite" / "pyproject.toml"
SNAPSHOT = ROOT / "docs" / "reference_help.json"

# The interactive menu has no argparse `--help`.
SKIP = {"stb-suite"}
WIDTH = "88"
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def capture(name, target):
    """Returns (help_text, None) for one command, or (None, reason)."""
    module, function = target.split(":")
    code = (f"import sys; sys.argv = [{name!r}, '--help']; "
            f"from {module} import {function} as main; main()")
    # PYTHONHASHSEED: some tools build argparse `choices` from a set, whose order
    # would otherwise change on every run (stb-translate's format lists do).
    env = {**os.environ, "COLUMNS": WIDTH, "PYTHONWARNINGS": "ignore", "PYTHONHASHSEED": "0",
           "NO_COLOR": "1", "PYTHONIOENCODING": "utf-8"}
    try:
        run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             env=env, stdin=subprocess.DEVNULL, timeout=300)
    except subprocess.TimeoutExpired:
        return None, "timed out"
    text = _ANSI.sub("", run.stdout)
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if line.startswith("usage:")), None)
    if run.returncode != 0 or start is None:
        tail = (run.stderr.strip().splitlines() or ["no output"])[-1]
        return None, f"exit {run.returncode}: {tail}"
    return "\n".join(line.rstrip() for line in lines[start:]).rstrip() + "\n", None


def collect():
    scripts = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))["project"]["scripts"]
    targets = {name: target for name, target in scripts.items() if name not in SKIP}
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = dict(zip(targets, pool.map(lambda kv: capture(*kv), targets.items())))
    data, failures = {}, []
    for name, (text, reason) in sorted(results.items()):
        if text is None:
            failures.append(f"  {name}: {reason}")
        else:
            data[name] = {"module": targets[name].split(":")[0], "help": text}
    return data, failures


def main():
    check = "--check" in sys.argv[1:]
    data, failures = collect()
    if failures:
        print("Could not capture --help for:\n" + "\n".join(failures), file=sys.stderr)
        return 1
    rendered = json.dumps(data, indent=1, sort_keys=True, ensure_ascii=False) + "\n"
    if check:
        current = SNAPSHOT.read_text(encoding="utf-8") if SNAPSHOT.is_file() else ""
        if current == rendered:
            print(f"{SNAPSHOT.name} is up to date ({len(data)} commands)")
            return 0
        old = json.loads(current) if current else {}
        changed = sorted(n for n in set(old) | set(data) if old.get(n) != data.get(n))
        print(f"{SNAPSHOT.name} is out of date; differing commands: {', '.join(changed)}\n"
              "Run: python docs/dump_help.py", file=sys.stderr)
        return 1
    SNAPSHOT.write_text(rendered, encoding="utf-8")
    print(f"Wrote {SNAPSHOT.relative_to(ROOT)}: {len(data)} commands, {len(rendered) // 1024} KiB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
