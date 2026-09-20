"""Facts about the suite that the documentation site is generated from.

Everything comes from the repository itself, read without importing the
package (importing it would need every scientific dependency, and the ML tools
refuse to import without PyTorch/MACE):

- the interactive menu (`stb_suite.py`'s six tool dictionaries), read with
  `ast`, including which `stb-*` command each menu code runs, found by following
  the `run_tool("stb-...")` calls in the code each entry points at;
- the console commands (`[project.scripts]` in `pyproject.toml`);
- every command's `--help` text, from the snapshot `docs/reference_help.json`
  that `docs/dump_help.py` writes (that one does need the package installed).

These come from this repository and, optionally, from extension repositories
(plugins for the menu, see `stb_suite._load_plugins`) named in `STB_DOCS_EXTRA_ROOTS`
(paths separated by `os.pathsep`). Each has a `docs_plugin.toml` saying where its
menu module, `pyproject.toml`, help snapshot and examples are; the site then covers
the extension's tools too, as if they were part of this repository.
"""

import ast
import json
import os
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MENU_SOURCE = ROOT / "stb-suite" / "src" / "stb" / "stb_suite.py"
PYPROJECT = ROOT / "stb-suite" / "pyproject.toml"
HELP_SNAPSHOT = ROOT / "docs" / "reference_help.json"

EXTRA_ROOTS_ENV = "STB_DOCS_EXTRA_ROOTS"


@dataclass(frozen=True)
class Source:
    """Where a set of menu entries, commands, `--help` texts and guides comes
    from: this repository, or an extension repository."""
    name: str
    menu: Path
    pyproject: Path
    help_snapshot: Path
    examples: Path
    repo_url: str = ""          # "" = this repository (gen_pages.REPO_URL)


PUBLIC = Source("stb", MENU_SOURCE, PYPROJECT, HELP_SNAPSHOT, ROOT / "examples")


def extra_sources():
    """The extension repositories named in STB_DOCS_EXTRA_ROOTS, if any."""
    found = []
    for entry in os.environ.get(EXTRA_ROOTS_ENV, "").split(os.pathsep):
        if not entry.strip():
            continue
        root = Path(entry).expanduser().resolve()
        manifest = root / "docs_plugin.toml"
        if not manifest.is_file():
            raise ValueError(f"{EXTRA_ROOTS_ENV}: {root} has no docs_plugin.toml")
        data = tomllib.loads(manifest.read_text(encoding="utf-8"))
        found.append(Source(name=data.get("name", root.name), menu=root / data["menu"],
                            pyproject=root / data["pyproject"],
                            help_snapshot=root / data["help_snapshot"],
                            examples=root / data.get("examples", "examples"),
                            repo_url=data.get("repo_url", "").rstrip("/")))
    return found


def sources():
    return [PUBLIC] + extra_sources()


# The six menu categories, in menu order: number -> (url slug, display name).
CATEGORIES = {
    1: ("inputs", "Inputs"),
    2: ("structures", "Structures"),
    3: ("analysis", "Analysis"),
    4: ("workflows", "Workflows"),
    5: ("ml-simulations", "ML Simulations"),
    6: ("utils", "Utils"),
}

# The tool dictionaries in stb_suite.py -> their category number.
_MENU_DICTS = {
    "INPUT_TOOLS": 1, "STRUCTURE_TOOLS": 2, "ANALYSIS_TOOLS": 3,
    "WORKFLOW_TOOLS": 4, "MLSIM_TOOLS": 5, "UTILITY_TOOLS": 6,
}

# The interactive menu itself: it has no argparse `--help` to snapshot.
INTERACTIVE_MENU = "stb-suite"


@dataclass
class Leaf:
    """One runnable menu code (`1.3`, or a workflow stage such as `4.1.2`)."""
    code: str
    title: str
    commands: list


@dataclass
class Entry:
    """One top-level menu item (`1.3`, `4.1`): a single tool, or a workflow
    whose stages are its leaves."""
    code: str
    cat: int
    item: int
    title: str
    description: str
    leaves: list = field(default_factory=list)
    source: Source = PUBLIC

    @property
    def commands(self):
        return [c for leaf in self.leaves for c in leaf.commands]

    @property
    def is_workflow(self):
        return len(self.leaves) > 1 or self.leaves[0].code != self.code

    @property
    def display_title(self):
        return clean_title(self.title)


def clean_title(title):
    """Drops the trailing "(stb-command)" the menu appends to some titles (one
    of them spells it with underscores)."""
    return re.sub(r"\s*\(stb[-_][^)]*\)\s*$", "", title)


def slugify(text):
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


# --- The interactive menu ------------------------------------------------------

def load_menu():
    """The menu as a list of `Entry`, in menu order: this repository's plus
    every extension's (a menu number used twice is an error)."""
    menu = []
    for source in sources():
        menu += _read_menu(source)
    seen = {}
    for entry in menu:
        if entry.code in seen:
            raise ValueError(f"menu code {entry.code} is defined by both '{seen[entry.code].source.name}' "
                             f"and '{entry.source.name}'")
        seen[entry.code] = entry
    return sorted(menu, key=lambda e: (e.cat, e.item))


def _read_menu(source):
    """The `Entry` list defined by one source's menu module."""
    tree = ast.parse(source.menu.read_text(encoding="utf-8"))
    functions = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}

    def commands_of(name, seen=()):
        """Every literal `run_tool("stb-...")` reachable from function `name`,
        following calls to other module-level functions."""
        if name in seen or name not in functions:
            return []
        found, callees = [], []
        for node in ast.walk(functions[name]):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id == "run_tool" and node.args \
                        and isinstance(node.args[0], ast.Constant):
                    found.append(node.args[0].value)
                elif node.func.id in functions and node.func.id != name:
                    callees.append(node.func.id)
        for callee in callees:
            found += commands_of(callee, seen + (name,))
        return list(dict.fromkeys(found))

    def literal(node):
        return node.value if isinstance(node, ast.Constant) else None

    def read_entry(node):
        entry = {}
        for key, value in zip(node.keys, node.values):
            name = literal(key)
            if name in ("title", "description"):
                entry[name] = literal(value)
            elif name == "func" and isinstance(value, ast.Name):
                entry["func"] = value.id
            elif name == "stages":
                entry["stages"] = {literal(k): read_entry(v)
                                   for k, v in zip(value.keys, value.values)}
        return entry

    menu = []
    for node in tree.body:
        if not (isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in _MENU_DICTS):
            continue
        cat = _MENU_DICTS[node.targets[0].id]
        for key, value in zip(node.value.keys, node.value.values):
            item, raw = literal(key), read_entry(value)
            entry = Entry(code=f"{cat}.{item}", cat=cat, item=item,
                          title=raw["title"], description=raw.get("description", ""), source=source)
            if "stages" in raw:
                for stage, stage_raw in raw["stages"].items():
                    entry.leaves.append(Leaf(f"{entry.code}.{stage}", stage_raw["title"],
                                             commands_of(stage_raw.get("func"))))
            else:
                entry.leaves.append(Leaf(entry.code, raw["title"], commands_of(raw.get("func"))))
            menu.append(entry)
    return menu


# --- Commands and their --help --------------------------------------------------

def load_command_sources():
    """Console commands with the source that defines each: name -> Source, in
    pyproject order, this repository first (a command defined twice is an error)."""
    owners = {}
    for source in sources():
        data = tomllib.loads(source.pyproject.read_text(encoding="utf-8"))
        for name in data["project"]["scripts"]:
            if name in owners:
                raise ValueError(f"command {name} is defined by both '{owners[name].name}' and '{source.name}'")
            owners[name] = source
    return owners


def load_scripts():
    """Console commands: name -> "module:function", in pyproject order (this
    repository's first, then each extension's)."""
    scripts = {}
    for source in sources():
        data = tomllib.loads(source.pyproject.read_text(encoding="utf-8"))
        scripts.update(data["project"]["scripts"])
    return scripts


def load_help():
    """The `--help` snapshots: command -> {"module": ..., "help": ...}."""
    helps = {}
    for source in sources():
        if source.help_snapshot.is_file():
            helps.update(json.loads(source.help_snapshot.read_text(encoding="utf-8")))
    return helps


_SECTION_HEADER = re.compile(r"^[A-Za-z][^:]*:$")


def help_summary(help_text):
    """The description paragraph argparse prints right after the usage block,
    as one line ("" if the command has none)."""
    lines = help_text.splitlines()
    i = 0
    while i < len(lines) and lines[i].strip():          # the usage block
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    paragraph = []
    while i < len(lines) and lines[i].strip() and not _SECTION_HEADER.match(lines[i].strip()):
        paragraph.append(lines[i].strip())
        i += 1
    return " ".join(paragraph)


def first_sentence(text, limit=170):
    """First sentence of `text`, shortened to `limit` characters."""
    sentence = re.split(r"(?<=[.!?])\s", text, maxsplit=1)[0]
    if len(sentence) > limit:
        sentence = sentence[:limit].rsplit(" ", 1)[0].rstrip(",;:") + " …"
    return sentence
