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
"""

import ast
import json
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MENU_SOURCE = ROOT / "stb-suite" / "src" / "stb" / "stb_suite.py"
PYPROJECT = ROOT / "stb-suite" / "pyproject.toml"
HELP_SNAPSHOT = ROOT / "docs" / "reference_help.json"

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
    """The menu as a list of `Entry`, in menu order."""
    tree = ast.parse(MENU_SOURCE.read_text(encoding="utf-8"))
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
                          title=raw["title"], description=raw.get("description", ""))
            if "stages" in raw:
                for stage, stage_raw in raw["stages"].items():
                    entry.leaves.append(Leaf(f"{entry.code}.{stage}", stage_raw["title"],
                                             commands_of(stage_raw.get("func"))))
            else:
                entry.leaves.append(Leaf(entry.code, raw["title"], commands_of(raw.get("func"))))
            menu.append(entry)
    return sorted(menu, key=lambda e: (e.cat, e.item))


# --- Commands and their --help --------------------------------------------------

def load_scripts():
    """Console commands: name -> "module:function", in pyproject order."""
    data = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    return dict(data["project"]["scripts"])


def load_help():
    """The `--help` snapshot: command -> {"module": ..., "help": ...}."""
    if not HELP_SNAPSHOT.is_file():
        return {}
    return json.loads(HELP_SNAPSHOT.read_text(encoding="utf-8"))


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
