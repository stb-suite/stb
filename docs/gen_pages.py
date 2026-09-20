"""Generates the documentation site's pages.

Two layers, both built at build time (mkdocs-gen-files) from the repository
itself -- nothing is copied into `docs/`:

Guides (`guides/`), one page per menu item:
- `examples/<code>-<name>/README.md` becomes that item's guide, with links
  between guides rewritten and a footer pointing back at the example folder on
  GitHub and forward to the command reference.
- A menu item with no example gets a short generated page (a "stub") that says
  so and points at its commands, so every tool in the suite has a page.
- `guides/index.md` is `examples/README.md` plus a table of the stubs.

Reference (`reference/`), one page per console command:
- Each page shows the command's `--help` text (from `docs/reference_help.json`,
  see `dump_help.py`), its menu code, and a link to its guide.
- `reference/stb-suite.md` is the map of every menu code.

Sidebars come from generated `SUMMARY.md` files (mkdocs-literate-nav), ordered
numerically by menu code. What the suite contains -- menu items, which command
each runs, console commands -- is read from `stb_suite.py` and `pyproject.toml`
by `stbdocs/catalog.py`, so a new tool shows up with no change here. Anything
inconsistent (a command with no `--help` snapshot, an example whose code is not
in the menu, ...) is logged as a warning, which fails the strict build.

The example READMEs are written for GitHub and rendered by Python-Markdown, so
each is adjusted on the way through: links between example folders become links
between the generated pages, and list structure is normalized
(`stbdocs/mdcompat.py`; the tests in `test/docs/` check that every page keeps the
list and code-block structure GitHub gives it).
"""

import logging
import posixpath
import re
import sys
from pathlib import Path

import mkdocs_gen_files

# mkdocs-gen-files runs this file with runpy.run_path, which does not put its
# folder on sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from stbdocs import catalog  # noqa: E402
from stbdocs.mdcompat import github_markdown_to_python_markdown  # noqa: E402

log = logging.getLogger("mkdocs.plugins.gen_pages")

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"

REPO_URL = "https://github.com/stb-suite/stb"
BRANCH = "main"

FOLDER_RE = re.compile(r"^(\d+)\.(\d+)-.+$")

# Three H1 styles exist across the examples:
#   "3.13 — `stb-sts`: Simulated STS Spectroscopy"
#   "1.3 — K-Grid Generator (`stb-kgrid`)"
#   "4.7 — Workflow: Hubbard U (Linear Response) (`stb-a` / `stb-b`)"
_H1_CMD_FIRST = re.compile(r"^(?P<code>\d+\.\d+) — `[^`]+`: (?P<title>.+)$")
_H1_CMD_LAST = re.compile(r"^(?P<code>\d+\.\d+) — (?P<title>.+?)(?: \(`[^)]*`\))?$")

# A link to another example folder: `](1.1-stb-inputfile/)` or
# `](../4.3-cohesive/)`, optionally with a `#anchor`.
_FOLDER_LINK = re.compile(r"\]\((?:\.\./)?(\d+\.\d+-[^/)#\s]+)/?(#[^)\s]*)?\)")


class Page:
    """One guide page in the sidebar: a real guide (from an example README) or a stub."""

    def __init__(self, entry, dest, label, text=None, edit_source=None, folder=None):
        self.entry, self.dest, self.label = entry, dest, label
        self.text, self.edit_source, self.folder = text, edit_source, folder

    @property
    def is_stub(self):
        return self.folder is None


def guide_title(folder, text):
    """Sidebar title of an example README, from its H1."""
    in_fence, h1 = False, None
    for line in text.splitlines():
        if line.lstrip().startswith(("```", "~~~")):
            in_fence = not in_fence
        elif not in_fence and line.startswith("# "):
            h1 = line[2:].strip()
            break
    if h1 is None:
        log.warning("examples/%s/README.md has no H1; using the folder name as its title", folder)
        return folder
    for pattern in (_H1_CMD_FIRST, _H1_CMD_LAST):
        m = pattern.match(h1)
        if m:
            return re.sub(r"^Workflow: ", "", m.group("title")).replace("`", "")
    log.warning("examples/%s/README.md: unrecognized H1 format %r; using it as-is", folder, h1)
    return h1.replace("`", "")


def collect_pages(menu):
    """Every guide page, real or stub, in menu order."""
    entries = {e.code: e for e in menu}
    pages, seen = [], set()
    for path in sorted(EXAMPLES.iterdir()):
        m = FOLDER_RE.match(path.name)
        if not (path.is_dir() and m):
            continue
        if not (path / "README.md").is_file():
            log.warning("examples/%s has no README.md; skipped", path.name)
            continue
        code = f"{m.group(1)}.{m.group(2)}"
        entry = entries.get(code)
        if entry is None:
            log.warning("examples/%s: no menu item has code %s (see stb_suite.py)", path.name, code)
            continue
        text = (path / "README.md").read_text(encoding="utf-8")
        slug = catalog.CATEGORIES[entry.cat][0]
        pages.append(Page(entry, f"guides/{entry.cat}-{slug}/{path.name}.md",
                          f"{code} {guide_title(path.name, text)}", text,
                          f"examples/{path.name}/README.md", path.name))
        seen.add(code)
    for entry in menu:
        if entry.code not in seen:
            slug = catalog.CATEGORIES[entry.cat][0]
            pages.append(Page(entry, f"guides/{entry.cat}-{slug}/{entry.code}-{catalog.slugify(entry.display_title)}.md",
                              f"{entry.code} {entry.display_title}"))
    return sorted(pages, key=lambda p: (p.entry.cat, p.entry.item))


# --- Small helpers --------------------------------------------------------------

def rel(dest, source_dest):
    """Relative link to page `dest`, written from page `source_dest`."""
    return posixpath.relpath(dest, posixpath.dirname(source_dest))


def ref_dest(command):
    return f"reference/{command}.md"


def md_escape(text):
    """Text that came from a --help description or a menu title, made safe to
    put in a Markdown paragraph or table cell (`*`, `_`, `<`, `|` ... are all
    common in these)."""
    text = text.replace("&", "&amp;").replace("<", "&lt;")
    return re.sub(r"([\\`*_\[\]|])", r"\\\1", text)


def write_page(dest, content, edit_source):
    with mkdocs_gen_files.open(dest, "w") as f:
        f.write(content)
    # edit_uri is rooted at docs/, so reach the repo root with "../"
    mkdocs_gen_files.set_edit_path(dest, f"../{edit_source}")


def write_nav(dest, lines):
    with mkdocs_gen_files.open(dest, "w") as f:
        f.writelines(line + "\n" for line in lines)


def rewrite_folder_links(text, source_page, dest_by_folder):
    """Turns links between example folders into links between the generated
    pages, relative to `source_page`. Fenced code blocks are left untouched."""

    def replace(match):
        folder, anchor = match.group(1), match.group(2) or ""
        dest = dest_by_folder.get(folder)
        if dest is None:
            log.warning("%s links to unknown example folder %r", source_page, folder)
            return match.group(0)
        return f"]({rel(dest, source_page)}{anchor})"

    out, in_fence = [], False
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith(("```", "~~~")):
            in_fence = not in_fence
        out.append(line if in_fence else _FOLDER_LINK.sub(replace, line))
    return "".join(out)


# --- Guides ---------------------------------------------------------------------

def guide_footer(page):
    tree = f"{REPO_URL}/tree/{BRANCH}/examples/{page.folder}"
    parts = [f"Example folder (input fixtures and guided script) on GitHub: [`examples/{page.folder}/`]({tree})"]
    scripts = sorted((EXAMPLES / page.folder).glob("example_*.sh"))
    if scripts:
        blob = f"{REPO_URL}/blob/{BRANCH}/examples/{page.folder}/{scripts[0].name}"
        parts.append(f"guided script [`{scripts[0].name}`]({blob})")
    footer = "\n\n---\n\n*" + " · ".join(parts) + "*\n"
    commands = page.entry.commands
    if commands:
        links = ", ".join(f"[`{c}`]({rel(ref_dest(c), page.dest)})" for c in commands)
        footer += f"\n*Command reference: {links}*\n"
    return footer


def stub_page(page, helps):
    entry = page.entry
    lines = [f"# {entry.code} — {entry.display_title}", ""]
    first = entry.commands[0] if entry.commands else None
    where = f" (the full option list is in the [command reference]({rel(ref_dest(first), page.dest)}))" if first else ""
    lines += [f"> There is no hands-on guide for this tool yet{where}.", ""]
    if entry.description:
        lines += [md_escape(entry.description), ""]
    if not entry.is_workflow:
        summary = catalog.help_summary(helps[first]["help"]) if first in helps else ""
        if summary and summary.rstrip(". ") != entry.description.rstrip(". "):
            lines += ["## What it does", "", md_escape(summary), ""]
        lines += ["## Run it", "",
                  f"- **Interactive menu:** run `stb-suite` and type `{entry.code}` at the main prompt.",
                  f"- **Direct command:** [`{first}`]({rel(ref_dest(first), page.dest)}) "
                  f"(add `--help` to print every option).", ""]
    else:
        lines += ["## Stages", "",
                  "Run the stages in order. In the interactive menu, type a stage's code at the main prompt.",
                  "", "| Stage | Command | What it does |", "|---|---|---|"]
        for leaf in entry.leaves:
            for command in leaf.commands:
                summary = catalog.first_sentence(catalog.help_summary(helps.get(command, {}).get("help", "")))
                lines.append(f"| `{leaf.code}` {md_escape(catalog.clean_title(leaf.title))} | "
                             f"[`{command}`]({rel(ref_dest(command), page.dest)}) | {md_escape(summary)} |")
        lines.append("")
    return "\n".join(lines)


def stubs_section(stubs, overview_dest):
    lines = ["## Tools without a hands-on guide yet", "",
             "Every tool has a page. These have only a short generated summary and the command "
             "reference, until a guide is written.", "",
             "| Code | Tool | What it does |", "|---|---|---|"]
    for page in stubs:
        lines.append(f"| `{page.entry.code}` | [{md_escape(page.entry.display_title)}]({rel(page.dest, overview_dest)}) "
                     f"| {md_escape(catalog.first_sentence(page.entry.description))} |")
    return "\n".join(lines) + "\n\n"


# --- Reference ------------------------------------------------------------------

def command_index(menu):
    """command -> (entry, leaf) for every command the menu runs."""
    return {c: (e, leaf) for e in menu for leaf in e.leaves for c in leaf.commands}


def reference_page(command, helps, in_menu, guide_dest_by_code):
    dest = ref_dest(command)
    # The template only adds a stylesheet that wraps the long --help lines.
    lines = ["---", "template: reference.html", "---", "", f"# `{command}`", ""]
    if command in in_menu:
        entry, leaf = in_menu[command]
        guide = f"[{entry.code} {entry.display_title}]({rel(guide_dest_by_code[entry.code], dest)})"
        if entry.is_workflow:
            lines.append(f"**{md_escape(catalog.clean_title(leaf.title))}** of the {guide} workflow · "
                         f"menu code `{leaf.code}`")
        else:
            lines.append(f"**{md_escape(entry.display_title)}** · menu code `{leaf.code}` · guide: {guide}")
    else:
        lines.append("Not in the `stb-suite` menu: run it directly from the command line.")
    lines += ["", "```text", helps[command]["help"].rstrip("\n"), "```", ""]
    return "\n".join(lines)


def suite_page(menu, guide_dest_by_code):
    dest = ref_dest(catalog.INTERACTIVE_MENU)
    lines = [f"# `{catalog.INTERACTIVE_MENU}`", "",
             "The interactive menu that wraps every tool. Run `stb-suite` with no arguments, then "
             "either browse the six categories or type a tool's code at the main prompt: `1.3` for the "
             "k-grid generator, `4.1.2` for the second stage of the stress-strain workflow.", "",
             "It asks the questions the command-line flags would answer and runs the same underlying "
             "command, so both ways give the same output. Every code below has a guide page and a "
             "command reference page.", ""]
    for cat, (_, name) in catalog.CATEGORIES.items():
        entries = [e for e in menu if e.cat == cat]
        if not entries:
            continue
        lines += [f"## {cat} · {name}", "", "| Code | Tool | Command |", "|---|---|---|"]
        for entry in entries:
            guide = rel(guide_dest_by_code[entry.code], dest)
            for leaf in entry.leaves:
                title = catalog.clean_title(leaf.title) if entry.is_workflow else entry.display_title
                shown = f"{entry.display_title}: {title}" if entry.is_workflow else title
                commands = ", ".join(f"[`{c}`]({rel(ref_dest(c), dest)})" for c in leaf.commands)
                lines.append(f"| `{leaf.code}` | [{md_escape(shown)}]({guide}) | {commands} |")
        lines.append("")
    return "\n".join(lines)


def reference_index(menu, scripts, helps, in_menu):
    dest = "reference/index.md"

    def row(command):
        summary = catalog.first_sentence(catalog.help_summary(helps[command]["help"]))
        code = f"`{in_menu[command][1].code}`" if command in in_menu else "—"
        return f"| [`{command}`]({rel(ref_dest(command), dest)}) | {code} | {md_escape(summary)} |"

    lines = ["# Command reference", "",
             "One page per console command, showing its `--help` text. Each links to the guide for "
             "its tool. [`stb-suite`](stb-suite.md) is the interactive menu; its page maps every "
             "menu code to a command.", ""]
    for cat, (_, name) in catalog.CATEGORIES.items():
        commands = [c for e in menu if e.cat == cat for c in e.commands if c in helps]
        if commands:
            lines += [f"## {cat} · {name}", "", "| Command | Menu code | What it does |", "|---|---|---|"]
            lines += [row(c) for c in commands] + [""]
    others = [c for c in scripts if c not in in_menu and c in helps]
    if others:
        lines += ["## Not in the menu", "", "| Command | Menu code | What it does |", "|---|---|---|"]
        lines += [row(c) for c in others] + [""]
    return "\n".join(lines)


def reference_nav(menu, scripts, helps, in_menu):
    lines = ["* [Overview](index.md)",
             f"* [{catalog.INTERACTIVE_MENU} (menu codes)]({catalog.INTERACTIVE_MENU}.md)"]
    for cat, (_, name) in catalog.CATEGORIES.items():
        commands = [c for e in menu if e.cat == cat for c in e.commands if c in helps]
        if commands:
            lines.append(f"* {cat} · {name}")
            lines += [f"    * [{c}]({c}.md)" for c in commands]
    others = [c for c in scripts if c not in in_menu and c in helps]
    if others:
        lines.append("* Not in the menu")
        lines += [f"    * [{c}]({c}.md)" for c in others]
    return lines


# --- Consistency checks (warnings fail the strict build) --------------------------

def check_consistency(menu, scripts, helps):
    for command in scripts:
        if command != catalog.INTERACTIVE_MENU and command not in helps:
            log.warning("%s has no --help snapshot: run `python docs/dump_help.py`", command)
    for command in helps:
        if command not in scripts:
            log.warning("docs/reference_help.json has %s, which is no longer a console command: "
                        "run `python docs/dump_help.py`", command)
    for entry in menu:
        for leaf in entry.leaves:
            if not leaf.commands:
                log.warning("menu item %s (%s) runs no stb-* command", leaf.code, leaf.title)
            for command in leaf.commands:
                if command not in scripts:
                    log.warning("menu item %s runs %s, which is not in [project.scripts]", leaf.code, command)


def main():
    menu = catalog.load_menu()
    scripts = catalog.load_scripts()
    helps = catalog.load_help()
    check_consistency(menu, scripts, helps)

    pages = collect_pages(menu)
    dest_by_folder = {p.folder: p.dest for p in pages if p.folder}
    guide_dest_by_code = {p.entry.code: p.dest for p in pages}
    in_menu = command_index(menu)
    stubs = [p for p in pages if p.is_stub]

    # Guides
    overview_dest = "guides/index.md"
    overview = (EXAMPLES / "README.md").read_text(encoding="utf-8")
    overview = re.sub(r"\A# Examples\b", "# Guides", overview, count=1)
    marker = "## Adding another example"
    if marker in overview:
        overview = overview.replace(marker, stubs_section(stubs, overview_dest) + marker, 1)
    else:
        overview += "\n" + stubs_section(stubs, overview_dest)
    write_page(overview_dest,
               rewrite_folder_links(github_markdown_to_python_markdown(overview), overview_dest, dest_by_folder),
               "examples/README.md")

    for page in pages:
        if page.is_stub:
            write_page(page.dest, stub_page(page, helps), "docs/gen_pages.py")
        else:
            content = github_markdown_to_python_markdown(page.text)
            write_page(page.dest, rewrite_folder_links(content, page.dest, dest_by_folder) + guide_footer(page),
                       page.edit_source)

    # Sidebar: categories in menu order, pages in numeric code order. The
    # nested-list indent is 4 spaces because Python-Markdown requires it.
    nav, current_cat = ["* [Overview](index.md)"], None
    for page in pages:
        if page.entry.cat != current_cat:
            current_cat = page.entry.cat
            nav.append(f"* {current_cat} · {catalog.CATEGORIES[current_cat][1]}")
        nav.append(f"    * [{page.label}]({posixpath.relpath(page.dest, 'guides')})")
    write_nav("guides/SUMMARY.md", nav)

    # Reference
    write_page("reference/index.md", reference_index(menu, scripts, helps, in_menu), "docs/gen_pages.py")
    write_page(ref_dest(catalog.INTERACTIVE_MENU), suite_page(menu, guide_dest_by_code), "stb-suite/src/stb/stb_suite.py")
    for command in scripts:
        if command in helps:
            module_file = "stb-suite/src/" + helps[command]["module"].replace(".", "/") + ".py"
            write_page(ref_dest(command), reference_page(command, helps, in_menu, guide_dest_by_code), module_file)
    write_nav("reference/SUMMARY.md", reference_nav(menu, scripts, helps, in_menu))


# mkdocs-gen-files runs this file with runpy.run_path (module name "<run_path>");
# importing it elsewhere (e.g. to test the helpers above) builds nothing.
if __name__ in ("__main__", "<run_path>"):
    main()
