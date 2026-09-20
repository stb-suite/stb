"""Generates the documentation site's guide pages from `examples/`.

`examples/` is the single source of truth: nothing is copied into `docs/`.
This script runs at build time (mkdocs-gen-files) and emits virtual pages:

- `guides/index.md`           <- examples/README.md
- `guides/<category>/<code>-<name>.md` <- examples/<code>-<name>/README.md
- `guides/SUMMARY.md`         <- the sidebar nav (mkdocs-literate-nav), ordered
                                 numerically by menu code (so 1.2 < 1.10)

Each example folder is named after the `stb-suite` menu's dotted code
(`<category>.<item>-<name>`), so the category and ordering come from the
folder name itself, and the sidebar label from the README's own H1. A new
`examples/<code>-<name>/README.md` shows up on the next build with no
change here.

The READMEs are written for GitHub and rendered here by Python-Markdown, so
each page is adjusted on the way through:

- Relative links between example folders (`../4.3-cohesive/`,
  `1.1-stb-inputfile/`) become links between the generated pages.
- List structure is normalized (`github_markdown_to_python_markdown`): the two
  renderers disagree on how deep a nested list is indented and on whether a
  list may follow a paragraph line directly. The tests in `test/docs/` check
  that every page keeps the list and code-block structure GitHub gives it.
- Each page gets a footer pointing back at its folder on GitHub, and an "edit
  this page" link to the README it came from.
"""

import logging
import posixpath
import re
from pathlib import Path

import mkdocs_gen_files

log = logging.getLogger("mkdocs.plugins.gen_pages")

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"

REPO_URL = "https://github.com/stb-suite/stb"
BRANCH = "main"

# Same 6 categories, in the same order, as stb_suite.py's main menu.
CATEGORIES = {
    1: ("inputs", "Inputs"),
    2: ("structures", "Structures"),
    3: ("analysis", "Analysis"),
    4: ("workflows", "Workflows"),
    5: ("ml-simulations", "ML Simulations"),
    6: ("utils", "Utils"),
}

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


class Guide:
    def __init__(self, folder):
        self.folder = folder
        m = FOLDER_RE.match(folder)
        self.cat, self.item = int(m.group(1)), int(m.group(2))
        self.code = f"{self.cat}.{self.item}"
        self.src = EXAMPLES / folder / "README.md"
        self.text = self.src.read_text(encoding="utf-8")
        slug, self.category_name = CATEGORIES.get(self.cat, (f"category-{self.cat}", f"Category {self.cat}"))
        self.dest = f"guides/{self.cat}-{slug}/{folder}.md"
        self.title = self._title()

    def _h1(self):
        in_fence = False
        for line in self.text.splitlines():
            if line.lstrip().startswith(("```", "~~~")):
                in_fence = not in_fence
            elif not in_fence and line.startswith("# "):
                return line[2:].strip()
        return None

    def _title(self):
        h1 = self._h1()
        if h1 is None:
            log.warning("%s has no H1; using the folder name as its title", self.src)
            return self.folder
        for pattern in (_H1_CMD_FIRST, _H1_CMD_LAST):
            m = pattern.match(h1)
            if m:
                title = re.sub(r"^Workflow: ", "", m.group("title"))
                return title.replace("`", "")
        log.warning("%s: unrecognized H1 format %r; using it as-is", self.src, h1)
        return h1.replace("`", "")

    @property
    def nav_label(self):
        return f"{self.code} {self.title}"


def collect_guides():
    guides = []
    for path in sorted(EXAMPLES.iterdir()):
        if path.is_dir() and FOLDER_RE.match(path.name):
            if (path / "README.md").is_file():
                guides.append(Guide(path.name))
            else:
                log.warning("examples/%s has no README.md; skipped", path.name)
    guides.sort(key=lambda g: (g.cat, g.item))
    return guides


def rewrite_folder_links(text, source_page, dest_by_folder):
    """Turns links between example folders into links between the generated
    pages, relative to `source_page`. Fenced code blocks are left untouched."""
    source_dir = posixpath.dirname(source_page)

    def replace(match):
        folder, anchor = match.group(1), match.group(2) or ""
        dest = dest_by_folder.get(folder)
        if dest is None:
            log.warning("%s links to unknown example folder %r", source_page, folder)
            return match.group(0)
        return f"]({posixpath.relpath(dest, source_dir)}{anchor})"

    out, in_fence = [], False
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith(("```", "~~~")):
            in_fence = not in_fence
        out.append(line if in_fence else _FOLDER_LINK.sub(replace, line))
    return "".join(out)


# A list item marker with the whitespace after it: `- `, `* `, `+ `, `1. `, `1) `.
_LIST_START = re.compile(r"^(\s*)([-*+]|\d+[.)])(\s+)\S")
# Lines that are never plain paragraph text: heading, table row, blockquote, HTML, rule.
_NOT_PARAGRAPH = re.compile(r"^\s*(?:#|\||>|<|(?:[-*_]\s*){3,}$)")
_QUOTE_LINE = re.compile(r"^ {0,3}>")
_FENCE = ("```", "~~~")


def _indent(line):
    return len(line) - len(line.lstrip())


def _content_offset(match):
    """Column where a list item's text starts (what its continuation lines and
    nested lists must be indented to, per CommonMark)."""
    lead, marker, gap = match.group(1), match.group(2), match.group(3)
    return len(lead) + len(marker) + (len(gap) if len(gap) <= 4 else 1)


def _normalize_lists(text):
    """Rewrites the list structure of GitHub-flavoured Markdown (CommonMark) so
    that Python-Markdown (MkDocs) renders it the same way. Two differences:

    - Nesting. CommonMark nests a list under an item when it is indented to
      the item's text (2 spaces under `- `, 3 under `1. `); Python-Markdown
      only nests at 4 spaces per level. Every list level, and every
      continuation line inside an item, is re-indented to 4 spaces per level.
    - Interrupting a paragraph. CommonMark lets a top-level list follow a
      paragraph line directly; Python-Markdown then folds it into the
      paragraph as literal text, so a blank line is inserted first. (An
      ordered list interrupts a paragraph only if it starts at 1.)

    Fenced code blocks are shifted along with the item they belong to and
    otherwise left alone."""
    out, stack = [], []          # stack: content offsets of the open list items
    prev_kind, prev_blank = None, True
    fence, fence_delta = None, 0
    for raw in text.splitlines(keepends=True):
        stripped = raw.strip()
        if fence is not None:
            if fence_delta >= 0:
                out.append(" " * fence_delta + raw if stripped else raw)
            else:
                out.append(raw[min(-fence_delta, _indent(raw)):] if stripped else raw)
            if stripped.startswith(fence):
                fence = None
                prev_kind, prev_blank = "special", False
            continue
        if not stripped:
            out.append(raw)
            prev_blank = True
            continue

        indent = _indent(raw)
        in_list_before = bool(stack)
        while stack and indent < stack[-1]:
            stack.pop()
        depth = len(stack)
        match = _LIST_START.match(raw)

        if match:
            interrupts = depth == 0 and prev_kind == "para" \
                and (not match.group(2)[0].isdigit() or match.group(2)[:-1] == "1")
            if interrupts:
                out.append("\n")
            out.append(" " * (4 * depth) + raw.lstrip())
            stack.append(_content_offset(match))
            prev_kind = "list"
        else:
            new_indent = indent if depth == 0 else 4 * depth + max(indent - stack[-1], 0)
            out.append(" " * new_indent + raw.lstrip())
            if stripped.startswith(_FENCE):
                fence, fence_delta = stripped[:3], new_indent - indent
                prev_kind = "special"
            elif depth > 0 or (in_list_before and not prev_blank):
                prev_kind = "list"    # item content, or a lazy continuation of one
            elif _NOT_PARAGRAPH.match(raw):
                prev_kind = "special"
            else:
                prev_kind = "para"
        prev_blank = False
    return "".join(out)


def github_markdown_to_python_markdown(text):
    """Applies `_normalize_lists` to the document and, separately, to the body
    of every blockquote (a list after a paragraph line inside a `>` block needs
    the same fix). Fenced code blocks, including any `>` inside them, are
    left untouched."""
    out, run, quoting, in_fence = [], [], False, False

    def flush():
        if quoting:
            inner = [re.sub(r"^ {0,3}> ?", "", line, count=1) for line in run]
            for line in _normalize_lists("".join(inner)).splitlines(keepends=True):
                out.append("> " + line if line.strip() else ">\n")
        else:
            out.append(_normalize_lists("".join(run)))

    for line in text.splitlines(keepends=True):
        is_quote = not in_fence and bool(_QUOTE_LINE.match(line))
        if run and is_quote != quoting:
            flush()
            run = []
        quoting = is_quote
        run.append(line)
        if not is_quote and line.lstrip().startswith(_FENCE):
            in_fence = not in_fence
    if run:
        flush()
    return "".join(out)


def footer(guide):
    tree = f"{REPO_URL}/tree/{BRANCH}/examples/{guide.folder}"
    parts = [f"Example folder (input fixtures and guided script) on GitHub: [`examples/{guide.folder}/`]({tree})"]
    scripts = sorted((EXAMPLES / guide.folder).glob("example_*.sh"))
    if scripts:
        blob = f"{REPO_URL}/blob/{BRANCH}/examples/{guide.folder}/{scripts[0].name}"
        parts.append(f"guided script [`{scripts[0].name}`]({blob})")
    return "\n\n---\n\n*" + " · ".join(parts) + "*\n"


def write_page(dest, content, edit_source):
    with mkdocs_gen_files.open(dest, "w") as f:
        f.write(content)
    # edit_uri is rooted at docs/, so reach the repo root with "../"
    mkdocs_gen_files.set_edit_path(dest, f"../{edit_source}")


def main():
    guides = collect_guides()
    dest_by_folder = {g.folder: g.dest for g in guides}

    # Overview page: examples/README.md (links to every folder + shared conventions).
    overview_dest = "guides/index.md"
    overview = (EXAMPLES / "README.md").read_text(encoding="utf-8")
    overview = re.sub(r"\A# Examples\b", "# Guides", overview, count=1)
    write_page(overview_dest,
               rewrite_folder_links(github_markdown_to_python_markdown(overview),
                                    overview_dest, dest_by_folder),
               "examples/README.md")

    for g in guides:
        content = rewrite_folder_links(github_markdown_to_python_markdown(g.text),
                                       g.dest, dest_by_folder) + footer(g)
        write_page(g.dest, content, f"examples/{g.folder}/README.md")

    # Sidebar: categories in menu order, pages in numeric code order. The
    # nested-list indent is 4 spaces because Python-Markdown requires it.
    nav = ["* [Overview](index.md)\n"]
    current_cat = None
    for g in guides:
        if g.cat != current_cat:
            current_cat = g.cat
            nav.append(f"* {g.cat} · {g.category_name}\n")
        nav.append(f"    * [{g.nav_label}]({posixpath.relpath(g.dest, 'guides')})\n")
    with mkdocs_gen_files.open("guides/SUMMARY.md", "w") as f:
        f.writelines(nav)


# mkdocs-gen-files runs this file with runpy.run_path (module name "<run_path>");
# importing it elsewhere (e.g. to test the helpers above) builds nothing.
if __name__ in ("__main__", "<run_path>"):
    main()
