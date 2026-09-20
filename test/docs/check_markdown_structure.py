#!/usr/bin/env python3
"""Checks that the docs' Markdown adjustments keep every example README looking
the way GitHub shows it.

The READMEs are written for GitHub, which renders CommonMark. Each one is
rendered twice -- once as-is with markdown-it-py (CommonMark, i.e. what GitHub
shows), once after stbdocs.mdcompat.github_markdown_to_python_markdown() with the
same Python-Markdown extensions mkdocs.yml uses -- and these must match:

- the sequence of <ul>/<ol>/<li> tags (a nested list that collapsed, or a list
  that fell into the paragraph above it as literal text, changes it);
- the tables' shape: number of tables, rows and cells (an unescaped or wrongly
  escaped `|` changes it);
- the number of code blocks;
- the text of every inline `code span` (a stray backslash changes it).

Usage: check_markdown_structure.py   (run from anywhere; exits 1 on any mismatch)
Needs: markdown, pymdown-extensions, markdown-it-py
"""
import html
import re
import sys
from pathlib import Path

import markdown
from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"

sys.path.insert(0, str(ROOT / "docs"))
from stbdocs.mdcompat import github_markdown_to_python_markdown  # noqa: E402

# Same extensions as mkdocs.yml (plus the ones MkDocs always enables).
EXTENSIONS = ["admonition", "tables", "sane_lists", "toc", "pymdownx.highlight",
              "pymdownx.inlinehilite", "pymdownx.superfences"]


def facts(rendered):
    """What must not change between the two renderings."""
    without_blocks = re.sub(r"<pre\b.*?</pre>", "", rendered, flags=re.S)
    return {
        "list tags": re.findall(r"</?(?:ul|ol|li)\b", rendered),
        "table tags": re.findall(r"</?(?:table|tr|th|td)\b", rendered),
        "code blocks": len(re.findall(r"<pre\b", rendered)),
        # whitespace is collapsed: CommonMark turns a line break inside a code
        # span into a space, Python-Markdown keeps it, and a browser shows both
        # the same way
        "inline code": [" ".join(html.unescape(re.sub(r"<[^>]+>", "", m)).split())
                        for m in re.findall(r"<code[^>]*>(.*?)</code>", without_blocks, flags=re.S)],
    }


def first_difference(a, b):
    if isinstance(a, int):
        return f"{a} vs {b}"
    i = next((k for k, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
    return f"{len(a)} (GitHub) vs {len(b)} (site), first difference at #{i}: " \
           f"{a[i] if i < len(a) else '-'!r} vs {b[i] if i < len(b) else '-'!r}"


def main():
    documents = [EXAMPLES / "README.md"] + sorted(EXAMPLES.glob("[0-9]*/README.md"))
    failures = 0
    for path in documents:
        original = path.read_text(encoding="utf-8")
        github = facts(MarkdownIt("commonmark").enable("table").render(original))
        site = facts(markdown.markdown(github_markdown_to_python_markdown(original), extensions=EXTENSIONS))
        for name in github:
            if github[name] != site[name]:
                failures += 1
                print(f"[FAIL] {path.relative_to(ROOT)}: {name}: {first_difference(github[name], site[name])}")
    print(f"{len(documents)} documents checked, {failures} difference(s) from GitHub's rendering")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
