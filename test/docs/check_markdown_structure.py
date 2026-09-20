#!/usr/bin/env python3
"""Checks that docs/gen_pages.py keeps every example README's list and
code-block structure intact when MkDocs (Python-Markdown) renders it.

The READMEs are written for GitHub, which renders CommonMark. Each one is
rendered twice -- once as-is with markdown-it-py (CommonMark, i.e. what GitHub
shows), once after gen_pages.github_markdown_to_python_markdown() with the same
Python-Markdown extensions mkdocs.yml uses -- and the sequence of <ul>/<ol>/<li>
tags plus the number of code blocks must match. A mismatch means a nested list
collapsed, or a list fell into the paragraph above it as literal text.

Usage: check_markdown_structure.py   (run from anywhere; exits 1 on any mismatch)
Needs: markdown, pymdown-extensions, markdown-it-py, mkdocs-gen-files
"""
import importlib.util
import re
import sys
from pathlib import Path

import markdown
from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"

# Same extensions as mkdocs.yml (plus the ones MkDocs always enables).
EXTENSIONS = ["tables", "sane_lists", "toc", "pymdownx.highlight",
              "pymdownx.inlinehilite", "pymdownx.superfences"]


def load_gen_pages():
    spec = importlib.util.spec_from_file_location("gen_pages_under_test", ROOT / "docs" / "gen_pages.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def list_tags(html):
    return re.findall(r"</?(?:ul|ol|li)\b", html)


def code_blocks(html):
    return len(re.findall(r"<pre\b", html))


def main():
    gen_pages = load_gen_pages()
    documents = [EXAMPLES / "README.md"] + sorted(EXAMPLES.glob("[0-9]*/README.md"))
    failures = 0
    for path in documents:
        original = path.read_text(encoding="utf-8")
        github = MarkdownIt("commonmark").enable("table").render(original)
        site = markdown.markdown(gen_pages.github_markdown_to_python_markdown(original),
                                 extensions=EXTENSIONS)
        tags_github, tags_site = list_tags(github), list_tags(site)
        if tags_github != tags_site or code_blocks(github) != code_blocks(site):
            failures += 1
            first = next((i for i, (a, b) in enumerate(zip(tags_github, tags_site)) if a != b),
                         min(len(tags_github), len(tags_site)))
            print(f"[FAIL] {path.relative_to(ROOT)}: list tags {len(tags_github)} (GitHub) vs "
                  f"{len(tags_site)} (site), code blocks {code_blocks(github)} vs {code_blocks(site)}, "
                  f"first difference at tag #{first}")
    print(f"{len(documents) - failures}/{len(documents)} documents keep GitHub's list/code structure")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
