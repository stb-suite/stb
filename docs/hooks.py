"""MkDocs hooks for the documentation site."""

from mkdocs.plugins import event_priority
from mkdocs.structure.files import InclusionLevel


@event_priority(-100)  # after mkdocs-literate-nav, which has already read the nav file
def on_nav(nav, config, files):
    """Keeps the generated sidebar file (guides/SUMMARY.md, see gen_pages.py)
    out of the built site. mkdocs-literate-nav only marks it "not in nav",
    which still renders it as a page and lists it in the sitemap and search."""
    for file in files:
        if file.src_uri.endswith("SUMMARY.md"):
            file.inclusion = InclusionLevel.EXCLUDED
