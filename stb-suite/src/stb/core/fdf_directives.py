"""Read-only view of the directives of an `.fdf` text, with SIESTA's own rules.

SIESTA's fdf reader ignores case and the characters `.`, `_` and `-` in a label,
and is *first-occurrence-wins* when a label is repeated. `parse_directives` mirrors
that for the plain `Label value` lines of one or several texts read in order
(e.g. `config_extra.fdf` then `calc.fdf`, which is how a `%include` on top of a
template behaves); the contents of `%block ... %endblock` and comments are
skipped. Used to *report* what a template says (stb-optical's Stage 1, stb-gqca's
Stage 2); never to rewrite it -- the suite edits templates only through a
`config_extra.fdf` sidecar.
"""

import re

_KGRID_LIST_RE = re.compile(r'kgrid[._\-]?monkhorstpack\s*\[\s*(\d+)\s+(\d+)\s+(\d+)\s*\]', re.IGNORECASE)
_KGRID_BLOCK_RE = re.compile(r'%block\s+kgrid[._\-]?monkhorstpack\s*\n(.*?)%endblock', re.IGNORECASE | re.DOTALL)


def normalize_label(label):
    """fdf labels ignore case and '.', '_', '-'."""
    return re.sub(r'[._\-]', '', label.lower())


def parse_directives(*texts):
    """{normalized label: value string} over `texts` in order, first occurrence wins."""
    values, in_block = {}, False
    for text in texts:
        for raw in text.splitlines():
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith('%block'):
                in_block = True
                continue
            if low.startswith('%endblock'):
                in_block = False
                continue
            if in_block or low.startswith('%'):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                values.setdefault(normalize_label(parts[0]), parts[1].strip())
    return values


def parse_kgrid(*texts):
    """The Monkhorst-Pack divisions as an (nx, ny, nz) tuple, from either the list
    form `kgrid.MonkhorstPack [a b c]` or the block form (first three rows, diagonal
    entries); None if the texts define none."""
    combined = "\n".join(texts)
    m = _KGRID_LIST_RE.search(combined)
    if m:
        return tuple(int(v) for v in m.groups())
    mb = _KGRID_BLOCK_RE.search(combined)
    if mb:
        try:
            rows = [r.split() for r in mb.group(1).strip().splitlines()][:3]
            return tuple(int(float(r[i])) for i, r in enumerate(rows))
        except (ValueError, IndexError):
            return None
    return None
