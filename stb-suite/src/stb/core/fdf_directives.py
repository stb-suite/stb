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


def filter_caution(*texts):
    """Warning text when the fdf texts (in SIESTA's reading order) turn the orbital
    filter on (FilterCutoff or FilterTol), else None. The filter changes the basis:
    with 'FilterCutoff 150 Ry' the fit residual of E(V) was 2.13 (Cu) and 1.43 (Al)
    meV/atom against 0.11 and 0.07 without it, and Cu's V0 moved from 11.904 to
    11.715 Ang^3 and B0 from 141 to 128 GPa (all-electron: 11.951, 141.3)."""
    values = parse_directives(*texts)
    active = [f"{name} {values[key]}" for name, key in (("FilterCutoff", "filtercutoff"), ("FilterTol", "filtertol"))
              if key in values]
    if not active:
        return None
    return (f"The orbital filter is ON ({', '.join(active)}): SIESTA filters the basis orbitals, which changes "
            "the basis itself and adds noise to energy differences between cells of different size -- E(V), "
            "elastic constants, convergence scans (Cu, FilterCutoff 150 Ry: fit residual 2.1 instead of 0.1 "
            "meV/atom, V0 -1.6 %, B0 -9 %; Al, Ag affected too; Si, graphene, hBN not). Remove the line and "
            "converge Mesh.Cutoff instead; 2.0.0's stb-inputfile/stb-cohesive templates wrote it by default.")


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


def describe_template_kgrid(*texts):
    """How the fdf texts set the SCF k-sampling, as a short string, or None when they set none (SIESTA then
    samples Gamma only): the Monkhorst-Pack divisions (list or block form), 'kgrid.Cutoff <value>' or
    'kgrid.File <name>'. A tool that forces a k-grid only "when the template has none" must test this, not
    parse_kgrid alone: stb-cohesive --calc once wrote a 10x10x10 block over a template's 'kgrid.Cutoff 15.0 Ang'
    (the block takes precedence) and SIESTA ran with 560 k-points."""
    grid = parse_kgrid(*texts)
    if grid is not None:
        return "x".join(str(v) for v in grid)
    values = parse_directives(*texts)
    for name, key in (("kgrid.File", "kgridfile"), ("kgrid.Cutoff", "kgridcutoff")):
        if key in values:
            return f"{name} {values[key]}"
    return None
