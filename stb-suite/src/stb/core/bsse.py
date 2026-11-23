"""Boys-Bernardi counterpoise (BSSE) ghost-fragment helpers, extracted from
adsorb_bsse.py once stackingfault_bsse.py became a second consumer of the
exact same "turn one fragment of an already-built structure into SIESTA
ghost atoms" need (same extract-on-second-use policy as the rest of
core/). stackingfault_bsse.py was later removed (BSSE correction doesn't
apply to a same-atom-count, same-basis lateral-shift comparison like a
gamma-surface scan -- there's no genuinely different-basis-size fragment
comparison to correct there, unlike an adsorption/interaction energy);
adsorb_bsse.py remains a consumer; species_block_overrides (per-species blocks of a
template re-keyed to one folder's species, ghosts included) is shared by stb-cohesive,
stb-adsorbBsse and stb-chargediffPrep. The underlying issue this module
corrects for, where it IS applicable: SIESTA's localized (PAO) basis set
lets two weakly-bound fragments artificially lower each other's energy by
borrowing basis functions across the gap -- a plane-wave code wouldn't
have this problem at all, but a localized-basis code always does for a
genuine two-different-fragment interaction/binding energy.
"""

import re

from stb.core import structure_io


def make_ghost_variant(base_structure, ghost_labels):
    """Returns a copy of `base_structure` (an FdfStructure) with every atom
    whose CURRENT label is a member of `ghost_labels` (a set of species
    label strings already declared in base_structure.species_meta) turned
    into a ghost species: '<label>_ghost', negative Z, no valence charge,
    same basis (from the same real pseudopotential file) as the real
    element -- SIESTA's standard ghost-atom convention, already used by
    cohesive_energy.py's BSSE ghost clusters there for "one atom's real
    local neighbors". Here it's applied to a whole fragment instead of a
    local neighbor shell -- the standard Boys-Bernardi counterpoise scheme
    for a 2-fragment interaction, exact by construction (no cutoff to
    truncate the correction, unlike cohesive_energy.py's --bsse-cutoff,
    since both fragments here are already complete/finite).

    Selecting by LABEL membership, not an [start, end) index range, is
    deliberate: adsorb_bsse.py's `base_structure` comes from a relaxed
    structure.fdf/.XV, where structure_io.write_fdf has already grouped
    atoms by species -- a slab/adsorbate fragment boundary is NOT
    guaranteed to fall at any particular index once the adsorbate shares an
    element with the slab (see adsorb.py's own '_slab'/'_ads' fragment
    labels, core/adsorption_sites.py::label_fragments -- this is exactly
    the label scheme that makes `ghost_labels` unambiguous here). Z is read
    straight from base_structure.species_meta[symbol] rather than via
    Element(symbol) -- `symbol` here can already be a compound fragment
    label (e.g. 'C_ads'), which pymatgen's Element() would reject.
    """
    species_meta = dict(base_structure.species_meta)
    new_atoms = []
    for symbol, pos in base_structure.atoms:
        if symbol in ghost_labels:
            label = f"{symbol}_ghost"
            if label not in species_meta:
                real_z = abs(species_meta[symbol]['Z'])
                used_ids = {str(info['id']) for info in species_meta.values()}
                next_id = 1
                while str(next_id) in used_ids:
                    next_id += 1
                species_meta[label] = {'id': str(next_id), 'Z': -real_z}
        else:
            label = symbol
        new_atoms.append((label, pos))

    species = list(dict.fromkeys(sym for sym, _ in new_atoms))
    return structure_io.FdfStructure(
        lattice=base_structure.lattice,
        lattice_constant=base_structure.lattice_constant,
        species=species,
        species_meta=species_meta,
        atoms=new_atoms,
        coord_format=base_structure.coord_format,
        raw_lines=[],
    )


def strip_config_extra_include(calc_text, config_extra_file="config_extra.fdf"):
    """Strips the '%include <config_extra_file>' sidecar line
    structure_io.prepend_include added when the prep stage first wrote
    this folder's calc.fdf, recovering the original --calc template text
    -- needed so a BSSE stage can prepend its OWN (different)
    config_extra.fdf (single-point + inherited D3, never the prep stage's
    own relaxation/spin/dipole blocks) without doubling the include line.
    A no-op if the file doesn't start with that exact prefix (defensive;
    e.g. a hand-edited calc.fdf).
    """
    prefix = f"%include {config_extra_file}\n\n"
    if calc_text.startswith(prefix):
        return calc_text[len(prefix):]
    return calc_text


# Blocks keyed by species label (first column, or a header line per species for
# PAO.Basis), compared on SIESTA's normalized name (case, '.', '_', '-' ignored).
SPECIES_BLOCKS = ("paobasissizes", "paopolarizationscheme", "paobasis")
_BLOCK_RE = re.compile(r"^[ \t]*%block[ \t]+(\S+)[ \t]*\n(.*?)^[ \t]*%endblock[^\n]*", re.I | re.S | re.M)


def _norm(label):
    return re.sub(r"[._\-]", "", label.lower())


def _block_entries(name, body):
    """{label: [lines]} of a per-species block, in order. PAO.Basis has one header line per
    species ('Li 2 0.0') followed by its shell lines ('n=2 0 2 P 1', radii, factors), whose
    first token is 'n=..' or a number; the other blocks have one 'label value' line each."""
    entries = {}
    current = None
    for raw in body.splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        first = line.split()[0]
        if _norm(name) == "paobasis":
            is_shell = first.lower().startswith("n=") or re.fullmatch(r"[-+0-9.eEdD]+", first)
            if is_shell and current is not None:
                entries[current].append(line)
                continue
        current = first
        entries.setdefault(current, []).append(line)
    return entries


def species_block_overrides(template_text, labels):
    """config_extra.fdf text that re-keys the template's per-species blocks (PAO.Basis,
    PAO.PolarizationScheme, PAO.BasisSizes) to the species `labels` of ONE folder, or ''
    when the template's blocks already fit it.

    The template is copied verbatim into folders whose species differ (an isolated F atom,
    an F atom with Li ghosts): an entry for an absent species stops SIESTA ('WRONG species
    symbol in PAO.PolarizationScheme: Li'), and a ghost label ('Li_ghost') with no entry of
    its own gets the default basis -- a BSSE ghost then carries a basis other than the
    real atom's, or its polarization orbital fails to generate. Each block is rewritten with
    the entries of the species present, a 'X_ghost' getting the entry of 'X'; config_extra.fdf
    is read first, and SIESTA keeps the first occurrence of a block."""
    labels = list(dict.fromkeys(labels))
    out = []
    for m in _BLOCK_RE.finditer(template_text):
        name, body = m.group(1), m.group(2)
        if _norm(name) not in SPECIES_BLOCKS or "<" in name:
            continue
        entries = _block_entries(name, body)
        new_lines = []
        for label in labels:
            source = label if label in entries else (
                label[:-len("_ghost")] if label.endswith("_ghost") and label[:-len("_ghost")] in entries else None)
            if source is None:
                continue
            lines = list(entries[source])
            head = lines[0].split(None, 1)
            lines[0] = label + (" " + head[1] if len(head) > 1 else "")
            new_lines.extend(lines)
        fits = set(entries) <= set(labels) and not any(
            lb.endswith("_ghost") and lb not in entries and lb[:-len("_ghost")] in entries for lb in labels)
        if fits:
            continue                      # the template's block already fits this folder exactly
        out.append(f"%block {name}\n" + "".join(f"  {ln.strip()}\n" for ln in new_lines)
                   + f"%endblock {name}\n")
    if not out:
        return ""
    return ("# Per-species blocks of the template, re-keyed to this folder's species (absent species\n"
            "# dropped, each X_ghost given the entry of X). Read before the template's own blocks.\n"
            + "".join(out))
