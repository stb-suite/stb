"""Boys-Bernardi counterpoise (BSSE) ghost-fragment helpers, extracted from
adsorb_bsse.py once stackingfault_bsse.py became a second consumer of the
exact same "turn one fragment of an already-built structure into SIESTA
ghost atoms" need (same extract-on-second-use policy as the rest of
core/). stackingfault_bsse.py was later removed (BSSE correction doesn't
apply to a same-atom-count, same-basis lateral-shift comparison like a
gamma-surface scan -- there's no genuinely different-basis-size fragment
comparison to correct there, unlike an adsorption/interaction energy);
adsorb_bsse.py remains the sole consumer. The underlying issue this module
corrects for, where it IS applicable: SIESTA's localized (PAO) basis set
lets two weakly-bound fragments artificially lower each other's energy by
borrowing basis functions across the gap -- a plane-wave code wouldn't
have this problem at all, but a localized-basis code always does for a
genuine two-different-fragment interaction/binding energy.
"""

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
