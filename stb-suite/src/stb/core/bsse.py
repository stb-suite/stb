"""Boys-Bernardi counterpoise (BSSE) ghost-fragment helpers, extracted from
adsorb_bsse.py once stackingfault_bsse.py became a second consumer of the
exact same "turn one fragment of an already-built structure into SIESTA
ghost atoms" need (same extract-on-second-use policy as the rest of
core/). Both consumers correct for the same underlying issue: SIESTA's
localized (PAO) basis set lets two weakly-bound fragments artificially
lower each other's energy by borrowing basis functions across the gap --
a plane-wave code wouldn't have this problem at all, but a localized-basis
code always does for any two-fragment interaction energy.
"""

from pymatgen.core.periodic_table import Element

from stb.core import structure_io


def make_ghost_variant(base_structure, ghost_start, ghost_end):
    """Returns a copy of `base_structure` (an FdfStructure, atoms in a KNOWN
    order such that [ghost_start, ghost_end) is exactly one physical
    fragment -- e.g. adsorb_bsse.py's slab-then-adsorbate order, or
    stackingfault_bsse.py's layer1-then-layer2 order from
    core/heterostructure.py::build_stacked_structure) with those atoms
    turned into ghost species: '<symbol>_ghost' label, negative Z, no
    valence charge, same basis (from the same real pseudopotential file)
    as the real element -- SIESTA's standard ghost-atom convention,
    already used by cohesive_energy.py's BSSE ghost clusters there for
    "one atom's real local neighbors". Here it's applied to a whole
    fragment instead of a local neighbor shell -- the standard
    Boys-Bernardi counterpoise scheme for a 2-fragment interaction, exact
    by construction (no cutoff to truncate the correction, unlike
    cohesive_energy.py's --bsse-cutoff, since both fragments here are
    already complete/finite).
    """
    species_meta = dict(base_structure.species_meta)
    new_atoms = []
    for i, (symbol, pos) in enumerate(base_structure.atoms):
        if ghost_start <= i < ghost_end:
            label = f"{symbol}_ghost"
            if label not in species_meta:
                real_z = Element(symbol).Z
                used_ids = {str(info['id']) for info in species_meta.values()}
                next_id = 1
                while str(next_id) in used_ids:
                    next_id += 1
                species_meta[label] = {'id': str(next_id), 'Z': -abs(real_z)}
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
