"""Single source of truth for reading/writing SIESTA .fdf structure files.

Consolidates 8 independent .fdf parsers that used to live in kgrid.py,
kpath.py, stacking2D.py, elastic_inputs.py, inputfile.py, cohesive_energy.py,
strain.py and translate.py, several of which diverged in how they handled
LatticeConstant and AtomicCoordinatesFormat. translate.py's own fractional/
cartesian coordinate-conversion engine (convert_coordinates, shared by every
format it supports, not just .fdf) stays in translate.py -- write_fdf() here
only knows how to serialize an already-resolved FdfStructure.

Note on units: `FdfStructure.lattice` is the physical lattice (LatticeConstant
already applied, in Angstrom) -- the right thing for tools that need the true
cell (k-grids, pymatgen Structure, symmetry...). Tools that instead rewrite an
existing .fdf's %block LatticeVectors while leaving its LatticeConstant line
untouched (elastic_inputs.py, strain.py) must use `raw_lattice_vectors()`
instead, or they'd silently double-apply the lattice constant on write.

Also hosts a small set of calc.fdf MD-directive helpers (`strip_fdf_comments`,
`read_effective_md_steps`, `read_md_state`, `is_siesta_true`,
`insert_include_after_structure`) shared by strain.py and elastic_inputs.py --
moved here once elastic_inputs.py became a second consumer of the exact same
"read the current MD.TypeOfRun/Steps/VariableCell state without rewriting
anything" need.
"""

from __future__ import annotations

import os
import re
import glob
from dataclasses import dataclass, field

import numpy as np
from pymatgen.core import Lattice, Structure

_SYSTEM_LABEL_RE = re.compile(r"^\s*SystemLabel\s+(\S+)", re.IGNORECASE | re.MULTILINE)


def get_system_label(fdf_path: str) -> str | None:
    """Best-effort SystemLabel from a SIESTA .fdf file: a tolerant regex
    scan, not a full read_fdf() parse -- usable on files read_fdf() would
    reject (missing lattice/species blocks, etc.), since this is only for
    identifying/labeling a run, not building a structure. None if the file
    can't be read or has no SystemLabel line.

    raman_analysis.py's detect_optical_label has its own identical copy of
    this same regex (a different lookup context -- globs a mode_*/ folder
    pattern first, then reads its calc.fdf); not consolidated onto this one
    here to avoid touching unrelated, already-working workflow code.
    """
    try:
        with open(fdf_path) as f:
            text = f.read()
    except OSError:
        return None
    m = _SYSTEM_LABEL_RE.search(text)
    return m.group(1) if m else None


def find_fdf_system_label(directory: str, explicit_label: str | None = None) -> str | None:
    """Returns `explicit_label` if given, else the SystemLabel from the
    first *.fdf (sorted by filename) in `directory` that has one. None if
    neither works. Shared by stb-status and stb-archive, which both need to
    answer "whose calculation is this" from a bare directory before they
    can find its <label>.out/.XV/etc.
    """
    if explicit_label:
        return explicit_label
    for fdf in sorted(glob.glob(os.path.join(directory, "*.fdf"))):
        label = get_system_label(fdf)
        if label:
            return label
    return None

_COORD_FORMAT_FRACTIONAL = {"fractional"}
_COORD_FORMAT_CARTESIAN = {"ang", "notscaledcartesianang"}


@dataclass
class FdfStructure:
    lattice: np.ndarray  # 3x3, Angstrom, LatticeConstant already applied
    lattice_constant: float  # as declared in the file (1.0 if absent)
    species: list[str]  # unique symbols, in %block ChemicalSpeciesLabel order
    species_meta: dict[str, dict]  # {symbol: {'id': str, 'Z': int}}, as declared in the file
    atoms: list[tuple[str, np.ndarray]]  # (symbol, position as given in file)
    coord_format: str  # 'fractional' or 'cartesian'
    raw_lines: list[str] = field(default_factory=list, repr=False)  # original file lines, if read from disk


def _strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def read_fdf(path: str) -> FdfStructure:
    """Parses a SIESTA .fdf file into a FdfStructure.

    Raises FileNotFoundError if the file doesn't exist, ValueError if a
    required block is missing/malformed or AtomicCoordinatesFormat isn't one
    of the values this parser understands (Fractional, Ang/NotScaledCartesianAng).
    """
    with open(path, "r") as f:
        raw_lines = f.readlines()

    lattice_constant = 1.0
    coord_format = None
    in_lattice_block = False
    in_species_block = False
    in_coords_block = False
    lattice_values: list[float] = []
    species_map: dict[str, str] = {}  # index -> symbol
    species_meta: dict[str, dict] = {}  # symbol -> {'id', 'Z'}
    raw_atoms: list[tuple[str, str, str, str]] = []  # (x, y, z, species_index)

    for line in raw_lines:
        cleaned = _strip_comment(line)
        if not cleaned:
            continue
        parts = cleaned.split()
        lower = cleaned.lower()

        if lower.startswith("latticeconstant"):
            try:
                lattice_constant = float(parts[1])
            except (IndexError, ValueError):
                raise ValueError(f"Malformed 'LatticeConstant' line in {path}: {line!r}")
            continue

        if lower.startswith("atomiccoordinatesformat"):
            try:
                fmt = parts[1].lower()
            except IndexError:
                raise ValueError(f"Malformed 'AtomicCoordinatesFormat' line in {path}: {line!r}")
            if fmt in _COORD_FORMAT_FRACTIONAL:
                coord_format = "fractional"
            elif fmt in _COORD_FORMAT_CARTESIAN:
                coord_format = "cartesian"
            else:
                raise ValueError(
                    f"Unsupported AtomicCoordinatesFormat '{parts[1]}' in {path}; "
                    "only Fractional and Ang/NotScaledCartesianAng are supported."
                )
            continue

        if lower == "%block latticevectors":
            in_lattice_block = True
            continue
        if lower == "%block chemicalspecieslabel":
            in_species_block = True
            continue
        if lower == "%block atomiccoordinatesandatomicspecies":
            in_coords_block = True
            continue
        if lower.startswith("%endblock"):
            in_lattice_block = in_species_block = in_coords_block = False
            continue

        if in_lattice_block:
            try:
                lattice_values.extend(float(x) for x in parts[:3])
            except ValueError:
                raise ValueError(f"Malformed line in %block LatticeVectors in {path}: {line!r}")
            continue

        if in_species_block:
            if len(parts) < 3:
                raise ValueError(f"Malformed line in %block ChemicalSpeciesLabel in {path}: {line!r}")
            index, z_str, symbol = parts[0], parts[1], parts[2]
            species_map[index] = symbol
            try:
                species_meta[symbol] = {"id": index, "Z": int(z_str)}
            except ValueError:
                raise ValueError(f"Malformed atomic number in %block ChemicalSpeciesLabel in {path}: {line!r}")
            continue

        if in_coords_block:
            if len(parts) < 4:
                raise ValueError(
                    f"Malformed line in %block AtomicCoordinatesAndAtomicSpecies in {path}: {line!r}"
                )
            raw_atoms.append((parts[0], parts[1], parts[2], parts[3]))
            continue

    if len(lattice_values) != 9:
        raise ValueError(
            f"Expected 9 values in %block LatticeVectors in {path}, found {len(lattice_values)}."
        )
    lattice = np.array(lattice_values).reshape(3, 3) * lattice_constant

    if not species_map:
        raise ValueError(f"%block ChemicalSpeciesLabel not found (or empty) in {path}.")
    if not raw_atoms:
        raise ValueError(f"%block AtomicCoordinatesAndAtomicSpecies not found (or empty) in {path}.")
    if coord_format is None:
        raise ValueError(
            f"'AtomicCoordinatesFormat' not found in {path}; cannot determine "
            "whether coordinates are fractional or cartesian."
        )

    try:
        species = list(dict.fromkeys(species_map[idx] for idx in sorted(species_map, key=int)))
    except ValueError:
        raise ValueError(f"Non-numeric species index in %block ChemicalSpeciesLabel in {path}.")

    atoms: list[tuple[str, np.ndarray]] = []
    for x, y, z, species_index in raw_atoms:
        if species_index not in species_map:
            raise ValueError(f"Species index '{species_index}' not declared in ChemicalSpeciesLabel ({path}).")
        atoms.append((species_map[species_index], np.array([float(x), float(y), float(z)])))

    return FdfStructure(
        lattice=lattice,
        lattice_constant=lattice_constant,
        species=species,
        species_meta=species_meta,
        atoms=atoms,
        coord_format=coord_format,
        raw_lines=raw_lines,
    )


def lattice_only(path_or_structure: str | FdfStructure) -> np.ndarray:
    """Physical lattice matrix (LatticeConstant applied), given a path or an already-read FdfStructure."""
    structure = read_fdf(path_or_structure) if isinstance(path_or_structure, str) else path_or_structure
    return structure.lattice


def raw_lattice_vectors(path_or_structure: str | FdfStructure) -> np.ndarray:
    """Lattice vectors as literally written in %block LatticeVectors (LatticeConstant NOT applied).

    Use this -- not lattice_only()/structure.lattice -- when the result will be
    written back into a file via rewrite_fdf_lattice() while that file's
    LatticeConstant line is left untouched, to avoid double-applying it.
    """
    structure = read_fdf(path_or_structure) if isinstance(path_or_structure, str) else path_or_structure
    return structure.lattice / structure.lattice_constant


def species_list(structure: FdfStructure) -> list[str]:
    """Unique species symbols, in declaration order."""
    return list(structure.species)


def species_dict(structure: FdfStructure) -> dict[str, dict]:
    """{symbol: {'id': declared species index (str), 'Z': declared atomic number (int)}}."""
    return dict(structure.species_meta)


def atom_counts(structure: FdfStructure) -> dict[str, int]:
    """{symbol: number of atoms actually placed in %block AtomicCoordinatesAndAtomicSpecies}.

    A symbol declared in %block ChemicalSpeciesLabel but never placed in the
    coordinates block (a legal, if unusual, SIESTA input -- e.g. leftover from
    editing a template) still shows up here, with a count of 0, rather than
    being silently omitted -- callers that need to KNOW about such a species
    (e.g. to display it) can, while callers that must never require data for
    it (e.g. cohesive_energy.py/cohesive_analysis.py deciding which isolated
    -atom folders/energies are actually needed) should filter on `count > 0`
    themselves. Moved here once cohesive_analysis.py's own copy of this exact
    logic became a second consumer (cohesive_energy.py's prep side needs the
    identical "which declared species are actually used" answer to avoid
    generating/requiring a DFT calculation for a species that never appears
    in the real structure).
    """
    from collections import Counter
    counts = Counter(symbol for symbol, _ in structure.atoms)
    for symbol in structure.species:
        counts.setdefault(symbol, 0)
    return dict(counts)


def to_pymatgen(structure: FdfStructure) -> Structure:
    """Builds a pymatgen Structure using the file's own coordinate format and physical lattice."""
    species = [symbol for symbol, _ in structure.atoms]
    coords = [pos for _, pos in structure.atoms]
    is_cartesian = structure.coord_format == "cartesian"
    return Structure(Lattice(structure.lattice), species, coords, coords_are_cartesian=is_cartesian)


def min_pairwise_distance(structure) -> float | None:
    """Returns the closest distance between any two sites in `structure` (a
    pymatgen Structure or Molecule), or None if it has fewer than 2 sites.

    `structure.distance_matrix` already does the right thing for either input
    type: minimum-image (periodic) distance for a Structure, plain Euclidean
    distance for a Molecule -- so callers don't need to special-case a
    non-periodic cluster housed in an artificial box (its own `distance_matrix`
    would otherwise pick up the box's periodic images as false neighbors).
    Shared by stb-crystalcast's post-build overlap check and (via
    core/structure_checks.py's run_malformation_checks) every tool's own
    [2]/[6] STRUCTURE VALIDATION atom-proximity row, including
    stb-crystalbuilder.
    """
    if len(structure) < 2:
        return None
    dm = structure.distance_matrix
    np.fill_diagonal(dm, np.inf)
    return dm.min()


def read_siesta_structure(path: str, fmt: str) -> Structure:
    """Returns a pymatgen Structure for one of this suite's own SIESTA
    formats: "fdf" (structure input, via read_fdf()/to_pymatgen() above) or
    "struct_out" (post-relaxation output, via ASE -- there's no shared .fdf
    -style reader for that format, it's ASE's own `struct_out` parser).
    Shared by stb-structural and stb-symmetry, the two tools whose analysis
    is meaningful on either a pre- or post-relaxation structure.
    """
    if fmt == "fdf":
        return to_pymatgen(read_fdf(path))
    from ase.io import read as ase_read
    from pymatgen.io.ase import AseAtomsAdaptor
    atoms = ase_read(path, format="struct_out")
    return AseAtomsAdaptor.get_structure(atoms)


def ensure_species_id(species_meta: dict[str, dict], symbol: str) -> dict[str, dict]:
    """Returns a copy of species_meta with `symbol` added, if not already present.

    The new entry gets the next free integer id (as a string, matching the
    rest of species_meta) and 'Z' looked up from pymatgen's periodic table.
    Used whenever a symbol not already declared in a file is introduced --
    e.g. from_pymatgen() below, or a tool adding/substituting an atom of a
    new element.
    """
    if symbol in species_meta:
        return species_meta

    from pymatgen.core.periodic_table import Element

    used_ids = {str(info["id"]) for info in species_meta.values()}
    next_id = 1
    while str(next_id) in used_ids:
        next_id += 1

    new_meta = dict(species_meta)
    new_meta[symbol] = {"id": str(next_id), "Z": Element(symbol).Z}
    return new_meta


def from_pymatgen(
    structure: Structure, species_meta: dict[str, dict] | None = None, coord_format: str = "fractional"
) -> FdfStructure:
    """Builds an FdfStructure from a pymatgen Structure (e.g. after make_supercell()).

    lattice_constant is always 1.0 (the physical lattice is written out as-is,
    matching write_fdf()'s expectations). If species_meta is given (typically
    the original structure's, via species_dict()), its 'id' assignments are
    reused for any symbol present in it, so species numbering stays stable
    across the transformation; any symbol not in it gets the next free id
    (see ensure_species_id()), with 'Z' looked up from the periodic table.
    """
    symbols = [site.specie.symbol for site in structure]
    species = list(dict.fromkeys(symbols))

    meta: dict[str, dict] = {}
    if species_meta:
        for symbol, info in species_meta.items():
            if symbol in species:
                meta[symbol] = dict(info)

    for symbol in species:
        if symbol not in meta:
            meta = ensure_species_id(meta, symbol)

    is_cartesian = coord_format == "cartesian"
    atoms = [
        (site.specie.symbol, np.array(site.coords if is_cartesian else site.frac_coords))
        for site in structure
    ]

    return FdfStructure(
        lattice=np.array(structure.lattice.matrix),
        lattice_constant=1.0,
        species=species,
        species_meta=meta,
        atoms=atoms,
        coord_format=coord_format,
        raw_lines=[],
    )


def write_fdf(structure: FdfStructure, path: str, header_comment: str | list[str] | None = None) -> None:
    """Writes a fresh .fdf file from scratch, built from an FdfStructure.

    Species with zero atoms in structure.atoms are omitted (mirrors the old
    per-caller "if count > 0" filtering). Positions are written in whatever
    structure.coord_format says (fractional or cartesian); this function does
    not convert between the two -- build the FdfStructure already in the
    desired format (e.g. via a fractional/cartesian conversion done by the
    caller) before calling this.

    `header_comment`, if given, overrides the default single-line
    "# automatic create using stb-translate..." header -- a string (one
    line) or a list of strings (one comment line each), for a caller that
    wants to record its own provenance (e.g. stb-mlrelax noting the
    structure was MACE-relaxed, with a couple of convergence details).
    Each line is prefixed with "# " unless it's already commented.
    """
    atoms_by_species: dict[str, list[np.ndarray]] = {}
    for symbol, pos in structure.atoms:
        atoms_by_species.setdefault(symbol, []).append(pos)

    species_with_atoms = [s for s in structure.species if atoms_by_species.get(s)]
    if not species_with_atoms:
        raise ValueError("Cannot write .fdf: no species with at least one atom.")

    raw_lattice = raw_lattice_vectors(structure)
    coord_label = "Fractional" if structure.coord_format == "fractional" else "Ang"
    total_atoms = sum(len(atoms_by_species[s]) for s in species_with_atoms)

    if header_comment is None:
        header_lines = ["# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)\n"]
    else:
        raw_header = [header_comment] if isinstance(header_comment, str) else header_comment
        header_lines = [line if line.startswith("#") else f"# {line}" for line in raw_header]
        header_lines = [f"{line}\n" for line in header_lines]

    lines = [
        *header_lines,
        "\n",
        f"NumberOfSpecies    {len(species_with_atoms)}\n",
        f"NumberofAtoms      {total_atoms}\n\n",
        "%block ChemicalSpeciesLabel\n",
    ]
    for symbol in species_with_atoms:
        meta = structure.species_meta[symbol]
        lines.append(f" {meta['id']}   {meta['Z']}   {symbol}\n")
    lines.append("%endblock ChemicalSpeciesLabel\n\n")
    lines.append(f"LatticeConstant {structure.lattice_constant} Ang\n\n")
    lines.append(f"AtomicCoordinatesFormat  {coord_label}\n\n")
    lines.append("%block LatticeVectors\n")
    for vec in raw_lattice:
        lines.append(f" {vec[0]:.8f}   {vec[1]:.8f}   {vec[2]:.8f}\n")
    lines.append("%endblock LatticeVectors\n\n")
    lines.append("%block AtomicCoordinatesAndAtomicSpecies\n")
    for symbol in species_with_atoms:
        species_id = structure.species_meta[symbol]["id"]
        for pos in atoms_by_species[symbol]:
            lines.append(f"  {pos[0]:.8f}   {pos[1]:.8f}   {pos[2]:.8f}   {species_id}\n")
    lines.append("%endblock AtomicCoordinatesAndAtomicSpecies\n")

    with open(path, "w") as f:
        f.writelines(lines)


def rewrite_fdf_lattice(source_path: str, new_lattice: np.ndarray, out_path: str) -> None:
    """Writes out_path as a copy of source_path with only %block LatticeVectors replaced.

    Everything else (species, positions, basis set, XC block, LatticeConstant
    line, comments...) is preserved verbatim -- new_lattice must already be in
    whatever basis source_path's LatticeConstant expects (see
    raw_lattice_vectors() docstring).
    """
    with open(source_path, "r") as f:
        lines = f.readlines()

    out_lines: list[str] = []
    in_lattice_block = False
    wrote_new_vectors = False
    for line in lines:
        lower = _strip_comment(line).lower()
        if lower == "%block latticevectors":
            in_lattice_block = True
            out_lines.append(line)
            continue
        if in_lattice_block and lower.startswith("%endblock"):
            in_lattice_block = False
            out_lines.append(line)
            continue
        if in_lattice_block:
            if not wrote_new_vectors:
                for vec in new_lattice:
                    out_lines.append(f"  {vec[0]:20.12f}  {vec[1]:20.12f}  {vec[2]:20.12f}\n")
                wrote_new_vectors = True
            continue
        out_lines.append(line)

    if not wrote_new_vectors:
        raise ValueError(f"%block LatticeVectors not found in {source_path}; nothing to rewrite.")

    with open(out_path, "w") as f:
        f.writelines(out_lines)


def rewrite_fdf_positions(source_path: str, new_positions: np.ndarray, out_path: str) -> None:
    """Writes out_path as a copy of source_path with only the atomic positions
    inside %block AtomicCoordinatesAndAtomicSpecies replaced, one-to-one with
    `new_positions` (same order, same count) -- everything else (lattice,
    species, basis, SCF settings, comments...) is preserved verbatim,
    including each atom's original species-index token (and anything after
    it on the line): only the leading x/y/z numbers are replaced.

    `new_positions` must already be in whatever coordinate system
    source_path's AtomicCoordinatesFormat declares (fractional or
    cartesian) -- this function does not convert between them, same
    contract as rewrite_fdf_lattice(). The positions analog of that
    function; added for stb-mlff, which needs to write many rattled/
    perturbed configurations of the SAME reference calculation (same
    basis/pseudos/SCF settings) without hand-reconstructing the whole
    .fdf from scratch (write_fdf() only knows how to write a bare-minimum
    file, not preserve an existing calculation's full setup).
    """
    with open(source_path, "r") as f:
        lines = f.readlines()

    out_lines: list[str] = []
    in_coords_block = False
    atom_i = 0
    for line in lines:
        lower = _strip_comment(line).lower()
        if lower == "%block atomiccoordinatesandatomicspecies":
            in_coords_block = True
            out_lines.append(line)
            continue
        if in_coords_block and lower.startswith("%endblock"):
            in_coords_block = False
            out_lines.append(line)
            continue
        if in_coords_block and _strip_comment(line):
            if atom_i >= len(new_positions):
                raise ValueError(
                    f"{source_path} has more atoms than new_positions ({len(new_positions)})."
                )
            parts = line.split(None, 3)
            suffix = parts[3] if len(parts) > 3 else "\n"
            pos = new_positions[atom_i]
            out_lines.append(f"  {pos[0]:20.12f}  {pos[1]:20.12f}  {pos[2]:20.12f}   {suffix}")
            atom_i += 1
            continue
        out_lines.append(line)

    if atom_i != len(new_positions):
        raise ValueError(
            f"{source_path} has {atom_i} atom(s) but new_positions has {len(new_positions)}."
        )

    with open(out_path, "w") as f:
        f.writelines(out_lines)


def alias_single_atom_species(source_path: str, out_path: str, species: str,
                               atom_index: int, alias_label: str) -> None:
    """Writes out_path as a copy of source_path with exactly one atom
    relabeled from `species` to a new species entry `alias_label` (same Z,
    a fresh species index), appended to %block ChemicalSpeciesLabel.
    Everything else -- lattice, every other atom, comments, basis blocks --
    is preserved verbatim, mirroring rewrite_fdf_lattice()'s approach.

    `atom_index` is 1-based and counts ALL atoms in %block
    AtomicCoordinatesAndAtomicSpecies in file order (the same convention
    stb-defect's --index already uses, see defect.py's parse_index_list) --
    not just atoms of `species`.

    Used by stb-hubbardu to isolate the Cococcioni-de Gironcoli perturbation
    onto a single atom when --species matches more than one atom: SIESTA's
    %block LDAU.proj perturbs an entire species label at once, so with no
    aliasing every atom of a multi-atom species would be perturbed
    simultaneously, measuring a collective response instead of the isolated
    single-site response the supercell method requires.

    Raises ValueError if `species` isn't declared, if `alias_label` already
    exists, if `atom_index` is out of range, or if the atom at `atom_index`
    isn't actually `species`.
    """
    with open(source_path, "r") as f:
        lines = f.readlines()

    species_ids: dict[str, str] = {}
    max_id = 0
    species_z: str | None = None
    in_species_block = False
    for line in lines:
        cleaned = _strip_comment(line)
        lower = cleaned.lower()
        if lower == "%block chemicalspecieslabel":
            in_species_block = True
            continue
        if in_species_block and lower.startswith("%endblock"):
            in_species_block = False
            continue
        if in_species_block and cleaned:
            parts = cleaned.split()
            if len(parts) >= 3:
                idx, z, symbol = parts[0], parts[1], parts[2]
                species_ids[symbol] = idx
                max_id = max(max_id, int(idx))
                if symbol == species:
                    species_z = z

    if species_z is None:
        raise ValueError(f"Species '{species}' not declared in %block ChemicalSpeciesLabel in {source_path}.")
    if alias_label in species_ids:
        raise ValueError(f"Alias species '{alias_label}' already exists in {source_path}; pick a different label.")

    alias_id = str(max_id + 1)
    target_species_id = species_ids[species]

    out_lines: list[str] = []
    in_species_block = False
    in_coords_block = False
    atom_counter = 0
    aliased = False
    for line in lines:
        cleaned = _strip_comment(line)
        lower = cleaned.lower()

        # NumberOfSpecies must match the (now one-longer) %block
        # ChemicalSpeciesLabel -- SIESTA uses it to size the species array
        # when reading that block, so a stale count would silently truncate
        # or misread the block instead of picking up the new alias species.
        if lower.startswith("numberofspecies"):
            parts = cleaned.split()
            try:
                count = int(parts[1])
                out_lines.append(line.replace(parts[1], str(count + 1), 1))
                continue
            except (IndexError, ValueError):
                pass

        if lower == "%block chemicalspecieslabel":
            in_species_block = True
            out_lines.append(line)
            continue
        if in_species_block and lower.startswith("%endblock"):
            out_lines.append(f" {alias_id}   {species_z}   {alias_label}\n")
            in_species_block = False
            out_lines.append(line)
            continue

        if lower == "%block atomiccoordinatesandatomicspecies":
            in_coords_block = True
            out_lines.append(line)
            continue
        if in_coords_block and lower.startswith("%endblock"):
            in_coords_block = False
            out_lines.append(line)
            continue
        if in_coords_block and cleaned:
            atom_counter += 1
            if atom_counter == atom_index:
                parts = cleaned.split()
                if len(parts) < 4 or parts[3] != target_species_id:
                    raise ValueError(
                        f"Atom {atom_index} in {source_path} is not species '{species}' "
                        f"(expected species index {target_species_id})."
                    )
                parts[3] = alias_id
                out_lines.append("  " + "   ".join(parts) + "\n")
                aliased = True
                continue

        out_lines.append(line)

    if not aliased:
        raise ValueError(
            f"Atom index {atom_index} is out of range in {source_path} "
            f"(found {atom_counter} atom(s))."
        )

    with open(out_path, "w") as f:
        f.writelines(out_lines)


# ==========================================
# calc.fdf MD-directive helpers -- shared by strain.py and elastic_inputs.py
# (moved here once elastic_inputs.py became a second consumer of the exact
# same "read the current MD.TypeOfRun/Steps/VariableCell state, without
# rewriting anything" need, same extract-on-second-use policy as the rest
# of this module).
# ==========================================

_MD_TYPEOFRUN_VALUE_RE = re.compile(r'MD\.TypeOfRun\s+(\S+)', re.IGNORECASE)
_MD_STEPS_VALUE_RE = re.compile(r'MD\.Steps\s+(\d+)', re.IGNORECASE)
_MD_NUMCGSTEPS_VALUE_RE = re.compile(r'MD\.NumCGsteps\s+(\d+)', re.IGNORECASE)
_MD_VARIABLECELL_VALUE_RE = re.compile(r'MD\.VariableCell\s+(\S+)', re.IGNORECASE)
_SIESTA_TRUE_VALUES = {'t', 'true', '.true.', 'yes'}


def strip_fdf_comments(calc_text: str) -> str:
    """Removes everything from the first '#' onward on every line (SIESTA's
    own fdf comment convention -- also used for trailing inline comments
    after a real value, e.g. 'DM.UseSaveDM .true.  #(...)'). Used before
    scanning calc_text for a directive's CURRENT value -- otherwise a
    plain-English comment merely mentioning a directive's name (e.g. '##
    MD.VariableCell should stay true', a real, caught-live bug while
    writing strain.py's own example) could be misread as the actual value
    if that comment happens to appear earlier in the file than the real
    directive line."""
    return "\n".join(line.split('#', 1)[0] for line in calc_text.splitlines())


def read_effective_md_steps(calc_text: str) -> tuple[str | None, str | None]:
    """Returns (key_name, value) for the relaxation step count -- prefers
    'MD.Steps' (what a real relaxation calc.fdf and stb-inputfile's own
    generated template both use), falling back to 'MD.NumCGsteps' (an
    older/alternate spelling some templates may still use) if 'MD.Steps'
    is absent. (None, None) if neither is present."""
    calc_text = strip_fdf_comments(calc_text)
    steps_match = _MD_STEPS_VALUE_RE.search(calc_text)
    if steps_match:
        return 'MD.Steps', steps_match.group(1)
    numcg_match = _MD_NUMCGSTEPS_VALUE_RE.search(calc_text)
    if numcg_match:
        return 'MD.NumCGsteps', numcg_match.group(1)
    return None, None


def read_md_state(calc_text: str) -> dict:
    """Read-only scan of calc_text's CURRENT MD.TypeOfRun/MD.Steps (or
    MD.NumCGsteps)/MD.VariableCell directives, for reporting purposes --
    never rewrites anything itself, purely informational/advisory. A value
    is None when the tag is absent (a normal template, not a malformed
    one)."""
    stripped = strip_fdf_comments(calc_text)
    run_match = _MD_TYPEOFRUN_VALUE_RE.search(stripped)
    varcell_match = _MD_VARIABLECELL_VALUE_RE.search(stripped)
    steps_key, steps_value = read_effective_md_steps(calc_text)
    return {
        'typeofrun': run_match.group(1) if run_match else None,
        'steps_key': steps_key,
        'steps': steps_value,
        'variablecell': varcell_match.group(1) if varcell_match else None,
    }


def is_siesta_true(value: str | None) -> bool:
    """True if value.lower() is one of SIESTA's boolean-true spellings."""
    return value is not None and value.lower() in _SIESTA_TRUE_VALUES


def insert_include_after_structure(calc_text: str, structure_basename: str, include_name: str) -> str:
    """Inserts '%include <include_name>' right after the existing
    '%include <structure_basename>' line -- the exact position the
    Geometry.Constraints block already occupies (commented out, as a
    template) in strain.py's own reference calc.fdf files. Falls back to
    prepending at the very top of the file if no such line is found.
    --calc itself is otherwise untouched -- no other directive is
    rewritten here.

    Only safe to use when `include_name`'s own directives do NOT already
    appear earlier in calc_text: SIESTA's fdf reader is first-occurrence
    -wins for duplicate labels (verified in hubbardu.py::write_run_folder),
    so an include inserted here would be silently ignored for any
    directive the base file already sets earlier in the file. strain.py's
    own config_extra.fdf (a %block Geometry.Constraints) is safe because
    that block name never pre-exists in a real calc.fdf; a caller wanting
    to override a directive (e.g. MD.TypeOfRun/MD.Steps) that commonly
    IS already present in a real template should prepend instead (see
    elastic_inputs.py's own local prepend_include)."""
    include_line = f"%include {include_name}"
    pattern = re.compile(
        r'^([ \t]*%include[ \t]+' + re.escape(structure_basename) + r'[ \t]*)$',
        re.IGNORECASE | re.MULTILINE)
    new_text, count = pattern.subn(r'\1\n' + include_line, calc_text, count=1)
    if count == 0:
        new_text = include_line + "\n\n" + calc_text
    return new_text
