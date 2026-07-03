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
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pymatgen.core import Lattice, Structure

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

    species = list(dict.fromkeys(species_map[idx] for idx in sorted(species_map, key=int)))

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


def to_pymatgen(structure: FdfStructure) -> Structure:
    """Builds a pymatgen Structure using the file's own coordinate format and physical lattice."""
    species = [symbol for symbol, _ in structure.atoms]
    coords = [pos for _, pos in structure.atoms]
    is_cartesian = structure.coord_format == "cartesian"
    return Structure(Lattice(structure.lattice), species, coords, coords_are_cartesian=is_cartesian)


def write_fdf(structure: FdfStructure, path: str) -> None:
    """Writes a fresh .fdf file from scratch, built from an FdfStructure.

    Species with zero atoms in structure.atoms are omitted (mirrors the old
    per-caller "if count > 0" filtering). Positions are written in whatever
    structure.coord_format says (fractional or cartesian); this function does
    not convert between the two -- build the FdfStructure already in the
    desired format (e.g. via a fractional/cartesian conversion done by the
    caller) before calling this.
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

    lines = [
        "# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)\n\n",
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
