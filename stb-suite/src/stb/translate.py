#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
from stb import __version__ as VERSION

import os
import sys
import glob
from datetime import datetime
import argparse
from typing import List, Dict
from ase.io import read as ase_read
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections
from stb.core.symmetry import crystal_system
from stb.core.ase_view import view_structure_interactive

# Import libraries:

# Supported formats
##### MODIFICADO #####
# Adicionado "fdf" aos formatos de entrada
INPUT_FORMATS = {"poscar", "cif", "siesta", "xyz", "fhi", "dftb", "xsf", "fdf"}
OUTPUT_FORMATS = {"cif","xyz", "poscar", "fdf", "dftb", "xsf", "fhi", "lammps"}

# Extension -> --in-format guess, used when -if/--in-format is omitted.
EXTENSION_FORMAT_MAP = {
    ".poscar": "poscar", ".vasp": "poscar",
    ".cif": "cif",
    ".xyz": "xyz",
    ".fdf": "fdf",
    ".xsf": "xsf",
    ".gen": "dftb",
    ".in": "fhi",
    ".struct_out": "siesta",
}
# Filenames with no (or a non-informative) extension that VASP/SIESTA use by convention.
BASENAME_FORMAT_MAP = {
    "poscar": "poscar", "contcar": "poscar",
}

# --out-format -> file extension used to auto-name each output in batch mode
# (--in-file as a glob pattern matching more than one file).
OUTPUT_EXTENSION_MAP = {
    "poscar": "poscar", "cif": "cif", "xyz": "xyz", "fdf": "fdf",
    "dftb": "gen", "xsf": "xsf", "fhi": "in", "lammps": "data",
}

# --out-file extension -> --out-format guess, used when -of/--out-format is omitted.
OUTPUT_EXTENSION_GUESS = {
    ".cif": "cif", ".xyz": "xyz", ".fdf": "fdf", ".xsf": "xsf",
    ".poscar": "poscar", ".vasp": "poscar",
    ".gen": "dftb", ".dftb": "dftb",
    ".in": "fhi", ".fhi": "fhi",
    ".data": "lammps", ".lmp": "lammps", ".lammps": "lammps",
}

# Report filename written alongside the output (--out-dir in batch mode, or
# --out-file's own directory otherwise) -- same print_dual-to-screen-and-file
# convention as the workflow analysis stages (e.g. stb-ramanAnalysis's
# raman_stage3.txt). Not written for --dry-run, which by definition has no
# disk side effects.
REPORT_FILE = "stb_translate_report.txt"
CIF_SYMPREC = 0.01          # --cif-symprec; None = P1, positions exactly as read


def guess_format_from_filename(path):
    """Best-effort --in-format guess from a filename; None if not recognized."""
    base = os.path.basename(path).lower()
    if base in BASENAME_FORMAT_MAP:
        return BASENAME_FORMAT_MAP[base]
    _, ext = os.path.splitext(base)
    return EXTENSION_FORMAT_MAP.get(ext)


# Dictionary with all periodic elements table
def periodic_table():
    element_atomicnumber = {
        "H": 1,     # Hidrogênio
        "He": 2,    # Hélio
        "Li": 3,    # Lítio
        "Be": 4,    # Berílio
        "B": 5,     # Boro
        "C": 6,     # Carbono
        "N": 7,     # Nitrogênio
        "O": 8,     # Oxigênio
        "F": 9,     # Flúor
        "Ne": 10,   # Neônio
        "Na": 11,   # Sódio
        "Mg": 12,   # Magnésio
        "Al": 13,   # Alumínio
        "Si": 14,   # Silício
        "P": 15,    # Fósforo
        "S": 16,    # Enxofre
        "Cl": 17,   # Cloro
        "Ar": 18,   # Argônio
        "K": 19,    # Potássio
        "Ca": 20,   # Cálcio
        "Sc": 21,   # Escândio
        "Ti": 22,   # Titânio
        "V": 23,    # Vanádio
        "Cr": 24,   # Cromo
        "Mn": 25,   # Manganês
        "Fe": 26,   # Ferro
        "Co": 27,   # Cobalto
        "Ni": 28,   # Níquel
        "Cu": 29,   # Cobre
        "Zn": 30,   # Zinco
        "Ga": 31,   # Gálio
        "Ge": 32,   # Germanio
        "As": 33,   # Arsênio
        "Se": 34,   # Selênio
        "Br": 35,   # Bromo
        "Kr": 36,   # Kriptônio
        "Rb": 37,   # Rubídio
        "Sr": 38,   # Estrôncio
        "Y": 39,    # Ítrio
        "Zr": 40,   # Zircônio
        "Nb": 41,   # Nióbio
        "Mo": 42,   # Molibdênio
        "Tc": 43,   # Tecnécio
        "Ru": 44,   # Rutênio
        "Rh": 45,   # Ródio
        "Pd": 46,   # Paládio
        "Ag": 47,   # Prata
        "Cd": 48,   # Cádmio
        "In": 49,   # Índio
        "Sn": 50,   # Estanho
        "Sb": 51,   # Antimônio
        "Te": 52,   # Telúrio
        "I": 53,    # Iodo
        "Xe": 54,   # Xenônio
        "Cs": 55,   # Césio
        "Ba": 56,   # Bário
        "La": 57,   # Lantânio
        "Ce": 58,   # Cério
        "Pr": 59,   # Praseodímio
        "Nd": 60,   # Neodímio
        "Pm": 61,   # Promécio
        "Sm": 62,   # Samário
        "Eu": 63,   # Európio
        "Gd": 64,   # Gadolínio
        "Tb": 65,   # Térbio
        "Dy": 66,   # Disprósio
        "Ho": 67,   # Holmium
        "Er": 68,   # Érbio
        "Tm": 69,   # Túlio
        "Yb": 70,   # Itérbio
        "Lu": 71,   # Lutécio
        "Hf": 72,   # Háfnio
        "Ta": 73,   # Tântalo
        "W": 74,    # Wolframio
        "Re": 75,   # Rênio
        "Os": 76,   # Ósmio
        "Ir": 77,   # Irídio
        "Pt": 78,   # Platina
        "Au": 79,   # Ouro
        "Hg": 80,   # Mercúrio
        "Tl": 81,   # Tálio
        "Pb": 82,   # Chumbo
        "Bi": 83,   # Bismuto
        "Po": 84,   # Polônio
        "At": 85,   # Astato
        "Rn": 86,   # Radônio
        "Fr": 87,   # Frâncio
        "Ra": 88,   # Radônio
        "Ac": 89,   # Actínio
        "Th": 90,   # Tório
        "Pa": 91,   # Protactínio
        "U": 92,    # Urânio
        "Np": 93,   # Netúnio
        "Pu": 94,   # Plutônio
        "Am": 95,   # Amerício
        "Cm": 96,   # Curió
        "Bk": 97,   # Berquélio
        "Cf": 98,   # Califórnio
        "Es": 99,   # Einstênio
        "Fm": 100,  # Férmio
        "Md": 101,  # Mendelevio
        "No": 102,  # Nobelio
        "Lr": 103,  # Laurêncio
        "Rf": 104,  # Rutherfórdio
        "Db": 105,  # Dúbnio
        "Sg": 106,  # Seabórgio
        "Bh": 107,  # Bóhrio
        "Hs": 108,  # Hassio
        "Mt": 109,  # Meitnério
        "Ds": 110,  # Darmstádio
        "Rg": 111,  # Roentgênio
        "Cn": 112,  # Copernício
        "Nh": 113,  # Nihônio
        "Fl": 114,  # Fleróvio
        "Mc": 115,  # Moscóvio
        "Lv": 116,  # Livermório
        "Ts": 117,  # Tenesso
        "Og": 118   # Oganessônio
    }

    # Dictionary with all atomic numbers table (invert Dictionary)
    atomicnumber_element = {str(v): k for k, v in element_atomicnumber.items()}
    return element_atomicnumber, atomicnumber_element


def readfile(filedata):
    with open(filedata, 'r') as fil:
        data = [line.split() for line in fil if line.strip()
                ]
    return data


def dic_atoms_position(atomsposition):
    dic_atomspos = {}
    for elt in atomsposition:
        if elt[0] not in dic_atomspos:
            dic_atomspos[elt[0]] = []
        dic_atomspos[elt[0]].append([elt[1], elt[2], elt[3]])
    return dic_atomspos


############# Extract Data Functions ######################

# This Function define the number of atoms types for xyz file
def getatomsandvectors_xyz(dataxyz, vectors):
    """vectors: the 3 lattice vectors as [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]
    (already Cartesian Angstrom strings) -- see --lattice-vectors in main()."""
    element, atomicnumber = periodic_table()
    xyz = readfile(dataxyz)[2:]
    latticeparameter = "1.0"
    typevectors = 'Cartesian'
    # Atoms position
    atoms = []
    j = 1
    getatoms = []
    for i in range(len(xyz)):
        if xyz[i][0] not in atoms:
            atoms.append(xyz[i][0])
            getatoms.append([j, element[xyz[i][0]], xyz[i][0]])
            j = j+1
    atomsposition = xyz
    for elem in getatoms:
        cont = 0
        for elxyz in xyz:
            if elem[2] == elxyz[0]:
                cont = cont+1
        elem.append(str(cont))
    atomic_position = dic_atoms_position(atomsposition)
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


# This Function define the number of atoms types for vasp file
def getatomsandvectors_vasp(poscar):
    element, atomicnumber = periodic_table()
    datavasp = readfile(poscar)
    latticeparameter = datavasp[1][0]
    vectors = datavasp[2:5]
    # a negative scale factor is VASP's way of giving the cell VOLUME: the factor is then (|V| / det(vectors))^(1/3)
    if float(latticeparameter) < 0:
        det = abs(np.linalg.det(np.array(vectors, dtype=float)))
        latticeparameter = f"{(abs(float(latticeparameter)) / det) ** (1 / 3):.10f}"
    # Optional "Selective dynamics" line (VASP only requires the first letter to
    # match) shifts the coordinate-type line and the start of positions by one.
    coord_line_idx = 8 if datavasp[7][0].lower().startswith('s') else 7
    # VASP only looks at the first letter of this line (c/C/k/K = Cartesian, anything
    # else = Direct), so "direct", "D", "Cart", ... are all valid. Normalize to the
    # canonical 'Direct'/'Cartesian' the rest of this module compares against.
    typevectors = 'Cartesian' if datavasp[coord_line_idx][0][0].lower() in ('c', 'k') else 'Direct'
    getatoms = []
    for i in range(len(datavasp[5])):
        getatoms.append([i+1, element[datavasp[5][i]],datavasp[5][i], datavasp[6][i]])
    atomsposition = []
    cont = coord_line_idx + 1
    for el in getatoms:
        for i in range(int(el[3])):
            atomsposition.append([el[2], datavasp[cont][0], datavasp[cont][1], datavasp[cont][2]])
            cont = cont+1
    atomic_position = dic_atoms_position(atomsposition)
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


# This Function define the data for cif file
def getatomsandvectors_fhi(input_fhi):
    element, atomicnumber = periodic_table()
    datafhi = readfile(input_fhi)
    vectors = []
    atomdata = {}
    latticeparameter = "1.00"
    for lines in datafhi:
        if lines[0] == "atom_frac":
            atomdata.setdefault(lines[4], [])
            atomdata[lines[4]].append([
                f"{float(lines[1]):.8f}",
                f"{float(lines[2]):.8f}",
                f"{float(lines[3]):.8f}"])
            typevectors = 'Direct'
        elif lines[0] == "atom":
            atomdata.setdefault(lines[4], [])
            atomdata[lines[4]].append([
                f"{float(lines[1]):.8f}",
                f"{float(lines[2]):.8f}",
                f"{float(lines[3]):.8f}"])
            typevectors = 'Cartesian'
        elif lines[0] == 'lattice_vector':
            vectors.append([lines[1], lines[2], lines[3]])
    getatoms = []
    icont = 1
    for el in atomdata:
        getatoms.append([str(icont), str(element[el]),
                        el, str(len(atomdata[el]))])
        icont = icont+1
    atomic_position = atomdata

    return typevectors, latticeparameter, vectors, getatoms, atomic_position


# This Function define the data for cif file
# This Function define the data for cif file
def getatomsandvectors_cif(input_cif):
    """
    Extrai dados de estrutura de um arquivo .cif usando a biblioteca ASE.
    Retorna os dados no formato padrão esperado pelo script.
    """
    
    # 1. Obter o dicionário de elementos definido localmente no script
    element, atomicnumber = periodic_table()
    
    # 2. Ler o arquivo CIF usando a função 'ase_read'
    # Esta função foi importada no topo do script como 'ase_read'
    # Adicionamos um try...except para o caso da ASE não estar instalada.
    # NOTE: raises ValueError (not sys.exit) on every failure below -- caught by
    # convert_one's own try/except, so one bad CIF in a batch run is reported
    # and skipped instead of killing the whole batch.
    try:
        structure = ase_read(input_cif)
    except Exception as e:
        raise ValueError(f"could not read CIF: {e}") from e

    if len(structure) == 0:
        raise ValueError("parsed with 0 atoms -- check the file")

    occupancy = structure.info.get('occupancy')
    if occupancy:
        partial_sites = [(tag, occ) for tag, occ in occupancy.items()
                         if any(frac < 0.999 for frac in occ.values())]
        if partial_sites:
            details = "; ".join(
                f"atom {tag} ({', '.join(f'{sym}:{frac:g}' for sym, frac in occ.items())})"
                for tag, occ in partial_sites)
            raise ValueError(
                f"has partial-occupancy/disordered site(s): {details} -- stb-translate "
                "does not resolve disorder, pre-process the CIF (pick one occupant per "
                "site) before converting")

    # 3. Extrair vetores de rede (vectors)
    # O método .get_cell() da ASE retorna os vetores de rede completos,
    # não escalonados (em Angstroms).
    # O formato de retorno deve ser: [['v1x', 'v1y', 'v1z'], [...], [...]]
    vectors = []
    for vec in structure.get_cell():
        vectors.append([f"{vec[0]:.8f}", f"{vec[1]:.8f}", f"{vec[2]:.8f}"])

    # 4. Extrair posições atômicas (atomic_position)
    # .get_positions() retorna coordenadas Cartesianas em Angstroms.
    # O formato de retorno deve ser: {'Simbolo': [['x1', 'y1', 'z1'], ['x2', 'y2', 'z2'], ...]}
    
    symbols = structure.get_chemical_symbols() # Lista de símbolos, ex: ['C', 'C', 'H']
    positions = structure.get_positions()      # Array numpy de posições
    
    # Usamos 'dicatoms' como na função original para construir o dicionário
    dicatoms = {} 
    for sym, pos in zip(symbols, positions):
        if sym not in dicatoms:
            dicatoms[sym] = []
        # Adiciona a posição formatada como strings
        dicatoms[sym].append([f"{pos[0]:.8f}", f"{pos[1]:.8f}", f"{pos[2]:.8f}"])
    
    atomic_position = dicatoms

    # 5. Extrair informações das espécies (getatoms)
    # O formato de retorno deve ser: [['indice', 'num_atomico', 'simbolo', 'contagem'], ...]
    
    getatoms = []
    icont = 1 # O índice da espécie deve começar em 1
    
    # Iteramos sobre o dicionário de átomos que acabámos de criar
    for sym in dicatoms: 
        getatoms.append([
            f"{icont}",            # Indice (ex: '1')
            f"{element[sym]}",     # Num. Atômico (usando a função local)
            f"{sym}",              # Símbolo (ex: 'C')
            f"{len(dicatoms[sym])}" # Contagem (ex: '12')
        ])
        icont += 1 # Incrementa o índice para a próxima espécie

    # 6. Definir tipo de vetores e parâmetro de rede
    # Como a ASE retorna coordenadas Cartesianas e vetores de rede completos,
    # definimos o formato como Cartesiano e o parâmetro de rede como 1.0.
    typevectors = 'Cartesian'
    latticeparameter = '1.00'

    # 7. Retornar os dados no formato padrão esperado
    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_siesta(input_siesta):
    """Reads a SIESTA STRUCT_OUT-style file: 3 lattice-vector lines, an
    atom-count line, then one `species_idx Z x y z` line per atom.

    **Real, verified bug fixed here**: STRUCT_OUT identifies each atom by
    its SPECIES INDEX (`pos[0]`), not by element symbol -- a structure
    with fragment-labeled species that share an element (e.g. stb-adsorb's
    own `C_slab`/`C_ads`, distinct SIESTA species indices, same Z) has 2+
    DISTINCT species indices for the same Z. Grouping is done correctly by
    species index (`dicatoms`, keyed by `pos[0]`) for `getatoms`'s own
    counts -- but an earlier version of this function built `atomic_
    position` keyed by ELEMENT SYMBOL instead (`atomicnumber[position[1]]`),
    silently MERGING every same-Z species-index group's positions into one
    shared list. Downstream, `build_ase_atoms`/`writefilefdf` look up each
    `getatoms` row's position list by that same symbol key -- so every
    same-Z group ended up pulling the SAME (merged) position list, in
    full, once per group: confirmed live on a real 84-atom, 4-species
    (2 distinct C indices + 2 distinct O indices) STRUCT_OUT, this
    produced a 168-atom (73+73+11+11) `.fdf`, whose OWN malformed,
    duplicate-labeled `%block ChemicalSpeciesLabel` then compounded the
    error AGAIN on the post-write round-trip re-read (168 -> 336).

    Fixed by keeping every species-index group distinct all the way
    through: a group is labeled with its bare element symbol only if that
    symbol hasn't been used by an earlier (lower-index) group, otherwise
    `f"{symbol}_spec{species_idx}"` -- unique per group, never colliding,
    and honestly NOT claiming a semantic fragment name (STRUCT_OUT itself
    carries no label text, only the numeric species index, so `_slabN`/
    `_adsN`-style naming would be a fabricated guess). This label is used
    consistently as both `getatoms`'s own species-column AND `atomic_
    position`'s dict key, restoring the 1:1 correspondence every other
    getatomsandvectors_* function already relies on. A STRUCT_OUT with no
    same-Z species-index collisions (the common case) is completely
    unaffected: every group's symbol is already unique, so no suffix is
    ever added.
    """
    element, atomicnumber = periodic_table()
    datasiesta = readfile(input_siesta)
    latticeparameter = "1.00"
    typevectors = "Direct"
    vectors = datasiesta[:3]
    dicatoms = {}
    for pos in datasiesta[4:]:
        if pos[0] not in dicatoms:
            dicatoms[pos[0]] = []
        dicatoms[pos[0]].append([pos[0], pos[1], atomicnumber[pos[1]]])

    # One unique label per species index, in STRUCT_OUT's own index order --
    # bare symbol for a Z seen for the first time, suffixed for a repeat.
    seen_symbols = set()
    label_by_species_idx = {}
    for species_idx in dicatoms:
        symbol = dicatoms[species_idx][0][2]
        label = symbol if symbol not in seen_symbols else f"{symbol}_spec{species_idx}"
        seen_symbols.add(symbol)
        label_by_species_idx[species_idx] = label

    getatoms = []
    for species_idx in dicatoms:
        getatoms.append([dicatoms[species_idx][0][0],
                        dicatoms[species_idx][0][1],
                        label_by_species_idx[species_idx],
                        str(len(dicatoms[species_idx]))])
    atomic_position = AtomPositions()
    atomic_position.order = []
    for position in datasiesta[4:]:
        label = label_by_species_idx[position[0]]
        if label not in atomic_position:
            atomic_position[label] = []
        atomic_position.order.append((label, len(atomic_position[label])))
        atomic_position[label].append(
            [position[2], position[3], position[4]])
    return typevectors, latticeparameter, vectors, getatoms, atomic_position


def getatomsandvectors_dftb(input_dftb):
    element, atomicnumber = periodic_table()
    datadftb = readfile(input_dftb)
    vectors = []
    getatoms = []
    atomic_position = {}
    dic_data = {}
    latticeparameter = "1.00"
    if datadftb[0][1] == 'S':
        typevectors = 'Cartesian'
    elif datadftb[0][1] == 'F':
        typevectors = 'Direct'
    else:
        # Raises (not sys.exit/exit()) so convert_one's try/except catches it --
        # one bad file in a batch run is reported and skipped, not fatal to the
        # whole batch.
        raise ValueError("coordinate type not defined on line 1 -- only 'S' or 'F' accepted")
    # the last three lines are a, b, c in that order (the origin line comes before them); reading them from the
    # end swapped a and c, and a non-cubic .gen came back as another structure
    for i in (3, 2, 1):
        vectors.append([f"{float(datadftb[-i][0]):.8f}",
                        f"{float(datadftb[-i][1]):.8f}",
                        f"{float(datadftb[-i][2]):.8f}"])
    icont = 0
    dic_data = {f"{i + 1}": elem for i, elem in enumerate(datadftb[1])}
    for line in datadftb[2:-4]:
        atomic_position.setdefault(dic_data[line[1]], [])
        atomic_position[dic_data[line[1]]].append(
            [f"{float(line[2]):.8f}", f"{float(line[3]):.8f}", f"{float(line[4]):.8f}"])
    icont = 1
    for el in atomic_position:
        getatoms.append([f"{icont}",
                         f"{element[el]}",
                         f"{el}",
                         f"{len(atomic_position[el])}"])
        icont = icont+1

    return typevectors, latticeparameter, vectors, getatoms, atomic_position


def getatomsandvectors_xsf(input_xsf):
    print("[WARNING] Only for PRIMVEC format")
    element, atomicnumber = periodic_table()
    dataxsf = readfile(input_xsf)
    vectors = []
    getatoms = []
    atomic_position = {}
    dic_data = {}
    latticeparameter = "1.00"
    typevectors = 'Cartesian'
    dataxsf = [v for v in dataxsf if not str(v[0]).startswith("#")]
    for j in range(len(dataxsf)):
        if dataxsf[j][0] == 'PRIMVEC':
            for i in range(1, 4):
                vectors.append([f"{float(dataxsf[j+i][0]):.8f}",
                                f"{float(dataxsf[j+i][1]):.8f}",
                                f"{float(dataxsf[j+i][2]):.8f}"])
    for j in range(len(dataxsf)):
        if dataxsf[j][0] == 'PRIMCOORD':
            na = dataxsf[j+1][0]
            for i in range(int(na)):
                atomic_position.setdefault(atomicnumber[dataxsf[j+i+2][0]], [])
                atomic_position[atomicnumber[dataxsf[j+i+2][0]]].append(
                    [f"{float(dataxsf[j+i+2][1]):.8f}",
                     f"{float(dataxsf[j+i+2][2]):.8f}",
                     f"{float(dataxsf[j+i+2][3]):.8f}"])
    icont = 1
    for el in atomic_position:
        getatoms.append([f"{icont}",
                         f"{element[el]}",
                         f"{el}",
                         f"{len(atomic_position[el])}"])
        icont = icont+1

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_fdf(input_fdf):
    """
    Extrai dados de estrutura de um arquivo .fdf do SIESTA.

    Delegates parsing to stb.core.structure_io.read_fdf and adapts the result
    to translate.py's own canonical 5-tuple format.
    """
    structure = structure_io.read_fdf(input_fdf)

    typevectors = 'Direct' if structure.coord_format == 'fractional' else 'Cartesian'
    latticeparameter = str(structure.lattice_constant)
    # Raw (unscaled) vectors: writefilefdf writes latticeparameter as its own
    # line, so the vectors themselves must not be pre-multiplied by it.
    vectors = [[f"{v:.8f}" for v in row] for row in structure_io.raw_lattice_vectors(structure)]

    atomic_position = AtomPositions()
    atomic_position.order = []
    for symbol, pos in structure.atoms:
        group = atomic_position.setdefault(symbol, [])
        atomic_position.order.append((symbol, len(group)))
        group.append([f"{pos[0]:.8f}", f"{pos[1]:.8f}", f"{pos[2]:.8f}"])

    getatoms = []
    for symbol in structure.species:
        meta = structure.species_meta[symbol]
        count = len(atomic_position.get(symbol, []))
        if count > 0:
            getatoms.append([meta['id'], str(meta['Z']), symbol, str(count)])

    return typevectors, latticeparameter, vectors, getatoms, atomic_position


###################### Write functions ################


def writefilefdf(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    """Writes a .fdf file, delegating serialization to stb.core.structure_io.write_fdf.

    translate.py keeps owning the coordinate conversion (convert_coordinates)
    since it's shared by every output format this tool writes, not just .fdf.
    """
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )

    species = [elem[2] for elem in getatoms]
    species_meta = {elem[2]: {'id': elem[0], 'Z': elem[1]} for elem in getatoms}
    order = getattr(final_positions, "order", None)
    if order is None or sorted(order) != sorted((e[2], i) for e in getatoms for i in range(len(final_positions[e[2]]))):
        order = [(e[2], i) for e in getatoms for i in range(len(final_positions[e[2]]))]   # grouped by species
    atoms = [
        (label, np.array([float(c) for c in final_positions[label][i][:3]]))
        for label, i in order
    ]
    lattice_constant = float(latticeparameter)
    raw_lattice = np.array(vectors, dtype=float)

    structure = structure_io.FdfStructure(
        lattice=raw_lattice * lattice_constant,
        lattice_constant=lattice_constant,
        species=species,
        species_meta=species_meta,
        atoms=atoms,
        coord_format='fractional' if final_type == 'Direct' else 'cartesian',
    )
    structure_io.write_fdf(structure, outfilename)


def writefileposcar(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    # --- NEW: Convert coordinates ---
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )
    # --- END NEW ---

    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    outfile.append(f"{latticeparameter}")
    for lin in vectors:
        outfile.append(f"{lin[0]}   {lin[1]}   {lin[2]} ")
    lineatoms = ""
    linenatoms = ""
    for atoms in getatoms:
        lineatoms = lineatoms + f"{atoms[2]}   "
        linenatoms = linenatoms+f"{atoms[3]}   "
    outfile.append(f"{lineatoms}")
    outfile.append(f"{linenatoms}")
    
    # --- MODIFIED: Use final_type ---
    if final_type == 'Direct':
        outfile.append("Direct")
    if final_type == 'Cartesian':
        outfile.append("Cartesian")
    # --- END MODIFIED ---

    # --- MODIFIED: Use final_positions ---
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"{position[0]}   {position[1]}   {position[2]}")
    # --- END MODIFIED ---

    np.savetxt(outfilename, outfile, fmt='%s')
    return

def writefilecif(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    """
    Writes a CIF file via pymatgen's CifWriter, which detects and writes the
    real space group/symmetry operations instead of always forcing P1.
    Note: CIF *always* uses fractional (Direct) coordinates.
    """
    from pymatgen.io.cif import CifWriter

    if coord_format and coord_format.lower() == 'cartesian':
        print("[WARNING] CIF format requires fractional (Direct) coordinates.")
        print("[INFO]    Ignoring '--coord-format cartesian' and converting to Direct.")

    # Force conversion to Direct, regardless of user input
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'direct'  # Force 'direct'
    )

    atoms = build_ase_atoms(final_type, latticeparameter, vectors, getatoms, final_positions)
    structure = AseAtomsAdaptor.get_structure(atoms)
    # CIF_SYMPREC (--cif-symprec): with a tolerance, CifWriter finds the space group, moves every atom onto its exact
    # symmetric position (by up to the tolerance) and writes the asymmetric unit in the standard setting; None writes
    # P1 with the positions exactly as read
    CifWriter(structure, symprec=CIF_SYMPREC).write_file(outfilename)




def writefilexyz(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    # --- NEW: Conversion logic for XYZ ---
    if coord_format and coord_format.lower() == 'direct':
        print("[WARNING] XYZ format requires Cartesian (Angstrom) coordinates.")
        print("[INFO]    Ignoring '--coord-format direct' and converting to Cartesian.")

    # Force conversion to Cartesian
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'cartesian' # Force 'cartesian'
    )
    # --- END NEW ---
    # NOTE: no leading "# automatic create..." attribution line here, unlike
    # every other writer -- XYZ's format is a strict 2-line header (atom
    # count, then a single comment line) with no room for a 3rd line; a stray
    # extra line here silently corrupts every stb-translate-written .xyz
    # (getatomsandvectors_xyz's fixed `[2:]` header skip would then treat the
    # real comment line as if it were the first atom) -- caught live by the
    # --dry-run/write round-trip check (verify_round_trip) reporting a
    # KeyError on the comment text when re-reading a just-written .xyz.
    outfile = []
    natoms_total = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        natoms_total = natoms_total+int(el[3])
    outfile.append(f"{natoms_total}")
    outfile.append(comment)

    # --- MODIFIED: Simplified to use final_positions ---
    # Since we know final_positions is 'Cartesian', the logic is simple
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"{elem[2]}   {position[0]}   {position[1]}   {position[2]}")
    # --- END MODIFIED ---

    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefiledftb(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    # --- NEW: Convert coordinates ---
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )
    # --- END NEW ---
    # NOTE: no leading "# automatic create..." attribution line here -- like
    # XYZ, the DFTB+ .gen format's first line is fixed-position ("N S"/"N F"),
    # and getatomsandvectors_dftb reads it as datadftb[0][1]; an extra line
    # above it silently corrupts every stb-translate-written .gen the same
    # way it did for XYZ (caught live by the write round-trip check).
    outfile = []
    natoms_total = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        natoms_total = natoms_total+int(el[3])

    # --- MODIFIED: Use final_type ---
    if final_type == 'Cartesian':
        outfile.append(f"{natoms_total}   S") # S = Cartesian
    elif final_type == 'Direct':
        outfile.append(f"{natoms_total}   F") # F = Fractional
    # --- END MODIFIED ---

    atoms = ""
    for i in range(len(getatoms)):
        atoms = atoms + f"{getatoms[i][2]}   "
    outfile.append(atoms)

    icont = 1
    # --- MODIFIED: Use final_positions ---
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"    {icont}  {elem[0]}  {position[0]}   {position[1]}   {position[2]}")
            icont = icont+1
    # --- END MODIFIED ---

    outfile.append(f"    0.00000000  0.00000000 0.00000000")
    vectors_np = np.array(vectors, dtype="float")*float(latticeparameter)
    for i in range(3):
        outfile.append(f"    {vectors_np[i][0]:.8f}   {vectors_np[i][1]:.8f}   {vectors_np[i][2]:.8f}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefilexsf(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    # --- NEW: Conversion logic for XSF ---
    if coord_format and coord_format.lower() == 'direct':
        print("[WARNING] XSF format (PRIMCOORD) requires Cartesian (Angstrom) coordinates.")
        print("[INFO]    Ignoring '--coord-format direct' and converting to Cartesian.")

    # Force conversion to Cartesian
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'cartesian' # Force 'cartesian'
    )
    # --- END NEW ---
    
    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    natoms_total = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        natoms_total = natoms_total+int(el[3])
    vectors_np = np.array(vectors, dtype="float")*float(latticeparameter)
    outfile.append('# create by stb-translate ')
    outfile.append(f'# {comment}\n')
    outfile.append('CRYSTAL')
    outfile.append(f"PRIMVEC")
    for i in range(3):
        outfile.append(f"    {vectors_np[i][0]:.8f}   {vectors_np[i][1]:.8f}   {vectors_np[i][2]:.8f}")
    
    # --- MODIFIED: Simplified to use final_positions ---
    # Since we know final_positions is 'Cartesian':
    outfile.append(f"PRIMCOORD")
    outfile.append(f"{natoms_total}  1")
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"    {elem[1]}   {position[0]}   {position[1]}   {position[2]}")
    # --- END MODIFIED ---
    
    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefilefhi(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    # --- NEW: Convert coordinates ---
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )
    # --- END NEW ---

    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    natoms_total = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        natoms_total = natoms_total+int(el[3])
    vectors_np = np.array(vectors, dtype="float")*float(latticeparameter)
    for vector in vectors_np:
        outfile.append(f"lattice_vector   {vector[0]:.8f}   {vector[1]:.8f}   {vector[2]:.8f}")
    outfile.append(" ")
    
    # --- MODIFIED: Use final_type and final_positions ---
    if final_type == 'Cartesian':
        for elem in getatoms:
            for position in final_positions[elem[2]]:
                outfile.append(
                    f"atom   {position[0]}   {position[1]}   {position[2]}   {elem[2]}")
    elif final_type == 'Direct':
        for elem in getatoms:
            for position in final_positions[elem[2]]:
                outfile.append(
                    f"atom_frac   {position[0]}   {position[1]}   {position[2]}   {elem[2]}")
    # --- END MODIFIED ---

    np.savetxt(outfilename, outfile, fmt='%s')
    return


def writefilelammps(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    """Writes a LAMMPS data file (atom_style 'atomic', 'metal' units --
    Angstrom/eV, matching a DFT-derived structure) via ASE's own writer, for
    handing a relaxed SIESTA structure off to classical MD/force-field
    simulations. Reuses build_ase_atoms (the same conversion this file's
    --view/--view-image use) rather than hand-rolling LAMMPS's triclinic
    box-bounds convention (tilt factors for a non-orthogonal cell need a
    basis rotation ASE's Prism already implements and tests)."""
    if coord_format and coord_format.lower() == 'direct':
        print("[WARNING] LAMMPS data format requires Cartesian (Angstrom) coordinates.")
        print("[INFO]    Ignoring '--coord-format direct' and converting to Cartesian.")

    atoms = build_ase_atoms(typevectors, latticeparameter, vectors, getatoms, atomsposition)
    from ase.io import write as ase_write
    ase_write(outfilename, atoms, format='lammps-data', masses=True)


##### NEW HELPER FUNCTION #####
class AtomPositions(dict):
    """{label: [[x, y, z], ...]} (every format's grouped-by-species layout) that also remembers the
    file's own atom order as `order`, a list of (label, index within that label). Formats that
    group atoms by species (POSCAR, DFTB+...) ignore it; .fdf output uses it so the atom order of
    an .fdf / STRUCT_OUT input survives (fragment-size and atom-index based tools depend on it)."""
    order = None


def _carry_order(src, dst):
    """dst (a plain dict from a coordinate conversion) gets src's atom order, if it has one"""
    order = getattr(src, "order", None)
    if order is None or dst is src:
        return dst
    out = AtomPositions(dst)
    out.order = order
    return out


def convert_coordinates(typevectors_in, latticeparameter, vectors, atomsposition_in, coord_format_out):
    final_type, final_positions = _convert_coordinates(
        typevectors_in, latticeparameter, vectors, atomsposition_in, coord_format_out)
    return final_type, _carry_order(atomsposition_in, final_positions)


def _convert_coordinates(typevectors_in: str,
                        latticeparameter: str, 
                        vectors: List[List[str]], 
                        atomsposition_in: Dict[str, List[List[str]]], 
                        coord_format_out: str) -> (str, Dict[str, List[List[str]]]):
    """
    Converts coordinates between Direct (Fractional) and Cartesian (Angstrom)
    based on the desired output format.

    Args:
        typevectors_in (str): The input format ('Direct' or 'Cartesian').
        latticeparameter (str): The lattice parameter (scaling factor).
        vectors (List[List[str]]): The lattice vectors.
        atomsposition_in (Dict): The input positions dictionary.
        coord_format_out (str): The desired output format ('direct' or 'cartesian', or None).

    Returns:
        Tuple[str, Dict]: (final_type, final_positions)
                         The final format and the converted positions dictionary.
    """
    
    # 1. Prepare the lattice matrix, scaled by the lattice parameter
    lattice_matrix = np.array(vectors, dtype=float) * float(latticeparameter)
    
    # 2. Determine input and output types
    # If the user did not specify a format, use the input format
    if coord_format_out is None:
        return typevectors_in, atomsposition_in

    # Normalize names to 'Direct' and 'Cartesian' for comparison
    type_in_norm = 'Direct' if typevectors_in.lower() == 'direct' else 'Cartesian'
    type_out_norm = 'Direct' if coord_format_out.lower() == 'direct' else 'Cartesian'

    # 3. Check if conversion is needed
    if type_in_norm == type_out_norm:
        # No conversion needed. Return original data, but with normalized type.
        return type_out_norm, atomsposition_in

    # 4. Perform conversion
    atomsposition_out = {}
    
    if type_in_norm == 'Direct' and type_out_norm == 'Cartesian':
        # --- Convert from Direct -> Cartesian ---
        for symbol, positions in atomsposition_in.items():
            atomsposition_out[symbol] = []
            for pos in positions:
                pos_np = np.array(pos, dtype=float)
                # v_cart = v_direct . M  (lattice_matrix rows are the lattice vectors)
                cart_pos = np.dot(pos_np, lattice_matrix)
                # Format back to list of strings
                atomsposition_out[symbol].append(
                    [f"{cart_pos[0]:.8f}", f"{cart_pos[1]:.8f}", f"{cart_pos[2]:.8f}"]
                )
        return 'Cartesian', atomsposition_out

    elif type_in_norm == 'Cartesian' and type_out_norm == 'Direct':
        # --- Convert from Cartesian -> Direct ---
        # v_direct = v_cart . M_inv
        inv_lattice = np.linalg.inv(lattice_matrix)
        for symbol, positions in atomsposition_in.items():
            atomsposition_out[symbol] = []
            for pos in positions:
                pos_np = np.array(pos, dtype=float)
                direct_pos = np.dot(pos_np, inv_lattice)
                # Format back to list of strings
                atomsposition_out[symbol].append(
                    [f"{direct_pos[0]:.8f}", f"{direct_pos[1]:.8f}", f"{direct_pos[2]:.8f}"]
                )
        return 'Direct', atomsposition_out
    
    # Fallback (should not happen)
    return typevectors_in, atomsposition_in
##### END NEW HELPER FUNCTION #####


##### 3D VIEWER HELPERS (--view / --view-image) #####
def build_ase_atoms(typevectors, latticeparameter, vectors, getatoms, atomsposition):
    """Builds an ase.Atoms from the 5-tuple every getatomsandvectors_* returns,
    for the optional --view/--view-image 3D visualization (also used for
    the pre-write structure summary/close-contact check and the post-write
    round-trip verification -- see convert_one/reread_for_check/
    verify_round_trip).

    **Real, verified bug fixed here**: `elem[2]` (a getatoms row's own
    species-column string) is not always a real element symbol -- any
    reader can legitimately return a compound, suffixed label there (this
    suite's own `_slab`/`_ads`/`_ghost` fragment-label convention for an
    fdf input, or getatomsandvectors_siesta's own `_specN` disambiguation
    for a STRUCT_OUT input with same-Z species groups, see that function's
    docstring). `ase.Atoms(symbols=[...])` looks up each string in ASE's
    own periodic table and raises `KeyError` on anything that isn't a bare
    element symbol -- confirmed live: `Atoms(symbols=['C_slab'], ...)`
    raises `KeyError: 'C_slab'` outright. `elem[1]` (every getatoms row's
    OWN atomic-number column, int or str depending on the reader --
    normalized via str() here) is unaffected by any such label and always
    resolves to the true element, so it -- not elem[2] -- is what should
    drive the ASE-facing chemical symbol.
    """
    from ase import Atoms
    _, atomicnumber = periodic_table()

    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'cartesian'
    )
    cell = np.array(vectors, dtype=float) * float(latticeparameter)

    symbols = []
    positions = []
    for elem in getatoms:
        label, real_symbol = elem[2], atomicnumber[str(elem[1])]
        for pos in final_positions[label]:
            symbols.append(real_symbol)
            positions.append([float(pos[0]), float(pos[1]), float(pos[2])])

    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)



def save_structure_image(atoms, path):
    """Renders the structure (atoms + unit cell edges) in 3D with matplotlib and
    saves it to `path`. Works headless (no display needed)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    cell = atoms.get_cell()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for species in sorted(set(symbols)):
        idx = [i for i, s in enumerate(symbols) if s == species]
        pts = positions[idx]
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], label=species, s=80, depthshade=True)

    if cell is not None and np.linalg.det(cell) != 0:
        import itertools
        coeffs_list = list(itertools.product([0, 1], repeat=3))
        corners = np.array([sum(c * v for c, v in zip(coeffs, cell)) for coeffs in coeffs_list])
        for i in range(8):
            for j in range(i + 1, 8):
                if sum(a != b for a, b in zip(coeffs_list[i], coeffs_list[j])) == 1:
                    p1, p2 = corners[i], corners[j]
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                            color='black', linewidth=0.8)

    ax.set_xlabel('X (Ang)')
    ax.set_ylabel('Y (Ang)')
    ax.set_zlabel('Z (Ang)')
    ax.set_box_aspect((1, 1, 1))
    ax.legend()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
##### END 3D VIEWER HELPERS #####


##### ANALYSIS, CHECKS AND REPORT #####
CLOSE_CONTACT = 0.5          # Ang: two atoms closer than this are almost certainly a unit/coordinate mix-up
SHORT_BOND_FACTOR = 0.75     # x the sum of covalent radii: a bond this short deserves a look
SYMPRECS = (1e-3, 1e-2, 1e-1)


def guess_output_format(path):
    """--out-format guess from --out-file's name; None if not recognized."""
    base = os.path.basename(path).lower()
    if base in BASENAME_FORMAT_MAP:
        return BASENAME_FORMAT_MAP[base]
    _, ext = os.path.splitext(base)
    return OUTPUT_EXTENSION_GUESS.get(ext)


def _how(given, guessed_from):
    return "given (-if/-of)" if given else f"auto, from the name '{guessed_from}'"


def _file_info(path):
    try:
        size = os.path.getsize(path)
        with open(path, errors="replace") as f:
            lines = sum(1 for _ in f)
    except OSError:
        return "-"
    unit = f"{size} bytes" if size < 10240 else f"{size / 1024:.1f} KiB"
    return f"{unit}, {lines} lines"


def read_structure(fmt, path, lattice_vectors=None):
    """The 5-tuple of the matching getatomsandvectors_* reader."""
    if fmt == "xyz":
        v = lattice_vectors
        vec = [[str(v[0]), str(v[1]), str(v[2])], [str(v[3]), str(v[4]), str(v[5])],
               [str(v[6]), str(v[7]), str(v[8])]]
        return getatomsandvectors_xyz(path, vec)
    readers = {"poscar": getatomsandvectors_vasp, "cif": getatomsandvectors_cif, "siesta": getatomsandvectors_siesta,
               "fhi": getatomsandvectors_fhi, "dftb": getatomsandvectors_dftb, "xsf": getatomsandvectors_xsf,
               "fdf": getatomsandvectors_fdf}
    return readers[fmt](path)


WRITERS = {"xyz": "writefilexyz", "poscar": "writefileposcar", "fdf": "writefilefdf", "dftb": "writefiledftb",
           "xsf": "writefilexsf", "fhi": "writefilefhi", "cif": "writefilecif", "lammps": "writefilelammps"}


def reread_for_check(out_format, out_file, vectors, latticeparameter):
    """Re-reads the written file with the matching reader: (ase.Atoms, coordinate type as stored). An xyz output
    carries no cell, so the input's is reused; lammps is output-only here and goes through ASE."""
    if out_format == "lammps":
        from ase.io import read as ase_read_atoms
        return ase_read_atoms(out_file, format='lammps-data'), "cartesian"
    t = read_structure(out_format, out_file) if out_format != "xyz" else getatomsandvectors_xyz(out_file, vectors)
    return build_ase_atoms(*t), t[0]


def written_coordinates(fmt, path):
    """How the written file stores the positions, read from the file itself."""
    fixed = {"cif": "fractional (CIF always)", "xyz": "Cartesian, Ang (XYZ always)",
             "siesta": "fractional (STRUCT_OUT always)",
             "xsf": "Cartesian, Ang (XSF always)", "lammps": "Cartesian, Ang (LAMMPS data)"}
    if fmt in fixed:
        return fixed[fmt]
    try:
        lines = [l.strip() for l in open(path, errors="replace")]
    except OSError:
        return "-"
    if fmt == "poscar":                       # line 8 (or 9 after 'Selective dynamics') holds the mode
        k = 8 if len(lines) > 7 and lines[7][:1].lower() == "s" else 7
        return "fractional (Direct)" if lines[k][:1].lower() == "d" else "Cartesian"
    if fmt == "fdf":
        for l in lines:
            if l.lower().startswith("atomiccoordinatesformat"):
                v = l.split()[1] if len(l.split()) > 1 else "?"
                return ("fractional" if v.lower().startswith("frac") else "Cartesian") + f" (AtomicCoordinatesFormat {v})"
    if fmt == "dftb":
        tok = lines[0].split()[1] if lines and len(lines[0].split()) > 1 else "?"
        return {"F": "fractional (F)", "S": "Cartesian, Ang (S)", "C": "Cartesian, cluster (C)"}.get(tok.upper(), tok)
    if fmt == "fhi":
        return "fractional (atom_frac)" if any(l.startswith("atom_frac") for l in lines) else "Cartesian (atom)"
    return "-"


def structure_facts(atoms, getatoms, vacuum_gap=10.0):
    """Everything the report says about one structure (no printing)."""
    from ase.data import atomic_masses, covalent_radii
    from stb.core.kspace import detect_vacuum_axes
    cell = np.array(atoms.cell)
    n = len(atoms)
    lengths, angles = atoms.cell.lengths(), atoms.cell.angles()
    volume = abs(np.linalg.det(cell))
    mass = float(sum(atomic_masses[z] for z in atoms.numbers))
    try:
        vac = detect_vacuum_axes(atoms.get_scaled_positions(wrap=True), cell, vacuum_gap)
    except Exception:
        vac = [False] * 3
    species = [(str(elem[2]), int(elem[1]), int(elem[3])) for elem in getatoms]
    facts = {"n": n, "lengths": lengths, "angles": angles, "volume": volume, "vol_atom": volume / max(n, 1),
             "density": mass / volume * 1.66053906660 if volume > 0 else None, "mass": mass, "vacuum": vac,
             "species": [(lab, z, cnt, cnt / n, float(atomic_masses[z])) for lab, z, cnt in species]}
    # nearest neighbour per species pair (minimum image), shortest bonds flagged
    pairs = {}
    if n > 1:
        d = atoms.get_all_distances(mic=True)
        np.fill_diagonal(d, np.inf)
        sym = atoms.get_chemical_symbols()
        for i in range(n):
            for j in range(i + 1, n):
                key = tuple(sorted((sym[i], sym[j])))
                if d[i, j] < pairs.get(key, (np.inf,))[0]:
                    pairs[key] = (float(d[i, j]), i, j)
        facts["dmin"] = float(d.min())
    else:
        facts["dmin"] = None
    facts["pairs"] = [(a, b, v[0], v[1], v[2],
                       v[0] / (covalent_radii[atoms.numbers[v[1]]] + covalent_radii[atoms.numbers[v[2]]]))
                      for (a, b), v in sorted(pairs.items(), key=lambda kv: kv[1][0])]
    return facts


def symmetry_facts(atoms):
    """Space group at several tolerances, and the primitive cell at 0.01 Ang (spglib)."""
    import spglib
    cell = (np.array(atoms.cell), atoms.get_scaled_positions(wrap=True), atoms.numbers)
    out = []
    for sp in SYMPRECS:
        try:
            ds = spglib.get_symmetry_dataset(cell, symprec=sp)
        except Exception:
            ds = None
        if ds is None:
            out.append((sp, None, None, None))
            continue
        g = getattr(ds, "international", None) or ds["international"]
        num = getattr(ds, "number", None) or ds["number"]
        rot = getattr(ds, "rotations", None)
        rot = rot if rot is not None else ds["rotations"]
        out.append((sp, g, num, len(rot)))
    prim = None
    try:
        p = spglib.standardize_cell(cell, to_primitive=True, symprec=1e-2)
        prim = len(p[2]) if p is not None else None
    except Exception:
        pass
    system, pg = None, None
    try:
        system, pg = crystal_system(AseAtomsAdaptor.get_structure(atoms))
    except Exception:
        pass
    return {"table": out, "primitive": prim, "system": system, "point_group": pg}


def compare_structures(a_in, a_out):
    """Round-trip comparison of the input and the re-read output."""
    from pymatgen.core import Composition
    c_in, c_out = Composition(a_in.get_chemical_formula()), Composition(a_out.get_chemical_formula())
    _, z_in = c_in.get_reduced_composition_and_factor()
    _, z_out = c_out.get_reduced_composition_and_factor()
    r = {"formula_in": c_in.reduced_formula, "formula_out": c_out.reduced_formula,
         "n_in": len(a_in), "n_out": len(a_out),
         "vfu_in": abs(np.linalg.det(a_in.cell)) / z_in, "vfu_out": abs(np.linalg.det(a_out.cell)) / z_out,
         "cell_diff": None, "pos_diff": None, "note": None}
    same_cell = np.allclose(np.array(a_in.cell), np.array(a_out.cell), atol=1e-4)
    r["cell_diff"] = float(np.abs(np.array(a_in.cell) - np.array(a_out.cell)).max())
    if len(a_in) == len(a_out) and same_cell:
        from ase.geometry import get_distances
        _, dm = get_distances(a_in.positions, a_out.positions, cell=a_in.cell, pbc=a_in.pbc)
        sym_in, sym_out = np.array(a_in.get_chemical_symbols()), np.array(a_out.get_chemical_symbols())
        dm[sym_in[:, None] != sym_out[None, :]] = np.inf
        r["pos_diff"] = float(dm.min(axis=1).max())         # each input atom to its nearest same-element output atom
    elif not same_cell:
        r["note"] = ("the writer re-expressed the cell (the standard setting of a symmetrized CIF): compared as the same "
                     "crystal up to a change of cell/origin (pymatgen StructureMatcher)")
        try:
            from pymatgen.analysis.structure_matcher import StructureMatcher
            si, so = AseAtomsAdaptor.get_structure(a_in), AseAtomsAdaptor.get_structure(a_out)
            sm = StructureMatcher(ltol=0.02, stol=0.05, angle_tol=1, primitive_cell=True, scale=False)
            r["match"] = bool(sm.fit(si, so))
            rms = sm.get_rms_dist(si, so) if r["match"] else None
            r["match_rms"] = float(rms[0]) * (abs(np.linalg.det(a_in.cell)) / len(a_in)) ** (1 / 3) if rms else None
        except Exception:
            r["match"] = None
    return r


def convert_one(args, in_file, in_format, out_file, out_format):
    """Reads, analyses and (unless --dry-run) writes one file; returns a dict for the report (never raises)."""
    rec = {"in_file": in_file, "in_format": in_format, "out_file": out_file, "out_format": out_format,
           "ok": False, "error": None, "atoms": None}
    try:
        t = read_structure(in_format, in_file, args.lattice_vectors)
        typevectors, latticeparameter, vectors, getatoms, atomsposition = t
        atoms = build_ase_atoms(*t)
    except FileNotFoundError as e:
        rec["error"] = f"file not found: {e}"
        return rec
    except Exception as e:
        rec["error"] = (f"could not parse '{in_file}' as '{in_format}': {e} -- check that the file matches the expected "
                        "structure for this format (-if)")
        return rec
    vals = np.array([float(c) for pos in atomsposition.values() for p in pos for c in p[:3]] or [0.0])
    rec["coord_range"] = (float(vals.min()), float(vals.max()))
    rec.update(atoms=atoms, typevectors=typevectors, latticeparameter=latticeparameter, getatoms=getatoms,
               in_size=_file_info(in_file), facts=structure_facts(atoms, getatoms), sym=symmetry_facts(atoms))
    if args.dry_run:
        rec["ok"] = True
        return rec
    try:
        globals()[WRITERS[out_format]](typevectors, latticeparameter, vectors, getatoms, atomsposition, out_file,
                                       args.coord_format)
    except Exception as e:
        rec["error"] = f"could not write '{out_file}' as '{out_format}': {e}"
        return rec
    rec["ok"] = True
    rec["out_size"] = _file_info(out_file)
    rec["coords_out"] = written_coordinates(out_format, out_file)
    try:
        atoms_out, _ = reread_for_check(out_format, out_file, vectors, latticeparameter)
        rec["check"] = compare_structures(atoms, atoms_out)
    except Exception as e:
        rec["check_error"] = str(e)
    return rec


def report_file(rec, sections, f_out, warn, tag=""):
    """The report sections of one converted file."""
    t = f" -- {os.path.basename(rec['in_file'])}" if tag else ""
    if not rec["ok"] or rec.get("facts") is None:
        sections("FAILED" + t)
        print_dual(color_text(f"[FAIL] {rec['error']}", 'red'), f_out)
        warn(f"{rec['in_file']}: {rec['error']}")
        return
    f, s, atoms = rec["facts"], rec["sym"], rec["atoms"]
    from pymatgen.core import Composition
    comp = Composition(atoms.get_chemical_formula())
    reduced, zfu = comp.get_reduced_formula_and_factor()

    sections("INPUT" + t)
    print_dual(f"File              : {rec['in_file']} ({rec['in_size']})", f_out)
    print_dual(f"Format            : {rec['in_format']} ({rec['in_how']})", f_out)
    print_dual(f"Coordinates       : {written_coordinates(rec['in_format'], rec['in_file'])}; lattice parameter "
               f"{rec['latticeparameter']} Ang", f_out)

    lo, hi = rec["coord_range"]
    shortest_axis = float(min(atoms.cell.lengths()))
    if str(rec["typevectors"]).lower().startswith("d") and (hi > 1.25 or lo < -0.25):
        print_dual(color_text(f"Coordinate range  : {lo:.3f} .. {hi:.3f} -- fractional coordinates far outside [0, 1]",
                              'yellow'), f_out)
        warn(f"{os.path.basename(rec['in_file'])}: the coordinates are read as fractional but span {lo:.2f} .. "
             f"{hi:.2f} -- a Cartesian file labelled Direct/Fractional? The structure below is then wrong.")
    elif str(rec["typevectors"]).lower().startswith("c") and hi <= 1.0 and lo >= 0.0 and shortest_axis > 2.0 \
            and len(atoms) > 1:
        print_dual(color_text(f"Coordinate range  : {lo:.3f} .. {hi:.3f} Ang -- every Cartesian coordinate inside "
                              "[0, 1] Ang", 'yellow'), f_out)
        warn(f"{os.path.basename(rec['in_file'])}: every Cartesian coordinate is between 0 and 1 Ang in a cell of "
             f"{shortest_axis:.1f} Ang or more -- fractional coordinates labelled Cartesian?")
    else:
        unit = "" if str(rec["typevectors"]).lower().startswith("d") else " Ang"
        print_dual(f"Coordinate range  : {lo:.4f} .. {hi:.4f}{unit} (as read)", f_out)

    sections("STRUCTURE" + t)
    print_dual(f"Formula           : {atoms.get_chemical_formula()} -- reduced {reduced}, Z = {zfu:g} formula units", f_out)
    print_dual(f"Atoms             : {f['n']} ({len(f['species'])} species)", f_out)
    rows = [([lab, chemical_symbol(z), str(z), str(cnt), f"{100 * frac:.1f}", f"{m:.3f}"], None)
            for lab, z, cnt, frac, m in f["species"]]
    print_table(["label", "element", "Z", "count", "fraction (%)", "mass (u)"], rows, f_out)
    a, b, c = f["lengths"]
    al, be, ga = f["angles"]
    print_dual(f"Cell              : a = {a:.5f}, b = {b:.5f}, c = {c:.5f} Ang; alpha = {al:.3f}, beta = {be:.3f}, "
               f"gamma = {ga:.3f} deg", f_out)
    print_dual("Lattice vectors (Ang):", f_out)
    for k, row in zip("abc", np.array(atoms.cell)):
        print_dual(f"    {k}: {row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}", f_out)
    dens = f"{f['density']:.4f} g/cm^3" if f["density"] else "-"
    print_dual(f"Volume            : {f['volume']:.4f} Ang^3 ({f['vol_atom']:.4f} Ang^3/atom); density {dens} "
               f"(mass {f['mass']:.3f} u)", f_out)
    nvac = sum(f["vacuum"])
    kind = {0: "3D bulk (no vacuum gap)", 1: "2D slab/layer", 2: "1D wire", 3: "0D molecule/cluster in a box"}[nvac]
    axes = ", ".join("abc"[i] for i in range(3) if f["vacuum"][i])
    print_dual(f"Periodicity       : {kind}" + (f" -- vacuum (>= 10 Ang gap) along {axes}" if axes else ""), f_out)
    if f["pairs"]:
        rows = []
        for e1, e2, d, i, j, ratio in f["pairs"]:
            flag = "TOO CLOSE" if d < CLOSE_CONTACT else ("short" if ratio < SHORT_BOND_FACTOR else "")
            rows.append(([f"{e1}-{e2}", f"{d:.4f}", f"{e1}{i + 1} - {e2}{j + 1}", f"{ratio:.2f}", flag],
                         "yellow" if flag else None))
        print_dual("Shortest distance per element pair (minimum image):", f_out)
        print_table(["pair", "distance (Ang)", "atoms", "/ covalent", "check"], rows, f_out)
        if f["dmin"] is not None and f["dmin"] < CLOSE_CONTACT:
            warn(f"{os.path.basename(rec['in_file'])}: two atoms are {f['dmin']:.3f} Ang apart -- almost certainly a "
                 "Direct/Cartesian mix-up or a wrong lattice scale (check the input, or -cf).")
        elif any(r[5] < SHORT_BOND_FACTOR for r in f["pairs"]):
            warn(f"{os.path.basename(rec['in_file'])}: a bond is shorter than {SHORT_BOND_FACTOR} x the sum of the "
                 "covalent radii -- check the geometry.")

    sections("SYMMETRY" + t)
    rows = [([f"{sp:g}", g or "-", str(num) if num else "-", str(nops) if nops else "-"], None)
            for sp, g, num, nops in s["table"]]
    print_table(["symprec (Ang)", "space group", "number", "operations"], rows, f_out)
    groups = {r[1] for r in s["table"] if r[1]}
    if len(groups) > 1:
        print_dual("Tolerance         : the space group depends on the tolerance -- the structure is close to a higher "
                   "symmetry (or slightly distorted from it)", f_out)
    print_dual(f"Crystal system    : {s['system'] or '-'}; point group {s['point_group'] or '-'}", f_out)
    if s["primitive"]:
        print_dual(f"Primitive cell    : {s['primitive']} atoms (symprec 0.01 Ang)"
                   + (f" -- the input cell is {f['n'] // s['primitive']}x larger" if s["primitive"] < f["n"] else
                      " -- the input cell is primitive"), f_out)

    if rec["out_file"] is None:
        return
    sections("OUTPUT" + t)
    print_dual(f"File              : {rec['out_file']} ({rec.get('out_size', '-')})", f_out)
    print_dual(f"Format            : {rec['out_format']} ({rec['out_how']})", f_out)
    print_dual(f"Coordinates       : {rec.get('coords_out', '-')} "
               + ("" if rec["out_format"] in ("cif", "xyz", "xsf", "lammps") else
                  "(requested, -cf)" if rec.get("cf") else "(as the input; -cf to choose)"), f_out)
    if rec["out_format"] == "cif":
        print_dual("CIF symmetry      : " + (f"found at {CIF_SYMPREC:g} Ang and written as the asymmetric unit in the "
                                             "standard setting -- atoms moved onto their exact symmetric positions; "
                                             "--cif-symprec 0 writes P1 exactly as read" if CIF_SYMPREC else
                                             "none (P1, --cif-symprec 0): cell and positions exactly as read"), f_out)
    if "check_error" in rec:
        print_dual(color_text(f"Verification      : not possible ({rec['check_error']})", 'yellow'), f_out)
        warn(f"{os.path.basename(rec['out_file'])}: the written file could not be re-read for verification: "
             f"{rec['check_error']}")
        return
    ck = rec["check"]
    rows = []
    ok_f = ck["formula_in"] == ck["formula_out"]
    rows.append((["composition", ck["formula_in"], ck["formula_out"], "same" if ok_f else "DIFFERENT"],
                 None if ok_f else "yellow"))
    rows.append((["atoms", str(ck["n_in"]), str(ck["n_out"]),
                  "same" if ck["n_in"] == ck["n_out"] else "re-expressed cell"], None))
    ok_v = np.isclose(ck["vfu_in"], ck["vfu_out"], rtol=1e-3)
    rows.append((["volume / formula unit (Ang^3)", f"{ck['vfu_in']:.4f}", f"{ck['vfu_out']:.4f}",
                  "same" if ok_v else "DIFFERENT"], None if ok_v else "yellow"))
    if ck["note"] is None:
        rows.append((["largest cell difference (Ang)", "", f"{ck['cell_diff']:.2e}", "ok"], None))
    symmetrized = (rec["out_format"] == "cif" and CIF_SYMPREC and ck["pos_diff"] is not None
                   and 1e-3 <= ck["pos_diff"] <= 2 * CIF_SYMPREC)
    if ck["pos_diff"] is not None:
        ok_p = ck["pos_diff"] < 1e-3
        rows.append((["largest atom displacement (Ang)", "", f"{ck['pos_diff']:.2e}",
                      "ok" if ok_p else ("symmetrized" if symmetrized else "CHECK")],
                     None if ok_p or symmetrized else "yellow"))
    if ck.get("match") is not None:
        rms = ck.get("match_rms")
        rows.append((["same crystal (StructureMatcher)", "", "yes" if ck["match"] else "NO",
                      (f"rms {rms:.1e} Ang" if rms is not None else "") if ck["match"] else "DIFFERENT"],
                     None if ck["match"] else "yellow"))
        if not ck["match"]:
            warn(f"{os.path.basename(rec['out_file'])}: the written file is not the same crystal as the input "
                 "(StructureMatcher) -- check the conversion.")
    print_dual("Verification (the written file read back):", f_out)
    print_table(["quantity", "input", "output", "verdict"], rows, f_out)
    if ck["note"]:
        print_dual(f"Note              : {ck['note']}", f_out)
    if not ok_f or not ok_v:
        warn(f"{os.path.basename(rec['out_file'])}: the written file does not reproduce the input "
             f"(composition {ck['formula_in']} -> {ck['formula_out']}, volume/f.u. {ck['vfu_in']:.3f} -> "
             f"{ck['vfu_out']:.3f} Ang^3) -- check the conversion.")
    elif symmetrized:
        print_dual(f"Symmetrized       : the CIF writer moved atoms by up to {ck['pos_diff']:.4f} Ang onto their exact "
                   f"symmetric positions (--cif-symprec {CIF_SYMPREC:g}); --cif-symprec 0 keeps them as read", f_out)
    elif ck["pos_diff"] is not None and ck["pos_diff"] >= 1e-3:
        warn(f"{os.path.basename(rec['out_file'])}: an atom moved by {ck['pos_diff']:.4f} Ang in the written file "
             "(precision of the format, or a wrapping/convention change).")


def chemical_symbol(z):
    from ase.data import chemical_symbols
    return chemical_symbols[z]


def build_parser():
    parser = argparse.ArgumentParser(
        description="Structure file converter: reads one structure (or a glob of them) in one format and writes it in "
                    "another, with a report of what was read and a check of what was written. Both formats are "
                    "detected from the file names when -if/-of are omitted.",
        epilog="Example usage:\n"
               "  stb-translate -i si.cif -o si.fdf                 # both formats from the names\n"
               "  stb-translate -i POSCAR -o structure.xyz -cf cartesian\n"
               "  stb-translate -i siesta.STRUCT_OUT -of cif -o relaxed.cif --save-report\n"
               "  stb-translate -i '*.cif' -of fdf --out-dir fdf/   # batch\n"
               "  stb-translate -i slab.fdf --dry-run                # report only\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-if", "--in-format", required=False, default=None, choices=sorted(INPUT_FORMATS),
                        help="Input format: poscar, cif, siesta [SIESTA STRUCT_OUT], xyz, fhi, dftb, xsf, fdf [SIESTA "
                             "input]. Default: from --in-file's name (.vasp/POSCAR/CONTCAR, .cif, .STRUCT_OUT, .xyz, "
                             ".in, .gen, .xsf, .fdf), per file in batch mode.")
    parser.add_argument("-i", "--in-file", required=True,
                        help="Input file, or a glob pattern (e.g. '*.cif') for batch conversion.")
    parser.add_argument("-of", "--out-format", required=False, default=None, choices=sorted(OUTPUT_FORMATS),
                        help="Output format: cif, xyz, poscar, fdf, dftb, xsf, fhi, lammps [LAMMPS data, atom_style "
                             "atomic, metal units]. Default: from --out-file's name (.cif, .xyz, .fdf, .vasp/.poscar/"
                             "POSCAR, .gen, .in, .xsf, .data/.lmp). Required in batch mode.")
    parser.add_argument("-o", "--out-file", required=False, default=None,
                        help="Output file. Required unless --dry-run, or in batch mode (--out-dir instead).")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory in batch mode: each output is named after its input with the output "
                             "format's extension.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Read and report only (structure, symmetry, checks); write nothing.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also write the report to {REPORT_FILE} (next to --out-file, or in --out-dir).")
    parser.add_argument("-cf", "--coord-format", choices=['direct', 'cartesian'], default=None,
                        help="Coordinates of the output: direct (fractional) or cartesian. Default: the input's, or the "
                             "output format's own.")
    parser.add_argument("--cif-symprec", type=float, default=0.01, metavar="ANG",
                        help="CIF output only: tolerance with which the space group is found and written (default: "
                             "0.01 Ang). The atoms are moved onto their exact symmetric positions (by up to this much) "
                             "and the asymmetric unit is written in the standard setting -- what a database or a "
                             "paper expects. 0 writes P1: the cell and the positions exactly as read.")
    parser.add_argument("--lattice-vectors", type=float, nargs=9, default=None,
                        metavar=("AX", "AY", "AZ", "BX", "BY", "BZ", "CX", "CY", "CZ"),
                        help="The 3 lattice vectors (Cartesian, Ang), needed for an XYZ input, e.g. 5.43 0 0 0 5.43 0 0 "
                             "0 5.43. Applied to every file in batch mode.")
    parser.add_argument("--view", action="store_true",
                        help="Open the structure in ASE's 3D viewer at the end (needs a display).")
    parser.add_argument("--view-image", default=None, metavar="PATH",
                        help="Save a 3D snapshot to PATH (a directory in batch mode, one PNG per file). No display "
                             "needed.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-translate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    global CIF_SYMPREC
    CIF_SYMPREC = args.cif_symprec if args.cif_symprec and args.cif_symprec > 0 else None

    if args.intro:
        show_intro(["Siesta ToolBox Suite - File Translator",
                    "Structure files between SIESTA, VASP, CIF, XYZ, XSF, FHI-aims, DFTB+ and LAMMPS",
                    f"Version {VERSION} | University of Brasilia - 2026",
                    "Developed by Dr. Carlos M. O. Bastos"])

    is_glob = any(ch in args.in_file for ch in '*?[')
    matches = sorted(glob.glob(args.in_file)) if is_glob else [args.in_file]
    if not matches:
        parser.error(f"No files matched '{args.in_file}'.")
    batch = len(matches) > 1

    # input format: given, else per file from its name
    in_formats = {}
    for m in matches:
        fmt = args.in_format or guess_format_from_filename(m)
        if fmt is None:
            parser.error(f"Could not infer the input format of '{m}' from its name; pass -if "
                         f"({', '.join(sorted(INPUT_FORMATS))}).")
        in_formats[m] = fmt
    if any(f == "xyz" for f in in_formats.values()) and not args.lattice_vectors:
        parser.error("An XYZ input needs --lattice-vectors (XYZ carries no cell).")

    # output format: given, else from --out-file's name
    out_how = "given (-of)"
    if batch:
        if args.out_file:
            parser.error(f"--out-file can't be used with {len(matches)} input files; use --out-dir.")
        if not args.dry_run and (args.out_format is None or args.out_dir is None):
            parser.error("In batch mode -of/--out-format and --out-dir are required (unless --dry-run).")
        if not args.dry_run:
            os.makedirs(args.out_dir, exist_ok=True)
    elif not args.dry_run:
        if args.out_file is None:
            parser.error("--out-file (-o) is required unless --dry-run.")
        if args.out_format is None:
            args.out_format = guess_output_format(args.out_file)
            if args.out_format is None:
                parser.error(f"Could not infer the output format from '{args.out_file}'; pass -of "
                             f"({', '.join(sorted(OUTPUT_FORMATS))}).")
            out_how = f"auto, from the name '{os.path.basename(args.out_file)}'"

    report_path = None
    if not args.dry_run and args.save_report:
        report_dir = args.out_dir if batch else (os.path.dirname(args.out_file) or ".")
        report_path = os.path.join(report_dir, REPORT_FILE)
    f_out = open(report_path, "w") if report_path else None
    warns = []
    sections = Sections(f_out)

    def warn(msg):
        warns.append(msg)

    records = []
    try:
        print("\n" + color_text("TRANSLATE:", 'bold'))
        print("-" * 60)
        print_dual(color_text("===== STB-TRANSLATE REPORT =====", 'magenta'), f_out)
        sections("RUN METADATA")
        print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
        try:
            import ase
            import spglib
            from importlib.metadata import version as _v
            print_dual(f"Software          : stb {VERSION}, ase {ase.__version__}, pymatgen {_v('pymatgen')}, "
                       f"spglib {spglib.__version__}", f_out)
        except Exception:
            pass
        if batch:
            fmts = sorted(set(in_formats.values()))
            print_dual(f"Input             : '{args.in_file}' -- {len(matches)} files, format "
                       f"{', '.join(fmts)} ({'given (-if)' if args.in_format else 'auto, per file name'})", f_out)
            print_dual("Output            : " + ("none (--dry-run)" if args.dry_run else
                                                 f"{args.out_dir}/ as {args.out_format} (given (-of))"), f_out)
        else:
            m = matches[0]
            print_dual(f"Input             : {m} -- {in_formats[m]} ({_how(args.in_format, os.path.basename(m))})", f_out)
            print_dual("Output            : " + ("none (--dry-run)" if args.dry_run else
                                                 f"{args.out_file} -- {args.out_format} ({out_how})"), f_out)
        if args.coord_format:
            print_dual(f"Coordinates       : {args.coord_format} (-cf)", f_out)
        if not args.dry_run and (args.out_format == "cif"):
            print_dual(f"CIF symmetry      : " + (f"{CIF_SYMPREC:g} Ang (--cif-symprec)" if CIF_SYMPREC else
                                                 "none, P1 (--cif-symprec 0)"), f_out)
        if args.lattice_vectors:
            print_dual("Lattice vectors   : from --lattice-vectors (XYZ input)", f_out)
        if report_path:
            print_dual(f"Report file       : {report_path}", f_out)

        for m in matches:
            if args.dry_run:
                out_file, out_fmt = None, None
            elif batch:
                base = os.path.splitext(os.path.basename(m))[0]
                out_file, out_fmt = os.path.join(args.out_dir, f"{base}.{OUTPUT_EXTENSION_MAP[args.out_format]}"), args.out_format
            else:
                out_file, out_fmt = args.out_file, args.out_format
            print(f"  ... reading {m}", flush=True)
            rec = convert_one(args, m, in_formats[m], out_file, out_fmt)
            rec["in_how"] = _how(args.in_format, os.path.basename(m))
            rec["out_how"] = "given (-of)" if batch else out_how
            rec["cf"] = args.coord_format
            records.append(rec)
            report_file(rec, sections, f_out, warn, tag=batch)

        ok = [r for r in records if r["ok"]]
        bad = [r for r in records if not r["ok"]]
        sections("SUMMARY & FILES")
        if batch:
            rows = []
            for r in records:
                if r["ok"]:
                    sg = next((g for sp, g, _n, _o in r["sym"]["table"] if sp == 1e-2 and g), "-")
                    rows.append(([os.path.basename(r["in_file"]), r["atoms"].get_chemical_formula(), str(len(r["atoms"])),
                                  sg, os.path.basename(r["out_file"]) if r["out_file"] else "-", "ok"], None))
                else:
                    rows.append(([os.path.basename(r["in_file"]), "-", "-", "-", "-", "FAILED"], "red"))
            print_table(["input", "formula", "atoms", "space group (0.01)", "output", "status"], rows, f_out)
        status = "OK" if not bad else (f"{len(bad)} of {len(records)} FAILED" if batch else "FAILED")
        print_dual(f"Status            : {status}" + (f" ({len(warns)} warning(s), see the last section)" if warns else ""),
                   f_out)
        if not batch and ok:
            r = ok[0]
            sg = next((g for sp, g, _n, _o in r["sym"]["table"] if sp == 1e-2 and g), "-")
            print_dual(f"Structure         : {r['atoms'].get_chemical_formula()}, {len(r['atoms'])} atoms, {sg}", f_out)
        written = [r["out_file"] for r in ok if r["out_file"]]
        if written:
            print_dual(f"Written           : {', '.join(written)}", f_out)
        if report_path:
            print_dual(f"Report            : {report_path}", f_out)
        not_asked = [fl for fl, on in (("--save-report", args.save_report), ("--view", args.view),
                                       ("--view-image", bool(args.view_image)), ("-cf", bool(args.coord_format)))
                     if not on]
        if not_asked and not args.dry_run:
            print_dual(f"Also available    : {' '.join(not_asked)}", f_out)
        sections("WARNINGS")
        if warns:
            for i, w in enumerate(warns, 1):
                print_dual(color_text(f"  {i}. {w}", 'yellow'), f_out)
        else:
            print_dual(color_text("  none", 'green'), f_out)
    finally:
        if f_out:
            f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Converting input files is 10% coding, 90% crying.\n", 'bold'))

    # 3D visualization: the very last step, so a blocking window never hides the report
    viewable = [(r["in_file"], r["atoms"]) for r in records if r["ok"] and r["atoms"] is not None]
    if (args.view or args.view_image) and viewable:
        print(color_text("3D VISUALIZATION", 'magenta'))
        print("-" * 60)
        if args.view_image and batch:
            os.makedirs(args.view_image, exist_ok=True)
        for in_file, atoms in viewable:
            label = f" ({in_file})" if batch else ""
            if args.view_image:
                image_path = args.view_image
                if batch:
                    image_path = os.path.join(args.view_image, f"{os.path.splitext(os.path.basename(in_file))[0]}.png")
                save_structure_image(atoms, image_path)
                print(f"[OK] Saved 3D structure image to {image_path}{label}")
            if args.view:
                print(f"[INFO] Opening interactive 3D viewer (ASE GUI){label} -- close the window to continue...")
                view_structure_interactive(atoms)

    if any(not r["ok"] for r in records):
        sys.exit(1)


if __name__ == "__main__":
    main()
