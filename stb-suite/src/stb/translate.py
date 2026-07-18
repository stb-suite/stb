#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.1"    

import os
import sys
import glob
import warnings
import subprocess
from time import sleep
from datetime import datetime
import argparse
import textwrap
from typing import List, Dict
import argparse
from ase.io import read as ase_read
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section
from stb.core.symmetry import crystal_system, space_group_label

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

# Report filename written alongside the output (--out-dir in batch mode, or
# --out-file's own directory otherwise) -- same print_dual-to-screen-and-file
# convention as the workflow analysis stages (e.g. stb-ramanAnalysis's
# raman_stage3.txt). Not written for --dry-run, which by definition has no
# disk side effects.
REPORT_FILE = "stb_translate_report.txt"


def guess_format_from_filename(path):
    """Best-effort --in-format guess from a filename; None if not recognized."""
    base = os.path.basename(path).lower()
    if base in BASENAME_FORMAT_MAP:
        return BASENAME_FORMAT_MAP[base]
    _, ext = os.path.splitext(base)
    return EXTENSION_FORMAT_MAP.get(ext)


def print_structure_summary(typevectors, latticeparameter, vectors, getatoms, atoms=None, f_out=None):
    """Prints a structure summary: formula, atom count, coordinate type,
    lattice vectors, cell volume, and (when `atoms` is given) the detected
    space group/crystal system/point group. Shown for every successful read
    (not just --dry-run), via print_dual so it's also captured in the run's
    persisted REPORT_FILE."""
    formula = " ".join(f"{elem[2]}{elem[3]}" for elem in getatoms)
    natoms_total = sum(int(elem[3]) for elem in getatoms)
    vectors_np = np.array(vectors, dtype=float) * float(latticeparameter)
    volume = abs(np.linalg.det(vectors_np))

    print_dual(f"Formula           : {formula}", f_out)
    print_dual(f"Total atoms       : {natoms_total}", f_out)
    print_dual(f"Coordinate type   : {typevectors}", f_out)
    print_dual(f"Lattice parameter : {latticeparameter}", f_out)
    print_dual("Lattice vectors (Ang):", f_out)
    for row in vectors_np:
        print_dual(f"    {row[0]:12.8f}  {row[1]:12.8f}  {row[2]:12.8f}", f_out)
    print_dual(f"Cell volume       : {volume:.4f} Ang^3", f_out)

    if atoms is not None:
        try:
            structure = AseAtomsAdaptor.get_structure(atoms)
            sg_label = space_group_label(structure)
            system, point_group = crystal_system(structure)
            print_dual(f"Space group       : {sg_label}", f_out)
            print_dual(f"Crystal system    : {system}", f_out)
            print_dual(f"Point group       : {point_group}", f_out)
        except Exception as e:
            print_dual(color_text(f"[WARNING] Could not determine symmetry: {e}", 'yellow'), f_out)


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
    # Optional "Selective dynamics" line (VASP only requires the first letter to
    # match) shifts the coordinate-type line and the start of positions by one.
    coord_line_idx = 8 if datavasp[7][0].lower().startswith('s') else 7
    typevectors = datavasp[coord_line_idx][0]
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
    getatoms = []
    for atoms in dicatoms:
        getatoms.append([dicatoms[atoms][0][0],
                        dicatoms[atoms][0][1],
                        dicatoms[atoms][0][2],
                        str(len(dicatoms[atoms]))])
    atomic_position = {}
    for position in datasiesta[4:]:
        if atomicnumber[position[1]] not in atomic_position:
            atomic_position[atomicnumber[position[1]]] = []
        atomic_position[atomicnumber[position[1]]].append(
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
    for i in range(1, 4):
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

    atomic_position = {}
    for symbol, pos in structure.atoms:
        atomic_position.setdefault(symbol, []).append([f"{pos[0]:.8f}", f"{pos[1]:.8f}", f"{pos[2]:.8f}"])

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
    atoms = [
        (elem[2], np.array([float(p[0]), float(p[1]), float(p[2])]))
        for elem in getatoms
        for p in final_positions[elem[2]]
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
    CifWriter(structure, symprec=1e-3).write_file(outfilename)




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
def convert_coordinates(typevectors_in: str,
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
        print("[INFO] Converting coordinates from Direct -> Cartesian...")
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
        print("[INFO] Converting coordinates from Cartesian -> Direct...")
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
    for the optional --view/--view-image 3D visualization."""
    from ase import Atoms

    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'cartesian'
    )
    cell = np.array(vectors, dtype=float) * float(latticeparameter)

    symbols = []
    positions = []
    for elem in getatoms:
        species = elem[2]
        for pos in final_positions[species]:
            symbols.append(species)
            positions.append([float(pos[0]), float(pos[1]), float(pos[2])])

    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)


def view_structure_interactive(atoms):
    """Opens ASE's interactive 3D structure viewer. Needs a display (local X11,
    or `ssh -X`/`-Y` to a remote machine) and a working Tk installation."""
    try:
        from ase.visualize import view
        # block=True: waits for the GUI subprocess so a failure (e.g. no display)
        # raises here and is reported, instead of failing silently in the background.
        view(atoms, block=True)
    except Exception as e:
        print(color_text(f"[FAIL] Could not open the interactive 3D viewer: {e}", 'red'))
        print(color_text(
            "       This needs a display (local X11, or `ssh -X`/`-Y` to a remote "
            "machine) and a working Tk installation. Try --view-image instead to "
            "save a static picture without needing a display.", 'yellow'))


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


##### INTEGRITY CHECKS (close contacts + post-conversion round-trip) #####
def check_close_contacts(atoms, threshold=0.5, f_out=None):
    """Reports on any pair of atoms closer than `threshold` Angstrom (mic
    -aware, via ASE's get_all_distances) -- a common symptom of mixing up
    Direct/Cartesian coordinates or the wrong lattice scale. Not fatal: an
    intentionally unrelaxed starting guess can legitimately have short
    contacts, so this only warns -- but always prints a result (clean or
    not), so the check is visibly part of the report rather than silent."""
    n = len(atoms)
    if n < 2:
        print_dual("[OK] Close-contact check skipped (fewer than 2 atoms).", f_out)
        return
    distances = atoms.get_all_distances(mic=True)
    symbols = atoms.get_chemical_symbols()
    close_pairs = [(i, j, distances[i, j])
                   for i in range(n) for j in range(i + 1, n)
                   if distances[i, j] < threshold]
    if close_pairs:
        print_dual(color_text(
            f"[WARNING] {len(close_pairs)} suspiciously close contact(s) found "
            f"(< {threshold} Ang) -- check for a Direct/Cartesian mixup or a wrong "
            f"lattice scale:", 'yellow'), f_out)
        for i, j, d in close_pairs:
            print_dual(color_text(f"           {symbols[i]}{i} -- {symbols[j]}{j}: {d:.4f} Ang", 'yellow'), f_out)
    else:
        print_dual(f"[OK] No suspiciously close contacts found (< {threshold} Ang).", f_out)


def reread_for_check(out_format, out_file, vectors, latticeparameter):
    """Re-parses the just-written output file with the matching reader and
    returns an ase.Atoms, for the post-conversion round-trip check below.
    `vectors`/`latticeparameter` are the ORIGINAL input's, reused for an xyz
    output re-read since an .xyz file itself carries no cell information.

    Every format except lammps goes through this file's own
    getatomsandvectors_* + build_ase_atoms (the same pair every other format
    conversion uses); lammps has no getatomsandvectors_lammps (it's an
    output-only format here) and ASE can already read a lammps-data file
    directly, so that one path is a straight ase.io.read call instead.
    """
    if out_format == "lammps":
        from ase.io import read as ase_read_atoms
        return ase_read_atoms(out_file, format='lammps-data')

    readers = {
        "poscar": lambda: getatomsandvectors_vasp(out_file),
        "cif": lambda: getatomsandvectors_cif(out_file),
        "xyz": lambda: getatomsandvectors_xyz(out_file, vectors),
        "fhi": lambda: getatomsandvectors_fhi(out_file),
        "dftb": lambda: getatomsandvectors_dftb(out_file),
        "xsf": lambda: getatomsandvectors_xsf(out_file),
        "fdf": lambda: getatomsandvectors_fdf(out_file),
    }
    return build_ase_atoms(*readers[out_format]())


def verify_round_trip(atoms_in, atoms_out, f_out=None):
    """Compares the pre-write structure against the just-written file
    (re-read via reread_for_check): reduced-formula stoichiometry, and cell
    volume PER FORMULA UNIT. Prints [OK] (confirms it verified) when they
    match, [WARNING] when not -- non-fatal, since the output file is already
    written by the time this runs.

    Deliberately normalized by formula unit (not raw atom count/cell volume):
    a symmetry-aware writer (writefilecif, via pymatgen's CifWriter) can
    legitimately re-express the same crystal in a differently-sized standard/
    conventional cell (e.g. a small hand-built cell expanded to the correct
    conventional cell for its space group) -- verified live against the
    suite's own Si test fixture, where CifWriter re-expressed a 2-atom cell
    as a 6-atom rhombohedral setting (exact 3x volume/atom-count): comparing
    raw totals flagged that as a false-positive mismatch; per-formula-unit
    comparison correctly reports it as verified."""
    from pymatgen.core import Composition

    comp_in = Composition(atoms_in.get_chemical_formula())
    comp_out = Composition(atoms_out.get_chemical_formula())
    formula_in, formula_out = comp_in.reduced_formula, comp_out.reduced_formula
    _, z_in = comp_in.get_reduced_composition_and_factor()
    _, z_out = comp_out.get_reduced_composition_and_factor()
    vol_in = abs(np.linalg.det(atoms_in.get_cell())) / z_in
    vol_out = abs(np.linalg.det(atoms_out.get_cell())) / z_out

    if formula_in == formula_out:
        print_dual(color_text(
            f"[OK] Stoichiometry verified: {formula_in} ({len(atoms_in)} -> {len(atoms_out)} "
            f"atoms) -- matches input.", 'green'), f_out)
    else:
        print_dual(color_text(
            f"[WARNING] Stoichiometry mismatch: input is {formula_in}, output is "
            f"{formula_out} -- check the conversion!", 'yellow'), f_out)

    if np.isclose(vol_in, vol_out, rtol=1e-3):
        print_dual(color_text(
            f"[OK] Cell volume verified: {vol_out:.4f} Ang^3 per formula unit -- "
            f"matches input.", 'green'), f_out)
    else:
        print_dual(color_text(
            f"[WARNING] Cell volume mismatch (per formula unit): input is {vol_in:.4f} "
            f"Ang^3, output is {vol_out:.4f} Ang^3 -- check the conversion!", 'yellow'), f_out)
##### END INTEGRITY CHECKS #####


def convert_one(args, in_file, out_file, f_out=None, header=None):
    """Runs the read -> structure/symmetry summary -> integrity checks ->
    (optional) write -> round-trip verification pipeline for a single file,
    printing a report block via print_dual (screen + the run's REPORT_FILE).
    `in_file`/`out_file` are passed explicitly (rather than read off `args`)
    since they vary per file in batch mode; every other option (in_format,
    out_format, coord_format, lattice_vectors, dry_run) is shared across the
    whole batch. `header`, if given, labels this file's block (batch mode).

    Does NOT do --view/--view-image -- that's deferred by main() to run once,
    after every file has been read/written/verified and reported, as the
    final step of the whole run.

    Returns (ok, atoms): ok is False on failure (prints its own [FAIL], never
    raises/exits, so a batch run can continue past one bad file); atoms is
    the parsed ase.Atoms on success (for the deferred view step), else None.
    """
    if header:
        print_dual(f"\n--- {header} ---", f_out)
    try:
        ##### MODIFICADO #####
        # Adicionado o 'case' para o novo formato 'fdf'
        match (args.in_format):
            case "poscar":
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_vasp(
                    in_file)
            case "cif":
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_cif(
                    in_file)
            case "siesta":
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_siesta(
                    in_file)
            case "xyz":
                v = args.lattice_vectors
                xyz_vectors = [[str(v[0]), str(v[1]), str(v[2])],
                               [str(v[3]), str(v[4]), str(v[5])],
                               [str(v[6]), str(v[7]), str(v[8])]]
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xyz(
                    in_file, xyz_vectors)
            case "fhi":
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fhi(
                    in_file)
            case "dftb":
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_dftb(
                    in_file)
            case "xsf":
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xsf(
                    in_file)
            case "fdf": ##### NOVO #####
                typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fdf(
                    in_file)


        print_dual(f"[OK] Read the file {in_file} ({args.in_format})", f_out)

        atoms = build_ase_atoms(typevectors, latticeparameter, vectors, getatoms, atomsposition)
        print_structure_summary(typevectors, latticeparameter, vectors, getatoms, atoms, f_out)
        check_close_contacts(atoms, f_out=f_out)

        if args.dry_run:
            print_dual("[INFO] Dry run -- no output file written.", f_out)
            return True, atoms

        ##### MODIFIED #####
        # All write functions now receive 'args.coord_format'
        match (args.out_format):
            case "xyz":
                writefilexyz(typevectors, latticeparameter, vectors,
                             getatoms, atomsposition, out_file, args.coord_format)
            case "poscar":
                writefileposcar(typevectors, latticeparameter, vectors,
                                getatoms, atomsposition, out_file, args.coord_format)
            case "fdf":
                writefilefdf(typevectors, latticeparameter, vectors,
                             getatoms, atomsposition, out_file, args.coord_format)
            case "dftb":
                writefiledftb(typevectors, latticeparameter, vectors,
                              getatoms, atomsposition, out_file, args.coord_format)
            case "xsf":
                writefilexsf(typevectors, latticeparameter, vectors,
                             getatoms, atomsposition, out_file, args.coord_format)
            case "fhi":
                writefilefhi(typevectors, latticeparameter, vectors,
                             getatoms, atomsposition, out_file, args.coord_format)

            case "cif":
                writefilecif(typevectors, latticeparameter, vectors,
                             getatoms, atomsposition, out_file, args.coord_format)
            case "lammps":
                writefilelammps(typevectors, latticeparameter, vectors,
                                getatoms, atomsposition, out_file, args.coord_format)
    except FileNotFoundError as e:
        print_dual(color_text(f"[FAIL] File not found: {e}", 'red'), f_out)
        return False, None
    except (IndexError, KeyError, ValueError) as e:
        print_dual(color_text(
            f"[FAIL] Could not parse '{in_file}' as '{args.in_format}': {e}",
            'red'), f_out)
        print_dual(color_text(
            "       Check that the file matches the expected structure for --in-format.",
            'yellow'), f_out)
        return False, None

    print_dual(f"[OK] Writing the file {out_file} ({args.out_format})", f_out)

    try:
        atoms_out = reread_for_check(args.out_format, out_file, vectors, latticeparameter)
        verify_round_trip(atoms, atoms_out, f_out=f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Could not verify output file integrity: {e}", 'yellow'), f_out)

    return True, atoms


def main():

    parser = argparse.ArgumentParser(
        description="File format converter using stb-translate."
    )

    ##### MODIFICADO #####
    # in-format agora é opcional: se omitido, é inferido da extensão do --in-file
    parser.add_argument("-if", "--in-format", required=False, default=None, choices=INPUT_FORMATS,
                        help="Input file format (options: poscar, cif, "
                             "siesta [SIESTA STRUCT_OUT-style output], xyz, fhi, dftb, xsf, "
                             "fdf [SIESTA .fdf input file]). If omitted, inferred from "
                             "--in-file's name/extension.")
    parser.add_argument("-i", "--in-file", required=True,
                        help="Path to the input file, or a glob pattern (e.g. '*.poscar') "
                             "matching several files for batch conversion.")
    parser.add_argument("-of", "--out-format", required=False, default=None, choices=OUTPUT_FORMATS,
                        help="Output file format (options: cif, xyz, poscar, fdf, dftb, xsf, fhi, "
                             "lammps [LAMMPS data file, atom_style atomic, metal units -- for "
                             "handing a relaxed structure off to classical MD/force-field runs]). "
                             "Required unless --dry-run.")
    parser.add_argument("-o", "--out-file", required=False, default=None,
                        help="Path to the output file. Required unless --dry-run, or in batch "
                             "mode (use --out-dir instead).")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory for batch conversion, used when --in-file is a "
                             "glob pattern matching more than one file. Each output is named "
                             "after its input's basename with the output format's extension. "
                             "Required in batch mode unless --dry-run; ignored otherwise.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Parse --in-file and print a summary (formula, atom count, "
                             "coordinate type, lattice vectors, cell volume, symmetry) without "
                             "writing any output file. --out-format/--out-file are not needed.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE} (next to --out-file, "
                             "or in --out-dir in batch mode), not just print it. Off by default; "
                             "ignored with --dry-run (no disk side effects there regardless).")

    ##### NEW ARGUMENT #####
    parser.add_argument(
        "-cf", "--coord-format",
        choices=['direct', 'cartesian'],
        default=None,
        help="Specify the output coordinate format (direct/fractional or cartesian). "
             "If not specified, uses the input format or the output format's default."
    )
    ##### END NEW ARGUMENT #####

    parser.add_argument(
        "--lattice-vectors", type=float, nargs=9, default=None,
        metavar=("AX", "AY", "AZ", "BX", "BY", "BZ", "CX", "CY", "CZ"),
        help="The 3 lattice vectors in Cartesian Angstrom, required only for XYZ "
             "input, e.g. --lattice-vectors 5.43 0 0 0 5.43 0 0 0 5.43. Applied to "
             "every matched file in batch mode.")

    ##### OPTIONAL 3D VIEWER #####
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D viewer (ASE GUI) for the parsed "
                             "structure. Needs a display (local X11, or ssh -X/-Y). Runs "
                             "as the very last step, once every file has been read, "
                             "written and verified.")
    parser.add_argument("--view-image", default=None, metavar="PATH",
                        help="Save a static 3D snapshot of the parsed structure to PATH "
                             "(e.g. structure.png). In batch mode, PATH is a directory "
                             "instead, and each file's snapshot is named after its stem. "
                             "Works without a display; can be used together with, or "
                             "instead of, --view. Runs as the very last step.")
    ##### END OPTIONAL 3D VIEWER #####

    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-translate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("TRANSLATE:", 'bold'))
    print("-"*60)

    # --in-file as a glob pattern (batch mode) vs. a plain path (single file)
    is_glob = any(c in args.in_file for c in '*?[')
    matches = sorted(glob.glob(args.in_file)) if is_glob else [args.in_file]
    if not matches:
        parser.error(f"No files matched '{args.in_file}'.")
    batch_mode = len(matches) > 1

    # Infer --in-format from the first matched file's name when omitted
    if args.in_format is None:
        args.in_format = guess_format_from_filename(matches[0])
        if args.in_format is None:
            parser.error(
                f"Could not infer --in-format from '{matches[0]}'; pass it explicitly "
                f"(options: {', '.join(sorted(INPUT_FORMATS))}).")

    if batch_mode:
        if args.out_file:
            parser.error(
                "--out-file can't be used when --in-file matches multiple files "
                f"({len(matches)} matched '{args.in_file}'); use --out-dir instead.")
        if not args.dry_run:
            if args.out_format is None or args.out_dir is None:
                parser.error(
                    "--out-format and --out-dir are required in batch mode, unless --dry-run.")
            os.makedirs(args.out_dir, exist_ok=True)
    else:
        if not args.dry_run and (args.out_format is None or args.out_file is None):
            parser.error(
                "--out-format and --out-file are required unless --dry-run is set.")

    # Validate lattice vectors requirement (shared across every matched file)
    if args.in_format == "xyz" and not args.lattice_vectors:
        parser.error(
            "The --lattice-vectors argument is required when input format is XYZ.")

    # Persisted report (print_dual'd alongside the screen output), same
    # convention as the workflow analysis stages (e.g. stb-ramanAnalysis's
    # raman_stage3.txt) -- opt-in via --save-report (off by default), and
    # never written for --dry-run regardless, which by definition has no
    # disk side effects.
    report_path = None
    if not args.dry_run and args.save_report:
        report_dir = args.out_dir if batch_mode else (os.path.dirname(args.out_file) or ".")
        report_path = os.path.join(report_dir, REPORT_FILE)
    f_out = open(report_path, "w") if report_path else None

    try:
        print_dual(color_text("===== STB-TRANSLATE REPORT =====", 'magenta'), f_out)

        print_section("[0] RUN METADATA", f_out)
        print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
        if batch_mode:
            print_dual(f"Input pattern     : '{args.in_file}' ({len(matches)} file(s) "
                        f"matched, {args.in_format})", f_out)
        else:
            print_dual(f"Input file        : {matches[0]} ({args.in_format})", f_out)
        if args.dry_run:
            print_dual("Output            : (dry run -- no file written)", f_out)
        elif batch_mode:
            print_dual(f"Output directory  : {args.out_dir} ({args.out_format})", f_out)
        else:
            print_dual(f"Output file       : {args.out_file} ({args.out_format})", f_out)
        if args.coord_format:
            print_dual(f"Coordinate format : {args.coord_format} (requested)", f_out)
        if args.in_format == "xyz":
            print_dual("Lattice vectors   : inline (--lattice-vectors)", f_out)
        if report_path:
            print_dual(f"Report file       : {report_path}", f_out)

        print_section("[1] READING, STRUCTURE & INTEGRITY CHECKS", f_out)
        results = []  # (in_file, ok, atoms, out_file)
        if not batch_mode:
            in_file = matches[0]
            ok, atoms = convert_one(args, in_file, args.out_file, f_out)
            results.append((in_file, ok, atoms, args.out_file))
        else:
            for in_file in matches:
                if args.dry_run:
                    out_file = None
                else:
                    base = os.path.splitext(os.path.basename(in_file))[0]
                    out_file = os.path.join(args.out_dir, f"{base}.{OUTPUT_EXTENSION_MAP[args.out_format]}")
                ok, atoms = convert_one(args, in_file, out_file, f_out, header=in_file)
                results.append((in_file, ok, atoms, out_file))

        ok_results = [r for r in results if r[1]]
        failed_results = [r for r in results if not r[1]]

        print_section("[2] SUMMARY & FILES", f_out)
        if batch_mode:
            print_dual(f"Status            : {len(ok_results)} ok, {len(failed_results)} failed", f_out)
            if failed_results:
                print_dual(color_text(
                    "Failed            : " + ", ".join(r[0] for r in failed_results), 'red'), f_out)
        else:
            print_dual(f"Status            : {'OK' if ok_results else 'FAILED'}", f_out)
        if not args.dry_run:
            written = [r[3] for r in ok_results]
            if written:
                print_dual(f"Files written     : {', '.join(written)}", f_out)
        if report_path:
            print_dual(f"Report            : {report_path}", f_out)
    finally:
        if f_out:
            f_out.close()

    print("\n" + "-"*60)
    print(color_text("Converting input files is 10% coding, 90% crying.\n", 'bold'))

    # [3] 3D VISUALIZATION -- deliberately the LAST step of the whole run,
    # after every file has been read, written and verified/reported above, so
    # a blocking --view window never delays or hides any check that follows it.
    viewable = [(in_file, atoms) for in_file, ok, atoms, _ in results if ok and atoms is not None]
    if (args.view or args.view_image) and viewable:
        print("\n" + color_text("[3] 3D VISUALIZATION", 'magenta'))
        print("-" * 60)
        if args.view_image and batch_mode:
            os.makedirs(args.view_image, exist_ok=True)
        for in_file, atoms in viewable:
            label = f" ({in_file})" if batch_mode else ""
            if args.view_image:
                image_path = args.view_image
                if batch_mode:
                    stem = os.path.splitext(os.path.basename(in_file))[0]
                    image_path = os.path.join(args.view_image, f"{stem}.png")
                save_structure_image(atoms, image_path)
                print(f"[OK] Saved 3D structure image to {image_path}{label}")
            if args.view:
                print(f"[INFO] Opening interactive 3D viewer (ASE GUI){label} "
                      "-- close the window to continue...")
                view_structure_interactive(atoms)

    if failed_results:
        sys.exit(1)

if __name__ == "__main__":
    main()
