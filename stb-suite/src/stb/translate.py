#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.8.0"    

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import argparse
from ase.io import read as ase_read
import numpy as np
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro

# Import libraries:

# Supported formats
##### MODIFICADO #####
# Adicionado "fdf" aos formatos de entrada
INPUT_FORMATS = {"poscar", "cif", "siesta", "xyz", "fhi", "dftb", "xsf", "fdf"} 
OUTPUT_FORMATS = {"cif","xyz", "poscar", "fdf", "dftb", "xsf", "fhi"}


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
def getatomsandvectors_xyz(dataxyz, latticedata):
    element, atomicnumber = periodic_table()
    xyz = readfile(dataxyz)[2:]
    lattice = readfile(latticedata)
    # Lattice
    latticeparameter = lattice[0][0]
    typevectors = 'Cartesian'
    vectors = [lattice[1], lattice[2], lattice[3]]
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
    typevectors = datavasp[7][0]
    getatoms = []
    for i in range(len(datavasp[5])):
        getatoms.append([i+1, element[datavasp[5][i]],datavasp[5][i], datavasp[6][i]])
    atomsposition = []
    cont = 8 
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
    try:
        structure = ase_read(input_cif)
    except NameError:
        # Este erro acontece se a importação no Passo 1 falhar
        print(color_text("[ERRO] A função 'ase_read' não foi encontrada.", 'red'))
        print(color_text("       Verifique se adicionou 'from ase.io import read as ase_read' no topo do script.", 'yellow'))
        sys.exit(1)
    except Exception as e:
        # Captura outros erros de leitura da ASE
        print(color_text(f"[ERRO] Falha ao ler o arquivo CIF com ASE: {e}", 'red'))
        print(color_text("       Verifique se a biblioteca 'ase' está instalada (pip install ase)", 'yellow'))
        sys.exit(1)


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
        getatoms.append([dicatoms[atoms[0]][0][0],
                        dicatoms[atoms[0]][0][1],
                        dicatoms[atoms[0]][0][2],
                        str(len(dicatoms[atoms[0]]))])
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
        print("[FAIL] Type of coordinate not define: Only S or F are accepted")
        exit()
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

    # Calcula ângulos entre vetores
def angle(u, v):
        cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
        

def writefilecif(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    """
    Writes a CIF file in the standard format.
    Note: CIF *always* uses fractional (Direct) coordinates.
    """
    
    # --- NEW: Conversion logic for CIF ---
    if coord_format and coord_format.lower() == 'cartesian':
        print("[WARNING] CIF format requires fractional (Direct) coordinates.")
        print("[INFO]    Ignoring '--coord-format cartesian' and converting to Direct.")
    
    # Force conversion to Direct, regardless of user input
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'direct' # Force 'direct'
    )
    # --- END NEW ---

    vectors_np = np.array(vectors, dtype=float) * float(latticeparameter)

    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    outfile.append("data_generated")
    outfile.append("_symmetry_space_group_name_H-M   'P 1'")
    outfile.append("_symmetry_Int_Tables_number      1")
    outfile.append("_cell_length_a    {:.8f}".format(np.linalg.norm(vectors_np[0])))
    outfile.append("_cell_length_b    {:.8f}".format(np.linalg.norm(vectors_np[1])))
    outfile.append("_cell_length_c    {:.8f}".format(np.linalg.norm(vectors_np[2])))

    alpha = angle(vectors_np[1], vectors_np[2])
    beta = angle(vectors_np[0], vectors_np[2])
    gamma = angle(vectors_np[0], vectors_np[1])

    outfile.append("_cell_angle_alpha  {:.8f}".format(alpha))
    outfile.append("_cell_angle_beta   {:.8f}".format(beta))
    outfile.append("_cell_angle_gamma  {:.8f}".format(gamma))
    outfile.append(" ")

    outfile.append("loop_")
    outfile.append("_symmetry_equiv_pos_as_xyz")
    outfile.append("  'x, y, z'")
    outfile.append(" ")

    outfile.append("loop_")
    outfile.append("_atom_site_label")
    outfile.append("_atom_site_type_symbol")
    outfile.append("_atom_site_fract_x")
    outfile.append("_atom_site_fract_y")
    outfile.append("_atom_site_fract_z")

    # --- MODIFIED: Simplified to use final_positions ---
    # Since we know final_positions is 'Direct', the logic is simple
    for elem in getatoms:
        for pos_str_list in final_positions[elem[2]]:
            outfile.append(
                f"{elem[2]}   {elem[2]}   {pos_str_list[0]}   {pos_str_list[1]}   {pos_str_list[2]}"
            )
    # --- END MODIFIED ---

    np.savetxt(outfilename, outfile, fmt='%s')
    return




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
    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    outfile.append(f"{sum}")
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
    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    
    # --- MODIFIED: Use final_type ---
    if final_type == 'Cartesian':
        outfile.append(f"{sum}   S") # S = Cartesian
    elif final_type == 'Direct':
        outfile.append(f"{sum}   F") # F = Fractional
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
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
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
    outfile.append(f"{sum}  1")
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
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
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
                # v_cart = M . v_direct
                cart_pos = np.dot(lattice_matrix, pos_np)
                # Format back to list of strings
                atomsposition_out[symbol].append(
                    [f"{cart_pos[0]:.8f}", f"{cart_pos[1]:.8f}", f"{cart_pos[2]:.8f}"]
                )
        return 'Cartesian', atomsposition_out

    elif type_in_norm == 'Cartesian' and type_out_norm == 'Direct':
        # --- Convert from Cartesian -> Direct ---
        print("[INFO] Converting coordinates from Cartesian -> Direct...")
        # v_direct = M_inv . v_cart
        inv_lattice = np.linalg.inv(lattice_matrix)
        for symbol, positions in atomsposition_in.items():
            atomsposition_out[symbol] = []
            for pos in positions:
                pos_np = np.array(pos, dtype=float)
                direct_pos = np.dot(inv_lattice, pos_np)
                # Format back to list of strings
                atomsposition_out[symbol].append(
                    [f"{direct_pos[0]:.8f}", f"{direct_pos[1]:.8f}", f"{direct_pos[2]:.8f}"]
                )
        return 'Direct', atomsposition_out
    
    # Fallback (should not happen)
    return typevectors_in, atomsposition_in
##### END NEW HELPER FUNCTION #####


def main():

    parser = argparse.ArgumentParser(
        description="File format converter using stb-translate."
    )

    ##### MODIFICADO #####
    # Atualizado o texto de ajuda para incluir 'fdf'
    parser.add_argument("-if", "--in-format", required=True, choices=INPUT_FORMATS,
                        help="Input file format (options: poscar, cif, siesta, xyz, fhi, dftb, xsf, fdf)")
    parser.add_argument("-i", "--in-file", required=True,
                        help="Path to the input file")
    parser.add_argument("-of", "--out-format", required=True, choices=OUTPUT_FORMATS,
                        help="Output file format (options: cif , xyz, poscar, fdf, dftb, xsf, fhi)")
    parser.add_argument("-o", "--out-file", required=True,
                        help="Path to the output file")
    
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
        "--lattice", help="Lattice vectors file, required only for XYZ output")
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




    # Validate lattice parameter requirement
    if args.in_format == "xyz" and not args.lattice:
        parser.error(
            "The --lattice argument is required when input format is XYZ.")

    print(f"\n[INFO] Converting {args.in_file} ({args.in_format}) to {args.out_file} ({args.out_format})...")

    # Add info about coordinate format
    if args.coord_format:
        print(f"[INFO] Requested output coordinate format: {args.coord_format}")


    if args.out_format == "xyz":
        print(f"[INFO] Lattice vector file: {args.lattice}")

    ##### MODIFICADO #####
    # Adicionado o 'case' para o novo formato 'fdf'
    match (args.in_format):
        case "poscar":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_vasp(
                args.in_file)
        case "cif":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_cif(
                args.in_file)
        case "siesta":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_siesta(
                args.in_file)
        case "xyz":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xyz(
                args.in_file, args.lattice)
        case "fhi":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fhi(
                args.in_file)
        case "dftb":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_dftb(
                args.in_file)
        case "xsf":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xsf(
                args.in_file)
        case "fdf": ##### NOVO #####
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fdf(
                args.in_file)
        

    print(f"[OK] Read the file {args.in_file} ({args.in_format})")

    ##### MODIFIED #####
    # All write functions now receive 'args.coord_format'
    match (args.out_format):
        case "xyz":
            writefilexyz(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)
        case "poscar":
            writefileposcar(typevectors, latticeparameter, vectors,
                            getatoms, atomsposition, args.out_file, args.coord_format)
        case "fdf":
            writefilefdf(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)
        case "dftb":
            writefiledftb(typevectors, latticeparameter, vectors,
                          getatoms, atomsposition, args.out_file, args.coord_format)
        case "xsf":
            writefilexsf(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)
        case "fhi":
            writefilefhi(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)

        case "cif":
            writefilecif(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)

    print(f"[OK] Writing the file {args.out_file} ({args.out_format})")
    
    print("[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Converting input files is 10% coding, 90% crying.\n\n", 'bold'))

if __name__ == "__main__":
    main()
