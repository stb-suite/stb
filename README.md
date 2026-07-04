# 🧰 STB-SUITE — Siesta Toolbox Suite

**A unified command-line toolkit for SIESTA DFT workflows**

![Version](https://img.shields.io/badge/version-1.9.1-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)
![Compatibility](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)

**Author:** Dr. Carlos M. O. Bastos — University of Brasília (UnB), 2025
🔗 [bastoscmo.github.io](https://bastoscmo.github.io)

---

## 📖 Overview

**STB-SUITE (Siesta Toolbox Suite)** is a comprehensive collection of command-line tools designed to assist users of the **SIESTA** DFT code through every step of their workflow — from input generation to post-processing and structural analysis.

It provides a unified, intuitive interface that simplifies and accelerates computational materials research.

---

## 🚀 Features

### 🧩 Calculation Preparation

* **`stb-inputfile`** – Generate FDF input files from structure files with suggested settings.
* **`stb-kgrid`** – Automatically suggest Monkhorst-Pack k-point grids based on a target density.
* **`stb-kpath`** – Generate high-symmetry paths for band-structure calculations.
* **`stb-strain`** – Create supercells with uniaxial or biaxial strain (Cartesian coordinates).
* **`stb-elasticInputs`** – Generate batches of deformed structures (strain tensor components) for elastic constant calculations.
* **`stb-cohesive`** – Set up cohesive energy workflows (bulk structure plus isolated-atom calculations).
* **`stb-phononsCreate`** – Generate displaced supercells for phonon calculations via Phonopy.
* **`stb-2Dstacking`** – Stack two monolayers into a heterostructure using the ZSL algorithm, with twist/shift control.

### 📊 Analysis & Post-Processing

* **`stb-bands`** – Analyze SIESTA band structures and calculate band gaps.
* **`stb-dos`** – Parse `PDOS.xml` files for total and projected density of states.
* **`stb-convdos`** – Apply Gaussian convolution to DOS data for smoothing.
* **`stb-structural`** – Compute lattice parameters, coordination numbers, and ECN values.
* **`stb-symmetry`** – Identify space group, point group, crystal system, and Wyckoff positions.
* **`stb-strainAnalysis`** – Extract stress-strain curves and mechanical properties (Young's modulus, UTS) from `strain_*` folders.
* **`stb-elasticAnalysis`** – Compute the stiffness matrix, elastic moduli, and Born stability criteria from `strain_*` folders.
* **`stb-cohesiveAnalysis`** – Calculate cohesive energy per atom from a completed cohesive-energy workflow.
* **`stb-bader`** – Perform Bader charge (AIM) analysis on SIESTA charge density grids.
* **`stb-cube`** – Convert SIESTA grid files (VT, VH, RHO) to Gaussian Cube format.
* **`stb-density`** – Export charge density as 2D slice maps or 3D point clouds for plotting.
* **`stb-workfunction`** – Calculate work function from the planar-averaged electrostatic potential.
* **`stb-phononsPos`** – Post-process phonon displacement calculations into thermal properties.

### ⚙️ Utilities & Interfaces

* **`stb-translate`** – Convert between structure formats (CIF, POSCAR, FDF, XYZ, XSF, FHI, DFTB).
* **`stb-siesta2wtb`** – Export SIESTA Hamiltonians to the **Wantibexos** tight-binding format.
* **`stb-clean`** – Remove unnecessary calculation files and clean directories.
* **`stb-suite`** – Unified terminal interface providing interactive access to all tools.

---

## 🧠 Requirements

* **Python ≥ 3.9**
* **Conda** (recommended)

### Python dependencies

```
numpy
matplotlib
ase
pymatgen
spglib
sisl
argparse (builtin)
```

---

## 📦 Installation

### 🔹 Conda (recommended)

```bash
conda install bastoscmo::stb_suite
```

### 🔹 Manual installation

```bash
git clone https://github.com/bastoscmo/stb-suite.git
cd stb-suite
pip install .
```

---

## ▶️ Usage

Launch the interactive interface:

```bash
stb_suite
```

Run individual tools directly:

```bash
stb-inputfile structure.fdf --type relax
stb-kgrid --file POSCAR --density 0.2
stb-symmetry --input struct.cif --filetype cif
```

---

## 🧾 License

Distributed under the **MIT License**.
© 2025 Dr. Carlos M. O. Bastos – University of Brasília (UnB)

---


