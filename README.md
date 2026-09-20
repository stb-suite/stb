# 🧰 STB-SUITE — Siesta Toolbox Suite

**A unified command-line toolkit for SIESTA DFT workflows**

![Version](https://img.shields.io/badge/version-1.9.1-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Python](https://img.shields.io/badge/python-3.9%20to%203.12-blue.svg)
![Compatibility](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)
[![Documentation](https://img.shields.io/badge/docs-stb--suite.github.io%2Fstb-informational.svg)](https://stb-suite.github.io/stb/)

**Author:** Dr. Carlos M. O. Bastos — University of Brasília (UnB), 2025
🔗 [bastoscmo.github.io](https://bastoscmo.github.io)

📚 **Documentation: <https://stb-suite.github.io/stb/>** — a guide for every tool and a reference page for every command.

---

## 📖 Overview

**STB-SUITE (Siesta Toolbox Suite)** is a collection of command-line tools that assist users of the **SIESTA** DFT code through every step of a project: generating inputs, building and converting structures, post-processing results, and running complete prepare-and-analyze workflows for a specific property.

Every tool is a standalone `stb-*` command, and `stb-suite` is an interactive menu that wraps all of them.

SIESTA itself is not bundled: the Inputs and Workflow tools write ready-to-run calculation folders, and the Analysis tools read what SIESTA produced.

---

## 🚀 What's inside

The menu groups the tools into six categories. Each tool has a guide with a runnable example in [`examples/`](examples/) (published on the [documentation site](https://stb-suite.github.io/stb/)).

| # | Category | What it's for | Tools |
|---|----------|---------------|-------|
| 1 | **Inputs** | Set up a SIESTA run | `stb-inputfile` `stb-pseudo` `stb-kgrid` `stb-kpath` `stb-dftu` `stb-fetch` `stb-mlrelax` |
| 2 | **Structures** | Build, generate or transform structure files | `stb-2Dstacking` `stb-supercell` `stb-slab` `stb-nanotube` `stb-defect` `stb-sqs` `stb-unitcell` `stb-crystalbuilder` `stb-passivate` `stb-molecule` `stb-amorphize` `stb-crystalcast` |
| 3 | **Analysis** | Analyze simulation results | `stb-bands` `stb-dos` `stb-convdos` `stb-structural` `stb-symmetry` `stb-bader` `stb-workfunction` `stb-density` `stb-xrd` `stb-fatbands` `stb-stm` `stb-wfdensity` `stb-sts` `stb-coop` `stb-ipr` `stb-effmass` `stb-spintexture` `stb-aimdAnalysis` |
| 4 | **Workflow** | Paired preparation + analysis pipelines for one property | Stress-strain, elastic constants, cohesive energy, phonons, convergence tests, structure solution (XRD), Hubbard U, adsorption, NEB, stacking fault, Raman, IR, HER, OER, GQCA, optical properties, custom ML force field, equation of state, charge density difference, Hirshfeld-I |
| 5 | **ML Simulations** | Run simulations with a MACE machine-learning potential instead of SIESTA | `stb-mlmd` `stb-mlphonons` `stb-mlelastic` `stb-mlsearch` `stb-mlmelting` `stb-mlconvergence` `stb-mlneb` `stb-mldiffusion` `stb-mlgcmc` `stb-mladsorb` `stb-mleos` |
| 6 | **Utils** | File management and format conversion | `stb-translate` `stb-clean` `stb-cube` `stb-siesta2wtb` `stb-ani2traj` `stb-status` `stb-archive` `stb-nativecharges` |

A workflow is several commands used in order (for example `stb-strain` then `stb-strainAnalysis`); the [documentation](https://stb-suite.github.io/stb/) lists every command of each one.

---

## 🧠 Requirements

* **Python 3.9 to 3.12**
* Installed automatically: `numpy`, `ase`, `matplotlib`, `pymatgen`, `spglib`, `sisl`, `pybader`, `scipy`, `phonopy`, `pandas`, `icet`, `pyxtal`
* Optional **`ml` extra** (`torch`, `mace-torch`, `torch-dftd`) for the ML tools — see below

---

## 📦 Installation

```bash
git clone https://github.com/stb-suite/stb.git
cd stb/stb-suite
pip install .
```

With the optional ML tools:

```bash
pip install ".[ml]"
```

If the suite is already installed, `pip install "stb_suite[ml]"` adds the ML extra to it.

There is no PyPI package: install from source as above.

---

## ▶️ Usage

Launch the interactive menu, then browse the categories or type a tool's code (for example `1.3`, or `4.1.2` for a workflow stage):

```bash
stb-suite
```

Or run any tool directly:

```bash
stb-inputfile structure.fdf --type relax
stb-kgrid -f structure.fdf -d 0.2
stb-symmetry --file structure.fdf --format fdf
```

Every tool prints its options with `--help` (for example `stb-kgrid --help`).

---

## 🗂️ Repository layout

* `stb-suite/` — the Python package (`src/stb/`) and its packaging metadata
* `examples/` — one tutorial per tool, with a runnable script and small input files
* `test/` — the developers' smoke tests
* `docs/`, `mkdocs.yml` — the documentation site, built from `examples/`

---

## 🧾 License

Distributed under the **MIT License**.
© 2025 Dr. Carlos M. O. Bastos – University of Brasília (UnB)
