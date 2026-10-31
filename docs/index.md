# STB-SUITE — Siesta Toolbox Suite

**A unified command-line toolkit for SIESTA DFT workflows.**

!!! warning dev-notice "Stable version 2.0.0, under testing"

    **This is the first stable release, and it is still being tested in real projects.**

    **Verify every result** before relying on it: check it independently, for
    example against a reference calculation or the literature.

[:material-file-pdf-box: Download the manual (PDF)](assets/stb_manual.pdf){ .md-button .md-button--primary }
&nbsp; every tool, workflow and option in one document, with the guides of this site.

STB-SUITE is a collection of independent command-line tools that assist users of
the [SIESTA](https://siesta-project.org/siesta/) DFT code through every step of a
project: generating inputs, building and converting structures, post-processing
results (bands, DOS, symmetry, charges, spectroscopy), and running complete
prepare-and-analyze workflows for a specific property (elastic constants,
phonons, adsorption, NEB barriers, ...).

Every tool is a `stb-*` command, and `stb-suite` is an interactive menu that wraps
all of them.

## Install

From source:

```bash
git clone https://github.com/stb-suite/stb.git
cd stb/stb-suite
pip install .
```

Python 3.9 to 3.12. The ML Simulations tools additionally need PyTorch and `mace-torch`,
which the optional `ml` extra installs:

```bash
pip install ".[ml]"                 # from stb/stb-suite, instead of the plain install
pip install "stb_suite[ml]"         # or, to add it to an existing install
```

There is no PyPI package: install from source as above.

## Recommended input layout

!!! tip "Three files per calculation (recommended, not required)"

    The suite works best when a SIESTA input is split in three files:

    | File | Holds | Written by |
    |------|-------|------------|
    | `structure.fdf` | the geometry only: atoms, species, cell (`%block LatticeVectors`, `%block AtomicCoordinatesAndAtomicSpecies`, ...) | you, or any structure tool (`stb-translate`, `stb-supercell`, `stb-slab`, ...) |
    | `calc.fdf` | the level of theory: `SystemLabel`, basis, XC functional, mesh cutoff, k-grid, SCF and MD settings, and the line `%include structure.fdf` | you, or `stb-inputfile` (1.1) |
    | `config_extra.fdf` | what a workflow stage must impose (a single point, a fixed cell, spin, dipole correction, D3, a k-grid matched to a supercell, an optical block, ...) | the workflow, never you |

    A workflow copies your `calc.fdf` unchanged into every folder it writes, rewrites only
    `structure.fdf`, and puts its own settings in `config_extra.fdf`, included at the **top** of
    the copy. SIESTA keeps the **first** value it reads for any label, so the forced settings
    win while your template stays exactly as you wrote it, and a later stage can read
    `config_extra.fdf` back to know what was really applied.

    Not mandatory: most tools also accept one self-contained `.fdf`. But with the split, the
    level of theory is defined once, by you, and every folder of a study shares it.

## Two ways to run any tool

**Directly, with flags** — scriptable, and what you use once you know a tool's options:

```bash
stb-kgrid -f structure.fdf -d 0.2
```

**Through the interactive menu** — guided questions instead of flags, a good way to
discover what a tool can do and what its defaults are:

```bash
stb-suite
```

At the menu's main prompt, either browse the categories or type a tool's code
(for example `1.3`) to jump straight to it. Both paths run the exact same underlying
tool and produce the same output.

## What's in the suite

The menu groups every tool into seven categories:

| # | Category | What it's for |
|---|----------|---------------|
| 1 | Inputs | Set up a SIESTA run: input file, pseudopotentials, k-grid, k-path |
| 2 | Analysis | Analyze simulation results: bands, DOS, structure, symmetry, charge density, charges |
| 3 | Utils | Helpers for file management and format conversion |
| 4 | Structures | Build, generate, or transform structure files: stacking, supercells, slabs, defects, SQS |
| 5 | Workflow | Complete prepare + analyze pipelines for one property: strain, elastic constants, phonons, adsorption, NEB |
| 6 | Catalysis | Catalysis workflows (HER, OER), provided by an extension package |
| 7 | ML Simulations | Train a MACE force field on SIESTA data, and run simulations with it or with the foundation model |

## Where to go next

The [**Guides**](guides/index.md) are hands-on tutorials, one per tool (or one per
workflow, for the multi-stage ones). Each explains what the tool does, the theory
behind it, and walks through a runnable example with both the direct command and the
menu path. Tools that do not have a written guide yet still get a short page listing
their commands.

The [**Reference**](reference/index.md) has one page per `stb-*` command, with its
full `--help` text, and a [map of every menu code](reference/stb-suite.md) for the
interactive menu.

## Source and license

The source code, the example inputs behind every guide, and the issue tracker are on
[GitHub](https://github.com/stb-suite/stb). Distributed under the MIT License.
Developed by Dr. Carlos M. O. Bastos, University of Brasília (UnB) —
[bastoscmo.github.io](https://bastoscmo.github.io).
