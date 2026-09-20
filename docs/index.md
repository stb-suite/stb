# STB-SUITE — Siesta Toolbox Suite

**A unified command-line toolkit for SIESTA DFT workflows.**

STB-SUITE is a collection of independent command-line tools that assist users of
the [SIESTA](https://siesta-project.org/siesta/) DFT code through every step of a
project: generating inputs, building and converting structures, post-processing
results (bands, DOS, symmetry, charges, spectroscopy), and running complete
prepare-and-analyze workflows for a specific property (elastic constants,
phonons, adsorption, NEB barriers, ...).

Every tool is a `stb-*` command, and `stb-suite` is an interactive menu that wraps
all of them.

## Install

From conda:

```bash
conda install bastoscmo::stb_suite
```

From source:

```bash
git clone https://github.com/stb-suite/stb.git
cd stb/stb-suite
pip install .
```

The ML Simulations tools additionally need PyTorch and `mace-torch`:

```bash
pip install "stb_suite[ml]"
```

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

The menu groups every tool into six categories:

| # | Category | What it's for |
|---|----------|---------------|
| 1 | Inputs | Set up a SIESTA run: input file, pseudopotentials, k-grid, k-path |
| 2 | Structures | Build, generate, or transform structure files: stacking, supercells, slabs, defects, SQS |
| 3 | Analysis | Analyze simulation results: bands, DOS, structure, symmetry, charge density |
| 4 | Workflow | Complete prepare + analyze pipelines for one property: strain, elastic constants, phonons, adsorption, NEB |
| 5 | ML Simulations | Run simulations with a MACE machine-learning potential instead of SIESTA |
| 6 | Utils | Helpers for file management and format conversion |

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
