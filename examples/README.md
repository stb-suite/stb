# Examples

Hands-on, readable tutorials aimed at people **learning to use** stb_suite —
as opposed to `test/`, which is the developers' smoke-test suite (edge
cases, expected errors, PASS/FAIL counts). Each folder here is a small
lesson: what the tool does, why you'd reach for it, some of the theory
behind why it works the way it does, and a commented script you run and
read at the same time.

This is an ongoing effort — new examples are added a few at a time, one per
tool, not all at once. Not every `stb-*` tool has one yet.

## How the suite itself is organized

`stb-suite`, the interactive front-end, groups every tool into 6
categories:

| # | Category         | What it's for                                                                                        |
|---|------------------|-------------------------------------------------------------------------------------------------------|
| 1 | Inputs           | Tools to set up a SIESTA run (input file, k-grid, k-path)                                             |
| 2 | Structures       | Tools to build, generate, or transform structure files (stacking, supercells, slabs, defects, SQS...) |
| 3 | Analysis         | Tools to analyze simulation results (bands, DOS, structures, charge density...)                       |
| 4 | Workflow         | Complete prep + analysis pipelines for a specific property (strain, elastic constants, cohesive energy, phonons) |
| 5 | ML Simulations   | Run MD (and more) with a trained MACE potential — the foundation model or your own fine-tuned one     |
| 6 | Utils            | Helper tools for file management and conversion                                                       |

Each example folder here is named after the **dotted code** the interactive
`stb-suite` menu uses to jump straight to a tool — `<category>.<item>`, and
one level deeper for `Workflow`'s prep/analysis pairs (e.g. `4.1.2`). So
`1.1-stb-inputfile/` corresponds to typing `1.1` at the `stb-suite` main
prompt (category `1` = Inputs, item `1` = Input File Generator), or running
`stb-inputfile` directly from the command line — both paths are shown side
by side in every example.

To see the full, current list of codes yourself (this file's own Index
below only lists the tools that already have a written example), run
`stb-suite` and either browse the category menus, or read the
`INPUT_TOOLS`/`STRUCTURE_TOOLS`/`ANALYSIS_TOOLS`/`WORKFLOW_TOOLS`/
`MLSIM_TOOLS`/`UTILITY_TOOLS` dictionaries directly in
`stb-suite/src/stb/stb_suite.py` — that source is the single source of
truth for the code assigned to any tool.

## Two ways to run any tool

Every example shows both:

1. **Direct CLI** — call the `stb-*` command yourself with flags. Faster,
   scriptable, and what you'd use in a script or once you already know the
   tool's flags by heart.
2. **Interactive `stb-suite` menu** — run `stb-suite` and either navigate
   the category menus or type the tool's dotted code (e.g. `1.1`) straight
   from the main prompt. It walks you through the same choices as guided
   questions instead of flags — a good way to discover what a tool can do,
   and what its defaults are, before memorizing its options.

Both paths call the exact same underlying tool and produce the exact same
output — every example script proves this directly, by driving both paths
with the same inputs and diffing (or otherwise comparing) the result.

## What's inside an example folder

Every example follows the same shape:

- **`README.md`** — what the tool does, when you'd reach for it, and any
  relevant theory (increasingly common for the newer ML-driven tools:
  what the underlying method/model actually does, not just how to invoke
  it). Ends with a flag reference and pointers to what to try next.
- **`example_<code>.sh`** — the guided walkthrough script. Not an automated
  test (see the matching folder under `test/` for that); it runs real
  commands, one self-contained case at a time, into its own
  `output/<case>/` folder, and prints the specific piece of output that
  proves what just happened. Pauses between sections (`[Press Enter to
  continue]`) so you can read before moving on — safe to re-run any time,
  it always starts by wiping its own `output/`.
- **Small input fixtures** (`*.fdf`, occasionally other formats) — checked
  into git, purpose-built for that example (not shared with `test/`'s own
  fixtures, even when conceptually similar), usually with a comment inside
  explaining what's deliberately unusual about them (a wrong lattice
  constant, a too-close pair of atoms, a specific symmetry, ...).
- **`.gitignore`** — excludes `output/` and any files the tool itself
  generates on every run (`references.bib`, a `stb_<tool>_report.txt`,
  convergence plots, ...), so a fresh clone only ever has the tutorial
  itself, not its own regenerable output.

## Conventions shared across examples

Several conventions recur across the newer input-generating tools (see
`1.1-stb-inputfile/`, `1.4-stb-dftu/`, `1.5-stb-fetch/`, `1.6-stb-mlrelax/`,
`2.1-stb-2Dstacking/`, `2.2-stb-supercell/`, `2.3-stb-slab/`) and are worth
knowing once, rather than re-discovering per example:

- **A numbered, sectioned report** (`[0] RUN METADATA`, `[1] ...`, ...)
  printed to the console and, with `--save-report`, also saved to
  `stb_<tool>_report.txt` — the same report either way, `--save-report`
  just also persists it to disk.
- **`references.bib`** — written automatically, no flag needed, whenever a
  tool's output depends on a citable method (SIESTA itself, a specific
  correction/potential, a database it queried, ...). Merges with whatever
  `references.bib` is already in the same folder instead of overwriting
  it, so running several tools into one working directory accumulates one
  combined bibliography.
- **Structure validation with an explicit checklist** — atom proximity,
  lattice handedness, atomic density (and a couple of tool-specific extras)
  each get their own reported `[OK]`/`[WARNING]`/`[SKIPPED]` row, not just
  a pass/fail summary, so it's always clear exactly what was checked.
- **`--view`** — an optional, off-by-default flag that opens the structure
  in ASE's interactive 3D viewer right before the tool exits (needs a local
  display, or `ssh -X`/`-Y`). Never exercised by the example scripts
  themselves (an unattended walkthrough shouldn't pop a GUI window), but
  every README mentions the exact command to try it yourself.

## Prerequisites

Install the package in editable mode before running any example:

```bash
cd stb-suite
pip install -e .
```

This puts every `stb-*` command (and `stb-suite` itself) on your PATH.

A few examples (currently `1.6-stb-mlrelax/`, anything under the future
`5.x-*` ML Simulations examples) additionally need the optional `ml` extra
— PyTorch + `mace-torch`:

```bash
pip install stb_suite[ml]
```

Each such example's own script checks for this up front and exits cleanly
with an install hint if it isn't available, rather than failing partway
through.

## Index

| Code | Tool            | Folder                                    |
|------|-----------------|--------------------------------------------|
| 1.1  | `stb-inputfile` | [`1.1-stb-inputfile/`](1.1-stb-inputfile/)  |
| 1.2  | `stb-kgrid`     | [`1.2-stb-kgrid/`](1.2-stb-kgrid/)          |
| 1.3  | `stb-kpath`     | [`1.3-stb-kpath/`](1.3-stb-kpath/)          |
| 1.4  | `stb-dftu`      | [`1.4-stb-dftu/`](1.4-stb-dftu/)            |
| 1.5  | `stb-fetch`     | [`1.5-stb-fetch/`](1.5-stb-fetch/)          |
| 1.6  | `stb-mlrelax`   | [`1.6-stb-mlrelax/`](1.6-stb-mlrelax/)      |
| 2.1  | `stb-2Dstacking`| [`2.1-stb-2Dstacking/`](2.1-stb-2Dstacking/)|
| 2.2  | `stb-supercell` | [`2.2-stb-supercell/`](2.2-stb-supercell/)  |
| 2.3  | `stb-slab`      | [`2.3-stb-slab/`](2.3-stb-slab/)            |

## Adding another example

If you're adding the next one: pick any tool without a folder yet here,
name it `<code>-<tool-name>/` after its `stb-suite` dotted code, and follow
the shape described in "What's inside an example folder" above — the
existing folders (especially the most recently added one) are the best
template to copy from.
