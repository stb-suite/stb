# 2.7 — Unit Cell Finder (`stb-unitcell`)

## What this tool does

Detects a structure's true crystal symmetry (via pymatgen's symmetry
analyzer, spglib) and rebuilds it as the **primitive** (smallest possible),
**conventional** (standardized, usually larger), or **refined**
(conventional-sized, positions snapped exactly onto the detected symmetry)
cell. Same **Structures** category (2) as `2.1-stb-2Dstacking/`,
`2.2-stb-supercell/`, ..., `2.6-stb-sqs/` — everything here builds,
generates, or transforms a structure file. Unlike most of its siblings,
this tool doesn't change the crystal at all — it only changes *which*
repeating cell describes it.

## Why this matters (a bit of theory)

### One crystal, several valid cells

A crystal's infinite periodicity can be described by more than one
repeating cell, and all of them describe the *exact same* structure —
they're just different choices of "the box that tiles to make it":

- **Primitive** — the smallest possible repeating unit. Fewer atoms means
  a cheaper DFT calculation: fewer basis functions, a smaller
  Hamiltonian, faster SCF — the natural choice for a real calculation,
  especially phonon or band-structure work where the cost scales with
  atom count per k-point.
- **Conventional** — a larger, standardized cell following the
  International Tables / IUCr convention. Easier for a human to read and
  compare against a database entry (most databases, e.g. Materials
  Project or ICSD, store this one) — at the cost of more atoms than
  strictly necessary for a calculation.
- **Refined** — conventional-sized, but with every atomic position
  snapped exactly onto the detected symmetry. Doesn't change which atoms
  are present or the cell size — it only cleans up numerical noise (e.g.
  from a DFT relaxation, or a hand-typed/CIF-derived structure) that
  keeps a structure just barely off its true, exact symmetry.

### A real reduction, and a real surprise: cell shape ≠ crystal system

`primitive-fcc/` below reduces FCC nickel's 4-atom conventional cubic
cell down to its true 1-atom primitive cell — and that primitive cell's
own lattice vectors turn out to have 60° angles between them (a
rhombohedral-*looking* cell), even though the crystal **system** is still
reported as cubic. This isn't a contradiction: the crystal system comes
from the full set of symmetry operations the space group has (cubic
Fm-3m, No. 225 either way), not from how any *one* particular choice of
repeating cell happens to look geometrically. `primitive-rocksalt/`
verifies the same reduction preserves a compound's exact stoichiometric
ratio too (NaCl's 4:4 Na:Cl → 1:1), not just a single-species metal's
atom count.

### Reduction is symmetry-*preserving* — the tool's own before/after table proves it

Unlike `2.2-stb-supercell/`, `2.5-stb-defect/`, or `2.6-stb-sqs/` (which
all deliberately *change* a structure's symmetry), `stb-unitcell`'s job
is the opposite: find a *different cell for the same symmetry*. Every run
now prints the identical `[7] SYMMETRY ANALYSIS (BEFORE / AFTER)` table
those other tools use — but here it's a **correctness check**: Before and
After are *expected* to match exactly, and the report says so explicitly
when they do (verified live in every case below).

### `--symprec`: a real, dramatic consequence, not just a numeric knob

`--symprec` sets the numerical tolerance spglib uses to decide whether
two positions count as symmetry-equivalent. `symprec-sensitivity/` below
runs the exact same (deliberately slightly noisy) structure through two
tolerances: the default (`1e-3`, looser than the ~2×10⁻⁵ Å noise) sees
straight through it to the true Fm-3m symmetry and reduces the cell; a
tolerance *tighter* than the noise itself (`1e-8`) sees only genuine,
literal asymmetry — space group P1, no symmetry at all — and correctly
refuses to reduce anything, since there's nothing it can prove is
redundant at that tolerance. Same atoms, same physical structure, two
completely different answers, purely from the tolerance.

**A real bug found and fixed while building this example**: the
structure-validation sections' own space-group line used to ignore
`--symprec` entirely (always defaulting to `1e-3` internally), while the
symmetry-detection and before/after sections correctly used whatever
`--symprec` was requested. At a tight `--symprec`, this meant the *same
run* could show `Fm-3m` in one section and `P1` in another — confusing,
and specifically misleading for this tool, since `--symprec` is not an
auxiliary setting here (as it is for `2.2`/`2.5`/`2.6`'s own before/after
tables) but the primary parameter controlling what "the primitive cell"
even means. Fixed by threading `--symprec` through the validation
sections too; verified live that all sections now agree.

### The vacuum-axis and origin/relabeling caveats

This tool operates on the literal 3D periodic cell as given — it has no
special vacuum-axis awareness the way `stb-slab`/`stb-supercell` do, so a
2D/slab-type input gets an explicit `[WARNING]` (`vacuum-warning/`
below). The vacuum thickness itself is preserved, but — as documented for
spglib-based tools elsewhere in this suite (`core/symmetry.py`) — the
in-plane lattice vectors and atom ordering are **not** guaranteed to
match the input in any mode: spglib rebuilds the cell from the detected
symmetry operations from scratch, free to pick any symmetry-equivalent
origin/ordering. `vacuum-warning/` shows this directly: the output's
in-plane vectors come back re-expressed differently from the input, even
though it's the exact same physical lattice.

## When you'd reach for it

- Switching to the primitive cell before a real DFT calculation, to cut
  the atom count (and cost) without changing the physics.
- Getting the standardized conventional cell for comparison against a
  database entry or a paper.
- Cleaning up a relaxed or hand-built/CIF-derived structure's numerical
  noise (`--mode refined`) before a symmetry-sensitive analysis (e.g.
  group-theory-based phonon selection rules).

## Two ways to run it

A — direct CLI:
```bash
stb-unitcell -f fcc_ni.fdf --mode primitive
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.7
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.7.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[10] SUMMARY & FILES` —
  `[5]` only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_unitcell_report.txt`. Earlier
  versions of this tool always wrote a separate `unitcell_report.txt`;
  that file is gone — everything it used to contain is now part of this
  same numbered report.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density — run once on the input structure and once on the reduced
  cell, both using the same `--symprec` as the reduction itself.
- **A before/after symmetry comparison table** — crystal system, 3D space
  group, layer group, point group, Hall symbol — expected (and verified)
  to match exactly, since cell reduction doesn't change the crystal.
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax` was used.
- **A provenance header** written into the output `.fdf`: the mode used,
  the detected space group, the reduction factor, and MACE
  convergence/energy detail if used.

## Optional (off by default)

- **`--mode {primitive,conventional,refined}`** — which cell to build
  (default `primitive`).
- **`--symprec`** / **`--angle-tolerance`** — the symmetry-detection
  tolerance, used both for the reduction itself and the before/after
  table.
- **`--ml-relax`** (+ `--ml-relax-cell`, `--model`/`--custom-model`) —
  pre-relax the reduced cell generically with MACE.
- **`--save-report`** — also persist the full report to
  `stb_unitcell_report.txt`.
- **`--view`** — opens the input structure and the reduced cell in ASE's
  interactive 3D viewer (`ase-gui`), as pageable frames. Needs a local
  display. Never exercised by `example_2.7.sh` itself.

## Files in this folder

- `fcc_ni.fdf` — bulk FCC nickel, the conventional 4-atom cell (`a` =
  3.52 Å, space group Fm-3m No. 225) — the classic textbook
  conventional-vs-primitive case.
- `nacl_rocksalt.fdf` — NaCl, rocksalt structure, conventional 8-atom
  cell (`a` = 5.64 Å, also Fm-3m No. 225 — a deliberate parallel with
  `fcc_ni.fdf`) — two species, used to verify the reduction preserves
  their exact ratio.
- `nacl_noisy.fdf` — the same NaCl crystal with a ~2×10⁻⁵ Å asymmetric
  perturbation on every position, for the `refined-noise/` and
  `symprec-sensitivity/` cases.
- `graphene.fdf` — 2D monolayer graphene (vacuum along `c`), for the
  `vacuum-warning/` case.
- `example_2.7.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_unitcell_report.txt`.

## Running the walkthrough

```bash
cd examples/2.7-stb-unitcell
./example_2.7.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Eight self-contained cases are generated (one
is skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                    | What it shows                                                        |
|-----------------------------|--------------------------------------------------------------------|
| `primitive-fcc/`            | FCC Ni's 4→1 reduction, and its 60° rhombohedral primitive cell    |
| `primitive-rocksalt/`       | NaCl's 8→2 reduction preserving the exact 1:1 species ratio        |
| `conventional-noop/`        | A legitimate no-op, clearly flagged instead of silently doing nothing |
| `refined-noise/`            | Noisy positions snapped back to their exact symmetry-consistent values |
| `symprec-sensitivity/`      | The same noisy structure at two tolerances: Fm-3m vs. P1           |
| `vacuum-warning/`           | The vacuum-axis caveat and the origin/relabeling caveat on graphene |
| `ml-relax/`                 | Generic MACE pre-relaxation of the reduced cell                    |
| `full-report/`              | `--save-report` + the validation checklist + `references.bib`      |

## Try it yourself

```bash
# Show the input structure and the reduced cell side by side
stb-unitcell -f fcc_ni.fdf --mode primitive --view

# The standardized conventional cell, for comparison against a database entry
stb-unitcell -f nacl_rocksalt.fdf --mode conventional

# Pre-relax with MACE, cell included, with your own fine-tuned model
stb-unitcell -f fcc_ni.fdf --ml-relax --ml-relax-cell --custom-model my_finetuned.model

# Clean up a relaxed structure's numerical noise before a symmetry-sensitive analysis
stb-unitcell -f your_relaxed_structure.fdf --mode refined
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-f/--file`          | Input structure file, `.fdf` (required)                                |
| `--mode`             | `primitive` (default) / `conventional` / `refined`                     |
| `--symprec`          | Symmetry precision, default `1e-3` (used for both the reduction and the before/after table) |
| `--angle-tolerance`  | Symmetry angle tolerance in degrees, default `5.0`                     |
| `--ml-relax`         | Pre-relax the reduced cell with a MACE potential                       |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the cell (vacuum axes stay fixed)        |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file                              |
| `-o/--output`        | Output `.fdf` file name (default `unitcell.fdf`)                       |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the input and reduced-cell structure in ASE's interactive viewer  |

## What's next

See `2.2-stb-supercell/` for the opposite direction (building a larger
cell from a smaller one), or `1.6-stb-mlrelax/` for a closer look at the
MACE pre-relaxation `--ml-relax` reuses here.
