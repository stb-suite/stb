# 2.2 — Supercell Builder (`stb-supercell`)

## What this tool does

Takes one structure file and tiles it into a larger periodic cell, following
an integer transformation matrix: either 3 numbers for a diagonal supercell
(`n_a n_b n_c`, the common case — repeat the cell `n_a`/`n_b`/`n_c` times
along each lattice vector) or 9 numbers for a full row-major 3×3 matrix
(needed for non-orthogonal transformations, see `reconstruction/` below).
Same **Structures** category (2) as `2.1-stb-2Dstacking/` — everything here
builds, generates, or transforms a structure file, as opposed to **Inputs**
(1), which configures an actual SIESTA run.

## Why this matters (a bit of theory)

### The transformation matrix

A supercell is described by a 3×3 integer matrix **M**: the new lattice
vectors are `M` applied to the original ones (`pymatgen`'s own
`Structure.make_supercell`). Two matrix shapes cover almost every real use:

- **Diagonal**, `diag(n_a, n_b, n_c)` — the everyday case: repeat the cell
  `n_a`×`n_b`×`n_c` times. The **determinant** (`n_a·n_b·n_c` here) is the
  atom-count multiplication factor — `stb-supercell` always reports it
  explicitly next to the matrix.
- **Full (non-diagonal)** — needed whenever the enlarged cell isn't simply
  "the same shape, bigger": a hexagonal 2D material's classic
  √3×√3 R30° reconstruction cell (`reconstruction/` below) needs off
  -diagonal entries, because the new lattice vectors point along different
  crystallographic directions than the original ones, not just longer
  copies of them.

A **negative determinant** is still geometrically valid — `mirrored/` below
shows it — but it is worth pausing on: mirroring a cell that has no
inversion/mirror symmetry of its own (a *chiral* crystal) produces a
genuinely different, enantiomeric structure, not the same one relabeled.
For a centrosymmetric crystal like diamond silicon (which already contains
inversion symmetry) the mirrored cell is physically indistinguishable from
the original — confirmed below, its full symmetry table comes back
identical, entry for entry.

### Why build a supercell at all

A few recurring reasons, each pointing at a different downstream tool in
this suite:

- **Diluting a point defect** (`stb-defect`) — a defect's periodic images
  interact with each other through the periodic boundary; a big enough
  supercell keeps that spurious self-interaction small and sets the
  effective defect *concentration* (`1 / N_atoms`).
- **Phonon force constants** (`stb-phononsCreate`/`stb-mlphonons`) — the
  finite-displacement method needs a real-space supercell large enough that
  a single atom's displacement has decayed to ~0 force on the supercell's
  own periodic image of itself.
- **Non-orthogonal reconstructions and Moiré-like commensurate cells** —
  the full 3×3 matrix path, see `reconstruction/`.
- **Matching a target k-point density** cheaply in real space instead of a
  denser Monkhorst-Pack grid, when a workflow specifically needs a bigger
  real-space cell (e.g. before `stb-sqs` places atoms across it for a
  disordered-alloy model).

### A genuine gotcha, verified live: Space Group survives, Point Group can look like it doesn't

A well-formed supercell doesn't create a new crystal, it re-describes the
same one in a bigger cell — so you'd expect the *whole* symmetry table to
come back unchanged. The **Space Group**, **Layer Group**, and **Hall
Symbol** columns do exactly that (verified below on both `basic/` and
`reconstruction/`): they're derived from `spglib`'s own internally
standardized cell, which is independent of the shape you actually handed
it.

The **Point Group** column can mislead you, though. `pymatgen`'s
`SpacegroupAnalyzer.get_point_group_symbol()` classifies the point group
from the raw rotation matrices found for the *exact* cell you gave it, in
that cell's own basis — not from the standardized space-group data. An
**anisotropic** diagonal supercell (e.g. doubling only along `c`, see
`nonuniform/`) genuinely restricts *which* of the crystal's true symmetry
operations still map that specific, now non-cubic-shaped cell onto itself
as an integer matrix. The physical crystal is still fully cubic; the
*given* cell just isn't shaped like a cube any more, so the reported point
group quietly drops from `m-3m` to `4/mmm`. Verified directly against
`spglib`'s own symmetry dataset: `len(dataset.rotations)` drops from 192 to
128 operations found for the exact same physical crystal, purely because of
the cell's shape — while `dataset.pointgroup` (a *different* spglib field,
computed from the standardized cell) still correctly reports `m-3m`.

**Takeaway**: after building an anisotropic supercell, trust the **Space
Group** row for "is this still the same crystal", not the **Point Group**
row — the latter reports the point symmetry of the specific cell shape you
built, which a non-uniform tiling can legitimately (not buggily) reduce.

## When you'd reach for it

- Preparing a supercell for `stb-defect` (a vacancy/substitution needs
  periodic images far enough apart) or `stb-passivate`/`stb-slab`.
- Setting up the finite-displacement supercell a phonon calculation needs.
- Building a 2D material's reconstruction/adsorbate-superlattice cell
  (`reconstruction/`).
- Getting a bigger, MACE-pre-relaxed starting cell before a real SIESTA
  relaxation of a defect or disordered structure (`ml-relax/`).

## Two ways to run it

A — direct CLI:
```bash
stb-supercell -f si_bulk.fdf -d 2 2 2
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.2
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.2.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[9] SUMMARY & FILES` — `[4]`
  only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_supercell_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density, each an explicit `[OK]`/`[WARNING]`/`[SKIPPED]` row — run once on
  the input structure and once on the final (possibly ML-relaxed)
  supercell, so a mirrored cell's left-handedness or a too-small density is
  never silently missed.
- **A before/after symmetry comparison table** — crystal system, the 3D
  space group, the layer group (for a 2D-periodic structure, same
  `core.symmetry.layer_group_label` `stb-fetch`/`stb-2dstacking` already
  use), point group, and Hall symbol, for the input structure vs. the final
  supercell.
- **`references.bib`** — SIESTA always (the output is a `.fdf`), plus the
  MACE papers if `--ml-relax` was used. Merges with whatever
  `references.bib` is already in the working directory.
- **A provenance header** written into the output `.fdf` itself: the input
  file, the transformation matrix and its determinant, the atom-count
  change, and (with `--ml-relax`) the MACE model used plus its convergence
  and energy change.

## Optional (off by default)

- **`-sp`/`--symprec`** — symmetry tolerance (Å) for the before/after table
  (default `0.01`).
- **`--ml-relax`** (needs the optional `ml` extra: `pip install
  stb_suite[ml]`) — pre-relaxes the *built* supercell with a MACE potential
  before writing it out. Positions only by default; add
  **`--ml-relax-cell`** to also relax the cell — any vacuum-padded axis
  always stays exactly fixed, same `core.mace_relax.build_cell_mask` masking
  `stb-mlrelax`/`stb-2dstacking` use. `--model small/medium/large` (default
  `small`) or `--custom-model PATH` picks which MACE potential to use. A
  fast heuristic, not a substitute for a real DFT relaxation.
- **`--save-report`** — also persists the full report (validation, the
  symmetry table, and, with `--ml-relax`, the full MACE simulation detail)
  to `stb_supercell_report.txt`.
- **`--view`** — opens the input structure and the final supercell side by
  side in ASE's interactive 3D viewer (`ase-gui`), as pageable frames. Needs
  a local display (or `ssh -X`/`-Y`). Never exercised by `example_2.2.sh`
  itself.

## Files in this folder

- `si_bulk.fdf` — bulk silicon, diamond cubic, the real conventional 8-atom
  cell (`a` = 5.431 Å, space group Fd-3m No. 227) — a clean, highly
  symmetric structure so the before/after symmetry table's lesson is
  unambiguous.
- `graphene.fdf` — a primitive 2-atom hexagonal graphene monolayer
  (`a` = 2.46 Å, vacuum along `c`) for the non-orthogonal reconstruction
  case.
- `example_2.2.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_supercell_report.txt`.

## Running the walkthrough

```bash
cd examples/2.2-stb-supercell
./example_2.2.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Six self-contained cases are generated (the
last one is skipped, with a hint, if the optional `ml` extra isn't
installed):

| Folder              | What it shows                                                        |
|----------------------|-----------------------------------------------------------------------|
| `basic/`             | Diagonal 2×2×2 supercell of bulk Si — the everyday case, symmetry table unchanged entry for entry |
| `nonuniform/`        | Diagonal but anisotropic 1×1×2 — the Space-Group-survives-but-Point-Group-drops gotcha, verified live |
| `reconstruction/`    | Full row-major 3×3 matrix — graphene's √3×√3 R30° reconstruction cell |
| `mirrored/`          | A negative-determinant (mirrored) cell — still valid, and physically identical here since Si is centrosymmetric |
| `ml-relax/`          | `--ml-relax`/`--ml-relax-cell` — MACE pre-relaxation of the built supercell |
| `full-report/`       | `--save-report` + the validation checklist + `references.bib`        |

## Try it yourself

```bash
# A general (non-diagonal) matrix, e.g. a different 2D reconstruction
stb-supercell -f graphene.fdf -d 2 0 0 0 2 0 0 0 1

# Open the result in ASE's interactive viewer, input vs. supercell
stb-supercell -f si_bulk.fdf -d 2 2 2 --view

# Pre-relax with MACE, cell included, with your own fine-tuned model
stb-supercell -f si_bulk.fdf -d 2 2 2 --ml-relax --ml-relax-cell --custom-model my_finetuned.model

# A tighter symmetry tolerance for a noisier input structure
stb-supercell -f si_bulk.fdf -d 2 2 2 -sp 0.1
```

## Flag reference

| Flag                | Meaning                                                                |
|----------------------|-------------------------------------------------------------------------|
| `-f/--file`          | Input structure file, `.fdf` (required)                                |
| `-d/--dim`           | 3 numbers (diagonal) or 9 numbers (full row-major matrix) (required)   |
| `-o/--output`        | Output `.fdf` file name (default `supercell.fdf`)                      |
| `-sp/--symprec`      | Symmetry tolerance, Å, for the before/after table (default `0.01`)     |
| `--ml-relax`         | Pre-relax the built supercell with a MACE potential before writing     |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the cell (vacuum-padded axes stay fixed) |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file for `--ml-relax`             |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the input structure and final supercell in ASE's interactive viewer |

## What's next

See `2.1-stb-2Dstacking/` for a different kind of cell-building tool (two
independent monolayers matched via the ZSL algorithm instead of one
structure tiled by an integer matrix), or `1.6-stb-mlrelax/` for a closer
look at the MACE pre-relaxation `--ml-relax` reuses here.
