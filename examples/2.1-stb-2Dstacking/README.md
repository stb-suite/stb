# 2.1 — 2D Monolayer Stacker (`stb-2Dstacking`)

## What this tool does

Takes two 2D monolayer structures (each its own `.fdf`, periodic in-plane
with vacuum along `c`) and stacks them into a single van der Waals
heterostructure: a bottom layer and a top layer, separated by a chosen
interlayer gap, sharing one commensurate supercell. This is the first tool
in the **Structures** category (2) — everything in this category builds,
generates, or transforms a structure file, as opposed to **Inputs** (1),
which configures an actual SIESTA run.

## Why this matters (a bit of theory)

Two independent 2D crystals almost never share an exact lattice. Even a
"close" pair like graphene (a = 2.46 Å) and monolayer h-BN (a = 2.504 Å)
differ by ~1.8% — enough that neither lattice tiles the other exactly at
any size. Two chemically unrelated 2D materials can differ far more.

The **ZSL algorithm** (Zur & McGill, 1984 — see `references.bib` after any
run) solves this by searching, up to a maximum supercell area
(`--max_area`), over integer supercell transformation matrices for layer 1
and layer 2, looking for a pair whose resulting lattice vectors match each
other to within a maximum strain fraction (`--max_strain`). Every match is
a trade-off: forcing a match in a small cell usually needs more strain;
allowing a larger cell gives the search more candidates to find a
lower-strain match, at the cost of more atoms in the final structure. The
search itself (the expensive part) runs once per invocation; building
several stacking configurations from the same match (different gaps,
shifts, twists) is cheap and reuses it.

### Two independent match metrics — not just one "strain" number

Every candidate ZSL match is scored on two genuinely different quantities,
both shown in the `-i`/`--interactive` candidate table and in the tool's
own `[2]` report section:

- **Strain (%)** — the lattice-constant mismatch left after forcing the
  two supercells to the same size. This is what `--max_strain` bounds.
- **Ang. Strain (°)** — how much layer 2 has to be rigidly rotated
  ("un-rotated") for its supercell vectors to line up with layer 1's. ZSL
  matches lattice vectors by **length only**, never absolute orientation
  — a candidate can have ~0% length strain and still need a large
  rotation to actually fit. That rotation is a real geometric distortion
  of the input, not a bookkeeping detail, and it matters most once
  `--twist` is involved (see below). Matches with *both* numbers under
  ~1 are highlighted in green in the interactive table — a fast visual
  cue for "physically clean" candidates.

`match_id=0` (the default, and what `-id 0` picks) is always the
**lowest-strain** match, ties broken by lowest angular strain — it does
**not** consider the two symmetrically. This is almost always what you
want (see `basic/` and `batch-sym/` below), with one important exception
covered in the `twist/` case: a self-twisted homobilayer generates a whole
family of exact ~0%-strain "trivial" matches regardless of angular strain,
which can bury the genuine low-angular-strain Moiré cell far down the
candidate list.

### 21.8° vs. the "magic angle" — a twist-angle caveat, verified live

`--twist` rotates layer 2 by a fixed angle before the ZSL search runs —
the same geometric idea behind twisted-bilayer-graphene research. The demo
below uses 21.8°, which is not an arbitrary choice: for a hexagonal
lattice it is an **exact** commensurate angle, `cos(θ) = (3p² + 3pq +
q²/2) / (3p² + 3pq + q²)` with the smallest nontrivial case `p = q = 1` —
a small coincidence cell exists at this angle by construction. This is
deliberately **not** the real ~1.1° "magic angle" of flat-band twisted
bilayer graphene physics (Bistritzer–MacDonald): that one needs a
supercell of several thousand atoms, far beyond what a quick CLI tutorial
should build. 21.8° is chosen here specifically because it stays small
enough to run in seconds while still being a genuine commensurate twist.

A real, verified gotcha follows directly from the previous point: because
`-id 0` sorts on strain first, and a self-twisted homobilayer's trivial
"un-rotate the whole twist away" matches also have ~0% strain, `-id 0`
at 21.8° silently returns an **untwisted** structure — verified live: its
symmetry table shows `p6/mmm`, identical to the plain input, at every
`--max_area` from 50 up to 200. The genuine Moiré cell is in the candidate
list (verified: match ID 70 of 233, at `--max_area 50`) — it's just sorted
far below where the interactive table's 15-row preview stops, since its
strain is ~0% but not *exactly* the same floating-point 0.0 the trivial
matches get. Once picked directly (`-id 70 -a 50`), the result is
genuinely twisted: symmetry drops to triclinic `p1`, a real, physically
distinct structure from the naive pick. `example_2.1.sh`'s `twist/` case
shows both side by side.

## When you'd reach for it

- Building a realistic van der Waals heterostructure (e.g. graphene/h-BN,
  a 2D heterojunction) as a starting point for a real SIESTA calculation.
- Screening several high-symmetry registries (`--batch_sym`) of the same
  pair before committing to expensive DFT on all of them.
- Generating a set of structures at different interlayer gaps
  (`--gap_range`) as the first step toward an interlayer binding-energy
  curve.
- Building a twisted bilayer (`--twist`) of two identical (or different)
  monolayers, e.g. for Moiré-superlattice physics.
- Sampling an arbitrary registry point (`-tx`/`-ty`, not one of
  `--batch_sym`'s 6 fixed ones) — e.g. one point of a full gamma-surface
  (stacking-fault energy) scan.

## Two ways to run it

A — direct CLI:
```bash
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.1
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.1.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[6] SUMMARY & FILES` — `[3]`
  only appears with `--ml-relax`, see below) printed to the console and,
  with `--save-report`, also saved to `stb_2dstacking_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density, each an explicit `[OK]`/`[WARNING]`/`[SKIPPED]` row — run for
  **every** heterostructure the invocation builds (not just one), even in
  `--batch_sym`/`--gap_range` batch mode.
- **A symmetry comparison table** for layer 1 / layer 2 / the
  heterostructure — crystal system, the 3D **space group** of the
  vacuum-padded cell, the **layer group** (the physically correct symmetry
  classification for a genuinely 2D-periodic structure — same
  `core.symmetry.layer_group_label` stb-fetch already uses), point group,
  and Hall symbol. Folded directly into the main report (no separate file
  — there used to be a `--sym_out symmetry_report.txt`, removed once this
  table moved into the report itself).
- **`references.bib`** — SIESTA (the output is a `.fdf`) plus the ZSL
  algorithm paper. Merges with whatever `references.bib` is already in the
  working directory instead of overwriting it.
- **A provenance header** written into the output `.fdf` itself: which two
  layer files were used, the stacking configuration name, twist angle,
  strain mode, the max applied linear strain for that structure, and (with
  `--ml-relax`) the MACE model used plus its convergence/energy.

## Optional (off by default)

- **`-i`/`--interactive`** — instead of silently taking the lowest-strain
  match (`-id 0`), prints the top 15 candidates (area, strain, angular
  strain, transformation matrix) and lets you pick by hand. Rows where
  *both* strain and angular strain are under ~1 print in green — a quick
  visual shortcut for "this one's physically clean" without reading every
  column. Mutually exclusive with `-id`.
- **`-tx`/`-ty` (`--shift_x`/`--shift_y`)** — place layer 2 at any
  fractional shift instead of one of `--batch_sym`'s 6 fixed registries
  (ignored if `--batch_sym` is also given). Useful for sampling a full
  gamma-surface (stacking-fault energy landscape) point by point. An
  arbitrary shift usually lands on **no** special symmetry point (`p1`,
  no symmetry) — see `custom-shift/` below for why `--batch_sym`'s 6
  points specifically are the ones a hexagonal bilayer's own symmetry can
  protect.
- **`--vacuum`** — overrides the vacuum thickness instead of inheriting
  layer 1's own.
- **`--save-report`** — also persists the full report (structure
  validation, the symmetry table, and, with `--ml-relax`, the full MACE
  simulation detail) to `stb_2dstacking_report.txt`.
- **`--view`** — opens every generated heterostructure in ASE's
  interactive 3D viewer (`ase-gui`), as separate, pageable frames. Needs a
  local display (or `ssh -X`/`-Y`). Never exercised by `example_2.1.sh`
  itself.
- **`--ml-relax`** (needs the optional `ml` extra: `pip install
  stb_suite[ml]`) — pre-relaxes each generated heterostructure with a MACE
  potential before writing it out. Positions only by default; add
  **`--ml-relax-cell`** to also relax the in-plane cell (`a`/`b`) — the
  vacuum axis (`c`) always stays exactly fixed, via the same
  `core.mace_relax.build_cell_mask` masking `stb-mlrelax` uses for a slab.
  `--model small/medium/large` (default `small`) or `--custom-model PATH`
  picks which MACE potential to use. The report shows the full simulation
  detail per structure: model parameter count/cutoff radius, steps used,
  convergence, wall time, and a **before/after table** (same convention as
  `1.6-stb-mlrelax/`'s own comparison table) — energy, max force, and, with
  `--ml-relax-cell`, lattice `a`/`b`/`gamma` alongside proof that `c` (the
  vacuum axis) never moves. A fast heuristic, not a substitute for a real
  DFT relaxation.

## Files in this folder

- `graphene.fdf` / `hbn.fdf` — a primitive 2-atom hexagonal graphene cell
  (a = 2.46 Å) and a primitive 2-atom hexagonal h-BN cell (a = 2.504 Å),
  a ~1.8% lattice-mismatched pair — the classic van der Waals
  heterostructure test case, also used for the `--twist` demo (with itself
  as both layers).
- `example_2.1.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_2dstacking_report.txt`.

## Running the walkthrough

```bash
cd examples/2.1-stb-2Dstacking
./example_2.1.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Eight self-contained cases are generated (the
last one is skipped, with a hint, if the optional `ml` extra isn't
installed):

| Folder              | What it shows                                             |
|---------------------|------------------------------------------------------------|
| `basic/`             | The mismatched graphene + h-BN pair, `-id 0`                |
| `batch-sym/`         | `--batch_sym` — all 6 high-symmetry registries in one run   |
| `gap-range/`         | `--gap_range` — the same registry at several interlayer gaps|
| `strain-modes/`      | `-sm top/bottom/sym` side by side                            |
| `twist/`             | `--twist` — a genuine 21.8° Moiré cell vs. the naive, silently-untwisted `-id 0` pick|
| `custom-shift/`      | `-tx`/`-ty` + `--vacuum` — an arbitrary registry point, contrasted with `--batch_sym`'s symmetry-protected ones|
| `full-report/`       | `--save-report` + the validation checklist + `references.bib`|
| `ml-relax/`          | `--ml-relax`/`--ml-relax-cell` — MACE pre-relaxation, vacuum axis fixed|

## Try it yourself

```bash
# Pick a match by hand instead of -id 0 (shows the full candidate table)
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -i

# Open the result in ASE's interactive viewer
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --view

# Combine batch stacking with a gap scan
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --batch_sym --gap_range 3.0 4.0 6

# Pre-relax with MACE before writing, cell shape included (vacuum axis stays fixed)
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --ml-relax --ml-relax-cell

# The genuine 21.8 deg Moire cell (match ID 0 alone silently un-rotates this away)
stb-2Dstacking -l1 graphene.fdf -l2 graphene.fdf -id 70 -a 50 -t 21.8
```

## Flag reference

| Flag                          | Meaning                                                          |
|-------------------------------|-------------------------------------------------------------------|
| `-l1/--layer1`, `-l2/--layer2`| Bottom / top monolayer `.fdf` (both required)                     |
| `-a/--max_area`               | Max ZSL supercell area, Å² (default 150.0)                        |
| `-s/--max_strain`             | Max ZSL strain fraction (default 0.05)                            |
| `-i/--interactive`            | Show the full match table and pick by hand                        |
| `-id/--match_id`              | Pick a specific match ID directly (default 0)                     |
| `--batch_sym`                 | Build all high-symmetry registries automatically                  |
| `-g/--gap`                    | Van der Waals gap in Å (default 3.2)                              |
| `--gap_range START END STEPS` | Build the same registry at several gaps (energy-curve mode)        |
| `--vacuum`                    | Target vacuum thickness, Å (default: inherits layer 1's)           |
| `-t/--twist`                  | Twist angle of layer 2, degrees, before the ZSL search             |
| `-tx/--shift_x`, `-ty/--shift_y` | Fractional shift of layer 2 (ignored with `--batch_sym`)        |
| `-sm/--strain_mode`           | `top` (default) / `bottom` / `sym` — which layer absorbs strain    |
| `-o/--output`                 | Output `.fdf` base name (default `stacked_structure.fdf`)          |
| `-sp/--symprec`               | Symmetry tolerance, Å (default 0.01)                               |
| `--save-report`               | Persist the full report (incl. symmetry table) to disk            |
| `--view`                      | Open every result in ASE's interactive 3D viewer                   |
| `--ml-relax`                  | Pre-relax with a MACE potential before writing (positions only)    |
| `--ml-relax-cell`             | With `--ml-relax`, also relax `a`/`b` (vacuum axis `c` stays fixed) |
| `--model`                     | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`)|
| `--custom-model`              | Path to a custom fine-tuned `.model` file for `--ml-relax`         |

## What's next

See `1.1-stb-inputfile/` for turning a stacked heterostructure into a full
SIESTA run (basis, pseudopotentials, k-grid), or `1.6-stb-mlrelax/` for a
fast MACE-based pre-relaxation before committing to real DFT.
