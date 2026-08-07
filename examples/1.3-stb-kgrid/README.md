# 1.3 — K-Grid Generator (`stb-kgrid`)

## What this tool does

SIESTA needs to know how finely to sample reciprocal space for a periodic
calculation — the `kgrid.MonkhorstPack` block in `calc.fdf`. `stb-kgrid`
computes a sensible one for you: given a structure and a target k-point
*density*, it works out how many k-points to use along each lattice
direction, correctly forcing a single k-point along any vacuum-padded
(non-periodic) direction. Instead of guessing a grid or hand-deriving one
from the reciprocal lattice, you answer one question — "how dense?" — and
get back a ready-to-use grid.

## Why this matters (a bit of theory)

A periodic crystal's electronic states are labeled by a wavevector **k**
(Bloch's theorem) ranging continuously over the **Brillouin zone** (BZ) —
reciprocal space's own unit cell. Total energy, forces, and everything else
a DFT code computes are, in principle, an *integral* over that zone. In
practice, that integral is approximated by a weighted sum over a finite set
of sampled k-points: too few, and the result is a poor approximation that
keeps changing as you add more; enough, and it converges to a stable value.

**Monkhorst-Pack** (Monkhorst & Pack, 1976 — the exact reference
`stb-kgrid` writes to `references.bib`, see below) is the standard
convention for *how* to lay out that finite set: a uniform mesh of
`N1 x N2 x N3` points along the three reciprocal lattice directions,
instead of picking points ad hoc.

`stb-kgrid` itself does **not** re-derive that paper's own
numerical-optimality argument for where exactly to place points — what it
computes is simpler and more practical: given a target k-point **density**
(in 1/Ang), it picks a mesh size per axis,

```
N_i = ceil( |b_i| / density )
```

where `b_i` is the length of reciprocal lattice vector `i`. A finer density
(smaller number) means more k-points — more accurate, slower. A coarser one
means fewer k-points — faster, less accurate. Metals (sharp features at the
Fermi surface) typically need a denser mesh than semiconductors/insulators
for the same accuracy — this is exactly what the density-recommendation
table (shown below, before you even pick a value) summarizes.

**Dimensionality**: Bloch's theorem, and the whole idea of a k-point, only
applies along a genuinely *periodic* direction. `stb-kgrid` detects a
vacuum-padded axis (the same `core/kspace.py::detect_vacuum_axes` heuristic
`stb-inputfile` also uses) and forces a single k-point there regardless of
density — there's no dispersion to sample along a direction with no real
periodicity, so "finer" sampling there wouldn't mean anything physically.

## When you'd reach for it

Whenever you need a k-grid for a periodic calculation and don't want to
guess one by hand. In practice you rarely run this standalone —
`stb-inputfile` (example `1.1`) already calls the same underlying
`core/kspace.py` logic to fill in `calc.fdf`'s own `kgrid.MonkhorstPack`
line automatically. `stb-kgrid` is for when you want to inspect, justify,
or experiment with that choice on its own (e.g. deciding how much accuracy
a given density buys you before committing to a full calculation).

## Two ways to run it

**A — direct CLI**:

```bash
stb-kgrid --file structure.fdf --density 0.2
```

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 1.3
```

`1.3` asks for the structure file, then — as of this session — shows the
density recommendation guide **before** asking for a density, so the choice
is informed instead of guessed, then runs the exact same `stb-kgrid`
command underneath. `example_1.3.sh` proves both paths agree.

## Files in this folder

- `structure.fdf` — bulk-silicon example structure (same 8-atom conventional
  cubic diamond cell used throughout `examples/1.1-stb-inputfile/`).
- `structure_molecule.fdf` / `structure_chain.fdf` / `structure_graphene.fdf` —
  an isolated CH4 molecule (0D), an infinite carbon chain (1D), and a
  graphene monolayer (2D).
- `example_1.3.sh` — the guided walkthrough (**not** an automated test — see
  `test/1-inputs/3-k-grid/test.sh` for that). Pauses between sections so you
  can read before moving on; safe to re-run.
- `output/` — created by `example_1.3.sh` when you run it (git-ignored, not
  checked in). See below.

## Running the walkthrough

```bash
./example_1.3.sh
```

Every case below is generated with `--save-report`, each into its own
folder under `output/` — `stb_kgrid_report.txt` (the exact console report)
plus `references.bib` (SIESTA + Monkhorst-Pack citations):

| Folder            | Structure                          | Dimensionality | Grid (density 0.2) |
|--------------------|-------------------------------------|:---:|:---------------------|
| `output/silicon/`  | `structure.fdf` (bulk silicon)     | 3D | `6  6  6`   (fully periodic)     |
| `output/molecule/` | `structure_molecule.fdf` (CH4)     | 0D | `1  1  1`   (isolated)           |
| `output/chain/`    | `structure_chain.fdf` (carbon chain)| 1D | `25  1  1`  (a wire)             |
| `output/graphene/` | `structure_graphene.fdf` (graphene)| 2D | `15  15  1` (a slab)             |

### The accuracy/cost tradeoff, made concrete

The script also scans the *same* silicon structure at three densities and
prints the resulting k-point *count* (the product of the grid divisions —
the actual driver of SCF cost) side by side:

| Density | Grid          | k-points | 
|---------|---------------|:--------:|
| 0.1     | `12 12 12`    | 1728     |
| 0.2     | `6  6  6`     | 216      |
| 0.3     | `4  4  4`     | 64       |

Going from 0.3 to 0.1 (3x finer) costs **27x** more k-points — for a 3D
solid, mesh size scales with the *cube* of the linear density change. This
is exactly why picking a sensible density (not just the smallest possible
one) matters in practice.

### The `--vacuum-gap` threshold is a real, tunable choice

Detecting "vacuum" isn't magic — it's a threshold (default 10 Ang): the
largest empty gap along an axis, wrapped periodically, must exceed it to be
treated as non-periodic. The script builds a structure with a 12 Ang gap
along `c` and shows it classified as 2D (grid `6 6 1`) at the default
threshold, but 3D (grid `6 6 2`) once `--vacuum-gap 15` is given — same
structure, different threshold, different physical interpretation *and*
different k-grid.

### Proof: CLI and the interactive menu agree

The script also drives the same silicon/density-0.2 case through the
interactive `stb-suite` menu (non-interactively, via a piped `printf`, no
folder kept for it), confirms the density guide prints before the density
prompt, and diffs the resulting grid line against `output/silicon/` —
proving the CLI and the menu produce identical output.

## Try it yourself

```bash
stb-kgrid --file structure.fdf --density 0.05   # high precision, expensive
stb-kgrid --file structure.fdf --density 0.5    # low precision, cheap
```

Compare the resulting grids and k-point counts against the density-scan
table above.

## Flag reference

```
stb-kgrid --file <structure.fdf> --density <value> [--vacuum-gap <Ang>] [--save-report]
```

- `-f/--file` (required): path to a SIESTA structure file (`.fdf` only).
- `-d/--density` (required): target k-point density, in 1/Ang.
- `--vacuum-gap` (optional, default `10.0`): minimum empty span (Ang) along
  an axis, wrapped periodically, to treat it as vacuum-padded.
- `--save-report` (optional): also persist the report to
  `stb_kgrid_report.txt`. Off by default.

`references.bib` (SIESTA + Monkhorst-Pack citations) is always written —
there's no flag for it.

Run `stb-kgrid --help` for the full list of options.

## What's next

The suggested grid is exactly what belongs in `calc.fdf`'s
`kgrid.MonkhorstPack` line — but you rarely need to copy it there by hand:
`stb-inputfile` (example `1.1`) already computes and writes it for you,
using this same underlying logic.
