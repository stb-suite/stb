# 2.5 — Point Defect Generator (`stb-defect`)

## What this tool does

Introduces a point defect — a **vacancy** (remove an atom), a
**substitution** (swap an atom's species), or an **interstitial** (add a
new atom) — into a structure, choosing the site by raw atom index, by the
position closest to a target, or automatically at *every* symmetrically
distinct site at once. Same **Structures** category (2) as
`2.1-stb-2Dstacking/`, `2.2-stb-supercell/`, `2.3-stb-slab/`,
`2.4-stb-nanotube/` — everything here builds, generates, or transforms a
structure file.

## Why this matters (a bit of theory)

### Three kinds of point defect, three different physical roles

- **Vacancy** — a missing atom. The simplest defect, and a real
  equilibrium feature of every crystal at finite temperature; also the
  starting point for diffusion (an atom hops into a vacancy) and for many
  color centers (e.g. diamond's nitrogen-vacancy center).
- **Substitution** — one atom's species swapped for another. The
  mechanism behind doping (Si:P, Si:B for n-/p-type semiconductors) and
  alloying.
- **Interstitial** — an extra atom squeezed in between existing ones.
  Common for small atoms (H, Li, N) diffusing through a host lattice —
  e.g. hydrogen embrittlement, or Li-ion battery intercalation.

### A point defect always breaks symmetry — verified, and it's dramatic

A perfect crystal's high symmetry comes from every atom's environment
being equivalent to many others under the crystal's own symmetry
operations. A single point defect, embedded in an otherwise-periodic cell,
breaks nearly all of that at once: the defect site is no longer equivalent
to its neighbors, and any symmetry operation that used to map it onto a
different (now unperturbed) atom no longer holds. `vacancy-index/` below
shows this directly and quantitatively: removing **one** atom from a
perfectly symmetric 8-atom cubic silicon cell (space group Fd-3m, No. 227,
48 symmetry operations) drops it straight to P-43m (No. 215, 24
operations) — losing inversion symmetry entirely, from removing a single
atom out of eight.

### Not every site is equivalent — `--all-inequivalent-sites`

In a crystal with more than one Wyckoff position for a given element, a
vacancy/substitution at one site is **not** the same defect as one at
another — different local coordination, different formation energy,
different physics. Picking a site by hand (`--index`) requires already
knowing which sites are distinct; `--all-inequivalent-sites` instead finds
every symmetrically distinct site automatically (via spglib) and writes
one output structure per site — see `all-inequivalent-sites/` below, using
magnetite (Fe₃O₄), whose Fe sits on two genuinely different sites
(tetrahedral and octahedral) while O sits on only one.

### `--ml-rank`: which site is actually more stable?

Knowing the sites is only half the story — `--ml-rank` (only valid with
`--all-inequivalent-sites`) relaxes each candidate with a MACE potential
and ranks them by energy, a fast pre-screen for which site is worth real
DFT time. `ml-rank/` below reproduces a genuine, physically sensible
result: an Fe vacancy at magnetite's octahedral site comes out **~0.06 eV
more stable** than at the tetrahedral site — a real, non-degenerate
answer, not a coin flip.

### `--ml-rank` vs. the generic `--ml-relax`

Both use the same MACE machinery, for two different jobs:
- **`--ml-rank`** (only with `--all-inequivalent-sites`) relaxes *every*
  candidate site and **ranks** them — the point is comparing sites.
- **`--ml-relax`** (any selection mode, including a single `--index`
  choice) just pre-relaxes whichever structure(s) you're already writing
  — the point is a better starting geometry, not a comparison. The two are
  mutually exclusive when both would apply to the same run (using
  `--ml-rank` already relaxes every candidate; there's nothing left for
  `--ml-relax` to add without picking a winner first).
`ml-relax/` below uses the generic form, on a single chosen vacancy.

### Diluting a defect: pair this tool with `stb-supercell`

A defect embedded in a small periodic cell interacts with its own
periodic images — usually not what you want. The fix is a big enough
supercell first (`2.2-stb-supercell/`): the defect *concentration* is
`1 / N_atoms`, so a bigger cell means a more dilute, more physically
realistic defect. This example's fixtures are deliberately small (fast,
readable output); a real calculation would typically run `stb-supercell`
first, then `stb-defect` on the result.

## When you'd reach for it

- Building a vacancy/substitution/interstitial structure for a defect
  formation-energy or diffusion-barrier DFT calculation.
- Screening every symmetrically distinct site of a multi-site element
  before committing to expensive DFT on all of them.
- Doping a semiconductor host at a specific, hand-picked site.
- Getting a MACE-pre-relaxed starting geometry for a defect calculation.

## Two ways to run it

A — direct CLI:
```bash
stb-defect -f si_bulk.fdf --type vacancy --index 1
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.5
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.5.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[7] SUMMARY & FILES` — `[4]`
  only appears with `--ml-relax`/`--ml-rank`) printed to the console and,
  with `--save-report`, also saved to `stb_defect_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density — run once on the input structure and once per output
  structure.
- **A before/after symmetry comparison table** per output structure —
  crystal system, 3D space group, layer group, point group, Hall symbol.
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax`/`--ml-rank` was used.
- **A provenance header** written into every output `.fdf`: the defect
  type, which atom(s)/site were involved, and MACE convergence/energy
  detail if used.

## Optional (off by default)

- **`--nearest X Y Z`** (+ `--nearest-format`, `--filter-species`) — pick
  the site closest to a target position instead of a raw index.
- **`--all-inequivalent-sites`** (+ `--filter-species`) — one output
  structure per symmetrically distinct site (vacancy/substitution only).
- **`--ml-rank`** (needs the optional `ml` extra: `pip install
  stb_suite[ml]`, only with `--all-inequivalent-sites`) — relax and rank
  every candidate site by MACE-relaxed energy.
- **`--ml-relax`** (+ `--ml-relax-cell`, `--model`/`--custom-model`) —
  pre-relax the output structure(s) generically, any selection mode.
  Mutually exclusive with `--ml-rank`.
- **`--save-report`** — also persists the full report to
  `stb_defect_report.txt`.
- **`--view`** — opens the input structure and every output structure in
  ASE's interactive 3D viewer (`ase-gui`), as pageable frames. Needs a
  local display. Never exercised by `example_2.5.sh` itself.

## Files in this folder

- `si_bulk.fdf` — bulk silicon, diamond cubic, the real conventional
  8-atom cell (`a` = 5.431 Å, space group Fd-3m No. 227) — highly
  symmetric on purpose, so the symmetry-breaking demo is unambiguous.
- `magnetite.fdf` — magnetite (Fe₃O₄), 14 atoms, space group Fd-3m
  (No. 227) — Fe on two genuinely distinct sites (tetrahedral,
  octahedral), O on one.
- `example_2.5.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_defect_report.txt`.

## Running the walkthrough

```bash
cd examples/2.5-stb-defect
./example_2.5.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Eight self-contained cases are generated (two
are skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                     | What it shows                                                     |
|------------------------------|----------------------------------------------------------------------|
| `vacancy-index/`             | A Si vacancy by `--index` — the symmetry-breaking demo, Fd-3m -> P-43m |
| `substitution/`               | Doping, Si -> Ge by `--index`                                       |
| `interstitial/`               | Adding an N interstitial at a fractional position                   |
| `nearest-position/`           | Selecting a site by `--nearest`, with periodic wraparound            |
| `all-inequivalent-sites/`     | Magnetite's 2 distinct Fe sites, found automatically                |
| `ml-rank/`                    | Ranking those 2 Fe sites by MACE-relaxed energy (octahedral wins)   |
| `ml-relax/`                   | Generic MACE pre-relaxation of a single chosen vacancy               |
| `full-report/`                | `--save-report` + the validation checklist + `references.bib`       |

## Try it yourself

```bash
# Show the input structure and the defect structure side by side
stb-defect -f si_bulk.fdf --type vacancy --index 1 --view

# Every inequivalent O site in magnetite, no MACE needed
stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species O

# Pre-relax with MACE, cell included, with your own fine-tuned model
stb-defect -f si_bulk.fdf --type vacancy --index 1 --ml-relax --ml-relax-cell --custom-model my_finetuned.model

# Dilute the defect first: a bigger cell, then a vacancy in it
stb-supercell -f si_bulk.fdf -d 2 2 2 -o si64.fdf --no-intro
stb-defect -f si64.fdf --type vacancy --index 1
```

## Flag reference

| Flag                     | Meaning                                                             |
|----------------------------|-------------------------------------------------------------------------|
| `-f/--file`                | Input structure file, `.fdf` (required)                               |
| `--type`                   | `vacancy` / `substitution` / `interstitial` (required)                |
| `--index`                  | Comma-separated 1-indexed atom(s), e.g. `3,7`                          |
| `--nearest X Y Z`          | Select the atom closest to this position                              |
| `--nearest-format`         | `fractional` (default) or `cartesian`, for `--nearest`                |
| `--all-inequivalent-sites` | One output structure per symmetrically distinct site                  |
| `--filter-species`         | Restrict `--nearest`/`--all-inequivalent-sites` to one element        |
| `--symprec`                | Symmetry precision for site-finding and the before/after table (default `1e-3`) |
| `--ml-rank`                | Relax and rank every candidate site (needs `--all-inequivalent-sites`) |
| `--new-species`            | Replacement element, required for `--type substitution`               |
| `--position X Y Z`         | New atom's position, required for `--type interstitial`               |
| `--position-format`        | `fractional` (default) or `cartesian`, for `--position`                |
| `--species`                 | New atom's element, required for `--type interstitial`                |
| `--ml-relax`               | Pre-relax the output structure(s) with a MACE potential               |
| `--ml-relax-cell`          | With `--ml-relax`, also relax the cell (vacuum axes stay fixed)       |
| `--model`                  | MACE-MP-0 size for `--ml-relax`/`--ml-rank`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`           | Path to a custom fine-tuned `.model` file                             |
| `-o/--output`               | Output `.fdf` file name (default `defect.fdf`)                        |
| `--save-report`            | Persist the full report (incl. symmetry table) to disk                |
| `--view`                    | Open the input and output structure(s) in ASE's interactive viewer    |

## What's next

See `2.2-stb-supercell/` for building a bigger, more dilute cell before
introducing a defect into it, or `1.6-stb-mlrelax/` for a closer look at
the MACE pre-relaxation `--ml-relax`/`--ml-rank` reuse here.
