# 2.12 — Random Crystal Generator (`stb-crystalcast`)

## What this tool does

`stb-crystalcast` wraps **pyxtal** to generate or transform structures via
symmetry groups, in 4 modes:

- **`generate`** (default, no mode flag) — give a symmetry group + a
  composition, pyxtal randomly assigns atoms to compatible Wyckoff
  positions and draws the free coordinates. This is the **inverse** of
  `stb-crystalbuilder` (2.8): there, you already know the exact Wyckoff
  sites and want to build; here, you only know the composition and want
  pyxtal to find valid candidate placements for you (e.g. as starting
  guesses for structure prediction).
- **`--substitute`** — swap elements in an existing structure while
  preserving the symmetry framework (every atom keeps its Wyckoff site).
- **`--subgroup`** / **`--supergroup`** — search for a related
  lower-/higher-symmetry structure of an existing one. The only
  symmetry-breaking/-restoring search anywhere in this suite.

**`--analyze` (Wyckoff decomposition of an existing structure) was removed**
in this pass — use `stb-symmetry` or `stb-unitcell` instead to inspect an
existing structure's own symmetry.

## How this differs from similar-sounding tools in the suite

Several other tools also involve "symmetry" or "generating a structure
from scratch" — none of them do what `stb-crystalcast` does, but it's
easy to confuse them at a glance. `example_2.12.sh`'s first case builds
the *exact same* real rock-salt NiO through both `stb-crystalbuilder` and
`stb-crystalcast` live, side by side, as direct proof of the first row
below:

| Tool | What it actually does | Why it's not `stb-crystalcast` |
|------|------------------------|----------------------------------|
| **`stb-crystalbuilder`** (2.8) | Builds a structure from a space group + the exact Wyckoff sites *you already know* | The **inverse** problem: you supply the placement, it expands it. It has no composition to check against — give it the wrong site and it silently builds a different, but still "valid", structure. `stb-crystalcast` only ever needs the composition; it works out (or randomizes) the placement itself. |
| **`stb-sqs`** (2.6) | Fits a **Special Quasirandom Structure** on a *fixed, already-existing* sublattice, matching real-alloy short-range-order statistics (via `icet`'s cluster expansion) | Solves a completely different problem: it never touches symmetry groups or Wyckoff positions at all — it's about how to arrange species on sites you *already have* so the *local correlations* look like a real random alloy, not about generating a new cell from a group. |
| **`stb-defect`** `--all-inequivalent-sites` | Enumerates symmetrically distinct **sites in an existing structure**, to choose where to place ONE vacancy/substitution/interstitial | Uses the same underlying idea (symmetry-equivalent orbits, via `SpacegroupAnalyzer`) as pyxtal's Wyckoff search, but the goal is picking a defect site in a structure you already have, not building/transforming a whole new one. |
| **`stb-symmetry`** / **`stb-unitcell`** | Analyze an *existing* structure's own symmetry (operations, equivalent atoms, primitive/conventional cell) | This is what `--analyze` used to do here before this session removed it. `stb-crystalcast` **generates or transforms**; it doesn't have an inspection-only mode of its own anymore. |
| **`stb-molecule`** | Builds **one** isolated reference molecule (from ASE's G2 database) in a vacuum box | `stb-crystalcast --molecular` instead packs whole rigid molecules (from *pyxtal's own* bundled collection, not G2) **into a periodic symmetry group** — a molecular crystal candidate, not a single isolated reference. `--dim 0` (bare atoms, point group, vacuum box) is the closer analog of `stb-molecule`'s own vacuum convention, but still composition-driven, not a G2 lookup by name. |
| **`stb-mlsearch`** | Basin-hopping / simulated-annealing search for the **best atomic arrangement at a fixed composition and cell** | A real optimization over configuration space at one fixed composition/cell. `stb-crystalcast --ml-rank` is much simpler: it generates independent random candidates first (across possibly many distinct symmetry groups/cells), then just relaxes and ranks each one once — a cheap pre-screen, not a search. |

## Why this matters (a bit of theory)

### Wyckoff positions: why "group + composition" is enough to build a structure

A space/layer/rod/point group's own symmetry operations partition all
possible atomic positions into **Wyckoff positions** — orbits of
symmetry-equivalent sites. Placing one atom at a Wyckoff site
automatically generates every symmetry-equivalent copy of it (the site's
*multiplicity*). `generate` mode's whole job is finding a combination of
Wyckoff sites whose multiplicities sum to your requested composition
(e.g. 4 Ni + 8 O for space group 225), then drawing random coordinates
for whichever sites still have free parameters (some Wyckoff sites are
fully fixed by symmetry alone — no randomness needed there at all).
`--sites` lets you pin specific Wyckoff labels yourself instead of
leaving every site random — the same information `stb-crystalbuilder`
needs as direct input.

### `--subgroup`/`--supergroup`: symmetry-breaking and symmetry-restoring searches

No other tool in this suite does this. `--subgroup` searches for related
structures with **lower** symmetry — a controlled symmetry-breaking
distortion (`--group-type t/k`, translationengleiche/klassengleiche
group-subgroup relations; `--eps` is the coordinate perturbation applied)
— useful for generating candidate low-symmetry phases of a known
high-symmetry parent (e.g. a real ferroelectric distortion of a cubic
perovskite). `--supergroup` does the reverse: given a (possibly
distorted) structure and a target *higher*-symmetry group, it searches
for the symmetric parent structure group-subgroup theory says it should
be a distortion of (`--d-tol` is the maximum displacement tolerated when
testing a candidate match). Verified live that this genuinely can go
either way for a real distortion: pyxtal sometimes finds a matching
parent, and sometimes correctly reports "not found" instead of a wrong
answer — both are legitimate outcomes of a real search, not a crash.

### A real, verified pymatgen quirk this tool works around

While standardizing this report, live testing surfaced a real
inconsistency in pymatgen's `SpacegroupAnalyzer`: for a `--subgroup`
candidate distorted with the default `--eps 0.05`, it can report a space
group that hasn't changed (`Fm-3m (225)`) *alongside* a point group that
implies a genuinely different, lower symmetry (`4/mmm`, not `m-3m`) — an
internally inconsistent pair from the exact same analyzer call. Since
`--subgroup`/`--supergroup` candidates already carry pyxtal's own
authoritative group assignment (literally what the search built towards),
`stb-crystalcast` now reports **that** label for "Space Group" in its
`[5] SYMMETRY ANALYSIS` table instead of re-detecting it via pymatgen —
verified live that this correctly shows `I4/mmm (No. 139)` for a targeted
`--subgroup` search, where the pymatgen re-detection alone would have
misleadingly still said `Fm-3m (225)`.

### `--ml-rank`: now with a choice of model

`--ml-rank` relaxes every generated candidate's positions with a MACE
potential and ranks them by relaxed energy — a fast, relative pre-screen
for which candidate(s) to prioritize for real DFT, not an absolute DFT
formation energy. Previously always used the small MACE-MP-0 foundation
checkpoint with no way to choose; `--model`/`--custom-model` (only valid
together with `--ml-rank`) now match the same convention every other
`--ml-rank`/`--ml-relax` feature in this suite already has
(`stb-crystalbuilder`, `stb-defect`, `stb-amorphize`, ...).

## When you'd reach for it

- A batch of symmetry-valid candidate structures for a given composition,
  when you don't already know (or want to hand-pick) the exact Wyckoff
  sites — use `stb-crystalbuilder` instead once you do.
- Screening candidate low-/high-symmetry phases related to a known
  structure (`--subgroup`/`--supergroup`), e.g. exploring possible
  ferroelectric distortions of a known high-symmetry parent.
- A quick element-substitution variant of an existing structure that
  keeps the same symmetry framework (`--substitute`).

## Two ways to run it

A — direct CLI:
```bash
stb-crystalcast --group 225 --species Ni O --num-ions 4 8
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.12
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.12.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[8] SUMMARY & FILES`)
  printed to the console and, with `--save-report`, also saved to
  `stb_crystalcast_report.txt`.
- **Structure validation** (atom proximity, lattice handedness, atomic
  density) on every candidate written.
- **A symmetry table** — before/after (`--substitute`/`--subgroup`/
  `--supergroup`) or requested-vs-detected (`generate`).
- **`references.bib`** — SIESTA + PyXtal always; + MACE/MACE-MP with
  `--ml-rank`. Previously never written at all in any mode.
- **A provenance header** written into every output `.fdf`: requested
  group/composition (or input file/substitutions/target group), detected
  symmetry, and the `--ml-rank` outcome if used.

## Optional (off by default)

- **`--count`** — generate/keep several independent candidates in one run
  (numbered `<output>_1.fdf`, `_2.fdf`, ...).
- **`--seed`** — reproducible batches (generation mode; with `--molecular`,
  only the lattice is reproduced this way, an upstream pyxtal limitation).
- **`--molecular`** — pack whole rigid molecules instead of bare atoms.
- **`--dim 0/1/2`** — point-group cluster / rod group / layer group,
  instead of the default `--dim 3` space group.
- **`--sites`** — pre-assign specific Wyckoff positions instead of fully
  random.
- **`--ml-rank`** (+ `--model`/`--custom-model`) — MACE relaxed-energy
  screening of generated candidates.
- **`--save-report`** / **`--view`** — same as every other stb-suite tool.

## Files in this folder

- `example_2.12.sh` — the guided walkthrough (**not** an automated test —
  see `test/2-structures/12-crystalcast/test.sh` for that). No committed
  input fixtures: like `2.8-stb-crystalbuilder/`, `generate`/`--molecular`/
  `--dim` build everything from flags; the one structure
  `--substitute`/`--subgroup`/`--supergroup` need is generated by this
  script's own first case, not checked in.
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_crystalcast_report.txt`.

## Running the walkthrough

```bash
cd examples/2.12-stb-crystalcast
./example_2.12.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Seven or eight self-contained cases are
generated (`ml-rank/` skipped, with a clear message, if the optional `ml`
extra isn't installed):

| Folder                      | What it shows                                                          |
|-------------------------------|-------------------------------------------------------------------------|
| `vs-crystalbuilder/`          | The exact same real NiO built two opposite ways — live proof of the comparison table above |
| `generate/`                   | Rock-salt NiO from `--group 225` + composition; `--save-report`, provenance header, `references.bib` |
| `count-and-seed/`              | `--count 3 --seed 7`: distinct candidates, reproducible batches         |
| `molecular-and-clusters/`      | `--molecular` (water in group 19) and `--dim 0` (isolated point-group cluster) |
| `substitute/`                  | `O:S` substitution on the NiO structure, symmetry framework preserved   |
| `subgroup-supergroup/`         | Lower-symmetry search, then a higher-symmetry search back               |
| `ml-rank/`                     | `--ml-rank --model small`: candidates ranked by MACE relaxed energy     |

## Try it yourself

```bash
# Pre-assign specific Wyckoff sites instead of leaving everything random
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 --sites 4a 4b

# See every molecule name bundled with pyxtal (usable with --molecular)
stb-crystalcast --list-molecules

# Screen candidates with your own fine-tuned MACE model
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --count 10 \
    --ml-rank --custom-model my_finetuned.model
```

## Flag reference

| Flag                  | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--group`             | Symmetry group (number or symbol); meaning depends on `--dim`           |
| `--species`/`--num-ions` | Composition (generation mode)                                        |
| `--dim`               | 3=space group (default), 2=layer, 1=rod, 0=point group                  |
| `--molecular`         | Pack whole rigid molecules instead of bare atoms                        |
| `--sites`             | Pre-assign specific Wyckoff positions                                   |
| `--lattice`           | Fix the cell instead of estimating it from `--volume-factor`            |
| `--count`/`--seed`    | Batch size / reproducibility (generation, `--subgroup`, `--supergroup`) |
| `--substitute`        | `OLD:NEW` element swap(s) on an existing structure (`-f`)                |
| `--subgroup`          | Lower-symmetry search on an existing structure (`-f`)                   |
| `--supergroup`        | Higher-symmetry search on an existing structure (`-f`, needs `--target-group`) |
| `--target-group`/`--group-type`/`--eps`/`--d-tol` | `--subgroup`/`--supergroup` search parameters |
| `--ml-rank`           | Rank generated candidates by MACE relaxed energy                        |
| `--model`/`--custom-model` | MACE model for `--ml-rank` (only valid together with it)          |
| `-o/--output`         | Output `.fdf` file name (default `crystalcast.fdf`)                     |
| `--save-report`       | Persist the full report to disk                                         |
| `--view`              | Open the written structure(s) in ASE's interactive viewer                |

## What's next

See `2.8-stb-crystalbuilder/` for the inverse workflow (you already know
the exact Wyckoff sites), or `stb-symmetry`/`stb-unitcell` to inspect an
existing structure's own symmetry/Wyckoff decomposition (what `--analyze`
used to do here). Take a promising `--ml-rank`-screened candidate into a
real SIESTA relaxation next.
