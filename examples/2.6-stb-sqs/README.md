# 2.6 — Special Quasirandom Structure Generator (`stb-sqs`)

## What this tool does

Takes a structure with one ordered species and disorders it into a target
substitutional alloy composition (e.g. Ni → 50:50 Ni:Fe), searching for the
one small periodic cell whose local statistics best approximate the true,
infinite random alloy — instead of just placing atoms at random. Same
**Structures** category (2) as `2.1-stb-2Dstacking/`, `2.2-stb-supercell/`,
`2.3-stb-slab/`, `2.4-stb-nanotube/`, `2.5-stb-defect/` — everything here
builds, generates, or transforms a structure file.

## Why this matters (a bit of theory)

### The problem: DFT needs periodicity, real alloys don't have any

A real substitutional alloy or solid solution (NiFe, SiGe, a doped oxide,
a high-entropy alloy, ...) has **no periodicity** — each site is occupied
by one species or another following the bulk composition, with no
repeating pattern. A DFT calculation, on the other hand, needs a finite,
periodic cell. The naive fix — randomly assign species to a small cell's
sites according to the target composition — doesn't actually solve this:
a handful of atoms is nowhere near enough to statistically reproduce a
truly random arrangement's **short-range order** (how often a given
species finds itself next to another, at each neighbor distance). Two
different random draws of the same small cell can behave quite
differently from each other, and from the real, infinite random alloy.

### Special Quasirandom Structures (SQS)

**Zunger, Wei, Ferreira & Bernard, Phys. Rev. Lett. 65, 353 (1990)**
proposed a different approach: instead of a literal random draw, *search*
a small periodic cell for the one atomic arrangement whose **correlation
functions** — how often pairs, triplets, etc. of sites at given distances
are occupied by the same vs. different species — best match those of the
ideal, infinitely large random alloy. The resulting cell isn't random at
all (it's a specific, deterministic arrangement found by the search), but
it's *quasi*-random in the sense that it statistically mimics randomness
far better than an actually-random draw of the same size.

`stb-sqs` wraps [icet](https://icet.materialsmodeling.org/)'s
implementation of this search (via pymatgen's `SQSTransformation`), using
the objective function described by **van de Walle et al., Calphad 42, 13
(2013)** — the same measure used internally by icet's own `generate_sqs`.
Concretely, icet minimizes

```
Q = -w*L + sum_alpha |correlation_alpha - target_correlation_alpha|
```

where the sum runs over every cluster correlation the search is asked to
match (see `--cluster-cutoffs` below), and `-w*L` is a bonus term that
rewards a *perfect* match extending out to a larger cluster radius `w`.
**More negative is better** — a lower (more negative) `Q` means a closer
match to the ideal random alloy's correlation functions. `Q` is what the
tool's `[5] SQS SEARCH` section reports as "Objective function".

### Two search methods, cross-checked live in this example

- `--method monte_carlo` (default) — simulated annealing: cheap, scales to
  large cells, but not guaranteed to find the true optimum (a stochastic
  search can get stuck in a local one).
- `--method enumeration` — exhaustively checks every distinct arrangement
  at the given cell size. Guaranteed to find the **true** global optimum,
  but only feasible for small cells (the number of arrangements grows
  combinatorially with size).

`enumeration/` in this walkthrough runs the exact same problem as
`monte-carlo/` and compares the two objective functions live — a genuine
cross-check that Monte Carlo really found the global optimum, not just a
locally-good one.

### Not every cell size fits a composition exactly

A target composition needs an **exact integer atom count per species** at
whatever cell size is searched. icet auto-detects the smallest
`--scaling` that allows this (`[4] SUPERCELL SIZING`, "Minimal valid
scaling") from its own auto-detected primitive cell; any `--scaling` that
isn't a multiple of it is rejected up front with a clear message, rather
than letting icet's own search silently find zero candidates.
`scaling-constraint/` below demonstrates both the rejection and the fix.

### `--cluster-cutoffs`: more correlations to match = a harder search

`--cluster-cutoffs` (default icet's own `2:3, 3:2, 4:1` — pair clusters
out to the 3rd neighbor shell, triplets out to the 2nd, quadruplets to the
1st) controls exactly which clusters the search tries to match. Asking
for a *wider* set of correlations, with the same small cell, is a
strictly harder target: `cluster-cutoffs/` below shows the objective
function get measurably worse (and the search take measurably longer)
when the default cutoffs are widened, on the exact same problem.

### A real disordering breaks symmetry substantially

Same theme as introducing a point defect (`2.5-stb-defect/`): a perfectly
ordered crystal's high symmetry comes from every site's environment being
equivalent under the crystal's own symmetry operations. Disordering even
a *single* species into two breaks nearly all of that — `[8] SYMMETRY
ANALYSIS` compares the ordered input structure against the final SQS
structure and shows this directly (verified live: FCC Ni's cubic Fm-3m,
No. 225, 48 operations, routinely drops to a low-symmetry monoclinic or
trigonal space group after disordering into Ni-Fe).

### `--sublattice` only disorders the species you name

A structure can already have more than one species (e.g. an oxide).
`--sublattice` targets exactly one of them; every other species keeps its
role in the crystal. `sublattice-mgo/` below verifies this on a real
(Mg,Ni)O rocksalt solid solution (a system studied for catalysis and
electrocatalysis): only the Mg sites are disordered with Ni, and although
icet rebuilds the cell from its own internally-detected primitive cell
(so the O atoms' exact coordinates/ordering can come out different from
the input — the same "equivalent, not identical" caveat noted for
spglib-based tools elsewhere in this suite, e.g. `core/symmetry.py`), the
underlying rocksalt bonding geometry survives: the cation-O
nearest-neighbor distance comes out *exactly* unchanged.

### `--ml-relax`: an SQS cell's local environment is inherently irregular

Unlike a perfect crystal, an SQS cell's atoms genuinely sit at slightly
different local equilibrium positions depending on which species their
neighbors happen to be. `--ml-relax` (+ `--ml-relax-cell`) pre-relaxes the
found SQS structure with a MACE potential before writing it out — the
same convention `stb-supercell`/`stb-slab`/`stb-nanotube`/`stb-defect`
already use.

### Known limitation: practically 2 species per sublattice

`--composition` with 3+ species crashes during the search itself
(`ASE Atoms only supports ordered structures`) — a **verified upstream
issue** in the installed pymatgen's own `SQSTransformation`/`IcetSQS`
wrapper (`pymatgen/io/icet.py`), which builds its internal reference
structure the same way this tool's own `minimal_scaling()` used to (fixed
in this same session — see `--help` for `--composition`), but pymatgen's
own copy of that pattern is out of `stb-sqs`'s reach to fix without
bypassing `SQSTransformation` entirely. Not exercised by this example.

## When you'd reach for it

- Building a DFT-ready cell for a substitutional alloy or solid solution
  (metallic alloys, doped oxides/semiconductors, ...) instead of a single
  hand-picked ordered arrangement.
- Comparing how a specific composition's local order changes the physics
  (relative to the pure end-member) in a subsequent DFT or `stb-ml*`
  calculation.
- Screening a composition ratio before committing to real DFT, by first
  checking `--ml-relax`'s MACE energy on the SQS cell.

## Two ways to run it

A — direct CLI:
```bash
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.6
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.6.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[11] SUMMARY & FILES` —
  `[6]` only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_sqs_report.txt`. Earlier versions of
  this tool always wrote a separate `sqs_report.txt`; that file is gone —
  everything it used to contain is now part of this same numbered report.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density — run once on the input structure and once on the final SQS
  structure.
- **A before/after symmetry comparison table** — crystal system, 3D space
  group, layer group, point group, Hall symbol — comparing the ordered
  input structure against the final (disordered) SQS structure.
- **`references.bib`** — SIESTA always, plus icet (the SQS search
  library itself), plus the MACE papers if `--ml-relax` was used.
- **A provenance header** written into the output `.fdf`: the disordered
  sublattice/composition, the scaling used, the SQS search method and
  objective function, and MACE convergence/energy detail if used.

## Optional (off by default)

- **`--scaling N`** — override icet's auto-detected minimal supercell
  size (must be an integer multiple of it).
- **`--method {monte_carlo,enumeration}`** — simulated annealing
  (default) or exhaustive search.
- **`--instances N`** — run N independent Monte Carlo searches in
  parallel and keep the best (a second global-optimum cross-check).
- **`--temperature`**, **`--mc-steps`** — Monte Carlo search settings.
- **`--cluster-cutoffs SIZE:SHELL,...`** — which correlation functions
  the search tries to match (default icet's own `2:3, 3:2, 4:1`).
- **`--ml-relax`** (+ `--ml-relax-cell`, `--model`/`--custom-model`) —
  pre-relax the found SQS structure generically with MACE.
- **`--save-report`** — also persist the full report to
  `stb_sqs_report.txt`.
- **`--view`** — opens the input structure and the final SQS structure in
  ASE's interactive 3D viewer (`ase-gui`), as pageable frames. Needs a
  local display. Never exercised by `example_2.6.sh` itself.

## Files in this folder

- `fcc_ni.fdf` — bulk FCC nickel, the conventional 4-atom cell (`a` =
  3.52 Å, space group Fm-3m No. 225) — a single-species metal, disordered
  into a binary Ni-Fe alloy by most cases below.
- `mgo_rocksalt.fdf` — MgO, rocksalt (NaCl-type) structure, conventional
  8-atom cell (`a` = 4.246 Å, also Fm-3m No. 225 — a deliberate parallel
  with `fcc_ni.fdf`) — already has two species, used to show
  `--sublattice` disordering only one of them.
- `example_2.6.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_sqs_report.txt`.

## Running the walkthrough

```bash
cd examples/2.6-stb-sqs
./example_2.6.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Eight self-contained cases are generated (one
is skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                | What it shows                                                          |
|------------------------|------------------------------------------------------------------------|
| `monte-carlo/`         | The basic SQS search and its objective function                        |
| `enumeration/`         | Exact search, cross-checked live against Monte Carlo's result          |
| `scaling-constraint/`  | A `--scaling` incompatible with the composition, rejected then fixed   |
| `cluster-cutoffs/`     | Wider correlation targets make the same search measurably harder       |
| `instances/`           | Parallel independent searches, a second global-optimum cross-check     |
| `sublattice-mgo/`      | (Mg,Ni)O: disordering only Mg, verifying O's bonding geometry survives |
| `ml-relax/`            | Generic MACE pre-relaxation of the found SQS structure                 |
| `full-report/`         | `--save-report` + the validation checklist + `references.bib`          |

## Try it yourself

```bash
# Show the input structure and the final SQS structure side by side
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --view

# A different composition ratio
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.25,Fe:0.75

# Pre-relax with MACE, cell included, with your own fine-tuned model
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --ml-relax --ml-relax-cell --custom-model my_finetuned.model

# A larger, more dilute SQS: build a bigger cell first
stb-supercell -f fcc_ni.fdf -d 2 2 2 -o ni32.fdf --no-intro
stb-sqs -f ni32.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-f/--file`          | Input structure file, `.fdf` (required)                                |
| `--sublattice`       | Existing species whose sites become the disordered alloy sublattice (required) |
| `--composition`      | Target composition on that sublattice, e.g. `Fe:0.5,Ni:0.5` (must sum to 1.0; practically 2 species, see limitation above) |
| `--scaling`          | icet supercell-size multiplier (default: auto-detect the minimal valid one) |
| `--method`           | `monte_carlo` (default) or `enumeration`                               |
| `--instances`        | Parallel search instances (default 1)                                  |
| `--temperature`      | Monte Carlo starting temperature (default 1.0)                         |
| `--mc-steps`         | Caps the number of Monte Carlo steps (`--method monte_carlo` only)     |
| `--cluster-cutoffs`  | `size:shell` pairs, e.g. `2:3,3:2` (default icet's own `2:3, 3:2, 4:1`) |
| `-sp/--symprec`      | Symmetry tolerance for the before/after table (default `0.01`)         |
| `--ml-relax`         | Pre-relax the SQS structure found with a MACE potential                |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the cell (vacuum axes stay fixed)        |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file                              |
| `-o/--output`        | Output `.fdf` file name (default `sqs.fdf`)                            |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the input and final SQS structure in ASE's interactive viewer     |

## What's next

See `2.2-stb-supercell/` for building a bigger, more dilute cell before
running an SQS search on it, `2.5-stb-defect/` for another example of a
single missing/substituted atom breaking symmetry substantially, or
`1.6-stb-mlrelax/` for a closer look at the MACE pre-relaxation
`--ml-relax` reuses here.
