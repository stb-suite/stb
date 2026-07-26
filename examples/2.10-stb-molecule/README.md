# 2.10 — Reference Molecule Builder (`stb-molecule`)

## What this tool does

Builds an isolated reference molecule (H2O, CO2, CH4, benzene, ...) from
ASE's bundled **G2 database** (162 small molecules), placed in a vacuum
box ready for a SIESTA calculation. Same **Structures** category (2) as
`2.1-stb-2Dstacking/`, ..., `2.9-stb-passivate/`, but the odd one out
geometrically: every other tool in this category ends up with some kind
of periodicity (a bulk crystal, a slab, a tube); this one deliberately
builds something with **none** — a single, isolated molecule.

## Why this matters (a bit of theory)

### SIESTA is always periodic — an "isolated" molecule needs a trick

SIESTA (like almost every plane-wave/localized-basis DFT code) solves the
electronic structure of a periodic crystal. There is no native "just one
molecule, nowhere else" mode. To model an isolated molecule anyway, you
put it in a large, empty box and repeat *that* periodically — the
molecule's periodic images are far enough apart that they don't
meaningfully interact, so the calculation behaves as if the molecule were
truly alone. `--vacuum` controls how big that empty margin is (in every
direction around the molecule); `stb-molecule` exists specifically to
build this box correctly, so you don't have to hand-place the atoms and
guess a cell size yourself.

### The G2 database: a standard benchmark set, not just a random list

ASE's G2 collection is the "G2/97 test set" from computational
thermochemistry (Curtiss et al.) — a standard, widely-used set of small
molecules with well-characterized reference geometries and energies,
originally assembled to benchmark quantum chemistry methods against
experiment. That's exactly why it's a good source for this tool: every
name in it (`--list` shows all 162) is a real, chemically sensible
molecule with a reasonable starting geometry, not an arbitrary guess.

### Point groups: a molecule's *own* symmetry, distinct from a crystal's

A crystal's symmetry is a **space group** (translations + rotations,
230 possibilities, see `2.8-stb-crystalbuilder/`); an isolated molecule
has no translational symmetry to speak of, so it's classified by its
**point group** instead — the set of rotations/reflections/inversions
that map the molecule onto itself, in the Schoenflies notation chemists
use (`C2v`, `Td`, `D6h`, ...). This is why `stb-molecule`'s own symmetry
detection (`core_symmetry.point_group_label`, via pymatgen's
`PointGroupAnalyzer`) is a completely different code path from the
space-group detection every periodic tool in this suite uses — there is
no periodic cell here for a space group to even apply to.
`point-group-gallery/` below verifies 4 real, textbook point groups
directly from the built geometries: water's bent `C2v`, methane's
tetrahedral `Td`, benzene's hexagonal `D6h`, and CO2's linear `D*h`
(pymatgen's notation for the infinite-order D∞h).

### Why symmetry is expected to survive a relaxation — and what it means if it doesn't

The G2 database's own reference geometries are already close to each
molecule's true equilibrium shape, so relaxing with MACE (`--ml-relax`)
is expected to preserve the point group exactly — a small, real
correction to the bond lengths/angles, not a change in shape. Every run's
own `[5] SYMMETRY ANALYSIS (BEFORE / AFTER)` table proves this instead of
just asserting it. If a real calculation on your OWN starting geometry
ever shows a point-group *change* here, that's a genuine, useful signal:
either the relaxation found a genuinely different (lower-symmetry, lower
-energy) minimum, or the starting geometry wasn't as symmetric as it
looked.

### A real pitfall: `--vacuum` too small fools the tool's own dimensionality check

`stb-molecule` always builds something physically 0D (fully isolated),
but the report's own "Dimensionality" line is computed generically (same
`core.kspace.detect_vacuum_axes` heuristic every periodic-aware tool in
this suite uses), by checking whether the gap between periodic images
exceeds a fixed 10 Ang threshold on each axis. The **default** `--vacuum
10.0` sits right at that threshold, correctly reading `0D (e.g., a
molecule)`. `vacuum-and-dimensionality/` below shows what happens if you
turn it down: at `--vacuum 3`, the box is real and the molecule is still
just as isolated physically, but the report's own dimensionality check
gets fooled into reporting `3D (bulk material)` — and even prints an
`Atomic density` validation row that means nothing for a lone molecule in
a box. The **point group**, by contrast, is completely unaffected (it
never looks at the box at all) and stays correct either way. Lesson:
trust the point group here, and keep `--vacuum` at its default (or
higher) for the dimensionality line to read correctly too.

## When you'd reach for it

- An isolated-molecule energy reference for a workflow that needs one —
  e.g. `stb-cohesive`'s isolated-atom references, or a gas-phase molecule
  (H2, O2, H2O) needed by `stb-her`/`stb-oer`'s computational hydrogen
  electrode descriptors.
- A quick, chemically sensible starting geometry for any small-molecule
  calculation, without hand-building the structure file yourself.

## Two ways to run it

A — direct CLI:
```bash
stb-molecule --name H2O
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.10
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.10.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[8] SUMMARY & FILES` —
  `[3]` only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_molecule_report.txt`.
- **Structure validation** — atom proximity, lattice handedness (atomic
  density skipped: always vacuum-padded on all 3 axes by construction) —
  run once before and once after any ML relax.
- **Point-group symmetry detection**, before and after, in a dedicated
  `[5] SYMMETRY ANALYSIS` table (see the theory above for why this is a
  point group, not a space group).
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax` was used.
- **A provenance header** written into the output `.fdf`: the molecule
  name, vacuum padding, detected point group, and MACE convergence/energy
  detail if used.

## Optional (off by default)

- **`--vacuum`** — padding in Ang around the molecule (default `10.0`;
  see the pitfall above before lowering it).
- **`--ml-relax`** (+ `--model`/`--custom-model`) — pre-relax the
  molecule's geometry with MACE. Positions only — there is no
  `--ml-relax-cell` option here at all (not just a no-op): an isolated
  molecule's vacuum box has no physically meaningful cell to relax.
- **`--save-report`** — also persist the full report to
  `stb_molecule_report.txt`.
- **`--view`** — opens the as-built and final molecule in ASE's
  interactive 3D viewer (`ase-gui`), as pageable frames. Needs a local
  display. Never exercised by `example_2.10.sh` itself.
- **`--list`** — prints all 162 available names and exits (names are
  case-sensitive).

## Files in this folder

- `example_2.10.sh` — the guided walkthrough (see below). No input
  fixtures needed — like `2.8-stb-crystalbuilder/`, everything is built
  from `--name`/`--vacuum` flags, not read from a structure file.
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_molecule_report.txt`.

## Running the walkthrough

```bash
cd examples/2.10-stb-molecule
./example_2.10.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Six self-contained cases are generated (one
is skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                        | What it shows                                                          |
|--------------------------------|-------------------------------------------------------------------------|
| `point-group-gallery/`         | 4 real point groups verified directly: water `C2v`, methane `Td`, benzene `D6h`, CO2 `D*h` |
| `unknown-name-suggestions/`    | `--list`, a case-insensitive "did you mean", and a genuine typo suggestion |
| `vacuum-and-dimensionality/`   | The real pitfall: `--vacuum` below the detection threshold fools the dimensionality check |
| `ml-relax-symmetry-preserved/` | MACE pre-relaxation of methane, with `Td` proven unchanged before/after |
| `full-report/`                 | `--save-report` + `references.bib`                                     |

## Try it yourself

```bash
# Any of the 162 G2 names -- see them all
stb-molecule --list

# A reference molecule for stb-her/stb-oer's computational hydrogen electrode workflow
stb-molecule --name H2

# Pre-relax with your own fine-tuned MACE model before a real SIESTA calculation
stb-molecule --name CH3OH --ml-relax --custom-model my_finetuned.model
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--name`             | Exact, case-sensitive G2 database name (required unless `--list`)      |
| `--list`             | Print all 162 available names and exit                                 |
| `--vacuum`           | Vacuum padding in Ang around the molecule (default `10.0`)              |
| `--ml-relax`         | Pre-relax the molecule's geometry with a MACE potential (positions only) |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file                              |
| `-o/--output`        | Output `.fdf` file name (default `molecule.fdf`)                       |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the as-built and final molecule in ASE's interactive viewer       |

## What's next

See `4.3-cohesive/` (Workflow) for a real use of an isolated reference
(atoms there, not molecules, but the same "periodic code needs an
artificial vacuum box" idea), or `2.8-stb-crystalbuilder/` for the
periodic, space-group side of the symmetry story this tool's point-group
detection deliberately steps outside of.
