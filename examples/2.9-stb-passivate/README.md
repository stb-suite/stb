# 2.9 — Surface Passivator (`stb-passivate`)

## What this tool does

Caps **dangling bonds** on a cut surface (typically a slab from `stb-slab`,
`2.3`) with a passivating atom (default H), placed purely from each
undercoordinated atom's own local coordination geometry — no bulk-lattice
assumptions hard-coded in. Same **Structures** category (2) as
`2.1-stb-2Dstacking/`, ..., `2.8-stb-crystalbuilder/`. Also available as a
same-step convenience flag on `stb-slab --passivate` — this tool is the
standalone version, usable on any structure, not just fresh-from-`stb-slab`
output.

## Why this matters (a bit of theory)

### Cutting a crystal creates unphysical dangling bonds

A bulk crystal's atoms are all fully coordinated — every bond an atom
"should" have (given its position in the lattice) is actually present.
Cutting a slab out of that bulk (to model a surface) necessarily severs
some of those bonds: the atoms right at the new surface are left with
fewer neighbors than they have in the bulk. In a real material, these
broken bonds get satisfied somehow — surface reconstruction, adsorbed
atoms/molecules, an oxide layer. A bare, idealized DFT slab has none of
that: the dangling bonds are simply unphysical half-bonds, which
introduce spurious mid-gap electronic states localized at the surface —
states that don't exist in the real material and can contaminate the very
electronic structure (band gap, DOS, work function, ...) you built the
slab to study in the first place. Capping each dangling bond with a
cheap, chemically simple passivant (H is the standard choice) removes
these spurious states, letting the slab's interior behave like the real
bulk while the surface you actually care about (the *other* side of the
slab, in an asymmetric-termination or adsorption study) stays real.

### The geometry: an exact result for one missing bond, a genuine ambiguity for two or more

For an atom missing exactly one bond, the direction is **exactly
determined**: the sum of the unit vectors to its existing neighbors, when
negated, points precisely at the ideal position for the missing one — a
direct geometric consequence of an atom's own equilibrium bonding
geometry (e.g. tetrahedral for Si), not an approximation.
`good-termination/` below proves this isn't just a plausible heuristic:
the recovered H position sits at *exactly* 109.47°, the textbook ideal
tetrahedral angle, from every one of Si's existing bonds.

For an atom missing **two or more** bonds, this breaks down: with only
one (or zero) existing bond vectors left to anchor the geometry, there
are infinitely many ways to place the missing atoms that still give the
right coordination number — nothing pins down the actual directions
without extra assumptions this tool deliberately doesn't make (it has no
hard-coded knowledge of what "tetrahedral" or any other target geometry
even means for the species involved). Guessing would risk silently
writing an unphysical structure; instead, these sites are reported and
left unpassivated for manual review. `bad-termination/` below shows this
directly, with a real, physical cause: cutting the *exact same* Si(111)
surface through a different atomic plane (the strongly-bonded interior of
a bilayer, instead of the weakly-bonded gap between bilayers) leaves each
surface atom missing 3 bonds instead of 1 — geometrically ambiguous, by
the reasoning above, not a limitation of this particular structure.

### Passivation should be symmetry-preserving, when it's the clean, single-missing-bond case

Since each capped bond is placed exactly along the direction the
structure's own local symmetry dictates, adding a well-determined
passivant never breaks the crystal's overall symmetry. Every run's own
`[6] SYMMETRY ANALYSIS (BEFORE / AFTER)` table proves this concretely,
not just asserts it — `symmetry-preserved/` below shows a real Si(111)
slab's space group, layer group, point group, and Hall symbol all coming
back *identical* before and after passivation.

### `--bond-length`/`--cutoff`: a physical approximation, and how to correct it

The default bond length is the sum of the two species' pymatgen atomic
radii — a fast, generic approximation (e.g. Si+H = 1.35 Ang), not the
real experimental value (Si-H is actually closer to ~1.48 Ang). Good
enough as a starting geometry, but `--bond-length` lets you dial in the
real value for your specific system before a real DFT relaxation — or
just let `--ml-relax` correct it automatically (see below).

## When you'd reach for it

- Right after `stb-slab` cuts a surface with an unphysical unsaturated
  termination, before handing the structure to SIESTA.
- On any structure (not just a fresh `stb-slab` output) with genuinely
  undercoordinated atoms you want to cap.
- As a fast pre-check for whether a given termination is even a good
  candidate: a termination that comes back with mostly `[WARNING]`
  unresolved sites is telling you it's a geometrically awkward cut (see
  `bad-termination/`) — maybe try a different termination index instead.

## Two ways to run it

A — direct CLI:
```bash
stb-passivate -f slab.fdf
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.9
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.9.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[9] SUMMARY & FILES` — `[4]`
  only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_passivate_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density (skipped for a vacuum-padded slab, since density isn't
  meaningful there) — run once before and once after passivation.
- **Full symmetry detection** via a before/after comparison table —
  crystal system, 3D space group, **layer group** (meaningful here, since
  a typical input is a genuinely 2D-periodic slab), point group, Hall
  symbol.
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax` was used.
- **A provenance header** written into the output `.fdf`: the passivant,
  how many dangling bonds were found/auto-passivated/left unresolved, and
  MACE convergence/energy detail if used.

## Optional (off by default)

- **`--passivant`** — element to cap dangling bonds with (default `H`).
- **`--cutoff`** — neighbor-search radius for coordination counting
  (default: auto, 1.25x the shortest pairwise distance).
- **`--bond-length`** — passivant bond length, Ang (default: auto, sum of
  atomic radii).
- **`--ml-relax`** (+ `--ml-relax-cell`, `--model`/`--custom-model`) —
  pre-relax the passivated structure with MACE. Vacuum-aware: any
  vacuum-padded axis (the typical slab case) always stays exactly fixed
  under `--ml-relax-cell`, same convention as `stb-unitcell`/`stb-slab`.
- **`--save-report`** — also persist the full report to
  `stb_passivate_report.txt`.
- **`--view`** — opens the input and passivated structures in ASE's
  interactive 3D viewer (`ase-gui`), as pageable frames. Needs a local
  display. Never exercised by `example_2.9.sh` itself.

## Files in this folder

- `example_2.9.sh` — the guided walkthrough (see below). No committed
  input fixtures: like `2.8-stb-crystalbuilder/`'s own
  `bulk-graphite-to-slab/` case, this example builds bulk Si with
  `stb-crystalbuilder` and cuts both Si(111) terminations with `stb-slab`
  itself, live, at the start of the script — the exact same 2-step
  pipeline (`2.8` -> `2.3`) a real user would run before ever reaching
  for `stb-passivate`.
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_passivate_report.txt`.

## Running the walkthrough

```bash
cd examples/2.9-stb-passivate
./example_2.9.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Six self-contained cases are generated (one
is skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                     | What it shows                                                          |
|-----------------------------|-------------------------------------------------------------------------|
| `good-termination/`         | The clean case: 8 deficit-1 sites, all auto-capped at the exact 109.47° tetrahedral angle |
| `bad-termination/`          | The *same* surface, cut one atomic plane over: 8 deficit-3 sites, all correctly left unresolved |
| `bond-length-override/`     | The default 1.35 Ang (atomic-radii-sum) approximation vs. a dialed-in real Si-H bond length |
| `symmetry-preserved/`       | The full report + before/after symmetry table: identical space/layer/point group |
| `ml-relax/`                 | MACE pre-relaxation of the passivated slab, vacuum axis provably fixed |
| `full-report/`              | `--save-report` + `references.bib`                                    |

## Try it yourself

```bash
# Build the bulk crystal and cut a slab first (see 2.8/2.3)
stb-crystalbuilder --spacegroup Fd-3m --a 5.43 --site Si 0 0 0 --reduce conventional -o si_bulk.fdf
stb-slab -f si_bulk.fdf --hkl 1 1 1 -o slab.fdf

# Passivate the dangling bonds, correcting the default bond-length approximation
stb-passivate -f slab.fdf --bond-length 1.48

# Cap with a different element entirely
stb-passivate -f slab.fdf --passivant F --bond-length 1.4

# Let MACE relax the passivated geometry before a real SIESTA calculation
stb-passivate -f slab.fdf --ml-relax --ml-relax-cell
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-f/--file`          | Input structure file, `.fdf` (required)                                |
| `--passivant`        | Element to cap dangling bonds with (default `H`)                       |
| `--cutoff`           | Neighbor-search radius, Ang (default: auto)                            |
| `--bond-length`      | Passivant bond length, Ang (default: auto, sum of atomic radii)        |
| `-sp/--symprec`      | Symmetry tolerance for the before/after table (default `0.01`)         |
| `--ml-relax`         | Pre-relax the passivated structure with a MACE potential               |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the cell (vacuum axis stays fixed)        |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file                              |
| `-o/--output`        | Output `.fdf` file name (default `passivated.fdf`)                     |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the input and passivated structures in ASE's interactive viewer   |

## What's next

See `2.3-stb-slab/` for cutting the slab this tool is usually applied to
(and its own `--passivate` convenience flag, for doing both in one step),
or `2.8-stb-crystalbuilder/` for building the bulk crystal this whole
pipeline starts from.
