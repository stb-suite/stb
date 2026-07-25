# 2.8 — Crystal Builder (`stb-crystalbuilder`)

## What this tool does

Builds a structure **from scratch**, given only a space group and the
minimal set of symmetrically-distinct atomic sites (one per Wyckoff
position) — pymatgen expands everything else using the space group's own
symmetry operations. Same **Structures** category (2) as
`2.1-stb-2Dstacking/`, ..., `2.7-stb-unitcell/`, but it's the **inverse**
of that last one: `stb-unitcell` takes a full structure and reduces it
down to its primitive/conventional cell; `stb-crystalbuilder` takes the
minimal symmetry description and builds the full structure up.

## Is it only for 3D systems? **Yes — and this matters.**

`--spacegroup` accepts one of the **230 international space groups**,
which are the complete classification of symmetry groups for **3D
-periodic crystals only**. There is no 2D (layer group) or 1D (rod group)
analog here — pymatgen's `Structure.from_spacegroup` (what this tool
wraps) is inherently a 3D-crystal builder, and every run's own `[2]
STRUCTURE VALIDATION` section now prints this explicitly:

```
Dimensionality : 3D (bulk material) -- every space group is inherently
3D-periodic; use stb-slab/stb-nanotube/stb-2Dstacking on the output to cut
a 2D/1D structure from it.
```

If you want a slab, a nanotube, a monolayer, or a wire, the workflow is
**two steps, not one**: build the 3D bulk crystal here, then cut it down
with the tool that actually understands vacuum/lower dimensionality:

| Want...                     | Build with `stb-crystalbuilder`, then cut with... |
|------------------------------|----------------------------------------------------|
| A surface slab               | `stb-slab` (`2.3`)                                  |
| An isolated 2D monolayer      | `stb-slab` with a thin enough `--min-slab-size`, or `stb-2Dstacking` (`2.1`) if you're also stacking two layers |
| A nanotube or nanoribbon      | `stb-nanotube` (`2.4`), which itself takes a 2D monolayer as input |

`bulk-graphite-to-slab/` below proves this concretely: it builds real
bulk graphite with `stb-crystalbuilder` (space group `P6_3/mmc`), then
pipes the result straight into `stb-slab` to cut a 2D-periodic slab out
of it — `stb-slab`'s own before/after table shows `Layer Group` going
from `N/A (not 2D-periodic)` on the bulk input to a real, detected
`p-3m1` on the cut slab, confirming the bulk crystal really was 3D-only.
(`stb-crystalbuilder`'s own report never shows a `Layer Group` row at
all — every space group it builds from is inherently 3D-periodic, so it
would always read `N/A`, never actual information.)

## Why this matters (a bit of theory)

### Space groups and Wyckoff positions

A crystal's full symmetry is one of exactly 230 **space groups** — every
combination of translations, rotations, mirrors, and screw/glide
operations that maps an infinite 3D lattice onto itself. Rather than
listing every atom's position by hand, a crystallographer only needs to
give the **minimal, symmetrically-distinct set** of atoms (their **Wyckoff
positions**) — the space group's own operations generate every other atom
in the cell automatically. This is exactly how real crystal structures
are published and stored (crystallographic databases like ICSD/COD list a
space group + a handful of Wyckoff sites, not a giant explicit atom list)
— `--site` mirrors that convention directly: give one entry per Wyckoff
letter, not every atom.

### Special positions: you can accidentally ask for *more* symmetry than requested

A Wyckoff site can land on a **general position** (no extra local
symmetry) or a **special position** (sitting exactly on a mirror, axis, or
inversion center of the space group). If your lattice/site combination
happens to have *even more* symmetry than the space group you asked for
allows, the actual structure snaps to that higher symmetry — the tool
detects and reports this explicitly rather than silently building
something subtly different from what you requested. `special-position
-polonium/` below shows a real, physical example of exactly this.

### Reduction after construction: primitive vs. conventional, revisited

Space-group construction naturally produces the **conventional** cell
(the one an international-tables/database entry actually lists) — often
larger than strictly necessary for a DFT calculation. `--reduce
{primitive,conventional,refined}` (new in this version) folds in
`stb-unitcell`'s own reduction directly, so you can go from "space group +
Wyckoff sites" straight to "the cheapest cell for a real calculation" in
one command — see `2.7-stb-unitcell/` for the full theory of what each of
these three modes means.

## When you'd reach for it

- You know a crystal's space group and Wyckoff positions (from a paper,
  International Tables, or a database entry) and want to build the actual
  atomic structure rather than typing every position by hand.
- Constructing a hypothetical/idealized structure (a specific
  substitution pattern, a textbook example) that isn't sitting in any
  database as a downloadable file.
- The starting point for a lower-dimensional structure: build the 3D bulk
  here, then hand it to `stb-slab`/`stb-nanotube`/`stb-2Dstacking` (see
  above).

## Two ways to run it

A — direct CLI:
```bash
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.8
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.8.sh` proves this directly at the end. Path B prints
the bulk-only warning immediately, before asking for the space group, so
it's the first thing you see even if you never read `--help`.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[10] SUMMARY & FILES` —
  `[4]` only with `--reduce`, `[5]` only with `--ml-relax`) printed to the
  console and, with `--save-report`, also saved to
  `stb_crystalbuilder_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density, plus the dimensionality note above — run once on the raw
  space-group-expanded structure and once on the final structure (after
  any `--reduce`/`--ml-relax`).
- **Full symmetry detection** — space group, point group, crystal system,
  Hall symbol — with an explicit `[NOTE]` whenever the detected space
  group differs from the one you requested (see "special positions"
  above).
- **A before/after symmetry comparison table** — crystal system, space
  group, point group, Hall symbol — expected (and verified) to match
  exactly, since reduction/relaxation never changes the crystal's actual
  symmetry. No `Layer Group` row here (unlike `stb-unitcell`/`stb-slab`):
  every space group this tool builds from is inherently 3D-periodic, so
  it would always read `N/A (not 2D-periodic)`, never real information.
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax` was used.
- **A provenance header** written into the output `.fdf`: requested vs.
  detected space group, lattice, sites, reduce mode and MACE convergence
  detail if used.

## Optional (off by default)

- **`--reduce {primitive,conventional,refined}`** — reduce the built
  structure to a smaller/standardized/symmetry-cleaned cell before
  writing (reuses `stb-unitcell`'s own reduction).
- **`--symprec`** / **`--angle-tolerance`** — the symmetry-detection
  tolerance, used for both the post-build symmetry check and `--reduce`.
- **`--ml-relax`** (+ `--ml-relax-cell`, `--model`/`--custom-model`) —
  pre-relax the final structure with MACE. Since every space group is
  inherently 3D-periodic, `--ml-relax-cell` here always relaxes all 3
  cell axes — no vacuum-axis masking is needed (unlike `stb-unitcell`,
  which must handle arbitrary user-supplied slabs/wires too).
- **`--save-report`** — also persist the full report to
  `stb_crystalbuilder_report.txt`.
- **`--view`** — opens the raw built structure and the final structure in
  ASE's interactive 3D viewer (`ase-gui`), as pageable frames. Needs a
  local display. Never exercised by `example_2.8.sh` itself.

## Files in this folder

- `example_2.8.sh` — the guided walkthrough (see below). No input
  fixtures are needed — unlike most `2-structures` tools,
  `stb-crystalbuilder` builds everything from `--spacegroup`/`--site`
  flags, not from a structure file.
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_crystalbuilder_report.txt`.

## Running the walkthrough

```bash
cd examples/2.8-stb-crystalbuilder
./example_2.8.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Nine self-contained cases are generated (one
is skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                        | What it shows                                                          |
|--------------------------------|-------------------------------------------------------------------------|
| `fcc-nickel/`                  | The basic case: 1 Wyckoff site -> 4 atoms, symbol/number equivalence   |
| `magnetite-spinel/`             | A real oxide with 2 distinct Fe sites + 1 O site -> 56 atoms           |
| `special-position-polonium/`    | Requesting less symmetry than the sites actually have -> the `[NOTE]`  |
| `reduce-to-primitive/`          | `--reduce primitive` on magnetite: 56 -> 14 atoms, one command         |
| `out-of-range-error/`           | The fixed bug: `--spacegroup 255` now gives a clear error message      |
| `ml-relax/`                     | MACE pre-relaxation of the built + reduced structure                   |
| `full-report/`                  | `--save-report` + the validation checklist + `references.bib`          |
| `bulk-graphite-to-slab/`        | **The 3D-only answer, proven**: bulk graphite -> a real 2D slab via `stb-slab` |

## Try it yourself

```bash
# Build a structure from a space group you found in a paper or database
stb-crystalbuilder --spacegroup Fd-3m --a 8.396 \
    --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 --site O 0.379 0.379 0.379

# Go straight to the cheapest cell for a real DFT calculation
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 --reduce primitive

# Build a 3D bulk crystal, then cut a 2D slab out of it (the 3D-only workaround)
stb-crystalbuilder --spacegroup "P6_3/mmc" --a 2.464 --c 6.711 --gamma 120 \
    --site C 0 0 0.25 --site C 0.333333 0.666667 0.25 -o graphite.fdf
stb-slab -f graphite.fdf --hkl 0 0 1 --min-slab-size 3.0 --min-vacuum-size 15
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--spacegroup`       | International symbol (e.g. `Fm-3m`) or number 1-230 (required)         |
| `--a/--b/--c`        | Lattice constants, Ang (`--b`/`--c` default to `--a`) (required `--a`) |
| `--alpha/--beta/--gamma` | Lattice angles, degrees (default 90)                                |
| `--site`             | `SYMBOL X Y Z`, repeatable, one per Wyckoff position (required)        |
| `--symprec`          | Symmetry precision, default `1e-3` (used for detection and `--reduce`) |
| `--angle-tolerance`  | Symmetry angle tolerance in degrees, default `5.0`                     |
| `--reduce`           | `primitive`/`conventional`/`refined`, off (no reduction) by default    |
| `--ml-relax`         | Pre-relax the final structure with a MACE potential                    |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the cell (all 3 axes)                    |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file                              |
| `-o/--output`        | Output `.fdf` file name (default `crystal.fdf`)                        |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the built and final structures in ASE's interactive viewer        |

## What's next

See `2.7-stb-unitcell/` for the inverse operation (reducing a full
structure down) and the full theory of primitive/conventional/refined
cells that `--reduce` reuses here, or `2.3-stb-slab`/`2.4-stb-nanotube` for
what to do next with a bulk crystal built here when you actually need a
lower-dimensional structure.
