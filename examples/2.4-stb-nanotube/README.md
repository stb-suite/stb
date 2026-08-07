# 2.4 — Nanotube/Nanoribbon Builder (`stb-nanotube`)

## What this tool does

Rolls a 2D monolayer into a nanotube (a cylinder) or cuts it into a
nanoribbon (a finite-width flat strip), given a chirality `(n, m)`. Same
**Structures** category (2) as `2.1-stb-2Dstacking/`, `2.2-stb-supercell/`,
`2.3-stb-slab/` — everything here builds, generates, or transforms a
structure file. Unlike those three, `stb-nanotube` **reduces**
dimensionality: a 2D monolayer (1 vacuum axis) always becomes a genuinely
1D-periodic tube or ribbon (2 vacuum axes) — see the symmetry table below
for what that means for the reported Layer Group.

Deliberately **not hexagonal-only**: the classic carbon-nanotube
`(n, m)` chirality math (gcd/`(n-m)%3d` shortcuts) only works for a
honeycomb lattice. `stb-nanotube` instead searches directly over the
monolayer's own lattice vectors, so it works for graphene, but equally for
h-BN, a buckled/corrugated sheet (silicene, phosphorene, ...), or any other
2D material — see `corrugated/` below.

## Why this matters (a bit of theory)

### The chiral vector and the translation vector

A chirality `(n, m)` defines the **chiral vector** `C = n*a1 + m*a2`, using
the monolayer's own in-plane lattice vectors. For **tube** mode, `C`
becomes the tube's circumference (rolled so that walking all the way around
`C` brings you back to an equivalent lattice point) — the diameter follows
directly, `2R = |C| / pi`. For **ribbon** mode, `C` instead sets the finite
**width** direction (tiled `--repeats` times), with no wrapping.

Both modes also need a second vector, the **translation vector `T`**, along
which the result is periodic (tube axis, or ribbon length). `T` must be a
real lattice vector of the monolayer, not parallel to `C` — but "not
parallel" alone isn't enough for a *tube*: `roll_to_tube` maps the flat
`(frac_C, frac_T)` cell onto a cylinder treating `C` and `T` as an
**orthogonal** (circumference, axial) frame. That only preserves real bond
lengths/connectivity if `T` genuinely has no component along `C` — i.e. `T`
must be as close to **perpendicular** to `C` as an integer lattice vector
allows, not merely non-parallel to it.

### A real, serious bug, found and fixed: `T` needs to be perpendicular, not just non-parallel

This is not a hypothetical concern — it was a genuine, verified bug in this
tool before this session, silently producing physically wrong tubes for
almost every chirality except plain zigzag. Verified live against
`ase.build.nanotube` (an independent, trusted reference implementation):
the old "just pick the shortest non-parallel `T`" logic gave graphene
`(6, 0)` a 12-atom cell where every carbon atom had only **1** real bonded
neighbor within 1.8 Ang, instead of the correct 3. The fix: every candidate
`T` is first reduced modulo `C` (subtracting the integer multiple of `C`
that minimizes its remaining along-`C` component — `T' = T + k*C` always
generates the *same* periodic cell, since `C x C = 0`), then candidates are
ranked by how close to perpendicular they land, length only as the
tie-break. `basic-tube/` below reproduces the corrected, physically
verified result directly.

### Zigzag, armchair, chiral — and a graphene-specific electronic hint

For a genuinely **hexagonal** lattice only (equal-length vectors, 60/120
degrees apart — checked automatically, never assumed), `(n, m)` also gets
the classic carbon-nanotube labels: **zigzag** (`m = 0`), **armchair**
(`n = m`), or **chiral** (anything else). For a single-species hexagonal
lattice with a 2-atom basis (i.e. graphene itself, not e.g. h-BN, where the
two sublattices are chemically different), the well-known rule
`metallic if (n - m) mod 3 == 0, else semiconducting` is reported too — see
`electronic-hint/` below, verified against several chiralities.

### `gcd(n, m) > 1`: a physically correct cell, just not always the smallest one

A related, subtler point, also found and documented while fixing the bug
above: for a chirality where `gcd(n, m) > 1`, the shortest *exactly
perpendicular* `T` this tool finds can require **more** repeating cells
than the classic literature (Dresselhaus/Saito) CNT unit cell — because
that classic construction additionally allows a combined
**rotation+translation** ("screw") symmetry between repeats, which
`roll_to_tube`/`build_ribbon` deliberately never use (screw symmetry is a
hexagonal-lattice-specific trick, and this tool works for any 2D lattice).
The cell `stb-nanotube` builds is always physically correct — fully,
uniformly bonded, verified live — just occasionally larger than the
absolute minimum. `basic-tube/` shows the tool's own in-report note about
this directly (graphene `(6, 0)` already has `gcd = 6`).

### Curvature strain: bond lengths compress on small-diameter tubes

Rolling a flat sheet without stretching along the circumference still
changes real 3D bond lengths: `curvature-trend/` below shows this directly,
measured, across increasing diameters — bond length rises monotonically
from a compressed ~1.39 Ang at `(4, 0)` toward flat graphene's real 1.42
Ang as the tube gets wider. `ml-relax/` then shows this strain is real and
relaxable: MACE finds a substantial energy drop and a large force-reduction
purely from letting the initial, curvature-compressed geometry (and the
axial period itself) relax.

### Ribbons have real, physical edges

A **ribbon** (finite in the `C` direction, periodic only along `T`) is a
genuinely different kind of object from a tube: it has two real physical
edges. A **tube**, being both circumferentially closed *and* axially
periodic (infinite in the idealized periodic cell), has none at all — every
atom already has its full coordination by construction. `passivate/` below
verifies both sides of this directly: 0 dangling bonds on a tube, real ones
correctly capped with H on a ribbon's two edges.

## When you'd reach for it

- Building a carbon nanotube (or a nanotube of any other 2D material) at a
  specific chirality for electronic-structure calculations.
- Building a graphene/h-BN/TMD nanoribbon, with a specific edge
  termination (`--passivate`) for a SIESTA calculation.
- Comparing metallic vs. semiconducting chiralities before committing to
  expensive DFT.
- Getting a MACE-pre-relaxed starting geometry for a small-diameter,
  strongly-curved tube before a real SIESTA relaxation.

## Two ways to run it

A — direct CLI:
```bash
stb-nanotube -f graphene.fdf --chirality 6 0
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.4
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.4.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[9] SUMMARY & FILES` — `[4]`
  only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_nanotube_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density — run once on the input monolayer and once on the output
  tube/ribbon.
- **A before/after symmetry comparison table** — crystal system, 3D space
  group, layer group, point group, Hall symbol. The monolayer (before)
  always has a real Layer Group; the tube/ribbon (after) never does (it's
  genuinely 1D-periodic, and spglib has no 1D rod-group classification) —
  the opposite contrast from `2.3-stb-slab/`, where it's the *bulk* input
  that lacks one.
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax` was used.
- **A provenance header** written into the output `.fdf`: chirality, mode,
  translation vector, `N_cells`, geometry (diameter/width, axial period),
  and passivation/MACE detail if used.

## Optional (off by default)

- **`--mode ribbon`** (default `tube`) — cut a finite-width strip instead
  of rolling a cylinder.
- **`--repeats N`** — axial repeats (tube) or width repeats along `C`
  (ribbon). Default 1.
- **`--passivate`** (+ `--passivant`/`--cutoff`/`--bond-length`) — cap
  dangling bonds. See `passivate/` below.
- **`--ml-relax`** (needs the optional `ml` extra: `pip install
  stb_suite[ml]`) — pre-relaxes the built tube/ribbon with a MACE
  potential. Positions only by default; add **`--ml-relax-cell`** to also
  relax the axial period — the vacuum-padded width/thickness axes always
  stay exactly fixed. `--model small/medium/large` (default `small`) or
  `--custom-model PATH` picks which MACE potential to use.
- **`--save-report`** — also persists the full report to
  `stb_nanotube_report.txt`.
- **`--view`** — opens the input monolayer and the final tube/ribbon in
  ASE's interactive 3D viewer (`ase-gui`), as pageable frames. Needs a
  local display. Never exercised by `example_2.4.sh` itself.

## Files in this folder

- `graphene.fdf` — a primitive 2-atom hexagonal graphene monolayer
  (`a` = 2.46 Å, vacuum along `c`) — the main fixture for every case below.
- `silicene_buckled.fdf` — a synthetic buckled honeycomb monolayer
  (`a` = 3.87 Å, 0.22 Å buckling between the 2 basis atoms) — for the
  `corrugated/` case.
- `example_2.4.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_nanotube_report.txt`.

## Running the walkthrough

```bash
cd examples/2.4-stb-nanotube
./example_2.4.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Nine self-contained cases are generated (two
are skipped, with a hint, if the optional `ml` extra isn't installed):

| Folder                | What it shows                                                          |
|------------------------|--------------------------------------------------------------------------|
| `basic-tube/`          | Zigzag `(6,0)` tube — full report, the corrected perpendicular-`T` fix, and the `gcd(n,m)>1` note |
| `electronic-hint/`     | Metallic `(6,0)` vs. semiconducting `(7,0)` vs. always-metallic armchair `(6,6)` |
| `chiral-index/`        | A genuinely chiral `(7,1)` tube (`gcd=1`, no cell-size caveat) — chiral angle, `CNT type: chiral` |
| `curvature-trend/`     | Measured bond-length compression vs. diameter, `(4,0)` through `(20,0)` |
| `ribbon/`              | Ribbon mode — finite width scaling linearly with `--repeats`           |
| `passivate/`           | A closed tube has 0 dangling bonds; a ribbon's two edges genuinely do  |
| `ml-relax/`            | MACE pre-relaxation measurably releasing curvature strain               |
| `corrugated/`          | A buckled (silicene-like) sheet rolls into a genuinely corrugated tube |
| `full-report/`         | `--save-report` + the validation checklist + `references.bib`          |

## Try it yourself

```bash
# Show the input monolayer and the final tube side by side
stb-nanotube -f graphene.fdf --chirality 6 0 --view

# A wide, nearly-strain-free tube
stb-nanotube -f graphene.fdf --chirality 20 0

# Pre-relax with MACE, cell included, with your own fine-tuned model
stb-nanotube -f graphene.fdf --chirality 6 0 --ml-relax --ml-relax-cell --custom-model my_finetuned.model

# A passivated nanoribbon, ready for a SIESTA edge-state calculation
stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 6 --passivate
```

## Flag reference

| Flag                | Meaning                                                                |
|----------------------|-------------------------------------------------------------------------|
| `-f/--file`          | Input 2D monolayer structure file, `.fdf` (required)                   |
| `--chirality N M`    | Chirality indices `(n, m)` (required)                                  |
| `--mode`             | `tube` (default) or `ribbon`                                           |
| `--repeats`          | Axial repeats (tube) / width repeats (ribbon), default `1`             |
| `--vacuum-gap`       | Gap, Å, to detect the input monolayer's vacuum axis (default `10.0`)   |
| `--min-vacuum-size`  | Vacuum padding around the output tube/ribbon, Å (default `15.0`)       |
| `--lattice-tol`      | Fractional-coordinate tolerance for cell membership (default `1e-6`)   |
| `--search-bound`     | Override the search range for `T` (default: auto)                      |
| `--symprec`          | Symmetry tolerance, Å, for the before/after table (default `0.01`)     |
| `--passivate`        | Cap dangling bonds on the built tube/ribbon                            |
| `--passivant`        | Element to passivate with, only with `--passivate` (default `H`)       |
| `--cutoff`           | Neighbor-search radius, Å, only with `--passivate` (default: auto)     |
| `--bond-length`      | Passivant bond length, Å, only with `--passivate` (default: auto)      |
| `-o/--output`        | Output `.fdf` file name (default `nanotube.fdf`)                       |
| `--ml-relax`         | Pre-relax the built tube/ribbon with a MACE potential before writing   |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the axial period (vacuum axes stay fixed) |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file for `--ml-relax`             |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the input monolayer and final tube/ribbon in ASE's interactive viewer |

## What's next

See `2.3-stb-slab/` for a different kind of cut (a finite Miller-index
slab instead of a rolled/tiled monolayer), or `1.7-stb-mlrelax/` for a
closer look at the MACE pre-relaxation `--ml-relax` reuses here.
