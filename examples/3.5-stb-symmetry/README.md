# 3.5 — Crystal Symmetry Analyzer (`stb-symmetry`)

## What this tool does

`stb-symmetry` reads a SIESTA structure file (`.fdf` or a post-relaxation
`.STRUCT_OUT`) and reports its symmetry, with no SIESTA re-run needed:

- The **3D space group**, Wyckoff positions/orbits, site symmetry, and the
  full list of symmetry operations.
- For a **vacuum-padded 2D slab** (exactly one vacuum axis): the **layer
  group** instead — the physically correct 2D classification.
- For a **vacuum-padded isolated molecule** (all 3 axes vacuum-padded): the
  **point group** instead.
- A per-atom **distortion** (Å): how far each atom sits from where the
  detected symmetry operations would place it exactly.
- Optionally: a **tolerance-sensitivity scan** (`--scan-symprec`), a
  **comparison** between two structures (`--compare-to`, e.g. before/after
  a SIESTA relaxation), and a **symmetry-refined structure**
  (`--write-refined`).

## Why this matters (a bit of theory)

### Space groups, Wyckoff positions, and symmetry operations

A crystal's **space group** is the complete set of symmetry operations
(rotations, reflections, screw axes, glide planes, and the translations
that repeat the unit cell) that map the structure onto itself. `stb-symmetry`
detects this via [spglib](https://spglib.readthedocs.io/), wrapped by
pymatgen's `SpacegroupAnalyzer` — the same detector `stb-unitcell` and
`stb-fetch --unitcell` already use in this suite (`core/symmetry.py`).

Atoms related to each other by a symmetry operation belong to the same
**Wyckoff orbit** (e.g. NaCl's 4 Na atoms all share orbit `4a`, its 4 Cl
atoms all share `4b` — 8 atoms, but only 2 physically distinct
environments). Each orbit also has a **site symmetry** — the subgroup of
operations that fix that particular position (e.g. `m-3m` for a
rock-salt ion sitting exactly on a high-symmetry point).

Symmetry **operations** are printed in compact `x,y,z` (Seitz) notation —
e.g. `-x, -y, -z` is an inversion through the origin, `x, y+1/2, z+1/2` is
a translation by half the cell along `b` and `c`.

### The distortion metric — and why it's *not* an absolute measure

Real structures are rarely perfectly symmetric (numerical noise from a
DFT relaxation, or a deliberately displaced defect). `stb-symmetry`
reports a per-atom **distortion** (Å): apply every *detected* symmetry
operation to an atom's position and measure the distance to the nearest
same-species atom that operation should map it onto exactly if the
structure had perfect symmetry; the atom's distortion is the largest such
distance over all non-identity operations.

**Limitation**: this is bounded by construction to roughly `<= --symprec`
(that's the definition of `symprec`: the tolerance within which spglib
accepted this operation set as valid in the first place) — it's a
finer-grained, per-atom breakdown *within* that tolerance budget, not an
absolute "how far from ideal" number independent of what tolerance you
asked for. An earlier attempt used `SpacegroupAnalyzer.get_refined_structure()`
directly, but that method empirically returns the input positions
essentially unchanged (verified: every atom's distance to its own
"refined" counterpart came out exactly 0.0, even for a deliberately noisy
structure) — recomputing the residuals directly from the operations
themselves, as done here, is what actually produces meaningful non-zero
values correlated with the injected noise.

### Why a vacuum-padded structure needs LAYER GROUP or POINT GROUP

spglib's 3D space-group detection has no concept of "vacuum" — it treats
an empty gap along one axis as just an unusually tall periodic cell. For
a slab (one vacuum axis), this means the reported 3D space group actually
depends on the arbitrary vacuum thickness you chose, **not** the real 2D
symmetry of the sheet itself. `stb-symmetry` calls this out explicitly
(see the `[1] DIMENSIONALITY & VACUUM AXES` section) and, when possible,
reports the physically meaningful classification instead:

- **1 vacuum axis** (a slab) → the **layer group**, via spglib's own
  dedicated `get_layergroup()` (needs spglib ≥ 2.1.0; ≥ 2.2.0 for a
  Hall-number bugfix) — a real, separate 2D-symmetry detector, not a
  heuristic filter of the 3D result.
- **3 vacuum axes** (an isolated molecule) → the **point group**, via
  pymatgen's separate, non-periodic `PointGroupAnalyzer` (it treats the
  structure as a plain molecule — Cartesian coordinates only, no
  lattice — the same way a quantum-chemistry package would). Its
  detection tolerance is a molecular, Å-scale parameter, **independent of
  `--symprec`** (a much tighter, crystallographic default) — pymatgen's
  own default (0.3 Å) is used unless overridden by `--scan-symprec`'s own
  point-group sweep.
- **2 vacuum axes** (a wire) → **neither section** — see Limitations
  below.

### `--symprec` / `--angle-tolerance` and `--scan-symprec`

`--symprec` (Å) and `--angle-tolerance` (degrees) are the tolerances
spglib uses to decide whether two atoms/angles are "the same" for
symmetry purposes. A structure fresh out of a DFT relaxation often has
small positional noise that hides its true symmetry at a tight tolerance
— `--scan-symprec` sweeps a range of tolerances (tight to loose) and
reports the first one at which the detected group changes, revealing
symmetry that's really there but momentarily obscured by noise.

### `--write-refined`

Writes a symmetry-refined structure: positions snapped to the detected
symmetry, same reduction `stb-unitcell --mode refined` performs. Refines
against whichever classification is physically meaningful for the
structure's dimensionality (3D space group / layer group / point group).
The file is written **inside `--output-dir`** (a plain filename, not a
path of its own) alongside `references.bib`/the report, and carries a
header documenting exactly how it was generated — source file/format,
`--symprec`/`--angle-tolerance`, and which group (space/layer/point) it
was refined against — so the file is self-describing even opened on its
own, without the console report next to it. **Known caveat** (documented
in `core/symmetry.py`): the refined
structure's atom order and coordinate origin are **never guaranteed** to
match the input, in any mode — spglib/pymatgen rebuild the cell from the
detected symmetry operations from scratch and are free to pick any
symmetry-equivalent origin, which even tiny input noise can flip for
highly symmetric structures. Not a bug; it's the same crystal, just
possibly relabeled.

## Limitations

- **No rod-group detection.** A wire (1D-periodic, 2 vacuum axes) has no
  fallback at all — spglib doesn't implement rod groups (the periodic
  analogue of a layer group for a 1D system). `stb-symmetry` reports this
  explicitly rather than silently falling back to the physically
  meaningless 3D space group.
- **Layer-group detection can fail even with exactly 1 vacuum axis** —
  the structure may not be genuinely 2D-periodic within `--symprec`, or
  the installed spglib may predate `get_layergroup()` (≥ 2.1.0).
- **The point-group tolerance is a separate parameter**, not tied to
  `--symprec` — don't expect `--symprec` to change point-group detection
  sensitivity; use `--scan-symprec` to see the point group's *own*
  tolerance sweep instead.
- **The distortion metric is relative to `--symprec`**, not an absolute
  measure of "how far from an ideal position" (see above).
- **`--write-refined`'s atom order/origin is not preserved** (see above).
- **A single isolated atom** (all 3 axes vacuum-padded, e.g.
  `stb-cohesive`'s own isolated-atom reference structures) is a
  geometrically degenerate case for pymatgen's `PointGroupAnalyzer` — it's
  caught and reported as "detection... failed" rather than crashing the
  whole run, but no point group is available for it either.

## The report: console output, `--save-report`, `--view`

Every run prints a numbered report to the console:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | date/time, input file/format, tolerances, output dir, active options |
| `[1] DIMENSIONALITY & VACUUM AXES` | which axes (if any) are vacuum-padded, and which section below is physically meaningful |
| `[2] SPACE GROUP` | the 3D space group, Hall symbol, point group, crystal system, Pearson symbol |
| `[3] LAYER GROUP` / `[3] POINT GROUP` | (conditional, exactly one or neither) the 2D or 0D classification |
| `[4] LATTICE` | cell parameters, volume, reduced formula, lattice vectors |
| `[5] TOLERANCE SENSITIVITY SCAN` | (conditional, `--scan-symprec`) how the detected group changes with tolerance |
| `[6] SYMMETRICALLY DISTINCT SITES` | the Wyckoff orbits |
| `[7] ATOMIC SITES` | per-atom fractional coordinates + distortion |
| `[8] SYMMETRY OPERATIONS` | (conditional, on unless `--no-operations`) the full x,y,z-notation list |
| `[9] SYMMETRY COMPARISON` | (conditional, `--compare-to`) space/layer/point group comparison between two structures |
| `[10] REFERENCES` | writes `references.bib` (SIESTA) |
| `[11] SUMMARY & FILES` | status and a recap of every file written |

- **`--save-report`** — also persists the full numbered report to
  `stb_symmetry_report.txt`. Off by default, so a plain run only ever
  writes `references.bib` — no text report file at all (the old,
  always-on `symmetry.dat` this tool used to write unconditionally is
  gone; a stale one left by an older run is cleaned up automatically
  instead of being mistaken for current output).
- **`--view`** — opens an interactive 3D view (via ASE) of the analyzed
  structure — or, if `--write-refined`/`--compare-to` is also given, both
  structures side by side (page through frames in ase-gui). Needs a
  display. Off by default.

## When you'd reach for it

- A quick sanity check on a structure's symmetry before committing to a
  real SIESTA run (e.g. confirming a hand-built structure has the space
  group you intended).
- Checking whether symmetry **survived a relaxation** (`--compare-to` a
  pre-relaxation `.fdf` against its post-relaxation `.STRUCT_OUT`).
- Recovering symmetry hidden by small numerical noise (`--scan-symprec`).
- Classifying a 2D slab or an isolated molecule correctly, instead of
  trusting a 3D space group that only describes the padded supercell.
- Generating a clean, symmetry-refined structure (`--write-refined`) for
  a subsequent calculation.

## Two ways to run it

A — direct CLI:
```bash
stb-symmetry --file nacl.fdf --format fdf
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.5
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.5.sh` proves this directly at the end.

## Files in this folder

- `nacl.fdf` — textbook rock-salt NaCl (space group Fm-3m, No. 225).
- `nacl_noisy.fdf` — the same NaCl with a small random atomic displacement,
  for the `--scan-symprec` demonstration.
- `graphene_slab.fdf` — a 2-atom graphene layer, 20 Å vacuum along `c`.
- `molecule.fdf` — an isolated water molecule, 10 Å vacuum on `a`/`b`/`c`.
- `wire.fdf` — a 1D-periodic chain, 10 Å vacuum on `a`/`b` — the
  documented "no rod group" limitation case.
- `example_3.5.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/5-symmetry/test.sh` for that).
- `.gitignore` — excludes `output/` and other generated files.

## Running the walkthrough

```bash
cd examples/3.5-stb-symmetry
./example_3.5.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder               | What it shows                                                              |
|-----------------------|-----------------------------------------------------------------------------|
| `bulk-nacl/`          | Space group, Wyckoff orbits, symmetry operations, zero distortion on textbook NaCl |
| `scan-symprec/`       | `--scan-symprec` recovering symmetry hidden by numerical noise |
| `layer-group/`        | Why a 2D slab needs LAYER GROUP instead of (or in addition to) SPACE GROUP |
| `point-group/`        | An isolated water molecule's POINT GROUP (C2v) |
| `limitation-wire/`    | The documented "no rod group" limitation for a 1D wire |
| `full-report/`        | Default (no report file) vs. `--save-report` (`stb_symmetry_report.txt`), `references.bib` |

## Try it yourself

```bash
# A structure you just built or relaxed
stb-symmetry --file my_calc.fdf --format fdf --save-report

# Check symmetry survived a relaxation
stb-symmetry --file my_calc.fdf --format fdf --compare-to my_relaxed.STRUCT_OUT --compare-format struct_out

# Recover symmetry hidden by small numerical noise
stb-symmetry --file my_relaxed.STRUCT_OUT --format struct_out --scan-symprec
```

## Flag reference

| Flag                 | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--file`             | Path to the structure file.                                            |
| `--format`           | `fdf` (SIESTA input) or `struct_out` (post-relaxation `.STRUCT_OUT`).  |
| `--symprec`          | Symmetry-detection tolerance in Å (default `1e-3`).                    |
| `--angle-tolerance`  | Symmetry-detection angle tolerance in degrees (default `5.0`).         |
| `--no-operations`    | Skip the full symmetry-operations list (keeps just the count).         |
| `--scan-symprec`     | Sweep a range of tolerances and report where the group changes.        |
| `--compare-to`       | Also analyze a second structure and compare symmetry.                  |
| `--compare-format`   | File format of `--compare-to` (required if given).                     |
| `--write-refined`    | Write the symmetry-refined structure as this filename, inside `--output-dir`, with a provenance header. |
| `-o/--output-dir`    | Where to write `references.bib` (and the report, with `--save-report`). |
| `--save-report`      | Persist the full numbered report to `stb_symmetry_report.txt`.         |
| `--view`             | Open an interactive 3D view (via ASE) after the report is written.     |

## What's next

`stb-structural` (3.4) covers a complementary question — a structure's
local coordination environment (bond lengths/angles, effective
coordination number) — rather than its global space/layer/point-group
symmetry. `stb-unitcell` reduces a structure to its primitive/conventional
cell using the same underlying `core/symmetry.py` machinery this tool's
`--write-refined` calls for a genuinely 3D bulk structure.
