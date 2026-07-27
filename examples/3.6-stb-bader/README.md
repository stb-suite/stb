# 3.6 — Bader Charge Analysis (`stb-bader`)

## What this tool does

`stb-bader` computes Bader (Atoms-in-Molecules) atomic charges from a
finished SIESTA calculation, with no extra SIESTA run needed:

- Reads `<label>.RHO` (the real-space charge-density grid) plus
  `<label>.XV` or `<label>.fdf` (the geometry), converts to a `.cube`
  file via [sisl](https://sisl.readthedocs.io/), and partitions the
  density into atomic basins with [PyBader](https://github.com/kerrigoon/pybader).
- Reports each atom's **population** (integrated electrons in its own
  basin), **net charge** (valence − population), and whether it's acting
  as a donor or acceptor.
- Auto-detects each species' **valence electron count** (`Z_val`) from
  `<label>.out` (or `--ref`), falling back to a hardcoded periodic-table
  guess per species when it can't.
- If the `.RHO` file is spin-polarized, also reports each atom's net spin
  (magnetic moment).
- Also reports each atom's **Bader volume** and minimum **surface
  distance**, a **per-species mean/std summary**, and cross-checks
  **symmetry-equivalent atoms** against each other.

## Why this matters (a bit of theory)

### Bader's Atoms-in-Molecules (AIM) theory

There is no operator in quantum mechanics for "how much charge belongs to
this atom" — a molecule or crystal's electron density is one continuous
field, not divided into per-atom pieces by anything in the underlying
physics. Bader's AIM theory (1990) gives a rigorous, *purely
density-based* way to divide it anyway: partition all space into regions
separated by **zero-flux surfaces** — surfaces where the gradient of the
electron density, `∇ρ(r)`, has zero component perpendicular to the
surface at every point. Each such region (a "basin") is assigned to
whichever nucleus it encloses, and integrating the density over a basin
gives that atom's Bader population.

This is attractive specifically *because* it needs nothing extra: no
fitted atomic radii, no basis-set-dependent partitioning (unlike e.g.
Mulliken charges) — just the density itself, which is in principle an
observable quantity. **Limitation**: it's still one specific, if
well-motivated, choice of how to divide up a field that isn't naturally
divided at all — "the charge on this atom" is not a direct experimental
observable, and different partitioning schemes (Bader, Hirshfeld,
Voronoi, ...) can disagree, sometimes substantially, especially for
strongly covalent or highly delocalized bonding.

### The grid-based algorithm — and its own approximation

Finding the *exact* zero-flux surfaces analytically is only tractable for
toy densities. PyBader instead uses the **grid-based** approach of Tang,
Sanville & Henkelman (*J. Phys.: Condens. Matter* **21**, 084204, 2009):
starting from each voxel of the charge-density grid, follow the steepest
ascent of the density (via a weighted average of neighboring voxels) until
you land on a density maximum (a nucleus) — every voxel that flows to the
same maximum belongs to that atom's basin. The "near-grid" refinement
step (this tool's default) corrects a known bias of the naive on-grid
method, where basin edges can align suspiciously with the grid axes
themselves instead of the true (grid-independent) zero-flux surface.

**Limitation**: this is fundamentally a **finite-resolution
approximation** — basin boundaries are only as accurate as the density
grid itself (SIESTA's real-space mesh, set by `MeshCutoff`/`PAO.MeshCutoff`
originally). `--speed fast` switches to the plain on-grid method (faster,
no edge-refinement pass) — useful for a quick look, but explicitly
**not** recommended for a number you intend to quote.

### Z_val detection and the unit-correction factor

A Bader *population* (raw integrated electron count in a basin) only
becomes a *net charge* once you know how many valence electrons that
species' pseudopotential started with (`Z_val` — net charge = `Z_val −
population`). `stb-bader` parses this straight out of the actual SIESTA
`.out` log for the pseudopotential really used (falling back to a
generic periodic-table guess, with a loud warning, for any species it
can't find there — see [Limitations](#limitations)).

Because the whole calculation is on a finite real-space grid, the total
integrated population across every atom can drift slightly from the
theoretical total (`sum(Z_val)`) — usually a small, uniform effect of
mesh resolution. If that drift exceeds 10%, `stb-bader` applies a single
global multiplicative **correction factor** so the reported total matches
the target. **Limitation**: this is one scalar applied uniformly to
every atom — it's the right fix if the deviation really is a uniform
mesh-resolution effect, but if the real cause is localized (one bad atom,
a coarse region of the mesh, or a mismatched `.RHO`/geometry pairing), the
correction just smears that same error evenly across every atom's
reported charge instead of fixing it. A large correction factor is a cue
to investigate further, not a guarantee the correction is physically
right.

### Symmetry cross-check

Atoms the structure's space group treats as equivalent (e.g. every O atom
on the same Wyckoff site) are *physically required* to carry the same
Bader charge. `stb-bader` runs pymatgen's `SpacegroupAnalyzer` on the
geometry (best-effort — silently skipped if symmetry detection fails,
e.g. for an already-relaxed structure with no exact symmetry left) and
flags any symmetry-equivalent group whose charges disagree by more than
`0.1 e-` — a red flag for numerical noise, a resolution problem localized
to part of the cell, or a genuinely mis-set-up structure, invisible from
looking at any single atom's row on its own.

## Limitations

- **Non-nuclear attractors.** A real Bader analysis can occasionally find
  a density maximum that isn't centered on any atom at all (a "non-nuclear
  attractor" — most common in ionic or metallic bonding). PyBader always
  folds such a basin into the nearest atom; `stb-bader` cannot detect or
  correct for this happening.
- **`Z_val` accuracy.** Detected values come from parsing the actual
  pseudopotential generation log in `.out` — accurate for the exact
  pseudopotential used. The hardcoded `FALLBACK_VALENCE` table (used only
  for a species the `.out` parse can't find) is a generic "standard
  valence" guess and can be wrong for a semicore-inclusive pseudopotential
  (e.g. one that treats normally-core `d`/`f` electrons as valence).
- **`--speed fast`** trades basin-boundary accuracy for speed (see above)
  — a quick preview, not a final number.
- **Near-zero population is a red flag, not a clean result.** An atom
  reporting essentially 0 electrons (below `0.01 e-`) usually means
  PyBader found no density feature of its own to anchor a basin on there
  (e.g. a pseudopotential with deep semicore states frozen into the core,
  leaving too little resolvable valence density near that nucleus) and
  its whole region got folded into a neighbor's basin instead — this
  fixture's own Sn atoms #1/#2 demonstrate exactly this (see the
  walkthrough below).
- **The unit-correction factor is a single global scalar** (see above) —
  it can mask a localized problem instead of fixing it.
- **A wrong `.RHO`/geometry pairing isn't always detectable.** The lattice
  and atom-count/order cross-checks catch a genuinely different
  calculation, but neither can catch a same-cell, same-atom-order
  mismatch — a `.RHO` grid carries no independent atomic-position record
  to compare against.

## The report: console output, `--save-report`

Every run prints a numbered report to the console:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | label, output dir, speed mode, `--ref`, vacuum tolerance, threads, cube/export flags |
| `[1] VALENCE (Z_val) SETUP` | Z_val source (detected vs. hardcoded fallback) per species |
| `[2] PYBADER CONFIGURATION` | method, refine method, threads, vacuum tolerance, spin-polarization |
| `[3] PER-ATOM BADER POPULATIONS` | the main per-atom table (population, Z_val, net charge, state, volume/surface distance, spin) |
| `[4] PER-SPECIES SUMMARY` | mean/std net charge per species |
| `[5] DIAGNOSTICS & WARNINGS` | unit-correction factor, near-zero-population atoms, symmetry-inconsistent groups, general limitations |
| `[6] REFERENCES` | writes `references.bib` (SIESTA); names the Bader/PyBader theory papers |
| `[7] SUMMARY & FILES` | status and a recap of every file written |

- **`--save-report`** — also persists the full numbered report to
  `stb_bader_report.txt`. Off by default, so a plain run only ever writes
  the `.cube` file(s) and `references.bib` — no text report file at all
  (the old, always-on `<label>_BADER.txt` this tool used to write
  unconditionally is gone; a stale one left by an older run is cleaned up
  automatically instead of being mistaken for current output).
- **`--output-dir`** — where the `.cube` file(s), `references.bib`, the
  report, and any `--export-volumes` files all land (default: current
  directory).

## When you'd reach for it

- A quick charge-transfer sanity check after a SIESTA SCF run (e.g.
  confirming an oxide's cations/anions show the expected donor/acceptor
  sign and roughly the expected magnitude).
- Comparing charge distribution across a series of related structures
  (a defect, a substitution, before/after adsorption).
- Visualizing individual atomic volumes in VESTA/VMD (`--export-volumes`).
- A magnetic-moment estimate per atom, for free, when the `.RHO` is
  spin-polarized.

## Two ways to run it

A — direct CLI:
```bash
stb-bader --label Sn3O4
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.6
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.6.sh` proves this directly at the end.

## Files in this folder

- `Sn3O4.out`, `Sn3O4.RHO`, `Sn3O4.XV` — a real, finished SIESTA
  calculation on a 14-atom Sn₃O₄ structure (copied from
  `test/3-analysis/6-bader/`). Chosen because it genuinely exercises
  several of the tool's diagnostics at once: two chemical species with
  different bonding character (donor Sn, acceptor O), a real
  near-zero-population case on 2 of the 6 Sn atoms (see the walkthrough),
  and enough symmetry for the equivalent-atom cross-check to run.
- `example_3.6.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/6-bader/test.sh` for that).
- `.gitignore` — excludes generated `.cube` files, `references.bib`,
  the report, and `output/`.

## Running the walkthrough

```bash
cd examples/3.6-stb-bader
./example_3.6.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated (each run
takes a little while — PyBader is partitioning a real, fairly dense
charge-density grid):

| Folder               | What it shows                                                              |
|-----------------------|-----------------------------------------------------------------------------|
| `basic/`              | Z_val detection, the per-atom population/net-charge table, per-species summary |
| `speed-fast/`         | `--speed fast`'s on-grid method vs. the default near-grid method |
| `export-volumes/`     | `--export-volumes`: one `.cube` per atom, for VESTA/VMD |
| `full-report/`        | Default (no report file) vs. `--save-report` (`stb_bader_report.txt`), `references.bib` |

## Try it yourself

```bash
# A finished SIESTA calculation of your own
stb-bader --label my_calc --save-report

# Point at a specific .out if <label>.out isn't the right one
stb-bader --label my_calc --ref relax/my_calc.out --save-report

# Visualize individual atomic volumes in VESTA/VMD
stb-bader --label my_calc --export-volumes -o bader_out
```

## Flag reference

| Flag                 | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-l/--label`         | SystemLabel used in SIESTA.                                            |
| `-o/--output-dir`    | Where the `.cube` file(s)/`references.bib` (and the report, with `--save-report`) land. |
| `--save-report`      | Persist the full numbered report to `stb_bader_report.txt`.            |
| `--ref`              | Path to a specific `.out` file to read `Z_val` from.                   |
| `--speed`            | `normal` (near-grid, default) or `fast` (on-grid, less precise edges). |
| `--threads`          | Worker threads for PyBader.                                            |
| `--vacuum-tol`       | Charge-density threshold below which a voxel is bucketed as vacuum (slabs/wires/molecules). |
| `--no-cube`          | Delete the intermediate `.cube` file(s) after the run.                 |
| `--export-volumes`   | Also write each atom's individual Bader volume as its own `.cube` file. |

## What's next

`stb-structural` (3.4) and `stb-symmetry` (3.5) both analyze a structure
file directly (coordination environment, space/layer/point-group
symmetry); `stb-bader` is this suite's only tool that instead analyzes a
finished **SCF charge density**, giving physically grounded atomic
charges rather than a purely geometric picture.
