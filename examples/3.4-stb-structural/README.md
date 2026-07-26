# 3.4 — Structure Analyzer (`stb-structural`)

## What this tool does

`stb-structural` reads a SIESTA structure file (`.fdf` or a post-relaxation
`.STRUCT_OUT`) and reports a broad set of local and extended structural
properties, with no SIESTA re-run needed:

- Lattice parameters, cell volume, density.
- **Effective Coordination Number (ECN)**, computed by 4 independent
  methods side by side.
- Bond distances and bond angles, pooled per species-pair / ligand-center
  -ligand triplet.
- Coordination-polyhedron **distortion**: bond-length distortion (BLD),
  bond-angle variance (BAV), and polyhedron volume.
- The shortest same-species contact distance in the structure.
- Coordination-polyhedron **connectivity** (corner-/edge-/face-sharing)
  between same-species centers.
- Optionally, the **radial distribution function** g(r).

`--mode mean` analyzes every atom in the structure, aggregated per
species. `--mode list --list 1,4,5` restricts the atom-level analyses
(ECN, distortion) to only the atoms you name — see
["`--mode mean` vs. `--mode list`"](#--mode-mean-vs---mode-list-what-stays-whole-structure)
below for exactly what does and doesn't get scoped down.

## Why this matters (a bit of theory)

There is no single, universally-agreed "coordination number" for a real
(often distorted, sometimes disordered) structure — every method below is
a different **operational definition** of "how many neighbors does this
atom have," built on a different criterion for where a neighbor shell
ends. That's exactly why this tool reports 4 of them side by side instead
of picking one: where they agree, you can trust the number; where they
disagree noticeably, that itself is informative (a strong sign the local
environment is genuinely ambiguous or distorted, not a bug in one of the
methods).

### Effective Coordination Number (ECN) — the 4 methods

All four are computed with **weighted** coordination numbers (a
continuous, "effective" count — e.g. `5.98` rather than a flat integer
`6`) rather than a hard neighbor cutoff, via
[pymatgen](https://pymatgen.org)'s `local_env` module.

- **MinDistNN** (`MinimumDistanceNN`) — the simplest of the four: find the
  closest neighbor at distance `d_min`, then include every neighbor within
  `(1 + tol) * d_min` (default `tol = 0.1`, i.e. within 10% of the closest
  bond). Each included neighbor's weight is `d_min / d_i`, so a neighbor
  right at the closest distance contributes weight 1.0 and one further
  away contributes proportionally less.
  **Limitation**: a single global relative tolerance either lumps a
  legitimately slightly-longer bond into the count or excludes it
  entirely — there's no adaptation to how spread out a given site's own
  bond lengths naturally are.

- **BrunnerNN** (`BrunnerNNRelative`) — sorts *every* neighbor distance up
  to a cutoff, then finds the largest **relative gap** (the biggest ratio
  `d_{i+1} / d_i` between two consecutive sorted distances). That gap is
  taken as the boundary between the first coordination shell and
  everything beyond it; every neighbor closer than the gap is included,
  weighted the same way as MinDistNN (`d_min / d_i`).
  **Limitation**: it assumes a clearly separated shell exists at all. In
  densely/continuously packed environments (common in metals or
  disordered structures) there may be no single obviously-largest gap, in
  which case the chosen cutoff is not physically meaningful, just
  whichever ratio happened to be biggest.

- **EconNN** (Hoppe, 1979 — the original ECoN, "Effective Coordination
  Number") — a smooth, self-consistent formula instead of any hard
  cutoff: each neighbor's weight is `w_i = exp(1 - (d_i / d_avg)^6)`,
  where `d_avg` is itself solved **iteratively** as the weighted mean bond
  length (start from the simple mean, recompute the weights, recompute
  the mean, repeat to convergence). The final `ECoN = sum(w_i)`. Because
  the weight decays smoothly rather than switching on/off at a cutoff, a
  neighbor doesn't have to be "in" or "out" — it can contribute a small
  fraction.
  **Limitation**: designed originally for ionic/ceramic bonding (Hoppe's
  full formula can use "fictive ionic radii" instead of raw bond
  distances); this tool uses raw bond distances only
  (`use_fictive_radius=False`), so the ionic-radius refinement isn't
  applied here.

- **CrystalNN** (Zimmermann & Jain) — the most elaborate of the four: a
  Voronoi tessellation gives each candidate neighbor a solid-angle-based
  weight, which is then refined by a smooth distance cutoff (based on the
  sum of covalent radii) and a preference for larger Pauling
  electronegativity differences (chemically motivated: a cation is more
  likely bonded to an anion than to another cation at a similar
  distance). Its own docstring notes the default parameters were
  benchmarked specifically for inorganic crystal structures.
  **Limitation**: the most "opinionated" of the four — its result depends
  on more internal parameters (electronegativity weighting, distance
  -cutoff smoothing) than a purely geometric method, so it can diverge
  most from the others for chemistries or bonding types (e.g. molecular
  crystals, MOFs) its defaults weren't tuned for.

**JmolNN was deliberately left out.** It was tested and consistently gave
the worst agreement with the actual local geometry of the five candidate
methods — its reference distance comes from a fixed Jmol bonding-radius
lookup table rather than the atom's own closest-neighbor distance, so its
weight isn't on a comparable scale to the other four.

### Bond distances and bond angles

For every atom in scope (all atoms in `mean` mode, just the requested
ones in `list` mode), every neighbor found by `CrystalNN` (unweighted
here — a plain neighbor list, not a coordination-number calculation) at
that site contributes one bond distance and one angle per neighbor pair
at that same center. Results are pooled **by species pair** (bonds) or
**ligand-center-ligand species triplet** (angles) — e.g. every `O-Sn`
bond in the structure is averaged together, every `O-Sn-O` angle
separately from every `Sn-O-Sn` angle — rather than one meaningless
structure-wide average mixing chemically distinct bond types.

In `mean` mode, the same physical bond can be counted from both of its
endpoints' neighbor lists if the neighbor relationship is symmetric
(`n=` in the report makes this visible, so it's never hidden). This
doesn't bias the pooled *average* value — it's the same number counted
twice — but does mean `n` is not simply "number of bonds in the
structure."

### Coordination-polyhedron distortion: BLD, BAV, Volume

- **BLD (Bond-Length Distortion, %)** — Baur's index: the mean absolute
  deviation of a site's own bond lengths from *that site's own* mean bond
  length, as a percentage:
  ```
  BLD = 100 * mean(|d_i - d_avg| / d_avg)
  ```
  **Limitation**: needs at least 2 bonds. It only measures how *unequal*
  the bond lengths are — a site with 6 identical-length bonds arranged in
  a wildly non-octahedral shape (e.g. all pointing into one hemisphere)
  would score BLD = 0%, since BLD says nothing about angles.

- **BAV (Bond-Angle Variance, deg²)** — based on Robinson et al. (1971),
  but **deliberately modified** here in two ways, both because a
  fractional/"effective" coordination number has no single well-defined
  ideal polyhedron to compare against:
  1. Uses the site's own **observed mean angle** as the reference,
     instead of a theoretical ideal angle (90°/109.5°/etc. for an
     idealized octahedron/tetrahedron).
  2. Includes **every** pairwise ligand-center-ligand angle at that site,
     not just the 12 "cis" (~90°) angles Robinson's original octahedral
     formula uses — there's no general way to separate "cis" from
     "trans" pairs outside an assumed 6-coordinate octahedral topology.

  For a real 6-coordinate site, including the ~180° "trans" angles too
  inflates BAV well above textbook octahedral values.
  **Limitation (important)**: because of this, **BAV values from this
  tool are not directly comparable to BAV values quoted in the
  literature** — they're still a valid, self-consistent distortion metric
  for comparing sites/species *within one report*, just not against
  external benchmarks. Needs at least 3 neighbors (2 angles).

- **Volume (Å³)** — the volume of the convex hull (`scipy.spatial.
  ConvexHull`) of the neighbor atoms' Cartesian positions.
  **Limitation**: needs at least 4 non-coplanar neighbors to bound a 3D
  volume (`N/A` otherwise — you'll see this for e.g. a 3-coordinate site
  in the example below). It also says nothing about *where* the central
  atom sits relative to that hull — an atom pushed far off-center inside
  an otherwise-normal-sized polyhedron still reports a normal volume.

### Same-species minimum distance

The single shortest distance between any two atoms of the same species
anywhere in the structure (e.g. the closest `O···O` contact). A quick
packing/steric sanity check — an implausibly short same-species distance
usually flags a bad structure (overlapping atoms) before you waste time
on a real SIESTA run — or, for metals, a genuine metal-metal bonding
distance.
**Limitation**: it's a single number, the global minimum — not a
distribution. One unusually close pair dominates the report; it says
nothing about the *typical* same-species spacing.

### Coordination-polyhedron connectivity

How two same-species coordination polyhedra (e.g. two `SnO_x` polyhedra)
are linked, classified by how many ligand atoms they have in common:

| Shared ligands | Classification  |
|----------------|------------------|
| 1              | corner-sharing   |
| 2              | edge-sharing     |
| 3+             | face-sharing     |

This is a standard way to describe extended structural topology in
crystallography (e.g. rutile's edge-sharing `TiO6` octahedral chains vs.
perovskite's corner-sharing `BO6` network).
**Limitation**: the neighbor list used here comes from a separate,
**unweighted** `CrystalNN` pass — a plain yes/no neighbor list, not the
same *weighted* coordination number reported in the ECN section above.
The same atom can therefore show e.g. `CrystalNN: 4.004` in the ECN table
while its connectivity classification was built from a discrete neighbor
set that isn't identical to "the nearest 4 atoms."

### Radial distribution function g(r)

The pair correlation function: `g(r)` is proportional to the probability
of finding an atom at distance `r` from a reference atom, relative to
what a uniform (ideal-gas-like) density would predict. `g(r) = 1`
everywhere would mean "no structure at all"; the real, sharp peaks you
see are the coordination shells.

Standard number-density normalization, per species pair (A as reference,
B as target):
```
g_AB(r) = <n_AB(r, r+dr)> / (4*pi*r^2*dr*rho_B)
```
averaged over every A-type atom as the reference. The **first peak**'s
position is the average nearest-neighbor A-B distance; **integrating
g(r) under that first peak recovers the coordination number**:
```
N_AB = 4*pi*rho_B * integral( r^2 * g_AB(r) dr )   [from 0 to the shell boundary]
```
This is a genuine independent cross-check against the ECN section above
(same physics, computed a completely different way) — the walkthrough
below does exactly this integral by hand and recovers the known O-Sn
coordination number.
**Limitations**: needs a large enough `--rdf-rmax` (and, implicitly, a
large enough cell/supercell) that periodic-image finite-size effects
don't distort the tail of the curve; the bin width `dr` (fixed at 0.05 Å
here) trades peak resolution against statistical noise; and `g(r)` is
purely radial — it carries no information about **angular** order (that's
what the bond-angle section above is for).

### `--mode mean` vs. `--mode list`: what stays whole-structure

`--mode list --list 1,7` restricts the **ECN** and **distortion**
(BLD/BAV/Volume) sections to only the atoms you name — computing (and
printing) results for the whole structure would be wasted work if you
only care about one or two sites.

The **same-species minimum distance** and **coordination-polyhedron
connectivity** sections, however, are always computed over the *entire*
structure regardless of `--mode`/`--list` — both describe properties of
the extended network (the closest contact anywhere; how polyhedra link
up across the whole cell), not a property of any single atom, so scoping
them down to a couple of named atoms wouldn't make physical sense.

## When you'd reach for it

- A quick coordination-environment sanity check on a structure before
  committing to a real SIESTA run (spotting an accidentally-collapsed or
  overlapping structure, an unphysical bond length, etc.).
- Comparing coordination number/bond length/distortion across several
  candidate structures (e.g. different space-group guesses from
  `stb-crystalcast`, or before/after a relaxation).
- Characterizing one or two specific sites of interest (`--mode list`) —
  e.g. a defect site, a substitution site, a surface atom — without
  wading through the whole structure's statistics.

## Two ways to run it

A — direct CLI:
```bash
stb-structural --file structure.fdf --format fdf --mode mean
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.4
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.4.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[11] SUMMARY & FILES`)
  printed to the console.
- **`rdf.dat`** — the full g(r) curve, unless `--no-rdf` is given.
- **`references.bib`** — SIESTA (every structure file analyzed here is
  SIESTA input/output).

## Optional (off by default)

- **`--save-report`** — also persists the full numbered report to
  `stb_structural_report.txt`. Without it, only `rdf.dat` and
  `references.bib` are written — no text report file at all (the old,
  always-on `structural_information.dat` this tool used to write
  unconditionally is gone; a stale one left by an older run is cleaned up
  automatically instead of being mistaken for current output).
- **`--no-rdf`** — skip the radial distribution function entirely.
- **`--rdf-rmax`** — g(r) cutoff radius in Å (default `10.0`).

## Files in this folder

- `structure.fdf` — a real, small (14-atom, Sn₃O₄-composition) SIESTA
  structure, copied from `test/3-analysis/4-structure/structure.fdf`.
  Small enough to be a lightweight example, but rich enough to exercise
  every section of the report: two chemical species with clearly
  different coordination environments (6-coordinate Sn, 3-coordinate O —
  including the `Volume: N/A` case for a site with fewer than 4
  neighbors), a genuine O-O and Sn-Sn same-species contact, and both
  corner- and edge-sharing polyhedra.
- `example_3.4.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/4-structure/test.sh` for that, which additionally
  covers the post-relaxation `.STRUCT_OUT` format on the same structure).
- `.gitignore` — excludes `output/`.

## Running the walkthrough

```bash
cd examples/3.4-stb-structural
./example_3.4.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder              | What it shows                                                              |
|----------------------|-----------------------------------------------------------------------------|
| `mode-mean/`          | Full-structure analysis: ECN (4 methods) per species, bond distances/angles, distortion, connectivity |
| `mode-list/`          | `--mode list --list 1,7`: ECN/distortion scoped to 2 atoms, same-species distance/connectivity still whole-structure |
| `rdf/`                | g(r), and hand-integrating its first O-Sn peak to independently recover the ECN section's coordination number |
| `full-report/`        | Default (no report file) vs. `--save-report` (`stb_structural_report.txt`), `references.bib` |

## Try it yourself

```bash
# Your own finished (or pre-relaxation) SIESTA structure
stb-structural --file my_calc.fdf --format fdf --mode mean --save-report

# Focus on specific atoms (e.g. a defect site you just built)
stb-structural --file my_calc.fdf --format fdf --mode list --list 12,45
```

## Flag reference

| Flag              | Meaning                                                                |
|--------------------|-------------------------------------------------------------------------|
| `--file`           | Path to the structure file.                                            |
| `--format`         | `fdf` (SIESTA input) or `struct_out` (post-relaxation `.STRUCT_OUT`).  |
| `--mode`           | `mean` (whole structure, per species) or `list` (specific atoms).      |
| `--list`           | Comma-separated 1-based atom indices. Required with `--mode list`.     |
| `-o/--output-dir`  | Where to write `rdf.dat`/`references.bib` (and the report, with `--save-report`). |
| `--no-rdf`         | Skip the radial distribution function (`rdf.dat` not written).         |
| `--rdf-rmax`       | g(r) cutoff radius, Å (default `10.0`).                                |
| `--save-report`    | Persist the full numbered report to `stb_structural_report.txt`.       |

## What's next

The ECN methods here (`pymatgen.analysis.local_env`) are used only by
this tool in the suite — no shared `core/` module for coordination
geometry exists yet (extract-on-second-use, per this project's own
convention, once a second tool needs the same machinery). `stb-symmetry`
(3.5) covers a complementary question — the structure's space-group/
point-group symmetry — rather than its local coordination environment.
