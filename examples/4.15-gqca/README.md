# 4.15 — Workflow: GQCA, Generalized Quasi-Chemical Approximation (`stb-gqca` / `stb-gqcaAnalysis`)

A 2-stage workflow that turns **3 relaxed DFT structures** into the full
mixing-thermodynamics surface of a substitutional alloy A(1-x)Bx — mixing
enthalpy/entropy/free energy and a short-range-order (SRO) parameter, at
*any* composition `x` and temperature `T`, on *any* lattice (2D or 3D).
This example builds a real, textbook semiconductor alloy (GaAs → AlGaAs)
through Stage 1 for real, then works through Stage 2 on hand-fabricated
(clearly labeled) cluster energies chosen to reproduce, and then
deliberately contradict, a well-known real-material fact — so you can see
exactly what the tool reports in both cases before trusting it on your
own alloy.

## 1. Theory

### 1.1 The pair-cluster (N=2) GQCA and the mass-action equilibrium

The Generalized Quasi-Chemical Approximation (GQCA — Sher, van
Schilfgaarde, Chen, Chen, *Phys. Rev. B* **36**, 4279, 1987; the classic
Guggenheim quasi-chemical approximation applied to a semiconductor alloy
sublattice) treats a disordered A(1-x)Bx sublattice as a lattice of
*independent* nearest-neighbor pairs, each in one of 3 local states:

```
n_j = 0  (AA)      n_j = 1  (AB)      n_j = 2  (BB)
```

with degeneracy `g = [1, 2, 1]` (an AB pair can be read in 2 orders, AA/BB
in 1). Their equilibrium populations `x_j(x, T)` follow a mass-action law
— literally the same algebra as a chemical equilibrium `2(AB) <-> AA +
BB` — driven by each pair type's own excess energy `Delta_j` (relative to
a linear interpolation between the pure end-members). Solve that
equilibrium once you have the 3 `Delta_j`, and you have `x_j` at *any*
composition/temperature — no per-composition DFT run needed.

### 1.2 Why exactly 3 structures, always at 50:50

A 2-site cluster has exactly 3 possible local compositions — there is no
"choose your x" knob at this stage. `stb-gqca` (Stage 1) builds:

- `cluster_n0/` — the pure A end-member (trivial: every sublattice site
  set to A).
- `cluster_n2/` — the pure B end-member (same, set to B).
- `cluster_n1/` — the ordered 50:50 structure whose *local pair
  environment* most purely realizes "every neighbor pair is unlike" (AB)
  — **not** a random 50:50 mixture, whose pairs would average to a mix of
  AA/AB/BB and contaminate `Delta_1` with disorder that isn't there in
  the idealized pair-cluster model.

Your real target composition (e.g. `Al(0.3)Ga(0.7)As`) is swept entirely
in Stage 2 by combining these 3 fixed reference energies — they are the
complete, sufficient input for the whole `x`/`T` surface, not a coarse
sample of it.

### 1.3 Mixing enthalpy, entropy, free energy, and SRO

Once `x_j(x,T)` is known: `H_mix = sum_j x_j * Delta_j`, `S_mix = -kB *
sum_j x_j*ln(x_j/g_j)` (the "gas of independent clusters" entropy),
`G_mix = H_mix - T*S_mix`, and an SRO parameter measuring how far `x_1`
(the AB/mixed-pair fraction) sits from its *random*-mixing value at the
same `x`. A **negative** `Delta_1`(AB pairing favored) pushes `x_1`
*above* random and gives an ordering-tendency (negative-by-convention)
SRO; a **positive** `Delta_1` pushes it below and favors clustering/
unmixing instead.

### 1.4 The convexity limitation: SRO yes, phase separation no

`stb-gqcaAnalysis` also runs a discrete-grid common-tangent (miscibility
-gap) search on `G_mix(x)`. **It essentially never finds one.** Within
this independent-cluster formalism, `G_mix(x)` is *provably convex* for
any `Delta_j` — a convex-analysis theorem, not a numerical accident (see
`core/gqca_solver.py::common_tangent_gaps`'s own docstring for the
argument). This model captures short-range order, **not** genuine phase
separation, which needs cluster-*cluster* correlations beyond this
independent-pair treatment. A "no gap" result is an intrinsic property of
the model — never report it as a physical miscibility claim about the
real material.

### 1.5 Frustrated sublattices: FCC/HCP/triangular can't reach −1.0

`cluster_n1`'s "every pair unlike" target is a pair-correlation value of
exactly −1.0. A **bipartite** (2-colorable) sublattice — simple cubic,
honeycomb — reaches it exactly. An **FCC** sublattice cannot: FCC's
close-packed planes contain triangles (an odd cycle), so no 2-coloring
makes every nearest-neighbor pair unlike — the identical geometric
frustration behind the classic FCC antiferromagnetic Ising model never
reaching perfect Néel order. Both **rocksalt's and zincblende's cation
sublattices are FCC** — this example's own GaAs structure hits exactly
this limit (Section 3.4). `stb-gqca` reports the best achievable ordering
at the smallest exact-composition cell and flags the deviation; it is a
real, unavoidable limit of the pair-cluster approximation on such
lattices, not a bug.

### 1.6 What's confirmed vs. provisional

The mixing-enthalpy and mass-action-law equations are confirmed against
a secondary source restating the original 1987 formalism. The
mixing-entropy formula and the SRO-parameter normalization are tagged
`[UNVERIFIED]` in the tool's own output — a fully-pinned closed form for
both could not be confirmed against the primary literature during this
tool's development (see `core/gqca_solver.py`'s docstrings). Every
number derived from an `[UNVERIFIED]` quantity carries that caveat
forward; treat `S_mix`/`SRO` as qualitatively meaningful, not exact.

## 2. Libraries and external dependencies used

- **`icet`** — `ClusterSpace`/`enumerate_structures` (Stage 1's pair-orbit
  detection and the Hart-Forcade search for the best AB-ordered
  structure at the smallest exact-composition cell).
- **`pymatgen`** — periodic neighbor search (sublattice nearest-neighbor
  distance) and the ASE/pymatgen structure adaptors icet needs.
- **`numpy`** — the mass-action solver and mixing-thermodynamics formulas
  (Stage 2).
- No network access, no ML/MACE dependency anywhere in this workflow —
  every stage is a pure, local, deterministic file operation. Stage 1
  needs no SIESTA either (it only *writes* input files); only reading
  Stage 1's own output back in Stage 2 needs a real `calc.out`.

## 3. Stage 1: cluster structure generation (`stb-gqca`, code `4.15.1`)

### 3.1 What it does

Given a structure (`structure.fdf`, 2D or 3D) and a `calc.fdf` template,
disorders one existing sublattice (`--sublattice`) with a second species
(`--species-b`) and writes `cluster_n{0,1,2}/`, each with its own
`structure.fdf` + `calc.fdf` (relaxation directives forced, see Section
2 of `stb-gqca --help`) + copied pseudopotentials + a density-based
k-grid — ready to relax with SIESTA, never run for you here.

### 3.2 The real system: GaAs → Al(x)Ga(1-x)As

`structure.fdf` is GaAs in the zincblende structure (conventional cubic
cell, `a = 5.6533 Ang`, the real experimental lattice constant).
`--sublattice Ga --species-b Al` disorders the **cation** sublattice,
modeling `Al(x)Ga(1-x)As` — the most widely used III-V heterostructure
alloy (laser diodes, HEMTs, solar cells) and, not coincidentally, the
textbook example of a **nearly ideal solid solution**: AlAs and GaAs
share almost exactly the same lattice constant (~0.1% mismatch), so real
AlGaAs mixes close to randomly at essentially any composition and growth
temperature. Section 5 checks whether Stage 2 reproduces exactly that
qualitative fact — and what it reports if you feed it energies for a
*hypothetically* non-ideal version of the same lattice instead.

### 3.3 Live output: pair orbit, cluster folders

```
$ stb-gqca -f structure.fdf --sublattice Ga --species-b Al -c calc.fdf -O output/algaas

[2] PAIR ORBIT
Shortest sublattice-sublattice distance : 3.9975 Ang
Orbit radius    : 1.9987 Ang
Multiplicity    : 6

[3] CLUSTER STRUCTURES
  [OK] cluster_n0 (8 atoms, k-grid [28, 28, 28], 3D full-cell relaxation)
  [OK] cluster_n2 (8 atoms, k-grid [28, 28, 28], 3D full-cell relaxation)
  [OK] cluster_n1 (4 atoms, k-grid [40, 40, 28], 3D full-cell relaxation, AB-ordered)
```

`cluster_n1` has *fewer* atoms than `cluster_n0`/`cluster_n2` — it comes
from `icet`'s own smallest-exact-composition primitive-cell enumeration,
not the 8-atom input cell. All 3 folders are written for a **real**
relaxation (`MD.TypeOfRun CG`, plus `MD.VariableCell T` here since this
input is 3D bulk) — the 3 structures generally have different equilibrium
lattice constants, and evaluating them at a fixed cell would contaminate
the mixing energy with spurious lattice-mismatch strain.

### 3.4 The FCC-frustration NOTE, live

```
[NOTE] The best achievable AB-ordered pair correlation at cell size
2x primitive deviates by 0.6667 from the ideal -1.0 (every-pair-unlike)
extreme -- expected on a geometrically frustrated sublattice (e.g.
FCC/HCP/triangular; see build_ab_structure's own docstring), not
necessarily a bug. Evaluated 2 candidate(s).

[NOTE] Delta_1 sensitivity check: a 4x-primitive cell does not
meaningfully improve the AB ordering (deviation 0.6667 vs. 0.6667 at
the size actually used, from 5 probed candidate(s)) -- the deviation
at the smallest cell already appears representative, not an artifact
of under-sizing.
```

Real, expected physics (Section 1.5) — the Ga sublattice is FCC. The
sensitivity check confirms this isn't just an under-sized search: a
larger cell doesn't meaningfully improve it either. `Delta_1` (and every
downstream `x_1`/`H_mix`/`SRO` value derived from it) is therefore a
somewhat less clean estimate than a bipartite sublattice would give —
qualitatively meaningful, not a rigorous quantitative bound. This is
carried forward as a `[NOTE]` in Stage 2's own report too.

### 3.5 Report structure

`gqca_stage1.txt` (also printed to the console) has 4 numbered sections
(`[0] RUN METADATA`, `[1] DIMENSIONALITY`, `[2] PAIR ORBIT`, `[3] CLUSTER
STRUCTURES`) plus a `[4] SUMMARY & NEXT STEPS`, and ends with 3 machine
-readable lines (`Sublattice`/`Species B`/`AB deviation`) that
`stb-gqcaAnalysis` reads back automatically — you never re-type
`--sublattice`/`--species-b` in Stage 2.

### 3.6 Running it both ways

```bash
stb-gqca -f structure.fdf --sublattice Ga --species-b Al -c calc.fdf -O gqca_study
```

or, at the `stb-suite` main prompt, type `4.15.1` and answer the guided
prompts (structure file, sublattice, species B, calc template,
pseudopotential source, output directory). Both paths call the exact
same underlying tool; the walkthrough script proves this directly by
diffing `cluster_n0`'s geometry from each path.

## 4. Stage 2: mass-action analysis (`stb-gqcaAnalysis`, code `4.15.2`)

### 4.1 What it does

Reads the 3 relaxed `cluster_n{0,1,2}/` folders, normalizes each to a
per-sublattice-site energy (the 3 folders generally have different cell
sizes — see `read_cluster_energy`'s own docstring), solves the
mass-action equilibrium over a composition/temperature grid, and writes
a long-format CSV (`x, T_K, eta, x_0, x_1, x_2, H_mix_eV, S_mix_eV_per_K,
G_mix_eV, SRO`) plus a `.dat`/`.gplot` pair (`H_mix`/`G_mix` vs. `x` at
one reference temperature) and, optionally (`--save-report`), the same
numbered report persisted to `gqca_stage2.txt`.

### 4.2 Composition/temperature grid

`--x-min`/`--x-max`/`--x-points` set the composition grid (default: the
full `[0,1]` range, 21 points); `--temp` fixes a single `T`, or
`--temp-min`/`--temp-max`/`--temp-points` sweeps a range — mutually
exclusive with `--temp`.

### 4.3 Sanity diagnostics and the miscibility-gap search

`[3] SANITY DIAGNOSTICS` always reconfirms the `x=0`/`x=1` pure-end
-member limits against the real data just read (not a standalone unit
test). `[4] MISCIBILITY GAP` runs the common-tangent search from Section
1.4 — read its `[STRUCTURAL LIMITATION]` note carefully before ever
quoting a "no gap" result as evidence of real miscibility.

### 4.4 Running it both ways

```bash
stb-gqcaAnalysis --directory gqca_study --temp 900
```

or type `4.15.2` at the `stb-suite` main prompt.

## 5. Worked example: comparing against a known real-material fact

Real SIESTA output isn't available in this walkthrough (no SIESTA binary
here) — exactly like `4.13-her`/`4.14-oer`'s own examples, Section 2 of
`example_4.15.sh` writes `calc.out` files with a **hand-chosen** FreeEng,
scaled so all 3 folders land on the intended per-site energy (`cluster_
n0`/`cluster_n2` have 4 Ga/Al sites each from the 8-atom input cell;
`cluster_n1` has 2, from `icet`'s smaller primitive-based cell).

### 5.1 Case A — the near-ideal AlGaAs limit

`e_0 = e_2 = e_1 = -50.0 eV/site` exactly → `Delta_1 = 0` for all `j`.
At `x=0.5`, `T=900 K` (a representative AlGaAs MBE/MOCVD growth
temperature):

```
x_j = [0.25, 0.50, 0.25]     H_mix = 0     SRO = 0
```

This **is** the exact random/binomial mixing distribution — and it is
the well-established textbook fact about real Al(x)Ga(1-x)As: an almost
perfectly ideal (Vegard's-law) solid solution, with no measurable
short-range order at any accessible growth temperature. Case A is chosen
specifically to land here, confirming the solver reduces to the correct
closed form when fed energies consistent with a genuinely ideal alloy.

### 5.2 Case B — a hypothetical, NOT-real-AlGaAs contrast

Same GaAs/Al lattice, but `cluster_n1`'s energy is shifted so `Delta_1 =
-0.1 eV/site` (AB pairing favored) — **not** a claim about real AlGaAs,
purely illustrative:

```
x_1: 0.50 -> 0.78      SRO: 0 -> -0.57
```

`x_1` moves above its random value and SRO goes negative (ordering
-favoring), exactly the expected direction for a negative `Delta_1`
(Section 1.3) — proof the solver tracks the sign/magnitude of whatever
cluster energies you actually give it, real or hypothetical. And even at
this artificially strong bias, the miscibility-gap search still reports
"fully miscible (no gap)" — the convexity limitation of Section 1.4 in
action, not a claim that this hypothetical alloy is really miscible.

### 5.3 Temperature sweep: order fading with heat

Sweeping `--temp-min 300 --temp-max 1500` on Case B's energies, SRO at
`x=0.5`:

```
T =  300 K   SRO = -0.9591
T =  600 K   SRO = -0.7474
T =  900 K   SRO = -0.5681
T = 1200 K   SRO = -0.4491
T = 1500 K   SRO = -0.3686
```

`|SRO|` shrinks monotonically as `T` rises — thermal disorder competing
against (here, hypothetical) ordering energy, the same qualitative order
-disorder temperature dependence real alloys show. This trend, unlike a
specific numeric enthalpy, is a general, model-independent consequence
of the mass-action equilibrium (Section 1.1) and is safe to expect
regardless of the actual material.

## 6. Known, deliberate limitations

- **No genuine phase separation** (Section 1.4) — this model captures
  short-range order only. Do not use a "no gap" result as evidence a
  real alloy is fully miscible.
- **Frustrated-sublattice `Delta_1`** (Section 1.5) is a less clean
  estimate than a bipartite sublattice's — qualitatively meaningful, not
  a rigorous bound. Affects any FCC/HCP/triangular sublattice (zincblende
  and rocksalt cation sublattices included).
- **2D inputs get positions-only relaxation** — `MD.VariableCell` is left
  off entirely for a vacuum-padded input (a confirmed SIESTA syntax for
  excluding just the vacuum axis from cell relaxation could not be
  found during this tool's development). A 2D alloy's in-plane lattice
  mismatch between end-members is therefore not captured. See
  `stb-gqca --help` / `force_full_relaxation`'s own docstring.
- **`S_mix`/SRO formulas are `[UNVERIFIED]`** (Section 1.6) — the
  mixing-enthalpy and mass-action-law equations are confirmed; these two
  are provisional.
- Only the pair (N=2) cluster is implemented — no larger-cluster GQCA
  extension.

## 7. Step-by-step: running this workflow on your own alloy

```bash
stb-gqca -f structure.fdf --sublattice <A> --species-b <B> -c calc.fdf -O gqca_study
# run SIESTA (a real CG/variable-cell relaxation) in every
# gqca_study/cluster_n{0,1,2}/ folder, then:
stb-gqcaAnalysis --directory gqca_study --temp <T_K>
# or --temp-min/--temp-max/--temp-points for a sweep.
```

Read every `[NOTE]`/`[WARNING]` in `gqca_stage1.txt` first, especially if
your sublattice is FCC/HCP/triangular (Section 1.5/3.4) — `Delta_1` is
then a less clean estimate, still qualitatively meaningful. `--save
-report` on Stage 2 persists the same numbered report to
`gqca_stage2.txt` alongside the CSV/plot files.

## 8. Files in this folder

| File | Purpose |
|---|---|
| `structure.fdf` | GaAs, zincblende, conventional 8-atom cubic cell (`a = 5.6533 Ang`) — a real semiconductor with a genuinely FCC (frustrated) cation sublattice. |
| `calc.fdf` | Shared template for the 3 cluster relaxations — `stb-gqca` itself rewrites the `MD.*` directives regardless of what this template says (Section 3.1). |
| `example_4.15.sh` | The guided walkthrough (**not** an automated test — see `test/4-workflow/15-gqca/{prep,analysis}/test.sh` for that, including its own bipartite/2D fixtures). Pauses between sections so you can read before moving on; safe to re-run. |
| `output/` | Created by `example_4.15.sh` when you run it (git-ignored, not checked in). See below. |

## 9. Running the script

```bash
./example_4.15.sh
```

| Case | Command(s) | What it shows |
|---|---|---|
| `output/algaas/` (Stage 1) | `stb-gqca --sublattice Ga --species-b Al` | Pair orbit, 3 cluster folders, the FCC-frustration `[NOTE]`s |
| `output/algaas/ideal.csv` (Stage 2, Case A) | `stb-gqcaAnalysis --temp 900` on `Delta_1=0` energies | Exact random-mixing recovery — the real AlGaAs limit |
| `output/algaas/ordering.csv` (Stage 2, Case B) | Same directory, `Delta_1=-0.1 eV/site` energies | Ordering-favoring shift, and the "no gap" structural limitation |
| `output/algaas/ordering_sweep.csv` | `stb-gqcaAnalysis --temp-min 300 --temp-max 1500 --temp-points 5` | SRO fading with temperature |
| *(no folder — a diff only)* | Stage 1 via `printf … \| stb-suite` | Proof the interactive menu (`4.15.1`) agrees with the CLI |

## What's next

- **`4.7-hubbardu`** — another 3-stage workflow built around reading a
  prior stage's results back through a persisted, machine-readable report
  footer, the same convention Stage 2 here uses for `Sublattice`/
  `Species B`/`AB deviation`.
- **`2.6-stb-sqs`** — a different (special-quasirandom-structure)
  approach to modeling a disordered alloy in a *single* finite cell,
  useful context for what GQCA's independent-cluster treatment leaves out.
- **`4.1-strain`** — the shortest two-stage workflow example in this
  suite, if you'd like a simpler starting point before this one's
  frustrated-sublattice/mass-action machinery.
- Every other `4.x` workflow example generates real candidate/perturbed
  geometries and expects a real SIESTA run in between stages, the same
  two-way (direct CLI / interactive menu) split, and the same
  numbered-report convention used here.
