# 4.3 — Workflow: Cohesive Energy (`stb-cohesive` / `stb-cohesiveAnalysis`)

This workflow has 2 stages: **Stage 1** (`stb-cohesive`, code `4.3.1`) takes
an already-relaxed structure and writes one ready-to-run, single-point
SIESTA folder for the full structure, one Gamma-only folder per chemical
species (each atom alone in a large vacuum box), and — by default — a second
set of isolated-atom-type folders where that same atom is surrounded by
"ghost" copies of its real neighbors, correcting for a systematic bias
localized-basis DFT has for this exact quantity. **Stage 2**
(`stb-cohesiveAnalysis`, code `4.3.2`) reads the finished SIESTA runs back
and reports the cohesive energy per atom — uncorrected and (if present)
BSSE-corrected — plus whether that correction has actually converged.
Together they answer one physical question neither `4.1-strain` nor
`4.2-elastic` does: not how the structure *responds* to a deformation, but
**how strongly bound it is in the first place** — the energy per atom needed
to pull the solid apart into free, isolated atoms.

Both stages live in this one folder and this one tutorial — not one folder
per tool — for the same reason `4.1-strain/` and `4.2-elastic/` do: Stage
1's output only exists to feed Stage 2, and Stage 2 only exists to interpret
Stage 1's output.

## 1. Theory

### 1.1 What "cohesive energy" means

For a structure with `N` atoms and total (Kohn-Sham) energy `E_bulk`, and
each chemical species' isolated-atom energy `E_atom(species)`:

```
E_coh = (E_bulk - sum_i E_atom(species_i)) / N
```

(the sum runs over every atom in the cell, so a species with 2 atoms
contributes `2 * E_atom(species)`). **Negative `E_coh` means a bound solid**
— it's the energy *released* by assembling `N` free atoms into the observed
structure, per atom; equivalently, `-E_coh` is the energy you'd have to
*supply* to pull the solid apart back into isolated atoms. A more negative
number means a more strongly bound material. This is exactly the arithmetic
`stb-cohesiveAnalysis`'s own `[3] COHESIVE ENERGY RESULTS` table performs,
with `E_bulk` from `structure/calc.out` and each `E_atom` from `atoms/<species>/calc.out`.

### 1.2 Why LCAO cohesive energies over-bind: Basis Set Superposition Error

SIESTA uses **localized atomic-like orbitals (LCAO/PAO)**, not plane waves —
a finite set of basis functions centered on each atom, not a smooth,
basis-independent grid filling all of space. This creates a subtle but real
bias specific to cohesive-energy-type calculations: when an atom sits alone
in `atoms/<species>/`, it only has *its own* basis functions to vary the
wavefunction with. When that *same* atom sits inside the real solid, it also
has variational freedom from *its neighbors'* basis functions — even though
those neighbors' own electrons aren't "given" to it, their orbital tails
still overlap and let the atom's own wavefunction relax into a
lower-energy shape than its own basis alone would allow. The isolated-atom
reference energy computed with only its own basis is therefore **too high**
(under-relaxed) — and since `E_coh` subtracts that reference, a
too-high isolated-atom energy makes `E_coh` **too negative**: the structure
looks more strongly bound than it really is. This is a basis-set artifact,
not real physics — it doesn't happen in plane-wave DFT, where every
calculation already uses the same basis (a uniform grid) regardless of
what's actually inside the box. Section 4 shows exactly how large this bias
is for a real material (spoiler: large — about 40%).

### 1.3 The counterpoise (ghost-atom) correction

The standard fix, adapted here from quantum chemistry's own
counterpoise-correction literature, is to give the isolated atom the *same*
basis-function environment it has inside the solid, **without adding any
real physics**: recompute its energy while surrounding it with "ghost"
copies of its actual nearest neighbors — same basis orbitals as the real
element, but (SIESTA's own ghost-atom convention) **zero nuclear charge and
zero valence electrons**. A ghost atom contributes basis functions for the
real atom to relax into, and nothing else — no extra electron, no extra
nuclear attraction, no new physics. `stb-cohesive` builds this automatically
from the real structure's own local geometry (`--bsse-correction`, ON by
default): for each species, it finds that species' real nearest neighbors
in `structure.fdf` within a chosen radius (section 1.4) and writes a small
cluster — the one real atom at the center, ghost copies of its real
neighbors placed at their real relative positions — into
`atoms_bsse/<species>/`. `stb-cohesiveAnalysis` then reports the cohesive
energy computed both ways side by side.

### 1.4 The cutoff radius (`--bsse-cutoff`): what it controls, how to read it

`--bsse-cutoff` (default `4.0` Ang) is the radius, in Angstrom, around each
species' real anchor atom within which *real* neighbors from the structure
become *ghost* atoms in its BSSE cluster. How many coordination shells this
actually reaches depends entirely on the material's bond length — for
short-bonded covalent solids (graphene, diamond, ~1.4-1.5 Ang bonds) even a
modest cutoff reaches several shells; for longer-bonded ionic/metallic
solids it may only reach the first. `stb-cohesive`'s own `[5] BSSE
(COUNTERPOISE) CORRECTION` table prints the actual **"Ghost nb."** count it
found for each species/site — this is the direct, checkable signal for
whether a chosen cutoff is doing anything at all (a count of `0` prints a
`[WARNING]`: the "corrected" reference would be identical to the
uncorrected one).

### 1.5 Convergence: why one cutoff is never enough on its own

Choosing a cutoff and trusting it blindly is not enough — the correction
must be shown to have **stabilized** as the cutoff grows, the same
methodological discipline `4.1-strain`/`4.2-elastic` apply to their own
strain range/mesh choices. `stb-cohesive` offers two ways to check this,
both gated behind `--bsse-convergence-check`:

- **A single before/after point** — `--bsse-convergence-increment` given
  *one* value (the default, `2.0`) generates exactly one extra reference at
  `--bsse-cutoff + increment`, written to `atoms_bsse_check/`.
  `stb-cohesiveAnalysis` reports the shift between the two as
  `BSSE cutoff convergence shift` in `[3]` — small means converged at this
  cutoff already; large means keep increasing `--bsse-cutoff`.
- **A full multi-cutoff scan** — give `--bsse-convergence-increment` *two or
  more* values (e.g. `--bsse-convergence-increment 2 4 6`) and
  `stb-cohesive` generates one extra reference per value, each in its own
  `atoms_bsse_check_<cutoff>/`. `stb-cohesiveAnalysis` auto-detects this and
  adds a `[3b] BSSE CUTOFF CONVERGENCE SCAN` section — a table of cohesive
  energy per atom at every scanned cutoff plus the change on the *last*
  step — and a third panel on `cohesive_correction.png` plotting the whole
  trend. This is the only way to actually *see* whether the correction has
  plateaued, rather than infer it from a single before/after gap; section 4
  tells the real story of why this scan feature exists at all.

### 1.6 The vacuum box, periodic images, and a real bug found while building this

SIESTA is **always periodic** — there is no true "isolated atom" mode.
`atoms/<species>/` and every BSSE cluster are really a sparse periodic
lattice of copies of themselves, spaced apart by a large empty
(`--vacuum`-sized) box, and simply relies on that spacing being large
enough that neighboring copies don't meaningfully interact through their own
basis orbitals' tails. If the box is too small for how far those orbitals
actually reach, you get a real, physical-looking but entirely spurious
energy shift from an atom "seeing" its own periodic image.

A real bug in this exact spot was found and fixed while building this
workflow: the ghost-cluster box used to be a **fixed** `--vacuum`-sized cube
*regardless* of `--bsse-cutoff` — so a larger convergence-check cutoff
(section 1.5) placed its outermost ghost atom farther from the box center
while the box itself stayed the same size, silently eating into the very
buffer meant to prevent periodic self-interaction. This is exactly backward:
the moment you deliberately widen the cutoff to *check* convergence is the
moment the box most needs to grow too. The fix: the box side is now
`2 * cutoff + --vacuum` for every BSSE-type reference — a **plain isolated
atom** (`atoms/<species>/`, cutoff `0`) is unaffected (box `= --vacuum`,
exactly as before), while a BSSE cluster or convergence-check cluster always
keeps a guaranteed `--vacuum / 2` buffer beyond its own outermost ghost
atom, no matter how large the cutoff gets. Section 4 walks through the real
numbers this bug produced (or, as it turned out, mostly didn't — read on).

### 1.7 Multi-site weighting

A species doesn't always sit in one local environment. When a species
occupies 2 or more *symmetrically distinct* crystallographic sites — real,
different local coordination, not just different unit-cell copies of the
same site — a single flat BSSE reference for that species would silently
average over 2 physically different situations. `--bsse-multi-site`
(default ON) detects this via space-group symmetry
(`core/symmetry.py::find_inequivalent_sites`) and builds one ghost cluster
**per distinct site**, in nested `atoms_bsse/<species>/site_<wyckoff>_x<mult>/`
folders. `stb-cohesiveAnalysis` then reads all of them back and computes a
single reference energy **weighted by each site's own multiplicity**
(`read_bsse_energy`'s weighted average) — so a species with, say, 6 atoms in
one environment and 2 in another gets a reference that's actually
representative of its real population, not whichever site happened to be
generated first. Section 2.5 demonstrates this on a real material where it
genuinely matters: AB-stacked bilayer graphene, where carbon has 2 real,
distinct local environments in the *same* structure.

### 1.8 Spin and k-grid conventions for isolated-atom-type calculations

Every isolated-atom-type calculation (plain atom, BSSE cluster, convergence
-check cluster) is **always spin-polarized and Gamma-only**, regardless of
the *full* structure's own `--spin` setting or k-point density — correct
physics for the overwhelming majority of cases (an isolated atom's ground
state is essentially always non-trivially spin-polarized, e.g. a lone
carbon atom has 2 unpaired 2p electrons; a large vacuum box has no
dispersion for a denser k-grid to sample). `stb-cohesive` prints an explicit
`[NOTE]` about this in `[3] FULL STRUCTURE SETUP` whenever the full
structure itself is run non-polarized: if the *bulk* material is expected to
order magnetically, comparing its non-polarized bulk energy against
polarized atomic references is an *inconsistent*, not just incomplete,
comparison — pass `--spin` on the full structure too if in doubt.

## 2. Stage 1: generating the calculations (`stb-cohesive`, code `4.3.1`)

### 2.1 What it does

Reads one structure, then writes (all under one wrapping `-O/--output-dir`,
default `cohesive_runs/`, matching `strain_runs/`/`elastic_runs/`):

- `structure/` — the full structure, one single-point SIESTA folder.
- `atoms/<species>/` — one Gamma-only isolated-atom folder per species that
  actually has atoms in the cell (a species declared in
  `ChemicalSpeciesLabel` but never placed contributes exactly 0 to the
  cohesive-energy sum regardless of its own energy, so it's correctly
  skipped — section 2.6).
- `atoms_bsse/<species>/[site_<wyckoff>_x<mult>/]` — the BSSE-corrected
  reference(s), unless `--no-bsse-correction`.
- `atoms_bsse_check/` (one extra cutoff) or `atoms_bsse_check_<cutoff>/`
  (2+ cutoffs, a scan) — only with `--bsse-convergence-check`.

### 2.2 No `calc.fdf` prerequisite — and why

Unlike `stb-strain`/`stb-elasticInputs`, `stb-cohesive` takes **no**
`-c/--calc-fdf` input at all — it always builds its own single-point SCF
`calc.fdf` internally (DZP basis, GGA-PBE, `MeshCutoff 320 Ry`,
`SCF.DM.Tolerance 1e-5 eV`), identically across the full structure and every
isolated-atom-type calculation. This is deliberate: a cohesive-energy
comparison is only meaningful if every energy in the sum comes from
*exactly* the same basis/functional/mesh — safer to guarantee that by
generating every `calc.fdf` from one template than to trust a hand-edited
file to stay in sync everywhere. `-k/--k-density` (default `0.2` 1/Ang) and
`--spin`/`--dispersion` are the only physics knobs Stage 1 exposes for the
full structure; everything else about the SCF setup is fixed on purpose.

### 2.3 The vacuum-buffer guarantee, live

```
$ stb-cohesive -s structure.fdf -p dojo --bsse-cutoff 4 \
    --bsse-convergence-check --bsse-convergence-increment 2 --no-intro

[5] BSSE (COUNTERPOISE) CORRECTION
------------------------------------------------------------
Species | Site | Ghost nb. | Cutoff  | Used | Folder
------------------------------------------------------------
Ga      | --   | 4         | 4.0 Ang | yes  | atoms_bsse/Ga/
As      | --   | 4         | 4.0 Ang | yes  | atoms_bsse/As/
...
Also generated 'atoms_bsse_check/' at cutoff 6.0 Ang (--bsse-cutoff +
--bsse-convergence-increment) -- stb-cohesiveAnalysis will report the shift
between the two as a convergence diagnostic. Its box is 32.0 Ang wide (vs.
28.0 Ang for the --bsse-cutoff reference) -- both keep the same 10.0 Ang
buffer beyond their respective outermost ghost atom, so this shift reflects
the correction's real cutoff dependence, not a box-size artifact.
```

`atoms_bsse/Ga/structure.fdf`'s own `%block LatticeVectors` is a `28.0 Ang`
cube (`2*4.0 + 20.0`); `atoms_bsse_check/Ga/structure.fdf`'s is `32.0 Ang`
(`2*6.0 + 20.0`) — both leave exactly `10.0 Ang` (`--vacuum/2`) of empty
space beyond their own outermost ghost atom. Case 2 of `example_4.3.sh`
reproduces this exact check on the real GaAs structure in this folder.

### 2.4 Single check vs. multi-cutoff scan, live

```
$ stb-cohesive -s structure.fdf -p dojo --bsse-convergence-check \
    --bsse-convergence-increment 2 4 6 --no-intro

BSSE convergence check: ON, scanning 3 cutoffs (6.0, 8.0, 10.0 Ang)
...
Also generated 3 BSSE convergence-check references -- a full cutoff scan
from --bsse-convergence-increment's multiple values -- so
stb-cohesiveAnalysis can report and plot the cohesive energy vs. cutoff
trend instead of a single before/after point. Every box keeps the same
10.0 Ang buffer beyond its own outermost ghost atom (box side = 2*cutoff +
--vacuum), regardless of cutoff, so the trend reflects the correction's
real cutoff dependence, not a box-size artifact.
Cutoff   | Box side | Folder
--------------------------------------------
6.0 Ang  | 32.0 Ang | atoms_bsse_check_6.0/
8.0 Ang  | 36.0 Ang | atoms_bsse_check_8.0/
10.0 Ang | 40.0 Ang | atoms_bsse_check_10.0/
```

With **exactly one** `--bsse-convergence-increment` value (the default),
naming/behavior is unchanged from before this scan feature existed —
`atoms_bsse_check/`, no suffix. Case 3 of `example_4.3.sh` runs this live.

### 2.5 Multi-site detection and per-site folders, live

Run against `structure_multisite.fdf` in this folder (AB-stacked bilayer
graphene, section 1.7):

```
[5] BSSE (COUNTERPOISE) CORRECTION
------------------------------------------------------------
[INFO] 2 symmetrically distinct site(s) detected for C: d (x2), c (x2)
Species | Site   | Ghost nb. | Cutoff  | Used | Folder
-----------------------------------------------------------------------
C       | d (x2) | 24        | 4.0 Ang | yes  | atoms_bsse/C/site_d_x2/
C       | c (x2) | 22        | 4.0 Ang | yes  | atoms_bsse/C/site_c_x2/
...
Space group(s) detected (per species, same structure -- should agree): C=P-3m1 (No. 164)
```

The 2 sites are the well-known "eclipsed"/dimer vs. "non-eclipsed" carbon
environments of Bernal-stacked bilayer graphene (section 1.7's own
structural explanation) — 2 real coordination environments for the exact
same element, correctly kept separate rather than averaged into one
reference. Case 4 of `example_4.3.sh` runs this live and shows the nested
folder layout on disk.

### 2.6 Phantom species (declared but never placed) are correctly skipped

A `structure.fdf` can declare a species in `%block ChemicalSpeciesLabel`
that never actually appears in `%block AtomicCoordinatesAndAtomicSpecies`
(e.g. hand-edited from a template, or a placeholder for a future
substitution). Since such a species contributes exactly `0 * E_atom` to the
cohesive-energy sum no matter what its own energy would be, `stb-cohesive`
correctly does **not** generate an `atoms/<species>/` (or BSSE) folder for
it, and `[1] INPUT STRUCTURE` prints an explicit
`[INFO] Declared but never placed (no isolated-atom calculation needed): <species>`
instead of silently requiring (and blocking on) a calculation nobody needs.
Case 5 of `example_4.3.sh` demonstrates this directly.

### 2.7 Output layout, and running it both ways

| Path | What it is |
|---|---|
| `cohesive_runs/structure/` | Full structure, single-point SCF |
| `cohesive_runs/atoms/<species>/` | Isolated atom, Gamma-only, spin-polarized |
| `cohesive_runs/atoms_bsse/<species>/[site_.../]` | BSSE-corrected reference(s) |
| `cohesive_runs/atoms_bsse_check[_<cutoff>]/...` | Convergence check/scan reference(s) |

Both CLI and interactive menu (`stb-suite` → `4.3.1`) produce identical
output for identical answers — Case 6 of `example_4.3.sh` proves this with
a real `diff -rq`.

## 3. Stage 2: analyzing the results (`stb-cohesiveAnalysis`, code `4.3.2`)

### 3.1 Report structure

```
[0] RUN METADATA
[1] INPUT STRUCTURE
[2] ENERGY EXTRACTION
[3] COHESIVE ENERGY RESULTS
[3b] BSSE CUTOFF CONVERGENCE SCAN     <- only when a multi-cutoff scan was found
[4] CORRECTION PLOT
[5] SUMMARY & FILES
```

### 3.2 `-d/--dir` auto-defaults to `cohesive_runs`

Matching Stage 1's own `-O/--output-dir` default, so
`stb-cohesiveAnalysis -o calc.out` works with **no `--dir` at all** in the
common case, run from the same parent directory Stage 1 was — no manual
`cd` required either way. An explicit
`[FAIL] Directory '<dir>' not found.` fires immediately if it's missing,
rather than failing indirectly later.

### 3.3 Numerical-quality diagnostics

`[2] ENERGY EXTRACTION` checks the **full structure's own** `calc.out` (not
the isolated-atom-type ones — an isolated atom has zero net force by
symmetry, and a BSSE cluster's real atom only feels basis-incompleteness
-driven Pulay forces from its ghosts, not a physical force worth "relaxing"
away): SCF convergence, and residual force against `--force-tolerance`
(default `0.05` eV/Ang) — a `[WARNING]` here means the cohesive energy may
reflect a strained/off-equilibrium geometry, not the true minimum. Advisory
only, never blocks the result.

### 3.4 Reading the `[2]`/`[3]` tables

`[2]`'s per-species table has a **Delta** column for each correction,
relative to the uncorrected reference — since BSSE removes an artificial
over-stabilization of the isolated atom (section 1.2), a **positive** Delta
here is the expected, physically correct direction (the isolated atom's
energy going *up*, closer to what it should be). `[3]`'s method table
mirrors this with a **Delta vs Uncorrected** column on the final cohesive
energy per atom itself — the single number that answers "how much did BSSE
actually change my answer".

### 3.5 `[3b] BSSE CUTOFF CONVERGENCE SCAN`

Only printed when a multi-cutoff scan (section 2.4) is found. One row per
scanned cutoff (plus the base `--bsse-cutoff` reference, labeled `(base)`,
recovered from a small `bsse_cutoff.txt` sidecar Stage 1 writes alongside
every reference — degrades gracefully to "value not recoverable" for an
older run generated before this file existed, still showing the scan points
themselves). The closing **"Last step change"** line is the number to
actually look at: small means the last cutoff increase barely moved the
answer (converged); still large means keep increasing
`--bsse-cutoff`/`--bsse-convergence-increment` further. Section 4 shows this
exact section with real, physically-verified numbers.

### 3.6 The matplotlib plot

`cohesive_correction.png` is **always saved** (no flag needed) whenever a
complete BSSE-corrected reference is available — 2 panels normally (per
-species correction shift relative to uncorrected; net effect on cohesive
energy per atom), a 3rd panel (cohesive energy per atom vs. cutoff) added
automatically whenever `[3b]`'s scan data exists. `--view` additionally
pops the same figure up interactively (a no-op under a non-interactive
backend like the `MPLBACKEND=Agg` this example's own script uses).

### 3.7 Running it both ways

Both CLI and `stb-suite` → `4.3.2` (which now also prompts for
`--save-report`/`--view`, not just the required filename/directory) produce
the identical report for identical answers.

## 4. Worked example: a real GaAs calculation, end-to-end

This is not illustrative data — it's the exact real investigation that
motivated sections 1.5/1.6/2.4 and the multi-cutoff scan feature itself,
run on the real `structure.fdf` in this folder (GaAs, zincblende, DZP/GGA
-PBE, `dojo` pseudopotentials, `-k 0.2`, `--bsse-cutoff 4.0`,
`--bsse-convergence-check`).

**First pass** (before the vacuum-buffer fix, section 1.6 — box was a fixed
20 Ang cube regardless of cutoff):

| Method | Per atom (eV/atom) | Delta vs Uncorrected |
|---|---|---|
| Uncorrected | -4.6876 | (reference) |
| BSSE-corrected (4.0 Ang) | -3.6479 | +1.0397 |
| BSSE check (6.0 Ang) | -3.1724 | +1.5152 |

Reading this the way sections 1.2/1.5 set up:

- **The uncorrected result is badly over-bound** — a commonly cited
  experimental/reference cohesive energy for GaAs is around `~3.3` eV/atom;
  `-4.6876` is roughly **42% too strong**, the textbook BSSE signature
  (section 1.2) at essentially textbook scale.
- **BSSE-corrected (`-3.6479`) lands much closer** to that reference value
  (~11% high, a very ordinary residual GGA-PBE error) — strong evidence the
  correction is doing exactly what section 1.3 says it should.
- **BSSE check (`-3.1724`) is a large, suspicious jump** — `+0.4755 eV/atom`
  between a 4.0 and 6.0 Ang cutoff is about 13% of the entire cohesive
  energy, far too large to treat the 4.0 Ang result as converged.

That large jump is what triggered the actual investigation: was it real
physics (the cutoff genuinely too small), or the vacuum-box artifact section
1.6 warns about? The evidence pointed at the box: with the *old* fixed 20
Ang cube, the 6.0 Ang-cutoff cluster's outermost ghost atom sat only `4.0`
Ang from the box edge (`half_box - cutoff = 10 - 6`) — and the real PAO
orbital cutoff radii read straight out of `calc.out` (`grep "rc ="`) showed
Ga's own basis reaching as far as **7.12 Ang**, comfortably farther than
that 4.0 Ang buffer. A real risk of periodic self-interaction, not a
theoretical one.

**Second pass**, after fixing the vacuum-buffer bug (section 1.6: box side
`= 2*cutoff + vacuum`, guaranteeing a constant `10.0 Ang` buffer at *any*
cutoff) and re-running for real:

| Method | Per atom (eV/atom), fixed box | Change vs. first pass |
|---|---|---|
| BSSE-corrected (4.0 Ang) | -3.6502 | +0.0023 eV/atom |
| BSSE check (6.0 Ang) | -3.1741 | +0.0017 eV/atom |

**The numbers barely moved at all** — about 1-3 meV, an order of magnitude
smaller than typical real-space-grid ("eggbox") noise, let alone a
correction of order 1 eV. If the old fixed box really had been contaminating
the result, tripling its buffer (4.0 → 10.0 Ang) should have changed the
answer visibly. It didn't.

**Conclusion**: the large 4.0→6.0 Ang shift is **real physics, not a box
-size artifact** — 4.0 Ang (4 ghost neighbors, the first coordination shell
only) simply isn't a large enough cluster for GaAs's BSSE correction to have
converged. The vacuum-buffer fix was still the right thing to do (it removes
a real risk that, for a different material or a larger requested cutoff,
absolutely could bite — sections 1.6/2.3), but it wasn't *this* material's
actual problem. This is exactly what directly motivated building the
multi-cutoff scan (section 1.5/2.4): a single before/after pair can't tell
you *whether* a shift is converging, only that it changed. **An honest note
for anyone reproducing this**: a real multi-point scan
(`--bsse-convergence-increment 4 8 12`, say) has not yet been run on this
real structure — a genuine next step, not a result already in hand.

One more small, secondary finding along the way, purely for illustration of
how much a numbered report surfaces without being asked: `stb-cohesive`
detected space group `R3m` (trigonal, No. 160) for this structure, not the
ideal cubic zincblende `F-43m` (No. 216) — traced to the full structure's
own SIESTA stress tensor showing a small residual shear (`-1.09` kbar
alongside a `-19.03` kbar hydrostatic component), i.e. this specific relaxed
geometry carries a small real trigonal-type distortion, not a symmetry
-detection bug. Zero practical effect here (each species has only 1 atom in
the cell, so site multiplicity is 1 regardless of point group), but exactly
the kind of thing worth noticing before trusting a multi-site weighting
result (section 1.7) on a less trivial structure.

## 5. Known, deliberate limitations

- **No `calc.fdf` customization** (section 2.2) — basis/functional/mesh are
  fixed by `stb-cohesive`'s own template; only `-k`, `--spin`, `--dispersion`
  are exposed, on purpose.
- **`--bsse-cutoff`/`--vacuum` are user-chosen, never auto-optimized** —
  `stb-cohesive` reports the ghost-neighbor count and (with
  `--bsse-convergence-check`) the convergence trend, but never picks a
  "safe" cutoff for you; section 4's own real result shows why that would
  be genuinely hard to automate reliably.
- **The convergence scan is diagnostic only** — it tells you *whether* the
  correction has stabilized, not what cutoff to use if it hasn't; you still
  have to re-run with a larger `--bsse-cutoff`/wider
  `--bsse-convergence-increment` values yourself.
- **Isolated-atom-type spin/k-grid are fixed** (section 1.8) — always
  spin-polarized, always Gamma-only, not independently configurable from
  the full structure's own settings.
- **Dimensionality/vacuum-gap detection** shares the exact same heuristic
  (and the same caveats) as `stb-kgrid`/`stb-strain`/`stb-elasticInputs`.

## 6. Step-by-step: running this workflow on your own structure

0. **Prerequisite**: an already-relaxed structure (`structure.fdf`,
   fractional coordinates). Unlike `4.1`/`4.2`, there is **no** separate
   `calc.fdf` to prepare (section 2.2) — `stb-cohesive` builds its own.

1. **Decide on `--bsse-correction`** (default ON, section 1.2/1.3): keep it
   on unless you specifically want the (systematically over-bound)
   uncorrected number for some other reason.

2. **Decide `--bsse-cutoff`** (section 1.4): start at the default `4.0` Ang,
   check the printed "Ghost nb." count per species, and plan to verify
   convergence (step 5) rather than trusting one value blindly.

3. **Run Stage 1** (`stb-cohesive`, code `4.3.1`):
   - CLI: `stb-cohesive -s structure.fdf -p <bank> --bsse-convergence-check
     --bsse-convergence-increment <inc1> [inc2 ...]` (one increment for a
     single check, several for a full scan — section 2.4).
   - Interactive: `stb-suite` → `4.3.1` — same questions instead of flags
     (the convergence-increment prompt also accepts several space-separated
     values now).
   - Writes everything under `cohesive_runs/` (section 2.1/2.7).

4. **Run SIESTA yourself in every generated folder.** Stage 1 never runs
   SIESTA — `cd` into `structure/`, each `atoms/<species>/`, each
   `atoms_bsse/...` and `atoms_bsse_check.../` folder and run it like any
   other SIESTA calculation.

5. **Run Stage 2** (`stb-cohesiveAnalysis`, code `4.3.2`) once every run
   finishes:
   - CLI: `stb-cohesiveAnalysis -o calc.out` (no `--dir` needed in the
     common case, section 3.2).
   - Interactive: `stb-suite` → `4.3.2`.

6. **Read `[2]`'s SCF/force diagnostics** (section 3.3) before trusting the
   numbers, then **`[3]`/`[3b]`** (sections 3.4/3.5) for the actual
   cohesive-energy result and whether the BSSE correction has converged —
   if `[3b]`'s last step is still large, go back to step 2 with a larger
   cutoff/wider scan.

7. **(Optional) persist or view**: `--save-report` writes
   `cohesive_results.dat`; `--view` opens `cohesive_correction.png`
   interactively (it's always saved to disk regardless, section 3.6).

## 7. Files in this folder

| File | What it is |
|---|---|
| `structure.fdf` | Real, relaxed GaAs, 2 atoms/cell — the material section 4's real numbers are for |
| `structure_multisite.fdf` | Purpose-built AB-stacked bilayer graphene, 4 C atoms/cell, 2 distinct sites — multi-site demo only (section 2.5) |
| `example_4.3.sh` | This walkthrough's runnable script (both stages, 8 cases) |

## 8. Running the script

```bash
bash example_4.3.sh
```

## 9. What's next

- **`4.1-strain`** / **`4.2-elastic`** / **`4.4-phonons`** — a different
  mechanical or vibrational question (how the structure *responds* to a
  deformation, or how it vibrates) about the same kind of already-relaxed
  structure this workflow starts from.
- There is currently **no ML-preview twin tool** for this workflow (unlike
  `4.2-elastic`'s `stb-mlelastic`) — a real, documented gap in the suite
  today, not an oversight in this tutorial.
