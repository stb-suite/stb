# 4.8 — Workflow: Adsorption (`stb-adsorb` / `stb-adsorbBsse` / `stb-adsorbAnalysis`)

Computing a real DFT adsorption energy means comparing three numbers —
the combined slab+adsorbate system, the bare slab, and the isolated
adsorbate — each from its own SIESTA calculation:

```
E_ads = E_site - E_clean_slab - E_adsorbate
```

Getting an honest number out of that comparison depends on a handful of
details that are easy to get subtly wrong by hand: whether the cell is
allowed to relax, whether the slab's own periodic images interact with the
adsorbate, and whether the basis-set mismatch between the three separate
calculations is corrected for. This item has **three** stages, not two:

- **Stage 1 (`stb-adsorb`, code `4.8.1`)** generates every folder the
  uncorrected `E_ads` needs (`clean_slab/`, `adsorbate*/`, one folder per
  candidate `sites/site_*/`), with the fixed-cell/self-interaction details
  below already handled correctly.
- **Stage 2 (`stb-adsorbBsse`, code `4.8.2`)** generates the
  Boys-Bernardi counterpoise (BSSE) correction — but only once Stage 1's
  own sites have actually finished relaxing in SIESTA, since the
  correction is only meaningful at the real, relaxed geometry (Section 2).
- **Stage 3 (`stb-adsorbAnalysis`, code `4.8.3`)** reads everything Stage
  1/2 produced and reports `E_ads`, with and without the BSSE correction,
  ranked most stable first.

`example_4.8.sh` is a guided, runnable walkthrough of all three stages —
**including a worked numeric example (Section 4) where the BSSE
correction changes which site you'd conclude is the most stable one.**
That is the single most important thing to take from this page: an
uncorrected adsorption-energy ranking is not just "a little too negative,"
it can point at the wrong site outright.

## 1. Why the cell stays fixed

Every folder `stb-adsorb`/`stb-adsorbBsse` write (`clean_slab/`,
`adsorbate*/`, every `sites/site_*/`, and the BSSE ghost-fragment folders)
gets a `config_extra.fdf` sidecar forcing `MD.VariableCell false`, `%include`-d
ahead of your own `--calc` template (same sidecar-override convention
`stb-hubbardu`/`stb-elasticInputs`/`stb-phononsCreate` already use, since
SIESTA's fdf reader is first-occurrence-wins for a duplicate directive —
the override has to physically precede anything your own template
already sets).

Two independent reasons this has to be fixed, not left to whatever your
`--calc` template happens to say:

- **The lattice constant is a bulk property.** The slab's in-plane lattice
  vectors come from the underlying bulk crystal's own equilibrium
  geometry — a single adsorbed atom or molecule has no physical business
  changing it. If the cell were allowed to relax, the "adsorption energy"
  you'd compute would be contaminated by a spurious change in the
  substrate's own lattice constant, not just the adsorbate's binding.
- **The vacuum gap is a numerical parameter, not a physical one.** The
  gap along `c` exists purely so the slab doesn't interact with its own
  image above/below it (Section 3.2 below is the in-plane analogue of
  this same idea). Letting `MD.VariableCell` touch that axis would change
  a number that was never meant to have a physical value in the first
  place.

Atomic positions are a different story: whether the slab, adsorbate, and
site should relax their **positions** is entirely up to your own `--calc`
template's `MD.TypeOfRun`/`MD.NumCGsteps` (a real adsorption study
normally does want that) — `stb-adsorb` never touches those. The one
exception is the BSSE ghost-fragment folders, next section.

## 2. BSSE — the counterpoise correction, and why it needs a stage of its own

### 2.1 The problem

SIESTA is a localized-basis (LCAO) code. A fragment computed **alone**
only has its own atoms' basis functions available to describe its
electrons. In the **combined** slab+adsorbate calculation, each fragment
additionally has the *other* fragment's basis functions sitting nearby,
even though those orbitals belong to atoms that, from that fragment's own
point of view, aren't really part of it. This extra, unearned
variational freedom systematically **lowers** the combined system's
energy relative to either fragment's own smaller-basis calculation —
Basis Set Superposition Error (BSSE) — making every LCAO adsorption
energy look more strongly bound than it really is, purely as a basis-set
artifact, not a physical effect. Section 4 below shows exactly how
misleading that artifact can be: not just "a bit too strong," but capable
of flipping which candidate site looks best.

### 2.2 The fix: Boys-Bernardi counterpoise, both directions

The standard correction (Boys & Bernardi, 1970) re-evaluates each
fragment's energy **in the full combined basis**, by adding "ghost"
atoms standing in for the atoms that are physically absent from that
fragment — same real basis functions and pseudopotential, but zero
nuclear/valence charge, so they add basis freedom without adding any real
physics. `stb-adsorbBsse` (Stage 2) writes both directions, at the exact
RELAXED geometry of the site they belong to:

- `bsse_slab/` — the real **substrate**, with the adsorbate replaced by
  ghosts → the substrate's own energy, evaluated in the full dimer basis.
- `bsse_adsorbate/` — the real **adsorbate**, with the substrate replaced
  by ghosts → the adsorbate's own energy, evaluated in the full dimer
  basis.

`stb-adsorbAnalysis` (Stage 3) then applies:

```
E_ads(BSSE-corrected) = E_site - E_bsse_slab - E_bsse_adsorbate
```

Ghost species use SIESTA's own ghost-atom convention: a new species label
(`<element>_ghost`) with a **negative** atomic number in `%block
ChemicalSpeciesLabel`, resolving to the *same* pseudopotential file as the
real element (`core/pseudopotentials.py::copy_pseudo`'s `dest_label`
parameter) but contributing zero charge to the calculation.

Because both fragments here — a finite slab, a finite adsorbate — are
already complete, self-contained systems, this counterpoise is **exact**:
there's no truncation choice to make in either direction. Contrast this
with `stb-cohesive`'s own BSSE correction for a *periodic bulk* crystal,
which has no natural "edge" to stop at and needs an arbitrary
neighbor-shell cutoff (`--bsse-cutoff`) instead.

### 2.3 Why BSSE needs its own stage, AFTER the sites relax

A counterpoise correction only means what it's supposed to mean if every
fragment is evaluated at the **same geometry** — the geometry of the
actual, combined, relaxed site. That geometry doesn't exist yet at Stage
1 (`stb-adsorb`) time: `sites/site_*/` is a real relaxation (whatever
`MD.TypeOfRun` your `--calc` template requests), and its atoms only
settle into their true adsorption geometry once SIESTA has actually
finished running there. An earlier version of this workflow wrote the
BSSE ghost-fragment folders directly from `stb-adsorb`'s own pre-
relaxation initial guess — which, for a real chemisorption bond, can be
substantially wrong: a live test case saw an O-Si bond shrink from a
2.0 Ang initial guess to ~1.62 Ang after real relaxation (plus ~0.46 Ang
of substrate reconstruction), and BSSE grows quickly as atoms get closer
together, so evaluating it at the wrong (too-separated) geometry
under-estimates the true correction.

`stb-adsorbBsse` (Stage 2) fixes this by running strictly AFTER Stage 1's
sites have finished relaxing: for each `sites/site_*/` with a finished
`siesta.XV`, it reads that RELAXED geometry directly (skipping, with a
clear report, any site that hasn't finished yet) and writes
`bsse/site_*/bsse_slab/`+`bsse/site_*/bsse_adsorbate/` at that exact
geometry — with `config_extra.fdf` additionally forcing single-point SCF
on top of the fixed cell (`MD.TypeOfRun CG`, `MD.Steps 0`,
`MD.NumCGsteps 0` — the same single-point-enforcement idea `stb-hubbardu`
uses for its own perturbation probes, for the analogous reason: a
counterpoise correction, like a linear-response slope, is only meaningful
when every point being compared shares the exact same geometry). If
`bsse_slab/`/`bsse_adsorbate/` were instead allowed to relax further,
each ghost-fragment calculation would drift to its *own* locally-optimal
geometry instead of staying at the site's, and the three energies being
subtracted would no longer refer to the same physical configuration.

## 3. Is the slab big enough? Two separate self-interaction checks

A periodic slab calculation has *two* independent directions in which a
badly-sized supercell can make the adsorbate spuriously interact with a
copy of itself, and `stb-adsorb` checks both. `example_4.8.sh` exercises
both live.

### 3.1 Out-of-plane: the isolated-adsorbate reference box

The isolated-adsorbate reference (`adsorbate/`) sits in a cubic
`--vacuum-box` (default 20 Ang). If the molecule's own extent exceeds
half that box, it can reach its own periodic image in *that* separate,
throwaway cell — `stb-adsorb` already warns about this
(`molecule_extent` vs. `--vacuum-box`).

### 3.2 In-plane: the slab's own lateral supercell

This is a different, easy-to-miss failure mode: the **site** folders use
the slab's *own* in-plane periodicity (vectors `a`/`b`), which is fixed by
whatever structure you gave `stb-adsorb` — not something this tool
controls. If that in-plane supercell is small (a bare primitive cell, for
instance), the adsorbate's periodic images in neighboring cells can sit
right next to each other, turning what was meant to be an isolated
adsorption event into an artificial 2D "monolayer" of adsorbates at full
coverage instead.

`stb-adsorb` checks this directly: `min_adsorbate_image_distance`
computes the minimum distance between the adsorbate and the nearest copy
of itself across the slab's `a`/`b` lattice vectors (translation-
invariant — it only depends on the lattice and the adsorbate's own
internal geometry, not which site it ends up on, so it's checked once per
adsorbate, in `[1] REFERENCE FOLDERS`, before any site is even placed).
Below **8 Ang** (a common rule-of-thumb separation for avoiding lateral
adsorbate-adsorbate interaction — dipole coupling, steric contact — in
the literature), it prints a `[WARNING]` recommending a larger in-plane
supercell (`stb-supercell`).

This is exactly the failure mode `structure.fdf`'s bare 2-atom graphene
primitive cell (lattice constant 2.46 Ang) demonstrates on purpose in
`example_4.8.sh`'s `output/stage1_lateral_warning/` section — *far* below
the 8 Ang recommendation, since it's only meant to exercise `stb-adsorb`'s
mechanics quickly, not to be a physically meaningful adsorption
calculation on its own. `output/stage1_fixed/` immediately follows it up
by tiling the same cell 4x4 in-plane with `stb-supercell` and re-running
`stb-adsorb`, confirmed live to make the warning disappear — that bigger
cell is what the rest of the walkthrough (Section 4) actually uses.

## 4. Libraries and external dependencies used

- **`pymatgen`** (`pymatgen.analysis.adsorption.AdsorbateSiteFinder`) —
  finds candidate ontop/bridge/hollow sites and reduces them by symmetry
  (`--symprec`). Note: for `structure.fdf`'s graphene basis, pymatgen's
  finder never returns a `hollow` candidate at any cell size tried —
  ontop and bridge are what the worked example below compares; this is a
  property of the site-finding algorithm on this exact basis, not a
  limitation `stb-adsorb` itself imposes.
- **`sisl`** — reads each site's relaxed `siesta.XV` in Stage 2
  (`stb-adsorbBsse`).
- **`matplotlib`** — the site-layout plot (`adsorption_sites.png`, Stage
  1) and the adsorption-energy ranking plot (`adsorption_ranking.png`,
  Stage 3, `--view-plots` to preview live).
- **`gnuplot`** (optional, always written by Stage 3) — `adsorption_curve.dat`/
  `.gplot`, a portable per-site `E_ads`/`E_ads_BSSE` scatter.
- **MACE** (optional, `pip install stb_suite[ml]`) — only touched by
  `--ml-prerelax` (Stage 1, pre-relax the isolated adsorbate) and
  `--ml-rank` (Stage 1, pre-screen candidate sites before committing to
  real SIESTA folders). Neither is used in this walkthrough.

## 5. Stage 1 mechanics, live (`stb-adsorb`, code `4.8.1`)

```
$ stb-adsorb --adsorbate O -O output/stage1 --no-intro
```

writes `clean_slab/`, `adsorbate/` (forced spin-polarized — many single
atoms/molecules, like O, N, or NO, have a non-zero ground-state spin), and
one `sites/site_1_ontop/` (the default site type/index). By default
(`--force-spin`, ON unless you pass `--no-force-spin`), every
`sites/site_*/` folder is ALSO forced spin-polarized — the physically
relevant quantity, `E_ads`, comes from the **combined** slab+adsorbate
calculation, and a single adsorbate atom bonded to a slab commonly leaves
that combined system with a net magnetic moment too (most simply, an odd
total valence-electron count, which a spin-restricted calculation cannot
represent at all). `clean_slab/` itself is left alone (no adsorbate
present, no clear universal reason to force it). Every folder also gets
`Slab.DipoleCorrection T` forced on `clean_slab/` and every `sites/site_*/`
(never the isolated `adsorbate/` reference, which is a molecule in an
all-around vacuum box, not a slab) — adsorbing on only one face breaks the
slab's inversion/mirror symmetry along the vacuum axis, giving the cell a
net dipole along a PERIODIC direction that would otherwise contaminate the
total energy with a spurious periodic-image field (same reasoning, and the
same fdf tag, `stb-her`/`stb-oer` already force unconditionally for their
own one-sided-adsorbate site folders). `clean_slab/`'s `config_extra.fdf`
is exactly:

```
# Auto-generated by stb-adsorb -- keeps the cell fixed (see the example's
# README.md for why).
MD.VariableCell false
```

while a `sites/site_*/`'s `config_extra.fdf` additionally has:

```
# Auto-generated -- forces spin-polarized SCF (see write_reference_folder's
# force_spin docstring for why).
Spin                polarized
```

Report structure (`[0]`–`[6]`, `--save-report` persists it to
`adsorption_prep_report.txt`):

```
[0] RUN METADATA
[1] REFERENCE FOLDERS
[2] ADSORPTION SITES: FINDING & COUNT
[3] ML PRE-SCREENING
[4] WRITING SITE FOLDERS
[5] SUMMARY & NEXT STEPS
[6] LIBRARY WARNINGS
```

`--site-type all --all-sites` writes one folder per symmetrically-distinct
site of every type instead of just one; `--both-sides` additionally mirrors
each selected site onto the bottom face for a free-standing 2D material
like this graphene fixture (verified: `site_1_ontop_bothsides/structure.fdf`
has `NumberofAtoms 4` — 2 substrate + 1 adsorbate per face); `--height-sweep
<min> <max> <step>` writes an approach curve instead of one fixed height.

**A caution worth knowing before you script this yourself:** calling
`stb-adsorb` twice into the *same* `-O` output directory with two
different single `--site-type` values (e.g. once for `ontop`, once for
`bridge`) overwrites `sites/adsorption_sites.txt` — the machine-readable
site table `stb-adsorbAnalysis` reads for each site's adsorbate/height —
silently dropping the first call's row. Use a single `--site-type all
--all-sites` call instead (as Section 4 below does) if you want more than
one site type in the same study.

Interactive path: `stb-suite`, type `4.8.1` at the main prompt — same
flags, same order, as the CLI above. `example_4.8.sh` proves the two agree.

## 6. Stage 2 mechanics, live (`stb-adsorbBsse`, code `4.8.2`)

```
$ stb-adsorbBsse --dir adsorption_run --save-report
```

Report structure (`[0]`–`[4]`, `--save-report` persists it to
`adsorption_bsse_report.txt`):

```
[0] RUN METADATA
[1] SITE SCAN
[2] WRITING BSSE FOLDERS
[3] SUMMARY & NEXT STEPS
[4] LIBRARY WARNINGS
```

`[1]` prints a table with one row per `sites/site_*/`: `ready` (has a
`siesta.XV`, gets a BSSE folder pair written), `SKIP (no .XV yet -- not
relaxed)` (not fatal — re-run this stage later once it's finished), or a
mismatched/stale-folder warning if a site's atom count no longer matches
its own `structure.fdf`. `[2]` writes `bsse/site_*/bsse_slab/` and
`bsse/site_*/bsse_adsorbate/` only for the `ready` sites, reusing that
site's own already-copied pseudopotentials directly (no `-p` flag on this
stage). Interactive path: `stb-suite`, type `4.8.2`.

## 7. Stage 3 mechanics, live (`stb-adsorbAnalysis`, code `4.8.3`)

```
$ stb-adsorbAnalysis --dir adsorption_run --save-report --apply production.fdf
```

Report structure (`[0]`–`[5]`, `--save-report` persists it to
`adsorption_report.txt`):

```
[0] RUN METADATA
[1] REFERENCE ENERGIES
[2] SITE RESULTS: CONFIGURATION COUNT & TABLE
[3] SUMMARY & PLOT
[4] APPLY
[5] LIBRARY WARNINGS
```

`[2]` prints one row per site (`E_ads`, and `E_ads_BSSE` once Stage 2 has
run there, plus an SCF-convergence/residual-force quality flag — a folder
missing `calc.out` entirely is `SKIP`ped, not fatal to the rest of the
batch) and a `BSSE coverage: complete N, incomplete N, absent N` summary.
`[3]` reports the most stable site **twice** — once by the uncorrected
ranking, once by the BSSE-corrected ranking — plus the BSSE correction's
exact magnitude at the corrected winner. `[4]` (`--apply <path>`) copies
the most stable site's `structure.fdf` to `<path>`, preferring the
BSSE-corrected ranking whenever it's available. Interactive path:
`stb-suite`, type `4.8.3`.

## 8. Worked example: watching the BSSE correction flip the ranking

This is the section `example_4.8.sh`'s `output/workflow/` walks through
line by line, and the reason this whole workflow needs a Stage 2. No real
SIESTA binary is available inside this walkthrough, so — same convention
`4.7-hubbardu`'s own worked example uses for its linear-response fit —
every `calc.out` below is a **hand-chosen, exact** number, picked in
advance to make one very concrete point.

Two candidate sites on the 4x4-supercell-fixed structure (Section 3.2),
`site_1_ontop` and `site_2_bridge`:

| Quantity | Value (eV) |
|---|---|
| `E_clean_slab` | `-200.000000` |
| `E_adsorbate` (isolated O) | `-13.500000` |
| `E_site` (ontop) | `-214.100000` |
| `E_site` (bridge) | `-214.150000` |
| `E_bsse_slab` (ontop) | `-200.050000` |
| `E_bsse_adsorbate` (ontop) | `-13.520000` |
| `E_bsse_slab` (bridge) | `-200.090000` |
| `E_bsse_adsorbate` (bridge) | `-13.560000` |

Uncorrected `E_ads = E_site - E_clean_slab - E_adsorbate`:

```
E_ads(ontop)  = -214.100000 - (-200.000000) - (-13.500000) = -0.600000 eV
E_ads(bridge) = -214.150000 - (-200.000000) - (-13.500000) = -0.650000 eV   <- more negative, "wins"
```

**Bridge looks more stable by 0.050 eV.** Now the BSSE-corrected
`E_ads(BSSE) = E_site - E_bsse_slab - E_bsse_adsorbate`:

```
E_ads_BSSE(ontop)  = -214.100000 - (-200.050000) - (-13.520000) = -0.530000 eV
E_ads_BSSE(bridge) = -214.150000 - (-200.090000) - (-13.560000) = -0.500000 eV   <- now LESS negative than ontop
```

**Ontop is now more stable by 0.030 eV — the ranking flips.** Both
corrections make the binding weaker (less negative), exactly the
direction Section 2.1 predicts (BSSE always over-binds); the bridge site
simply happened to have the larger correction (`+0.150` eV vs. ontop's
`+0.070` eV — plausible physically, since a bridge adsorbate typically
sits closer to more substrate atoms at once, borrowing more basis
functions). `stb-adsorbAnalysis`'s own report says it directly:

```
Most stable site (uncorrected):   site_2_bridge  (E_ads = -0.650000 eV, exothermic (favorable))
Most stable site (BSSE-corrected): site_1_ontop  (E_ads = -0.530000 eV, exothermic (favorable))
BSSE correction at that site: +0.070000 eV (uncorrected LCAO adsorption energies
systematically over-bind -- expect this to make the energy less negative)
```

`--apply best_production.fdf` in this same run copies **`site_1_ontop`**
(the BSSE-corrected winner) forward — proof the tool itself treats the
correction as authoritative once it's available, not just an FYI number
next to the "real" one.

The remaining two candidate sites from Stage 1 (`site_3_bridge`,
`site_4_bridge`) are deliberately left without a `siesta.XV`/`calc.out` in
this walkthrough, at no extra scripting cost demonstrating that neither
Stage 2 nor Stage 3 treats an unfinished site as fatal — both simply skip
it and report exactly which one, and why (`SKIP (no .XV yet -- not
relaxed)` / `SKIP (missing calc.out)`).

## 9. Known, deliberate limitations

- **No SIESTA run is ever performed by this suite.** All three stages
  generate/read input and output files; you run SIESTA yourself in
  `clean_slab/`, every `adsorbate*/`, every `sites/site_*/`, and every
  `bsse/site_*/bsse_slab/`+`bsse_adsorbate/`, in between stages.
- **`--no-force-spin` leaves the combined slab+adsorbate calculation's
  `Spin` setting up to you.** By default (Section 5) both the isolated
  reference AND every `sites/site_*/` folder are forced spin-polarized; if
  you disable that, verify by hand whether `Spin polarized` is still needed
  for an open-shell adsorbate — adsorption often (not always) quenches the
  adsorbate's own spin.
- **The BSSE correction assumes each finished `sites/site_*/`'s
  `siesta.XV` really is the converged, relaxed geometry** — Stage 2 has no
  way to distinguish a genuinely converged run from one that stopped early
  for other reasons; check your own SIESTA convergence before trusting the
  correction.
- **`stb-adsorb`'s own site-finder is pymatgen's `AdsorbateSiteFinder`, not
  a physically-informed search** — as Section 4/8 notes, it can miss a
  real candidate site type entirely (no `hollow` site found for this
  graphene basis at any tested cell size); always sanity-check
  `adsorption_sites.png` against what you'd expect chemically.
- **Calling `stb-adsorb` more than once into the same `-O` directory with
  different single `--site-type` values silently drops the earlier call's
  site-table row** (Section 5) — use one `--site-type all --all-sites`
  call, or separate output directories, instead.

## 10. Step-by-step: running this workflow on your own structure

1. Have a slab/2D structure (`.fdf`, vacuum along `c`) and a working
   `calc.fdf` template.
2. Run Stage 1: `stb-adsorb -s <structure.fdf> -c <calc.fdf> --adsorbate
   <El-or-G2-molecule> --site-type all --all-sites`. Read `[1]`'s
   self-interaction checks (Section 3) before proceeding — tile with
   `stb-supercell` first if either warns.
3. Run SIESTA yourself in `adsorption_run/clean_slab/`, every
   `adsorption_run/adsorbate*/`, and every `adsorption_run/sites/site_*/`.
4. Run Stage 2: `stb-adsorbBsse --dir adsorption_run`. Re-run it later for
   any site reported `SKIP (no .XV yet -- not relaxed)` once that site
   finishes.
5. Run SIESTA yourself in every generated `bsse/site_*/bsse_slab/` and
   `bsse/site_*/bsse_adsorbate/`.
6. Run Stage 3: `stb-adsorbAnalysis --dir adsorption_run --save-report
   --apply production.fdf`. Read `[3]` for **both** rankings before
   trusting a number — if they disagree, trust the BSSE-corrected one
   (Section 8).

## 11. Files in this folder

| File | Purpose |
|---|---|
| `structure.fdf` | Bare 2-atom graphene primitive cell, free-standing (10 Ang vacuum along `c`) — small and fast, deliberately laterally too small to exercise Section 3.2's self-interaction check. |
| `calc.fdf` | Shared `calc.fdf` template (non-polarized, Gamma-only, modest basis/mesh — mechanics only, not physically converged). |
| `example_4.8.sh` | The guided walkthrough (**not** an automated test — see `test/4-workflow/8-adsorption/{prep,bsse,analysis}/test.sh` for that). Pauses between sections so you can read before moving on; safe to re-run. |
| `output/` | Created by `example_4.8.sh` when you run it (git-ignored, not checked in). See below. |

## 12. Running the script

```bash
./example_4.8.sh
```

| Case | Command(s) | What it shows |
|---|---|---|
| `output/stage1/` | `stb-adsorb --adsorbate O -O output/stage1` | Stage 1 mechanics: fixed cell, spin-polarization note |
| `output/stage1_lateral_warning/` | same command | Section 3.2's lateral self-interaction `[WARNING]`, live |
| `output/stage1_fixed/` | `stb-supercell -d 4 4 1` then `stb-adsorb --site-type all --all-sites` | the warning fixed by tiling; the `[2]` candidate-count table |
| `output/stage1_bothsides/` | `stb-adsorb --adsorbate H --both-sides` | 2D free-standing material, both faces |
| `output/workflow/` | Stage 1 → (fabricated `.XV`/`calc.out`) → Stage 2 → (fabricated BSSE `calc.out`) → Stage 3 | **the worked example (Section 8): BSSE flips the ranking** |
| *(no folder — a diff only)* | Stage 1 via `printf … \| stb-suite` | proof the interactive menu (`4.8.1`) agrees with the CLI |

## What's next

- **`stb-supercell`** — tile a small cell (like the graphene primitive
  cell here) into a laterally larger one before running `stb-adsorb`, if
  Section 3.2's check warns you.
- **`stb-cohesive`** — the sibling BSSE correction for a periodic bulk
  crystal (needs a neighbor-shell cutoff, unlike this exact one).
- **`stb-mladsorb`** — a pure-MACE, no-SIESTA extraction of the same
  site-screening idea, for a fast pre-DFT ranking of candidate sites
  before committing to the real workflow above.
