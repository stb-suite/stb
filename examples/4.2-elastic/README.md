# 4.2 — Workflow: Elastic Constants (`stb-elasticInputs` / `stb-elasticAnalysis`)

This workflow has 2 stages: **Stage 1** (`stb-elasticInputs`, code `4.2.1`)
applies canonical Voigt strains to a relaxed structure and writes one
ready-to-run, single-point SIESTA folder per strain value; **Stage 2**
(`stb-elasticAnalysis`, code `4.2.2`) reads the finished SIESTA runs back and
fits the full symmetric stiffness matrix `C_ij`, then checks whether the
material is even mechanically stable. Together they answer a question
`4.1-strain` deliberately does **not**: not "how stiff is this one direction",
but **the complete elastic response of the material** — every normal and
shear stiffness constant, at once, exploiting the crystal's own point-group
symmetry to keep the DFT cost down.

Both stages live in this one folder and this one tutorial — not one folder
per tool — for the same reason `4.1-strain/` does: Stage 1's output only
exists to feed Stage 2, and Stage 2 only exists to interpret Stage 1's
output.

## 1. Theory

### 1.1 From a single modulus to the stiffness tensor `C_ij`

`4.1-strain` fits one number per direction: an "initial slope", `stress =
modulus × strain`, in the spirit of Young's modulus. The real, rigorous
object behind that number is the full **stiffness tensor**, and Hooke's law
in its general form:

```
sigma = C . epsilon
```

`sigma` (stress) and `epsilon` (strain) are symmetric 3×3 tensors; `C` is a
4th-rank tensor with, in general, 21 independent components. Everything in
this workflow (and the rest of the suite) uses **Voigt notation** to collapse
this into linear algebra: the 6 independent stress/strain components indexed
`1..6 = XX/YY/ZZ/YZ/XZ/XY`, so `C` becomes an ordinary symmetric 6×6 matrix
(3×3 for a 2D material's in-plane response — `zz`/`yz`/`xz` are physically
undefined when that axis is vacuum-padded) and `sigma = C . epsilon` becomes
an honest matrix-vector product. This is the exact same Voigt convention
`4.1-strain` already introduced for its own `%block Geometry.Constraints`.

Point-group symmetry is what makes this tractable: for a *general* triclinic
crystal you'd need all 21 constants, each from an independent DFT
calculation. For anything with real symmetry — cubic, hexagonal, tetragonal,
... — most of those 21 numbers are either zero or exact multiples of each
other, and the whole point of `--symmetry-method` (section 1.3) is
exploiting that instead of brute-forcing every component.

### 1.2 Two independent physical routes to the same tensor

`stb-elasticInputs`/`stb-elasticAnalysis` support 2 structurally different
ways to measure `C`, both selected via `--method`:

- **`stress`** (the default) — apply a small strain `epsilon`, read SIESTA's
  own stress tensor back, and fit `sigma = C . epsilon` — a **linear** fit.
  One strain series (several magnitudes of the *same* Voigt direction) gives
  a whole *column* of `C` at once, since `sigma = C . epsilon` for a
  single-component `epsilon` is just `C`'s corresponding column scaled by
  the strain.
- **`energy`** — apply the same kind of strain, read SIESTA's *total energy*
  back instead, and fit the **parabolic** form `E = E0 + (V0/2) * gamma^T C
  gamma` (`gamma` = the Voigt engineering-strain vector, `V0` = the
  unstrained cell's volume/area/length). A **pure** single-component strain's
  energy curvature only ever determines that component's own *diagonal*
  constant (`C_mm`) — off-diagonal constants (`C_12`, ...) need a **combined**
  pattern instead, e.g. `'xx+yy'`, straining 2 Voigt components at once so
  their cross term `2*C_12*gamma_1*gamma_2` shows up in the curvature.
  `--method energy` on the Stage 1 side auto-selects exactly the pure+combined
  pattern set the detected point group actually needs — never a manual list.

These are **genuinely independent observables** — one is a code's stress
tensor, the other its total energy — computed from the *same* DFT method, but
with different numerical sensitivities (section 1.4). Comparing them is a
real physical cross-check, not just re-running the same calculation twice;
Stage 2's mandatory eggbox cross-check (section 3.3) does exactly this
automatically, for free, using data a single `--method stress` run already
collected.

### 1.3 Symmetry reduces how many directions you actually need — `basic` vs `full`

`--symmetry-method` (Stage 1 *and* Stage 2, and they must agree) controls
**how** `--dirs all`'s 6 canonical directions get reduced to the minimal set
actually run:

- **`basic`** — only pools 2 directions when a **single** point-group
  operation maps one exactly onto the other (e.g. a 4-fold rotation mapping
  `xx` onto `yy` in a cubic crystal). This is a *pairwise* test: it cannot
  see reductions that only follow from the *whole* tensor's symmetry at
  once. Hexagonal's `C11 = C22` is exactly such a case (section 4 shows this
  concretely: `basic` needs all 3 in-plane directions on a real hexagonal
  -type 2D material where `full` needs only 1).
- **`full`** — fits the *entire* stiffness matrix in one global least
  -squares pass, constrained a priori to the point group's own
  symmetry-allowed subspace (however many truly independent constants that
  point group has — reported directly, e.g. `"Point group -6m2 has 5
  independent elastic constant(s)"`). This catches every crystal system's
  correct reduction uniformly (cubic, hexagonal, trigonal, ...), at a real
  cost: because every measured direction's data is folded jointly into one
  global fit, `C_12` (say) is no longer independently measurable from 2
  different directions' own data — it's true by construction, not something
  2 separate DFT calculations can be shown to *agree* on. That's exactly why
  Stage 2's tensor-symmetry check (section 3.3) only exists for `basic` —
  see section 4 for the real trade-off in DFT cost vs. this lost
  cross-check.

`basic` is the CLI default (matching `stb-elasticInputs --help`); the
interactive menu (`stb-suite` → `4.2.1`) also defaults there but offers
`full` as option `2`.

### 1.4 Dimensionality changes the physical quantity, not just the unit

Same convention `4.1-strain` already uses, extended to a full tensor instead
of one modulus, auto-detected by Stage 2 from `reference_structure.fdf`'s own
vacuum padding (section 2.5):

- **3D** — `C_ij` in **GPa**, a genuine stress-strain stiffness.
- **2D** (vacuum along `z`) — `C_ij` in **N/m** (a force per unit length, not
  per area) — SIESTA's own volumetric stress (eV/Å³) un-diluted by the cell's
  arbitrary vacuum height `Lz`, the same convention 2D-materials literature
  uses (e.g. graphene's ~340 N/m in-plane stiffness).
- **1D** (vacuum on 2 axes, periodic wire/tube along `z`) — `C_ij` in **nN**,
  a raw axial force, un-diluted by the wire's arbitrary vacuum
  cross-section, no assumed physical cross-section area.

`--method energy`'s `V0` in section 1.2's formula is this same
dimensionality-aware quantity: cell volume (3D), in-plane area (2D), or axial
length (1D) — read from `reference_structure.fdf`, never guessed.

### 1.5 The eggbox effect, and why Stage 2 checks for it automatically

SIESTA evaluates everything on a finite real-space grid (`Mesh.CutOff`). As
atoms/cells move continuously across that fixed grid, quantities computed
from *derivatives* of the energy — forces, and especially the **stress**
tensor — pick up small, spurious, non-physical oscillations tied to the grid
spacing itself: the "eggbox effect". Total **energy** is comparatively smooth
(it's the thing being differentiated, not the derivative), so it's
structurally less susceptible. This is precisely why section 1.2's 2
independent methods (stress vs. energy) make a *physically meaningful* cross
-check, not a redundant one: in the exact continuum limit `sigma(epsilon) =
dE/d(epsilon) / V0`, so for a perfectly converged calculation both methods
must agree exactly — any real disagreement is a genuine hint of eggbox (or
other numerical) contamination, concentrated on the stress side. Stage 2's
`[3] NUMERICAL QUALITY DIAGNOSTICS` (section 3.3) runs exactly this check,
automatically, whenever `--method stress` is used — reusing the *same*
`calc.out` files, no new DFT calculations.

## 2. Stage 1: generating strained structures (`stb-elasticInputs`, code `4.2.1`)

### 2.1 What it does

`stb-elasticInputs` takes a relaxed structure and a reference `calc.fdf`,
applies a canonical Voigt strain (or, for `--method energy`, a pure or
combined pattern) to the lattice vectors for a range of strain values, and
writes one ready-to-run folder per value: the deformed structure, a copy of
`calc.fdf` (with the single-point override added, section 2.2), a new
`config_extra.fdf`, and any linked pseudopotentials. It also writes ONE
`reference_structure.fdf` — the *undeformed* structure — at the top of
`--output-dir` (section 2.4). It does **not** run SIESTA itself — you run
each folder yourself, then hand the results to Stage 2 (section 3).

Unlike `stb-strain`, there is **no `--relax-mode`**: every generated folder
is unconditionally forced to a pure single-point SCF (section 2.2) — the
imposed strain must stay exactly at the sampled geometry for the
stress/energy fit to be physically meaningful, so there's no "which
convention do you want the cell to relax under" choice to make in the first
place.

### 2.2 The single-point-SCF precedence proof, live

`--calc` (`calc.fdf` in this folder) is the SAME file used for
`structure.fdf`'s own real relaxation — DZP basis, a 12×12×1 Monkhorst-Pack
grid (only 1 k-point along the vacuum axis), 320 Ry mesh cutoff, PBE-GGA,
DFT-D3 dispersion, and a REAL, still-active relaxation block:
`MD.TypeOfRun CG`, `MD.VariableCell true`, `MD.Steps 200`. Unlike
`stb-strain` (which only *checks and warns* about a bad relaxation block,
section 2.2 of `4.1-strain/README.md`), `stb-elasticInputs` never trusts
`--calc` to already be single-point — it **forces** one, unconditionally,
regardless of what `--calc` contains:

```
$ stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault --dirs all --symmetry-method full --steps 5 --no-intro
[4] SINGLE-POINT SCF ENFORCEMENT
------------------------------------------------------------
Calc template (current state, before forcing):
  MD.TypeOfRun=CG  Steps: MD.Steps=200  MD.VariableCell=true
Every generated folder is forced to a pure single-point SCF -- NO ionic relaxation, NO cell relaxation -- regardless of the state above, since the imposed strain must stay exactly at the sampled geometry for the stress/energy fit to be meaningful. Forced via '%include config_extra.fdf' PREPENDED at the very top of every generated calc.fdf (before your own directives, including the structure %include) -- SIESTA's fdf reader is first-occurrence-wins for duplicate labels, so this ordering guarantees the forced values win even if your own template already sets any of them (the common case for a real relaxation calc.fdf).

config_extra.fdf (written into every generated folder):
  # Auto-generated by stb-elasticInputs -- forces a pure single-point SCF (no ionic or
  # cell relaxation) at the imposed-strain geometry, regardless of --calc's own settings.
  MD.TypeOfRun       CG
  MD.Steps           0
  MD.NumCGsteps      0
  MD.VariableCell    false
```

The mechanism (`prepend_include`, not `stb-strain`'s own
`insert_include_after_structure`) matters here: `config_extra.fdf` is
`%include`-d at the **very top** of the generated `calc.fdf`, before your own
directives (including the structure `%include`). SIESTA's fdf reader is
**first-occurrence-wins** for duplicate labels — so the forced values win
even though `calc.fdf`'s own `MD.Steps 200`/`MD.VariableCell true` are still
sitting right there, further down, completely untouched:

```
First line of the generated calc.fdf:
%include config_extra.fdf

...calc.fdf's own MD block, still there, further down, silently overridden:
MD.TypeOfRun            CG
MD.VariableCell         true
MD.MaxForceTol          0.01 eV/Ang
MD.Steps                200
```

If the override were instead inserted *after* the structure `%include` (the
convention `stb-strain` itself uses for its own `%block Geometry.Constraints`
— safe there only because that block name never pre-exists in a real
`calc.fdf`), it would be silently ignored whenever `--calc` already sets any
of `MD.TypeOfRun`/`MD.Steps`/`MD.VariableCell` — the common case for a real
relaxation template, exactly like this one.

### 2.3 `--method stress` vs. `energy` at generation time

```
$ stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault --method energy --steps 5 --no-intro
[3] DEFORMATION DIRECTIONS / SYMMETRY REDUCTION
------------------------------------------------------------
[INFO] Point group -6m2 has 5 independent elastic constant(s); 2 strain pattern(s) needed to determine the ones reachable without leaving the periodic plane/axis (of 6 physically valid candidate(s), 21 in total).

  DEFORMATION DIRECTIONS (symmetry-method energy)
------------------------------------------------------------
  Detected symmetry : point group -6m2 -- 12 operation(s) ...
------------------------------------------------------------
  Direction   Status        Reconstructed from
  --------------------------------------------------------
  xx          RUN           --
  yy          SUPPRESSED    symmetry-allowed fit
  xy          RUN           --
  xx+yy       SUPPRESSED    symmetry-allowed fit
  xx+xy       SUPPRESSED    symmetry-allowed fit
  yy+xy       SUPPRESSED    symmetry-allowed fit
  --------------------------------------------------------
  2 of 6 run -- 4 suppressed by symmetry (exactly reconstructible by stb-elasticAnalysis)

Modes: ['xx', 'xy']
```

`--method energy` always auto-selects its own pure+combined pattern set from
the detected point group — there's no `--dirs`/`--symmetry-method` choice to
make for it (`--dirs` stays at its `all` default; `--symmetry-method` is read
by `--method stress` only). Note that for *this specific* point group, the
2 patterns needed both happen to be pure (`xx`, `xy`) — no `'xx+yy'`-style
combined pattern was actually necessary here; the auto-selection is a rank
-selection over every physically valid candidate (`21` combined+pure
patterns total for this material, only `2` needed), not a fixed rule.

### 2.4 `reference_structure.fdf` — what it is, and why Stage 2 needs it

Every `stb-elasticInputs` run writes exactly ONE `reference_structure.fdf` at
the top of `--output-dir` (`elastic_runs/` by default) — a plain copy of the
**undeformed** input structure, never one from inside a `strain_*/` folder.
Stage 2 reads it back for 3 things:

1. **Dimensionality auto-detection** (section 1.4) — from its vacuum
   padding, same `--vacuum-gap` convention as Stage 1's own `[2]` section.
2. **Symmetry/point-group detection** — for `--symmetry-method
   basic`/`full`'s own pooling/reconstruction (section 1.3).
3. **`--method energy`'s reference volume/area/length** (`V0` in section
   1.2's formula) — needed to convert a fitted energy curvature into an
   actual stiffness constant.

It must genuinely be the **undeformed** structure: a single strained
direction genuinely *lowers* the crystal's symmetry (e.g. this example's
`-6m2` structure, strained along `xx` alone, is no longer `-6m2` — it's
exactly tetragonal-like in-plane) — reading symmetry from inside a strained
folder would silently misdetect it, the same reasoning `4.1-strain` gives for
never straining atomic positions directly (only fractional coordinates, which
scale for free with the deformed cell).

### 2.5 Output layout, and running it both ways

`--output-dir <dir>/<direction>/strain_<direction>_<pct>/` — same nested
convention `stb-strain` itself uses, with `reference_structure.fdf` at the
top level (section 2.4) so `cd <output-dir> && stb-elasticAnalysis` (or, from
the parent directory, `stb-elasticAnalysis --dir <output-dir>`) finds
everything with no extra flags:

```bash
# Direct CLI
stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault \
    --dirs all --symmetry-method full --steps 5 --no-intro

# Interactive
stb-suite   # then type 4.2.1, or navigate Workflow -> Elastic Constants -> Stage 1
```

Both call the exact same underlying tool and produce byte-identical output —
`example_4.2.sh`'s own Case 6 proves this directly, driving both paths with
the same inputs and `diff -rq`-ing the result.

## 3. Stage 2: analyzing the results (`stb-elasticAnalysis`, code `4.2.2`)

### 3.1 Report structure

```
[0] RUN METADATA                      -- run dir, method, dimensionality, tolerances
[1] <N>D STIFFNESS MATRIX             -- the fitted C_ij (or C33 alone for 1D)
[1b] <N>D COMPLIANCE MATRIX           -- S = C^-1 (eps = S . sigma), when invertible
[2] (MAIN|RELEVANT) CONSTANTS         -- the labeled diagonal/off-diagonal values
[2b] COUPLING CONSTANTS               -- normal-shear cross terms, 3D only, shown
                                          only when numerically significant
[3] NUMERICAL QUALITY DIAGNOSTICS     -- eggbox cross-check, fit quality, tensor
                                          symmetry (--method stress only; see 3.4)
[4] STABILITY AND PROPERTIES          -- Young's/bulk modulus, Poisson ratio,
                                          anisotropy index, Born stability VERDICT
[5] SUMMARY & FILES                   -- one glanceable results table + file paths
```

`--save-report` persists the exact same report to `elastic_stage2.txt`
(matching `elastic_inputs.py`'s own `elastic_stage1.txt` naming logic).

### 3.2 `-d/--dir` and dimensionality — both auto-detected, no `cd` required

`-d/--dir` (default `elastic_runs`, matching Stage 1's own `--output-dir`
default) accepts either the nested `<direction>/strain_.../` layout Stage 1
writes, or a flat directory of `strain_*/` folders — and, crucially, can be
pointed at from the **parent** directory, no `cd` needed:

```
$ stb-elasticAnalysis --dir elastic_runs --no-intro
[INFO] Detected dimensionality: 2D (e.g., a slab or surface)
[INFO] Scanning 'elastic_runs' for 'strain_*' folders...
```

`--dimensionality auto` (the default) reads `reference_structure.fdf`
(section 2.4) the exact same way Stage 1's own `[2]` section does; pass
`--dimensionality 3d`/`2d`/`1d` to override and skip detection entirely.

### 3.3 `[3] NUMERICAL QUALITY DIAGNOSTICS` — the eggbox cross-check, and 2 more

Whenever `--method stress` is used, 3 checks run automatically — zero new
DFT calculations, all reusing data a single run already collected:

- **Eggbox cross-check** (the theory from section 1.5) — every diagonal
  constant actually measured gets refit a *second*, independent way from the
  *same* `calc.out` files: the usual stress-linear fit, and a parabolic fit
  of the total energy already sitting in the same files. `--eggbox-tolerance`
  (default 5.0%) sets how much disagreement is flagged.
- **Per-direction fit quality** — `R²` and residual stress at zero strain,
  from the *same* stress-strain fit each direction's column already needed.
  `--fit-quality-tolerance` (default `R²=0.995`) flags a noisy or
  non-linear (too-large-`--max`) series, or a non-negligible residual stress
  hinting the reference structure isn't fully relaxed.
- **Tensor symmetry (`C` vs. `C^T`)** — `--symmetry-method basic` **only**
  (section 1.3 explains why): 2 independently-run directions can genuinely
  disagree on the value they each imply for the same off-diagonal constant
  (`C_12` from `xx`'s own column vs. from `yy`'s own column) — a real,
  useful internal-consistency check `full`/`energy` cannot offer, since both
  force `C_12 = C_21` by construction and would report an uninformative
  exact 0 every time.

A real, clean run (this folder's own synthetic data, Case 4 of
`example_4.2.sh`, `--symmetry-method basic`) with everything passing:

```
[3] NUMERICAL QUALITY DIAGNOSTICS
------------------------------------------------------------
Tensor symmetry (C vs C^T): max |C_ij - C_ji| = 0.000 N/m (i=1, j=1)  [OK]
Direction fit quality (stress vs. strain):
  xx: R^2 = 1.000000 | residual stress =    0.00 N/m  [OK]
  xy: R^2 = 1.000000 | residual stress =    0.00 N/m  [OK]
  yy: R^2 = 1.000000 | residual stress =    0.00 N/m  [OK]
Eggbox cross-check (stress vs. energy, same data):
  xx: stress =   200.00 N/m | energy =   200.00 N/m | diff =  0.00%  [OK]
  yy: stress =   200.00 N/m | energy =   200.00 N/m | diff =  0.00%  [OK]
  xy: stress =    80.00 N/m | energy =    80.00 N/m | diff =  0.00%  [OK]
[INFO] All checked direction(s) agree within tolerance -- no indication of eggbox contamination in the stress data.
```

The exact same data analyzed with `--symmetry-method full` instead omits the
`Tensor symmetry` line entirely — not a `0.000` value, no line at all — the
concrete proof this workflow's own `example_4.2.sh` (Case 5) runs live.

**A real bug, found and fixed on this exact material's data during this
session's development**: `eggbox_cross_check`'s energy-side conversion used
to also multiply by the *stress*-side geometric dilution factor (the cell
height `Lz` for 2D, a cross-section for 1D) — a factor the energy side
should never need (total energy isn't diluted by vacuum the way SIESTA's raw
volumetric stress is), silently inflating the reported "energy" value by an
extra factor of `Lz` (a real ~20× error was observed live on this material,
a spurious `~1900%` false "possible eggbox contamination" warning where the
two methods actually agreed to ~4%, section 4). Fixed in `elastic_analysis.py`
— 3D was never affected (its own conversion constant has no such geometric
factor to begin with).

### 3.4 Plotting: gnuplot and matplotlib, both opt-in

Same convention as `4.1-strain`'s own Stage 2: `--save-gnuplot` writes
`<mode>_fit.dat`/`.gplot` per direction/pattern plus a combined overview into
`--plot-dir` (default `elastic_plots/`, render with `cd elastic_plots &&
gnuplot <script>.gplot`); `--view` shows the same fits as an interactive
matplotlib plot instead, points-and-line, no file written. Both off by
default — the `[5]` section reminds you if you skipped them.

### 3.5 Running it both ways

```bash
# Direct CLI
stb-elasticAnalysis --dir elastic_runs --symmetry-method full --no-intro

# Interactive
stb-suite   # then type 4.2.2
```

## 4. Worked example: a real 2D material, end-to-end

Everything above was verified live against a real, fully-converged SIESTA
dataset — the SAME material `structure.fdf`/`calc.fdf` in this folder are
taken from (a nice continuity with `4.1-strain/`'s own worked example,
section 4 there — also a real 2D SiC monolayer): 1 Si + 1 C atom, point group
`-6m2`, run all 3 ways, `+/-2.0%, 5 steps` each.

| Run | Real DFT calcs | Directions run | C11 | C22 | C12 | C66 | Eggbox check | Tensor-symmetry check |
|---|---|---|---|---|---|---|---|---|
| `--method stress --symmetry-method full` | **5** | `xx` only | 200.35 | 200.34 | 37.18 | 81.58 N/m | `xx`: 4.10% [OK] | not run (exact by construction) |
| `--method stress --symmetry-method basic` | **15** | `xx`, `yy`, `xy` | 200.35 | 200.43 | 37.04 | 81.51 N/m | `xx`/`yy`/`xy`: 4.10/4.05/0.06% [OK] | max\|C_ij-C_ji\|=0.280 N/m (i=1,j=2) [OK] |
| `--method energy` | **10** | `xx`, `xy` patterns | 192.13 | 192.12 | 29.01 | 81.56 N/m | n/a (energy method never runs it) | n/a |

Reading this the way sections 1 and 3 set up:

- **Both `stress` runs are internally consistent**: `C66 = (C11-C12)/2`
  closes exactly in both (`(200.35-37.18)/2 = 81.585` vs. the reported
  `81.58`; `(200.35-37.04)/2 = 81.655` vs. `81.51` — small residual
  fit-noise, still `R^2 > 0.998` in both directions per `[3]`) — a genuine
  physical-consistency check the tool never enforces, just reports.
- **`full` really is 3× cheaper here** (section 1.3, section 2.3): 5 real
  SIESTA calculations vs. 15 for the exact same physical answer (`C11`
  agrees to 4 significant figures between the two) — the real-cost version
  of `example_4.2.sh`'s own Case 2.
- **The ~4% gap between the `stress` runs and the `energy` run** traces to a
  real, correctly-flagged residual, not noise or a bug: every run reports
  `[WARNING] Pattern 'xx': linear energy term (3.09e-01 eV) is not negligible
  next to the quadratic term (5.00e+01 eV)` — the reference relaxation
  converged forces well (~2×10⁻⁴ eV/Å, well under its own `MD.MaxForceTol
  0.005`), but still carries a small residual cell stress that shows up
  proportionally larger at this workflow's small (±2%) strain range. A
  useful reminder that "converged forces" and "zero residual stress at this
  specific strain range" are related but distinct claims.
- **The eggbox cross-check (`basic`/`full` rows) independently confirms**
  the same ~4% gap *within* the stress-method data alone (`xx`'s own
  stress-vs-energy diff is 4.10% in both stress runs) — exactly the physical
  observable section 1.5 says it should catch, and exactly why it's
  mandatory rather than opt-in.
- **`basic`'s own tensor-symmetry check passes cleanly** (`0.280 N/m` out of
  `~37 N/m`, well inside noise) — an actual, independently-computed
  agreement between `C12` read off 2 different directions' own data, not
  assumed by construction the way `full`'s own `C12` is.

**When is each method actually the better choice?**

- **`stress + full`** — the cheapest default whenever you trust the detected
  point group (common for well-characterized materials): fewest DFT
  calculations for the complete tensor.
- **`stress + basic`** — costs more DFT time (here, 3×) but buys a real,
  independent internal-consistency check (the tensor-symmetry line) that
  `full` cannot offer by construction. Worth it for an unfamiliar or
  low-symmetry system, or whenever you don't fully trust the auto-detected
  point group and want the data itself to prove the assumed symmetry holds.
- **`--method energy`** — a genuinely independent physical observable
  (energy curvature, not stress), most valuable as a fully separate
  cross-check against an entire `stress` result. The mandatory eggbox check
  (section 3.3) already gets you most of this value *for free* from a single
  `stress` run's own data — reach for a full separate `--method energy` run
  when you want the strongest possible independent confirmation (a
  completely uncorrelated dataset, not just a second fit of the same
  `calc.out` files), or when eggbox flags a real disagreement you want to
  pin down.

## 5. Known, deliberate limitations

- **No `--relax-mode`** (section 2.1) — every generated folder is always a
  pure single-point SCF; there is no clamped-cell-vs-relaxed-transverse-stress
  choice the way `stb-strain` offers, since the imposed strain must stay
  exactly at the sampled geometry for either fitting method to be valid.
- **`--method energy`'s DFT cost scales with symmetry, not always
  favorably** — a genuinely low-symmetry (triclinic/monoclinic) structure
  can need many more combined patterns than either `stress` variant needs
  canonical directions.
- **`[3]` quality checks are `--method stress`-only** — `--method energy`
  never runs the eggbox/fit-quality/tensor-symmetry section at all (it has
  no analogous "2 independent fits of the same data" structure to check).
- **Auto-dimensionality/symmetry both need a readable
  `reference_structure.fdf`** (section 2.4) — without it, dimensionality
  falls back to 3D with a `[WARNING]`, and `--symmetry-method full`/`--method
  energy` hard-fail (they cannot function without symmetry detection);
  `--symmetry-method basic` degrades gracefully to independent
  per-direction fitting instead.

## 6. Step-by-step: running this workflow on your own structure

0. **Prerequisite**: an already-relaxed structure (`structure.fdf`,
   fractional coordinates) and the exact `calc.fdf` used for that relaxation
   — this workflow strains an existing equilibrium structure, it doesn't
   find one for you.

1. **Decide `--method`** (section 1.2): `stress` (default, cheaper, needs
   the mandatory `[3]` checks to trust) or `energy` (independent
   observable, own DFT cost). If `stress`, also **decide
   `--symmetry-method`** (section 1.3): `basic` (default, more DFT calls,
   gets you the tensor-symmetry cross-check) or `full` (fewer DFT calls,
   loses that check).

2. **Run Stage 1** (`stb-elasticInputs`, code `4.2.1`):
   - CLI: `stb-elasticInputs -s structure.fdf -c calc.fdf -p <bank> --dirs
     all --symmetry-method <basic|full> --max <pct> --steps <n>` (or
     `--method energy` instead of `--symmetry-method`).
   - Interactive: `stb-suite` → `4.2.1` — same questions instead of flags.
   - Either way, this writes `elastic_runs/<direction>/strain_<direction>_<pct>/`
     folders (section 2.5) plus one top-level `reference_structure.fdf`
     (section 2.4).

3. **Run SIESTA yourself in every generated folder.** Stage 1 never runs
   SIESTA — `cd` into each `strain_<direction>_<pct>/` and run it like any
   other SIESTA calculation (`siesta calc.fdf > calc.out`, or through your
   own batch/queue script). Repeat for every `<direction>/` subfolder.

4. **Run Stage 2** (`stb-elasticAnalysis`, code `4.2.2`) once every run
   finishes:
   - CLI: `stb-elasticAnalysis --dir elastic_runs --symmetry-method
     <basic|full>` (must match Stage 1's own choice; omit for `--method
     energy`, which always uses its own symmetry-allowed-subspace fit).
   - Interactive: `stb-suite` → `4.2.2`.
   - `--dir` and dimensionality are both auto-detected (sections 2.4/3.2) —
     nothing to set manually in the common case, and no `cd` required.

5. **Read `[3] NUMERICAL QUALITY DIAGNOSTICS`** (section 3.3) before
   trusting the numbers: any eggbox/fit-quality `[WARNING]` is a real,
   specific hint (contaminated stress data, too-large `--max`, an
   incompletely relaxed reference) about exactly which direction to
   re-examine, not a generic "something's wrong somewhere".

6. **Read `[4] STABILITY AND PROPERTIES`** — the Born mechanical-stability
   criteria and the final `VERDICT: STABLE`/`UNSTABLE`; **`[5] SUMMARY &
   FILES`** for the same numbers in one glanceable table plus where
   everything ended up on disk.

7. **(Optional) persist or plot**: `--save-report` writes
   `elastic_stage2.txt`; `--save-gnuplot` writes `<mode>_fit.dat`/`.gplot`
   into `--plot-dir` (section 3.4); `--view` opens the same fits
   interactively in matplotlib instead.

8. **(Optional) go further**: `stb-mlelastic` (ML Simulations) previews the
   exact same physics with a MACE potential, no SIESTA/`strain_*` folders
   needed — a fast sanity check before committing to this real DFT workflow,
   or to decide between `basic`/`full`/`energy` cheaply first (section 9).

## 7. Files in this folder

| File              | What it is                                                          |
|-------------------|----------------------------------------------------------------------|
| `structure.fdf`   | Real, relaxed 2D SiC monolayer, 2 atoms/cell, point group `-6m2`     |
| `calc.fdf`        | The real relaxation `calc.fdf` for `structure.fdf` (DZP/12x12x1/PBE-D3) |
| `example_4.2.sh`  | This walkthrough's runnable script (both stages, 6 cases)            |

## 8. Running the script

```bash
bash example_4.2.sh
```

## 9. What's next

- **`4.1-strain`** (Workflow `4.1`) — the single-direction modulus this
  workflow's full tensor is deliberately not needed for when you only care
  about one direction's stiffness/strength/toughness, not the complete `C_ij`.
- **`stb-mlelastic`** (ML Simulations, code `5.3`) — a fast MACE-based
  preview of the exact same stiffness-tensor physics, no SIESTA/`strain_*`
  folders needed; a good sanity check (or a way to cheaply try `basic` vs.
  `full` vs. `energy` before committing DFT time) before the real workflow
  here.
- **`stb-eosInputs`/`stb-eosAnalysis`** (Workflow `4.18`) / **`stb-mleos`**
  (ML Simulations, code `5.11`) — bulk modulus from the curvature of the
  energy-volume curve instead, an independent cross-check on a 3D
  structure's own `B` (section 3, `[4]`) against this workflow's
  stress-strain-derived value.
- **`stb-cohesive`/`stb_cohesive_analysis`** (Workflow `4.3`) — a different
  mechanical-stability question (cohesive energy per atom) about the same
  kind of structure.
