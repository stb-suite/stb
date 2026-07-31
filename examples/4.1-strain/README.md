# 4.1 — Workflow: Stress-Strain (`stb-strain` / `stb-strainAnalysis`)

This workflow has 2 stages: **Stage 1** (`stb-strain`, code `4.1.1`)
deforms a relaxed structure by a small Cartesian strain and writes
ready-to-run SIESTA folders, one per strain value; **Stage 2**
(`stb-strainAnalysis`, code `4.1.2`) reads the finished SIESTA runs back
and fits a stress-strain curve — initial modulus, peak strength, critical
strain, toughness, and optionally yield. Together they answer one
physical question: **how stiff and how strong is this material along a
given direction, and how far can it be deformed before it breaks?**

Both stages live in this one folder and this one tutorial — not one
folder per tool — because neither stage means much without the other:
Stage 1's output only exists to feed Stage 2, and Stage 2 only exists to
interpret Stage 1's output.

## 1. Theory

### 1.1 Deforming a crystal: the strain tensor

A small, uniform deformation of a crystal lattice is described by a
**deformation gradient** `F = I + ε`, where `ε` is the (small) **strain
tensor** and `I` is the identity — every lattice vector `a` maps to
`F . a`. This is the standard small-strain (linear elasticity)
approximation: valid as long as `ε` stays small enough that the crystal's
response is still reversible and the higher-order terms of a full
nonlinear deformation theory are negligible. `stb-strain` builds `F`
exactly this way (`apply_cartesian_strain`) and applies it to the lattice
vectors only — atomic positions, given in **fractional** coordinates,
scale along with the cell for free, which is exactly why `stb-strain`
requires fractional input and refuses Cartesian coordinates outright (a
Cartesian position would NOT follow the deformed cell automatically).

`stb-strain` supports two families of deformation:

- **Uniaxial** (`x`, `y`, `z`) — `ε` has exactly one nonzero diagonal
  component. Physically: stretching (or compressing) the crystal along
  one Cartesian axis only.
- **Biaxial** (`xy`, `xz`, `yz`) — `ε` has two equal nonzero diagonal
  components, the two named axes strained together by the same amount
  (e.g. `xy` means `ε_xx = ε_yy = ε`, `ε_zz = 0`). This is the natural
  choice for probing a 2D material's **in-plane** response, since a
  single in-plane axis alone doesn't respect the material's own in-plane
  symmetry the way straining both together does.

Stress and strain components throughout this workflow (and the rest of
this suite) use the standard **Voigt notation**: the 6 independent
components of a symmetric 3×3 tensor collapsed into one index,
`1..6 = XX/YY/ZZ/YZ/XZ/XY` — the same numbering SIESTA itself uses for
`%block Geometry.Constraints` and its own `Stress tensor Voigt` output
line.

### 1.2 What the "initial slope" measures — and what it doesn't

Stage 2 fits a straight line to the stress-strain data at small strain
and reports its slope as the direction's **initial modulus** — a
stiffness, in the same spirit as Young's modulus (`stress = modulus ×
strain`, Hooke's law, in the linear-elastic small-strain limit). But it
is **not**, in general, the same number as a rigorous elastic constant
`C_ij` from the full elastic tensor:

- With `--relax-mode cell-fixed`, every OTHER Voigt component is held
  exactly fixed while only the imposed one is strained — a "clamped
  -cell" measurement.
- With `--relax-mode stress-constrained`, the other periodic components
  are left free to relax to zero stress instead — a "relaxed-transverse
  -stress" measurement.

Both are legitimate, standard DFT conventions (see `stb-strain`'s own
`--relax-mode` docs), and both differ from the *fully coupled* response a
true `C_ij` requires — that needs the sibling workflow
**`stb-elasticInputs`/`stb-elasticAnalysis`** (elastic constants), which
fits a whole system of coupled strain directions simultaneously.
`stb-strainAnalysis`'s own report says so explicitly:
`[INFO] Single-direction slope under a clamped transverse cell -- for the
rigorous small-strain elastic-tensor Young's Modulus, use
stb-elasticAnalysis.`

Think of this workflow's own modulus as **"how stiff is this material if
you pull on it exactly like this, under exactly this transverse
boundary condition"** — directly meaningful and reproducible in its own
right, just not automatically equal to a textbook single-crystal elastic
constant unless the two conventions happen to coincide (they do, for the
right symmetry and relax-mode combination).

### 1.3 Reading a stress-strain curve: modulus, UTS, critical strain, toughness, yield

A typical stress-strain curve for a solid rises roughly linearly at small
strain, curves over, reaches a peak, and then falls (or plateaus) as the
material's bonding network is pushed past its ability to sustain load.
Stage 2 extracts 4 numbers (and a 5th, optional one) that summarize this
shape:

| Quantity | What it is | Physical meaning |
|---|---|---|
| **Initial Slope** | linear fit at small strain (section 1.2) | stiffness — how much stress it takes to produce a given small strain |
| **Peak Stress/Force (UTS)** | the single largest value on the curve | the "ultimate tensile strength" — the most load-bearing capacity this direction ever exhibits in the tested range |
| **Critical Strain** | the strain at which the peak occurs | how far you can deform the material before it starts to soften/fail |
| **Toughness** | the area under the whole curve (`∫ stress d(strain)`) | energy absorbed per unit volume (3D) / area (2D) / length (1D) before failure — a fracture-resistance measure, distinct from stiffness or strength |
| **Yield (0.2% offset)**, opt-in via `--yield` | where the curve departs from the initial-slope line by a strain offset of 0.002 | a **macroscopic metallurgy concept** (onset of dislocation-glide plasticity) — reported for comparison, but a perfect, defect-free periodic crystal under affine DFT strain has no dislocation-nucleation mechanism to actually yield through; interpret with care |

If the peak sits at either edge of the strain range you actually tested
(`--stmin`/`--stmax`), Stage 2 flags it explicitly — the true peak may
lie beyond what you scanned, and the reported UTS/critical-strain values
are only a **lower bound** in that case.

### 1.4 Dimensionality changes the physical quantity, not just the unit

SIESTA always reports a genuine 3D stress (force per unit area), computed
from the **full simulation cell**, vacuum included. That's the right
quantity as-is for a bulk (3D) structure — but not for a 2D sheet or a 1D
wire, where the vacuum padding is an arbitrary simulation-box choice, not
a physical thickness:

- **3D (bulk)** — Stress in **GPa** (`1 kBar = 0.1 GPa`), used directly.
- **2D (a slab/sheet, vacuum along `z`)** — SIESTA's own stress is
  diluted by the *entire* cell height, vacuum included, which is
  arbitrary (you could always add more vacuum and get a smaller number
  for the same physical sheet). Multiplying back by that same cell
  height cancels the arbitrary dilution and recovers the physically
  meaningful quantity: **force per unit width, in N/m** — exactly the
  convention the 2D-materials literature itself uses (e.g. graphene's
  own in-plane stiffness is conventionally quoted in N/m, around
  ~340 N/m, not as a "3D stress in GPa" that depends on an arbitrarily
  chosen layer thickness).
- **1D (a wire/tube, vacuum along `x` and `y`)** — same problem, one
  dimension worse: there's no natural "cross-section" at all without
  assuming a wire wall thickness. Instead of a made-up stress, Stage 2
  reports the raw **axial Force, in nN** (`stress × cell cross-section`,
  independent of how much vacuum padding was chosen). If you *do* know
  the wire's real physical cross-section, `--cross-section` additionally
  back-derives a "conventional" GPa-equivalent stress for comparison.

Stage 2 auto-detects which of these three applies (section 3.2) instead
of asking you to already know it.

### 1.5 Symmetry reduces how many directions you actually need

By Neumann's principle, if some symmetry operation of the crystal's point
group maps one direction onto another, straining either one must give
the identical mechanical response — computing both would just repeat the
same physics twice. `stb-strain`'s `[3] AXIS SYMMETRY ADVISORY` (uniaxial)
and the interactive menu's full direction-picker (both uniaxial and
biaxial, section 2.5) detect this automatically from the structure's own
point group, so you only ever need to actually run the symmetry-distinct
directions.

## 2. Stage 1: generating strained structures (`stb-strain`, code `4.1.1`)

### 2.1 What it does

`stb-strain` takes a relaxed structure and a reference `calc.fdf`, applies
the deformation from section 1.1 to the lattice vectors for a range of
strain values, and writes one ready-to-run folder per value: the
deformed structure, a copy of `calc.fdf`, a new `config_extra.fdf` (see
2.3), and any linked pseudopotentials. It does **not** run SIESTA itself
— you run each folder yourself, then hand the results to Stage 2
(section 3).

### 2.2 The `calc.fdf` you pass in — step by step

**Golden rule: `--calc` must be the exact same file used for the
structure's original relaxation** — same basis set, k-points, XC
functional, SCF settings. `stb-strain` doesn't regenerate any of that; it
only strains the lattice and adds one `%include` line (section 2.3), so
the physics of the strained evaluation stays consistent with whatever
produced `structure.fdf` in the first place. This example's `calc.fdf` is
exactly that: a plausible relaxation input for `structure.fdf`'s 8-atom
bulk Si cell (DZP basis, an 8x8x8 Monkhorst-Pack grid, PBE, standard SCF
tolerances).

Two settings in that file matter specifically for Stage 1, and
`stb-strain`'s own `[4]` report section checks (and warns about) both,
**without ever forcing either one**:

1. **`MD.VariableCell` stays `true`.** `stb-strain` itself never disables
   cell relaxation — that job belongs to the `%block Geometry.Constraints`
   it writes into the new `config_extra.fdf` (section 2.3), which only
   has *any* effect while SIESTA's cell is actually variable in the first
   place. Turn `MD.VariableCell` off (or leave it out) and the constraints
   block becomes a no-op — SIESTA never even looks at which stress
   components are fixed, since the cell was never going to move anyway.
2. **`MD.Steps` stays non-zero** (this is the real SIESTA keyword — the
   older `MD.NumCGsteps` spelling is also recognized as a fallback). This
   is the ionic-relaxation step count; with `0`, nothing relaxes no matter
   which `--relax-mode` you pick.

`calc.fdf` in this folder has both right (`MD.VariableCell true`,
`MD.Steps 100`) — open it and read the comment above its `STRUCTURE
RELAXATION` block. `calc_missing_relax.fdf`, also in this folder, is the
same file with both **deliberately wrong** (no `MD.VariableCell`,
`MD.Steps 0`), used only to demonstrate the warnings below live — never
use it as a real template.

#### The gotcha, live

```
$ stb-strain -s structure.fdf -c calc_missing_relax.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro
...
[4] RELAXATION MODE & CELL CONSTRAINTS
Calc template (current state, read-only -- nothing below is forced):
  MD.TypeOfRun=CG  Steps: MD.Steps=0  MD.VariableCell=(absent)
[WARNING] MD.VariableCell is not enabled in this calc.fdf (or is absent).
The %block Geometry.Constraints written below only has an effect while
the cell is actually variable -- add 'MD.VariableCell true' to your
--calc yourself; this tool does not force it.
[WARNING] Relaxation step count is MD.Steps=0 -- with 0 (or no) steps,
atomic/cell relaxation is a no-op regardless of --relax-mode. Set
MD.Steps (or MD.NumCGsteps) > 0 in your --calc for this mode to actually
relax anything.
```

Both warnings fire together, and the run still completes (they're
warnings, not hard errors — `stb-strain` can't force a directive that
lives entirely in *your* file, so it tells you what's wrong instead).

### 2.3 `--relax-mode`: two conventions, one mechanism

Both `--relax-mode` choices (introduced physically in section 1.2) are
expressed through the exact same mechanism: a `%block Geometry.Constraints`
with `stress N` lines (Voigt numbering, section 1.1), written into a new
`config_extra.fdf` file and included via `%include config_extra.fdf`
inserted right after `%include structure.fdf` in every generated
`calc.fdf` copy. Nothing else in `calc.fdf` is touched. The 2 modes only
differ in **which** components that block fixes:

- **`cell-fixed`** — fixes all 6 components. The cell stays exactly at
  the imposed strain (the same practical effect as disabling
  `MD.VariableCell` entirely, just expressed as constraints so both modes
  share one mechanism); only the ions relax, via `calc.fdf`'s own
  `MD.Steps`.
- **`stress-constrained`** — fixes only the imposed strain direction's own
  component(s); every other periodic direction is left free to relax to
  zero stress. (If the structure had a vacuum-padded axis — a 2D slab, say
  — that axis's own component would *also* be fixed automatically, to
  protect it from spurious relaxation; this bulk 3D example has no vacuum
  axis, so that rule never triggers here.)

#### Side by side

```
$ stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x ...
config_extra.fdf:
  %block Geometry.Constraints
    stress 1  # Fixes XX
    stress 2  # Fixes YY
    stress 3  # Fixes ZZ
    stress 4  # Fixes YZ
    stress 5  # Fixes XZ
    stress 6  # Fixes XY
  %endblock Geometry.Constraints

$ stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir x ...
config_extra.fdf:
  %block Geometry.Constraints
    stress 1  # Fixes XX
  %endblock Geometry.Constraints
```

`stress-constrained` is the physically softer constraint (only 1 of the 6
components fixed here) — the standard "relaxed-ion, relaxed-transverse
-stress" method for elastic-constant work; `cell-fixed` is the stricter
"clamped-cell, relaxed-ion" method. Both are legitimate; which one you
want depends on which elastic-constant convention you're computing.

### 2.4 Symmetry reduces the direction list (CLI advisory)

`structure.fdf` is cubic bulk silicon (point group `m-3m`) — one of the
few point groups where every Cartesian axis really is mechanically
equivalent (section 1.5). `stb-strain`'s `[3] AXIS SYMMETRY ADVISORY`
section detects this automatically and tells you so, before you waste a
calculation straining `y` or `z` too:

```
[3] AXIS SYMMETRY ADVISORY
  Detected symmetry : point group m-3m -- 48 operation(s) ...
  Axis    Status        Equivalent to
  x       REQUESTED     --
  y       EQUIVALENT    x
  z       EQUIVALENT    x
  y, z are equivalent to 'x' by symmetry -- straining them should give
  the same mechanical response; you may not need to compute both.
```

Purely advisory — it never blocks the run, and biaxial directions are out
of scope for this check (see `core/symmetry.py`'s own docstring for why).

### 2.5 The interactive menu's full direction-picker

The `[3] AXIS SYMMETRY ADVISORY` above is CLI-only, uniaxial-only, and
purely informational — it tells you `y`/`z` are redundant but doesn't stop
you from running them anyway, and it has nothing to say about biaxial
directions at all. The interactive menu (`stb-suite` → `4.1.1`) does the
equivalent check for **both** uniaxial and biaxial directions, together
with vacuum-axis detection, and turns it into an actual choice instead of
just a warning: it detects everything up front, then only offers the
directions that are physically valid (non-vacuum) and not
symmetry-redundant.

Biaxial equivalence is derived from the exact same axis-symmetry groups
the CLI's own advisory uses, just keyed by each biaxial direction's
*excluded* axis instead of its own axes — `xy` excludes `z`, so it's
equivalent to another biaxial direction whenever a point-group operation
maps their excluded axes onto each other. No separate biaxial symmetry
computation is needed.

On this folder's own `structure.fdf` (cubic Si, `m-3m`, no vacuum), the
menu collapses from 6 possible directions down to 2 independent ones:

```
------------------------------------------------------------
            STRAIN DIRECTION SELECTION
------------------------------------------------------------
  Detected symmetry : point group m-3m -- 48 operation(s) ...
------------------------------------------------------------
  Direction  Type      Status        Notes
  --------------------------------------------------------
  x          uniaxial  INDEPENDENT   --
  y          uniaxial  REDUNDANT     equivalent to 'x'
  z          uniaxial  REDUNDANT     equivalent to 'x'
  xy         biaxial   INDEPENDENT   --
  xz         biaxial   REDUNDANT     equivalent to 'xy'
  yz         biaxial   REDUNDANT     equivalent to 'xy'
  --------------------------------------------------------
  2 independent, non-vacuum direction(s): x, xy
------------------------------------------------------------

  1) x
  2) xy
  3) ALL independent UNIAXIAL directions (x)
  4) ALL independent BIAXIAL directions (xy)
  5) ALL independent directions -- uniaxial + biaxial (2 total)

Select a direction [1-5]:
```

Picking `5` runs Stage 1 once per independent direction (`x` and `xy`
here), each into its own `strain_runs/<direction>/` subfolder (section
2.6) — with no separate y/n follow-up prompt needed, and without ever
wasting a real calculation on `y`, `z`, `xz`, or `yz`. `3` and `4` are the
same idea scoped to just one strain type — handy when a structure has
several independent directions of each kind and you only want, say,
every uniaxial one.

On a vacuum-padded structure (a 2D slab, say), any direction touching the
vacuum axis is marked `VACUUM` instead and dropped from the menu entirely
— `stb-strain`'s own CLI would otherwise hard-error if you asked for one
of those by hand. If the structure can't be read, or the symmetry
pre-check fails for any reason, the menu falls back to the old free-text
prompt (`x`, `y`, `z`, `xy`, `xz`, `yz`) with no filtering.

### 2.6 Output layout: one subfolder per direction

Every direction gets its own `<output-dir>/<direction>/` subfolder —
`strain_runs/x/strain_x_0.00/`, `strain_runs/x/strain_x_2.00/`, and so on
— rather than everything flat directly under `strain_runs/`. Run a second
direction (`y`) into the same `--output-dir` and it gets its own
`strain_runs/y/` alongside, with no collision and no risk of one
direction's `--save-report` file overwriting another's. This is exactly
the layout Stage 2 (section 3) expects: point `--dir` at `strain_runs`
itself and every direction found underneath is read and compared
automatically.

### 2.7 Running it both ways

Every command below also works through the interactive menu —
`stb-suite` → `4.1.1` (or type `4.1.1` directly from the main prompt) asks
the same questions (structure file, calc file, relax mode, direction,
strain range, ...) instead of flags, and calls the exact same
`stb-strain` underneath. `example_4.1.sh` proves this directly by running
both paths with the same inputs and diffing the generated folders.

```bash
stb-strain -s structure.fdf -c calc.fdf \
    --relax-mode cell-fixed --stdir x --stmin 0 --stmax 2 --step 2
```

## 3. Stage 2: analyzing the results (`stb-strainAnalysis`, code `4.1.2`)

### 3.1 Report structure

Once you've run SIESTA in every `strain_runs/<direction>/strain_<direction>_<pct>/`
folder from Stage 1, Stage 2 reads the finished `.out` files back and
computes everything from section 1.3 (modulus, peak, critical strain,
toughness, optional yield). It reports in the same numbered `[N]` section
style as Stage 1's own report:

```
[0] RUN METADATA          -- the requested directory/file/dimensionality/output settings
[1] INPUT DATA            -- which strain_* folders were found, grouped by direction
[2] DIMENSIONALITY DETECTION
[3] MECHANICAL PROPERTIES -- a quick-reference table (if >1 direction) + one full block per direction
[4] OUTPUT FILES          -- what --save-gnuplot would write / did write
[5] SUMMARY & NEXT STEPS
```

### 3.2 Dimensionality auto-detection

`[2]` looks for a real `.fdf` structure file inside the first `strain_*`
folder found (Stage 1 always writes one — the deformed structure, same
basename you gave `--structure`) and classifies it via vacuum-axis
detection, applying section 1.4's rules exactly like Stage 1's own `[2]`
section and `stb-mlelastic`'s `--dimensionality auto`: vacuum along `z`
only → 2D (`N/m`); vacuum along `x`+`y` → 1D (`nN`); no vacuum → 3D
(`GPa`, this folder's own case). If no `.fdf` is found at all (a
hand-built folder holding only a bare SIESTA `.out`) it falls back to 3D
with a `[NOTE]`, not an error — pass `--dimensionality` explicitly
(`auto`/`3d`/`2d`/`1d`) to override either way.

### 3.3 Multiple directions are compared automatically, no flag needed

Point `--dir` at Stage 1's top-level output directory (`strain_runs`,
holding `x/`, `xy/`, ... subfolders — e.g. from picking Stage 1's menu
"ALL" option) and `[1]`/`[3]` report every direction found, side by side,
the same way `--dir strain_runs/x` alone still reports just that one
direction on its own.

### 3.4 The initial-slope fit window — and a real bug fixed here

`[3]` reports not just the modulus but the fit window itself (how many
points, what strain range) and the fit's standard error, e.g.
`Fit window : 3 point(s), strain 0.0000% to 4.0000%`. The fit (section
1.2/1.3) prefers every point within ±2% strain; when a sweep's own step
size is coarser than that (very common for a large-deformation sweep,
e.g. this workflow's own default 25% range at a few-percent step), fewer
than 3 points land inside ±2% — the window widens to the nearest few
points instead, with an explicit `[NOTE]` saying so.

**A real bug was found and fixed here, verified live**: the *previous*
version of this fallback fit the *entire* sweep (0% all the way to
whatever `--stmax` was) whenever too few points fell in ±2%, instead of
widening to just the nearest few — silently averaging the genuine
small-strain elastic response together with large-strain, manifestly
nonlinear/plastic behavior into one number. On a real 2D SiC monolayer
dataset (3 directions, 21 steps each, 0-40% strain, 2%/step — so only the
0%/2% points ever fell in ±2%), that gave a "modulus" with R² as low as
~0.01 (statistically meaningless) instead of the correct near-zero-strain
slope (R² > 0.99 once the window widens to the nearest 3 points instead;
see the worked example, section 4). The fit line drawn on both the
gnuplot and matplotlib output is now also clipped to the actual fit
window (never drawn across the full strain range, which would otherwise
visually claim the material follows that straight line all the way to
its peak/plastic regime).

### 3.5 Plotting: gnuplot and matplotlib, both opt-in, points-and-line

`--save-gnuplot` writes `<direction>_curve.dat`/`.gplot` (or
`comparison_curve.dat`/`.gplot` for multiple directions) — the same
`.dat`+`.gplot` convention used throughout this suite; the `.gplot` script
references its own data/output files by bare filename, so run it from
*inside* `--output-dir` (`cd <output-dir> && gnuplot <script>.gplot` —
Stage 2's own `[4]` section prints this exact command). `--view` opens an
interactive matplotlib plot of the same curve instead — shown on screen
only, never written to disk. Both draw the stress/force-strain data as
connected points-and-line (not bare scatter points), the small-strain fit
line (dashed, clipped to its own fit window — section 3.4), and the peak
point. `--save-report` persists the numbered report itself to
`stb_strainAnalysis_report.txt`. All three default off.

```
$ stb-strainAnalysis --file calc.out --dir strain_runs --save-report --save-gnuplot
...
[3] MECHANICAL PROPERTIES
------------------------------------------------------------
2 direction(s) found -- comparing automatically (no flag needed). Full
detail per direction follows; quick-reference table first:
Direction | Slope (GPa) | R²     | Peak (GPa) | Crit. Strain (%) | Toughness (GJ/m^3) | Notes
-----------------------------------------------------------------------------------------------
X         | ...         | ...    | ...        | ...              | ...                | --
XY        | ...         | ...    | ...        | ...              | ...                | --

Direction       : X
...
  Fit window    : 3 point(s), strain 0.0000% to 4.0000%
...
```

### 3.6 Running it both ways

Same as Stage 1 (section 2.7): `stb-suite` → `4.1.2` asks the same
questions (directory, output filename, dimensionality — `auto`/`3d`/`2d`/
`1d`, advanced settings, yield, output directory, save-report/
save-gnuplot/view) instead of flags, and calls the exact same
`stb-strainAnalysis` underneath.

## 4. Worked example: a real 2D material, end-to-end

Everything above was verified live against a real, fully-converged SIESTA
dataset: a 2D SiC monolayer (fetched via `stb-fetch` from an OPTIMADE
provider, then run through this exact workflow), 3 directions (`x`, `y`,
`xy`), 21 strain steps each, 0-40% at 2%/step, all 63 SIESTA runs
converged normally.

```
$ stb-strainAnalysis --file calc.out --dir strain_runs --save-gnuplot
...
[2] DIMENSIONALITY DETECTION
[INFO] Auto-detected from a real structure file -- dimensionality: 2D
Dimensionality    : 2D (auto-detected)
[INFO] Detected cell Z-height: 20.7697 Ang (used for N/m conversion)

[3] MECHANICAL PROPERTIES
Direction | Slope (N/m) | R²     | Peak (N/m) | Crit. Strain (%)
X         | 154.3192    | 0.9975 | 19.5423    | 20.0000
XY        | 199.1299    | 0.9968 | 15.7829    | 16.0000
Y         | 167.3003    | 0.9999 | 19.7770    | 20.0000
```

Reading this the way section 1 sets up:

- **Dimensionality auto-detected correctly** (section 3.2) — 2D, purely
  from the real `structure.fdf`'s vacuum gap along `z`, no
  `--dimensionality` flag needed. The Z-height (20.7697 Å) is this
  structure's own vacuum-inclusive cell height, auto-read from the SIESTA
  output to un-dilute the reported stress into N/m (section 1.4).
- **All 3 directions compared automatically** (section 3.3) — no
  `--compare`-style flag exists or was needed.
- **The fit-window fix in action** (section 3.4): every direction here
  only has 2 points within the literal ±2% window (0% and 2%), so the fit
  widened to the nearest 3 points (0%, 2%, 4%) instead of the old, buggy
  full 0-40% fit — recovering R² of 0.9968-0.9999 (excellent linear fits)
  where the old code would have given R² as low as ~0.01.
- **Physical sanity check**: 154-199 N/m is the right *order of
  magnitude* for a covalent 2D monolayer's in-plane stiffness — smaller
  than graphene's own ~340 N/m (expected: SiC's Si-C bonds are weaker and
  longer than graphene's C-C bonds), and the modest anisotropy between
  `X`/`Y`/`XY` (not all equal) is itself informative: it means this
  structure's in-plane point group is lower than a fully hexagonal
  (isotropic) one — exactly the kind of fact this workflow is built to
  surface, not something you'd know just from looking at the lattice
  vectors.
- **Peaks around 16-20% critical strain**, well past the small-strain
  elastic regime the modulus describes — a reminder that "Initial Slope"
  and "Peak" describe two different physical regimes of the *same* curve
  (section 1.3), not the same number scaled differently.

`example_4.1.sh`'s own Stage 2 cases demonstrate the same mechanics
against small, hand-built `calc.out` data (since this folder's own Stage
1 run doesn't invoke real SIESTA) — the numbers there are illustrative,
not a real material's properties; this section is the genuine, physically
verified result.

## 5. Known, deliberate limitations

- **Not a full elastic tensor** (section 1.2) — a single-direction,
  clamped-or-relaxed-transverse modulus, not a rigorous `C_ij`. Use
  `stb-elasticInputs`/`stb-elasticAnalysis` for that (or `stb-mlelastic`
  for a fast MACE-based preview of the same physics).
- **Biaxial "modulus" is not a true biaxial modulus** — Stage 2 reports
  only the stress component along the biaxial direction's first axis
  (e.g. `xy` → the `xx` component), flagged with a `[WARNING]`.
- **Yield strength (0.2% offset)** is a macroscopic metallurgy concept
  (dislocation-glide onset) — a defect-free periodic crystal under
  affine DFT strain has no such mechanism; report with care (section
  1.3).
- **Peak/critical strain at the edge of the tested range** is only a
  lower bound on the true peak — widen `--stmin`/`--stmax` if you see
  this flag.
- **Auto-dimensionality-detection** needs a real, parseable `.fdf` in the
  scanned folder; a hand-built/legacy folder holding only a bare `.out`
  falls back to 3D (section 3.2).

## 6. Files in this folder

| File                       | What it is                                                        |
|-----------------------------|--------------------------------------------------------------------|
| `structure.fdf`             | Bulk Si, 8-atom conventional cubic cell, fractional coordinates    |
| `calc.fdf`                  | The (correct) relaxation calc.fdf for `structure.fdf`              |
| `calc_missing_relax.fdf`    | Same file, deliberately misconfigured — for the gotcha demo only   |
| `example_4.1.sh`            | This walkthrough's runnable script (both stages)                   |

## 7. Running the script

```bash
bash example_4.1.sh
```

## 8. What's next

- **`stb-elasticInputs`/`stb-elasticAnalysis`** (Workflow `4.2`) — the
  rigorous full elastic tensor `C_ij` this workflow's own modulus is
  deliberately NOT a substitute for (section 1.2/5).
- **`stb-mlelastic`** (ML Simulations) — a fast MACE-based preview of the
  same stiffness-tensor physics, no SIESTA/`strain_*` folders needed; a
  good sanity check before committing to the real DFT workflow here.
- **`stb-cohesive`** (Workflow `4.3`) — a different mechanical-stability
  question (cohesive energy per atom) about the same kind of structure.
- **`stb-mleos`**/**`stb-eosInputs`+`stb-eosAnalysis`** — bulk modulus
  from the curvature of the energy-volume curve instead, an independent
  cross-check on a 3D structure's stiffness against this workflow's own
  stress-strain-derived value.
