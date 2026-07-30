# 4.1 — Workflow: Stress-Strain (`stb-strain` / `stb-strainAnalysis`)

This workflow has 2 stages: **Stage 1** (`stb-strain`, code `4.1.1`)
applies a Cartesian strain to a structure and writes ready-to-run SIESTA
folders; **Stage 2** (`stb-strainAnalysis`, code `4.1.2`) aggregates the
finished SIESTA runs back into a stress-strain curve and mechanical
properties. Both stages live in this one folder (one workflow, one
tutorial) instead of one folder per tool. **This part only covers Stage
1** — Stage 2 will be added to this same `README.md`/`example_4.1.sh`
in a follow-up update.

## 1. What Stage 1 does

`stb-strain` takes a relaxed structure and a reference `calc.fdf`, applies
a small Cartesian deformation (`F = I + strain_tensor`, standard small
-strain elasticity) to the lattice vectors for a range of strain values,
and writes one ready-to-run folder per value: the deformed structure, a
copy of `calc.fdf`, a new `config_extra.fdf` (see section 3), and any
linked pseudopotentials. It does **not** run SIESTA itself — you run each
folder yourself, then hand the results to Stage 2.

The deformation can be uniaxial (`x`, `y`, `z`) or biaxial (`xy`, `xz`,
`yz` — both listed axes strained together, e.g. for an in-plane 2D
scan); this example only exercises uniaxial `x`.

## 2. The `calc.fdf` you pass in — step by step

**Golden rule: `--calc` must be the exact same file used for the
structure's original relaxation** — same basis set, k-points, XC
functional, SCF settings. `stb-strain` doesn't regenerate any of that; it
only strains the lattice and adds one `%include` line (section 3), so the
physics of the strained evaluation stays consistent with whatever produced
`structure.fdf` in the first place. This example's `calc.fdf` is exactly
that: a plausible relaxation input for `structure.fdf`'s 8-atom bulk Si
cell (DZP basis, an 8x8x8 Monkhorst-Pack grid, PBE, standard SCF
tolerances).

Two settings in that file matter specifically for Stage 1, and
`stb-strain`'s own `[4]` report section checks (and warns about) both,
**without ever forcing either one**:

1. **`MD.VariableCell` stays `true`.** `stb-strain` itself never disables
   cell relaxation — that job belongs to the `%block Geometry.Constraints`
   it writes into the new `config_extra.fdf` (section 3), which only has
   *any* effect while SIESTA's cell is actually variable in the first
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

### The gotcha, live

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

## 3. `--relax-mode`: 2 modes, one mechanism

Both `--relax-mode` choices are expressed through the exact same
mechanism: a `%block Geometry.Constraints` with `stress N` lines (SIESTA's
own Voigt numbering — `1..6` = `XX/YY/ZZ/YZ/XZ/XY`), written into a new
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

### Side by side

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

## 4. Symmetry: why `x`, `y`, `z` are equivalent here

`structure.fdf` is cubic bulk silicon (point group `m-3m`) — one of the
few point groups where every Cartesian axis really is mechanically
equivalent. `stb-strain`'s `[3] AXIS SYMMETRY ADVISORY` section detects
this automatically and tells you so, before you waste a calculation
straining `y` or `z` too:

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

## 5. The interactive menu goes further: a full direction-picker (`stb-suite` → `4.1.1`)

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
here), each into its own `strain_runs/<direction>/` subfolder (section 6)
— with no separate y/n follow-up prompt needed, and without ever wasting a
real calculation on `y`, `z`, `xz`, or `yz`. `3` and `4` are the same idea
scoped to just one strain type — handy when a structure has several
independent directions of each kind and you only want, say, every
uniaxial one.

On a vacuum-padded structure (a 2D slab, say), any direction touching the
vacuum axis is marked `VACUUM` instead and dropped from the menu entirely
— `stb-strain`'s own CLI would otherwise hard-error if you asked for one
of those by hand. If the structure can't be read, or the symmetry
pre-check fails for any reason, the menu falls back to the old free-text
prompt (`x`, `y`, `z`, `xy`, `xz`, `yz`) with no filtering.

## 6. Output layout: one subfolder per direction

Every direction gets its own `<output-dir>/<direction>/` subfolder —
`strain_runs/x/strain_x_0.00/`, `strain_runs/x/strain_x_2.00/`, and so on
— rather than everything flat directly under `strain_runs/`. Run a second
direction (`y`) into the same `--output-dir` and it gets its own
`strain_runs/y/` alongside, with no collision and no risk of one
direction's `--save-report` file overwriting another's. This also sets up
Stage 2 to read every direction found under `strain_runs/` at once — more
on that once Stage 2 is added here.

## 7. Running it both ways

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

## 8. Files in this folder

| File                       | What it is                                                        |
|-----------------------------|--------------------------------------------------------------------|
| `structure.fdf`             | Bulk Si, 8-atom conventional cubic cell, fractional coordinates    |
| `calc.fdf`                  | The (correct) relaxation calc.fdf for `structure.fdf`              |
| `calc_missing_relax.fdf`    | Same file, deliberately misconfigured — for the gotcha demo only   |
| `example_4.1.sh`            | This walkthrough's runnable script (Stage 1 only, for now)         |

## 9. Running the script

```bash
bash example_4.1.sh
```

## 10. What's next

Stage 2 (`stb-strainAnalysis`, code `4.1.2`) — reading the `strain_runs/`
folders back after you've run SIESTA in each one, fitting the
stress-strain curve, and reporting mechanical properties — will be added
to this same `README.md` and `example_4.1.sh` in a follow-up update.
