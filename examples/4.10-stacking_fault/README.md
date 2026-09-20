# 4.10 — Workflow: 2D Stacking Fault / Gamma-Surface (`stb-stackingfault` / `stb-stackingfaultAnalysis`)

Slide one layer of a van der Waals bilayer (or heterostructure) laterally
over the other, recompute the total energy at a grid of offsets, and you
get a 2D energy landscape — the **generalized stacking-fault energy
surface**, or **gamma surface**. This item has **two** stages:

- **Stage 1 (`stb-stackingfault`, code `4.10.1`)** rigidly slides layer 2
  across a grid (or a 1D line) of lateral offsets and writes one SIESTA
  folder per point, with a choice of 3 strategies for the interlayer gap.
- **Stage 2 (`stb-stackingfaultAnalysis`, code `4.10.2`)** reads every
  finished energy back and reports the equilibrium stacking, the
  highest-energy registry, the corrugation, and the full gamma-surface
  map.

`example_4.10.sh` is a guided, runnable walkthrough of all of this,
including a real 30-point energy landscape (Case 5 / Section 4 below) and
a live proof that the CLI and the interactive `stb-suite` menu produce
byte-identical output (Case 7).

**A 3rd stage used to exist here** (`stb-stackingfaultBsse`, a BSSE/
counterpoise correction). It was removed — Section 1.2 and Section 5
explain why in detail, since it's the single most important thing to
understand about this workflow if you're coming from `stb-adsorbBsse` or
a cohesive-energy calculation and expect a similar correction step here.

## 1. Theory

### 1.1 The gamma-surface and its two headline numbers

For every sampled lateral offset `(shift_x, shift_y)` (fractional
coordinates of layer 2's in-plane lattice vectors), Stage 1 writes a
single structure with layer 1 fixed and layer 2 rigidly translated, and
Stage 2 reads back its total energy `E(shift_x, shift_y)`. Two numbers
summarize the resulting landscape:

- **Equilibrium stacking**: the offset with the lowest energy — e.g.
  graphene's `AB` (Bernal) registry, one of the 2 local minima a
  honeycomb bilayer has (the other, symmetric one is `BA`).
- **Corrugation**: `max(E) - min(E)` over the sampled grid — the energy
  barrier resisting interlayer sliding/shear. Directly relevant to
  interlayer friction, and (for a 3D crystal's actual stacking fault, as
  opposed to a 2D vdW bilayer) to dislocation dissociation energetics.

Reported in **meV/Å²** (millielectronvolt per square Ångström of
interface area) — the standard, cell-size-independent unit for this kind
of areal energy density, computed as `(E - E_min) * 1000 / area`, where
`area` is the in-plane cell area (`|a × b|` of the structure's own
lattice vectors — see Section 3.2).

### 1.2 Why there is no BSSE-correction stage here (read this if you know `stb-adsorbBsse`)

Basis Set Superposition Error (BSSE) correction is necessary when you
compare energies of systems that have **genuinely different atom counts
and basis sizes** — e.g. adsorption energy (slab+adsorbate vs. slab alone
vs. adsorbate alone) or cohesive energy (bulk crystal vs. isolated atom).
In those cases, the "combined" system's SCF benefits from more basis
functions than either isolated reference does on its own, artificially
lowering it — counterpoise correction (real+ghost single-point
references at the same geometry) fixes that asymmetry.

**A stacking-fault scan has no such asymmetry.** Every single point in
the whole grid has the *exact same atoms, the exact same basis set* —
only positions change. Write `E(shift) = E_layer1_isolated +
E_layer2_isolated + E_interaction(shift)`. Since neither layer's own
internal geometry changes with `shift` (only layer 2's *position*
changes, not its internal structure), `E_layer1_isolated` and
`E_layer2_isolated` are **constants**, identical at every grid point. They
cancel exactly in `max(E) - min(E)` — the raw total energy's corrugation
**is already, exactly**, the true interaction-energy corrugation. No
correction needed, or even meaningful, to apply here.

Applying the `stb-adsorbBsse`-style formula anyway
(`E_corrected = E - e_layer1(ghost) - e_layer2(ghost)`) doesn't just do
nothing useful — it makes things *worse*: `e_layer1(ghost)`/
`e_layer2(ghost)` (real+ghost single-points) are **lower** than the true
isolated references by exactly the BSSE amount, so subtracting them
*adds the BSSE signal back into the result*
(`E_corrected(shift) = E_interaction(shift) + BSSE_total(shift)`,
provably, see Section 5). This was verified live, on a real 30-point
scan, before the correction stage was removed from this workflow: the
"corrected" corrugation curve showed spurious kinks the raw energy curve
didn't have, traced to exactly this effect.

## 2. Stage 1: generating the grid (`stb-stackingfault`, code `4.10.1`)

### 2.1 Sweep shape: `--scan surface` (default) vs. `x`/`y`/`xy`

`--scan surface` (default) samples the full 2D grid, `-nx`/`-ny` points
along each in-plane lattice vector independently (an anisotropic lattice
often warrants a different density per axis). `--scan x`/`y`/`xy` trade
that for a cheaper 1D line — along `shift_x` only, `shift_y` only, or the
`shift_x = shift_y` diagonal — `-n`/`--scan-points` points. Same
`shift_II_JJ` folder-naming/manifest scheme either way (one index just
stays `00`) — Stage 2 auto-detects which shape it's looking at from
`sf_manifest.json`, no flag needed.

`example_4.10.sh`'s Case 1 runs a real `5x5 --scan surface` grid; Case 4
runs all 3 of `--scan x/y/xy`.

### 2.2 The 3 interlayer-gap strategies (`--mode`)

| Mode | Strategy | Cost | Accuracy |
|---|---|---|---|
| 1 | SIESTA relaxes z for real: a restricted CG relaxation per point (`x`,`y` frozen for every atom via `%block Geometry.Constraints`, only `z` free, `--relax-z-steps` steps) | highest (a real relaxation per point) | highest (same level of theory as the reported energy) |
| 2 | MACE-MP-0 relaxes z first (same x/y-frozen definition, cheap ML), then a plain SIESTA single-point at the result | middle | middle (mixes a cheap ML pre-relaxation into an otherwise-DFT number) |
| 3 | fixed `--gap` (default 3.2 Å), plain SIESTA single-point, no relaxation at all | lowest | lowest (every point uses the same interlayer distance, right for none of them exactly) |

`example_4.10.sh`'s Case 2 generates the same 2-point sweep 3 times, once
per mode, and shows each mode's own `config_extra.fdf` side by side —
mode 1's `Geometry.Constraints` block, mode 2/3's plain single-point
(`MD.Steps 0`). Mode 2 additionally needs the optional `ml` extra (`pip
install stb_suite[ml]`); if `--d3` is also on (the default), MACE-MP-0's
own D3(BJ) dispersion needs the separate optional `torch-dftd` package
too — the script degrades gracefully (falls back to `--no-d3` for that
one demonstration, or skips mode 2 entirely) if either is missing, same
as it will for you.

**Known caveat, worth reading before trusting a `--mode 1` production
run**: see Section 5.1.

### 2.3 ZSL matching — same-file sliding vs. real heterostructures

Pass the **same file** for `--layer1`/`--layer2` for the canonical case:
a material sliding against itself (e.g. graphite `ABA` vs. `ABC`
stacking). This is a trivial `1×1`, 0-strain match by construction — no
supercell search needed.

Pass **different files** to study an interlayer/heterostructure sliding
landscape instead — `stb-stackingfault` then runs the Zur & McGill
Zero-Strain-Layer (ZSL) algorithm to find the lowest-strain **commensurate
supercell** that fits both lattices, exactly the same matching logic
`stb-2Dstacking` uses. `--max_area`/`--max_strain`/`--match_id` control
the search; `-sm/--strain_mode` (`top`/`bottom`/`sym`) controls how the
mismatch strain is distributed between the two layers; `-t/--twist`
applies a *fixed* twist angle (the same physical system throughout the
whole grid — not swept).

`example_4.10.sh`'s Case 3 runs `graphene.fdf` (a=2.46 Å) against
`hbn.fdf` (a=2.504 Å, ~1.8% mismatched) and finds a real, non-trivial
match: **26.20 Å², 1.76% strain, 20 atoms per grid point** (10 per
layer) — vs. Case 1's trivial 5.24 Å², 0% strain, 4-atom homobilayer
case.

### 2.4 D3 dispersion, ML pre-relax/preview, and `sf_manifest.json`

`--d3` (default **on**) forces the Grimme DFT-D3 dispersion correction in
every generated `config_extra.fdf` — interlayer binding here is van der
Waals-dominated, which plain GGA misses. `--ml-prerelax-layers`
optionally relaxes each monolayer's own positions with MACE-MP-0 before
stacking (a safety net if you're not fully sure your inputs are already
at equilibrium); `--ml-preview` evaluates the whole grid's single-point
energy on MACE-MP-0 first (no SIESTA) and writes a quick PNG heatmap/line
sanity check of the gamma-surface's shape before committing to the real
DFT grid. Both need the optional `ml` extra.

Every run always writes `sf_manifest.json` (grid shape + per-point
shift/gap, machine-readable) under `sf_run/` — the primary source Stage
2 reads, independent of whether `--save-report`'s narrative
`stackingfault_setup.txt` was also requested.

## 3. Stage 2: analyzing the results (`stb-stackingfaultAnalysis`, code `4.10.2`)

### 3.1 Report structure

```
[0] RUN METADATA        -- directory, grid shape, mode
[1] GRID ENERGIES        -- one row per point: E(eV), dE(meV/A2), SCF, MaxF, GapFinal (mode 1 only)
[2] STACKING FAULT ANALYSIS  -- equilibrium/highest-energy registry, corrugation
[3] SUMMARY               -- gnuplot data (if --save-gnuplot), the always-on animation
[4] APPLY                 -- only if --apply was given
```

`--dir` (default `sf_run`, auto-detected as `.` if you're already inside
it) points at Stage 1's own run folder; `--file` (default `calc.out`) is
the SIESTA output filename expected inside every `positions/shift_II_JJ/`.

### 3.2 The meV/Å² scale

`[2]`'s leading `[INFO] Interface area` line reports the in-plane cell
area used for the conversion — always computed fresh from whichever grid
point's `structure.fdf` was read first, `|a × b|` of its own lattice
vectors (Å²). Every relative energy in the report (`[1]`'s `dE(meV/A2)`
column, `[2]`'s equilibrium/highest/corrugation lines) is `(E - E_min) *
1000 / area` — the standard "millielectronvolt per square Ångström"
convention used throughout the 2D-materials literature for this kind of
areal energy density. The absolute per-point `E(eV)` column stays in raw
eV (diagnostic only — an absolute SCF total energy in meV would just be
an unwieldy 6-7 digit number).

### 3.3 `--save-gnuplot` and `plot/`, `--view`, and the always-on animation

Gnuplot `.dat`+`.gplot` output is **opt-in** (`--save-gnuplot`, off by
default), written under `<dir>/plot/` — a `pm3d` heatmap for `--scan
surface`, a line plot for `--scan x/y/xy`. A gnuplot `pm3d` map needs
**every** point of the grid present to render at all (a single missing
point renders the entire map blank, verified live) — an incomplete grid
prints a `[WARNING]` and skips the plot files instead of writing a broken
one; the numeric analysis in `[1]`/`[2]` still works fine on a partial
grid regardless.

`--view` shows the *same* landscape as an interactive matplotlib preview
(a heatmap or line plot, `plt.show()`) — independent of, and more
forgiving than, `--save-gnuplot`'s file output: it renders even on an
incomplete grid (a `NaN` gap in the array, not a blank plot). Following
this suite's own convention (`WORKFLOW_TOOLS` writes gnuplot pairs;
matplotlib here is preview-only), it is **never** written to disk as a
PNG.

An **extended-XYZ animation** (`stackingfault_animation.xyz`, one frame
per analyzed point, `Lattice=`/per-frame `.info` with `label`,
`shift_x`/`shift_y`, `energy_eV`, `dE_meV_per_A2`) is **always** written,
unconditionally — viewable directly in VESTA/OVITO/ASE-GUI, or
interactively via `--view-animation` (ASE's own 3D
viewer, needs a display).

### 3.4 `--apply`

Copies the equilibrium (lowest-energy) point's `structure.fdf` to a path
of your choosing — e.g. as the starting geometry for a follow-up
calculation at that specific registry.

## 4. Worked example: a real 30-point energy landscape

`example_4.10.sh`'s Case 5 fabricates `calc.out` for a real `--scan x -n
30` sweep using **real `siesta: FreeEng` values** from an actual
bilayer-graphene DFT calculation (`PAO.BasisSize DZP`, `XC.Functional
GGA`/`PBE`, real SIESTA 5.4.2 single-points) — not invented numbers,
just written onto this folder's own small, unrelaxed `graphene.fdf` pair
(area 5.2408 Å², vs. the original run's own relaxed lattice) so the exact
meV/Å² figure below is this fixture's own.

```
$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 30 --mode 3 --no-intro
$ stb-stackingfaultAnalysis --dir sf_run --save-gnuplot --view --view-animation --no-intro
```

```
[1] GRID ENERGIES
------------------------------------------------------------
Point           Shift(x,y)          E(eV)           dE(meV/A2)    SCF   MaxF(eV/A)
------------------------------------------------------------------------------------
shift_00_00     (0.000,0.000)       -654.900928     3.188         OK    0.0100
shift_01_00     (0.033,0.000)       -654.901136     3.148         OK    0.0100
shift_02_00     (0.067,0.000)       -654.901574     3.065         OK    0.0100
shift_03_00     (0.100,0.000)       -654.902258     2.934         OK    0.0100
shift_04_00     (0.133,0.000)       -654.903136     2.767         OK    0.0100
shift_05_00     (0.167,0.000)       -654.904736     2.461         OK    0.0100
shift_06_00     (0.200,0.000)       -654.906246     2.173         OK    0.0100
shift_07_00     (0.233,0.000)       -654.907901     1.857         OK    0.0100
shift_08_00     (0.267,0.000)       -654.909639     1.526         OK    0.0100
shift_09_00     (0.300,0.000)       -654.911400     1.190         OK    0.0100
shift_10_00     (0.333,0.000)       -654.913058     0.873         OK    0.0100
shift_11_00     (0.367,0.000)       -654.914659     0.568         OK    0.0100
shift_12_00     (0.400,0.000)       -654.915805     0.349         OK    0.0100
shift_13_00     (0.433,0.000)       -654.916873     0.145         OK    0.0100
shift_14_00     (0.467,0.000)       -654.917482     0.029         OK    0.0100
shift_15_00     (0.500,0.000)       -654.917635     0.000         OK    0.0100
shift_16_00     (0.533,0.000)       -654.917482     0.029         OK    0.0100
shift_17_00     (0.567,0.000)       -654.916873     0.145         OK    0.0100
shift_18_00     (0.600,0.000)       -654.915805     0.349         OK    0.0100
shift_19_00     (0.633,0.000)       -654.914659     0.568         OK    0.0100
shift_20_00     (0.667,0.000)       -654.913058     0.873         OK    0.0100
shift_21_00     (0.700,0.000)       -654.911400     1.190         OK    0.0100
shift_22_00     (0.733,0.000)       -654.909639     1.526         OK    0.0100
shift_23_00     (0.767,0.000)       -654.907901     1.857         OK    0.0100
shift_24_00     (0.800,0.000)       -654.906246     2.173         OK    0.0100
shift_25_00     (0.833,0.000)       -654.904736     2.461         OK    0.0100
shift_26_00     (0.867,0.000)       -654.903136     2.767         OK    0.0100
shift_27_00     (0.900,0.000)       -654.902258     2.934         OK    0.0100
shift_28_00     (0.933,0.000)       -654.901574     3.065         OK    0.0100
shift_29_00     (0.967,0.000)       -654.901136     3.148         OK    0.0100
------------------------------------------------------------------------------------

[2] STACKING FAULT ANALYSIS
[INFO] Interface area (in-plane): 5.2408 Ang^2
Equilibrium stacking (min)    : shift_15_00  (shift 0.5000, 0.0000)  (E = -654.917635 eV, dE = 0.000 meV/Ang^2)
Highest-energy registry (max) : shift_00_00  (shift 0.0000, 0.0000)  (E = -654.900928 eV, dE = 3.188 meV/Ang^2)
Corrugation (stacking-fault) energy : 3.188 meV/Ang^2 (max - min over the sampled grid; 0.016707 eV over a 5.2408 Ang^2 cell)
```

**Reading this physically**: the energy decreases smoothly and
**perfectly symmetrically** from `shift_x=0` (eclipsed/highest energy)
down to `shift_x=0.5` (equilibrium, lowest energy) and back up — exactly
the shape expected for a graphene bilayer sliding along one lattice
direction, with the true minimum-energy `AB`-like registry sitting at
the midpoint of this particular path. No BSSE column anywhere in this
table (Section 1.2) — the raw energy already is the physically correct
comparison.

## 5. Known, deliberate limitations

### 5.1 `--mode 1`'s CG relaxation can converge inconsistently across the grid

Verified live, on a real 30-point `--mode 1` scan: with a loose
`SCF.DM.Tolerance`/`MD.MaxForceTol` (e.g. `0.01 eV/Å`), SIESTA's CG
relaxation can converge in as few as **2-4 steps** — satisfying the force
tolerance, but with very different margins from point to point (some
points land at a residual force of `~0.002 eV/Å`, others at `~0.009
eV/Å`, both technically "converged"). The resulting small interlayer-gap
inconsistency (a few hundredths of an Å) is nearly invisible in the raw
total energy (which varies slowly/quadratically near an already-close
geometry) but can produce a visible kink in `GapFinal(A)` between
neighboring grid points. **Fix**: tighten `MD.MaxForceTol` in your
`calc.fdf` (e.g. `0.001`-`0.005 eV/Å`) if you see a non-smooth
`GapFinal(A)` progression in `[1]`'s table — don't just increase
`--relax-z-steps`, since the step *count* usually isn't the limiting
factor (2-4 steps is typically enough to reach even a tight tolerance on
a system this small).

### 5.2 `--mode 3`'s fixed gap is an approximation

Every point uses the *same* nominal `--gap`, regardless of whether that's
close to the true equilibrium distance at that specific registry (an
eclipsed/high-energy point's true relaxed gap is typically *larger* than
an offset/low-energy point's). `[2]`'s report always prints a `[NOTE]`
reminding you the true, fully-relaxed corrugation may be somewhat lower —
`--mode 1`/`2` exist specifically to remove this approximation, at
increasing cost.

### 5.3 No BSSE correction — by design, not by omission

Covered in full in Section 1.2 — repeated here because it's the most
likely point of confusion for anyone who has used `stb-adsorbBsse` or a
cohesive-energy workflow and expects an analogous stage here. There isn't
one, and there shouldn't be: the raw energy's corrugation is already the
exact, correct physical answer for this specific kind of same-atom-count
comparison.

## 6. Step-by-step: running this workflow on your own structure(s)

1. Prepare 1 or 2 monolayer `.fdf` files (same file for a material
   sliding against itself; 2 different files for a heterostructure) and a
   `calc.fdf` template (basis, k-grid, XC — Stage 1 forces the
   single-point/z-relaxation behavior itself via `config_extra.fdf`,
   independent of whatever `MD.TypeOfRun`/`MD.Steps` your template has).
2. `stb-stackingfault -l1 <layer1.fdf> -l2 <layer2.fdf> -c calc.fdf
   --mode 1` (most accurate — start here unless you have a specific
   reason to trade accuracy for speed) — pick `-nx`/`-ny` (surface) or
   `--scan x/y/xy -n` (1D) depending on whether you need the full 2D map
   or a cheaper profile first.
3. Run SIESTA in every `sf_run/positions/shift_II_JJ/` folder.
4. `stb-stackingfaultAnalysis --dir sf_run --save-report --save-gnuplot
   --view` — inspect `[1]`'s `GapFinal(A)` column for the smoothness
   issue in Section 5.1 before trusting the corrugation number if you
   used `--mode 1`.

## 7. Files in this folder

```
graphene.fdf        Primitive 2-atom graphene monolayer (a=2.46 Ang)
hbn.fdf              Primitive 2-atom h-BN monolayer (a=2.504 Ang, ~1.8% mismatched)
calc.fdf             Shared calc.fdf template (DZP, GGA-PBE) -- Stage 1 overrides its
                      MD block per-point regardless (Section 2.2)
example_4.10.sh       This walkthrough
README.md             This file
.gitignore            Ignores output/ (regenerated by the script every run)
```

## 8. Running the script

```
cd examples/4.10-stacking_fault
./example_4.10.sh
```

Needs no real SIESTA binary for Stage 1 (writes `.fdf` files only, plus a
real MACE-MP-0 z-relaxation for `--mode 2` if the optional `ml` extra is
installed). Stage 2's worked example (Case 5) uses fabricated `calc.out`
data (Section 4) since this walkthrough doesn't invoke real SIESTA
either. Every case writes into its own `output/case*/` subfolder;
`output/` itself is gitignored and safe to delete/regenerate any time.

## 9. What's next

Once you have a real gamma-surface for your own material:
- The equilibrium `structure.fdf` (via `--apply`) is a natural starting
  point for a follow-up calculation at that specific registry — an NEB
  path between two symmetry-equivalent minima (workflow `4.9`), or a
  phonon calculation (workflow `4.4`) to check its dynamical stability.
- The corrugation energy is directly comparable to literature values for
  interlayer sliding/friction barriers in the same material family.
