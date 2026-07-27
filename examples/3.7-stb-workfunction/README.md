# 3.7 — Work Function Calculator (`stb-workfunction`)

## What this tool does

`stb-workfunction` computes a material's **work function** — the minimum
energy needed to remove an electron from the material and place it at rest
just outside its surface, in vacuum — straight from a finished SIESTA
calculation, with no extra SIESTA run needed:

- Reads `<label>.VT` (the total electrostatic potential, on SIESTA's
  real-space grid) and planar-averages it along the direction normal to
  the surface (`--axis`).
- Reads the Fermi energy from `<label>.out` (or `--fermi` to override).
- Automatically finds the flat **vacuum plateau(s)** in that averaged
  profile and computes `Φ = E_vacuum - E_F` for each one.
- Handles an **asymmetric slab** (different terminations on each side)
  correctly: two genuinely different vacuum levels are reported
  separately, not averaged into one physically meaningless number.
- Auto-detects the vacuum-padded axis from the structure's own geometry
  (`<label>.XV`/`.fdf`) when available.

## Why this matters (a bit of theory)

### The work function itself

The work function `Φ` is one of the most basic, experimentally
measurable properties of a surface — it sets the energy barrier for
electron emission (thermionic/photoemission), governs Schottky-barrier
heights at metal-semiconductor junctions, and determines band alignment
at heterointerfaces. In DFT, it's computed from two ingredients you
already have after any slab calculation:

```
Φ = E_vacuum - E_F
```

`E_F` is the Fermi energy (the "top of the electron sea" inside the
material); `E_vacuum` is the electrostatic potential energy an electron
would have sitting at rest far outside the surface, in the vacuum region
of the simulation cell.

### Why this only makes sense for a slab, not a 3D bulk crystal

A genuinely 3D-periodic crystal has **no vacuum** anywhere in its cell —
`E_vacuum` simply isn't defined. `stb-workfunction` needs a structure with
a real vacuum gap along one axis (a slab, a 2D monolayer, a wire) so that
axis's planar-averaged potential has a region *flat enough* to call
"vacuum." This is exactly why `--axis` auto-detection reuses
`core.kspace.detect_vacuum_axes` — the same vacuum-gap detector used by
`stb-kgrid`/`stb-mlrelax` elsewhere in this suite — and why the tool warns
explicitly if you point it at an axis that doesn't actually look
vacuum-padded.

### The planar average, and why it "flattens" in vacuum

`V(x,y,z)` on the raw 3D grid is full of atomic-scale wiggles (it peaks
sharply near every nucleus). Averaging it over the two in-plane
directions at each point along the surface normal — `V(z) = <V(x,y,z)>_xy`
— washes out those in-plane wiggles while keeping the physically
meaningful variation *along* the surface normal. Deep in the vacuum, far
from any atom, there's nothing left to average over that varies at all —
the profile goes flat, and that flat value **is** `E_vacuum`.
`stb-workfunction` finds this flat region automatically (the lowest
-gradient contiguous run of the averaged profile — see
`find_vacuum_plateaus`), rather than asking you to eyeball a plot and
pick a point.

**A concrete distinction worth knowing**: the reported "vacuum size" is
the width of this genuinely-flat region only — not the whole geometric
vacuum gap in your cell. Near the slab-vacuum interface the potential is
still transitioning (not yet flat), so the flat plateau is always
*narrower* than `cell length - slab thickness`. In the CrS example below,
the real vacuum gap is about 20 Å (a 23.35 Å cell minus a ~3.1 Å-thick
layer), but the genuinely flat plateau `stb-workfunction` finds and
averages over is only ~3.4 Å wide — both numbers are correct, they just
answer different questions ("how much empty space is there" vs. "how
much of it is flat enough to trust as *the* vacuum level").

### Asymmetric slabs and the surface dipole

If a slab's two surfaces are different (different terminations, an
adsorbate on only one side, a deliberately asymmetric cut), the two
vacuum regions genuinely sit at **different** potential levels — there
are two distinct work functions, and averaging them would throw away
real physics. The difference between the two vacuum levels, `dV`, is
directly proportional to the slab's own surface dipole moment (`dV ~
4π·μ/A`, in Gaussian-like units). `stb-workfunction` reports `dV` itself
rather than converting it to an absolute dipole moment in `e·Å`, since
that conversion needs a convention-specific prefactor (see e.g.
Bengtsson, *Phys. Rev. B* **59**, 12301 (1999) for the full derivation).

### A systematic slope = a missing dipole correction

A genuine vacuum plateau should be **flat**, not sloped. A slow,
systematic ramp across an otherwise-flat vacuum region (as opposed to
just noisy scatter) is the classic numerical symptom of an uncancelled
periodic-image electric field on an asymmetric slab — SIESTA's own fix
for this is `SlabDipoleCorrection`. `stb-workfunction` fits a line to
each plateau and flags this automatically (`detect_plateau_slope`),
without needing to parse whether the flag was actually set in your `.fdf`
— it's detected from the physical symptom itself, so it works regardless
of SIESTA version or how exactly that setting is logged.

## Limitations

- **The chosen axis must genuinely have vacuum.** If none of the 3 axes
  has a real vacuum gap (a 3D bulk crystal), there's nothing to compute —
  the tool fails cleanly rather than returning a meaningless number (see
  the "bulk, no vacuum" case in `test/3-analysis/7-workfunction/test.sh`).
- **A "noisy" plateau is a warning, not an automatic correction.** If the
  detected vacuum region's standard deviation is above a small tolerance,
  the reported vacuum level/work function should be treated as
  approximate.
- **`dV` is not converted to an absolute dipole moment** (see above) —
  it needs a convention-specific prefactor this tool doesn't assume for
  you.
- **An atom projecting into the "vacuum" plateau may mean it isn't real
  vacuum** — e.g. an isolated adsorbate on one side of an otherwise
  symmetric slab. `stb-workfunction` flags this rather than silently
  trusting the plateau.
- **The geometry cross-check can't catch every mismatch.** Comparing the
  `.XV`/`.fdf` cell against the `.VT` grid's own cell catches a
  genuinely different calculation being paired by accident, but a
  same-cell, different-atomic-position mismatch isn't detectable from
  these files alone.

## The report: console output, `--save-report`, `--save-gnuplot`, `--view`

Every run prints a numbered report to the console:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | label, Fermi/grid file, axis (requested or auto-detected), tolerances, output dir, active options |
| `[1] FERMI ENERGY` | value and its source (detected `.out` vs. `--fermi` override) |
| `[2] VACUUM AXIS DETECTION` | axis used, and the geometry/grid cross-check result |
| `[3] VACUUM PLATEAU DETECTION` | table of detected plateaus (position, mean potential, noise), with warnings (atoms projecting in, systematic slope) |
| `[4] WORK FUNCTION RESULTS` | vacuum size/slab thickness, Φ per plateau, average + `dV` if asymmetric |
| `[5] OUTPUT DATA & PLOTS` | whether `workfunction_data.dat`/`.gplot` were written (only with `--save-gnuplot`) |
| `[6] REFERENCES` | writes `references.bib` (SIESTA) |
| `[7] SUMMARY & FILES` | status and a recap of every file written |

- **`--save-report`** — also persists the full numbered report to
  `stb_workfunction_report.txt`. Off by default.
- **`--save-gnuplot`** — also writes `workfunction_data.dat` (the
  planar-averaged potential profile, with the detected vacuum level(s)/
  work function(s) noted in its header) and `workfunction.gplot` (an
  annotated gnuplot script — Fermi level and each vacuum level drawn as
  dashed reference lines, with the work function itself marked between
  them). Off by default — this tool used to write both **unconditionally**
  on every run; that's gone, they're opt-in now, same convention as
  `--save-gnuplot` on `stb-dos`/`stb-bader`/etc.
- **`--view`** — shows an interactive matplotlib preview of the same
  profile before exiting. Off by default (this **replaces** the old
  `--no-plot` flag, which showed a plot by default and only let you turn
  it *off* — inverted here to match every other tool in this suite:
  opt-in, not opt-out).

## When you'd reach for it

- Right after any slab/surface SIESTA calculation, to get the work
  function without hand-plotting the `.VT` file yourself.
- Comparing work functions across a series of surfaces/terminations/
  adsorbate coverages.
- Checking for a missing `SlabDipoleCorrection` on an asymmetric slab
  (the systematic-slope warning) before trusting any surface-related
  result from that calculation.

## Two ways to run it

A — direct CLI:
```bash
stb-workfunction --label siesta --file calc.out
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.7
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.7.sh` proves this directly at the end.

## Files in this folder

- `calc.out`, `siesta.VT`, `siesta.XV` — a real, finished, converged
  SIESTA calculation (copied from `test/3-analysis/7-workfunction/`) on a
  CrS monolayer (SystemLabel `siesta`) — a genuinely 2D material fetched
  from the `twodmatpedia` OPTIMADE database via `stb-fetch`, with a large
  vacuum gap along `c`, making it a real, physically meaningful
  work-function case (not a synthetic slab). `--file calc.out` is needed
  because this run's own log file isn't named `siesta.out`.
- `example_3.7.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/7-workfunction/test.sh` for that, which uses small
  synthetic grids instead, to exercise every code path — asymmetric
  slabs, sloped vacuum, mismatched geometry, etc. — without committing
  large binary fixtures).
- `.gitignore` — excludes `output/` and other generated files.

## Running the walkthrough

```bash
cd examples/3.7-stb-workfunction
./example_3.7.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder            | What it shows                                                              |
|--------------------|-----------------------------------------------------------------------------|
| `basic/`           | The full report on a real, converged DFT calculation: real Fermi energy, real work function |
| `save-gnuplot/`    | `--save-gnuplot`: `workfunction_data.dat` + the annotated `workfunction.gplot` |
| `full-report/`     | Default (no report file) vs. `--save-report` (`stb_workfunction_report.txt`), `references.bib` |

## Try it yourself

```bash
# Right after a slab calculation
stb-workfunction --label my_slab --save-report

# Also get the annotated gnuplot script for a publication-ready plot
stb-workfunction --label my_slab --save-gnuplot

# Force the Fermi energy if your .out file isn't the default <label>.out
stb-workfunction --label my_slab --file relax/my_slab.out
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-l/--label`        | SystemLabel used in SIESTA.                                            |
| `-f/--file`         | `.out` file to read the Fermi energy from (default `<label>.out`).     |
| `-g/--grid`         | Potential grid file (default `<label>.VT`).                            |
| `-z/--axis`         | Axis normal to the surface, `0`/`1`/`2` (default: auto-detected).      |
| `--fermi`           | Manually force the Fermi energy (eV), instead of reading `.out`.       |
| `--asymmetric-tol`  | Vacuum levels within this many eV are treated as one plateau (default `0.05`). |
| `-o/--output-dir`   | Where `references.bib` (and the report/data+gnuplot, with their flags) land. |
| `--save-report`     | Persist the full numbered report to `stb_workfunction_report.txt`.     |
| `--save-gnuplot`    | Also write `workfunction_data.dat`/`workfunction.gplot`.               |
| `--view`            | Show an interactive matplotlib preview of the potential profile.       |

## What's next

`stb-density` (3.8) reads the same class of SIESTA real-space grid files
(`.RHO` instead of `.VT`) for charge-density maps/clouds rather than a
single planar-averaged profile; `stb-bader` (3.6) partitions the charge
density into atomic basins instead of averaging the electrostatic
potential. All three read a finished SCF calculation's own real-space
output, no re-run needed.
