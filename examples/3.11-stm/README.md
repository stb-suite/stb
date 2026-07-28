# 3.11 — STM Simulator (`stb-stm`)

## What this tool does

`stb-stm` turns a SIESTA energy-integrated local density of states grid
(`<label>.LDOS`) into a simulated STM image, under the **Tersoff-Hamann
approximation**: tunneling current is proportional to the LDOS at the tip
position. Two modes:

- `--mode current` (**default**) — the classic *constant-current* STM
  image: for every `(x, y)`, finds the height where the LDOS first drops
  to `--iso` coming in from vacuum — the topography a real tip follows to
  hold tunneling current constant.
- `--mode height` — the simpler *constant-height* image: a flat 2D LDOS
  map at a fixed height `--z` above the surface.

## Why this matters (a bit of theory)

### The Tersoff-Hamann approximation

A real STM tip is a complicated, atomically-sharp, chemically-specific
object — modeling it exactly is a hard many-body tunneling problem. Tersoff
and Hamann (1985) showed that if you approximate the tip as a structureless
**s-wave point probe**, Bardeen's tunneling theory collapses to a strikingly
simple result:

```
I(x, y, z) is proportional to LDOS(x, y, z; E_F)
```

the tunneling current at tip position `(x, y, z)` is proportional to the
**sample's own local density of states, evaluated at the tip position**,
integrated over the electrons that can actually tunnel at the applied bias
(see below). This is *why* an STM image can be simulated directly from a
DFT calculation's electron density alone, with no separate tip calculation
at all — and why it is only a **proportionality**, not a calibrated,
absolute tunneling current in Amperes: `--iso` (or the LDOS value under
the tip in `--mode height`) is in the grid's own units (e/Bohr^3), tuned
for sensible image contrast, not read off a real experiment's setpoint
current.

### What a `.LDOS` file actually holds, and why the bias WINDOW matters

`%block LocalDensityOfStates Emin Emax` tells SIESTA to integrate the
local density of states over an energy window `[Emin, Emax]` and write the
result as a real-space grid — the same binary format as a charge-density
`.RHO` file, but holding

```
LDOS(r) = integral from Emin to Emax of  rho(r, E) dE
```

instead of the density integrated over *all* occupied states. **This
bias window is a property of the SIESTA run itself, baked into the grid
file — `stb-stm` cannot see or change it after the fact.** The energies
can be given as absolute values, or relative to the Fermi level using the
`EF` keyword as the first token — this is the setting used to generate
this example's own `siesta.LDOS`:

```
%block LocalDensityOfStates
EF -3.50 0.00  eV
%endblock LocalDensityOfStates
```

`EF -3.50 0.00 eV` integrates from 3.50 eV *below* the Fermi level up to
the Fermi level itself — i.e., **occupied states only**. This is the STM
analog of a negative sample bias in a real experiment (electrons tunnel
*from the sample to the tip*): the image shows where the *filled* states
live. A positive window (e.g. `EF 0.00 3.50 eV`) would instead image
**empty** states, and in general a different bias window can show visibly
different orbital character for the same material (a real, well-known STM
effect, not a simulation artifact) — always know which window you (or
whoever ran the calculation) actually used before interpreting an image.

**Note on the exact fdf syntax**: verified live while building this
example that SIESTA's parser requires `EF Emin Emax` and the unit
(`eV`) on the **same line** — splitting them across two lines (`EF -3.50
0.00` then `eV` on the next line) raises a hard `ERROR in
LocalDensityOfStates block!` and aborts the run. Always write it as one
line, exactly as shown above.

### `--mode height` vs. `--mode current`: two different physical pictures

- `--mode height` is the cheaper, simpler case: just read off the LDOS
  grid at one fixed height everywhere. It is NOT what a real constant
  -current STM does, and can be misleading for a corrugated surface — the
  tip would physically crash into a protruding atom at a height that's
  perfectly safe elsewhere. It also has no built-in safeguard against a
  poorly chosen `--z`: verified live in this walkthrough that the default
  `3.0` Ang gives an entirely flat, near-zero image for this example's own
  occupied-states LDOS (already decayed by that height), while `--mode
  current`'s feedback loop always finds a height with real signal (or
  clearly reports when it can't) on the exact same grid.
- `--mode current` reproduces the real feedback loop: at every `(x, y)`,
  the tip retreats from far above the surface until the LDOS first reaches
  `--iso` — exactly mimicking a real STM's constant-current topography
  scan, and the mode that produces the textbook "STM image" with visible
  atomic corrugation. `stb-stm` searches **outside-in** (farthest height
  first) so it finds the same point a real retreating tip would stop at,
  not an arbitrary crossing.

### Why `--iso`/`--z` now have defaults, and what they mean physically

Both used to be required with no default. `--iso` defaults to `0.001`
e/Bohr^3, `--z` to `3.0` Ang — reasonable starting points for a typical
slab, not a universally "correct" value for every material/bias window;
tune them if the default image looks saturated (too many points crash
into the surface or never reach `--iso`) or empty (nothing does — see
`iso-sensitivity/` below for exactly this on a real structure).

### Locating the surface: a real bug found using this exact fixture

Every height reported by this tool is **relative to the topmost atom**
along the detected (or `--axis`-forced) surface-normal direction — never
the raw simulation cell origin. This example's own `siesta.XV` (a real
CrS monolayer, fetched from the twodmatpedia database — the same
structure already used in `examples/3.7-stb-workfunction/` and
`examples/3.9-stb-xrd/`) is exactly what caught a real, verified bug while
building this walkthrough: its 4 atoms sit at fractional z = `0, 0,
0.066, 0.934` (2 Cr at the cell's stored origin, 2 S buckled slightly to
either side) — a naive "topmost atom = `max(z)`" picks the `z=0.934` S
atom and searches **upward** from there, straight into the tiny ~7% sliver
before wrapping back to the Cr atoms — collapsing the whole search window
to ~1.5 Ang instead of the genuine ~20 Ang vacuum region that actually
separates this thin monolayer's two faces. Fixed with a new
`core/kspace.py::find_surface_reference` helper: it reuses the same
gap-finding logic that already decides *whether* an axis is vacuum-padded
(`detect_vacuum_axes`), now also used to find *where* the real vacuum
starts, correctly identifying the `z=0.066` S atom as the reference
instead.

A second, related subtlety this same fixture exposed: a periodic cell's
stacking axis is topologically a ring, so **any** single compact atomic
region surrounded by vacuum has TWO faces exposed to that same gap (a
"top" and a "bottom"), not one. Searching the *entire* ~20 Ang gap gave
wildly unphysical "corrugation" values (18-20 Ang — real STM corrugation
is sub-Angstrom to a few Ang) at several `--iso` values tried live,
because the outside-in scan eventually crossed into the *far* face's own
LDOS tail. `find_surface_reference` caps the search at **half** the
identified gap — matching this tool's own pre-existing documented
limitation ("only images the surface exposed in the +axis direction...
a slab with two exposed faces only gets its 'top' one imaged this way"),
just now enforced numerically instead of relying on the atom arrangement
to accidentally already imply it (as a perfectly centered structure, like
this suite's own graphene test fixture, happens to). Verified live:
identical numbers as before on every existing (conventional) fixture —
this fix only changes anything for a structure like this one — and a
physically sensible ~1.0-1.7 Ang corrugation on this real CrS monolayer
once fixed.

`stb-stm` refuses to run rather than guess if it can't find exactly one
vacuum-padded axis at all (unlike `stb-workfunction`, which can degrade to
an assumed axis — an STM image without a genuine, confirmed surface is
meaningless, not just less precise).

## The report: console output, `--save-report`, `--save-gnuplot`, `--view`

Every run prints a numbered report, the same `[0]...[5]` style every newer
tool in this suite uses:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | input files, mode, `--iso`/`--z` (whichever applies), output dir, active options |
| `[1] INPUT DATA` | grid shape, LDOS range, surface-normal axis, vacuum-padded axes, topmost atom height, axis length |
| `[2] STM IMAGE` | mode-specific results: actual height used + LDOS stats (`height`), or iso/coverage/height-map stats + corrugation (`current`) |
| `[3] OUTPUT DATA & PLOTS` | whether `stm_<mode>.dat`/`.gplot` were written |
| `[4] REFERENCES` | writes `references.bib` (SIESTA + Tersoff-Hamann) |
| `[5] SUMMARY & FILES` | status and a recap of every file written |

- **`--save-report`** — persist the report to `stb_stm_report.txt`. Off by
  default.
- **`--save-gnuplot`** — write `stm_<mode>.dat` (a pm3d-blocked X/Y/value
  file) + `stm_<mode>.gplot` together. Off by default — this tool used to
  write both **unconditionally** on every run. The `.gplot` script's own
  `set output`/`splot` filenames are always bare basenames (`stm_current.pdf`,
  never `some_dir/stm_current.pdf`) — a real, verified bug found while
  rewriting this tool: the old `set output` line kept the `--output-dir`
  prefix (inconsistently with `splot`, which already stripped it), so
  running `gnuplot stm_current.gplot` from *inside* the output directory
  itself (the intended, documented usage) tried to write to a nested,
  usually-nonexistent `output_dir/output_dir/...` path.
- **`--view`** — show an interactive matplotlib preview (a filled 2D color
  map, mirroring the `.gplot`'s own `pm3d map`). Off by default — this
  tool previously had **no** matplotlib option at all.

## When you'd reach for it

- Simulating what a real STM experiment would see on a finished SIESTA
  slab calculation, before (or instead of) a costly real measurement.
- Comparing occupied- vs. empty-states images (different `%block
  LocalDensityOfStates` bias windows) to understand bias-dependent
  contrast on a real material.
- Getting a quick, physically grounded corrugation estimate (`--mode
  current`) to compare against a real STM measurement's own topography.

## Two ways to run it

A — direct CLI:
```bash
stb-stm --file siesta.LDOS --geometry-file siesta.XV --mode current
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.11
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.11.sh` proves this directly at the end. The menu
defaults `--iso`/`--z` the same way the CLI does (just press Enter), plus
separate `y/N` prompts for `--save-report`/`--save-gnuplot`/`--view`.

## Files in this folder

- `siesta.LDOS` / `siesta.XV` / `calc.fdf` / `structure.fdf` — a real
  SIESTA calculation of a CrS monolayer (fetched via `stb-fetch` from the
  twodmatpedia OPTIMADE database, id `2dm-2617` — the same structure
  already used in `examples/3.7-stb-workfunction/`/`examples/3.9-stb-xrd/`),
  with `%block LocalDensityOfStates EF -3.50 0.00 eV` (occupied states
  only, see the theory section above) producing `siesta.LDOS`. This exact
  fixture is what caught the real vacuum-side bug described above.
- `example_3.11.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/11-stm/test.sh` for that, which uses its own
  separate graphene fixture with a conventional, centered vacuum axis).
- `.gitignore` — excludes `output/` and every file this tool itself
  generates on a run.

## Running the walkthrough

```bash
cd examples/3.11-stm
./example_3.11.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder                  | What it shows                                                              |
|--------------------------|-----------------------------------------------------------------------------|
| `constant-current/`      | The default mode with default `--iso`: the real STM feedback-loop image, corrugation reported |
| `constant-height/`       | `--mode height` at the default `--z` going flat (LDOS already decayed to ~0 there) vs. a closer `--z` recovering real contrast |
| `iso-sensitivity/`       | The same structure at three `--iso` values -- too high (mostly NaN), sensible (default), and very low |
| `full-report/`           | Default (no report/data files) vs. `--save-report --save-gnuplot`, `references.bib` |

## Try it yourself

```bash
# Your own finished SIESTA slab calculation (needs %block LocalDensityOfStates -- see above)
stb-stm --label my_slab --mode current --save-report --save-gnuplot

# A shallower/deeper occupied-states window, or an empty-states one, needs a NEW SIESTA run
# (the bias window is baked into the .LDOS file itself, not a stb-stm option)

# Constant-height at a custom tip distance, with a live preview
stb-stm --label my_slab --mode height --z 4.0 --view
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--label`           | SIESTA label; shorthand for auto-detected `.LDOS`/`.STM.LDOS` + `.XV`/`.fdf`. |
| `--file`/`--geometry-file` | Explicit `.LDOS`/geometry paths (alternative to `--label`).      |
| `--axis`            | Surface-normal axis (0/1/2). Auto-detected from vacuum padding if omitted. |
| `--mode`            | `current` (default) or `height`.                                       |
| `--z`               | Height above the topmost atom for `--mode height` (default `3.0` Ang). |
| `--iso`             | LDOS threshold for `--mode current` (default `0.001` e/Bohr^3).       |
| `--z-max`           | Upper bound of the search/plot window (default: half the detected vacuum gap). |
| `-o/--output-dir`   | Where all generated files (and `references.bib`) land.                 |
| `--save-report`     | Persist the full report to `stb_stm_report.txt`.                       |
| `--save-gnuplot`    | Also write `stm_<mode>.dat`/`stm_<mode>.gplot`.                        |
| `--view`            | Show an interactive matplotlib preview.                                |

## What's next

`core/kspace.py`'s vacuum-axis detection (and now, the surface-reference
gap-finding this tool's own fix added) is shared with
`stb-workfunction`/`stb-mlrelax`/`stb-kgrid` and every other tool that
needs to tell a genuinely periodic axis apart from a vacuum-padded one.
`stb-density` covers the same real-space-grid family (`.RHO`, total
charge/spin density) with a richer set of cut modes (slice/3D/profile)
for a quantity that isn't tied to a surface the way an STM image is.
