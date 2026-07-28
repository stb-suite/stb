# 3.9 — Powder XRD Simulator (`stb-xrd`)

## What this tool does

`stb-xrd` simulates a powder X-ray diffraction pattern (peak list +
optional plot) from a structure, using pyxtal's diffraction engine. It
reports:

- **Peak positions and intensities** — every reflection (`hkl`) that falls
  inside a scanned 2-theta range, sorted by intensity — saved to
  `xrd_pattern.dat` (`--save-gnuplot`), not printed to the console/report.
- **Structure/symmetry info** — reduced formula, space group, crystal
  system, point group, lattice parameters, cell volume, density.
- Optionally, a **similarity score** against an experimental pattern
  (`--compare-to`), and an interactive matplotlib preview (`--view`).

## Why this matters (a bit of theory)

### Bragg's law

A crystal is a periodic stack of atomic planes. X-rays reflecting off
successive parallel planes interfere constructively only at specific
angles — Bragg's law:

```
n * lambda = 2 * d_hkl * sin(theta)
```

where `lambda` is the X-ray wavelength, `d_hkl` is the spacing between
`(hkl)` planes, and `theta` is half the scattering angle (`2*theta` is
what's actually measured). A powder sample contains crystallites in every
possible orientation, so every set of planes that satisfies Bragg's law
for *some* orientation contributes a peak at its own `2*theta` — the
"powder pattern" is the full peak table as a function of angle.

### The structure factor and systematic absences

Not every geometrically allowed reflection actually shows up: the
**structure factor** `F_hkl` (a sum over every atom's scattering
contribution, weighted by a phase factor depending on its position) sets
each peak's intensity, and a crystal's symmetry can force `F_hkl = 0`
*exactly* for certain `(hkl)` combinations — "systematic absences". This
is why a peak table is a fingerprint of a structure's symmetry, not just
its cell size: two structures with the same lattice parameters but
different space groups can have visibly different missing/present peaks.

### Wavelength choice

Different X-ray sources (Cu-Kalpha, Mo-Kalpha, ...) have different
wavelengths, which rescales the whole pattern via Bragg's law — a shorter
wavelength (e.g. Mo-Kalpha, ~0.71 Ang) reaches higher resolution (smaller
`d_hkl`) within the same 2-theta range than a longer one (e.g. Cu-Kalpha,
~1.54 Ang). `--wavelength` accepts a named source or a plain number in Ang.

### The similarity score (`--compare-to`)

Given an experimental pattern (2-theta, intensity pairs), `--compare-to`
interpolates both the simulated and experimental patterns onto a common
grid and reports pyxtal's cosine-weighted similarity (0-1, higher is more
similar) — a quick, single-number check of "does this candidate structure
plausibly match the measured pattern", useful e.g. for narrowing down
candidate phases before a full Rietveld refinement.

**Required column format for the `--compare-to` file**: 2-theta (degrees)
in the **first** column, intensity in the **last** column — whitespace- or
comma-separated, blank lines and `#` comments skipped. Any columns in
between are ignored, so this accepts both a plain 2-column file AND
stb-xrd's own `--save-gnuplot` output (`xrd_pattern.dat`, 6 columns:
`2theta d h k l intensity`) directly, without needing to strip the
`d`/`h`/`k`/`l` columns out yourself first.

**Peak lists are auto-broadened before comparing.** A real, continuously
-scanned diffractogram (many closely and roughly evenly spaced points) is
compared as-is. But a *discrete peak list* — e.g. a literature-reported
indexed peak table, or stb-xrd's own stick-pattern `xrd_pattern.dat` — has
no peak width at all, while the simulated pattern `--compare-to` compares
it against is always a Gaussian-broadened continuous profile (FWHM=0.1
deg). Comparing the two representations directly is apples-to-oranges:
pyxtal's cubic interpolation draws a spurious curve through the empty gaps
between sparse peaks, which does NOT match the broadened profile's true
near-zero baseline there — even for the exact same structure. stb-xrd
detects this (a peak list has a highly uneven spacing between consecutive
2-theta points — measured via the coefficient of variation of that
spacing — vs. a real scan's near-constant step) and automatically
Gaussian-broadens a detected peak list the same way before computing the
similarity. `--raw-experimental` disables this and compares the file
exactly as read, for the rare case you want that instead.

**Real bug found and fixed (verified live)**: feeding stb-xrd's own
`xrd_pattern.dat` for a structure back into `--compare-to` **for that same
structure** used to score as low as 0.32-0.46 similarity instead of ~1.0,
from two compounding issues: (1) `--compare-to` used to always read
intensity from column index 1, which in a 6-column `xrd_pattern.dat` is
`d` (the d-spacing in Ang), not intensity; and (2) even after fixing that,
comparing the raw, un-broadened stick pattern against the broadened
simulated profile still only scored ~0.32, for the reason explained above.
Fixed by reading intensity from the last column always, and by
auto-broadening peak-list-shaped input before comparing — verified live
that the exact same self-comparison now scores 1.0000, that a genuinely
different structure's pattern still scores well below 1.0 (so the fix
doesn't just make everything score 1.0), and that a genuine dense/uniform
continuous scan is correctly left un-broadened.

## Limitations

- **This is an ideal, simulated powder pattern** — no peak broadening, no
  instrumental profile, no preferred orientation, no temperature (Debye-
  Waller) factors. Good for phase identification/comparing candidate
  structures, not for Rietveld-quality quantitative refinement.
- **The similarity score needs >= 4 experimental points** (pyxtal's cubic
  interpolation requirement) — fewer raises a clear error instead of a
  cryptic one from deep inside pyxtal.
- **A vacuum-padded structure (a slab/2D monolayer) is still simulated as
  fully 3D-periodic** — this tool has no concept of "this axis is just
  vacuum". A huge artificial lattice parameter along the vacuum direction
  produces spurious LOW-ANGLE `(00l)` peaks with very short real-space
  periodicity that don't correspond to anything a real bulk powder sample
  would show (there's no actual repeating unit at that spacing — it's the
  simulation cell's vacuum gap, not a physical lattice plane). The
  walkthrough below hits this directly: CrS's own real geometry has a huge
  vacuum gap along `c`, and the strongest simulated peak is exactly this
  kind of `(001)` low-angle artifact. Worth knowing before reading too much
  into a low-angle peak from a slab/monolayer structure.

## The report: console output, `--save-report`, `--save-gnuplot`, `--view`

Every run prints a numbered report to the console:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | input file/format, wavelength, 2-theta range, `--top` setting, output dir, active options |
| `[1] STRUCTURE` | formula, site count, space group/crystal system/point group/Hall symbol/layer group, lattice parameters, cell volume, density |
| `[2] DIFFRACTION PATTERN` | peak count, strongest peak, resolution (min d-spacing) — a compact summary only, no peak-by-peak table |
| `[3] EXPERIMENTAL COMPARISON` | (conditional — `--compare-to`) file used, point count, whether a peak list was detected/auto-broadened, similarity score |
| `[4] OUTPUT DATA & PLOTS` | whether `xrd_pattern.dat`/`.gplot` were written |
| `[5] REFERENCES` | writes `references.bib` (SIESTA + pyxtal) |
| `[6] SUMMARY & FILES` | status and a recap of every file written |

**The full peak-by-peak table is not printed anywhere** (console or
`--save-report`'s file) — it used to be, in full, which could run to
hundreds of lines for a low-symmetry structure over a wide 2-theta range.
It only ever lives in `xrd_pattern.dat`.

- **`--save-report`** — persist the (table-free) report to `stb_xrd_report.txt`.
- **`--save-gnuplot`** — write `xrd_pattern.dat` (the actual peak list —
  2theta, d, h, k, l, intensity — with a complete, human-readable header:
  structure/formula/space group/wavelength/range, not just pyxtal's own
  terse one-liner) + `xrd_pattern.gplot` (a stick-pattern plot script)
  together — off by default (this tool used to write the data file
  **unconditionally** on every run; that's no longer the case). `--top`
  controls how many peaks land in this file (it used to just trim the
  console table, which no longer exists).
- **`--view`** — show an interactive matplotlib preview of the simulated
  stick pattern (or, with `--compare-to`, an overlay of both patterns).
  Renamed from the old `--plot`, same off-by-default behavior — but fixed:
  the old code delegated to a non-blocking call (`fig.show()`, or pyxtal's
  own `plot_pxrd()`, which was checked directly and also only ever calls
  the non-blocking `fig.show()` internally), so the window would disappear
  as soon as the script exited. `--view` now builds its own plot and calls
  the blocking `plt.show()`, the same convention already used correctly by
  `stb-workfunction`/`stb-density`'s own `--view`.

## When you'd reach for it

- A quick, cheap check of what a candidate structure's powder pattern
  should look like, before running a real diffraction experiment.
- Narrowing down candidate phases against a measured pattern
  (`--compare-to`).
- Sanity-checking a structure's symmetry: does the saved peak list show
  the systematic absences you'd expect from its space group?

## Two ways to run it

A — direct CLI:
```bash
stb-xrd --file structure.fdf --format fdf
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.9
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.9.sh` proves this directly at the end.

## Files in this folder

- `siesta.XV` — a real, converged SIESTA calculation's geometry (copied
  from `examples/3.7-stb-workfunction/`) on a CrS 2D monolayer (fetched
  from the twodmatpedia OPTIMADE database via `stb-fetch`).
- `mock_experimental.dat` — a **synthetic** "experimental" pattern (small
  angle jitter + intensity noise added to this structure's own simulated
  peaks), used to exercise `--compare-to` — not independent lab data.
- `example_3.9.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/9-xrd/test.sh` for that). Converts `siesta.XV` to
  `crs_structure.fdf` itself, up front (see the note below).
- `.gitignore` — excludes `crs_structure.fdf` (regenerated by the script),
  generated `.dat`/`.gplot`/`.pdf` files, `references.bib`, the report,
  and `output/`.

**Why the script converts `siesta.XV` via `sisl` directly, not
`stb-translate`**: `stb-xrd` only accepts `.fdf`/`.STRUCT_OUT` structure
files. While preparing this example, `stb-translate`'s own `siesta` input
reader (`getatomsandvectors_siesta` in `translate.py`) was found to
mis-handle this real file — it treats the `.XV`'s Cartesian (Bohr)
positions as fractional coordinates without converting either the
positions or the cell from Bohr to Angstrom, which then fails downstream
with a cryptic ASE cell-shape error. This is a separate, pre-existing bug
in `translate.py`, out of scope for `stb-xrd`'s own migration — reported
rather than silently worked around. `sisl` itself reads this exact file
correctly (confirmed directly, and it's the same mechanism
`stb-workfunction`'s own vacuum-axis detection already uses successfully
on it), so the example script uses `sisl` + `core.structure_io.write_fdf`
directly instead.

## Running the walkthrough

```bash
cd examples/3.9-stb-xrd
./example_3.9.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder               | What it shows                                                              |
|-----------------------|-----------------------------------------------------------------------------|
| `basic/`              | A real structure's (CrS) space group/lattice info and pattern summary |
| `save-gnuplot/`       | `--save-gnuplot`: the stick-pattern `.dat`/`.gplot` pair, only with the flag |
| `compare/`            | `--compare-to` a synthetic mock-experimental pattern, similarity score |
| `self-check/`         | `--compare-to` the structure's OWN `xrd_pattern.dat` — the fixed column/peak-list-broadening bug, and `--raw-experimental` to see the old (wrong) score |
| `full-report/`        | Default (no report/gnuplot files) vs. `--save-report`, `references.bib` |

## Try it yourself

```bash
# A quick look at a candidate structure's pattern
stb-xrd --file structure.fdf --format fdf --save-gnuplot --view

# A different X-ray source, narrower range, top peaks only
stb-xrd --file structure.fdf --format fdf --wavelength MoKa --two-theta-range 10 80 --top 10

# Compare against a measured pattern
stb-xrd --file structure.fdf --format fdf --compare-to experimental.dat --view
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--file`/`--format` | Structure file and its format (`fdf` or `struct_out`).                 |
| `--wavelength`      | X-ray source name or a wavelength in Ang (default: CuKa).              |
| `--two-theta-range` | 2-theta range to scan, in degrees (default: 0-90).                    |
| `--top`             | Only include the N strongest peaks in the saved data file (`--save-gnuplot`). |
| `--compare-to`      | Experimental pattern file (2-theta first column, intensity last column); prints a similarity score. |
| `--raw-experimental` | With `--compare-to`, never auto-broaden a peak-list-looking experimental file — compare it exactly as read. |
| `-o/--output-dir`   | Where all generated files (and `references.bib`) land.                 |
| `--save-report`     | Persist the full report to `stb_xrd_report.txt`.                       |
| `--save-gnuplot`    | Also write `xrd_pattern.dat`/`xrd_pattern.gplot`.                      |
| `--view`            | Show an interactive matplotlib preview.                                |

## What's next

`stb-xrdsearch`/`stb-xrdrank` (Workflow) use this same simulated-pattern
machinery to guide a structure search against a real experimental
pattern, rather than just reporting one structure's pattern in isolation.
