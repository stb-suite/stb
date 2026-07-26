# 3.1 — Bands Analyzer (`stb-bands`)

## What this tool does

`stb-bands` parses a SIESTA `.bands` file (band-structure eigenvalues
along a high-symmetry k-path) and reports the Valence Band Maximum
(VBM), Conduction Band Minimum (CBM), whether the gap is direct or
indirect, and whether the material is metallic — no SIESTA re-run
needed. Optionally cross-checks that gap against the full SCF k-mesh
(`--eig-file`, a `.EIG` file) in case the path missed the true extrema,
and reports each spin channel (+ half-metallic character) separately for
spin-polarized calculations.

**This is the first example for the Analysis category (3.x)** — the
same guided-walkthrough format 2.x already uses, applied to a tool that
*analyzes* an existing SIESTA output rather than building/transforming a
structure.

## Manual code + physics review (this session)

Before touching anything, `bands.py` and its shared parsing/physics
module (`core/siesta_bands.py`) were read in full and checked line by
line. **No bugs were found** — the parsing, VBM/CBM search, direct/
indirect classification, metal detection, and spin-channel handling are
all correct (see below for exactly what was checked). The one real gap
found was **presentation**, not physics: unlike almost every other
`*_analysis.py` tool in the suite, `bands.py` used to write
`bands_analysis.txt` *unconditionally*, with no `--save-report` opt-in
and no numbered `[0]...[N]` report — both fixed in this pass, with the
underlying numbers completely unchanged.

## Why this matters (a bit of theory)

### The `.bands` file format

Line 1 is the Fermi energy; line 4 is `nbands nspin nk`; then one block
per k-point (the k value, followed by `nbands*nspin` energies — wrapped
onto extra lines if there are many bands); a footer lists which of those
k-points are the high-symmetry ones, with labels (e.g. `'GAMMA'`, `'X'`).
Spin-up and spin-down values are **not interleaved per band** — the
first `nbands` values at each k are the full spin-up channel, the next
`nbands` (if `nspin=2`) are spin-down (confirmed against `sisl`'s own
reader).

### VBM, CBM, and direct vs. indirect gaps

At **each** k-point, the highest energy at or below the Fermi level is
that k-point's local VBM; the lowest energy above it is the local CBM.
- The **indirect (fundamental) gap** is `CBM - VBM` using the *global*
  extrema — the true VBM might be at a completely different k-point than
  the true CBM.
- The **direct gap** is the smallest `CBM(k) - VBM(k)` restricted to the
  *same* k-point.

Since the direct gap is the same search with an extra constraint
(same k), **indirect gap ≤ direct gap always** — `example_3.1.sh`'s first
case shows this on a real 21-band structure (indirect 1.2457 eV, direct
1.4484 eV, so `Gap type: Indirect`).

### Why a "gap" is never negative, and what `--gap-tol` actually decides

If bands cross the Fermi level, the naive `CBM - VBM` can come out
negative (global CBM below global VBM) — `stb-bands` clamps this to
exactly `0.0` rather than ever reporting a negative gap. Very small
*positive* gaps and true metals are otherwise indistinguishable, so
`--gap-tol` (default `0.01` eV) is the threshold: an indirect gap below
it is classified `Metallic`. `example_3.1.sh`'s
`metallic-threshold/` case proves this is a genuine threshold effect —
the *same* real ~1.25 eV gap is classified `Indirect` at a tight
`--gap-tol` and `Metallic` at a loose one; the physical gap itself never
moved.

### Half-metals: one spin channel with a gap, the other without

For a spin-polarized calculation, `stb-bands` runs the exact same
VBM/CBM search independently on each spin channel. A **half-metallic**
material has one channel metallic and the other with a real gap — it
conducts current of only one spin orientation, which is exactly why
half-metals matter for spintronics. `example_3.1.sh`'s `half_metal.bands`
fixture is a purpose-built synthetic example: spin-up channel gap
~0.0068 eV (Metallic, below the default `--gap-tol`), spin-down channel
gap exactly 2.0 eV (Direct) — a clean, unambiguous half-metal.

### Why a k-mesh comparison can only go one direction

A `.bands` file only samples a 1-D **path** through the Brillouin zone —
the true VBM/CBM could sit anywhere else in the zone, entirely off that
path. A `.EIG` file (the full SCF k-**mesh**, `--eig-file`) samples much
more of the zone, so it can catch this. The key asymmetry: a denser mesh
can only find a fundamental gap **smaller than or equal to** any 1-D
subset of it (the path *is* a subset of the mesh in the ideal case),
never larger. So:
- mesh gap **smaller** than the line gap → expected and informative: the
  path missed the true extrema, trust the mesh.
- mesh gap **larger** than the line gap → a red flag: the mesh itself is
  too coarse to have captured what the path already found.

`example_3.1.sh`'s `mesh-vs-line/` case demonstrates the first (smaller)
case directly, with a synthetic 3-point mesh whose extrema are known by
hand (`mesh.EIG`/`mesh.KP`).

## When you'd reach for it

- Getting the band gap and its direct/indirect character from a finished
  SIESTA band-structure calculation, without opening a plotting tool.
- Sanity-checking whether your k-path actually captured the true gap
  (`--eig-file`), before quoting a gap value from a path alone.
- Checking for half-metallic character in a spin-polarized calculation.

## Two ways to run it

A — direct CLI:
```bash
stb-bands --file semiconductor.bands --shift fermi
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.1
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.1.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[7] SUMMARY & FILES`)
  printed to the console.
- **`bands_gnuplot.dat` + `bands.gplot`** — the actual plot data (run
  `gnuplot bands.gplot` for a PDF).
- **`references.bib`** — SIESTA (every `.bands` file is SIESTA output).

## Optional (off by default)

- **`--save-report`** — also persists the full numbered report to
  `stb_bands_report.txt`. The old, always-on `bands_analysis.txt` is
  gone entirely — no flag brings it back under that name.
- **`--eig-file`** (+ `--kp-file`) — k-mesh vs. k-path gap comparison.
- **`--gap-tol`** — the Metallic/Indirect-or-Direct threshold (default
  `0.01` eV).
- **`--label`** — shorthand for `--file <label>.bands`, auto-detecting a
  matching `<label>.EIG`/`<label>.KP` if present.

## Files in this folder

- `semiconductor.bands` — a real 21-band structure, 4 k-points
  (Γ→X→M→Γ), with a genuine indirect gap. Copied from
  `test/3-analysis/1-bands/siesta_edge_nbands21.bands` (kept small
  deliberately — a real `.bands` file for a decent-sized cell is easily
  1+ MB, far too heavy for a lightweight example).
- `half_metal.bands` — synthetic spin-polarized fixture, purpose-built
  for the half-metallic demo (from `siesta_spin2.bands`).
- `mesh.EIG` / `mesh.KP` — a tiny synthetic 3-point k-mesh with
  known-by-hand extrema (from `mesh_kxkyz.{EIG,KP}`).
- `example_3.1.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/1-bands/test.sh` for that, which additionally
  covers the real, ~1.2 MB reference fixture this example deliberately
  skips).
- `.gitignore` — excludes `output/` and every file this tool itself
  generates on a run (`bands_gnuplot.dat`, `bands.gplot`,
  `references.bib`, `stb_bands_report.txt`).

## Running the walkthrough

```bash
cd examples/3.1-stb-bands
./example_3.1.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Five self-contained cases are generated:

| Folder                  | What it shows                                                          |
|--------------------------|-------------------------------------------------------------------------|
| `gap-analysis/`          | VBM/CBM and direct-vs-indirect gap on a real (tiny) indirect-gap structure |
| `metallic-threshold/`    | The same real gap classified Indirect vs. Metallic, only `--gap-tol` changes |
| `half-metallic/`         | One spin channel metallic, the other with a 2 eV gap                   |
| `mesh-vs-line/`          | `--eig-file`/`--kp-file`: the k-mesh catching what the k-path missed    |
| `full-report/`           | Default (no report file) vs. `--save-report` (`stb_bands_report.txt`), `references.bib` |

## Try it yourself

```bash
# Your own finished SIESTA calculation
stb-bands --label my_calc --shift fermi --save-report

# Cross-check the path-based gap against the full SCF mesh
stb-bands --file my_calc.bands --eig-file my_calc.EIG --shift fermi
```

## Flag reference

| Flag              | Meaning                                                                |
|--------------------|-------------------------------------------------------------------------|
| `--label`          | SIESTA label; shorthand for `--file <label>.bands` + auto-detected `.EIG`/`.KP` |
| `--file`           | Explicit path to the `.bands` file                                     |
| `--shift`          | Energy reference for the plot: `vbm`/`cbm`/`fermi`/`manual`             |
| `--manual-value`   | Custom shift value (required with `--shift manual`)                    |
| `--gap-tol`        | Metallic-classification threshold, eV (default `0.01`)                 |
| `--eig-file`       | Full SCF k-mesh eigenvalues, for the mesh-vs-line gap comparison        |
| `--kp-file`        | Cartesian k-points matching `--eig-file` (requires `--eig-file`)        |
| `-o/--output-dir`  | Where to write the plot data (and report, with `--save-report`)        |
| `--save-report`    | Persist the full numbered report to `stb_bands_report.txt`             |

## What's next

`core/siesta_bands.py`'s VBM/CBM machinery is shared with `stb-fatbands`
(orbital-projected bands) and, via `select_band_vbm_cbm`, with
`stb-effmass`/`stb-wfdensity` (finding the VBM/CBM band directly from a
`.WFSX` with no `.bands` file at all). `stb-dos` covers the same
occupied/empty-states question from the density-of-states side instead
of individual bands.
