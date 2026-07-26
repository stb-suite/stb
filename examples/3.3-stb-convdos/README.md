# 3.3 — DOS Processor / Convolution (`stb-convdos`)

## What this tool does

`stb-convdos` applies Gaussian broadening (convolution) to every DOS
column of a whitespace-separated data file (energy first, one or more DOS
columns after — typically `stb-dos`'s own output). No structure, no
SIESTA re-run — it's pure post-processing of an already-computed density
of states.

Two input modes:

- **`--file <file> --out <file>`** — broaden one file.
- **`--dir <folder>`** — broaden **every** `.dat` file found recursively
  under a folder (e.g. a whole `stb-dos` output tree: `dos_total.dat` +
  `dos_per_atom/*.dat` + `dos_per_species/*.dat`) with the exact same
  `--sigma`/`--fwhm`/`--size`, writing a broadened mirror of the entire
  tree — same subfolders, same filenames — under `--output-dir` (default:
  `<folder>_filtered`, a sibling directory). A file that fails validation
  (e.g. accidentally not a DOS file) is skipped with a warning instead of
  aborting the whole batch.

## Why this matters (a bit of theory)

### Why broaden a DOS at all?

A DFT calculation gives you the density of states at a finite, discrete
set of k-points. Especially with a coarse k-mesh, the "raw" DOS looks
like a comb of sharp, almost delta-function-like spikes rather than a
smooth curve — an artifact of finite sampling, not the real physical
density of states (which, for a bulk crystal, genuinely is smooth).
Convolving with a Gaussian kernel simulates the effect that any real
measurement or a denser calculation would show: finite instrumental/
thermal broadening, or just enough smoothing to see the underlying trend
instead of sampling noise. `example_3.3.sh`'s first case makes this
completely explicit with a synthetic "stick spectrum" (`dos_total.dat`
below), watching sharp spikes turn into smooth peaks.

### `sigma` vs. `FWHM`

Both describe the *width* of the Gaussian broadening kernel, just in two
conventions:

- **`sigma`** (standard deviation) — the natural parameter of a Gaussian,
  `exp(-(E-E0)^2 / (2*sigma^2))`.
- **`FWHM`** (Full Width at Half Maximum) — the width of the peak measured
  at half its height, the value more commonly quoted in spectroscopy and
  experimental resolution specs.

They're related by an exact constant:
```
FWHM = 2*sqrt(2*ln(2)) * sigma  ≈  2.3548 * sigma
```
`--fwhm` is simply converted to the equivalent `sigma` internally
(`sigma = fwhm / 2.3548`) before anything else happens — `--sigma 50` and
`--fwhm 117.741` are the exact same broadening, just expressed
differently (`example_3.3.sh` proves this directly).

### Why `--sigma`/`--fwhm` are in meV, not "samples"

The convolution kernel itself operates on **grid samples** — a "width of
5" only means something once you know how far apart the file's own
energy points actually are. `stb-convdos` instead asks for `--sigma`/
`--fwhm` in **meV**, a physical, file-independent unit, and converts it
internally:
```
sigma_samples = (sigma_meV / 1000) / d_energy
```
where `d_energy` is the **median** spacing between consecutive energy
points in that specific file (median, not mean — robust to the rare
irregular point). This is why the same `--sigma 50` gives the same
physical broadening whether the input file has 100 or 10,000 energy
points — and why, in `--dir` mode, files with different native energy
grids (e.g. a `stb-dos` run with a different grid resolution) can each
get a different `sigma_samples`/kernel size internally while still
representing the identical 50 meV physical width (see `[2] BROADENING
PARAMETERS` in the report, which shows this per-file when it happens).

### The kernel size, and why it must be odd

By default, the kernel width (in samples) is sized automatically to
cover about **3 standard deviations** on each side of center — over 99%
of a Gaussian's mass:
```
half_width = ceil(3 * sigma_samples)
size = 2 * half_width + 1
```
It's forced to be **odd** so there's a single, unambiguous center sample
— an even-sized kernel would shift the filtered curve by half a bin.
`--size` overrides this if you want a specific width; `stb-convdos`
rejects an even or non-positive value outright rather than silently
producing a shifted or garbage result.

### The conservation check

A properly normalized Gaussian kernel (this one sums to exactly 1) should
barely change the **total integrated DOS** of each column — broadening
redistributes weight in energy, it doesn't create or destroy states.
`[3] CONSERVATION CHECK` integrates (trapezoidal rule) each column before
and after broadening and reports both numbers side by side. A **large**
drift almost always means the energy window is too narrow for the
requested broadening — weight simply falls off the zero-padded edges of
the window and is lost. This is a genuine runtime physics sanity check,
not just a formatting nicety: if you ever see the "before" and "after"
numbers diverge by more than a percent or so, trust the warning, not the
filtered file.

### `--dir`: broadening a whole `stb-dos` output tree at once

`stb-dos` (3.2) can write dozens of files for a real structure —
`dos_total.dat`, one `dos_per_atom/<species>_<index>.dat` per atom, one
`dos_per_species/dos_<species>.dat` per species. Running `stb-convdos`
by hand on every one of them individually doesn't scale. `--dir` walks
the whole tree recursively, applies the identical broadening to every
`.dat` file found, and writes a **mirrored** tree (same relative paths,
same filenames) under `--output-dir` — a drop-in "broadened copy" of the
original `stb-dos` output, ready to plot the same way.

Only the **plot** behaves differently in `--dir` mode: with dozens of
files, showing a before/after figure for every single one would be
impractical (and likely to overwhelm the display). Only the **first**
file processed gets its plot shown; the rest are still broadened and
written normally, just without a plot window. `--no-plot` turns off that
one remaining plot too.

## When you'd reach for it

- Making a jagged, under-sampled DOS from a quick/coarse SIESTA
  calculation presentable, without re-running with a denser k-mesh.
- Matching a calculated DOS's resolution to an experimental spectrum's
  known instrumental broadening (if you know its FWHM).
- Broadening an entire `stb-dos` output folder (total + per-atom + per
  -species) in one command instead of one `stb-convdos` call per file.

## Two ways to run it

A — direct CLI:
```bash
stb-convdos --file dos_total.dat --sigma 50 --out dos_filtered.dat
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.3
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.3.sh` proves this directly at the end. The menu adds
one extra first choice (single file vs. a whole folder) before the usual
`--sigma`/`--fwhm`/`--size`/plot prompts.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[6] SUMMARY & FILES`)
  printed to the console.
- **The filtered `.dat` file(s)** — one for `--file`, a mirrored tree for
  `--dir`.
- **`references.bib`** — SIESTA (every DOS file processed here, directly
  or via `stb-dos`, is ultimately derived from a SIESTA calculation).

## Optional (off by default)

- **`--save-report`** — also persists the full numbered report to
  `stb_convdos_report.txt`.
- **`--size`** — explicit kernel width in samples (must be odd), instead
  of the automatic 3-sigma sizing.
- **`--output-dir`** (`--dir` mode only) — where the broadened tree goes;
  default `<dir>_filtered`.
- **`--no-plot`** — skip the before/after plot entirely (`--dir`: skips
  the one plot it would otherwise show, for the first file).

## Files in this folder

- `dos_total.dat` — a small, hand-built synthetic "stick spectrum": 41
  energy points (-2 to +2 eV, 0.1 eV spacing), two columns (`s`, `p`)
  that are exactly zero except for 2 sharp spikes each — deliberately
  built to make the broadening effect (spike -> smooth peak) and the
  conservation check (`before`/`after` areas match exactly) easy to see
  and verify by hand.
- `multi.PDOS.xml` — a small synthetic PDOS (3 atoms/2 species), the same
  fixture `examples/3.2-stb-dos/` uses — processed live by
  `example_3.3.sh` via a real `stb-dos` run to generate a genuine
  `dos_total.dat` + `dos_per_atom/*.dat` + `dos_per_species/*.dat` tree,
  used to demonstrate `--dir` on the exact use case it's built for.
- `example_3.3.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/3-dos-convolution/test.sh` for that).
- `.gitignore` — excludes `output/`.

## Running the walkthrough

```bash
cd examples/3.3-stb-convdos
./example_3.3.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder                | What it shows                                                              |
|-------------------------|-----------------------------------------------------------------------------|
| `basic-broadening/`    | Sharp spikes -> smooth peaks, the numbered report, the conservation check |
| `sigma-vs-fwhm/`       | `--sigma 50` and the equivalent `--fwhm 117.741` give the identical broadening |
| `dir-mode/`            | A real `stb-dos` output tree broadened in one `--dir` command, mirrored structure, one skipped-file example |
| `full-report/`         | Default (no report file) vs. `--save-report` (`stb_convdos_report.txt`), `references.bib` |

## Try it yourself

```bash
# Broaden one file from a real stb-dos run
stb-convdos --file dos_total.dat --sigma 50 --out dos_total_filtered.dat

# Broaden an entire stb-dos output folder in one command
stb-convdos --dir my_calc_dos/ --sigma 50 --save-report
```

## Flag reference

| Flag              | Meaning                                                                |
|--------------------|-------------------------------------------------------------------------|
| `--file`           | Input file to broaden. Requires `--out`. Alternative to `--dir`.       |
| `--dir`            | Input folder to broaden recursively (every `.dat` found). Alternative to `--file`. |
| `--sigma`          | Gaussian standard deviation, meV. Alternative to `--fwhm`.             |
| `--fwhm`           | Gaussian full width at half maximum, meV (`= 2.3548 * sigma`). Alternative to `--sigma`. |
| `--size`           | Explicit kernel width in samples (must be odd). Default: auto (3-sigma rule). |
| `--out`            | Output file. Required with `--file`, invalid with `--dir`.             |
| `-o/--output-dir`  | Output tree root for `--dir` (default `<dir>_filtered`). Invalid with `--file`. |
| `--no-plot`        | Skip the before/after plot.                                            |
| `--save-report`    | Persist the full numbered report to `stb_convdos_report.txt`.          |

## What's next

`stb-convdos` is a pure signal-processing post-processor — it has no
opinion about *where* its input DOS files came from, only that they're
whitespace-separated columns with energy first. `stb-dos` (3.2) is its
natural upstream partner (and the `--dir` mode above is built specifically
around that pairing), but any DOS-like data file in the same column
convention works the same way.
