# 3.2 — PDOS XML Parser (`stb-dos`)

## What this tool does

`stb-dos` reads a SIESTA `<label>.PDOS.xml` file (the projected density of
states) and writes it out as plain-text, per-orbital data ready to plot:

- `dos_total.dat` — the DOS summed over every atom and orbital.
- `dos_per_atom/<species>_<index>.dat` — one file per atom.
- `dos_per_species/dos_<species>.dat` — one file per chemical species.

Each file has one column per orbital angular momentum (`s`, `p`, `d`, `f`
with `--projection l`, the default; or the individual `m` sub-orbitals —
`px`, `py`, `pz`, `dxy`, ... — with `--projection ml`), split into
spin-up/spin-down columns for a spin-polarized (`nspin=2`) calculation.
`--type` picks which of the three output kinds to write (default: all
three).

## Why this matters (a bit of theory)

### What a PDOS.xml actually holds

Unlike a `.bands`/`.EIG` file (discrete eigenvalues at specific k-points),
a `PDOS.xml` is already a **broadened, continuous curve**: SIESTA applies
a Gaussian/Lorentzian smearing to the discrete eigenvalues and evaluates
the result on a shared energy grid, once per orbital of every atom. This
is what makes it directly plottable as a DOS curve, but it also means
there is no exact, discrete VBM/CBM sitting in the file the way there is
in a `.bands`/`.EIG` — see `--shift vbm/cbm` below.

### `--shift`: which energy sits at 0 eV in the output

- `fermi` (default) — the calculation's own Fermi energy, from
  `<fermi_energy>` in the XML.
- `vbm` / `cbm` — the Valence Band Maximum / Conduction Band Minimum (see
  next section).
- a manual number — any value you choose.

This only changes where the written energy axis's zero sits — it never
changes the shape of the DOS curve itself.

### `--shift vbm/cbm`: reusing `stb-bands`' own VBM/CBM machinery

Since a PDOS.xml alone has no exact band-edge value, `--shift vbm/cbm`
looks for a **companion file** next to it and computes the VBM/CBM from
that, using the exact same `core/siesta_bands.py` code `stb-bands` itself
uses — so the number reported here always matches what `stb-bands` would
report for the same file, never a second, independently-computed value.
One `--label` covers all three files at once (`<label>.PDOS.xml` +
`<label>.bands`/`<label>.EIG`); if you pass a `PDOS.xml` path directly
instead, the label is still derived automatically from its own filename.

The lookup follows a strict hierarchy, most-precise source first:

1. **`<label>.bands`** — the k-**path** file `stb-bands` defaults to.
   Tried first.
2. **`<label>.EIG`** — the full SCF k-**mesh**, read via `sisl`. Used only
   if no `.bands` is present.
3. **An estimate from the DOS itself** (`--estimate-from-dos`) — the last
   resort, only used if *neither* of the above exists, and only when
   explicitly opted into.

A `.bands` k-path and a `.EIG` k-mesh are, in general, different samplings
of the Brillouin zone (the same distinction `stb-bands`' own `--eig-file`
mesh-vs-line comparison is built around) — they can legitimately report
different VBM/CBM values for the same calculation, which is exactly why
the hierarchy always prefers a real eigenvalue file over the DOS-based
guess, and why a mismatch between the two files' own Fermi energies (or
`nspin`) is reported as a warning, not silently ignored.

Without either file, and without `--estimate-from-dos`, `stb-dos` refuses
outright rather than silently substituting a different-physics
approximation for an exact value.

### The DOS-estimate fallback (`--estimate-from-dos`)

When neither a `.bands` nor a `.EIG` is available, `--estimate-from-dos`
walks outward from the Fermi energy on each side of the (summed) total
DOS until it first rises above `--dos-threshold-frac * peak_DOS` (default
`0.01`, i.e. 1% of the tallest point in the curve) — the energy where the
curve leaves the near-zero gap region. This is explicitly an
**approximation**: the true band edge is blurred by the PDOS's own
broadening, so the estimated value can differ meaningfully from the exact
one a `.bands`/`.EIG` would give. `stb-dos` always prints an explicit
warning when this path is used, rather than presenting it as exact.

## The report: console output, `--save-report`, `--save-gnuplot`, `--view`

Every run prints a numbered report — the same `[0]...[5]` style every
newer tool in this suite uses:

- `[0] RUN METADATA` — input file, resolved label, requested `--type`(s),
  projection mode, shift mode, output directory.
- `[1] INPUT DATA` — Fermi energy, spin channels, energy grid, orbitals/
  atoms/species found, plus any non-fatal parsing warnings.
- `[2] ENERGY SHIFT` — which reference was used and why (the `--shift
  vbm/cbm` hierarchy messages above live here).
- `[3] WRITING OUTPUT FILES` — one `[OK]` line per `.dat`/directory
  actually written.
- `[4] REFERENCES` — SIESTA citations, written to `references.bib`.
- `[5] SUMMARY & FILES` — every file this run produced.

`--save-report` additionally persists that exact report to
`stb_dos_report.txt` — off by default, so a plain run only ever writes the
`.dat` files and `references.bib`, no text report file.

`--save-gnuplot` writes one `.gplot` script **per output category** — off
by default:

- `dos_total.gplot` (next to `dos_total.dat`) — plots every orbital/spin
  column of the total DOS, using gnuplot's own `columnheader(i)` to read
  column names straight from the `.dat` file's own
  `#Energy(eV)  s  p  ...` header line. Since there's only ever one
  `total` file, this is the most useful single-file view.
- `dos_per_atom/dos_per_atom.gplot` — ONE script overlaying every atom's
  own total DOS (that atom's orbital/spin columns collapsed into a
  single curve via gnuplot's `sum [i=2:N] column(i)`), so atoms can be
  compared against each other in one plot instead of opening one file
  per atom.
- `dos_per_species/dos_per_species.gplot` — the same idea, one curve per
  chemical species.

Run `gnuplot <name>.gplot` (from inside the same folder as the script)
for a PDF.

`--view` opens an interactive matplotlib preview right before the tool
exits — off by default, one figure per output category actually written
(`total`/`atom`/`species`, mirroring `--save-gnuplot`'s own grouping),
shown together only after every report line and file has already been
written, so a blocking preview window never delays or hides them.

## When you'd reach for it

- Plotting a projected density of states from a finished SIESTA
  calculation, split by atom/species/orbital.
- Referencing the DOS plot to the VBM or CBM instead of the Fermi energy
  — e.g. to line it up visually with a `stb-bands` plot shifted the same
  way.
- A quick, no-`.bands`/`.EIG`-available estimate of where the band edges
  sit, when only a `PDOS.xml` exists.

## Two ways to run it

A — direct CLI:
```bash
stb-dos --label example --shift vbm
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.2
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.2.sh` proves this directly at the end. The menu
additionally offers `--shift` as a **numbered menu** (`1` Fermi level,
`2` VBM, `3` CBM, `4` custom value — re-prompting on an invalid choice)
instead of typing the keyword, plus separate `y/N` prompts for
`--save-report`/`--save-gnuplot`/`--view`, each mapped straight onto the
matching CLI flag.

## Files in this folder

- `example.PDOS.xml` — a small, hand-built synthetic fixture: 1 atom
  (Si), 1 orbital (`s`), 13 energy points from -3 to +3 eV, Fermi energy
  0.0 eV. Its DOS drops to near-zero between -1 and +1 eV — a deliberate
  "gap" shape, used to demonstrate `--estimate-from-dos`.
- `example.bands` / `example.EIG` — a matching "same calculation" trio
  (same Fermi energy), with known-by-hand VBM/CBM at each level of the
  `--shift vbm/cbm` hierarchy.
- `multi.PDOS.xml` — a second fixture, 3 atoms (`Si_1`, `Si_2`, `O_3`,
  2 species), same energy grid/"gap" shape as `example.PDOS.xml` but each
  atom's DOS scaled differently (`x1`/`x2`/`x0.5`) so the three curves
  are visibly distinct — used to demonstrate `--save-gnuplot`/`--view`
  actually overlaying multiple atoms/species in one script/figure
  (`example.PDOS.xml`'s own single atom can't show that).
- `example_3.2.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/2-dos/test.sh` for that, which additionally covers
  a real ~22 MB reference PDOS.xml this example deliberately skips).
- `.gitignore` — excludes `output/`.

## Running the walkthrough

```bash
cd examples/3.2-stb-dos
./example_3.2.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Five self-contained cases are generated:

| Folder                | What it shows                                                              |
|-------------------------|-----------------------------------------------------------------------------|
| `basic-run/`           | The three output types (`dos_total.dat`, `dos_per_atom/`, `dos_per_species/`) |
| `vbm-hierarchy/`       | `--shift vbm`'s `.bands` → `.EIG` → `--estimate-from-dos` hierarchy, live  |
| `label-shorthand/`     | `--label` resolving both the `PDOS.xml` and its `.bands` in one go        |
| `multi-atom-plots/`    | `--save-gnuplot`/`--view` overlaying 3 atoms/2 species in one script/figure each, rendered with a real `gnuplot` |
| `full-report/`         | Default (no report/`.gplot` files) vs. `--save-report --save-gnuplot`     |

## Try it yourself

```bash
# Your own finished SIESTA calculation
stb-dos --label my_calc --shift vbm

# No .bands/.EIG available for my_calc -- opt into the DOS-based estimate
stb-dos --label my_calc --shift cbm --estimate-from-dos
```

## Flag reference

| Flag                   | Meaning                                                                    |
|-------------------------|-----------------------------------------------------------------------------|
| `filename`              | Path to the `.PDOS.xml` file (positional). Optional if `--label` is given. |
| `--label`               | SIESTA label; shorthand for `filename='<label>.PDOS.xml'`, and (with `--shift vbm/cbm`) used to locate `<label>.bands`/`<label>.EIG`. |
| `--type`                | Which output(s) to write: `total`, `atom`, `species` (default: all three). |
| `--shift`               | Energy reference: `fermi` (default), `vbm`, `cbm`, or a manual number.    |
| `--estimate-from-dos`   | With `--shift vbm/cbm`: fall back to a DOS-based estimate if no `.bands`/`.EIG` is found. Off by default. |
| `--dos-threshold-frac`  | Threshold (fraction of peak DOS) used only by `--estimate-from-dos` (default `0.01`). |
| `--projection`          | Orbital detail: `l` (s/p/d/f, default) or `ml` (individual m sub-orbitals). |
| `-o/--output-dir`       | Where to write the output files (default: current directory).            |
| `--save-report`         | Persist the full numbered report to `stb_dos_report.txt`. Off by default. |
| `--save-gnuplot`        | Write one `.gplot` script per output category (total/atom/species). Off by default. |
| `--view`                | Show an interactive matplotlib preview, one figure per output category. Off by default. |

## What's next

The VBM/CBM machinery this tool reuses (`core/siesta_bands.py`) is shared
with `stb-bands`, `stb-fatbands`, and — via `select_band_vbm_cbm` — with
`stb-effmass`/`stb-wfdensity`. `stb-bands` covers the same occupied/empty
-states question from the discrete-eigenvalue side instead of the DOS.
