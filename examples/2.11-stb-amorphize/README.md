# 2.11 — Amorphous Structure Generator (`stb-amorphize`)

## What this tool does

`stb-amorphize` runs a melt-quench protocol with the **MACE-MP-0**
machine-learned potential (or a custom model you fine-tuned yourself with
`stb-mlffAnalysis`) to turn a crystalline structure into a fast, heuristic
starting guess for an amorphous/glassy one — meant to give SIESTA a much
better starting point than an ad-hoc random-displacement structure. **This
is a heuristic, not a substitute for DFT.** NOT a substitute for a slower,
production-quality quench (more MD steps, a slower cool-down) either.
Bulk (3D periodic) structures only — melting a vacuum-padded slab/wire/
molecule is physically meaningless, and the input is rejected outright if
a vacuum-padded axis is detected.

## Why this matters (a bit of theory)

### Melt-quench: why heating something up is how you make it *lose* order

An amorphous/glassy structure (real short-range chemistry — roughly the
right bond lengths/angles/coordination — but no long-range crystalline
order) is hard to build directly: you can't just "un-order" atomic
positions at random without creating physically absurd overlaps. The
classic laboratory trick is to melt the material, then cool it faster
than it can re-crystallize. It works because a liquid has already
forgotten its crystalline arrangement while keeping sensible local
chemistry (atoms still push each other apart at the right distance); cool
it fast enough and there isn't time for atoms to find their way back into
an ordered lattice before the structure "freezes" kinetically — a glass,
not a crystal. `stb-amorphize` reproduces exactly this in silico: melt a
crystalline structure well above its real melting point (so it melts
even given the tiny simulated time available — real quenches take
seconds; this one takes picoseconds), hold it there long enough to
forget the crystal, then ramp the temperature back down.

### NPT-Berendsen has a *finite response time* — and that has two real consequences

The MD engine underneath is ASE's `NPTBerendsen`: not plain
constant-energy dynamics, but a thermostat (`--taut`) **and** a barostat
(`--taup`) that actively drive the system toward a target temperature and
0 GPa pressure — a liquid/quenched structure should relax its own volume
freely, not stay pinned at the crystal's original cell size. Both
`--taut`/`--taup` are **response times** (in fs), not instantaneous
snaps: Berendsen coupling relaxes the current temperature toward the
target roughly exponentially, with `--taut` as the time constant
(default 50 fs). That single fact drives two things worth knowing before
you pick step counts:

1. **`--melt-steps` needs to comfortably exceed `--taut`** for the system
   to actually reach the melt temperature, not just start drifting toward
   it. Measured live on the same 8-atom `si8.fdf` (`--seed 7` fixes the
   initial velocities so both runs are directly comparable — see
   `output/too-short-vs-adequate-melt/` below):

   | Protocol | `--melt-steps` (vs. `--taut`=50 fs) | Bond-angle std after MD |
   |---|---|---:|
   | Too short | 5 (0.1x `--taut`) | 3.07 deg |
   | Adequate  | 100 (2x `--taut`) | 8.64 deg |

   Both runs' bond-angle **mean** stays close to the crystal's own
   ~109.5 deg tetrahedral angle (short-range chemistry survives either
   way) — only the **std** (long-range disorder) tracks how much actual
   melting happened, exactly the diagnostic's intended signature.

2. **A `--quench-steps` ramp faster than `--taut` leaves the actual,
   instantaneous temperature lagging behind the falling target** — real
   residual thermal energy/force, not a clean local minimum. This is
   *why the final static relax exists at all* (see below) — measured
   live (`--melt-steps 30 --quench-steps 60 --seed 7`, `output/
   final-relax-comparison/`): after the quench ramp nominally targets
   300 K, the actual instantaneous temperature is still **760 K**.

### Why the final static relax is on by default

Following directly from the lag above: the quench MD alone leaves real
thermal motion, not a converged geometry. `stb-amorphize` always runs a
final position+cell relax afterward (`--no-final-relax` to skip it).
Measured live on the same run as above:

| Quantity | Right after quench MD | After final relax | Change |
|---|---:|---:|---:|
| Energy (eV) | -41.094170 | -42.966005 | -1.871835 (-0.233979/atom) |
| Max force (eV/Ang) | 1.5434 | 0.0093 | -1.5340 |

A max force of ~1.5 eV/Ang is real thermal noise (atoms still
mid-vibration when the MD stopped) — a handful of FIRE steps bring it
down to the same ~0.01-0.05 eV/Ang order of magnitude every other
stb-suite relaxation targets.

### The bond-angle mean/std check: how the tool proves amorphization happened

A crystal has a sharp, well-defined bond-angle distribution (std near
0 deg — every tetrahedral Si-Si-Si angle in `si8.fdf` is exactly
109.47 deg, since it's a perfect diamond-cubic cell). A genuine
melt-quench should **broaden** that distribution while keeping roughly
the same **mean** (short-range coordination survives; long-range order
doesn't) — the table above is this signature measured directly, not
just asserted.

### Symmetry before/after: expected to differ, unlike most other stb-suite tables

`[6] SYMMETRY ANALYSIS (BEFORE / AFTER)` is deliberately the *opposite*
kind of check from `stb-unitcell`'s own before/after table: there, an
unexpected symmetry change is a red flag; here, losing the crystal's
symmetry is the entire point. A genuinely melt-quenched structure often
collapses all the way down to `P1`/`P-1` (essentially no symmetry
operations left besides identity, or identity+inversion). Measured live
(`--melt-steps 100 --quench-steps 200 --seed 2`, `output/full-report/`):
`Fd-3m (227)` (diamond-cubic Si) &rarr; `P-1 (2)`.

### Why an 8-atom cell? (and why every number above is a qualitative illustration)

This example deliberately reuses one tiny 8-atom bulk Si cell
throughout so every comparison finishes in seconds on a CPU. Real
amorphization needs a much bigger supercell (tens to hundreds of atoms)
— with only 8 atoms, a single MD trajectory is a genuinely chaotic,
small-number-statistics system: re-running the exact same command
*without* `--seed` gives visibly different bond-angle numbers each time,
and even `--seed` only fixes the *initial* velocities, not every source
of run-to-run floating-point noise. Every comparison above still shows a
real, physically sensible **trend** — that's the point being
demonstrated — just don't expect the precise digits to reproduce
bit-for-bit on your machine, or to mean much on their own for a cell
this small.

### `--save-data`: keeping the MD trace that used to be thrown away

Until this feature, every intermediate MD step was discarded — only the
final frame was ever written anywhere, and the only trace of
temperature/energy during the run was a `\r`-overwritten progress line
on stderr that vanished the moment the run ended. `--save-data` writes
`<stem>_md_diagnostics.dat` (step, time, temperature,
E_pot/E_kin/E_total, cell volume) with **one continuous step/time axis
spanning both the melt and quench stages**, sampled every `--stride`
steps, using gnuplot's own `index` block convention — the melt stage is
`index 0`, the quench stage is `index 1` **in the same file**, so either
stage can be plotted alone, or both together via `index 0:1`, without
re-running anything. The companion `.gplot` (same `.dat`+`.gplot`
convention used throughout the suite, e.g. `4.1-strain/analysis/`)
renders both energy and temperature vs. step into one PDF, with the
melt&rarr;quench transition marked and each stage's target temperature
drawn as a reference line.

### `--save-traj`: watching the disorder happen, not just trusting a number

`--save-traj` writes a real multi-frame trajectory of the whole MD run
(same sampling as `--save-data`) — the same 3-format choice as
`stb-ani2traj`/`stb-mlmd`: `xsf` (OVITO/VMD-native, default), `pdb`
(VMD's own default), or `xyz` (OVITO-native). The bond-angle std number
proves disorder happened; actually opening this trajectory in OVITO or
VMD lets you *see* the crystal melt and re-quench, frame by frame.

## When you'd reach for it

- A better starting guess than a random-displacement structure for an
  amorphous/glassy SIESTA calculation (a-Si, a-SiO2, glassy alloys, ...).
- `--save-data`/`--save-traj` specifically: diagnosing whether a given
  `--melt-steps`/`--quench-steps`/`--melt-temp` protocol is actually doing
  what you think, before committing to a much bigger production cell.

## Two ways to run it

A — direct CLI:
```bash
stb-amorphize -f si8.fdf --save-data --save-traj
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.11
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.11.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[9] SUMMARY & FILES`) printed
  to the console and, with `--save-report`, also saved to
  `stb_amorphize_report.txt`.
- **Structure validation + bond-angle mean/std**, before and after the MD.
- **Symmetry analysis before/after** (`[6]`) — expected to differ (often
  down to `P1`/`P-1`), unlike most other stb-suite before/after tables.
- **`references.bib`** — SIESTA + MACE (+ MACE-MP if a foundation model
  was used, not a custom fine-tuned one).
- **A provenance header** written into the output `.fdf`: input file,
  melt/quench protocol, bond-angle before/after, final-relax summary.

## Optional (off by default)

- **`--save-data`** — writes `<stem>_md_diagnostics.dat` + a companion
  `.gplot` (step, time, temperature, E_pot/E_kin/E_total, cell volume) for
  the whole melt-quench process, one gnuplot `index` block per stage.
  Render it with `gnuplot <stem>_md_diagnostics.gplot`.
- **`--save-traj`** (+ `--traj-format xsf|pdb|xyz`, default `xsf`) — a
  multi-frame trajectory of the whole MD run, for viewing in OVITO/VMD.
- **`--stride`** — sampling interval (in MD steps) shared by both exports
  above (default `10`).
- **`--no-final-relax`** — skip the final static (position+cell) relax.
- **`--save-report`** / **`--view`** — same as every other stb-suite tool.

## Files in this folder

- `si8.fdf` — an 8-atom bulk Si cell (the same fixture
  `test/2-structures/11-amorphize/` uses), small enough that the whole
  walkthrough finishes in well under a minute on CPU. See "Why an 8-atom
  cell?" above for what that trades away.
- `example_2.11.sh` — the guided walkthrough (**not** an automated test —
  see `test/2-structures/11-amorphize/test.sh` for that). Pauses between
  sections so you can read before moving on; safe to re-run any time.
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_amorphize_report.txt`.

## Running the walkthrough

```bash
cd examples/2.11-stb-amorphize
./example_2.11.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. The whole script is skipped, with a clear
message, if the optional `ml` extra isn't installed. Six self-contained
cases are generated:

| Folder                          | What it shows                                                          |
|-----------------------------------|-------------------------------------------------------------------------|
| `too-short-vs-adequate-melt/`     | `--melt-steps` vs. `--taut`: the bond-angle-std table above, measured live |
| `final-relax-comparison/`         | Why the final relax is on by default: the energy/force table above    |
| `vacuum-rejected/`                | Rejected outright on a `stb-slab`-built (100) slab (bulk-only)         |
| `save-data/`                      | `--save-data`: the `.dat`+`.gplot` pair, rendered to a real PDF if gnuplot is installed |
| `save-traj/`                      | `--save-traj`: all 3 formats (xsf/pdb/xyz), read back and verified with ASE |
| `full-report/`                    | `--save-report`, the symmetry `Fd-3m -> P-1` collapse, `references.bib` |

## Try it yourself

```bash
# Diagnose a protocol before committing to a bigger production cell
stb-amorphize -f si8.fdf --melt-steps 2000 --quench-steps 4000 \
    --save-data --save-traj --stride 20

# Watch the melt-quench in OVITO/VMD directly
stb-amorphize -f si8.fdf --save-traj --traj-format pdb

# Pre-relax with your own fine-tuned MACE model, keeping the full trace
stb-amorphize -f si8.fdf --custom-model my_finetuned.model --save-data
```

## Flag reference

| Flag                  | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-f/--file`           | Input crystalline structure file (`.fdf`, bulk only)                    |
| `--melt-temp`         | Melt temperature, K (default `3000.0`)                                  |
| `--melt-steps`        | MD steps held at `--melt-temp` (default `500`) — see `--taut` above     |
| `--quench-temp`       | Final target temperature, K (default `300.0`)                           |
| `--quench-steps`      | MD steps for the linear cooling ramp (default `1000`)                   |
| `--timestep`          | MD timestep, fs (default `1.0`)                                         |
| `--taut`/`--taup`     | Berendsen thermostat/barostat response times, fs (default `50.0`/`200.0`) |
| `--compressibility`   | Barostat compressibility, eV/Ang^3 (default `4.57e-5`, water's value — a generic placeholder) |
| `--stride`            | Sample every N MD steps for `--save-data`/`--save-traj` (default `10`)  |
| `--save-data`         | Write `<stem>_md_diagnostics.dat` + `.gplot`                            |
| `--save-traj`         | Write a multi-frame MD trajectory                                       |
| `--traj-format`       | Trajectory format for `--save-traj`: `xsf`/`pdb`/`xyz` (default `xsf`)  |
| `--no-final-relax`    | Skip the final static (position+cell) relax                            |
| `--model`             | MACE-MP-0 size: `small`/`medium`/`large` (default `small`)              |
| `--custom-model`      | Path to a custom fine-tuned `.model` file                              |
| `-sp/--symprec`       | Symmetry tolerance, Ang, for the before/after table (default `0.01`)    |
| `-o/--output`         | Output `.fdf` file name (default `amorphous.fdf`)                       |
| `--save-report`       | Persist the full report (incl. symmetry table) to disk                 |
| `--view`              | Open the input and final structure in ASE's interactive viewer          |

## What's next

Take the amorphized `.fdf` into a real SIESTA relaxation/calculation (e.g.
via `stb-inputfile`, example `1.1`) — remember this is a fast heuristic
starting guess, not a DFT-verified amorphous structure. See
`2.10-stb-molecule/` for another Structures-category tool with the same
numbered-report/provenance-header conventions, or `1.7-stb-mlrelax/` for
another MACE-driven tool in this suite (including the same FIRE/BFGS/LBFGS
convergence theory `stb-amorphize`'s own final relax step also relies on).
If you have your own SIESTA training data, `stb-mlffAnalysis` (Workflow
menu) can fine-tune a custom MACE model for `--custom-model` that's more
accurate for your specific material than the generic foundation checkpoint.
