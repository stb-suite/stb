# 3.18 — `stb-aimdAnalysis`: AIMD Trajectory Analysis

`stb-aimdAnalysis` turns a raw AIMD/MD trajectory -- a SIESTA
`<label>.ANI` (+ `.out` for the true per-step cell) or a generic
ASE-readable file (e.g. one written by `stb-mlmd`) -- into physics:
structure (RDF), transport (MSD/diffusion), vibrations (VACF-derived
VDOS), thermodynamics (energy/temperature/volume/pressure, from SIESTA's
own `.MDE` file), and per-atom displacement/atom-pair relative-distance
tracking. It complements `stb-ani2traj`, which only converts a
trajectory's file *format* for external viewers and computes nothing.

## 1. Theory

### 1.1 Radial distribution function g(r)

For every frame, the minimum-image distance between every requested pair
of atoms is histogrammed and normalized against the ideal-gas shell
density -- the standard structural fingerprint of a liquid/amorphous/gas
system (bond lengths show up as sharp peaks). `--pair A-B` restricts
this to one species pair (e.g. `O-H`); the default is the total RDF over
every pair. Correctly re-derived per frame from that frame's OWN cell, so
a variable-cell (NPT/Parrinello-Rahman) run stays correct.

### 1.2 Mean-squared displacement (MSD) and diffusion coefficient

`MSD(t) = <|r(t) - r(0)|^2>`, averaged over multiple sliding time
origins (not just `t=0`, the standard way to reduce statistical noise)
and over every atom of a given species. A linear fit of the diffusive
-regime part of the curve (default: the second half of the trajectory)
gives the diffusion coefficient via the 3D Einstein relation `D =
slope/6`, converted from Ang^2/fs to the conventional cm^2/s.

### 1.3 Vibrational density of states (VDOS)

`.ANI`/most generic trajectory formats carry only POSITIONS, not
velocities -- velocities are estimated by central finite difference,
`v(t) = (x(t+dt) - x(t-dt))/(2 dt)`. The velocity autocorrelation
function (VACF) is computed via FFT (fast, `O(N log N)`) and its
Fourier transform (Wiener-Khinchin theorem) gives the VDOS. This
amplifies positional noise into the derivative, so the result is
explicitly qualitative, not benchmark-grade -- useful for "where's the
dominant vibrational frequency", not a publication-grade phonon DOS.

### 1.4 Single-atom displacement tracking

`--track-atom N` (0-based atom index) reports ONE specific atom's own
Cartesian displacement from its initial position, frame by frame --
built on the same PBC-unwrapped trajectory the MSD already uses, so a
real diffusing/hopping atom's displacement stays physically continuous
instead of resetting every time it crosses a periodic boundary. Useful
for watching one diffusing adatom/dopant/defect directly, instead of
only the species-averaged MSD.

### 1.5 Atom-pair relative-distance tracking

`--track-pair I-J` (two 0-based atom indices, e.g. `0-5`) reports the
relative separation between two SPECIFIC atoms every frame, using the
SAME minimum-image convention as the RDF -- deliberately NOT the
unwrapped trajectory `--track-atom` uses. What matters physically for a
bond length, a hydrogen-bond distance, or a reaction coordinate is the
atoms' TRUE instantaneous separation; unwrapping the two atoms
independently would not give this in general, since each atom can
accumulate its own separate drift/unwrap path that doesn't cancel back
down to the true short separation.

### 1.6 Thermodynamic time series (energy/temperature/volume/pressure)

`[7] THERMODYNAMIC TIME SERIES` reports and plots four quantities in one
4-panel figure. Volume is always available (computed directly from each
frame's own cell, any input source). For a `--label` run, Energy/
Temperature/Pressure are read straight from SIESTA's own dedicated
`<label>.MDE` file (`core.siesta_log.get_mde_trajectory`) -- a small,
clean per-step table (`Step T(K) E_KS(eV) E_tot(eV) Vol(Ang^3)
P(kBar)`) SIESTA already writes for exactly this purpose, rather than
re-scraping it back out of scattered `.out` log lines. Unlike the INPUT
`.fdf`, `.MDE` is always named after `SystemLabel` (like `.XV`/`.ANI`),
so no `--geometry-file`-style override is needed for it.

Both `E_tot` (total energy, kinetic+potential) and `E_KS` (the
electronic/potential-like energy alone) are plotted together on the
Energy panel when both are available -- `E_tot` trades between kinetic
and potential energy by design, so only the total is meaningful to
check for stability/drift (the same distinction `stb-mlmd`'s own energy
tracking relies on). For a generic `--trajectory` input, Energy/
Temperature fall back to each frame's own embedded `Epot`/`Temp` info
if present (e.g. from `stb-mlmd --out-format xyz`); Pressure is never
available there (no known source).

### 1.7 Why energy is plotted per atom, not as an absolute total

Energy is an **extensive** quantity -- it scales with how many atoms
are in the simulation cell. `-1029.53 eV` for an 8-atom cell and
`-1029.53 eV` for a 64-atom cell of the SAME material would be wildly
different physical states (roughly 8x less bound per atom in the
second case), yet the absolute number alone doesn't tell you that.
Dividing by the atom count gives an **intensive** quantity (`eV/atom`)
-- independent of system size, and the natural axis for judging
stability/drift or comparing two different-sized calculations of the
same material. This is why the Energy panel (both the matplotlib
`--view` figure and the gnuplot `<stem>_energy.gplot`/`<stem>_thermo
.gplot` scripts) plots `E_tot`/`E_KS` in `eV/atom` -- the raw absolute
totals are still written as extra columns in `<stem>_energy.dat` and
still printed in the `[7]` report table, just not what's on the y-axis.

## 2. Cross-check: `--track-pair` against the RDF's own first peak

`aimd.ANI`/`.XV`/`.fdf`/`.out` is a real SIESTA AIMD run: a 5-step
Verlet trajectory of an isolated O2 dimer in vacuum. With only 2 atoms,
`--track-pair 0-1` IS the O-O bond -- and it independently agrees with
the RDF's own first peak, two completely separate code paths (a
statistical histogram over all frames vs. one specific pair's per-frame
minimum-image distance) landing on the same answer:

```
$ stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 --no-intro
...
First peak        : r = 1.087 Ang (g(r) = 269.148)
...
Tracked pair      : atom 0 (O) -- atom 1 (O)
Initial distance (Ang)          | 1.0864
Final distance (Ang)            | 1.2284
```

The bond visibly stretches over the 5-step trajectory (1.086 -> 1.228
Ang) -- too short a run to mean anything thermodynamically, but exactly
the kind of event `--track-pair` exists to let you watch directly
instead of only seeing it blurred into a population-averaged RDF.

## 3. Features at a glance

- **`--geometry-file <path>`**: an explicit, `--label`-decoupled path to
  the real SIESTA input `.fdf` (needed for the true MD timestep) --
  the same `--sts`/`--coop`/`--ipr`/`--effmass`/`--spintexture` pattern
  this suite already uses for `.fdf`/`.HSX` inputs. Details and a live
  reproduction of the real-world filename mismatch it fixes: section 4.
- **`--list-atoms`**: a standalone early-exit mode -- prints every
  atom's 0-based index, species, and Cartesian coordinates (from the
  first frame only, so it's fast regardless of trajectory length), then
  exits immediately without running RDF/MSD/VDOS/anything else.
  Deliberately a separate opt-in mode rather than a table always
  printed in the report -- a real structure can have hundreds of atoms.
  The interactive `stb-suite` menu asks "List every atom's index/
  species/coordinates?" (y/N) right before prompting for
  `--track-atom`/`--track-pair`.
- **`--track-atom N`** (single-atom displacement) and **`--track-pair
  I-J`** (atom-pair relative distance) -- report sections `[5]`/`[6]`,
  always printed (say "Not requested" if the flag is omitted, same
  convention as `[8]`'s gnuplot section). Theory: sections 1.4/1.5.
- **`[7] THERMODYNAMIC TIME SERIES`** -- energy/temperature/volume/
  pressure in one 4-panel figure (matplotlib via `--view`, gnuplot via
  `--save-gnuplot`), reading Energy/Temperature/Pressure from SIESTA's
  own dedicated `<label>.MDE` file for a `--label` run (section 1.6).
  The Energy panel's y-axis plots `eV/atom`, not the raw absolute total
  (an intensive, system-size-independent quantity -- section 1.7);
  `[7]`'s report table and `<stem>_energy.dat`'s extra columns still
  carry the absolute eV values too. Verified on a real 500-step SIESTA
  NVT run: section 5.
- **Numbered `[0]`...`[10]` report** (`RUN METADATA`, `INPUT DATA`,
  `RDF`, `MSD & DIFFUSION`, `VACF/VDOS`, the two tracking sections, the
  thermodynamic section, `WRITING OUTPUT FILES`, `REFERENCES`,
  `SUMMARY & FILES`).
- **`--save-gnuplot`** writes a real `.dat` + `.gplot` pair for every
  computed quantity (RDF, MSD, VACF, VDOS, volume, energy, temperature,
  pressure, and the atom-displacement/pair-distance series if used)
  plus one combined 4-panel `<stem>_thermo.gplot`.
- **`--view`** shows the same figures as an interactive matplotlib
  preview, off by default.
- **`-o`/`--output-dir`** selects where every output file is written
  (default: the current directory).
- **`[9] REFERENCES`** writes `references.bib` citing SIESTA for a
  `--label` run; correctly prints a note instead for a generic
  `--trajectory` input, since that path isn't guaranteed to be a SIESTA
  run at all (section 6).

## 4. `--geometry-file`: the real-world `.fdf` filename mismatch, in detail

`--label` mode reads two things from a `.fdf`: the true MD timestep
(`MD.InitialTimeStep`/`MD.LengthTimeStep`) and, only as a LAST-resort
lattice fallback (after `<label>.out`'s per-step cell and `<label>.XV`),
the geometry itself. SIESTA always names `.XV`/`.ANI`/`.HSX`/`.WFSX`
after `SystemLabel`, but the INPUT `.fdf`'s own filename is chosen by
the user and is very often different (e.g. `SystemLabel siesta` with
the real input called `calc.fdf`) -- without an explicit override, this
silently degrades to "assume 1 fs per frame" with only an easy-to-miss
warning. Renaming this example's own `aimd.fdf` to `calc.fdf`
reproduces the scenario and its fix directly:

```
$ stb-aimdAnalysis --label aimd --no-intro   # no aimd.fdf present -- only calc.fdf
[WARNING] MD.LengthTimeStep not found in 'aimd.fdf (auto-detected)' --
assuming 1 fs per (strided) frame for MSD/VDOS time axes. If the real
input file isn't named 'aimd.fdf', pass --geometry-file <path>.

$ stb-aimdAnalysis --label aimd --geometry-file calc.fdf --no-intro
Geometry file     : calc.fdf
Timestep (post-stride) | 1.0000 fs
```

`--geometry-file` only applies to `--label` mode (a `--trajectory` input
already carries its own per-frame cell, no `.fdf` involved at all) -- and
is validated to exist up front, same as every other explicit-path flag
in this suite.

## 5. A real 500-step SIESTA NVT run: `.MDE` + `--geometry-file` together

`siesta.ANI`/`.XV`/`.MDE` + `calc.fdf`/`structure.fdf` are a SECOND,
real SIESTA run bundled with this example (distinct from the small
5-step `aimd.*` fixture above) -- 500 MD steps, an 8-atom SiC supercell,
Nose thermostat at a 500 K target. `SystemLabel siesta` but the real
input is `calc.fdf` -- reproducing the exact `--geometry-file` scenario
above on a richer run, AND giving `[7] THERMODYNAMIC TIME SERIES` real
`.MDE` data to read:

```
$ stb-aimdAnalysis --label siesta --geometry-file calc.fdf --save-gnuplot --no-intro
[7] THERMODYNAMIC TIME SERIES (ENERGY / TEMPERATURE / VOLUME / PRESSURE)
Quantity                    | Mean        | Std dev  | Min         | Max
Volume (Ang^3)              | 691.4382    | 0.0000   | 691.4382    | 691.4382
Energy, total (eV)          | -1029.5308  | 0.0029   | -1029.5339  | -1029.5154
Energy, potential (eV)      | -1030.4177  | 0.4158   | -1031.0612  | -1029.5250
Temperature (K)             | 503.0835    | 251.5614 | 137.9600    | 1250.3500
Pressure (kBar)             | -4.7330     | 13.6222  | -23.0030    | 29.1470
Energy, total (eV/atom)     | -128.691355 | 0.000363 | -128.691739 | -128.689429
Energy, potential (eV/atom) | -128.802211 | 0.051971 | -128.882655 | -128.690628
```

`-1029.5308 eV / 8 atoms = -128.691355 eV/atom` -- the per-atom rows are
exactly what `siesta_energy.dat`'s plotted columns (2/3) and
`siesta_energy.gplot`/`siesta_thermo.gplot`'s y-axis carry (section 1.7
above), while these absolute-eV rows and their own extra `.dat` columns
stay available for reference.

Three independent physical sanity checks all pass on this real data:
- **Volume is exactly constant** (std = 0.0000 Ang^3) -- correct for a
  fixed-cell NVT run (Nose thermostat controls temperature, not cell
  shape/volume).
- **Temperature fluctuates around the 500 K Nose target** (mean 503.1 K)
  -- expected for a small 8-atom system's canonical-ensemble
  fluctuations, not clamped the way a barostat/thermostat's instantaneous
  value never is.
- **`E_tot` is ~140x more stable than `E_KS`** (std 0.0029 eV vs.
  0.4158 eV) -- confirming `E_tot`, not `E_KS`, is the physically
  appropriate "conserved-ish" quantity to judge stability by (section
  1.6 above).

## 6. `--trajectory`: independent of SIESTA

Not every trajectory this tool can read comes from SIESTA -- `stb-mlmd`
(a MACE-driven MD run) writes the same xsf/pdb/xyz formats
`stb-aimdAnalysis --trajectory` accepts. Because that input isn't
guaranteed to be a SIESTA calculation, `[9] REFERENCES` correctly skips
the SIESTA citation for it instead of assuming one:

```
$ stb-aimdAnalysis --trajectory synthetic_traj.xyz --track-pair 0-1 --no-intro
[9] REFERENCES
No SIESTA-specific references for a generic --trajectory input -- cite
whichever tool produced it instead (e.g. stb-mlmd's own MACE/foundation
-model references).
```

An extended-xyz trajectory written with per-frame `Time` info (e.g.
`stb-mlmd --out-format xyz`) auto-detects `dt`; `xsf`/`pdb` carry no such
data, so `--dt` is required for those.

## 7. Known, deliberate limitations

- The VDOS is qualitative, not benchmark-grade (section 1.3).
- This example's small `aimd.*` fixture (5 frames) is far too short for
  the MSD/diffusion or VDOS numbers to carry real physical meaning --
  it only demonstrates that the pipeline runs correctly end to end.
- `--track-pair`'s minimum-image convention assumes no atom moves more
  than half a cell length between frames (the same assumption the RDF
  and PBC-unwrapping already make elsewhere in this suite).

## 8. Files in this folder

| File | What it is |
|------|------------|
| `aimd.ANI`/`.XV`/`.fdf`/`.out` | The real 5-step SIESTA AIMD run (O2 dimer in vacuum) |
| `siesta.ANI`/`.XV`/`.MDE`, `calc.fdf`, `structure.fdf` | The real 500-step SIESTA NVT run (8-atom SiC supercell, Nose thermostat) -- `SystemLabel siesta`, real input `calc.fdf` |
| `example_3.18.sh` | This walkthrough's runnable script |

## 9. Running it

```bash
bash example_3.18.sh
```

Or try the tool directly:

```bash
stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 \
    --save-report --save-gnuplot --view
```

## 10. What's next

- `stb-ani2traj` (Utils) -- converts the same `.ANI` trajectory's
  *format* for external viewers (OVITO/VMD), computing nothing; run it
  alongside this tool to actually watch the motion this tool quantifies.
- `stb-mlmd` (ML Simulations) -- a MACE-driven MD run whose own
  xsf/pdb/xyz output is exactly what `--trajectory` reads in section 6
  above.
