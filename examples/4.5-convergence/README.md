# 4.5 — Workflow: Convergence Tests (`stb-convergence` / `stb-convergenceAnalysis`)

Every SIESTA calculation depends on three numerical "dials" that have no
universally correct setting: how fine the real-space integration grid is
(`Mesh.CutOff`), how spatially complete the localized basis is
(`PAO.EnergyShift`), and how densely the Brillouin zone is sampled
(the k-grid). This workflow has **2 stages**: **Stage 1**
(`stb-convergence`, code `4.5.1`) sweeps one or more of those three
parameters, writing one **fully relaxed** (positions + cell) SIESTA folder
per swept value. **Stage 2** (`stb-convergenceAnalysis`, code `4.5.2`) reads
the finished runs back and reports, per parameter, the value at which
**both** the total energy **and** the relaxed geometry itself stop changing
beyond a chosen tolerance — not just one or the other.

Both stages live in this one folder and this one tutorial — not one folder
per tool — for the same reason `4.1-strain`/`4.2-elastic`/`4.3-cohesive` do:
Stage 1's output only exists to feed Stage 2, and Stage 2 only exists to
interpret Stage 1's output. This workflow is also the natural **first**
stop before any of those others: its entire purpose is making sure the
`calc.fdf` you're about to use for a real strain/elastic/cohesive/phonon
calculation isn't quietly wrong in a way that has nothing to do with
physics.

## 1. Theory

### 1.1 What "convergence testing" means, and why every calculation needs one first

SIESTA (like every real DFT code) never solves the exact many-electron
problem — it approximates it using a **finite basis** of localized atomic
orbitals, a **finite real-space grid** for the electron density, and a
**finite sampling** of the Brillouin zone. Each of these approximations has
a tunable "fineness" parameter, and — critically — **the computed total
energy, forces, and relaxed geometry all depend on the chosen fineness, not
just on the physics you're trying to compute.** A calculation run at a
too-coarse setting doesn't fail outright; it silently returns a number that
looks like an answer but contains numerical error that can dwarf the real
physical effect you're after (e.g. a strain-energy difference of a few meV,
Section 4.1/4.2's whole subject). **Convergence testing** is the systematic
fix: sweep a parameter, watch where the answer stops changing beyond some
tolerance, and use *that* value — not a "reasonable-sounding" textbook
default — for every production calculation downstream.

### 1.2 The three dials this tool turns, and the physics behind each

| Parameter | fdf keyword | Unit | What it actually controls | Direction that improves accuracy |
|---|---|---|---|---|
| `meshcutoff` | `Mesh.CutOff` | Ry | Fineness of the real-space grid SIESTA uses to compute the Hartree/XC energy and potential — the LCAO analogue of a plane-wave cutoff | **larger** |
| `energyshift` | `PAO.EnergyShift` | Ry | The energy increase, relative to the free orbital, at which each basis orbital's radial tail is truncated — smaller means the orbital is allowed to extend farther out, making the basis more spatially complete | **smaller** |
| `kgrid` | `kgrid.MonkhorstPack` (density-driven) | 1/Ang | Brillouin-zone sampling density — matters most for metals/small-gap systems where bands vary rapidly near the Fermi level, far less for large-gap insulators or molecules in a big vacuum box | **smaller** density value (= a denser grid) |

`--parameter` selects one, several, or `all` three; each is swept
**independently** in its own subfolder (Section 1.6). `kgrid` is swept as a
**density**, in 1/Ang — the same convention `stb-kgrid` itself uses — not as
a raw `N x N x N` grid directly; Stage 1 converts it via
`core/kspace.compute_monkhorts` (Section 2.5 shows this live).

### 1.3 Why every sweep point is a *full relaxation*, not a fixed-geometry single point — the eggbox effect

A finite real-space grid is not perfectly translationally invariant: moving
an atom by a fraction of a grid spacing can change the computed energy by a
small, spurious amount, purely because of where the atom sits *relative to
the grid points* — SIESTA's own documentation calls this the **"eggbox
effect."** This matters enormously for convergence testing specifically:
if you compared energies at a single **fixed, arbitrary geometry** across a
`Mesh.CutOff` sweep, two different settings could both look
self-consistent and "energy-converged" simply because they're evaluating
the same (possibly not-yet-relaxed) geometry — while the geometry SIESTA
would actually *relax to* at each setting is still quietly drifting. A
fixed-point comparison cannot see that drift at all.

This is exactly why Stage 1 forces `MD.TypeOfRun CG` + `MD.Steps` +
`MD.VariableCell true` into **every** generated folder, regardless of what
your own `--calc` template says (verified live in Section 2.2's captured
`[3] RELAXATION & PARAMETER ENFORCEMENT` table) — full positions-and-cell
relaxation at every single swept value. It's more expensive than a
single-point sweep, but it's the only honest way to test this, and it sets
up the entire reason Stage 2 needs *two* separate convergence criteria
(Section 1.4).

### 1.4 Two independent convergence criteria: energy AND relaxed structure

Because Section 1.3's eggbox effect can make total energy *look* converged
at a geometry that itself hasn't settled, Stage 2 tracks **two criteria
separately, per parameter**:

- **Energy**: `|Delta E/atom|` between consecutive swept values, against
  `--tolerance` (default `0.01` eV/atom).
- **Structure**: `|Delta V/V|`, the percent change in the *relaxed* cell
  volume between consecutive values, against `--volume-tolerance` (default
  `0.5` %).

Both are reported side by side in `[2.N]`'s table, and `[4]`'s summary
shows exactly where each one first stabilizes. Section 4's worked example
shows a real, code-verified case where **they disagree** — energy looking
converged two sweep steps before the relaxed geometry actually does. This
is the entire reason Stage 1 always does a full relaxation (Section 1.3)
instead of a cheaper single-point sweep: a single-criterion, fixed-geometry
tool structurally *cannot* catch this.

### 1.5 Why the plateau scan checks the *entire remaining tail*, not just the next point

A single small delta between two adjacent swept values could be a
coincidence — two points landing close together by chance, or ordinary SCF
noise — not a genuine plateau. Stage 2's internal scan (`_scan_plateau`)
only declares a value "converged" if **every later point in the sweep**
also stays under tolerance (and has confirmed SCF convergence), not just
the very next one. This avoids a false "converged" reading from one lucky
adjacent pair while the underlying trend hasn't actually settled yet — the
same discipline `4.3-cohesive`'s own BSSE-cutoff convergence check applies
(one before/after point is a hint, not proof; a scan is proof).

### 1.6 Independent, not joint, sweeps — a real limitation worth understanding before it bites

Each of the three parameters is swept **alone**, holding the other two
fixed at whatever your `--calc` template already says. The tool never
performs a joint/coupled optimization across all three at once. This
matters in practice whenever two parameters are physically coupled: for
example, a very coarse k-grid in a metallic system can itself introduce
enough energy noise that a `Mesh.CutOff` scan run *at that same coarse
k-grid* looks "converged" earlier than it really would at a denser,
properly-converged grid — the k-grid-driven noise floor can simply be
larger than the mesh-cutoff-driven changes you're trying to measure.
**Practical recommendation**: converge `kgrid` first, then re-verify
`meshcutoff`/`energyshift` with that denser k-grid now held fixed in your
`--calc` template, rather than trusting one pass's independent sweeps in
isolation for a system where this coupling is likely to matter (metals,
small cells).

### 1.7 Data-quality diagnostics are not a convergence criterion

`--force-tolerance` (default `0.05` eV/Ang) and the `[2.N]` table's
**Quality** column (`SCF?`, `F>tol`, `steps`) flag a *specific run's own*
numerical health — did SCF actually confirm converged, is the residual
force still large (meaning the geometry might not really be at a
relaxed minimum), did the relaxation exhaust its declared `MD.Steps`
budget (meaning it may have been cut off before finishing). This is
**completely separate machinery** from the tolerance-based plateau scan of
Section 1.4/1.5 — a flagged row never blocks or alters the convergence
scan itself, it's advisory only. Worth stating plainly because the two are
easy to conflate: `--tolerance`/`--volume-tolerance` ask "has the *trend*
settled," `--force-tolerance` asks "is *this one run* actually trustworthy
in the first place." Section 3.5 demonstrates a real flagged row.

## 2. Stage 1: generating the calculations (`stb-convergence`, code `4.5.1`)

### 2.1 What it does

Reads one structure and one `calc.fdf` template, then writes, under one
wrapping `-o/--output-dir` (default `convergence_runs/`):

```
<output-dir>/<parameter>/convergence_<parameter>_<value>/
    structure.fdf       (copied verbatim from -s)
    calc.fdf             (your -c template, with '%include config_extra.fdf' prepended)
    config_extra.fdf      (forces this folder's swept value + full relaxation)
    <pseudopotential files>   (only if -p/--pseudo-dir given)
```

One such subtree per selected parameter (Section 2.4), each an
**independent** sweep (Section 1.6).

### 2.2 What gets forced into every folder, and why `stb-convergence` needs a real `calc.fdf`

Unlike `4.3-cohesive` (which builds its own internal single-point template),
`stb-convergence` takes your **own** `-c/--calc` template and reuses it
almost unchanged — because the entire point is testing the numerical
settings of the calculation you're actually about to run, not some fixed
internal reference. Only `config_extra.fdf`, prepended via
`%include` (fdf is first-occurrence-wins, so this always overrides whatever
your template says), forces the swept parameter's value plus full
relaxation. Captured live, real tool output (`Mesh.CutOff = 250 Ry` case):

```
# Auto-generated by stb-convergence -- forces this folder's own swept parameter
# value plus full relaxation (positions + cell), regardless of --calc's own settings.
Mesh.CutOff           250.0000  Ry

MD.TypeOfRun       CG
MD.Steps           100
MD.VariableCell    true
```

`--relax-steps` (default `100`) controls the forced `MD.Steps`.

### 2.3 Single-parameter default sweep, live

```
$ stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter meshcutoff --no-intro

[2] PARAMETER SELECTION
------------------------------------------------------------
Parameter  | Range source      | Sweep range              | Convergence direction
-----------------------------------------------------------------------------------
meshcutoff | suggested default | 100 to 400, step 50 (Ry) | larger = more converged

[5] GENERATED CONVERGENCE FOLDERS
------------------------------------------------------------
Under 'convergence_runs/<parameter>/':
Folder                          | Parameter  | Value       | Files
------------------------------------------------------------------------------------------------------
convergence_meshcutoff_100.0000 | meshcutoff | 100.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
convergence_meshcutoff_150.0000 | meshcutoff | 150.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
convergence_meshcutoff_200.0000 | meshcutoff | 200.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
convergence_meshcutoff_250.0000 | meshcutoff | 250.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
convergence_meshcutoff_300.0000 | meshcutoff | 300.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
convergence_meshcutoff_350.0000 | meshcutoff | 350.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
convergence_meshcutoff_400.0000 | meshcutoff | 400.0000 Ry | structure.fdf, calc.fdf, config_extra.fdf
```

Real, verbatim output from this folder's own `structure.fdf`/`calc.fdf`
(bulk silicon, diamond cubic, DZP/GGA-PBE — Section 7). Case 1 of
`example_4.5.sh` reproduces this exactly.

### 2.4 `--parameter all`: three independent sweeps in one call

```
$ stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter all --no-intro

[2] PARAMETER SELECTION
------------------------------------------------------------
Parameter   | Range source      | Sweep range                    | Convergence direction
-------------------------------------------------------------------------------------------
meshcutoff  | suggested default | 100 to 400, step 50 (Ry)       | larger = more converged
energyshift | suggested default | 0.001 to 0.05, step 0.01 (Ry)  | smaller = more converged
kgrid       | suggested default | 0.05 to 0.3, step 0.05 (1/Ang) | smaller = more converged
```

producing 3 independent subtrees on disk:

```
convergence_runs/meshcutoff/convergence_meshcutoff_{100,150,200,250,300,350,400}.0000/
convergence_runs/energyshift/convergence_energyshift_{0.0010,0.0110,0.0210,0.0310,0.0410,0.0510}/
convergence_runs/kgrid/convergence_kgrid_{0.0500,0.1000,0.1500,0.2000,0.2500,0.3000}/
```

Case 2 of `example_4.5.sh` runs this live.

### 2.5 Custom ranges, and the k-grid density convention

`--meshcutoff-range`/`--energyshift-range`/`--kgrid-range MIN MAX STEP`
override the suggested default range for that parameter only (harmless to
pass for a parameter not selected in `--parameter` — Section 2.6 shows what
happens if you do, though: it's actually an error, not a silent no-op).
Remember `kgrid` sweeps a **density** (1/Ang), not a raw grid — real,
verbatim `config_extra.fdf` at density `0.05` on this folder's bulk-Si
structure:

```
# Auto-generated by stb-convergence -- forces this folder's own swept parameter
# value plus full relaxation (positions + cell), regardless of --calc's own settings.
kgrid.MonkhorstPack   [24  24  24]

MD.TypeOfRun       CG
```

`0.05` (1/Ang, a *dense* sampling — remember: smaller density number =
denser grid, Section 1.2) resolved to a real `24x24x24` Monkhorst-Pack
grid for this 2-atom, `5.43` Ang cubic cell, via the exact same
`core/kspace.compute_monkhorts` machinery `stb-kgrid` itself uses.

### 2.6 Validation: what `stb-convergence` correctly refuses

Real, verbatim error messages (all exit with a non-zero status, same
`argparse`-style convention as every other tool in the suite):

```
$ stb-convergence -s structure.fdf -c calc.fdf --parameter bogus --no-intro
stb-convergence: error: --parameter: unrecognized value 'bogus' (choose from meshcutoff, energyshift, kgrid, or 'all').

$ stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --energyshift-range 0.01 0.05 0.01 --no-intro
stb-convergence: error: --energyshift-range was given but 'energyshift' is not in --parameter.

$ stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --meshcutoff-range 100 400 0 --no-intro
stb-convergence: error: --meshcutoff-range: step must be > 0.
```

Case 3 of `example_4.5.sh` reproduces all three.

### 2.7 Output layout, and running it both ways

| Path | What it is |
|---|---|
| `convergence_runs/<parameter>/convergence_<parameter>_<value>/` | One fully-relaxed SIESTA folder per swept value |
| `convergence_stage1.txt` | Persisted report, only with `--save-report` |

Both CLI and interactive menu (`stb-suite` -> `4.5.1`) produce **byte
-identical** output for identical answers — verified directly (`diff -rq`,
real check, not assumed) and reproduced by Case 4 of `example_4.5.sh`.

## 3. Stage 2: analyzing the results (`stb-convergenceAnalysis`, code `4.5.2`)

### 3.1 Report structure

```
[0] RUN METADATA
[1] DISCOVERED PARAMETER SWEEPS
[2.N] <PARAMETER> SWEEP        <- one per discovered parameter, e.g. [2.1] MESHCUTOFF SWEEP
[3] OUTPUT FILES
[4] SUMMARY & NEXT STEPS
[5] APPLY                       <- only printed with --apply
```

### 3.2 `--dir` auto-detects both folder layouts — including a real, verified quirk in Stage 1's own message

`--dir` (default `convergence_runs`, matching Stage 1's own
`-o/--output-dir` default) auto-detects **both** the nested layout
(`<dir>/<parameter>/convergence_<parameter>_<value>/`, what Stage 1
actually writes) and a flat layout (`convergence_<parameter>_<value>/`
directly under `<dir>`, e.g. this folder's own test fixtures). **Worth
flagging explicitly**: Stage 1's own `[6] SUMMARY & NEXT STEPS` prints a
`[NOTE]` claiming *"stb-convergenceAnalysis (Stage 2) doesn't yet support
this per-parameter layout -- point --dir at one parameter's own
subfolder"* — but running `stb-convergenceAnalysis --dir convergence_runs`
directly against Stage 1's own nested output (verified live, not assumed)
discovers the parameter correctly regardless:

```
[1] DISCOVERED PARAMETER SWEEPS
------------------------------------------------------------
Parameter  | Folders found
--------------------------
meshcutoff | 7
```

This is a **stale message left in Stage 1's own report text**, not a real
functional limitation — Stage 2's actual discovery code supports both
layouts already. Point `--dir` at whichever level is convenient; both work.

### 3.3 Reading the `[2.N]` table

Real, verbatim output (Section 4's worked example):

```
[2.1] MESHCUTOFF SWEEP
------------------------------------------------------------
Folders found: 4 (skipped: 0 unreadable/incomplete)  |  Increasing meshcutoff improves accuracy: True
Value    | E/atom(eV)  | |Delta E/atom| | Volume(Ang^3) | |Delta V/V|(%) | Quality
----------------------------------------------------------------------------------
150.0000 | -250.450000 | --             | 153.1304      | --             | OK
200.0000 | -250.425000 | 0.025000       | 159.2201      | 3.9768         | OK
250.0000 | -250.422500 | 0.002500       | 161.8786      | 1.6697         | OK
300.0000 | -250.422000 | 0.000500       | 161.9677      | 0.0551         | OK
```

`E/atom` and `Volume` are absolute values at that swept value; `|Delta ...|`
columns are the **change from the previous row** (Section 1.4/1.5's
plateau-scan inputs) — blank (`--`) on the first row since there's no
previous point to compare against.

### 3.4 How the "recommended value" is chosen

`[4]`'s summary table reports **both** criteria's own converged value side
by side, plus a single **Recommended** value — whichever of the two is
more accurate-ward (i.e. the *tighter* of the two), so both criteria end up
satisfied simultaneously:

```
[4] SUMMARY & NEXT STEPS
------------------------------------------------------------
Parameter  | Energy-converged | Structure-converged | Recommended
-----------------------------------------------------------------
meshcutoff | 250.0000         | 300.0000            | 300.0000
```

When the two disagree, an explicit `[NOTE]` is printed (Section 4 shows it
verbatim) — this is exactly Section 1.4's "energy can look converged before
the relaxed geometry does" scenario, made concrete.

### 3.5 The Quality column, live: `SCF?`, `F>tol`, `steps`

Real, verbatim output — one deliberately poor 50 Ry run added to the same
sweep (Section 1.7):

```
Value    | E/atom(eV)  | |Delta E/atom| | Volume(Ang^3) | |Delta V/V|(%) | Quality
-------------------------------------------------------------------------------------
50.0000  | -249.100000 | --             | 132.6510      | --             | SCF?,F>tol
150.0000 | -250.450000 | 1.350000       | 153.1304      | 15.4385        | OK
...
[WARNING] 1/5 run(s) flagged (SCF? = SCF convergence not confirmed, F>tol = residual force above --force-tolerance 0.05 eV/Ang, steps = relaxation used its full declared step budget) -- their geometry may not be a genuine relaxed minimum.
```

`SCF?` means the SCF cycle never printed a confirmed-converged line;
`F>tol` means the residual force is still above `--force-tolerance`
(default `0.05` eV/Ang) — both flag that `50.0000`'s own row is
untrustworthy on its own terms, entirely separately from whether it fits
the overall convergence trend (Section 1.7). Case 6 of `example_4.5.sh`
reproduces this.

### 3.6 `--apply`: writing the recommended value back into your own `calc.fdf`

```
$ stb-convergenceAnalysis --dir . --no-intro --apply calc.fdf

[5] APPLY
------------------------------------------------------------
[Applied] meshcutoff = 300.0000 (300.0000 Ry) -> calc.fdf
```

verified: the target file's own `Mesh.CutOff` line is rewritten in place
(`Mesh.CutOff           300.0000  Ry`) — same tag-substitution machinery
Stage 1 itself uses internally, applied here to *your* production
`calc.fdf` instead of a generated one. A parameter with no recommended
value (e.g. still not converged anywhere in the scanned range) is skipped,
not guessed.

### 3.7 The plots: `--save-gnuplot` / `--view`

`--save-gnuplot` writes, per discovered parameter, a `.dat` + companion
`.gplot` pair (same convention as `elastic_analysis.py`/`strain_analysis.py`
— see `CLAUDE.md`): `<parameter>_convergence.dat` /
`<parameter>_convergence.gplot`, a 2-panel `pdfcairo` plot (energy per atom
on top, relaxed cell volume on the bottom — the exact two quantities
Section 1.4's dual criterion tracks). `--view` additionally pops the same
2 panels up interactively via matplotlib (a no-op under a non-interactive
backend, e.g. this example's own `MPLBACKEND=Agg`). Neither is written
unless requested.

### 3.8 Running it both ways

Both CLI and `stb-suite` -> `4.5.2` produce the identical report for
identical answers.

## 4. Worked example: energy convergence lies about geometry convergence

Bulk silicon (diamond cubic, this folder's own `structure.fdf`/`calc.fdf`
— DZP basis, GGA-PBE), a `Mesh.CutOff` sweep at 150/200/250/300 Ry. **An
honest note on the data**: these are small, hand-built `calc.out` files
(the same convention this suite's own test fixtures use,
`test/4-workflow/5-convergence/analysis/`) with realistic, physically
-plausible bulk-silicon-like numbers — not a completed real ab initio
calculation. Every number and table shown below, however, is **genuine,
verbatim `stb-convergenceAnalysis` output**, generated by actually running
the real tool against them — nothing here is hand-derived or approximated.

```
Value    | E/atom(eV)  | |Delta E/atom| | Volume(Ang^3) | |Delta V/V|(%) | Quality
----------------------------------------------------------------------------------
150.0000 | -250.450000 | --             | 153.1304      | --             | OK
200.0000 | -250.425000 | 0.025000       | 159.2201      | 3.9768         | OK
250.0000 | -250.422500 | 0.002500       | 161.8786      | 1.6697         | OK
300.0000 | -250.422000 | 0.000500       | 161.9677      | 0.0551         | OK

Energy converged at: meshcutoff = 250.0000 (or higher, |Delta E/atom| < 0.01 eV)
Structure converged at: meshcutoff = 300.0000 (or higher, |Delta V/V| < 0.5%, relaxed cell)
[NOTE] Energy and relaxed structure disagree on where meshcutoff converges -- using the more accurate-ward of the two (300.0000) so BOTH criteria are satisfied. This is exactly the scenario stb-convergence's own docstring warns about: the total energy can look converged before the relaxed geometry actually is.
```

Reading this the way Section 1.4 set up:

- **Energy alone would say 200 Ry is already "close enough"** — the delta
  from 200->250 Ry is only `0.0025` eV/atom, well under the `0.01`
  eV/atom default tolerance, and it only gets smaller from there
  (`0.0005` eV/atom by 300 Ry). If this tool only tracked energy, you'd
  walk away thinking `200` Ry (or even less) was a safe production
  setting.
- **The relaxed cell volume tells a different story.** The cube side grows
  `5.35 -> 5.42 -> 5.45 -> 5.451` Ang across the same 4 points — `Delta V/V`
  drops `3.98% -> 1.67% -> 0.055%`. It only falls below the `0.5%` default
  tolerance on the **last** step, `250 -> 300` Ry. At `200` or even `250`
  Ry, the relaxed geometry itself is still measurably drifting even though
  the energy already looks flat.
- **The tool's own recommendation, `300.0000` Ry, is the honest answer** —
  the smallest value where *both* criteria are simultaneously satisfied,
  not just the cheaper one two criteria would have picked in isolation.
  Using `200` Ry for a real strain/elastic-constant calculation downstream
  (`4.1-strain`/`4.2-elastic`) would risk comparing relaxed geometries that
  are still shifting with the mesh — numerical noise that could easily be
  mistaken for a real physical trend at the sub-percent strain scale those
  workflows care about.

This is precisely the failure mode Section 1.3's eggbox-effect discussion
warns about, made concrete with real tool output rather than left
abstract: **a coarser Mesh.CutOff can absolutely look energy-converged
while the geometry it would relax atoms to is still moving.**

## 5. Known, deliberate limitations

- **Independent, not joint, sweeps** (Section 1.6) — each parameter is
  varied alone; the true joint optimum across all three can differ,
  especially when two parameters are physically coupled (a coarse k-grid
  masking a real mesh-cutoff dependence in a metal, for example).
- **Default ranges are a fixed, evenly-spaced approximation** of standard
  SIESTA convergence-test scales, not the literal non-uniform point sets
  some tutorials quote, and are **not derived from your actual structure**
  — heavier elements or higher-accuracy work may need a wider
  `Mesh.CutOff` range than the default `100-400` Ry.
- **Full variable-cell relaxation at every single swept point, always**
  (Section 1.3) — the only honest way to catch the eggbox effect, but
  correspondingly more expensive per point than a fixed-geometry
  single-point sweep.
- **`--force-tolerance`/the Quality column are advisory only** (Section
  1.7) — they never gate or alter the tolerance-based plateau scan itself.
- **`--apply` for `kgrid` re-derives the Monkhorst-Pack grid from only the
  *first* discovered folder's own `structure.fdf`** for that parameter,
  assuming the same structure/vacuum layout holds across the whole sweep.
- **Stage 1's own `[6]` NOTE about Stage 2's layout support is stale**
  (Section 3.2, verified live) — a minor, real documentation inconsistency
  in the tool's own printed text, not a functional limitation of Stage 2
  itself.

## 6. Step-by-step: running this workflow on your own structure

0. **Prerequisite**: a `structure.fdf` **and** a real `calc.fdf` — the
   production settings you actually intend to use. Unlike `4.3-cohesive`,
   Stage 1 does **not** build its own template (Section 2.2); it tests
   *your* settings.
1. **Decide which parameter(s) to test** (Section 1.2): `meshcutoff`,
   `energyshift`, `kgrid`, or `all` three.
2. **Run Stage 1** (`stb-convergence`, code `4.5.1`) — CLI or `stb-suite`
   -> `4.5.1`. Writes everything under `convergence_runs/` (or your own
   `-o`).
3. **Run SIESTA yourself** in every generated
   `convergence_runs/<parameter>/convergence_*/` folder.
4. **Run Stage 2** (`stb-convergenceAnalysis`, code `4.5.2`) once they
   finish — CLI or `stb-suite` -> `4.5.2`.
5. **Read `[2.N]`'s Quality column first** (Section 3.5) — don't trust a
   flagged row's number until you know why it's flagged.
6. **Read `[4]`'s Energy-converged / Structure-converged / Recommended**
   columns (Section 3.4) — if the first two disagree, trust Recommended,
   not whichever one you happened to look at first.
7. **(Optional)**: `--apply` to write the recommended value straight into
   your own production `calc.fdf` (Section 3.6); `--save-gnuplot`/`--view`
   for the 2-panel energy/volume plot (Section 3.7).

## 7. Files in this folder

| File | What it is |
|---|---|
| `structure.fdf` | Bulk silicon, diamond cubic, 2 atoms/cell — the standard SIESTA convergence-test fixture (same file `test/4-workflow/5-convergence/` itself uses) |
| `calc.fdf` | DZP/GGA-PBE template — the actual settings this example convergence-tests |
| `example_4.5.sh` | This walkthrough's runnable script (both stages, 7 cases) |

## 8. Running the script

```bash
bash example_4.5.sh
```

## 9. What's next

- Once a parameter is converged here and applied (`--apply`) into your
  `calc.fdf`, feed that **same** `calc.fdf` into `4.1-strain`,
  `4.2-elastic`, `4.3-cohesive`, or `4.4-phonons` — this workflow's entire
  purpose is making sure those results aren't numerical-parameter noise in
  disguise.
- **`stb-mlconvergence`** (5-ML Simulations menu) sounds similar but tests a
  completely different axis: MACE-MP-0 **model size** (small/medium/large),
  never any DFT numerical parameter, and never compared against real DFT at
  all — don't confuse the two.
