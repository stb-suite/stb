# 4.4 — Workflow: Phonons (`stb-phononsCreate` / `stb-phononsPos`)

This workflow has 2 stages: **Stage 1** (`stb-phononsCreate`, code `4.4.1`)
takes an already-relaxed structure and its own `calc.fdf`, displaces each
symmetry-inequivalent atom by a tiny amount, and writes one ready-to-run,
single-point SIESTA folder per displacement. **Stage 2** (`stb-phononsPos`,
code `4.4.2`) reads the finished SIESTA runs (or an ML-computed force-constant
file) back and reports whether the structure is dynamically stable, its
phonon band structure, density of states, and the thermal properties
(free energy, entropy, heat capacity) derived from it. Together they answer
a question neither `4.1-strain`, `4.2-elastic`, nor `4.3-cohesive` does: not
how the structure *responds* to a deformation or how strongly it's bound,
but **how it vibrates**, and whether the relaxed geometry you started from
is actually a genuine energy minimum at all (a structure can be perfectly
converged in energy/forces and still be dynamically unstable).

Both stages live in this one folder and this one tutorial — same reasoning
as `4.1`/`4.2`/`4.3`: Stage 1's output only exists to feed Stage 2.

## 1. Theory

### 1.1 The harmonic approximation and the finite-displacement method

Near its equilibrium geometry, a crystal's potential energy surface is
approximated as a quadratic (harmonic) function of atomic displacements —
the same approximation behind a classical mass-and-spring lattice model.
The curvature of that quadratic form, the **force-constant matrix**
`Φ_ij = -dF_i/du_j` (how much the force on atom `i` changes per unit
displacement of atom `j`), is exactly what this workflow computes, via
finite differences: displace one atom a small, known amount `u`
(`--distance`, default `0.01` Ang), compute the resulting forces on every
*other* atom with a single-point SCF, and `Φ ≈ -F/u`. Diagonalizing the
resulting dynamical matrix (the mass-weighted Fourier transform of `Φ`) at
each q-point gives the phonon frequencies `ω(q)` — the actual physical
output of this workflow. This is a genuinely different kind of calculation
from `4.1-strain`/`4.2-elastic` (both static, zero-frequency responses to a
finite deformation) — phonons are the *dynamic*, frequency-dependent
response of the lattice.

### 1.2 Symmetry-reduced displacements: why AlP only needs 2 of 96

A brute-force calculation would displace every atom in the supercell along
every Cartesian direction, both `+` and `-` (Phonopy's own default when no
symmetry is assumed, since a real force response isn't perfectly linear)
— `N_atoms * 3 * 2` displacements: `16 * 3 * 2 = 96` for AlP's
2x2x2 = 16-atom supercell, exactly the "of 96" Section 3.3's live report
quotes. But if the crystal has real symmetry
(zincblende AlP: space group `F-43m`, 192 symmetry operations), most of
those displacements are related to each other by a symmetry operation and
give *redundant* information — Phonopy (via spglib) detects this and
generates only the symmetry-inequivalent subset, reconstructing the full
force-constant tensor afterward by applying the detected symmetry
operations to the computed subset. For AlP this is a **97.9% reduction**:
2 real SIESTA calculations instead of 96 (Section 3.3, live).

### 1.3 From force constants to a band structure: the q-path

A "band structure" plot needs a specific path through reciprocal space
connecting high-symmetry points (Γ, X, W, ...) appropriate to the crystal's
own Bravais lattice. `stb-phononsPos` auto-detects this using **ASE's own
Bravais-lattice machinery** (`Cell.get_bravais_lattice`/`Cell.bandpath`,
the same approach `stb-kpath` already uses for electronic band structures)
— deliberately **not** Phonopy's own `seekpath`-based `auto_band_structure`,
since `seekpath` isn't a dependency of this suite at all (Section 2). The
path itself (e.g. AlP's `Γ-X-W-K-Γ-L-U-W-L-K | U-X`) is standard for an
FCC-derived lattice and is printed in every Stage 2 report that uses
`--bands`.

### 1.4 Dynamical stability: what a negative frequency means, and why a hard zero is the wrong test

Because `ω(q)` comes from diagonalizing the dynamical matrix, its square
`ω²` is the actual eigenvalue — a **negative** eigenvalue gives an
**imaginary** frequency, conventionally *displayed* as a negative number
(e.g. `-5.8 THz`) even though "negative frequency" isn't physically
meaningful on its own. What it means: the harmonic energy actually
*decreases* along that specific atomic displacement pattern at that
specific q-point — the structure is not sitting in a true local energy
minimum along that direction. Causes range from a genuinely unstable phase
(a real soft mode, e.g. approaching a phase transition), an
insufficiently-relaxed input geometry, to pure numerical noise.

That last case is exactly why **a hard `< 0` test is the wrong tool**:
IEEE-754 floating point represents `-0.0` (negative zero) as distinct from,
but numerically equal to, `+0.0` — so `-0.0 < 0` evaluates to **False** in
Python, meaning a structure whose true acoustic-branch frequency is exactly
zero can print as `"-0.0000 THz"` while a hard `< 0` check calls it
*stable* — technically not wrong (`-0.0` really is zero), but with zero
visible reasoning, which reads as an accident rather than a deliberate
rule. Worse in the other direction: a real but tiny residual (e.g.
`-0.0003 THz`, ordinary acoustic-sum-rule numerical noise, see Section 4.2)
would trip a hard `< 0` test and print a scary `[WARNING]` for a structure
that is, in every physical sense, stable. `stb-phononsPos` instead:

1. Calls `phonon.symmetrize_force_constants()` right after loading —
   enforces the acoustic sum rule (translational invariance: the sum of
   forces on all atoms from a rigid whole-crystal translation must be
   exactly zero) — removing most of this numerical noise at the source.
2. Compares against an explicit, documented tolerance,
   `IMAGINARY_MODE_TOL_THZ = -0.01 THz` (`core/phonon_workflow.py`,
   shared with `stb-mlphonons`), instead of a bare `0` — anything less
   negative than this is treated as "zero, within numerical tolerance."

The report's own wording makes the distinction explicit rather than silent
— compare `"No imaginary modes found on the sampled mesh."` (clean,
positive/exact-zero result) against
`"No imaginary modes found along the band path (within numerical
tolerance)."` (a small negative residual, correctly classified stable, but
now honestly labeled as such) — Section 4.2 shows this live, on this
tutorial's own real data.

### 1.5 Thermal properties from the phonon density of states

Once `ω(q)` is known over a dense sampling mesh (`-m/--mesh`, default
`20 20 20`), each mode is an independent quantum harmonic oscillator, and
standard statistical mechanics gives the vibrational contribution to the
free energy, entropy, and constant-volume heat capacity as sums over every
mode at every temperature — this is what `[3] THERMAL PROPERTIES SUMMARY`
reports. Two textbook cross-checks any real result should pass, both
verified live on this tutorial's own data in Section 5:

- **Third law**: entropy `S -> 0` as `T -> 0 K` (a perfect, defect-free
  harmonic crystal has one unique, fully-ordered ground state).
- **Dulong-Petit limit**: heat capacity `Cv -> 3*N_atoms*R` (the classical
  equipartition limit, no longer quantum-mechanically suppressed) as
  `T -> infinity` — for AlP's 2-atom cell, `3*2*8.314 = 49.88 J/K/mol`.

`[3b] DENSITY OF STATES` additionally fits a **Debye model** (a free
-particle approximation, `g(omega) ~ omega^2`) to the low-frequency part of
the real DOS to extract a Debye frequency/temperature — a single-number
summary of the "typical" phonon energy scale, and the origin of the DOS
plot's own blue curve (Section 4.4).

### 1.6 What this workflow simplifies — explicitly

- **Harmonic only.** No phonon-phonon (anharmonic) coupling — no intrinsic
  phonon lifetime/linewidth, no direct thermal-conductivity calculation.
  Thermal *expansion* specifically needs a beyond-harmonic treatment
  (quasi-harmonic approximation, QHA: computing the harmonic phonon
  spectrum at several *volumes* and combining with an equation of state).
  **`stb-mlphonons` has a `--qha` mode; this DFT-based Stage 1/2 pair does
  not** — a real, current gap in the suite, not an oversight in this
  tutorial (Section 6).
- **No non-analytic term correction (NAC).** Polar materials (AlP and GaAs
  both qualify — heteropolar III-V semiconductors) have long-range
  dipole-dipole (Born-effective-charge-mediated) contributions to the
  dynamical matrix that cause LO-TO splitting (a discontinuity in the
  optical branches exactly at Γ) — this workflow has no `BORN` file
  support, so that splitting is simply absent from every band structure
  and DOS this tutorial shows (Section 6).
- **Symmetry-based displacement reduction assumes the input structure
  obeys its own detected symmetry almost exactly.** Real DFT-relaxed
  structures always carry *some* residual numerical noise; Section 5's
  GaAs case is a real, worked demonstration of this assumption breaking
  down in practice, not just a theoretical caveat.
- **Supercell size and `--symprec` are user-chosen, never auto-optimized
  or auto-converged.** Unlike `stb-mlphonons`'s own `--check-convergence`
  (affordable there because MACE evaluations are essentially free) —
  real DFT makes an automatic multi-supercell-size convergence check far
  too expensive to offer here (Section 6).

## 2. Libraries and external dependencies used

- **[Phonopy](https://phonopy.github.io/phonopy/)** (`4.3.1` in this
  environment) is the actual engine behind almost everything in this
  workflow: displacement generation with symmetry reduction, force-constant
  construction from `FORCE_SETS`, the dynamical matrix and its
  diagonalization, band structure/DOS/thermal-property calculations, and
  acoustic-sum-rule symmetrization. `core/phonon_workflow.py` wraps the
  pieces shared between Stage 1/Stage 2 (and 4 other tools:
  `stb-ramanModes`, `stb-irModes`, `stb-herRefs`, `stb-oerRefs`); the rest
  is direct Phonopy API calls in `phonons_create.py`/`phonons_pos.py`.
- **[ASE](https://wiki.fysik.dtu.dk/ase/)** (`3.29.0`) supplies the
  Bravais-lattice classification and q-path generation
  (`Cell.get_bravais_lattice`/`Cell.bandpath`, Section 1.3) — a
  deliberate choice over Phonopy's own `seekpath`-based path finder, since
  `seekpath` is not a dependency of this suite.
- **[spglib](https://spglib.readthedocs.io/)** (`2.7.0`) does the actual
  space-group/symmetry-operation detection, called internally by both
  Phonopy (Stage 1's displacement reduction, Stage 2's `primitive_matrix=
  'auto'` primitive-cell detection) and, elsewhere in the suite, by
  pymatgen's `SpacegroupAnalyzer` (the same underlying library, different
  Python wrapper).
- **numpy** / **scipy** — general array math; `scipy.constants` supplies
  the Planck/Boltzmann constants used for the Debye-temperature conversion
  (`h * freq_Debye / kB`).
- **matplotlib** (`3.11.0`) — the always-saved PNGs (`phonon_bands.png`,
  `thermal_properties.png`) and `--view`'s live preview, all via Phonopy's
  own `plot_band_structure()`/`plot_total_dos()`/`plot_projected_dos()`
  convenience methods, not hand-built figures.
- **gnuplot** (`6.0`, external, not a Python dependency) — only touched
  when `--save-gnuplot` is given; `stb-phononsPos` writes `.dat` + `.gplot`
  file pairs, gnuplot itself is never invoked by the suite.
- **MACE** (via `stb-mlphonons`) — used **only** by this tutorial's own
  `example_4.4.sh`, as a shortcut to get real force constants for Stage 2's
  demo cases without a SIESTA wait. The real DFT-based Stage 1/Stage 2 pair
  this tutorial documents never touches MACE or any ML dependency.

## 3. Stage 1: generating displacements (`stb-phononsCreate`, code `4.4.1`)

### 3.1 What it does

Reads a structure + `calc.fdf`, then writes into `phonon_runs/` (fixed
name, no `-o/--output-dir` flag — unlike `strain_runs`/`elastic_runs`/
`cohesive_runs`):

- `disp-001/`, `disp-002/`, ... — one single-point SIESTA folder per
  symmetry-inequivalent displacement, each with its own `calc.fdf`
  (forced to single-point SCF, Section 3.2) and a copy of the structure
  file with that one displacement applied.
- `phonopy_disp.yaml` — Phonopy's own metadata (supercell, primitive
  matrix, the exact displacement vectors used) — Stage 2 needs this file
  to reconstruct the full force-constant tensor from the finished runs.

### 3.2 Single-point SCF enforcement, live

Every `calc.fdf` — even one written for a full ionic+cell relaxation, like
this tutorial's own (`MD.TypeOfRun CG`, `MD.VariableCell true`,
`MD.Steps 150`) — is forced to a pure single-point SCF in every `disp-*/`
folder, since a finite-difference phonon calculation needs the force at
*exactly* the displaced geometry; any further relaxation during the run
would contaminate the very forces the whole method depends on. This uses
the **same `%include config_extra.fdf` prepend mechanism** as
`stb-elasticInputs` (`core/structure_io.py::prepend_include`) — SIESTA's
fdf reader is first-occurrence-wins for duplicate labels, so prepending the
override at the very top of the generated `calc.fdf` guarantees it wins
regardless of what the template itself says later:

```
$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --no-intro

[1] SINGLE-POINT SCF ENFORCEMENT
------------------------------------------------------------
Calc template (before forcing): MD.TypeOfRun=CG  Steps: MD.Steps=150  MD.VariableCell=true
Forced to single-point SCF via '%include config_extra.fdf'.

Supercell k-grid  : 5 5 5 (auto-suggested, density=0.2 1/Ang)

config_extra.fdf (written into every generated folder):
  # Auto-generated by stb-phononsCreate.
  # Supercell k-grid (auto-suggested, density=0.2 1/Ang) -- takes precedence over --calc's own k-grid,
  # which was tuned for the smaller unit cell.
  kgrid.MonkhorstPack   [5  5  5]

  # Forces a pure single-point SCF (no ionic or cell relaxation) at each displaced
  # supercell's exact geometry, regardless of --calc's own settings.
  MD.TypeOfRun       CG
  MD.Steps           0
  MD.NumCGsteps      0
  MD.VariableCell    false
```

### 3.3 Symmetry reduction and `--symprec`, live

```
[3] SYMMETRY REDUCTION
------------------------------------------------------------
Generating supercell [2, 2, 2] with 0.01 Ang displacements ...

>>> Detected space group: F-43m (216)  (symprec=0.01 Ang) <<<
Quantity                     | Value
-------------------------------------------------------------------
Symmetry precision (symprec) | 0.01 Ang
Space group                  | F-43m (216)
Point group                  | -43m
Symmetry operations          | 192
Displacements needed         | 2 (of 96 without symmetry reduction)
Reduction                    | 97.9% fewer SIESTA runs
```

`--symprec` (default `0.01` Ang) matters far more than it might look —
this session's own suite-wide fix (several tools used to default to
`1e-3`, tighter than pymatgen's real own default of `0.01` and far tighter
than typical DFT-relaxation noise) directly targeted this exact structure.
Reproduced live, with `--symprec 1e-5` (Phonopy's own raw, unexposed
internal default) instead:

```
$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --symprec 1e-5 --no-intro

[3] SYMMETRY REDUCTION
------------------------------------------------------------
>>> Detected space group: R3m (160)  (symprec=1e-05 Ang) <<<
Symmetry precision (symprec) | 1e-05 Ang
Space group                  | R3m (160)
Point group                  | 3m
Symmetry operations          | 48
Displacements needed         | 4 (of 96 without symmetry reduction)
Reduction                    | 95.8% fewer SIESTA runs
```

The *true* zincblende symmetry (`F-43m`, cubic) gets misdetected as the
lower-symmetry `R3m` (trigonal) purely from the input structure's own
tiny relaxation noise (its lattice vectors carry a real residual of
`~8.5e-6 Ang` in the position that should be exactly `0` — smaller than
`1e-5` but this shows the same failure mode at an even tighter tolerance
for contrast) — not just a report-line cosmetic difference: it doubles the
real SIESTA calculations needed (2 -> 4) AND changes which physical
symmetry operations the whole downstream force-constant reconstruction at
Stage 2 relies on.

### 3.4 The supercell k-grid: auto-suggested or explicit

`--calc`'s own k-grid was tuned for the small reference unit cell — using
it unmodified on the much larger phonon supercell would badly *over*-sample
reciprocal space (wasted computation) for no benefit. `stb-phononsCreate`
either auto-suggests a grid from a target k-point **density**
(`--kgrid-density`, default `0.2` 1/Ang — same convention as
`stb-kgrid`/`stb-strain`), or accepts an explicit `--kgrid X Y Z`
override — written into `config_extra.fdf` *first*, ahead of the
single-point-SCF block, so it's the very first directive SIESTA reads:

```
$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --kgrid-density 0.15 --no-intro
Supercell k-grid  : 7 7 7 (auto-suggested, density=0.15 1/Ang)

$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --kgrid 7 7 7 --no-intro
Supercell k-grid  : 7 7 7 (explicit --kgrid)
```

The interactive `stb-suite` menu (`4.4.1`) previews the auto-suggested
grid directly as the bracketed default at its own k-grid prompt
(`Supercell k-grid [5 5 5]: `), so blank == accept, matching this
session's own UX convention rather than a silent "blank means auto"
prompt.

### 3.5 Running it both ways

Both CLI and `stb-suite` -> `4.4.1` produce identical `phonon_runs/`
output for identical answers (Case 9 of `example_4.4.sh` demonstrates the
interactive path directly).

## 4. Stage 2: analyzing the results (`stb-phononsPos`, code `4.4.2`)

### 4.1 Report structure

```
[0] RUN METADATA
[1] FORCE EXTRACTION
[2] DYNAMICAL STABILITY
[2b] BAND STRUCTURE        <- only with --bands
[2c] MODE FREEZE           <- only with --freeze-unstable-mode
[3] THERMAL PROPERTIES SUMMARY
[3b] DENSITY OF STATES     <- only with --dos/--pdos
[3c] THERMAL DISPLACEMENTS <- only with --thermal-displacements
[4] SUMMARY & FILES
```

`-dir/--directory` defaults to `phonon_runs`, matching Stage 1's own fixed
output name — no `cd` or extra flag needed in the common case. `[1] FORCE
EXTRACTION` transparently handles two different input kinds: real SIESTA
`.FA` force files (built into `FORCE_SETS` via `phonopy-init --siesta -f`),
or a `phonopy_disp.yaml` with force constants **already embedded** (an
ML-sourced file, e.g. from `stb-mlphonons` — this tutorial's own Stage 2
shortcut, Section 2) — the latter skips SIESTA-specific extraction
entirely.

### 4.2 Dynamical stability, and reading "(within numerical tolerance)"

Live, on this tutorial's own real (MACE-computed, Section 2) AlP force
constants:

```
$ stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --no-intro

[2] DYNAMICAL STABILITY
------------------------------------------------------------
Supercell used         : 2 x 2 x 2 (4 displacements, 0.0100 Ang each)
Minimum mesh frequency : 0.2572 THz
No imaginary modes found on the sampled mesh.
[INFO] Building auto-detected high-symmetry q-path for band structure ...

[2b] BAND STRUCTURE
------------------------------------------------------------
Bravais lattice   : face-centred cubic (FCC)
Path              : Γ-X-W-K-Γ-L-U-W-L-K | U-X
Minimum band-path frequency : -0.0000 THz
No imaginary modes found along the band path (within numerical tolerance).
```

This is a **live, real reproduction of the exact "-0.00000 looks unstable
but reads stable" scenario** that motivated Section 1.4's fix: the
band-path minimum genuinely displays as `-0.0000 THz` (a tiny negative
residual, not the sampled mesh's own clean positive `0.2572 THz` result —
different q-points, different numerical noise level), and the report now
says so explicitly instead of leaving the reader to wonder whether that's
an oversight. Both are, physically, the same thing: an acoustic branch
correctly going to (numerically) zero at Γ.

### 4.3 Band structure — the zero-line, explained

`phonon_bands.png` (always saved when `--bands` is given) has a thin
dotted **blue** horizontal line at `Frequency = 0` running across every
panel. Traced directly to Phonopy's own source
(`phonopy.phonon.band_structure.BandPlot.decorate()`):

```python
ax.plot([spts[0], spts[-1]], [0, 0], linestyle=":", linewidth=0.5, color="b")
```

This is Phonopy's own **built-in reference guide**, drawn on every band
-structure plot it produces regardless of what tool calls it — stb-suite
adds nothing here. It exists purely so a reader can visually confirm the
acoustic branches actually touch zero at Γ; it is not a plotted mode, not
a warning, not a sign of instability. (This is exactly the finding from
this session's real investigation into the user's own production AlP
results before this tutorial was written.)

### 4.4 DOS/PDOS — the Debye-fit blue curve, and species projection

```
$ stb-phononsPos -dir phonon_runs -m 20 20 20 --dos --view --no-intro

[3b] DENSITY OF STATES
------------------------------------------------------------
Total DOS peak    : 10.8712 THz
Debye frequency   : 8.2113 THz (Debye temperature 394.08 K)
```

With `--dos --view`, Phonopy's own `phonon.plot_total_dos()` overlays a
**second, blue** curve on top of the real (red) DOS — traced to
`phonopy.phonon.dos.plot_total_dos`'s own `color_Debye: str = "blue"`
default argument. This is the Debye quadratic fit
(`g(omega) = Debye_fit_coef * omega^2`) Section 1.5 introduces, drawn from
`omega=0` up to the fitted Debye frequency, then dropping to zero — the
**other** "strange blue line" investigated this session, and unlike the
band-structure one, this one carries real information: it shows how well
the actual DOS matches ideal Debye (free-particle, low-frequency) behavior.
The correct reading for a real 2-atom-basis crystal: the fit and the real
DOS agree well only in the lowest few THz (the purely acoustic region) and
then visibly diverge once optical branches contribute — expected, textbook
behavior, not an artifact.

`--pdos` splits the same DOS per chemical species instead of (or alongside)
the total:

```
$ stb-phononsPos -dir phonon_runs -m 20 20 20 --pdos --thermal-displacements --no-intro

[3b] DENSITY OF STATES
------------------------------------------------------------
  Al PDOS peak    : 10.8712 THz
  P PDOS peak    : 10.8712 THz
```

### 4.5 Thermal displacements (ADP / `U_iso`)

`--thermal-displacements` computes each atom's isotropic mean-square
displacement (the Debye-Waller factor, `U_iso`), the same quantity X-ray
crystallography reports as atomic displacement parameters, and writes a
`tdispmat.cif` for direct visualization:

```
[3c] THERMAL DISPLACEMENTS
------------------------------------------------------------
ADP CIF (@ 1000.0 K) : tdispmat.cif
  Al U_iso @ 1000.0 K : 0.030916 Ang^2
  P U_iso @ 1000.0 K : 0.028382 Ang^2
```

(Computed on a denser, symmetry-**unreduced** mesh with eigenvectors,
reused for `--pdos` too when both are requested together — one extra pass,
not two.)

### 4.6 `--symprec` at Stage 2 — two effects, and a real interaction to respect

`--symprec` here controls **two independent things**: the tolerance
Phonopy uses to re-symmetrize the loaded force constants (Section 1.4),
and the tolerance ASE uses to classify the Bravais lattice for the q-path
(Section 1.3) — both threaded from the same one flag. On this tutorial's
own (already very close to ideal, MACE-relaxed) structure, loosening it
doesn't visibly change anything:

```
$ stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --symprec 0.5 --no-intro
Symmetry precision: 0.5 Ang (force-constant symmetrization + band-path detection)
Bravais lattice   : face-centred cubic (FCC)
```

**A real, found-and-verified interaction**: a Stage-2 `--symprec`
substantially **tighter** than what Stage 1 used to reduce the original
displacement set can fail to reconstruct the force constants at all —
verified live this session against this suite's own `test/4-workflow/
4-phonons/analysis/test.sh` fixture (a real SIESTA `.FA`-based dataset,
Stage 1 run at the default `symprec=0.01`, Stage 2 re-run with
`--symprec 1e-5`):

```
[1] FORCE EXTRACTION
------------------------------------------------------------
Force files found : 42 (label 'Sn3O4')
FORCE_SETS       : generated successfully
[INFO] Initializing Phonopy API and loading FORCE_SETS ...
[ERROR] Could not load phonopy data: Input forces are not enough to
calculate force constants,or something wrong (e.g. crystal structure does
not match).
```

The mechanism: force constants for a symmetry-reduced displacement set are
reconstructed by applying the *detected* symmetry operations to the
computed subset — if Stage 2 detects a **different, lower** symmetry (a
tighter tolerance rejects some of the operations a looser one accepted),
it can no longer explain all the needed tensor elements from what's
actually on disk. **Rule of thumb: loosen, never tighten blindly, if in
doubt** — the CLI help text and the interactive prompt both say this
explicitly now.

### 4.7 `--save-gnuplot` / `--view` — both opt-in this session

```
$ stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --dos --save-gnuplot --no-intro
Files               : phonon_runs/thermal_properties.png, phonon_runs/thermal_properties.dat, phonon_runs/phonon_bands.png, phonon_runs/band.yaml, phonon_runs/phonon_plots/ (7 .dat/.gplot pairs)
```

Gnuplot `.dat`/`.gplot` output used to be written unconditionally on every
run; this session made it opt-in via `--save-gnuplot` (default off),
matching the exact convention `aimd_analysis.py`/`effmass.py`/`coop.py`
already established — matplotlib PNGs (`phonon_bands.png`,
`thermal_properties.png`) are **unaffected**, still always saved either
way. `--view` (also new this session) additionally pops up the computed
figures live via matplotlib (bands, thermal properties, DOS, PDOS) — a
safe no-op under a non-interactive backend like `MPLBACKEND=Agg` (this
tutorial's own script uses it for exactly that reason). Note:
`--view` does **not** cover `--thermal-displacements` — Phonopy's
`plot_thermal_displacements()` needs a *different* internal calculation
(`run_thermal_displacements`, isotropic-only) than the ADP-tensor one this
tool actually runs (`run_thermal_displacement_matrices`), and there is no
matrices-plotting equivalent in Phonopy's own API (Section 6).

The interactive `stb-suite` menu (`4.4.2`) also gained an **`"all"`**
shortcut for its additional-analyses list this session — expands to
`bands+dos+pdos+thermal`, deliberately **excluding** mode-freeze (a
different kind of action, writing a structure file, not just a
plot/report):

```
Select by number (1-5), space/comma-separated, 'all' for 1-4, or blank for none: all
Selected: bands, dos, pdos, thermal
```

### 4.8 Running it both ways

Both CLI and `stb-suite` -> `4.4.2` produce the identical report for
identical answers (Case 9 of `example_4.4.sh`).

## 5. Worked example: two real materials, end-to-end (AlP clean, GaAs broken)

This is not illustrative data — both investigations happened this same
session, on real production folders the user actually ran (`stb
-phononsCreate` + `stb-phononsPos`, real SIESTA, DZP/GGA-PBE, `dojo`
pseudopotentials).

### 5.1 AlP: clean, stable, and a direct confirmation of this session's own fix

Stage 1's real, saved report on this exact structure (`structure.fdf`/
`calc.fdf` in this folder) already shows the correct result at the fixed
default: `F-43m (216)`, 192 symmetry operations, 2 of 96 displacements
needed (Section 3.3). Stage 2, run on the real SIESTA-computed force
constants (not this tutorial's own MACE shortcut):

- **Every one of the 6060 sampled band-path frequencies checked directly**
  (not just visually from the plot): minimum `+9.95e-08 THz` — positive,
  essentially exact zero at Γ, **zero negative frequencies anywhere**.
  Genuinely, cleanly dynamically stable.
- **Heat capacity at 1000 K: 48.97 J/K/mol**, vs. the Dulong-Petit
  classical limit of `49.88 J/K/mol` (Section 1.5) — 98.2% of the limit,
  exactly the plateau the theory predicts.
- **Entropy at 0 K: exactly 0** — Third Law respected.
- **DOS peak at 12.89 THz**, matching the flat optical-branch pileup
  visible in the band structure near X and along L-U-W-L.
- **Lattice constant**: this structure's own relaxed cell implies a
  conventional cubic lattice constant of `~5.53 Ang`, vs. AlP's real
  experimental value of `~5.46 Ang` — a `~1.2%` overestimate, exactly the
  expected direction/magnitude of PBE-GGA's own well-known systematic
  lattice-constant bias, not a red flag.

Every one of these is an independent physical cross-check, and every one
passes — strong, real evidence that both the physics *and* this session's
`--symprec` fix (Section 3.3) are working correctly end-to-end, not just
in a synthetic test fixture.

### 5.2 GaAs: a genuine, forensically-traced instability

The real GaAs production folder (same setup, `--symprec 0.05` explicitly
at Stage 1) tells a very different story. Stage 2's own saved report:

```
[2] DYNAMICAL STABILITY
------------------------------------------------------------
Minimum mesh frequency : -5.8340 THz
[WARNING] Negative (imaginary) phonon frequencies found ...

[2b] BAND STRUCTURE
------------------------------------------------------------
Minimum band-path frequency : -5.8363 THz
[WARNING] Negative (imaginary) frequency found along the band path ...
```

This is **not** the kind of near-zero noise Section 1.4/4.2 tolerates —
it's large, real, and visible directly on the band plot: the 3 branches
that should be acoustic (rising from 0 at Γ) instead dive to deeply
negative values away from Γ, almost mirroring the optical branches. Both
SIESTA displacement runs completed normally (`0_NORMAL_EXIT`, no failed
calculation feeding bad data in).

**Hypothesis 1, tested and ruled out live**: a Stage-1/Stage-2
`--symprec` mismatch (Section 4.6's real, documented risk). Re-running
Stage 2 with `--symprec 0.05` (matching Stage 1 exactly, instead of the
tighter `0.01` default) barely changes anything (`-5.816 THz` vs. the
original `-5.834 THz`) — this specific instability is **not** caused by a
tolerance mismatch.

**Hypothesis 2, supported by direct evidence**: checking the dynamical
matrix eigenvalues at Γ directly (not through the CLI, via Phonopy's
Python API) shows the acoustic branches are **essentially perfect there**
(`~1e-8 THz`, and the acoustic-sum-rule residual before any correction is
already tiny, `~1.5e-8`) — the force constants are **not** globally broken;
the problem is specific to q-points away from Γ, which depend on how
accurately the symmetry-based expansion (Section 1.2/4.6) reconstructed
the *full* force-constant tensor from just 2 real displacements. Comparing
the two structures' own residual numerical noise (the tiny near-zero term
in each `%block LatticeVectors`):

| Structure | Residual noise | 
|---|---|
| AlP (clean, stable) | `0.00000852 Ang` |
| GaAs (unstable) | `0.00395996 Ang` |

**GaAs carries ~465x more symmetry-breaking noise than AlP** in the same
structural position — small enough that `spglib` still correctly
*identifies* the space group as `F-43m`, but evidently not small enough
for the symmetry-operation-based force-constant *expansion* (which applies
those operations to the real, slightly-noisy atomic positions) to remain
accurate for atom pairs far apart in the supercell — exactly consistent
with "correct at Γ, broken away from it."

**Conclusion and recommendation**: this is a genuine data-quality/
convergence issue with the *input structure and/or supercell size* for
this specific material, not a stb-suite bug — the code is working
correctly and already flagged the problem with its own `[WARNING]`. Two
concrete next steps, neither of which could be verified further without
running new DFT (a real limitation of this tutorial, not glossed over):
(1) symmetry-refine the structure to the exact detected space group before
Stage 1 (e.g. `stb-unitcell --mode refined`) to remove the `0.004 Ang`
-scale noise at the source, and/or (2) increase the supercell size (e.g.
`3x3x3`) to check whether GaAs's real interatomic force range is simply
longer than what a 16-atom supercell can capture — the same
size-convergence assumption AlP happened to satisfy but GaAs, evidently,
does not.

## 6. Known, deliberate limitations

- **No non-analytic term correction (NAC/Born effective charges)** —
  LO-TO splitting at Γ is absent from every band structure/DOS this
  workflow produces, for any polar material (Section 1.6).
- **Harmonic only** — no anharmonic phonon lifetimes/linewidths, no
  thermal expansion for this DFT-based pair (`stb-mlphonons --qha` exists
  for the ML path only, Section 1.6).
- **Supercell size and `--symprec` are user-chosen, never auto-converged**
  — Section 5.2's GaAs case is the concrete demonstration of why this
  matters in practice, not just a caveat.
- **A Stage-2 `--symprec` tighter than Stage 1's can break force-constant
  reconstruction outright** (Section 4.6) — a real, sharp-edged risk, not
  a gentle accuracy tradeoff.
- **`--view` has no `--thermal-displacements` coverage** — Phonopy's
  `plot_thermal_displacements()` needs a different internal calculation
  than the ADP-tensor one this tool actually runs; there's no
  matrices-plotting equivalent in Phonopy's own API today (Section 4.7).
- **Symmetry-based displacement reduction assumes near-exact input
  symmetry** — real residual relaxation noise degrades the reconstructed
  force constants for distant neighbors specifically (Section 5.2), with
  no automatic detection of when this has actually happened (the mesh/
  band-path stability check is the only downstream signal).
- **`--kgrid-density`'s 0.2 1/Ang default, `--distance`'s `0.01` Ang
  default, and the supercell dimensions themselves are never
  auto-converged** — same "user-chosen, not automated" limitation as
  `4.1`/`4.2`/`4.3`'s own strain range/mesh/cutoff choices.

## 7. Step-by-step: running this workflow on your own structure

0. **Prerequisite**: an already-relaxed structure (`structure.fdf`) and the
   **same `calc.fdf`** used for that relaxation (Section 3.2 — Stage 1
   reads it only to know the basis/functional/mesh, then forces
   single-point SCF itself).

1. **Decide the supercell size** (`-dim`, e.g. `2 2 2`). Bigger captures
   longer-ranged interatomic forces more accurately but costs more SIESTA
   runs even after symmetry reduction — Section 5.2's GaAs case is a real
   example of a size that turned out to be too small for that specific
   material.

2. **Decide `--symprec`** (default `0.01`, Section 3.3) — start at the
   default; if the input structure is unusually noisy (freshly relaxed
   with loose force tolerances, say), consider checking with
   `stb-symmetry --scan-symprec` first (same tool `4.3`'s own tutorial
   points to).

3. **Run Stage 1** (`stb-phononsCreate`, code `4.4.1`):
   - CLI: `stb-phononsCreate -s structure.fdf -c calc.fdf -p <bank> -dim
     <nx> <ny> <nz>`.
   - Interactive: `stb-suite` -> `4.4.1`.
   - Writes `phonon_runs/disp-*/` (Section 3.1).

4. **Run SIESTA yourself in every `disp-*/` folder.** Stage 1 never runs
   SIESTA.

5. **Run Stage 2** (`stb-phononsPos`, code `4.4.2`) once every displacement
   finishes:
   - CLI: `stb-phononsPos -dir phonon_runs -m <mx> <my> <mz> --bands --dos
     --pdos --thermal-displacements` (mix and match the extra-analysis
     flags you want).
   - Interactive: `stb-suite` -> `4.4.2` (the `"all"` shortcut picks
     bands+dos+pdos+thermal in one go, Section 4.7).

6. **Read `[2]`/`[2b]` first** (Section 4.2) — if you see a `[WARNING]`
   that isn't explained by "(within numerical tolerance)", the geometry
   may genuinely not be a true energy minimum, or (Section 5.2) the
   supercell/symmetry assumptions may not hold for this material — don't
   trust `[3]`'s thermal properties until this is resolved (Phonopy itself
   silently drops imaginary modes from those sums, so they'd be wrong, not
   just absent).

7. **(Optional) persist or view**: `--save-report` writes
   `phonon_properties.txt`; `--save-gnuplot` writes `.dat`/`.gplot` pairs
   (off by default, Section 4.7); `--view` pops up the computed plots
   live.

## 8. Files in this folder

| File | What it is |
|---|---|
| `structure.fdf` | Real, relaxed AlP, zincblende, 2 atoms/cell — the material Section 5.1's real numbers are for |
| `calc.fdf` | The matching real calc template (DZP, PBE-GGA, 320 Ry, 11x11x11 k-grid) |
| `example_4.4.sh` | This walkthrough's runnable script (both stages, 9 cases) |

GaAs's own structure (Section 5.2) is not checked into this folder — that
investigation used the real production data directly, not a tutorial
fixture; the numbers quoted above are real and independently reproducible
against that structure, just not shipped here.

## 9. Running the script

```bash
bash example_4.4.sh
```

## 10. What's next

- **`4.1-strain`** / **`4.2-elastic`** / **`4.3-cohesive`** — the other
  mechanical/energetic questions about the same kind of already-relaxed
  structure this workflow starts from (deformation response, stiffness
  tensor, binding energy) — phonons answer a genuinely different one
  (vibrational stability and thermal behavior).
- **`stb-mlphonons`** (ML Simulations, code `5.2`) is this workflow's own
  ML-preview twin — same physics, MACE-MP-0 instead of SIESTA, no SIESTA
  wait, plus features this DFT-based pair doesn't have (`--qha` thermal
  expansion, `--check-convergence` supercell-size scanning, foundation
  -model comparison) — exactly the shortcut this tutorial's own script
  used to get real Stage 2 data without waiting on SIESTA.
