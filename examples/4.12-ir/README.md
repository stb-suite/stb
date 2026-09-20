# 4.12 — Workflow: IR Spectrum (`stb-ir` / `stb-irModes` / `stb-irAnalysis`)

This workflow has 3 stages: **Stage 1** (`stb-ir`, code `4.12.1`) takes an
already-relaxed structure and its own `calc.fdf`, builds a supercell,
symmetry-reduces which atoms actually need to be displaced, and writes one
ready-to-run, single-point SIESTA folder per displacement — the same
finite-displacement machinery `4.4-phonons` and `4.11-raman` share.
**Stage 2** (`stb-irModes`, code `4.12.2`) reads those forces back, finds
the vibrational modes at the Γ point, classifies which ones are IR-active
by symmetry, and writes whichever *new* SIESTA folder(s) the intensity
calculation needs — the choice of folder depends entirely on how periodic
the structure is (Section 1.3). **Stage 3** (`stb-irAnalysis`, code
`4.12.3`) reads those back and produces per-mode intensities, a synthetic
spectrum, and (optionally) an extended-XYZ animation of each mode.

All three stages live in this one folder and this one tutorial — Stage
1's output only exists to feed Stage 2, and Stage 2's only to feed Stage 3.

## 1. Theory

### 1.1 The harmonic approximation, briefly (full derivation: `4.4-phonons`)

Every phonon workflow in this suite (`4.4-phonons`, `4.11-raman`, this one)
starts the same way: displace each atom by a tiny amount `d` (Å), measure
the resulting force on every other atom, and build the force-constant
matrix `Φ_ij = -∂F_i/∂u_j` from the finite differences. Diagonalizing the
mass-weighted force-constant matrix at the Γ point (`q=0`) gives `3N`
frequencies and eigenvectors — 3 of them acoustic (translations, ω=0), the
rest optical. Symmetry cuts the number of *independent* displacements
needed dramatically (Section 3.3) — see `4.4-phonons` README Section 1.2
for the group-theory argument in full; this tutorial doesn't repeat it.

### 1.2 What makes a mode IR-active: `dμ/dQ`

A phonon mode is Raman-active if the crystal's *polarizability* changes
along it (`4.11-raman`'s own theory, second-order in the atomic
displacement). It's **IR**-active if the crystal's **dipole moment**
changes along it instead — first order, and governed by each atom's own
**Born effective charge tensor**, `Z*_κ,αβ` (`α`=field direction,
`β`=displacement direction):

```
dμ_α/dQ  =  Σ_κ Σ_β  Z*_κ,αβ · e_κ,β
```

where `e_κ,β` is mode `Q`'s eigendisplacement on atom `κ`, normalized to
unit norm across the whole cell. Intensity is `|dμ/dQ|²` (times a
frequency-dependent prefactor for the spectrum's Lorentzian weighting,
Section 5.3). Symmetry alone can already tell you which modes are
*allowed* to be IR-active (Section 1.4) — a centrosymmetric crystal's
`u`-parity modes (Raman-active) and `g`-parity modes (IR-active) are
mutually exclusive (Mutual Exclusion Rule), the reason NaCl's own single
optical branch (Section 6) is IR-active but Raman-silent.

### 1.3 Three regimes for `Z*`/`dμ/dQ`, and why they're genuinely different physics

`Z*` (equivalently, `dμ/dQ`) is **not** computed the same way for every
structure — how periodic the structure is along a given Cartesian
direction changes what "the dipole moment" even *means* there:

- **0D (molecule) / 1D (wire)** — along any non-periodic direction, SIESTA
  prints the total dipole moment directly (a well-defined, finite
  quantity for a finite or quasi-1D system). Stage 2 writes a `+delta`/
  `-delta` displaced pair per selected mode; Stage 3 central-differences
  the two dipole values. Simple, and exact.
- **3D bulk** — a fully periodic crystal has **no well-defined dipole
  moment** in the naive sense (moving every atom's position by one lattice
  vector is a gauge choice, not a physical difference, so `Σ q_i r_i`
  is ill-defined for a periodic charge distribution). The physically
  correct quantity is the **Berry-phase modern theory of polarization**
  (King-Smith & Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993)) — SIESTA
  automates this via `BornCharge T` + `MD.TypeOfRun FC` +
  `%block PolarizationGrids`, giving `Z*` **directly** from ONE
  equilibrium run, valid for *every* mode at once (no per-mode
  displacement needed at all — Section 4.2).
- **2D slab (one vacuum-padded axis)** — a genuine hybrid: the two
  in-plane (periodic) axes need the same Berry-phase machinery as bulk;
  the vacuum axis is non-periodic, so a naive dipole moment IS
  well-defined there and the 0D/1D method applies. **This is not a
  theoretical nicety** — it was root-caused as a real, silent bug this
  session: monolayer h-BN's dominant in-plane IR mode (E′, literature
  ~1420 cm⁻¹) came back with **exactly zero** computed intensity, because
  SIESTA's own raw dipole printout is an exact zero for any in-plane
  displacement of a periodic slab (confirmed directly in `calc.out`, not a
  parsing bug). `stb-irModes`' HYBRID path (Section 4.2) fixes this by
  computing the in-plane component via Born charges and the out-of-plane
  component via the dipole-difference method, on the SAME structure,
  combined per-mode. Empirically validated (not just theoretically
  argued) before being coded: the corrected h-BN E′ mode intensity came
  back nonzero and physically sane (~15, vs. the A2″ out-of-plane mode's
  ~0.8), and the frequency (after a separate, unrelated pseudopotential
  fix) matched literature almost exactly (1372 cm⁻¹ computed vs. 1372
  cm⁻¹, Geick, Perry & Rupprecht, *Phys. Rev.* **146**, 543 (1966)).

This workflow auto-detects which of the three applies from the number of
vacuum-padded axes (`kspace.detect_vacuum_axes`, `--vacuum-gap` threshold,
default 10 Å) — you never choose the path yourself.

### 1.4 Symmetry, Mulliken labels, and degenerate modes

Stage 2 always runs a full space-group/point-group analysis (via
`spglib`, through Phonopy) on the Γ-point eigenvectors and reports, per
mode: its irreducible representation (Mulliken label — `T1u`, `A2″`,
`E′`, ...) and whether that irrep is IR-active *by symmetry* (not just
"happens to have nonzero computed intensity" — a genuinely silent mode
should show ~0 intensity as a cross-check of the symmetry classification
itself). A **degenerate** irrep (dimension > 1, e.g. NaCl's `T1u`, h-BN's
`E′`) means several modes share exactly one frequency by symmetry — but
Phonopy's own choice of eigenvector *basis* within that degenerate
subspace is arbitrary (any orthogonal rotation within it is equally
valid). This has a direct, non-obvious consequence for IR intensity:
**each individual partner's own `dμ/dQ` *direction* is basis-dependent
and not itself physically meaningful** — only the **group's summed
intensity** (`Σ|dμ/dQ|²` over the whole degenerate subspace) is
basis-independent (an orthogonal change of basis preserves the sum of
squared projections). Stage 3 computes and reports this combined
quantity explicitly (Section 5.2) — it's the number to compare against a
single experimental peak, not any one partner's own row.

### 1.5 What this workflow simplifies — explicitly

- **No non-analytic (LO-TO) correction for polar bulk crystals** — the
  single most consequential limitation of this workflow, root-caused (not
  fixed). Full explanation: Section 6.4. `4.4-phonons`
  documents the identical limitation for its own (intensity-free) phonon
  band structures/DOS — this is a shared, suite-wide gap, not specific to
  IR.
- **Harmonic only** — same as every other phonon tool in this suite; no
  anharmonic linewidths/lifetimes, no temperature dependence of the
  frequency itself.
- **Supercell size, `--symprec`, `-d`/`--distance`, and `--kgrid-density`
  are all user-chosen, never auto-converged** — same "you decide,
  nothing auto-converges" limitation as `4.4-phonons`/`4.11-raman`.
- **1D wires (2 vacuum axes) stay on the plain 0D/1D dipole-difference
  path** — not the HYBRID treatment, deliberately out of scope for now
  (no 1D polar test case existed when the HYBRID path was built; the
  same underlying issue Section 1.3 describes for 2D almost certainly
  also applies to a periodic 1D polar wire, just unverified).

## 2. Libraries and external dependencies used

- **Phonopy** (`phonopy.interface.siesta`) — supercell generation,
  symmetry-reduced displacements, FORCE_SETS construction, Γ-point
  eigenvector/eigenvalue extraction. Same shared machinery as
  `4.4-phonons`/`4.11-raman`, via `core/phonon_workflow.py`.
- **spglib** (through Phonopy) — space-group/point-group detection and
  irrep classification for the symmetry analysis (Section 1.4).
- **SIESTA** itself — every actual force/dipole/Born-charge number comes
  from real DFT; nothing in this workflow is machine-learned or
  approximated (contrast with `5.x` ML Simulations' `stb-mlphonons`,
  which reuses the same Phonopy machinery but MACE-MP-0 forces instead).
- **ASE** (`ase.io.write`, extended XYZ) — Stage 3's mode-vibration
  animations (Section 5.4).
- **matplotlib** — Stage 3's optional `--view` spectrum preview, and this
  README/script's own literature-comparison plot (Section 6.3).

## 3. Stage 1: generating displacements (`stb-ir`, code `4.12.1`)

### 3.1 What it does

Reads `structure.fdf` + `calc.fdf`, builds a `-dim a b c` supercell,
symmetry-reduces the displacement list (Section 3.3), and writes one
`ir_study/phonon_disp/disp-NNN/` folder per required displacement — each a
complete, ready-to-run SIESTA input (structure + your `--calc` template +
an auto-written `config_extra.fdf`, Section 3.2). Also prints a `[1] INPUT
STRUCTURE` table (composition, atom count, cell volume, dimensionality,
supercell size) and an advisory `[NOTE]` telling you, in advance, which of
Section 1.3's three Stage-2 paths this structure will take.

### 3.2 Single-point SCF enforcement, live

Every generated folder is forced to a pure single-point SCF via
`%include config_extra.fdf` prepended at the top of `calc.fdf` (SIESTA's
fdf reader is first-occurrence-wins for a duplicate label, so this always
wins even over a relaxation template's own `MD.TypeOfRun CG`) — a
displaced supercell's force must be measured at *exactly* the displaced
geometry, or the finite-difference force constants (and every frequency/
intensity derived from them) are silently corrupted:

```
$ stb-ir -s structure.fdf -c calc.fdf -p dojo -dim 3 3 3 --kgrid-density 0.15 --symprec 0.01 --no-intro
```
```
[2] SINGLE-POINT SCF ENFORCEMENT
------------------------------------------------------------
Calc template (current state, before forcing):
  MD.TypeOfRun=CG  Steps: MD.Steps=150  MD.VariableCell=true
...
Supercell k-grid  : 5 5 5 (auto-suggested, density=0.15 1/Ang)

config_extra.fdf (written into every generated folder):
  kgrid.MonkhorstPack   [5  5  5]
  MD.TypeOfRun       CG
  MD.Steps           0
  MD.VariableCell    false
```

### 3.3 Symmetry reduction and `--symprec`

Phonopy's own raw `symprec` default (`1e-5`) is far tighter than any real
DFT-relaxed structure's residual numerical noise and can silently
**misdetect a lower symmetry than the crystal actually has** — `stb-ir`
uses `0.01` (pymatgen's own default) instead, matching `4.4-phonons`/
`4.11-raman`'s identical fix. The NaCl structure this example ships is
numerically very clean (built via `stb-crystalbuilder` at the exact
Fm-3m Wyckoff positions, then relaxed), so tightening `--symprec` here
does not, in practice, change anything (Case 2 of the script shows this
honestly — 2 displacements either way). For a structure where this DOES
bite in practice, see `4.4-phonons` README Section 5.2 (a real, forensically
-traced GaAs misdetection).

### 3.4 The supercell k-grid: auto-suggested or explicit

`--kgrid-density` (default `0.2` 1/Å) auto-suggests a Monkhorst-Pack grid
for the supercell's own single-point SCF — this is what the phonon
frequencies are actually built from, so an under-converged k-grid can move
the reported frequency significantly (dramatically so for a metal/
semimetal with a Kohn anomaly, verified live on graphene's own G-band this
session while investigating `stb-raman`: >100 cm⁻¹ shift from doubling the
density alone — see `stb-raman --help`'s own epilog Notes for the same
finding written into the tool itself).
`--kgrid` overrides it with an explicit grid. NaCl is a wide-gap insulator
and far less sensitive than graphene, but the mechanism is identical.

### 3.5 Running it both ways

CLI: `stb-ir -s structure.fdf -c calc.fdf -p dojo -dim 3 3 3 --kgrid-density 0.15 --symprec 0.01`.
Interactive: `stb-suite` → `4.12.1`, same questions in the same order.
Case 1/2/3/7 in the script demonstrate both.

## 4. Stage 2: modes & displacement folders (`stb-irModes`, code `4.12.2`)

### 4.1 Report structure

`[1] PHONON MODES AT GAMMA` (frequencies, symmetry labels, IR activity,
Cartesian-polarization character e.g. "98% z") → `[1b] SYMMETRY ANALYSIS`
(space group, point group, degenerate groups) → `[2] PSEUDOPOTENTIALS` →
`[3] IR DISPLACEMENT FOLDERS` (what gets written, Section 4.2) → `[4]
SUMMARY & NEXT STEPS`.

### 4.2 The three intensity paths, concretely

The `[Path]` line names which of Section 1.3's regimes applies and what
gets written:

| Path | Detected when | Folders written |
|---|---|---|
| `BULK` | 0 vacuum axes | ONE `born_charge_disp/equilibrium/` (covers every mode) |
| `HYBRID` | exactly 1 vacuum axis | the same shared equilibrium folder (in-plane) PLUS `dipole_disp/mode_NN_plus\|minus/` per mode (vacuum axis) |
| `NONBULK` | 2-3 vacuum axes | `dipole_disp/mode_NN_plus\|minus/` per mode only |

NaCl (this example) is `BULK`:

```
[Path] Bulk (3D periodic) -- one equilibrium Born-effective-charge SIESTA run
(MD.TypeOfRun FC + BornCharge T + PolarizationGrids) covers every selected
mode's IR intensity.
Non-acoustic Gamma modes : 3
  mode   1 (band   3) :     5.5653 THz  (T1u) [IR-active]  [50% z, 50% x]
  mode   2 (band   4) :     5.5653 THz  (T1u) [IR-active]  [50% z, 50% y]
  mode   3 (band   5) :     5.5653 THz  (T1u) [IR-active]  [50% x, 50% y]
```

Only **one** SIESTA folder is written no matter how many modes are
selected — the Born-charge tensor `Z*` doesn't depend on which mode you
later combine it with (Section 1.2's formula makes this explicit: `Z*` is
a per-atom property of the *equilibrium* structure, `e_κ` is what varies
per mode).

### 4.3 A real, known SIESTA bug you will hit running this for real

If you have SIESTA installed and run the `BornCharge`/`Optical`-module
equilibrium calculation yourself with `mpirun -np N` for `N>1`, it fails
with `cdiag: Error in Cholesky factorisation` — a real, reproducible bug
in SIESTA's own parallel diagonalization specific to the polarization/
Optical module (confirmed directly: `-np 1` works, `-np 2`
reproduces the failure, `-np 6` also fails; supercell force runs and
dipole-difference runs elsewhere in the SAME workflow are unaffected and
can use any rank count). The fix is simply `-np 1` for that one folder —
`config_extra.fdf`/`calc.fdf` themselves need no change, this is a SIESTA
binary limitation, not something `stb-irModes` writes wrong. The script
(Case 4) documents the exact commands.

### 4.4 Running it both ways

CLI: `stb-irModes -dir ir_study -p dojo -c calc.fdf --symprec 0.01`.
Interactive: `stb-suite` → `4.12.2`, one question per flag (`--modes`/
`--freq-min`/`--freq-max` mode selection, `--use-symmetry`,
`--skip-degenerate`, `--export-animations`). Not scripted turn-by-turn in
`example_4.12.sh` (Section 1.3's 3-way path makes the exact prompt count
depend on what Stage 1 found) — run it live yourself against this
example's own `output/workflow/ir_study/`.

## 5. Stage 3: analysis (`stb-irAnalysis`, code `4.12.3`)

### 5.1 Report structure

`[1] READING RUNS` (per-mode `dμ/dQ`, raw) → `[2] IR-ACTIVE MODES SUMMARY`
(table + degenerate-group combining, Section 5.2) → `[3] SPECTRUM` (peak
list, optional `[3b] EXPERIMENTAL COMPARISON`) → `[3c] MODE VIBRATIONS`
(always-on animation export) → `[4] SUMMARY & FILES`.

### 5.2 Combining degenerate modes

Section 1.4 already explained *why*: within NaCl's one degenerate `T1u`
group, each individual mode's own `dμ/dQ` direction is an arbitrary basis
choice —

```
Mode  THz       cm^-1     dmu/dQ_x    dmu/dQ_y    dmu/dQ_z    Intensity
1     5.5653    185.64    0.743228    0.024350    0.743569    1.105875
2     5.5653    185.64    0.024356    0.743331    0.743331    1.105676
3     5.5653    185.64    0.743435    0.743094    0.024358    1.105477

[Degenerate groups] modes 1, 2, 3 (5.5653 THz / 185.64 cm^-1) -> combined intensity 3.317027
```

Note the three individual intensities (~1.10 each) are NOT what should be
compared to a single experimental peak — their sum, **3.317027**, is. Note
also each row's own `dμ/dQ` picks out a different pair of Cartesian axes
(x,z / y,z / x,y respectively) purely as an artifact of Phonopy's
arbitrary eigenvector basis inside the degenerate subspace — physically,
NaCl's cubic symmetry makes the mode isotropic, and only the rotationally
-invariant combined intensity reflects that.

### 5.3 Spectrum, peaks, and `--experimental`

`[3] SPECTRUM` builds a Lorentzian-broadened synthetic spectrum (default
linewidth `10 cm⁻¹`) from every mode's frequency+intensity and always
lists detected peaks (`scipy.signal.find_peaks`, prominence-gated). NaCl
gives exactly one, at `185.59 cm⁻¹`, matching the single triply-degenerate
mode. `--experimental FILE` (2-column: cm⁻¹, intensity) overlays a real
measured spectrum and attempts peak-matching — this needs an actual
digitized/broadened experimental curve with several points to work (the
peak-finder needs neighboring samples to define prominence), **not** a
bare single-frequency literature value, which is why this example
compares against the literature TO frequency via a plain reference line on
a plot instead (Section 6.3), the same technique used for h-BN's own
literature comparison.

### 5.4 Mode vibrations / animations

`[3c] MODE VIBRATIONS` always (no flag needed) reloads Stage 1's force
constants and writes one extended-XYZ animation per analyzed mode under
`<directory>/mode_animations/` — open any of them in a molecular viewer,
or pass `--view-modes` to preview them directly via ASE's interactive
viewer during the run.

### 5.5 Running it both ways

CLI: `stb-irAnalysis -dir ir_study --save-gnuplot`. Interactive:
`stb-suite` → `4.12.3`, same `--save-gnuplot`/`--view`/`--view-modes`/
`--animation-frames`/`--experimental` questions.

## 6. Worked example: NaCl end-to-end, and a real, currently unfixed limitation

This is not illustrative data — this exact investigation happened this
session, on a real production folder (`stb-crystalbuilder` →
`stb-unitcell` → a real SIESTA relaxation → this workflow, all `dojo`
pseudopotentials, DZP/GGA-PBE).

### 6.1 Building and relaxing NaCl

`stb-crystalbuilder` (space group `Fm-3m`, Wyckoff sites Na `(0,0,0)`,
Cl `(0.5,0.5,0.5)`, starting lattice constant `a=5.6402 Å`, the accepted
experimental value) → `stb-unitcell` (conventional → primitive, 2 atoms)
→ a real SIESTA CG relaxation (`kgrid.MonkhorstPack [10 10 10]`,
`Mesh.CutOff 320 Ry`) converged to a primitive lattice vector modulus of
`3.988224 Å` (conventional `a=5.6395 Å` — a `0.01%` change from the
starting guess; PBE and experiment happen to agree almost exactly for
NaCl's lattice constant, unlike the systematic `~1%` PBE overestimate seen
for AlP in `4.4-phonons` Section 5.1).

### 6.2 Stage 1–3 results

- **Stage 1**: `-dim 3 3 3` (54-atom supercell), `2` of `324` possible
  displacements needed (Fm-3m symmetry), `--kgrid-density 0.15` → auto
  `5×5×5` supercell k-grid.
- **Stage 2**: space group `Fm-3m (225)`, point group `m-3m`, `1296`
  symmetry operations, **1** non-acoustic Γ mode, triply degenerate,
  irrep `T1u`, IR-active, frequency **5.5653 THz = 185.64 cm⁻¹**. Born
  effective charges: Na `+1.0529`, Cl `−1.0544` (isotropic, as cubic
  symmetry requires; acoustic sum rule `+1.0529 − 1.0544 = −0.0015`,
  essentially exact).
- **Stage 3**: combined mode intensity **3.317027**, single spectral peak
  at **185.59 cm⁻¹**.

### 6.3 Comparison to literature: a real, honest ~13% gap

NaCl is *the* textbook example for the Lyddane-Sachs-Teller relation
(`ω_LO²/ω_TO² = ε(0)/ε(∞)`). Inelastic-neutron-scattering literature:

> Raunio, Almqvist & Stedman, "Phonon Dispersion Relations in NaCl",
> *Phys. Rev.* **178**, 1496 (1969) — the standard reference for NaCl's
> zone-center TO/LO phonon frequencies, `ω_TO ≈ 164 cm⁻¹`,
> `ω_LO ≈ 264 cm⁻¹`.

Cross-checked independently via the VASP Wiki's own worked NaCl dielectric
example (`ε(0)=5.9`, `ε(∞)=2.34`, `ω_TO=0.02 eV=161.3 cm⁻¹`), giving
`ω_LO = ω_TO·√(ε(0)/ε(∞)) ≈ 256 cm⁻¹` — consistent with Raunio et al.
within a few percent.

Our computed **185.64 cm⁻¹** sits **~13% above the TO value**, and —
tellingly — almost exactly **midway** between the literature TO and LO
values. This is the physics content of Section 6.4, not a convergence
problem: k-grid density (Section 3.4) and supercell size were both within
normal ranges for this system, and this exact ~13% gap reproduces
identically regardless.

### 6.4 Root cause: no non-analytic (LO-TO) correction for a polar crystal

`stb-ir`'s Γ-point frequency comes **purely** from real-space
finite-displacement force constants (Section 1.1) — the Born charges
computed in Stage 2's `BULK` path are used **only** for IR intensity
(Section 1.2), never fed back into the frequency/dynamical-matrix
calculation. For a **polar** crystal (nonzero `Z*`, like NaCl's ionic
Na⁺/Cl⁻ bonding), this is a real, well-documented gap in the underlying
method, not specific to this suite:

> "An approach that the phonopy code employs is so-called non-analytical
> term correction (NAC), where dipole–dipole interaction is calculated
> from static dielectric constant tensor and Born effective charges...
> To apply this approach to the supercell force constants, **a part of
> the dipole–dipole interaction included at the supercell force constants
> calculation is subtracted from the supercell force constants to avoid
> double counting the contribution**." — Togo & Tanaka, "First-principles
> Phonon Calculations with Phonopy and Phono3py", *J. Phys. Soc. Jpn.*
> **92**, 012001 (2023).

In plain terms: the electrostatic dipole-dipole interaction has `1/r³`
range and does **not** decay within any finite periodic supercell — a
displaced atom's induced dipole spuriously couples to its own periodic
images across the supercell boundary. Without explicitly subtracting this
spurious, supercell-size-dependent contribution and replacing it with the
analytically correct (infinite-crystal, Ewald-summed) term — which needs
both `Z*` (already computed) **and** the high-frequency dielectric tensor
`ε(∞)` (not currently computed by any `stb` tool) — the reported
Γ-frequency for a polar material's IR-active mode is contaminated, and can
land anywhere between the true TO and LO values depending on supercell
size. NaCl's `185.6 cm⁻¹`, landing almost exactly midway between `164`
and `264 cm⁻¹`, is the textbook signature of exactly this (Gonze & Lee,
*Phys. Rev. B* **55**, 10355 (1997), is the original formal treatment;
`4.4-phonons` README Section 6 documents the identical limitation for its
own phonon band structures, for the same underlying reason).

This does **not** affect non-polar bulk crystals (`Z*≈0`: Si, diamond,
graphene) — nor, apparently, h-BN's out-of-plane A2″ mode as strongly
(Section 1.3), plausibly because the relevant dipole-dipole coupling is
more effectively screened by the slab's own vacuum padding; NaCl's fully
3D-periodic ionic lattice is the cleanest, strongest case for seeing this
limitation in this suite today.

**This limitation is real, understood, and currently unfixed in `stb-ir`
— every report from a `BULK`-path run now prints an explicit
`[LIMITATION]` block explaining it (Sections 4.2/5.1 above show it live).
A proper fix is planned for a future version**: it needs (1) a way to
obtain `ε(∞)` from SIESTA — likely its `Optical` module, not currently
automated by any `stb` tool, and (2) a dynamical-matrix-level change in
`core/phonon_workflow.py` implementing the Gonze-Lee subtraction before
diagonalization. Until then, treat a `BULK`-path polar crystal's reported
Γ frequency as a rough `(TO, LO)` bracket midpoint — useful for symmetry,
activity, and relative-mode-ordering comparisons, but not a precise
quantitative TO (or LO) value.

## 7. Known, deliberate limitations

- **No non-analytic term correction (NAC) for `BULK`-path polar
  crystals** — Section 6.4, the main finding of this tutorial's own
  worked example. Not yet implemented; planned for a future version.
- **1D wires stay on the plain dipole-difference path**, not a HYBRID
  treatment (Section 1.5) — likely has the same unverified limitation a
  2D slab's in-plane axes had before the HYBRID fix.
- **`--experimental` needs a genuine broadened/digitized spectrum**, not a
  bare literature peak position (Section 5.3) — this example's own
  literature comparison (Section 6.3) uses a plain reference-line plot
  instead, not the `--experimental` flag.
- **Harmonic only; no anharmonic linewidths/lifetimes/thermal shift** —
  same as every phonon tool in this suite.
- **Supercell size, `--symprec`, `-d`, `--kgrid-density` are all
  user-chosen, never auto-converged** — same limitation `4.4-phonons`/
  `4.11-raman` document for themselves.
- **A too-tight Stage-2 `--symprec` (tighter than Stage 1's) can degrade
  the force-constant reconstruction** — same sharp-edged risk documented
  in `4.4-phonons` README Section 4.6; pass the SAME `--symprec` to both
  stages.

## 8. Step-by-step: running this workflow on your own structure

0. **Prerequisite**: an already-relaxed structure (`structure.fdf`) and
   the SIESTA `calc.fdf` template you relaxed it with (any template is
   fine — it gets forced to single-point per displacement automatically,
   Section 3.2).
1. `stb-ir -s structure.fdf -c calc.fdf -p dojo -dim <a> <b> <c>` — pick a
   supercell large enough that atoms stop interacting across periodic
   images (same rule of thumb as `4.4-phonons`).
2. Run SIESTA in every `ir_study/phonon_disp/disp-*/` folder.
3. `stb-irModes -dir ir_study -p dojo -c calc.fdf` — reads `[NOTE]` from
   Stage 1 (or the `[Path]` line here) to know which folder(s) come next.
4. Run SIESTA in whatever Stage 2 wrote (`born_charge_disp/equilibrium/`
   and/or `dipole_disp/mode_NN_*/`, depending on the path) — remember
   `-np 1` if a `BULK`/`HYBRID` equilibrium run is involved (Section 4.3).
5. `stb-irAnalysis -dir ir_study --save-gnuplot` — read the `[LIMITATION]`
   block if your structure is a polar `BULK`-path crystal (Section 6.4)
   before treating the reported frequency as precise.

## 9. Files in this folder

- `structure.fdf` / `calc.fdf` — the real, relaxed NaCl primitive cell and
  SIESTA calc template used throughout this README and the script.
- `experimental_nacl_TO.dat` — the literature TO reference point (Section
  6.3), NOT used via `--experimental` (Section 5.3 explains why) but kept
  here as the citable, versioned source of the `164.0` value the script's
  own plot draws.
- `real_data/` — genuine SIESTA output (2× `siesta.FA` force files, 1×
  `siesta.BC` + `calc.out` from the Born-charge equilibrium run) from
  running this exact structure this same session, verified
  byte-for-byte reproducible against a fresh Stage 1/2 run (same
  structure + same flags ⟹ identical displaced geometries). The script
  copies these into freshly-generated folders instead of requiring a live
  SIESTA install to run the tutorial (Section "Injecting REAL..." in the
  script itself explains the exact commands you'd run instead, including
  the `-np 1` Cholesky workaround, Section 4.3).
- `example_4.12.sh` — the guided walkthrough script.
- `.gitignore` — excludes `output/` and other regenerable files.

## 10. Running the script

```bash
cd examples/4.12-ir
./example_4.12.sh
```

Safe to re-run any time — it wipes and rebuilds its own `output/` first.
Needs a real `stb_suite` install (`pip install -e .` from `stb-suite/`,
see the top-level `examples/README.md`) and `dojo` pseudopotentials
available (bundled). Does **not** need a working SIESTA binary — Stages
2-3 run against `real_data/`'s genuine, pre-computed SIESTA output
(Section "Files in this folder" above).

## 11. What's next

- Try `--use-symmetry` at Stage 2 to have `stb-irModes` skip probing
  symmetry-forbidden modes entirely, or `--skip-degenerate` to compute
  only one representative per degenerate group (Section 1.4) and derive
  the rest by symmetry, saving real SIESTA time on a larger structure.
- Try a non-polar bulk crystal (Si, diamond) to see the `BULK`-path
  frequency come back with NO Section 6.4 caveat needed (`Z*≈0` there).
- Try a 2D slab (h-BN, graphene) to exercise the `HYBRID` path directly
  (Section 1.3) — no dedicated example folder exists for it yet, but the
  same `stb-ir`/`stb-irModes`/`stb-irAnalysis` commands apply unchanged.
- Run `stb-raman`/`stb-ramanModes` on this SAME `structure.fdf` (no
  dedicated `4.11-raman` example folder exists yet, but the tool works
  unchanged) — Raman and IR are mutually exclusive by symmetry for a
  centrosymmetric crystal (Section 1.2), so NaCl's own `T1u` mode should
  come back Raman-**silent**, a genuine, checkable cross-validation of both
  workflows' symmetry classification against each other.
