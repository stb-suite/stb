# Examples

Hands-on, readable tutorials aimed at people **learning to use** stb_suite —
as opposed to `test/`, which is the developers' smoke-test suite (edge
cases, expected errors, PASS/FAIL counts). Each folder here is a small
lesson: what the tool does, why you'd reach for it, some of the theory
behind why it works the way it does, and a commented script you run and
read at the same time.

This is an ongoing effort — new examples are added a few at a time, one per
tool, not all at once. Not every `stb-*` tool has one yet.

## How the suite itself is organized

`stb-suite`, the interactive front-end, groups every tool into 7
categories:

| # | Category         | What it's for                                                                                        |
|---|------------------|-------------------------------------------------------------------------------------------------------|
| 1 | Inputs           | Tools to set up a SIESTA run (input file, pseudopotentials, k-grid, k-path, DFT+U, fetch, ML pre-relax) |
| 2 | Analysis         | Tools to analyze simulation results (bands, DOS, structures, charge density, charges...)              |
| 3 | Utils            | Helper tools for file management and conversion (translate, clean, cube, status, archive...)          |
| 4 | Structures       | Tools to build, generate, or transform structure files (stacking, supercells, slabs, defects, SQS...) |
| 5 | Workflow         | Complete prep + analysis pipelines for a specific property (strain, elastic constants, phonons, ...)  |
| 6 | Catalysis        | Catalysis workflows (HER, OER); their tools come from an extension package, not from this repository |
| 7 | ML Simulations   | Train a MACE force field on SIESTA data (7.1), and run MD and more with it or the foundation model    |

Each example folder here is named after the **dotted code** the interactive
`stb-suite` menu uses to jump straight to a tool — `<category>.<item>`, and
one level deeper for a workflow's stages (e.g. `5.1.2`). So
`1.1-stb-inputfile/` corresponds to typing `1.1` at the `stb-suite` main
prompt (category `1` = Inputs, item `1` = Input File Generator), or running
`stb-inputfile` directly from the command line — both paths are shown side
by side in every example.

To see the full, current list of codes yourself (this file's own Index
below only lists the tools that already have a written example), run
`stb-suite` and either browse the category menus, or read the
`INPUT_TOOLS`/`ANALYSIS_TOOLS`/`UTILITY_TOOLS`/`STRUCTURE_TOOLS`/
`WORKFLOW_TOOLS`/`CATALYSIS_TOOLS`/`MLSIM_TOOLS` dictionaries directly in
`stb-suite/src/stb/stb_suite.py` — that source is the single source of
truth for the code assigned to any tool.

## Recommended input layout: `calc.fdf`, `structure.fdf`, `config_extra.fdf`

> **Recommended, not required.** Most tools also accept one self-contained `.fdf`, but the
> examples, and every workflow, are built around this split:
>
> | File | Holds | Written by |
> |------|-------|------------|
> | `structure.fdf` | the geometry only: atoms, species, cell | you, or a structure tool (`stb-translate`, `stb-supercell`, `stb-slab`, ...) |
> | `calc.fdf` | the level of theory (`SystemLabel`, basis, XC, mesh cutoff, k-grid, SCF/MD) and `%include structure.fdf` | you, or `stb-inputfile` (1.1) |
> | `config_extra.fdf` | what a workflow stage must impose: single point, fixed cell, spin, dipole, D3, a matched k-grid, ... | the workflow, never you |
>
> A workflow copies `calc.fdf` unchanged into each folder, rewrites only `structure.fdf`, and
> includes `config_extra.fdf` at the **top** of the copy: SIESTA keeps the first value it reads for
> a label, so the forced settings win and your template stays as you wrote it. The level of theory
> is defined once, by you, for every folder of a study.

The full manual, with every tool, workflow and option and these guides, is also available as a
[PDF](https://stb-suite.github.io/stb/assets/stb_manual.pdf).

## Two ways to run any tool

Every example shows both:

1. **Direct CLI** — call the `stb-*` command yourself with flags. Faster,
   scriptable, and what you'd use in a script or once you already know the
   tool's flags by heart.
2. **Interactive `stb-suite` menu** — run `stb-suite` and either navigate
   the category menus or type the tool's dotted code (e.g. `1.1`) straight
   from the main prompt. It walks you through the same choices as guided
   questions instead of flags — a good way to discover what a tool can do,
   and what its defaults are, before memorizing its options.

Both paths call the exact same underlying tool and produce the exact same
output — every example script proves this directly, by driving both paths
with the same inputs and diffing (or otherwise comparing) the result.

## What's inside an example folder

Every example follows the same shape:

- **`README.md`** — what the tool does, when you'd reach for it, and any
  relevant theory (increasingly common for the newer ML-driven tools:
  what the underlying method/model actually does, not just how to invoke
  it). Ends with a flag reference and pointers to what to try next.
- **`example_<code>.sh`** — the guided walkthrough script. Not an automated
  test (see the matching folder under `test/` for that); it runs real
  commands, one self-contained case at a time, into its own
  `output/<case>/` folder, and prints the specific piece of output that
  proves what just happened. Pauses between sections (`[Press Enter to
  continue]`) so you can read before moving on — safe to re-run any time,
  it always starts by wiping its own `output/`.
- **Small input fixtures** (`*.fdf`, occasionally other formats) — checked
  into git, purpose-built for that example (not shared with `test/`'s own
  fixtures, even when conceptually similar), usually with a comment inside
  explaining what's deliberately unusual about them (a wrong lattice
  constant, a too-close pair of atoms, a specific symmetry, ...).
- **`.gitignore`** — excludes `output/` and any files the tool itself
  generates on every run (`references.bib`, a `stb_<tool>_report.txt`,
  convergence plots, ...), so a fresh clone only ever has the tutorial
  itself, not its own regenerable output.

## Conventions shared across examples

Several conventions recur across the newer input-generating tools (see
`1.1-stb-inputfile/`, `1.2-stb-pseudo/`, `1.5-stb-dftu/`, `1.6-stb-fetch/`, `1.7-stb-mlrelax/`,
`4.1-stb-2Dstacking/`, `4.2-stb-supercell/`, `4.3-stb-slab/`,
`4.4-stb-nanotube/`, `4.5-stb-defect/`, `4.6-stb-sqs/`, `4.7-stb-unitcell/`,
`4.8-stb-crystalbuilder/`, `4.9-stb-passivate/`, `4.10-stb-molecule/`,
`4.11-stb-amorphize/`, `4.12-stb-crystalcast/`, `2.1-stb-bands/`, `2.2-stb-dos/`,
`2.3-stb-convdos/`, `2.4-stb-structural/`) and are worth knowing once, rather than
re-discovering per example:

- **A numbered, sectioned report** (`[0] RUN METADATA`, `[1] ...`, ...)
  printed to the console and, with `--save-report`, also saved to
  `stb_<tool>_report.txt` — the same report either way, `--save-report`
  just also persists it to disk.
- **`references.bib`** — written automatically, no flag needed, whenever a
  tool's output depends on a citable method (SIESTA itself, a specific
  correction/potential, a database it queried, ...). Merges with whatever
  `references.bib` is already in the same folder instead of overwriting
  it, so running several tools into one working directory accumulates one
  combined bibliography.
- **Structure validation with an explicit checklist** — atom proximity,
  lattice handedness, atomic density (and a couple of tool-specific extras)
  each get their own reported `[OK]`/`[WARNING]`/`[SKIPPED]` row, not just
  a pass/fail summary, so it's always clear exactly what was checked.
- **`--view`** — an optional, off-by-default flag that opens the structure
  in ASE's interactive 3D viewer right before the tool exits (needs a local
  display, or `ssh -X`/`-Y`). Never exercised by the example scripts
  themselves (an unattended walkthrough shouldn't pop a GUI window), but
  every README mentions the exact command to try it yourself.

## Prerequisites

Install the package in editable mode before running any example:

```bash
cd stb-suite
pip install -e .
```

This puts every `stb-*` command (and `stb-suite` itself) on your PATH.

A few examples (`1.7-stb-mlrelax/`, `4.11-stb-amorphize/`, `4.12-stb-crystalcast/` and the
`7.x-*` ML Simulations examples) additionally need the optional `ml` extra
— PyTorch + `mace-torch`:

```bash
pip install stb_suite[ml]
```

Each such example's own script checks for this up front and exits cleanly
with an install hint if it isn't available, rather than failing partway
through.

## Index

| Code | Tool | Folder |
|------|------|--------|
| 1.1 | `stb-inputfile` | [`1.1-stb-inputfile/`](1.1-stb-inputfile/) |
| 1.2 | `stb-pseudo` | [`1.2-stb-pseudo/`](1.2-stb-pseudo/) |
| 1.3 | `stb-kgrid` | [`1.3-stb-kgrid/`](1.3-stb-kgrid/) |
| 1.4 | `stb-kpath` | [`1.4-stb-kpath/`](1.4-stb-kpath/) |
| 1.5 | `stb-dftu` | [`1.5-stb-dftu/`](1.5-stb-dftu/) |
| 1.6 | `stb-fetch` | [`1.6-stb-fetch/`](1.6-stb-fetch/) |
| 1.7 | `stb-mlrelax` | [`1.7-stb-mlrelax/`](1.7-stb-mlrelax/) |
| 2.1 | `stb-bands` | [`2.1-stb-bands/`](2.1-stb-bands/) |
| 2.2 | `stb-dos` | [`2.2-stb-dos/`](2.2-stb-dos/) |
| 2.3 | `stb-convdos` | [`2.3-stb-convdos/`](2.3-stb-convdos/) |
| 2.4 | `stb-structural` | [`2.4-stb-structural/`](2.4-stb-structural/) |
| 2.5 | `stb-symmetry` | [`2.5-stb-symmetry/`](2.5-stb-symmetry/) |
| 2.6 | `stb-bader` | [`2.6-stb-bader/`](2.6-stb-bader/) |
| 2.7 | `stb-workfunction` | [`2.7-stb-workfunction/`](2.7-stb-workfunction/) |
| 2.8 | `stb-density` | [`2.8-stb-density/`](2.8-stb-density/) |
| 2.9 | `stb-xrd` | [`2.9-stb-xrd/`](2.9-stb-xrd/) |
| 2.10 | `stb-fatbands` | [`2.10-stb-fatbands/`](2.10-stb-fatbands/) |
| 2.11 | `stb-stm` | [`2.11-stm/`](2.11-stm/) |
| 2.12 | `stb-wfdensity` | [`2.12-wfdensity/`](2.12-wfdensity/) |
| 2.13 | `stb-sts` | [`2.13-sts/`](2.13-sts/) |
| 2.14 | `stb-coop` | [`2.14-coop/`](2.14-coop/) |
| 2.15 | `stb-ipr` | [`2.15-ipr/`](2.15-ipr/) |
| 2.16 | `stb-effmass` | [`2.16-effmass/`](2.16-effmass/) |
| 2.17 | `stb-spintexture` | [`2.17-spintexture/`](2.17-spintexture/) |
| 2.18 | `stb-aimdAnalysis` | [`2.18-aimdanalysis/`](2.18-aimdanalysis/) |
| 2.19 | `stb-nativecharges` | [`2.19-stb-nativecharges/`](2.19-stb-nativecharges/) |
| 3.1 | `stb-translate` | [`3.1-stb-translate/`](3.1-stb-translate/) |
| 3.2 | `stb-clean` | [`3.2-stb-clean/`](3.2-stb-clean/) |
| 3.3 | `stb-cube` | [`3.3-stb-cube/`](3.3-stb-cube/) |
| 3.4 | `stb-siesta2wtb` | [`3.4-stb-siesta2wtb/`](3.4-stb-siesta2wtb/) |
| 3.5 | `stb-ani2traj` | [`3.5-stb-ani2traj/`](3.5-stb-ani2traj/) |
| 3.6 | `stb-status` | [`3.6-stb-status/`](3.6-stb-status/) |
| 3.7 | `stb-archive` | [`3.7-stb-archive/`](3.7-stb-archive/) |
| 4.1 | `stb-2Dstacking` | [`4.1-stb-2Dstacking/`](4.1-stb-2Dstacking/) |
| 4.2 | `stb-supercell` | [`4.2-stb-supercell/`](4.2-stb-supercell/) |
| 4.3 | `stb-slab` | [`4.3-stb-slab/`](4.3-stb-slab/) |
| 4.4 | `stb-nanotube` | [`4.4-stb-nanotube/`](4.4-stb-nanotube/) |
| 4.5 | `stb-defect` | [`4.5-stb-defect/`](4.5-stb-defect/) |
| 4.6 | `stb-sqs` | [`4.6-stb-sqs/`](4.6-stb-sqs/) |
| 4.7 | `stb-unitcell` | [`4.7-stb-unitcell/`](4.7-stb-unitcell/) |
| 4.8 | `stb-crystalbuilder` | [`4.8-stb-crystalbuilder/`](4.8-stb-crystalbuilder/) |
| 4.9 | `stb-passivate` | [`4.9-stb-passivate/`](4.9-stb-passivate/) |
| 4.10 | `stb-molecule` | [`4.10-stb-molecule/`](4.10-stb-molecule/) |
| 4.11 | `stb-amorphize` | [`4.11-stb-amorphize/`](4.11-stb-amorphize/) |
| 4.12 | `stb-crystalcast` | [`4.12-stb-crystalcast/`](4.12-stb-crystalcast/) |
| 5.1 | `stb-strain` / `stb-strainAnalysis` | [`5.1-strain/`](5.1-strain/) |
| 5.2 | `stb-elasticInputs` / `stb-elasticAnalysis` | [`5.2-elastic/`](5.2-elastic/) |
| 5.3 | `stb-cohesive` / `stb-cohesiveAnalysis` | [`5.3-cohesive/`](5.3-cohesive/) |
| 5.4 | `stb-phononsCreate` / `stb-phononsPos` | [`5.4-phonons/`](5.4-phonons/) |
| 5.5 | `stb-convergence` / `stb-convergenceAnalysis` | [`5.5-convergence/`](5.5-convergence/) |
| 5.6 | `stb-xrdsearch` / `stb-xrdrank` | [`5.6-xrd/`](5.6-xrd/) |
| 5.7 | `stb-hubbardu` / `stb-hubbarduAlphas` / `stb-hubbarduAnalysis` | [`5.7-hubbardu/`](5.7-hubbardu/) |
| 5.8 | `stb-adsorb` / `stb-adsorbBsse` / `stb-adsorbAnalysis` | [`5.8-adsorption/`](5.8-adsorption/) |
| 5.9 | `stb-neb` / `stb-nebCycle` / `stb-nebAnalysis` | [`5.9-neb/`](5.9-neb/) |
| 5.12 | `stb-ir` / `stb-irModes` / `stb-irAnalysis` | [`5.12-ir/`](5.12-ir/) |
| 5.13 | `stb-gqca` / `stb-gqcaAnalysis` | [`5.13-gqca/`](5.13-gqca/) |
| 5.14 | `stb-optical` / `stb-opticalAnalysis` | [`5.14-optical/`](5.14-optical/) |
| 5.15 | `stb-eosInputs` / `stb-eosAnalysis` | [`5.15-eos/`](5.15-eos/) |
| 5.16 | `stb-chargediffPrep` / `stb-chargediffAnalysis` | [`5.16-chargediff/`](5.16-chargediff/) |
| 5.17 | `stb-hirshfeldPrep` / `stb-hirshfeldIons` / `stb-hirshfeldAnalysis` | [`5.17-hirshfeld/`](5.17-hirshfeld/) |
| 7.1 | `stb-mlff` / `stb-mlffAnalysis` / `stb-mlffValidate` / `stb-mlffActiveLearning` | [`7.1-mlff/`](7.1-mlff/) |
| 7.2 | `stb-mlmd` | [`7.2-mlmd/`](7.2-mlmd/) |
| 7.3 | `stb-mlphonons` | [`7.3-mlphonons/`](7.3-mlphonons/) |
| 7.4 | `stb-mlelastic` | [`7.4-mlelastic/`](7.4-mlelastic/) |
| 7.5 | `stb-mlneb` | [`7.5-mlneb/`](7.5-mlneb/) |
| 7.6 | `stb-mladsorb` | [`7.6-mladsorb/`](7.6-mladsorb/) |

Workflow items (category `5`, and the MLFF workflow `7.1`) get ONE folder per workflow, covering every
stage (prep [+ perturbation] + analysis) together, instead of one folder
per tool -- `5.1-strain/` covers both stages (`stb-strain`, code `5.1.1`,
and `stb-strainAnalysis`, code `5.1.2`) as a single, theory-first tutorial,
`5.2-elastic/` does the same for `stb-elasticInputs` (`5.2.1`) and
`stb-elasticAnalysis` (`5.2.2`), `5.3-cohesive/` does the same for
`stb-cohesive` (`5.3.1`) and `stb-cohesiveAnalysis` (`5.3.2`),
`5.4-phonons/` does the same for `stb-phononsCreate` (`5.4.1`) and
`stb-phononsPos` (`5.4.2`), `5.5-convergence/` does the same for
`stb-convergence` (`5.5.1`, one fully-relaxed folder per swept
Mesh.CutOff/PAO.EnergyShift/k-grid-density value) and
`stb-convergenceAnalysis` (`5.5.2`, the energy- and relaxed-structure
-convergence report), `5.6-xrd/` does the same for `stb-xrdsearch`
(`5.6.1`) and `stb-xrdrank` (`5.6.2`), `5.7-hubbardu/` extends the
pattern across **three** stages instead of two: `stb-hubbardu` (`5.7.1`,
reference prep), `stb-hubbarduAlphas` (`5.7.2`, perturbation prep), and
`stb-hubbarduAnalysis` (`5.7.3`, the linear-response fit and final
computed `U`), and `5.8-adsorption/` does the same for `stb-adsorb`
(`5.8.1`, site prep), `stb-adsorbBsse` (`5.8.2`, BSSE counterpoise prep,
generated only once Stage 1's sites have actually relaxed), and
`stb-adsorbAnalysis` (`5.8.3`, the adsorption-energy ranking, with and
without the BSSE correction). `5.9-neb/` covers `stb-neb` (`5.9.1`, path
interpolation/prep) and `stb-nebAnalysis` (`5.9.2`, the barrier estimate),
plus `stb-nebCycle`, a CLI-only tool (deliberately not in the interactive
menu, meant for a cluster submission script's loop) that sits between
them for modes 2/3's real-DFT refinement cycles. `5.12-ir/` extends the
pattern across **three** stages like `5.7-hubbardu/`: `stb-ir` (`5.12.1`,
phonon displacement prep), `stb-irModes` (`5.12.2`, Γ-point mode/symmetry
analysis plus whichever dipole/Born-charge displacement folder(s) the
structure's dimensionality needs), and `stb-irAnalysis` (`5.12.3`,
per-mode intensities, degenerate-group combining, and the spectrum).
`5.13-gqca/` covers the 2-stage Generalized
Quasi-Chemical Approximation alloy-thermodynamics workflow: `stb-gqca`
(`5.13.1`, the 3 pair-cluster structures) and `stb-gqcaAnalysis`
(`5.13.2`, the mass-action solve over composition/temperature). Unlike
the earlier workflow examples it ships **real, finished SIESTA runs**
(`reference_results/`: Ga(1-x)In(x)As and Ga(1-x)Al(x)As, DZP and TZP) so its
Stage 2 runs genuine data, and its README includes a convergence
investigation of the mixing energy.
`5.14-optical/` covers the 2-stage Optical Properties workflow: `stb-optical`
(`5.14.1`, one SIESTA folder per polarization direction, with the SCF k-grid and the
single-point guards forced through `config_extra.fdf`) and `stb-opticalAnalysis`
(`5.14.2`, the dielectric function and everything derived from it, with sum-rule,
symmetry and Kohn-Sham-gap checks, a scissor what-if and a dashboard). It ships real
finished SIESTA runs of bulk silicon (plus an 11-run convergence study) and of
graphene, and its README covers the theory, the numerical parameters and the
limitations of independent-particle optics.
`7.1-mlff/` covers the 4-stage Custom ML Force Field workflow: `stb-mlff` (`7.1.1`,
training and held-out test configurations, with a foundation-model screen and optional
DFT-D3), `stb-mlffAnalysis` (`7.1.2`, a committee of MACE models fine-tuned on the
SIESTA labels), `stb-mlffValidate` (`7.1.3`, test-set errors with error bars, phonons
against SIESTA, an MD check and a verdict) and `stb-mlffActiveLearning` (`7.1.4`, the
committee picks the next SIESTA runs from MD with the model), linked by
`mlff_workflow.json`. It ships real finished SIESTA runs of graphene (a 3x3 supercell:
two rounds that reach the target, the phonon reference, a k-grid study of the Kohn
anomaly) and of a primitive-cell run where active learning sampled too hard and the
model got worse; its README has the flowchart of the loop.
`5.15-eos/` covers the 2-stage Equation of State workflow: `stb-eosInputs` (`5.15.1`, one SIESTA
folder per volume, the cell scaled isotropically, your `calc.fdf` untouched with the forced
fixed-cell directives in `config_extra.fdf`, the pseudopotentials) and `stb-eosAnalysis`
(`5.15.2`, the fit of `V0`/`E0`/`B0`/`B0'` with leave-one-out uncertainties, SIESTA's own pressure
against `-dE/dV`, a second bulk modulus from the stress, a checks table, a dashboard, gnuplot
files). It ships real finished SIESTA runs of bulk silicon -- 14 setups changing k-points, mesh
and basis one at a time -- and compares them with the all-electron Delta-test reference, which
shows why `B0` alone is a poor figure of merit for a basis set.
`5.16-chargediff/` covers the 2-stage Charge Density Difference workflow: `stb-chargediffPrep`
(`5.16.1`, one SIESTA folder for the combined system and one per fragment, every fragment with the
others as ghost atoms, the same mesh forced in `config_extra.fdf`) and `stb-chargediffAnalysis`
(`5.16.2`, `Δρ` with consistency checks, the counterpoise interaction energy, charge per fragment and
atom, the charge across the plane between fragments, the induced dipole checked against SIESTA's own,
a dashboard and a cube). It ships real SIESTA runs of a hydrogen bond (water dimer), a dative bond
(ammonia borane) and Li on graphene, with the basis-set superposition error measured and the
water dimer's basis, mesh and distance studies.
`5.17-hirshfeld/` covers the 3-stage Hirshfeld-I workflow: `stb-hirshfeldPrep` (`5.17.1`, the system
as a single point with SaveRho plus one isolated neutral atom per element), `stb-hirshfeldIons` (`5.17.2`,
the cation and anion of every element, after the neutral atoms ran) and `stb-hirshfeldAnalysis`
(`5.17.3`, the iterative stockholder charges next to Voronoi, SIESTA's Hirshfeld and simple Hirshfeld,
checks of every reference atom, the error still left when the iteration stops, a dipole check and a
dashboard). It ships real SIESTA runs of H2O, N2O, NaCl and MgO -- two nitrogens of opposite sign, a
crystal with its periodic images, an oxide past +-1 -- and the water molecule's mesh, box and basis studies.
`7.2-mlmd/` covers ML molecular dynamics (`stb-mlmd`): the theory of MD from Verlet to thermostats,
fluctuations, block-averaged errors, `g(r)`/MSD/Lindemann and what classical nuclei leave out; the report
section by section, the model card and the trained-range check. It ships the fine-tuned silicon model of the
`7.1` example with its card, and real GPU runs of 216-512 silicon atoms -- energy conservation vs timestep,
equipartition, finite size, thermal expansion against experiment, a melt far outside the training data --
with SIESTA single points on MD frames to say where the model can be trusted.
`7.3-mlphonons/` covers ML phonons (`stb-mlphonons`): the harmonic crystal, the dynamical matrix, finite
displacements and what a supercell can see, sound velocities as elastic constants, quantum thermodynamics and the
quasi-harmonic thermal expansion with its Grueneisen parameters; the report section by section. It uses the `7.2`
example's fine-tuned silicon model and checks it against neutron, elastic, calorimetric and thermal-expansion data,
against three foundation models (all ~30 % too soft), against the SIESTA phonons of `7.1` -- whose imaginary modes
are traced to a single bad force -- and against the MD of `7.2`; plus supercell, displacement and pre-relax studies.
`7.4-mlelastic/` covers ML elastic constants (`stb-mlelastic`): Hooke's law and the Voigt notation, the stress-strain
method, internal relaxation (clamped vs relaxed ions), the symmetry of the tensor and its independent constants, Born
stability, Voigt/Reuss/Hill moduli and the usual ratios, Christoffel sound velocities, directional moduli and 2D
sheets. With the `7.2` silicon model on a GPU: constants against experiment and against the `7.3` phonon slopes (1 %),
the force threshold of the internal relaxation, the strain amplitude, three foundation models (~60 % soft in shear),
a 216-atom cell, a structure off its minimum, and hcp Mg and graphene with the foundation model. SIESTA runs at the
model's own level of theory (`5.2` workflow) separate the model's error from that of PBE/DZP: shear constants 12–18 %
low, `C12` 15–21 % high, most of the `C44` error in the internal relaxation.
`7.5-mlneb/` covers the ML nudged elastic band (`stb-mlneb`): the potential energy surface and transition-state
theory, the elastic band and why it must be nudged, the improved tangent, the climbing image, IDPP, the
force-consistent spline -- every step computed on the Mueller-Brown surface -- and what a SIESTA single point on a
model's path measures. Its worked example is a real loop: a vacancy hop in 63-atom silicon with the `7.2` model gives
0.72 eV, inside every trained range and with a committee that agrees; seven SIESTA single points on the same path give
0.25 eV; the 14 configurations exported with `--to-mlff` train a round-3 model with 0.31 eV (+23 %), which is also
better on the perfect crystal.

## Adding another example

If you're adding the next one: pick any tool without a folder yet here,
name it `<code>-<tool-name>/` after its `stb-suite` dotted code, and follow
the shape described in "What's inside an example folder" above — the
existing folders (especially the most recently added one) are the best
template to copy from.
