# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

STB-SUITE (Siesta Toolbox Suite) is a collection of independent Python command-line
tools that assist users of the **SIESTA** DFT code (structure format conversion,
k-point/k-path generation, band/DOS/symmetry/elastic/phonon post-processing, Bader
charge analysis, etc.). It is distributed as a single pip/conda package, `stb_suite`,
whose source lives under `stb-suite/src/stb/`.

## Repository layout

- `stb-suite/src/stb/` — the actual package source. Every `*.py` file here is one
  standalone CLI tool.
- `stb-suite/build/` — a `setuptools` build artifact, gitignored/untracked. It is
  **not** regenerated on every `pip install -e .` and can drift from `src/stb` between
  full rebuilds (it has been out of date since at least the phase-7 reorg). Never edit
  here and never treat it as source of truth.
- `stb-suite/pyproject.toml` — package metadata and the `[project.scripts]` table that
  maps each `stb-*` shell command to a module's `main()`.
- `stb-suite/meta.yaml` — conda-forge recipe (mirrors `pyproject.toml` dependencies).
- `test/` — example input files and manual smoke scripts organized by category,
  mirroring the `stb-suite` main menu 1:1 (`1-inputs`, `2-structures`,
  `3-analysis`, `4-workflow`, `5-mlsimulations`, `6-utils`), each mirroring one or more tools.
  `4-workflow/<property>/` has `prep/` + `analysis/` subfolders for the 7
  paired prep+analysis properties (strain, elastic, cohesive, phonons,
  convergence, XRD structure solution, Hubbard U linear response). This is
  **not** an automated test suite (no pytest/unittest, no CI config exists).
  `test/6-utils/1-translator/test.sh` is representative: it generates sample
  structure files in every supported format, runs the built CLI against them, and
  greps the output for expected content.

## Build / install / "test" commands

```bash
# editable install for development (run from stb-suite/)
cd stb-suite
pip install -e .

# after install, tools are on PATH, e.g.:
stb-kgrid --file POSCAR --density 0.2 --type poscar
stb-translate --no-intro -if poscar -i si.poscar -of cif -o si.cif
stb-suite   # interactive menu wrapping all tools
```

There is no unit test suite to run. To validate a change to a tool, exercise it
against the fixtures under `test/` (or the generator scripts, e.g.
`bash test/6-utils/1-translator/test.sh`) and inspect the produced output files,
following the same "generate sample input → run the CLI → check the output file"
pattern used in `test.sh`.

## Architecture

**Each module in `stb-suite/src/stb/` is a thin CLI script** (argparse `main()`,
registered as a console script in `stb-suite/pyproject.toml` under
`[project.scripts]`, e.g. `stb-kgrid = "stb.kgrid:main"`) that delegates shared logic
to `stb-suite/src/stb/core/`:

- `core/structure_io.py` — the only `.fdf` reader/writer (`read_fdf`, `write_fdf`,
  `to_pymatgen`, `rewrite_fdf_lattice`, `rewrite_fdf_positions`, plus thin accessors
  like `lattice_only`, `species_list`/`species_dict`). Use `raw_lattice_vectors()`
  instead of `lattice_only()`/`.lattice` when you're about to write back into an
  existing file's `%block LatticeVectors` while leaving its `LatticeConstant` line
  untouched (see the module docstring) — using the wrong one double-applies the
  lattice constant. `rewrite_fdf_positions(source_path, new_positions, out_path)` is
  the positions analog of `rewrite_fdf_lattice`: replaces only `%block
  AtomicCoordinatesAndAtomicSpecies`, preserving everything else (basis, SCF,
  pseudopotential blocks) verbatim — added for `stb-mlff`, which needs many rattled
  copies of the same reference calculation without hand-reconstructing the whole
  `.fdf` (`write_fdf` only writes a bare-minimum file, losing the original
  calculation's setup).
- `core/siesta_log.py` — parsers for SIESTA `.out` logs: `get_fermi_energy`,
  `get_cell_height`, `get_stress_tensor` (matrix block + Voigt fallback, eV/Å³),
  `get_stress_voigt_kbar` (Voigt-only, raw kBar — kept separate from
  `get_stress_tensor` because their two callers do incompatible downstream math),
  `get_free_energy`, `parse_strain_folder_name`. All return `None` on
  not-found/parse-error rather than raising (these are used in loops that scan many
  `strain_*` folders and tolerate incomplete ones).
- `core/kspace.py` — `compute_monkhorts` (Monkhorst-Pack grid from lattice + target
  density). Raises `ValueError` on zero cell volume.
- `core/cli.py` — `COLORS`, `color_text()`, `show_intro(lines, delay=0.2)` (banner
  content/title is still per-tool, passed in as `lines`), `get_input`/
  `get_float_input`/`get_int_input`.
- `core/deps.py` — `require_sisl()`, the shared `import sisl` guard (prints a
  consistent install hint and exits if missing) used by every tool that reads SIESTA
  grid/Hamiltonian files via `sisl` (`cube.py`, `density.py`, `workfunction.py`,
  `bader.py`, `wantibexos.py`).
- `core/symmetry.py` — `reduce_to_unitcell(structure, mode, symprec, angle_tolerance)`,
  wrapping pymatgen's `SpacegroupAnalyzer` to reduce a structure to its primitive cell,
  conventional cell, or a symmetry-refined version of the input cell (positions snapped
  to the detected symmetry, same cell size). Shared by `stb-unitcell` and `stb-fetch`'s
  `--unitcell` option — extracted here once the second consumer needed it, same
  extract-on-second-use policy as `structure_io.py` below. Note: the output's atom order
  and coordinate origin are never guaranteed to match the input, in any mode — spglib
  rebuilds the cell from the detected symmetry operations from scratch and is free to
  pick any symmetry-equivalent origin, which even tiny input noise can flip for highly
  symmetric structures. Not a bug; the crystal is the same, just relabeled.
- `core/passivation.py` — `passivate_dangling_bonds(structure, passivant, cutoff,
  bond_length)`, caps undercoordinated atoms (e.g. a cut slab's surface) with a
  passivating atom along the missing-bond direction, determined purely from local
  coordination geometry (vector-sum of existing-neighbor directions, negated — no
  bulk-lattice assumptions hard-coded in; verified against a real Si(111) slab to
  recover the exact 109.5-degree tetrahedral angle). Only auto-caps single-missing
  -bond atoms; atoms missing 2+ bonds are geometrically underdetermined from local
  coordination alone and are reported instead of guessed. Shared by `stb-passivate`
  and `stb-slab`'s `--passivate` option.
- `core/mace_relax.py` — `build_cell_mask(vacuum_axes)`, `get_calculator(model,
  device)`, `relax(atoms, calc, cell_mask, ...)`: the MACE-MP-0 load/relax logic
  shared by `stb-mlrelax` and `stb-defect`'s `--ml-rank`. Callers must call
  `core.deps.require_mace()` themselves first — this module only imports
  `mace`/`ase.optimize`/`ase.filters` lazily inside each function, never at module
  level, so merely importing it (e.g. from `defect.py`, which most users run without
  ever touching `--ml-rank`) doesn't force the heavy optional dependency to load.
  `get_calculator`'s `model` argument accepts a custom `.model` file path (e.g. one
  fine-tuned by `stb-mlffAnalysis`) in addition to `"small"/"medium"/"large"` —
  autodetected via `os.path.isfile`, loading `mace.calculators.MACECalculator`
  instead of the `mace_mp` foundation loader. Currently only `stb-mlrelax` exposes
  this as a `--custom-model` CLI flag; the other consumers can gain one the same way
  once there's a real need.
- `core/md_traj.py` — `read_static_lattice`/`read_frame_lattices` (per-MD-step cell
  from `.out`, falling back to `.XV`/`.fdf`), `read_md_timestep_fs` (real MD time
  step from `.fdf`), `unwrap_trajectory` (minimum-image PBC unwrap). Extracted from
  `ani2traj.py` once `stb-aimdAnalysis` became a second consumer needing the exact
  same per-frame lattice/timestep/unwrap logic (same extract-on-second-use policy as
  the rest of `core/`); `ani2traj.py` now imports these instead of defining them
  locally, no behavior change there.

New format-specific structure readers/writers (POSCAR, CIF, XYZ, XSF, FHI, DFTB) still
live only in `translate.py` — it's the sole consumer of those formats, so there's
nothing to extract yet. If a second tool needs one of them, move it into
`core/structure_io.py` first rather than copying `translate.py`'s implementation.

When adding a new tool, mirror this pattern: a new standalone `stb/<name>.py` that
imports what it needs from `stb.core`, plus a new entry in `[project.scripts]` in
`pyproject.toml` (and update `test/` with a fixture if practical).

Newer `2-structures` tools worth knowing about specifically:
- `stb-unitcell` — reduces a structure to its primitive/conventional cell, or refines
  noisy positions to exact symmetry (`--mode primitive|conventional|refined`).
- `stb-crystalbuilder` — the inverse: builds a full structure from a space group +
  the minimal symmetrically-distinct Wyckoff sites, via pymatgen's
  `Structure.from_spacegroup`.
- `stb-defect --all-inequivalent-sites` — auto-enumerates symmetrically distinct sites
  (via `SpacegroupAnalyzer`) instead of requiring the user to pick atom indices by
  hand; writes one output structure per site. `--ml-rank` (needs the optional `ml`
  extra) additionally relaxes each candidate's local geometry with MACE-MP-0
  (positions only, via `core/mace_relax.py`) and prints them ranked by relaxed
  energy — a fast pre-screen for which site to prioritize for real DFT, verified
  against magnetite's 2 distinct Fe sites (tetrahedral vs octahedral): the
  octahedral-site vacancy ranked ~0.06 eV more stable, a physically sensible,
  non-degenerate result.
- `stb-fetch` — the suite's **first network-dependent tool** (everything else is pure
  local file processing). Fetches a structure by exact id or formula search from COD
  (anonymous REST, no key), Materials Project (pymatgen's native `MPRester`, needs a
  free `PMG_MAPI_KEY`), or any OPTIMADE-compliant database (pymatgen's
  `OptimadeRester`, generic client for ~30 known provider aliases or a raw base URL —
  chosen over one-off per-database clients since OPTIMADE is a standardized API many
  databases already implement, e.g. AFLOW, JARVIS, OQMD, and 2D-materials databases
  like `twodmatpedia`). Fetched structures are often disordered even for common
  materials (e.g. oxidation-state-split sites); same-element disorder is collapsed
  automatically, genuine multi-element disorder is rejected with a clear error rather
  than silently producing a wrong structure.
- `stb-passivate` / `stb-slab --passivate` — H-terminates (or any element) dangling
  bonds on a cut surface; see `core/passivation.py` above. Available both as its own
  tool (works on any structure, not just fresh-from-`stb-slab` output) and as a
  same-step convenience flag on `stb-slab`.
- `stb-molecule` — builds an isolated reference molecule (H2O, CO2, benzene, ...) in a
  vacuum box from ASE's bundled G2 database (`ase.build.molecule`, 162 entries, names
  are case-sensitive — `--list` prints them all).
- `stb-mlrelax` — the suite's **first tool with a heavy optional dependency**
  (`pip install stb_suite[ml]`, PyTorch + `mace-torch`; everything else installs with
  the core `dependencies` list alone). Fast pre-relaxation with the MACE-MP-0
  foundation potential before a real SIESTA relaxation -- like `stb-fetch`, needs
  network access on first use (model download, then cached under `~/.cache/mace/`).
  A heuristic, not a DFT replacement. `--relax-cell` auto-adapts to the structure's
  periodicity (reuses `core/kspace.py::detect_vacuum_axes`, same as `stb-kgrid`): a
  vacuum-padded axis (2D slab, 1D wire) is masked out of `FrechetCellFilter`'s Voigt
  strain (`core/mace_relax.py::build_cell_mask`) so only the genuinely periodic
  direction(s) relax and the vacuum thickness stays exactly fixed -- verified against
  a graphene slab, in-plane `a` relaxed to the real ~2.46-2.50 Ang value while the
  vacuum axis changed by exactly 0. A fully isolated structure (0D, vacuum on all 3
  axes) has nothing to relax; `--relax-cell` is a no-op there, positions-only still
  runs.
- `stb-amorphize` — melt-quench amorphous structure generator, the third
  `core/mace_relax.py` consumer: heats a crystalline structure (`ase.md.nptberendsen.
  NPTBerendsen`) above its melting point to erase crystalline memory, then ramps the
  temperature back down, giving a fast heuristic amorphous starting guess. Bulk (3D
  periodic) only — rejects any vacuum-padded axis (melting/NPT-relaxing a slab/wire/
  molecule is physically meaningless). Uses float32 for the MD stages (MACE's own
  guidance: faster, recommended for MD) and a separate float64 calculator for the
  optional final static relax (geometry optimization wants float64). Verified live on
  bulk Si: bond-angle std went from 0.00 deg (perfect crystal) to ~18.6 deg after
  melt-quench while the mean stayed near the tetrahedral angle — the expected
  signature of losing long-range order while keeping short-range coordination,
  matching real amorphous-Si physics.

**`stb_suite.py`** (`stb-suite` command) is the interactive front-end / dispatcher. It
shows a menu with 6 categories, in this order: **1 Inputs, 2 Structures, 3 Analysis,
4 Workflow, 5 ML Simulations, 6 Utils** — backed by six dicts, `INPUT_TOOLS`,
`STRUCTURE_TOOLS`, `ANALYSIS_TOOLS`, `WORKFLOW_TOOLS`, `MLSIM_TOOLS`,
`UTILITY_TOOLS`, each keyed by menu number to `{'title', 'description', 'func'}`.
`INPUT_TOOLS` (category 1) holds only the 3 tools that configure an actual SIESTA
run (input file, k-grid, k-path); `STRUCTURE_TOOLS` (category 2, e.g. `stb-slab`,
`stb-supercell`, `stb-defect`, `stb-crystalcast`) holds everything that builds/
generates/transforms a structure file, split out once `INPUT_TOOLS` grew to 17
entries mixing the two concerns. `MLSIM_TOOLS` (category 5, added alongside
`stb-mlmd`) holds tools that RUN a simulation using a trained MACE potential
(the foundation model or a custom one fine-tuned via `stb-mlffAnalysis`) instead
of driving a real SIESTA calculation — distinct from `WORKFLOW_TOOLS` (generating/
analyzing DFT calculations) and from the one-shot ML structure-preprocessing tools
already in `STRUCTURE_TOOLS` (`stb-mlrelax`, `stb-amorphize`). `WORKFLOW_TOOLS` is
one level deeper (each property has a `'stages'` dict of 2+ entries instead of a
`'func'`); `run_sub_menu()` recurses into `entry['stages']` automatically whenever
it finds one instead of a `'func'`, so no separate menu function was needed for the
extra level. `_flatten_tool_codes()` builds a flat `{"1.1": func, ..., "2.1": func,
..., "4.1.2": func, ..., "5.1": func, ..., "6.4": func}` lookup from these same
dicts at import time (`TOOL_CODES`), letting the main menu's prompt accept a dotted
code (e.g. `4.1.2`) to jump straight to a tool instead of navigating level by level
— regenerate nothing by hand here, it derives entirely from the 6 dicts above.
Every `run_*` function builds an `args: List[str]`
from interactive prompts (`get_input`/`get_float_input`/`get_int_input`) and dispatches
via `run_tool(tool_name, args)`, which shells out to the **installed** console command
by name (`subprocess.run`, must be on `PATH`) and centralizes error handling plus the
"Press Enter to continue" pause. Don't reintroduce a second dispatch path (resolving a
sibling script file via `__file__` and invoking it with `sys.executable`) — that used to
coexist for about a third of the tools for no documented reason and was consolidated
onto `run_tool()` alone.

`stb-ani2traj` (Utils, `ani2traj.py`) converts a SIESTA AIMD trajectory (`<label>.ANI`)
into a multi-frame format that also carries the lattice, for viewing in OVITO/VMD/etc.
`.ANI` itself is a bare multi-frame XYZ (via sisl's `aniSileSiesta`, a subclass of its
`xyzSile`) with no cell information at all — the tool reads all frames via
`sisl.get_sile(ani_file).read_geometry[::stride]()` (an sisl multi-frame slicing API,
not just a plain method call). Output format is one of `xsf` (multi-frame AXSF,
default — read natively by OVITO and VMD's xsf plugin), `pdb` (multi-model with
`CRYST1`, VMD's own default format), or `xyz` (extended XYZ with a `Lattice=` tag,
OVITO-native but VMD's plain xyz reader ignores that tag, so no PBC there) — all three
written via ASE's native multi-frame writers (`ase.io.write` given a list of `Atoms`),
not hand-rolled, since ASE already handles e.g. the triclinic-cell conventions each
format needs. Lattice source, in preference order: per-MD-step `outcell:` blocks in
`<label>.out` (via `core.siesta_log.get_md_trajectory`, added for this — correct even
for variable-cell/NPT/Parrinello-Rahman runs, since it reads the actual cell at every
step instead of assuming one fixed cell), then `<label>.XV`, then `<label>.fdf`. When
the `.out` source is used and `--out-format xyz`, each frame's `E_KS`/`Temp_ion` (also
from `get_md_trajectory`) and real simulation time (`MD.InitialTimeStep`/
`MD.LengthTimeStep`, parsed from `<label>.fdf`) are embedded as extended-XYZ header
properties too — only `xyz` has a slot for this (ASE's xsf/pdb writers don't). `--unwrap`
removes periodic-boundary "jumps" (minimum-image displacement accumulated frame to
frame, same convention as `Atoms.get_all_distances(mic=True)` elsewhere in this suite)
so a molecule that legitimately drifts across the box edge reads as continuous motion
instead of teleporting.

`stb-translate` also writes a `lammps` output format (`writefilelammps`, via
`ase.io.write(..., format='lammps-data')`, `atom_style` implicitly `atomic`) for
handing a relaxed SIESTA structure off to classical MD/force-field tools. It's
output-only (never a valid `--in-format`) since ASE's own LAMMPS reader is for the
dump/trajectory format, not `lammps-data`; deliberately not extended into
`stb-ani2traj`'s format list either, since ASE can only *write* `lammps-data` as a
single frame, not a multi-frame trajectory. Because there's no
`getatomsandvectors_lammps` input reader, the post-write round-trip integrity check
(`reread_for_check`) special-cases `out_format == "lammps"` to re-read the just
-written file directly via `ase.io.read(..., format='lammps-data')` instead of
going through the `getatomsandvectors_* -> build_ase_atoms` path every other format
uses.

`stb-aimdAnalysis` (Analysis, `aimd_analysis.py`) extracts physical quantities from a
SIESTA AIMD trajectory (`<label>.ANI` + `.out`) that `stb-ani2traj` doesn't: radial
distribution function g(r) (vectorized numpy over all frames at once, not one
`ase.Atoms`/frame — matters for long trajectories), mean-squared displacement (MSD)
and a diffusion coefficient (Einstein relation, `D = slope/6` from a linear fit over
a configurable time window), and a vibrational density of states (VDOS) from the
Fourier transform of the velocity autocorrelation function (VACF, Wiener-Khinchin
theorem). `.ANI` carries no velocities, only positions, so velocities are estimated
by central finite difference (`v(t) = (x(t+dt) - x(t-dt)) / (2*dt)`) — the VDOS is
explicitly documented as qualitative, not benchmark-grade, since this amplifies
positional noise. Reuses `core/md_traj.py` for per-frame lattice/timestep/unwrap
(same functions `stb-ani2traj` uses) rather than duplicating that parsing. Analysis
category (not Workflow), since it post-processes one existing trajectory rather than
generating new geometry or aggregating multiple folders (same reasoning as
`stb-bands`/`stb-dos` etc.).

**`stb-mlff` / `stb-mlffAnalysis`** (Workflow item 17, `mlff.py`/`mlff_analysis.py`) is
the suite's first tool that trains a genuinely custom ML potential, rather than just
using the generic MACE-MP-0 foundation model as-is (`stb-mlrelax`, `stb-defect
--ml-rank`, `stb-amorphize`, `stb-adsorb`/`stb-neb --ml-*`, `stb-mlphonons`). Stage 1
(`stb-mlff`) generates `mlff_config_NNN/` training configurations from a reference
`calc.fdf` — rattled displacements (`ase.Atoms.rattle`, one or more `--stdev`
amplitudes) and/or frames sampled from an existing AIMD trajectory via `--from-aimd`
(reusing `core/md_traj.py`) — each folder keeping the reference's exact basis/mesh
cutoff/pseudopotentials via the new `rewrite_fdf_positions`. The user runs SIESTA in
each folder. Stage 2 (`stb-mlffAnalysis`) aggregates the finished calculations
(energy via `core.siesta_log.get_free_energy`, forces via a new all-atom `.FA`
reader generalizing `her_analysis.py::read_fa_force`'s single-atom version, stress
best-effort via `core.siesta_log.get_stress_tensor`). The train/validation split is
done in Python (`split_train_valid`, seeded, writes `train_set.xyz`/`valid_set.xyz`)
rather than left to `mace_run_train`'s own `--valid_fraction`, so the exact held-out
configurations are known afterwards — needed for the parity plot and foundation
-model comparison below (`mace_run_train`'s own docs confirm `--valid_fraction` is
simply ignored once `--valid_file` is given). Then fine-tunes a MACE-MP checkpoint
on it via `mace_run_train` (subprocess, same "wrap the external tool" pattern as
`core/phonon_workflow.py`'s phonopy calls). Deliberately uses plain fine-tuning with
`--E0s average` (per-element atomic
reference energies fit from the SIESTA data itself), NOT
`--multiheads_finetuning`/`--E0s foundation` (mace_run_train's own recommended
anti-catastrophic-forgetting mode) — verified live that `foundation` E0s leaves a
~424 eV/atom constant offset baked into every energy target, because SIESTA's
absolute total energy (its own pseudopotentials/core treatment) is not on the same
scale as the foundation model's own reference DFT code; `average` fixes this since
it's fit from the actual training data, but mace_run_train's `multiheads_finetuning`
asserts against `average` outright, so the choice is: correct energy scale, or
anti-forgetting replay, not both. Verified end-to-end against a real 8-configuration
SIESTA dataset (O2 in vacuum): correct-scale RMSE convergence, a loadable
`.model` file, and a working `mace.calculators.MACECalculator` energy evaluation.
`--export-lammps` additionally compiles the fine-tuned model into a TorchScript file
(`<model>-lammps.pt`, via `mace_create_lammps_model`, also part of `mace-torch`)
loadable by LAMMPS's `pair_style mace` (needs a LAMMPS build with the MACE pair
style, e.g. https://github.com/ACEsuit/lammps -- this only produces the file, it
doesn't need or use a LAMMPS installation itself). Verified live: converts either of
the two `.model` files `mace_run_train` writes (plain and TorchScript-`_compiled`)
non-interactively for a single-head model (multi-head would otherwise prompt for a
head to export).

After training, `stb-mlffAnalysis` evaluates the fine-tuned model on the held-out
`valid_set.xyz` (`evaluate_on_configs`) and, unless `--skip-foundation-comparison`,
also evaluates the raw (non-fine-tuned) foundation model on the exact same
configurations via `core.mace_relax.get_calculator` — reporting RMSE energy/force
for both side by side and overlaying both on one `<name>_parity.png` (predicted vs
DFT, `plot_parity`, y=x reference line), so the improvement from fine-tuning is
directly visible rather than just a number in the training log. `--save-data` writes
the underlying predicted/DFT pairs as `.dat` files.

`stb-mlffActiveLearning` (Stage 3 of the same Workflow item, `mlff_active_learning.py`)
uses an already fine-tuned model to screen a fresh batch of candidate configurations
and writes the `--top-k` ones with the largest *predicted* atomic force as new
`mlff_config_NNN/` folders (continuing the numbering of any already present), for
the user to label with real SIESTA and fold back into another `stb-mlffAnalysis`
round. Two candidate-generation strategies, `--sampling-method`:
- `rattle` (default, cheap): independent static snapshots via `mlff.
  build_rattled_configs`, typically at a larger `--stdev` than the original
  training set.
- `md`: one short NVT trajectory driven by the model itself, reusing `stb-mlmd`'s
  own `build_dynamics` **directly** (not a reimplementation — same physically
  -verified fs/GPa unit conversions) via `sample_md_candidates`. More physically
  representative than rattling (visits the configuration space the system would
  actually explore thermally, instead of a blind perturbation around the static
  reference) — deliberately run at a HIGHER temperature than any real production
  use (`--md-temperature`, default 800 K) purely to cover more configuration space
  per step, a common MLIP training-set-generation trick, not a claim about the
  material's real operating conditions. The predicted max|F| score comes for free
  from each MD step (already computed by the integrator), no extra evaluation per
  frame needed. `stb-mlmd` itself stays unaware of any of this — it only knows how
  to run MD; `stb-mlffActiveLearning` decides to reuse that as a sampling engine,
  matching the same direct-import-reuse pattern already used with `mlff.py`.

Explicitly documented as NOT rigorous uncertainty quantification (that needs a
committee of independently-trained models) — a single-model heuristic: a
configuration where the potential itself predicts an unusually large local force is
a reasonable, cheap signal that the region is poorly represented in the current
training set, not a proof.

**`stb-mlmd`** (ML Simulations, `mlmd.py`) is the first tool in the new ML Simulations
category (5) — tools that RUN a simulation with a trained MACE potential rather than
generating/analyzing a real SIESTA calculation. Runs NVE (`ase.md.verlet.
VelocityVerlet`), NVT (`ase.md.langevin.Langevin`), or NPT (`ase.md.nptberendsen.
NPTBerendsen`, bulk-only — rejects a vacuum-padded axis, same reasoning as
`stb-amorphize`) molecular dynamics with either the MACE-MP-0 foundation model or a
custom fine-tuned one (`--custom-model`, same flag name/semantics as `stb-mlrelax`).
Physical units (`--timestep` in fs, `--friction` in fs⁻¹, `--pressure` in GPa) are
explicitly converted to ASE's internal unit system (`* ase.units.fs`, `/
ase.units.fs`, `* ase.units.GPa`) before being passed to the ASE MD classes, which
all expect ASE's own internal time/pressure units, not fs/GPa directly — verified
against ASE's own docstrings (`timestep: The time step in ASE time units`) after
finding a real, verified ~10x timestep bug in `stb-amorphize`'s own MD call
(`timestep=args.timestep * 1.0` instead of `* units.fs`, silently running its
melt-quench at ~10.18x the documented timestep since `ase.units.fs ≈ 0.0982`) —
fixed there too once found. Trajectory output reuses `stb-ani2traj`'s exact 3-format
choice (xsf/pdb/xyz) for viewer consistency; the final frame is also written as a
bare `.fdf` (same `from_pymatgen`+`write_fdf` pattern already used by `stb-mlrelax`/
`stb-amorphize`, not `rewrite_fdf_positions` — that one's for the mlff workflow's
different need of preserving many configs' full calc setup).

Manually reviewed + physics-verified live (real MACE-MP-0 runs on a 64-atom bulk Si
cell) after shipping, leading to 4 follow-up additions, all confirmed against that
same live data:
- **Total energy tracking**: the report only ever printed `E_pot`, which is NOT the
  NVE-conserved quantity (`E_pot` trades with `E_kin` by design) — verified this was
  misleading by hand-reconstructing `E_total = E_pot + E_kin` from a saved trajectory
  and finding ~2e-6 relative drift over 1000 steps, invisible from the tool's own
  original output. Now tracks and reports `E_total` directly every run (drift %,
  flagged if >1%, for NVE; mean/std `T` — and mean/std cell volume for NPT —
  instead of just the final instantaneous value, for NVT/NPT).
- **`--equilibration-steps`**: verified live that a freshly Maxwell-Boltzmann-seeded
  structure's `E_pot` relaxes for ~100-200 steps before settling into steady
  dynamics; this transient was previously included in the saved trajectory/
  statistics/final structure. Runs and discards N steps before production starts.
- **`<stem>_md_diagnostics.png`**: always-generated 2-panel plot (`E_total` and `T`
  vs time, matplotlib, same visual convention as the rest of the suite) plus
  `--save-data` for the raw numbers. The temperature panel's reference line is
  labeled "initial (NVE has no target)" for the nve ensemble specifically — NVE has
  no thermostat, so `--temperature` there is only the initial-velocity seed, not
  something the dynamics maintains; labeling it "target" (as nvt/npt correctly are)
  would misrepresent that.
- **`--taut`/`--taup`/`--compressibility`**: previously hardcoded in `build_dynamics`
  (100 fs / 500 fs / water's `4.57e-5` eV/Å³ — the same generic placeholder
  `stb-amorphize` uses). Exposed as CLI overrides for a better-calibrated NPT
  barostat response on a real, non-water-like material.

Also closed a real analysis gap found during this review: `stb-aimdAnalysis` could
only read a SIESTA `.ANI`, so an `stb-mlmd` trajectory had no path into the suite's
RDF/MSD/VDOS analysis at all. `stb-aimdAnalysis` now accepts `--trajectory <path>`
(mutually exclusive with `--label`) to read any ASE-readable multi-frame file
directly (`ase.io.read(path, index=':')`) instead of going through `sisl`/SIESTA at
all — every `compute_rdf`/`compute_msd`/`compute_vacf_vdos` function already worked
on plain `frac_positions`/`cells`/`symbols`/`dt_fs` arrays, so only the "read frames
and determine `dt_fs`" input path needed a second branch, not the analysis math
itself. `--dt` (fs between saved frames) is auto-detected from consecutive frames'
`Time` info when reading an extended-xyz written with it (e.g. `stb-mlmd
--out-format xyz`); required explicitly otherwise (xsf/pdb carry no such data).
Output filenames key off a generic `stem` (the SystemLabel, or the trajectory
file's basename) instead of `args.label` directly, so naming works the same in
either input mode.

**`stb-mlphonons`** (ML Simulations, `mlphonons.py`, v2.1.0) is the second tool in the
ML Simulations category — a standalone phonon calculation (displacements, force
constants, band structure, DOS, thermal properties, optionally QHA) driven entirely by
a MACE potential, no SIESTA and no separate `stb-phononsPos` analysis step needed. It
originally coexisted with two SIESTA-Phonons-workflow ML variants,
`stb-phononsML`/`stb-phononsQHA`, which only accepted `small`/`medium`/`large`
foundation checkpoints (no custom fine-tuned model). Once `stb-mlphonons` gained
foundation-model comparison and its own QHA mode (making it a strict superset), the
user had those two **deleted outright** (`phonons_ml.py`/`phonons_qha.py`, their
`pyproject.toml` entries, and their Phonons-workflow menu stages all removed) rather
than kept as parallel paths — `stb-mlphonons` is now the only ML-driven phonon tool in
the suite. Reuses, rather than duplicates, the algorithmic pieces originally extracted
for those two tools: `core.mace_phonons.generate_ml_displacements`/
`compute_force_constants` (this is now their sole consumer — the module docstring
explains why it still lives in `core/` despite that) for the physics, and
`phonons_pos.build_band_path`/`band_path_to_phonopy_format`/`band_tick_positions`/
`pretty_label` for the q-path (ASE's own Bravais-lattice/`Cell.bandpath` machinery,
deliberately NOT phonopy's seekpath-based `auto_band_structure` — `seekpath` isn't
installed/a dependency of this suite; `phonons_pos.py`'s own docstring already explains
this exact tradeoff). Plots (bands/DOS/thermal, matplotlib PNGs) are new code, not
reused from `phonons_pos.py` (which writes gnuplot `.dat`+`.gplot` pairs, the older
convention) — matching the newer matplotlib convention this session's tools
(`aimd_analysis.py`/`mlff_analysis.py`/`mlmd.py`) already use.

Four features on top of the original single-volume calculation:
- **Foundation-model comparison**: when `--custom-model <path>` is given, also runs
  the same pipeline with the MACE-MP-0 foundation model (`--skip-foundation-comparison`
  to disable) and overlays both on the same bands/DOS/thermal plots — an honest,
  at-a-glance check of whether fine-tuning actually changed the phonon physics.
  Verified live: an 8-atom Si model fine-tuned this session showed a visibly different
  top optical branch (~18.2 THz) vs. the foundation model (~11.7 THz) on the same
  structure.
- **Acoustic sum rule (ASR) correction**: calls `phonon.symmetrize_force_constants
  (show_drift=False)` before computing bands/DOS, eliminating spurious near-zero
  negative frequencies at Γ that are numerical noise, not a genuine instability. Even
  after this correction a residual like `-0.0000 THz` can still print, so the
  imaginary-mode check uses a tolerance (`_IMAGINARY_MODE_TOL_THZ = -0.01`, a module
  constant) instead of a strict `< 0`, avoiding a false-positive `[WARNING]` — verified
  live (a real run went from a spurious `-0.0003 THz` warning to a clean "No imaginary
  modes found").
- **Species-projected DOS + symmetry report**: `phonon.run_mesh(..., with_eigenvectors=
  True)` enables both `run_total_dos`/`run_projected_dos`, the latter summed per
  species (`phonon.primitive.symbols`) into `dos.dat`'s extra columns; a symmetry table
  (space group, point group, and the displacement-count reduction from symmetry) is
  printed alongside, adapted from the deleted `phonons_ml.py`. Verified live: a Ge
  substitutional defect in bulk Si showed a physically correct low-frequency PDOS bump
  from the heavier, more weakly-bonded Ge atom.
- **QHA (`--qha`)**: a self-contained volume scan (`--n-volumes`, `--strain-range`,
  `--eos`), adapted directly from the deleted `phonons_qha.py`'s `phonopy.api_qha.
  PhonopyQHA`-based fitting — isotropic strain via `np.linspace`, one full
  displacement+force-constant+thermal-properties pass per volume, then bulk modulus/
  thermal expansion/equilibrium volume vs. temperature. Requires `--n-volumes >= 4`
  (enforced via `parser.error`) since `PhonopyQHA`'s derivatives need at least that
  many points. Mutually exclusive with foundation-model comparison (no QHA-vs-QHA
  overlay yet — documented limitation, not a bug). Verified live with
  `--n-volumes 5 --strain-range 3.0`: qualitatively sensible EOS/thermal-expansion/
  bulk-modulus curves (values themselves noisy, as expected from deliberately toy
  scan parameters, not a code defect).

Verified live on a 64-atom bulk Si:Ge cell for the base pipeline: physically sensible
dispersion (correct acoustic branches to ~0 at Γ, max frequency in the right ballpark)
and textbook Debye-like thermal properties (entropy → 0 as T → 0, heat capacity
saturating at high T).

Three more features on a second pass, plus a real bug fix shared with the SIESTA path:
- **`--freeze-unstable-mode`**: mirrors `stb-phononsPos`'s own feature of the same name
  -- if the DOS mesh has a genuine imaginary mode (below `_IMAGINARY_MODE_TOL_THZ`),
  displaces the supercell along that mode's eigenvector (`phonon.run_modulations`,
  calibrated to the requested `--freeze-amplitude` in Angstrom via a unit-amplitude
  probe displacement first, since phonopy's own "amplitude" parameter is an internal
  eigenvector-normalization scale, not itself in Angstrom) and writes the result as a
  new `.fdf`. Unlike `phonons_pos.py`, no bohr conversion is needed on write since this
  tool's Phonopy object is already built in Angstrom (see `core/mace_phonons.py`).
- **`--animate-mode`**: writes a multi-frame trajectory (xsf/pdb/xyz, same convention as
  `stb-ani2traj`/`stb-mlmd`) animating one normal mode's vibration at a given q-point
  (`--animate-qpoint`, default Gamma) by sweeping `phonon.run_modulations`'s phase over
  a full cycle -- same modulation machinery `--freeze-unstable-mode` uses for a single
  snapshot. `--animate-dim` controls the modulation supercell size and must make the
  q-point commensurate (only the default Gamma point works with the default `1 1 1`).
- **`--check-convergence`**: runs the full single-volume pipeline at multiple
  `--convergence-dims` supercell sizes (default: the given `--dim` and `--dim+1` on
  every non-vacuum axis) and reports how the max frequency / heat capacity / entropy
  change between them -- MACE is cheap enough that this is an actually affordable
  convergence check, unlike the real DFT-based workflow. Mutually exclusive with
  `--qha`/`--freeze-unstable-mode`/`--animate-mode`.

**Real bug found and fixed while verifying `--freeze-unstable-mode` live**: the
amplitude-calibration probe (`mod.modulated_supercells[0].positions -
mod.supercell.positions`, comparing raw Cartesian coordinates) is wrong under periodic
wraparound -- an atom sitting near a cell face can have its wrapped fractional
coordinate flip from ~0.999 to ~0.001 between the unmodulated and modulated supercells
for a real displacement of a fraction of an Angstrom, which reads as a jump of nearly a
full lattice vector. Caught live on a real compressed (deliberately made unstable)
8-atom Si cell: a requested `--freeze-amplitude 0.2` measured as ~74 Ang achieved,
instead of the correct 0.2. Fixed with a minimum-image-convention helper
(`_pbc_atomic_shift`, fractional-difference-then-wrap-then-back-to-Cartesian, same
convention as `core/md_traj.py`'s trajectory unwrapping elsewhere in the suite) used for
both the freeze probe/achieved measurement and `--animate-mode`'s own amplitude
calibration (which needed the exact same probe-and-rescale treatment -- verified live
that passing the requested Angstrom value straight through as phonopy's "amplitude"
gave a real optical mode only ~0.094 Ang of motion for a nominal "1.0", not 1 Ang). The
identical latent bug existed in `phonons_pos.py`'s pre-existing `--freeze-unstable-mode`
(same raw-Cartesian-difference pattern) and was fixed there too with the same helper.

**`stb-mlelastic`** (ML Simulations, `mlelastic.py`, v1.1.0) is the third tool in the ML
Simulations category -- a standalone stiffness-matrix (stress-strain) calculation driven
entirely by a MACE potential, the same relationship `stb-mlphonons` has to
`stb-phononsCreate`/`stb-phononsPos`. Unlike that pair, it's a single self-contained
tool with no `strain_*/` folders to hand off to SIESTA and reload later: MACE stress is
cheap enough to evaluate in-memory for every canonical Voigt direction in one run, so
there's no DFT-cost reason for `stb-elasticInputs`' symmetry-based direction-reduction,
and no `--method energy` path either (that method exists in the DFT tool specifically to
sidestep SIESTA's real-space-grid "eggbox effect" on stress, which MACE has no analog
of). Reuses, rather than duplicates, `elastic_analysis.py`'s post-fit numerics --
`compute_stiffness_matrix`, `direction_fit_diagnostics`, `tensor_symmetry_check`,
`check_stability_and_report`, and `emit_elastic_report` (renamed from the original
`_emit_elastic_report` once this became a second consumer; it only reads a handful of
attribute names off any caller's argparse `Namespace`, so it doesn't care whether the
data came from SIESTA or MACE) -- since a stress-strain elastic-constant fit is the same
linear algebra regardless of data source; only the data-mining step (reading `calc.out`
per folder there) differs, replaced here by `compute_ml_elastic_data`'s in-memory MACE
evaluation loop. Also reuses `elastic_inputs.py`'s `get_strain_matrix`/`direction_axes`
for the identical deformation-matrix/vacuum-blocking convention, so a strain labeled
`xy` means the same deformation in both tools. Plots (`plot_stress_strain`) are new
matplotlib code, not `elastic_analysis.py`'s gnuplot `.dat`+`.gplot` writer -- matching
the newer matplotlib convention this session's ML tools already use. A **real bug** was
found and fixed while wiring the report reuse: `emit_elastic_report`'s closing "Full
report" line read `elastic_analysis.py`'s own hardcoded `REPORT_FILE` module constant
directly, so it always printed `mechanical_properties.txt` regardless of what
`stb-mlelastic` actually wrote -- fixed by adding a `report_path` parameter (defaulting
to the original constant, so `elastic_analysis.py`'s own 2 call sites are unaffected).
Verified live: an 8-atom bulk Si cell recovered cubic symmetry (C11=C22=C33,
C12=C13=C23, C44=C55=C66, all independently fit with no symmetry constraint imposed)
and passed the Born stability check; a 2-atom graphene monolayer (vacuum along c)
correctly dropped the vacuum-blocked zz/yz/xz directions, applied the same Lz dilution
correction `stb-elasticAnalysis` needs for SIESTA's own volumetric stress, and gave a
layer modulus (~264-318 N/m) in the right physical ballpark against the experimental
~340 N/m.

A second pass added 4 more features:
- **Foundation-model comparison**: when `--custom-model` is given, also runs the same
  pipeline with the MACE-MP-0 foundation model (`--skip-foundation-comparison` to
  disable) and overlays both on the stress-strain plot plus a side-by-side diagonal
  -constants table -- same pattern as `stb-mlphonons`/`stb-mlffAnalysis`.
- **Sound velocities + Debye temperature** (3D only): Vp/Vs/Vm and theta_D derived from
  the Voigt-Reuss-Hill bulk/shear moduli `check_stability_and_report` already computes,
  plus the structure's own mass density -- zero extra MACE cost. Verified live on bulk
  Si: density came out 2285 kg/m^3 (real Si: 2329), and the ~391 K Debye temperature
  (real Si: ~645 K) is proportionally consistent with MACE-MP-0's own known
  underestimate of Si's elastic constants (~95.6/57.8/26.8 GPa vs. the ~166/64/80 GPa
  experimental C11/C12/C44) -- an internal-consistency check, not an independent error.
- **Directional Young's-modulus anisotropy surface** (3D only, `--skip-anisotropy-surface`
  to disable): builds the full compliance tensor from the fitted C_sym
  (`compliance_voigt_to_full`, its own small helper -- NOT `core/symmetry.py`'s
  `voigt_to_full_tensor`, which uses the STIFFNESS engineering-strain factor convention,
  not compliance's) and evaluates `1/E(n) = s_ijkl n_i n_j n_k n_l` over a spherical grid,
  plotted as a 3D "material property surface" (ELATE-style). Verified numerically
  against a synthetic isotropic tensor (C44=(C11-C12)/2): E(n) constant to machine
  precision across arbitrary directions and matching the closed-form isotropic formula
  E=(C11-C12)(C11+2*C12)/(C11+C12) exactly.
- **`--check-convergence`**: re-fits at multiple `--max` strain magnitudes
  (`--convergence-strains`) and reports how each diagonal constant changes -- checks the
  linear-regime assumption itself, mutually exclusive with foundation comparison.

**Real bug found and fixed while testing `--dimensionality` manual override**: forcing
a lower-dimensional structure (e.g. 2D graphene) to report as 3D leaves several
canonical directions (zz/yz/xz) with no data regardless of the override (direction
selection is vacuum-based, unaffected by `--dimensionality`), making C_sym singular --
`plot_anisotropy_surface`'s `np.linalg.inv` crashed uncaught, and the sound-velocity
section would have silently returned a physically meaningless number (its Reuss-average
inversion already had a fallback, masking the same underlying rank-deficiency). Fixed
by checking `missing_directions` before attempting either section (prints an honest
`[WARNING] Skipped` instead) and additionally wrapping `plot_anisotropy_surface`'s
inversion in `try/except LinAlgError` as defense in depth.

**`stb-mlsearch`** (ML Simulations, `mlsearch.py`) is the fourth tool in the ML
Simulations category, and the first WITHOUT a direct analog in the SIESTA-based
Workflow tools: unlike `stb-mlphonons`/`stb-mlelastic` (fast ML previews of an existing
DFT prep+analysis pair), it explores configuration space itself (which atomic
arrangement, at a fixed cell, has the lowest energy) -- a genuinely new capability, not
a cheaper version of something the suite already does. Wraps `ase.optimize.basin.
BasinHopping` directly rather than reimplementing the algorithm: perturb (uniform
random atomic displacement, `--dr`), locally relax with MACE (`--fmax`), accept/reject
via the Metropolis criterion (probability `exp(-dE/kT)`, `--temperature` in Kelvin,
converted to eV via `ase.units.kB`), track the best locally-relaxed structure found
across `--steps` attempts. ASE's own implementation only perturbs atomic POSITIONS,
never the cell shape/volume (written originally for clusters, not crystal structure
prediction) -- appropriate for this tool's primary use case (finding the best atomic
arrangement at a fixed cell: defect/vacancy site occupation, disordered/alloy
configurations, local amorphous-like rearrangements), not full crystal structure
prediction. An optional final full relax (`--no-final-cell-relax` to disable) polishes
the winning configuration's cell shape too, via `core/mace_relax.py`'s
`build_cell_mask` (same vacuum-aware convention as `stb-mlrelax`/`stb-amorphize`).
`parse_basin_log` reads numbers straight out of `BasinHopping`'s own logfile (a stable,
fixed-format text line per step) for the `search_history.png` diagnostic plot, rather
than subclassing/monkeypatching the class to intercept them. The full search
trajectory (every attempted, locally-relaxed candidate, from `BasinHopping`'s own
`local_minima_trajectory`) is converted to the same xsf/pdb/xyz convention `stb-
ani2traj`/`stb-mlmd`/`stb-mlphonons` already use. Explicitly NOT rigorous global
optimization (a heuristic search, no proof of finding the true global minimum) and
explicitly NOT crystal structure prediction (cell shape untouched during the search
itself). Verified live: `--seed` gives bit-identical results across repeated runs; a
heavily rattled (0.3 Ang/atom stdev) 8-atom Si cell always converges back to
essentially the same low energy a plain local relax alone finds (diamond Si has one
dominant global minimum at this perturbation scale) -- the correct, expected result for
this simple single-basin test system, not a sign the search did nothing.

A second `--algorithm simulated-annealing` mode was added: one continuous NVT MD
trajectory (Langevin, reusing `mlmd.py`'s own `build_dynamics` -- same import-don't
-duplicate pattern `stb-mlffActiveLearning`/`stb-mlmelting` already use) whose target
temperature is cooled from `--temp-start` to `--temp-end` over `--steps` total MD steps
(`--schedule exponential`, the classic SA convention, updated on the fly via
`Langevin.set_temperature`, or `linear`), snapshotting and locally relaxing the running
configuration every `--snapshot-interval` steps to find candidate minima along the
cooling path -- the annealing analog of what `BasinHopping`'s own per-step Metropolis
+relax loop does for random perturbations, just driven by a physically continuous
cooling trajectory instead of independent random jumps. `run_basin_hopping`/
`run_simulated_annealing` share the same `(best_energy, n_new_minima, traj_out)` return
shape so `main()`'s final-polish/summary reporting tail needs no per-algorithm branching.
Verified live: the exponential schedule's exact per-step temperature was confirmed by
hand against the closed-form geometric-decay formula (`T_start * (T_end/T_start)
**frac`); both algorithms converge to essentially the same low energy on the same
heavily-rattled 8-atom Si test case.

**`stb-mlmelting`** (ML Simulations, `mlmelting.py`) is the fifth tool in the ML
Simulations category -- brackets a bulk material's melting point via a sequence of
short MACE-driven NVT MD runs (one per candidate temperature, reusing `mlmd.py`'s own
`build_dynamics` directly, same import-don't-duplicate pattern `stb-mlffActiveLearning`
already uses), each starting fresh from the same MACE-relaxed solid reference. Tracks
the Lindemann index (RMS atomic displacement from each atom's own time-averaged
position -- not a fixed lattice site, so rigid drift doesn't register as melting --
divided by the nearest-neighbor distance) vs. temperature: stays small and roughly flat
for a vibrating solid, rises sharply once atoms are no longer confined to a lattice
site. The classic Lindemann (1910) melting criterion, with the literature's usual
~0.10 threshold (`--lindemann-threshold`) used to bracket the transition by linear
interpolation between the two bracketing scanned temperatures. Also reports the
self-diffusion coefficient D(T) as a cross-check, reusing `aimd_analysis.py`'s own
`compute_msd`/`fit_diffusion_coefficient` (the identical MSD/Einstein-relation fit
already used for a SIESTA AIMD trajectory or an `stb-mlmd` run) -- D is clipped to 0 if
the linear fit lands slightly negative, a real but understood statistical artifact of a
short/noisy MSD fit deep in the near-zero-true-diffusion solid regime, not a sign of
unphysical backward diffusion. Bulk (3D periodic) only, same reasoning as
`stb-amorphize`/`stb-mlmd`'s NPT mode. Explicitly documented (both in the module
docstring and printed with the estimate itself) as a coarse "one-phase" heating method,
known to SUPERHEAT past the true thermodynamic melting point since a small, defect-free
periodic cell has no free surface to nucleate melting from -- verified live on a 4-atom
Al fcc cell (real Al: ~933 K): the tool's own estimate landed at ~3095 K, a large but
expected overestimate given the tiny cell and short runs used, not a bug.

**`stb-mlconvergence`** (ML Simulations, `mlconvergence.py`) is the sixth tool in the ML
Simulations category -- a model-size convergence check, not a physics-property
calculation: runs the exact same reference calculation (a full relax, cell + positions,
via `core/mace_relax.py`'s `relax`/`build_cell_mask`) independently with each requested
MACE-MP-0 size (`--models`, default `small medium large`) and any `--custom-models`
given, then reports how much the answer (energy/atom, cell volume, wall-clock time)
actually changes between consecutive sizes -- flagged NOT CONVERGED if a pair differs
by more than `--energy-tolerance`/`--lattice-tolerance`. Helps decide whether a bigger,
much slower model is worth it for a given structure, instead of guessing; explicitly
NOT a check against real DFT, only whether MACE-MP-0's own answer has converged with
size. Verified live on an 8-atom bulk Si cell: small vs. medium differ by ~29 meV/atom
and ~0.6% in cell volume (correctly flagged NOT CONVERGED at the tight 5 meV/atom, 1%
defaults; correctly flagged OK and recommending the smaller/faster model at looser
tolerances) -- both differences are physically small but real, exactly the kind of
signal this tool exists to surface before committing to a size for a larger workflow.

**`stb-mlneb`** (ML Simulations, `mlneb.py`) is the seventh tool in the ML Simulations
category -- a standalone climbing-image NEB (nudged elastic band). `stb-neb` already
has a `--ml-neb` preview mode (a real climbing-image NEB on MACE-MP-0 before committing
to the real SIESTA-based `image_NN/` folders that tool generates), but that mode is
bolted onto a DFT-oriented workflow: it always needs `--calc`/pseudopotentials even
though the ML path never touches either, and its `--ml-model` only accepts
`small`/`medium`/`large` (no `--custom-model`, unlike every other ML tool in this
suite). `stb-mlneb` is the pure-MACE extraction of that same capability for users who
only want the ML-level barrier estimate -- e.g. defect/vacancy migration, an adsorbate
hopping between sites -- with no DFT step at all. Deliberately reuses (not duplicates)
`neb.py`'s own pure-geometry helpers (`check_composition_match`, `wrap_into_cell`,
`resolve_lattice_mismatch`, `linear_interpolate_images`, `idpp_refine_images`,
`compute_frozen_indices`, `cumulative_reaction_coordinates`, `check_path_quality`,
`write_path_trajectory`, `write_ml_preview_plot`) -- interpolating/quality-checking a
reaction path is exactly the same problem whether the images end up going to SIESTA or
straight into a MACE NEB; only `write_image_folder` (SIESTA input/pseudopotential
-specific) isn't reused. The NEB relaxation itself reuses `core/mace_relax.py`'s
`relax_neb` (the same function `stb-neb --ml-neb` already calls). `--save-images`
writes each relaxed image as a bare `image_NN.fdf` (no calc/pseudopotentials), e.g. to
hand off to a real SIESTA NEB later. Verified live with a genuine single-atom vacancy
-migration hop in 7-atom bulk Si (both endpoints hand-built with matching atom
correspondence, not relying on `--autosort-tol`): a ~0.59 eV barrier (right order of
magnitude for the real, commonly-cited ~0.2-0.5 eV DFT literature value) and a reaction
energy of essentially 0 eV -- correct by symmetry, since both endpoints are the same
vacancy-in-bulk-Si situation at two symmetric-equivalent sites. The core pipeline
(endpoint pre-relax + interpolate/IDPP + climbing-image NEB) is itself extracted into
`run_single_neb`, once `stb-mldiffusion` needed to call it many times (once per
candidate defect hop) instead of duplicating it -- `main()`'s own single-pair report
just calls it once, folding what were previously 3 separate report sections into one
`[2] NEB PIPELINE (MACE)` section.

**`stb-mldiffusion`** (ML Simulations, `mldiffusion.py`) is the eighth tool in the ML
Simulations category -- vacancy migration-barrier screening. Where `stb-mlneb` needs
two already-built endpoint structures, this tool builds them itself: given a PERFECT
bulk reference structure and the index of one atom to remove (`--vacancy-index`), it
auto-detects every neighbor within the first coordination shell (`find_neighbor_shell`,
using each site's own local nearest-neighbor distance, not the structure's global
minimum -- matters for a non-uniform local environment), builds the vacancy-hop
endpoint pair for each candidate (`build_vacancy_hop_endpoints` -- same atom
-correspondence-preserving construction verified for `stb-mlneb`'s own test fixture:
every atom except the hopping one stays at the exact same list index/position), and
runs a real climbing-image NEB on every one of them via `stb-mlneb`'s own
`run_single_neb` (imported directly, not duplicated). Distinct from
`stb-mlmd`/`stb-aimdAnalysis`'s indirect, MSD-based diffusion coefficient: this gives
the actual migration barrier directly, useful in the rare-event/low-temperature regime
where an MD trajectory would need an infeasibly long run to see even one real hop.
Verified live on 8-atom bulk Si (vacancy at a site with exactly 4 tetrahedral nearest
neighbors -- the correct Si coordination number): all 4 auto-detected candidate hops
were found, and all 4 gave the exact same ~0.59 eV barrier -- expected by cubic
symmetry, since all 4 nearest-neighbor directions in diamond Si are crystallographically
equivalent, a strong physical-correctness signal for both the neighbor auto-detection
and the barrier computation.

**`stb-mlgcmc`** (ML Simulations, `mlgcmc.py`) is the ninth tool in the ML Simulations
category -- canonical/grand-canonical Monte Carlo adsorption for a single-atom
adsorbate species in a rigid host framework, driven by MACE. `--ensemble canonical`
keeps a fixed adsorbate count (displacement moves only); `--ensemble grand-canonical`
additionally inserts/deletes at a given chemical potential (`--mu`), following the
standard textbook (Frenkel & Smit) acceptance formulas `acc(insert) = min(1, V*z/(N+1)
* exp(-beta*dU))` / `acc(delete) = min(1, N/(V*z) * exp(-beta*dU))`, where `z =
exp(beta*mu)/Lambda^3` is the absolute activity and `Lambda` the adsorbate's own
thermal de Broglie wavelength (computed automatically from its mass and T). These
formulas were verified correct BEFORE ever touching MACE: a toy non-interacting
(`dU=0` always) ideal-gas Metropolis simulation, run for 2*10^5 steps, reproduced the
known-exact grand-canonical ideal-gas result `<N> = V*z` with `variance(N) ~= mean(N)`
(the Poisson-distribution signature) to within Monte Carlo noise.

**A real, important physics finding from live verification** (not a code bug -- the
formula was already proven exact above): `mu` is extremely sensitive. `beta =
1/(kB*T)` in eV^-1 is ~O(40) at room temperature, so the ideal-gas reference loading
`V*z` changes by roughly 50x per 0.1 eV change in `mu`. Verified live on a 105 Ang^3
graphene-monolayer cell with Ar at 300 K: `V*z` went from ~1e-4 (no adsorption) at
`mu=-0.5 eV` to ~2.6e4 (severe, unphysical overcrowding, 30-70+ Ar atoms crammed into
the cell) at `mu=0.0 eV` -- jumping straight past any sensible loading in between,
which is exactly why an initial coarse scan (`-1.0, -0.5, 0.0, 0.5, 1.0`) produced
wildly non-monotonic, overcrowded results. `run_gcmc` now prints this ideal-gas
reference loading up front (before the expensive interacting MC loop), flagging a
`[WARNING]` when it's far from O(0.1-10) -- catching a bad `--mu` choice cheaply
instead of after a long, wasted run. Once properly calibrated (solving for `mu` where
`V*z ~ O(0.1-1)`, here `mu` around -0.34 to -0.28 eV for this same cell), the isotherm
behaves physically sensibly: mean loading 0.0 (`mu=-1.0`, correctly far below the
transition) -> 5.0 (`mu=-0.28`), a clean, robust two-point increase used for
verification instead of a finer-grained scan (noisy at practical smoke-test step
counts, though still directionally correct: 0.49 -> 4.73 -> 4.30 across `-0.34, -0.30,
-0.26`, plausibly a real packing-limited plateau given the tiny graphene primitive
cell's lateral footprint, not just noise). Single-atom adsorbates only in this
version, and the host framework is always treated as rigid (never relaxed) -- both
documented limitations, exposed further only once there's a genuine need (same
"expose on first genuine use" policy as the rest of this suite).

**Real bug found via the interactive menu wrapper**: `run_mlgcmc_generator` never
prompted for `--equilibration-steps`, silently leaving it at its default (500) --
picking a modest `--steps` interactively (e.g. 150, well below 500) meant EVERY step
was discarded as equilibration, so the collected N/energy history arrays were always
empty (`<E>` printed as `nan`, no `mc_history.png`/`.dat` ever written). Fixed two ways:
`main()` now validates `--equilibration-steps < --steps` up front (`parser.error`
instead of a silent empty-statistics run), and the interactive wrapper now prompts for
`--equilibration-steps` explicitly instead of leaving it implicit.

**`stb-mladsorb`** (ML Simulations, `mladsorb.py`) is the tenth tool in the ML
Simulations category -- standalone adsorption-site screening. `stb-adsorb` already has
a `--ml-rank` mode (enumerate candidate ontop/bridge/hollow sites, relax each with
MACE-MP-0 substrate-fixed, rank by energy) used to prioritize which site(s) to send to
SIESTA, but it's bolted onto a DFT-oriented workflow: always needs `--calc`/
pseudopotentials even though the ML path never touches either, only reports a
RELATIVE energy (dE from the best site), and its `--ml-model` only accepts
`small`/`medium`/`large` (no `--custom-model`, the same gap `stb-mlneb` already found
and fixed for `stb-neb`). `stb-mladsorb` is the pure-MACE extraction of that
capability, additionally computing a genuine ML adsorption energy per site (`E_ads =
E_site - E_slab - E_isolated_adsorbate`, all three MACE energies -- the isolated
adsorbate reference is itself relaxed, matching real adsorption-energy convention)
rather than just a relative ranking number. Deliberately reuses (not duplicates)
`adsorb.py`'s own pure-geometry helpers (`resolve_slab_orientation`,
`parse_adsorbates`, `isolated_adsorbate_structure`, `molecule_extent`,
`write_site_plot`, `min_adsorbate_slab_distance`, `build_sweep_values`) -- finding
candidate sites/building the isolated-adsorbate reference is exactly the same problem
whether the result goes to SIESTA or straight into a MACE relaxation; only
`write_reference_folder`/`write_bsse_folders` (SIESTA-input-specific) aren't reused.
Verified live on H adsorption on a 2-atom graphene monolayer: the most stable site
(bridge, of the 4 candidates a 2-atom basis gives) came out at `E_ads ~ -0.97 eV`, the
right physical magnitude and sign for a real H-carbon chemisorption bond (MACE's
site-type preference can differ from DFT's, an expected fast-surrogate-potential
discrepancy, not a bug -- the point of this tool is prioritizing candidates for real
DFT, not replacing it).

A second pass (v1.1.0) added 4 more features, refactoring the single-model screening
loop into `screen_sites` (returns a sorted `scored` list) so both the main model and an
optional foundation-model comparison can share it without duplicating the relax loop:
- **Foundation-model comparison**: when `--custom-model` is given, also screens every
  site with the raw MACE-MP-0 foundation model (`--skip-foundation-comparison` to
  disable) and reports whether the same site ranks best for both -- same pattern as
  `stb-mlphonons`/`stb-mlelastic`. Verified live (wiring only, since no committed
  carbon+hydrogen SIESTA dataset exists to fine-tune against): a synthetic single-atom
  Si "slab" + Si adsorbate, using the same quick-fine-tuned Si model
  `test/5-mlsimulations/2-mlphonons`'s real SIESTA fixtures already produce for other
  tools' tests, correctly ran both models and compared their independent rankings.
- **`--n-orientations`**: for a multi-atom adsorbate, cheaply screens several random
  3D rotations per site (`pick_best_orientation`, single-point energy only, no relax)
  before the expensive full relax, keeping the best-looking one --
  `AdsorbateSiteFinder.add_adsorbate` always re-aligns the molecule's own z-axis to the
  surface normal regardless of input orientation (per its own docstring), so this
  samples the physically meaningful remaining degree of freedom (rotation about that
  normal, and any out-of-plane tilt). Verified live: H2O on graphene at a deliberately
  awkward height improved from `E_ads = 0.57 eV` (default G2 orientation) to `0.19 eV`
  with 8 sampled orientations -- a real, meaningful improvement, not noise.
- **Binding-energy-vs-height curve** (`plot_height_curves`, automatic whenever
  `--height-sweep` is used): plots E_ads vs. height per site instead of just flatly
  ranking every (site, height) pair together. Verified live that different starting
  heights for the SAME site can converge to genuinely different local minima after
  relaxation (a real, expected limitation of a local-relaxation screen, not a bug) --
  exactly the kind of thing this curve is meant to surface.
- **`--diffusion-barrier`**: estimates the adsorbate's own surface-diffusion barrier
  between its 2 most stable, distinct sites via a real climbing-image NEB, reusing
  `stb-mlneb`'s `run_single_neb` directly (`run_diffusion_barrier`) -- distinct from
  `stb-mldiffusion` (vacancy migration in the bulk), this is migration of an ADSORBATE
  across a surface. `freeze_substrate=True` is always used: the two site endpoints
  already share bit-identical substrate positions (both independently relaxed with the
  substrate held fixed by `screen_sites`), so `run_single_neb`'s own
  `compute_frozen_indices` auto-detects and freezes every substrate atom with zero
  extra bookkeeping. Verified live on H/graphene: the NEB's independently-computed
  reaction energy (0.2606 eV) matched the site-ranking energy difference between the
  same two sites (-0.7054 − (-0.9659) = 0.2605 eV) to within 0.0001 eV -- a rigorous
  cross-check between two completely separate code paths (the relax-based ranking and
  the NEB), not just a plausibility check.

**`stb-mleos`** (ML Simulations, `mleos.py`) is the eleventh tool in the ML Simulations
category -- a standalone equation-of-state (E vs V) calculation driven entirely by a
MACE potential: scans `--n-volumes` isotropically-scaled cell volumes around the
(pre-relaxed) equilibrium, relaxes atomic positions at each fixed scaled volume, and
fits an E(V) equation of state (`ase.eos.EquationOfState`, Birch-Murnaghan by default,
also Vinet/Murnaghan/stabilized-jellium) to get V0/E0/B0 (bulk modulus, GPa, via
`elastic_analysis.py`'s own `CONV_EVA3_TO_GPA`)/B0'. No phonons or force constants at
all -- unlike `stb-mlphonons --qha`, which needs the full volume-scan-plus-force
-constants machinery for its thermal-property-vs-volume curves, this tool only wants
the static E(V) curve itself, so it's much cheaper. Bulk (3D periodic) only, same
reasoning as `stb-mlmelting`/`stb-amorphize` (an equation of state assumes a
well-defined, bounded cell volume). Reuses the exact isotropic-scaling convention
`stb-mlphonons`'s `run_qha` already uses (`factor = (1 + strain_pct/100)**(1/3)`
applied to the cell vectors, so `--strain-range` means the same "percent volume
change" in both tools) rather than inventing a second one.

Deliberately distinct from `stb-mlelastic`'s own bulk modulus: that one comes from a
linear stress-strain fit around a single reference cell (the elastic STIFFNESS
TENSOR); this one comes from the curvature of the total-energy-vs-volume curve
itself. The two are independent methods and independent MACE evaluation loops on the
same structure -- comparing them is a genuine physical-consistency check on the MACE
potential, not a duplicate feature, and the report prints an explicit `[NOTE]`
pointing the user at this cross-check. Verified live on the same 8-atom bulk Si used
elsewhere in this session (`test/5-mlsimulations/3-mlelastic/si8.fdf`): `stb-mleos`
gives V0=163.19 Ang^3, E0=-42.967 eV (matching `stb-mlsearch`'s independently-found
energy minimum for the same structure, ~-42.966 eV), B0=71.22 GPa, B0'=4.32,
R^2=0.999993 (excellent fit) -- and `stb-mlelastic` on the exact same structure
reports Bulk Modulus (B) = 70.41 GPa, agreeing with `stb-mleos`'s B0 to ~1%, strong
evidence both independent methods are extracting the same real physics from the MACE
potential. `--custom-model` (+ `--skip-foundation-comparison` to opt out, on by
default) overlays a second E(V) curve and V0/E0/B0/B0' table fit with the raw MACE-MP-0
foundation model on the same volumes, same pattern as `stb-mlelastic`/`stb-mlphonons`.

**`stb-eosInputs`/`stb-eosAnalysis`** (Workflow item 18, `eos_inputs.py`/
`eos_analysis.py`) is the real-SIESTA analog of `stb-mleos` -- same isotropic-volume
-scan idea (`factor = (1 + strain_pct/100)**(1/3)` applied to the cell vectors, same
`--strain-range` meaning as `stb-mleos`/`stb-mlphonons`'s `--qha`), but generating
real DFT input folders instead of an in-memory MACE loop, mirroring the existing
`stb-strain`/`stb-strainAnalysis` prep+analysis pair. Stage 1 (`stb-eosInputs`) reads
a bulk (3D periodic) structure, rejects any vacuum-padded axis (same reasoning as
`stb-mleos`), and writes one `vol_<pct>/<input_basename>` folder per scanned volume
via `core/structure_io.py`'s `rewrite_fdf_lattice` -- only `%block LatticeVectors` is
replaced (scaled by the isotropic factor, `LatticeConstant` line untouched, same
`raw_lattice_vectors`/`rewrite_fdf_lattice` pattern `stb-strain` itself uses), so
positions (expected fractional) scale for free with the cell and everything else
(species, basis, pseudopotentials, SCF block) is preserved verbatim. The user runs
SIESTA in each folder with a FIXED cell (positions may relax; letting the cell relax
too would defeat the scan). Stage 2 (`stb-eosAnalysis`) aggregates every `vol_*`
folder's `calc.out` -- energy via `core.siesta_log.get_free_energy`, volume via
`core.siesta_log.get_outcell` (the actual cell SIESTA used for that run, read directly
from the `.out` file's last `outcell:` block, not re-derived from the folder's own
`.fdf` -- one fewer file dependency, same convention `stb-strainAnalysis` already
uses) -- and fits an equation of state. The fit itself (`ase.eos.EquationOfState`
wrapper, GPa conversion, R^2) is NOT duplicated from `stb-mleos`: both tools import it
from a new shared `core/eos_fit.py` (`fit_eos`/`normalize_eos_string`), extracted
specifically because `mleos.py` calls `core.deps.require_mace()` at import time --
merely importing it from a plain DFT-workflow tool that never touches MACE/torch
would have forced the heavy optional `ml` dependency to load just to fit a curve.
`stb-mleos` was refactored (v1.0.0 -> v1.0.1, no behavior change, confirmed via a
regression run reproducing its exact previous numbers) to use the same shared module.
Plots use this workflow category's own gnuplot `.dat`+`.gplot` convention (matching
`stb-strainAnalysis`/`stb-elasticAnalysis`/`stb-convergenceAnalysis`), not the newer
matplotlib convention the ML Simulations tools use -- `stb-eosAnalysis` is the first
new WORKFLOW_TOOLS analysis stage added this session, so it follows its own
category's established sibling convention instead. The report's closing `[NOTE]`
points at TWO independent cross-checks: `stb-elasticAnalysis`'s own stress-strain
-derived bulk modulus (same DFT setup, different method), and `stb-mleos` (the MACE
analog of this exact curvature-based method, for a DFT-vs-ML sanity check).

Verified end-to-end against a known synthetic ground truth rather than just "runs
without crashing": generated a real 7-point `vol_*` sweep from an 8-atom cubic-Si
fixture via `stb-eosInputs`, then injected synthetic `calc.out` files whose energies
follow an EXACT Birch-Murnaghan curve (`ase.eos.birchmurnaghan` itself, E0=-300 eV,
B0=90 GPa, B0'=4.3, V0=160.103 Ang^3) at each real scanned volume -- `stb-eosAnalysis`
recovered all four parameters exactly (R^2=1.000000), confirming the full pipeline
(folder scan -> `get_outcell` volume -> `get_free_energy` -> `fit_eos`) end-to-end.
Also verified the gnuplot `.dat`/`.gplot` pair actually renders (a real `gnuplot` run
produced a valid PDF from the two-block data file). `stb-eosInputs`'s own generation
was checked directly: `vol_0.00`'s raw lattice block is unchanged (scale factor
exactly 1 at zero strain), `vol_5.00`'s is scaled by exactly `(1.05)**(1/3)`, and a
vacuum-padded (2D) structure is correctly refused with the same message wording as
`stb-mleos`'s own bulk-only check.

A second pass added 3 more features to `stb-eosAnalysis` (v1.0.0 -> v1.1.0), all
operating on `core/eos_fit.py`'s shared fit -- no new SIESTA data-mining needed:
- **`--eos all`**: fits every supported form (birchmurnaghan/vinet/murnaghan/sj) on
  the exact same (volume, energy) data and reports V0/E0/B0/B0'/R^2 side by side in
  a `[3] EQUATION-OF-STATE FIT (ALL FORMS)` table -- a robustness check on the fit
  itself (do different functional forms agree?), not a second data-mining pass. The
  plotted curve still uses birchmurnaghan specifically. Mutually exclusive with
  `--target-pressure` (`parser.error` if both given -- inverting needs one specific
  form). Verified live against the same known-synthetic Birch-Murnaghan dataset used
  for the base tool's own verification: all 4 forms independently recovered B0 within
  ~0.2 GPa of the true 90 GPa.
- **`--target-pressure`**: inverts the fitted EOS (`P = -dE/dV`, via `core/
  eos_fit.py`'s new `invert_pressure` -- a numerical derivative, `np.gradient`, of
  the already-fitted curve, one implementation working for every EOS form rather
  than an analytic per-form derivative) to report the predicted equilibrium volume
  and linear lattice-parameter scale factor at one or more target external
  pressures (GPa) in a new `[3b] TARGET PRESSURE` section. Flags `[EXTRAPOLATED]`
  whenever the REQUESTED PRESSURE (not the returned volume) falls outside the range
  of pressures the fitted curve actually spans over the scanned volumes -- checking
  the returned volume instead was a real bug caught during design: `np.interp`
  clamps to the nearest volume at the scan boundary for an out-of-range pressure
  rather than extrapolating, so a wildly out-of-range request (e.g. 50 GPa against a
  scan only spanning about -4.6 to +6.4 GPa) would otherwise silently return an
  in-range boundary volume and look like a normal interpolation. Verified against
  the analytic Birch-Murnaghan pressure formula: inverting the same synthetic
  ground-truth fit at 5.0 GPa returned V=152.311 Ang^3, and evaluating the closed
  -form BM pressure at that volume gives back 4.9997 GPa; a 50 GPa request was
  correctly flagged `[EXTRAPOLATED]`.
- **v0-outside-scanned-range warning**: `core/eos_fit.py`'s `fit_eos` now calls
  `eos.fit(warn=False)` and checks `v0 outside [min(v), max(v)]` explicitly itself
  (`v0_outside_range` in the returned dict), instead of relying on `ase.eos`'s own
  `warnings.warn()` call for the same condition, which bypasses every caller's
  `print_dual`-formatted report entirely (easy to miss, wrong place/format if it
  shows up at all). Both `stb-eosAnalysis` and `stb-mleos` (bumped to v1.0.2) now
  print a clean in-report `[WARNING]` instead -- a real, if latent, gap in both
  tools that this pass fixed once noticed while implementing `--target-pressure`'s
  own more careful extrapolation check.

`stb-status` (Utils, `status.py`) prints a quick per-folder summary of a SIESTA calc:
run type (single-point/relaxation/AIMD, via `core.siesta_log.get_dynamics_type` +
`categorize_dynamics`), SCF convergence, max force, final energy, final ionic
temperature (AIMD only), and which key output files (`.XV`, `.bands`, `.DOS`, `.RHO`,
`.ANI`, ...) are present — "can I run stb-bands/-dos/-cube/-ani2traj on this folder
right now". `--path` accepts a glob (e.g. `strain_*`) to scan several folders at
once, printing a summary table at the end. Every field degrades to `None`/`unknown`
instead of raising so a batch scan always finishes and reports what it found per
folder, rather than aborting on the first folder missing a `.out`/`.fdf`.

`stb-archive` (Utils, `archive.py`) packages a finished calculation (inputs +
essential outputs, `DEFAULT_INCLUDE_EXTS`) into a single `.tar.gz` with an embedded
`MANIFEST.txt` (SystemLabel, run type, SCF convergence, final energy — the same
`core.siesta_log` fields `stb-status` reports) — the intentional complement of
`stb-clean`, which deletes that same class of files to prep for a restart instead of
preserving them for sharing/reproducibility. Default include set deliberately
excludes large regenerable intermediates (`.HSX`, `.DM`, `.RHO`, `.VT`, `.VH`,
extensionless logs like `MESSAGES`) — `--include`/`--exclude` adjust the set
per-run. `--path` glob support and per-directory failure handling mirror
`stb-status`. Shares `core.siesta_log.find_out_file` with `stb-status` (extracted
there once `archive.py` became a second consumer needing to locate a generically
-named log like `calc.out`, not just `<label>.out`) and
`core.structure_io.find_fdf_system_label` (same "extract on second use" pattern).

**`stb-xrd`** (Analysis, item 3.9, `xrd.py`/`core/xrd.py`, v2.1.0) simulates a powder
X-ray diffraction pattern from a structure using pyxtal's diffraction engine (not
pymatgen's — pyxtal is already a hard dependency for `stb-crystalcast`, and its
`Similarity` class has no pymatgen equivalent). Reports space group/crystal
system/lattice info and a diffraction-pattern summary (peak count, strongest peak,
resolution); the full peak-by-peak table (2theta, d, h, k, l, intensity) is never
printed to the console/report, only ever written to `xrd_pattern.dat`
(`--save-gnuplot`, with a complete header) alongside a stick-pattern `.gplot` script.
`--compare-to <path>` scores the simulated pattern against an experimental one via
`pyxtal.XRD.Similarity` (cosine-weighted, 0-1). `stb-xrdsearch`/`stb-xrdrank`
(Workflow item 6, XRD structure solution) reuse this same `core/xrd.py` pattern
-computation/experimental-reading machinery to drive a structure search against a
real experimental pattern, rather than just reporting one structure's own pattern.

**Two real, compounding bugs found and fixed in `--compare-to` (verified live)**:
comparing a structure's own simulated pattern against ITSELF — a natural
self-consistency sanity check — used to score only 0.32-0.46 similarity instead of
the expected ~1.0.
- **Bug 1 (silent column misread)**: `read_experimental_pattern` always read
  intensity from column index 1. That's correct for a plain 2-column file, but
  `stb-xrd`'s own `--save-gnuplot` output (`xrd_pattern.dat`) has 6 columns (`2theta
  d h k l intensity`), so column 1 there is `d` (the d-spacing in Ang) — silently
  fed into the similarity metric as if it were intensity, no error. Fixed by always
  reading intensity from the LAST column instead (identical behavior for a 2-column
  file, correct for `xrd_pattern.dat`'s 6 columns, and for any other convention that
  keeps intensity last).
- **Bug 2 (broadening mismatch, the deeper physics issue)**: even after fixing the
  column, the score was still only ~0.32 for the identical structure. The simulated
  side of the comparison (`xrd.get_profile()`) is always a Gaussian-broadened
  CONTINUOUS profile (FWHM=0.1 deg); a raw peak list (stb-xrd's own stick-pattern
  output, or any literature-reported indexed peak table) has no peak width at all.
  `pyxtal.XRD.Similarity`'s cubic interpolation draws a spurious curve through the
  empty gaps between sparse peaks, which does not match the broadened profile's true
  near-zero baseline there — comparing the two representations directly is
  apples-to-oranges regardless of which column is read. Fixed with a new
  `core/xrd.py::looks_like_peak_list` heuristic (coefficient of variation of the
  spacing between sorted 2-theta points: ~0 for a real, evenly-stepped continuous
  scan, ~1.7 for a genuine sparse peak list — verified on both a real 386-peak
  pattern and a synthetic uniform 0.02-deg-step scan) that auto-detects a peak-list
  -shaped `--compare-to` file and Gaussian-broadens it the same way
  (`core/xrd.py::broaden_peak_list`, same FWHM=0.1/res=0.01 defaults as the
  simulated side) before computing the similarity. `--raw-experimental` opts back
  out of this (compares exactly as read) for anyone who explicitly wants the old
  behavior. Verified live on Sn3O4: self-comparison via its own `xrd_pattern.dat`
  went from 0.4583 (both bugs) / 0.3204 (column fixed, still unbroadened) to
  1.0000 (both fixed); a genuinely different structure (NaCl/KBr) compared against
  the same pattern still scored well below 1.0, confirming the fix doesn't just
  collapse every comparison to 1.0; and a genuinely continuous, densely/evenly
  -sampled scan is correctly left un-broadened (still scores 1.0 against itself).

**`stb-adsorb`/`stb-adsorbBsse`/`stb-adsorbAnalysis`** (Workflow item 8, "Adsorption") was
restructured from a 2-stage prep+analysis pair into 3 stages after a real, verified physics
bug in the BSSE (Boys-Bernardi counterpoise) correction was found while reviewing a genuine
user SIESTA calculation: `stb-adsorb` (Stage 1) used to write the `bsse/site_*/
{bsse_slab,bsse_adsorbate}` ghost-fragment folders at PREP time, frozen at the same
pre-relaxation initial-guess geometry used for `sites/site_*/` itself -- before SIESTA had
ever actually run there. Read directly off that real calculation's `siesta.XV`: the O-Si
bond at one site relaxed from a 2.0 Ang initial guess down to ~1.62 Ang (plus ~0.46 Ang of
substrate puckering) -- BSSE grows quickly as atoms get closer, so evaluating the ghost
fragments at the far-too-separated initial geometry systematically UNDER-estimates the true
correction, undermining the "same geometry as the real, relaxed site" requirement
`examples/4.8-adsorption/README.md` (Section 2.3) already documented as the whole point of
a counterpoise correction. Per this suite's own established workflow-placement rule (an
"Analysis" stage must never generate new geometry/folders, only aggregate finished results
-- the same constraint `stb-hubbardu` already solved with its own 3-stage prep/alphas/
analysis split), the fix was a new dedicated Stage 2, `stb-adsorbBsse`: it reads each
`sites/site_*/`'s finished `siesta.XV` (skipping, with a clear per-site report, any site
that hasn't relaxed yet or whose atom count doesn't match its `structure.fdf` -- never a
hard failure for the whole batch) and writes the ghost-fragment folders
(`make_ghost_variant`/`write_bsse_folders`, moved here unchanged from `adsorb.py`) at THAT
relaxed geometry instead. Reuses each site's own already-copied pseudopotential files
directly as its ghost-species source (`pp_path=site_dir`), so this stage needs no
`-p`/`--pseudo-dir` flag of its own. `stb-adsorb` (Stage 1) no longer generates any `bsse/`
folder at all, and its `--bsse-correction`/`--no-bsse-correction` flags were removed
outright (not deprecated -- the old behavior was never physically correct, so keeping it as
a legacy opt-in would only perpetuate the bug). `stb-adsorbAnalysis` (renumbered Stage 2 ->
Stage 3) needed no aggregation-logic changes -- it already read `bsse/site_*/
{bsse_slab,bsse_adsorbate}/calc.out` exactly as `stb-adsorbBsse` now populates it; only its
explanatory report text was updated to point at the new stage. Verified via a synthetic
fixture (`test/4-workflow/8-adsorption/bsse/test.sh`): a fabricated `siesta.XV` (via sisl)
with the adsorbate deliberately moved 0.5 Ang from `stb-adsorb`'s own initial guess produces
ghost-fragment `structure.fdf` files whose written position matches that RELAXED geometry
exactly, not the original guess; a site with no `siesta.XV` yet is reported and skipped (not
fatal to the batch); a mismatched atom count (stale/hand-edited folder) is caught and
reported explicitly rather than silently misread.

**Second real, verified physics bug in `stb-adsorbBsse`, found from a real user SIESTA
calculation** (O adsorbed on a SiC slab): the BSSE-corrected `E_ads` came out MORE
negative (stronger binding) than the uncorrected value -- physically impossible for a
genuine Boys-Bernardi counterpoise correction, which can only ever weaken (never
strengthen) a binding energy, since adding ghost basis functions is a strictly
variational improvement. Root cause: `write_bsse_folders`'s `config_extra.fdf`
(`BSSE_SINGLE_POINT_BLOCK`) only ever forced `MD.VariableCell false` +
single-point-SCF -- it never carried over the site's own `Spin polarized`/
`Slab.DipoleCorrection` overrides (`read_original_calc_text` deliberately strips the
site's `%include config_extra.fdf` line to recover the user's raw `--calc` template, and
the BSSE stage's own replacement config never re-added them). So `bsse_slab`/
`bsse_adsorbate` silently fell back to the template's own `Spin non-polarized`, while
`E_site` (and the isolated-adsorbate reference) were genuinely spin-polarized -- for an
open-shell adsorbate like O, the resulting spin-restriction penalty on the ghost
fragment (~1.1 eV, verified against the real run's own `calc.out` files) was large
enough to swamp the genuine (~0.01-0.2 eV) BSSE effect and flip its sign. Fixed with a
new `read_site_theory_flags(site_dir)`: reads the site's own already-written
`config_extra.fdf` (rather than re-deriving the answer from CLI flags, so it's correct
regardless of which `stb-adsorb` version/flags produced that site) for whether `Spin`/
`Slab.DipoleCorrection` are present, and `write_bsse_folders` now composes its
`config_extra.fdf` as `BSSE_SINGLE_POINT_BLOCK` + (`core/adsorption_sites.py`'s own
`SPIN_POLARIZED_BLOCK`/`DIPOLE_CORRECTION_BLOCK` when the site actually used them) --
same shared constants `write_reference_folder` already uses, imported rather than
duplicated. `write_bsse_folders` now returns `(spin_polarized, dipole_corrected)` so
`main()`'s `[2] WRITING BSSE FOLDERS` section reports, per site, whether each was
inherited. Verified live against the real repro (`--no-force-spin` fixture in
`test/4-workflow/8-adsorption/bsse/test.sh`, 58 checks): a site built with
`--force-spin` (the default) propagates both `Spin polarized` and
`Slab.DipoleCorrection` into both ghost fragments; a site built with `--no-force-spin`
correctly gets neither `Spin` line (dipole correction is unconditional at the site
level, so it still propagates) -- confirming the fix reads the site's actual level of
theory rather than always forcing spin on.

**`stb-nebSites`** (Workflow item 9, "NEB / Reaction Path", `neb_sites.py`) is a new Stage 1
inserted ahead of the existing `stb-neb`/`stb-nebAnalysis` pair, which renumbered from
`4.9.1`/`4.9.2` to `4.9.2`/`4.9.3` -- the same insert-a-stage-and-renumber pattern already
used when `stb-adsorbBsse` became the new Stage 2 of item 8 above. Where `stb-neb` needs two
already-built, already-relaxed endpoint `.fdf` files, `stb-nebSites` builds the CANDIDATES for
that pair itself: given a slab and one adsorbate/molecule, it enumerates the symmetrically
distinct ontop/bridge/hollow sites (the exact same `AdsorbateSiteFinder` + `--symprec`
machinery `stb-adsorb` uses) and writes two endpoint folders, `site_A/`/`site_B/`, for
whichever two candidates you pick (`--site-a`/`--site-b`, 0-based indices into a single
combined table spanning every requested site type -- a hop can legitimately go ontop ->
hollow, not just between two sites of the same type). If either index is omitted, the tool
prints the candidate table (plus a "closest candidate pair" suggestion, the site-level
analog of `stb-mldiffusion`'s atom-level `find_neighbor_shell`) and exits, the same
ambiguous-selection convention `stb-hubbardu`'s `--atom-index` already uses -- never guesses.
Distinct from `stb-adsorb`: no isolated-adsorbate reference, no `E_ads` ranking, no
`--all-sites`/`--ml-rank`/`--height-sweep`/`--both-sides` -- this tool only ever cares about
two concrete points, not ranking many of them. Distinct from `stb-mldiffusion`: that one
builds a vacancy-hop pair in BULK via a nearest-neighbor-ATOM shell search and runs the NEB
itself with MACE; this tool builds an ADSORBATE-hop pair on a SLAB via symmetry-reduced site
enumeration and only writes real SIESTA input folders -- running the actual NEB stays
`stb-neb`/`stb-nebAnalysis`'s job.

Two refactors, both following this suite's own extract-on-(further)-use policy:
- **`core/adsorption_sites.py`** (new) -- `resolve_slab_orientation`/`reorient_vacuum_to_c`,
  `resolve_adsorbate`, `write_reference_folder`, `write_site_plot`, and
  `min_adsorbate_image_distance` moved out of `adsorb.py` once `stb-nebSites` became a second
  consumer; `adsorb.py` now imports them instead of defining them locally (`resolve_slab_
  orientation` was also changed to raise `ValueError`/return `(structure, relabel_note)`
  instead of printing-and-`sys.exit`-ing itself, so each of the two callers can report the
  error/info line in its own report's exact style). `parse_adsorbates`/`force_gamma_kgrid`/
  `force_spin_polarized`/`isolated_adsorbate_structure`/`min_adsorbate_slab_distance` stay in
  `adsorb.py`, specific to the isolated-`E_ads`-reference machinery `stb-nebSites` never
  touches (its own copy of the small, generic `min_adsorbate_slab_distance` distance-matrix
  check is a deliberate duplicate, not worth a cross-import -- same "plain generic helper"
  bar `stb-adsorb`'s own `build_sweep_values` duplicate already sets). Verified via full
  regression: `test/4-workflow/8-adsorption/{prep,bsse}/test.sh` both still pass unchanged
  (111 and 46 checks) after the extraction.
- **`core/structure_io.py::read_relaxed_or_input(path)`** (new) -- resolves a bare `.fdf` file
  (unchanged `read_fdf` behavior) or a directory (reads its `structure.fdf` as the base, and
  prefers a `*.XV` inside the same directory if one exists, same atom-count-mismatch guard as
  before) into `(FdfStructure, used_relaxed)`. The third consumer of "prefer a finished `.XV`
  over `structure.fdf`" (after `adsorb_bsse.py::read_relaxed_site_fdf` and
  `adsorb_analysis.py::read_site_geometry_atoms`) is what triggered the extraction:
  `stb-neb`'s own `--initial`/`--final` now accept a directory too (e.g. a `stb-nebSites`
  `site_A/`/`site_B/` folder) via this same function, printing a `[WARNING]` in `[0] RUN
  METADATA` (not fatal -- still builds the band from the unrelaxed guess, useful to preview
  the path) when no `.XV` is found yet. `adsorb_bsse.py::read_relaxed_site_fdf` was refactored
  into a thin wrapper around the shared function (its own "not relaxed yet" contract folds
  the shared `(structure, used_relaxed)` tuple into a plain `None`, since BSSE -- unlike
  `stb-neb` -- must skip an unrelaxed site outright, never fall back to its guess). Verified:
  `test/4-workflow/8-adsorption/bsse/test.sh` (46 checks) unchanged after the refactor; live,
  a fabricated `site_A/siesta.XV` (adsorbate shifted 0.3 Ang closer to the substrate, same
  sisl-write recipe as the BSSE test fixture) was confirmed to override `site_A/structure.fdf`
  in the generated `image_00/`, while an un-fabricated `site_B/` correctly still printed the
  `NOT YET RELAXED` warning and fell back to its own `structure.fdf` guess.

Verified live end-to-end on the same bare 2-atom graphene primitive cell `test/
4-workflow/8-adsorption/prep/structure.fdf` uses: `stb-nebSites --adsorbate H` (no
`--site-a`/`--site-b`) prints the same lateral-self-interaction `[WARNING]` Section 3.2 of
`examples/4.8-adsorption/README.md` documents for this exact fixture, plus a 4-row candidate
table (1 ontop + 3 bridge -- this basis exposes no `hollow` candidate to pymatgen's finder,
same observation already made for `stb-adsorb` on this fixture) and a "closest candidate
pair" suggestion; re-running with `--site-a 0 --site-b 3` writes `site_A/`
(`NumberofAtoms 3` -- 2 substrate + H) and `site_B/` sharing the exact same 2 substrate-atom
positions (trivial atom correspondence for the NEB step that follows, unlike
`stb-mldiffusion`'s vacancy case, which needs an explicit index-matching construction);
piping those two folders straight into `stb-neb --initial site_A --final site_B` produces a
5-image band with the `NOT YET RELAXED` warning on both, exactly as expected before any real
SIESTA relaxation has run.

A follow-up pass added three visualization outputs, matching the conventions `stb-adsorb`
already established rather than inventing new ones. All three share the same rotation trick:
`plot_slab(..., adsorption_sites=False)` draws the substrate, and every hand-placed marker's
raw (unrotated) Cartesian coordinate is first passed through
`pymatgen.analysis.adsorption.get_rot(pmg_structure)` -- the exact same rotation `plot_slab`
itself applies internally to its own adsorption-site markers (confirmed by reading its
source) -- before plotting, so every marker lands exactly where the substrate atoms are
drawn instead of being silently offset for any slab whose surface normal isn't already
parallel to Cartesian z:
- **`candidate_sites.png`** (`[2]`, `write_candidate_sites_plot`, whenever at least one
  candidate is found -- including the `< 2 candidates` error case, for context) -- every
  enumerated candidate marked with a small dot and its 0-based index as an adjacent text
  label (the exact numbering `--site-a`/`--site-b` and the printed table use), dot color
  encoding the site type (`_SITE_TYPE_COLORS`: ontop=red, bridge=blue, hollow=green) --
  deliberately small, uniform dots (not large/shaped markers) so the number label next to
  each one stays legible and the plot doesn't get visually noisy on a cell with many
  candidates. Deliberately NOT `core/adsorption_sites.py::write_site_plot` (still
  used by `stb-adsorb`'s own generic overview) -- that one draws pymatgen's OWN
  `find_adsorption_sites()` candidates at ITS OWN default settings, not necessarily this
  run's `--height`/`--symprec`/`--site-type`, and carries no index numbers at all; for a tool
  whose whole point is picking two specific numbered candidates, a plot showing a different,
  unnumbered site set would be actively misleading, not just unhelpful. Verified visually
  (rendered PNG inspected directly) against the same bare graphene fixture: index 0 (ontop,
  red) sits exactly on the substrate atom, indices 1-3 (bridge, blue) sit exactly at bond
  midpoints -- matching `[2]`'s printed table (`0=ontop, 1..3=bridge`) exactly.

  **Real bug found and fixed (reported live: dots/numbers invisible on a real, bigger
  structure, only ever showing up correctly on the tiny 2-atom test fixture)**: both plotting
  functions originally used a plain `zorder=100`/`101` for their own markers/labels, but
  `plot_slab` draws its own repeated (default `repeat=5`) substrate-atom circles at
  `zorder=2*idx`/`2*idx+1`, where `idx` runs over EVERY drawn copy -- up to `25*N` for an
  `N`-atom cell, so `zorder` alone can reach into the hundreds for anything much bigger than
  the 2-atom toy fixture, silently drawing our own markers UNDER the substrate instead of
  over it. Fixed with a shared `_FOREGROUND_ZORDER = 10_000` constant (comfortably above any
  realistic substrate-circle zorder) applied to every scatter/annotate call in both
  functions. Verified live on the 32-atom 4x4 graphene supercell (`stb-supercell -d 4 4 1`)
  that originally triggered the bug report: dots and index labels now render clearly on top.
- **`chosen_sites.png`** (`[4]`, `write_chosen_sites_plot`, once `--site-a`/`--site-b`
  resolve) -- marks ONLY the two chosen endpoints, labeled 'A'/'B'. Verified visually: site A
  (ontop) landed exactly on a substrate atom, site B (bridge) landed exactly at the midpoint
  between two neighboring substrate atoms -- both physically correct placements. Both plots'
  legends use `bbox_to_anchor` below the axes (not `loc='upper right'`) after a first attempt
  showed the legend clipping against the image edge for a site near a corner of the plotted
  cell.
- **`--view`** (needs `--site-a`/`--site-b` to have resolved) -- opens `clean_slab`,
  `site_A`, `site_B` as a 3-frame browser in ASE's interactive viewer
  (`core/ase_view.py::view_structure_interactive`, the same shared helper/graceful-
  no-display-failure convention `stb-adsorb`'s own `--view` already uses), a separate flag
  from `--view-plots` (matplotlib) for the same reason `stb-adsorb` keeps them separate.
  Verified headless (`MPLBACKEND=Agg DISPLAY=`, same convention as every other `--view` test
  in this suite): prints the 3 frame labels, then a graceful `[FAIL] Could not open the
  interactive 3D viewer` instead of hanging.

**`stb-neb` redesigned into 4 modes (`--mode 1-4`, default 3), plus a new CLI-only
`stb-nebCycle`** -- item `4.9`'s middle stage no longer has a single fixed path (plain
interpolation, optionally pre-shaped by a full MACE-MP-0 climbing-image NEB, always ending
in single-point SIESTA `image_NN/` folders). It now covers the whole cost/precision
spectrum between "100% MACE" and "100% real DFT" explicitly:

| Mode | Path engine | Output | SIESTA's role |
|---|---|---|---|
| 1 | 100% MACE-MP-0 (full climbing-image NEB) | one `neb_mace_result.json` | none at all |
| 2 | 100% MACE-MP-0, then one read | `image_NN/` (single-point) | one round, energy only |
| 3 (default) | MACE-MP-0 + a few real-DFT refinement cycles | `cycle_00/image_NN/` + printed loop | several rounds, real forces |
| 4 | 100% real-DFT NEB, no MACE at all | `cycle_00/image_NN/` + printed loop | every round, from cycle 0 |

`--mode` **replaces** the old `--ml-neb`/`--no-ml-neb` boolean outright (a deliberate
breaking CLI change, confirmed with the user rather than kept as a silently-coexisting
second way to ask for the same thing) -- modes 1/2/3 all still run through the exact same
`mace_relax.relax_neb` call the old `--ml-neb` did (`--ml-k`/`--ml-fmax`/`--ml-max-steps`/
`--ml-freeze-substrate`/`--ml-freeze-threshold` are unchanged, just now gated on `mode in
(1, 2, 3)` instead of a single flag); `--idpp`/`--ml-prerelax-endpoints` stay orthogonal to
`--mode`, available in all four. Modes 3/4 write to `cycle_00/` (not the output root
directly) and print a `[4] CLUSTER SUBMISSION` section with a ready-to-paste `sub.sh` loop
-- `--climb-after` defaults suggested per mode (`0` for mode 3, since MACE already
converged a well-shaped climbing-image band before handoff; `5` for mode 4, since climbing
too early on a raw interpolated path risks locking onto the wrong image as the saddle, the
same reasoning `relax_neb`'s own two-stage `climb=False`-then-`True` approach exists for).

**`stb-nebCycle`** (`neb_cycle.py`) is the new piece modes 3/4 depend on: a genuinely
CLI-only tool, deliberately **not** registered in `stb_suite.py`'s `WORKFLOW_TOOLS` (no
interactive-menu entry at all, per explicit request) since its whole reason to exist is
being called from *outside* any live Python process -- inside a cluster submission
script's loop, alternating with real (possibly hours-later, queued) SIESTA runs. One call
= one real-DFT NEB step:

1. Finds the highest-numbered `cycle_NN/` under `--dir`, reads every `image_*`'s
   `structure.fdf` (positions) + `calc.out` (energy, `core.siesta_log.get_free_energy`) +
   `<SystemLabel>.FA` (forces, `core.siesta_log.read_fa_forces` -- see below).
2. Builds one `ase.Atoms` per image, each wrapped in an
   `ase.calculators.singlepoint.SinglePointCalculator(atoms, energy=E, forces=F)` --
   injecting an already-finished result as if it came from a live calculator.
3. Builds `ase.mep.neb.NEB(images, k=..., climb=<cycle >= --climb-after>,
   method="improvedtangent")` -- same construction `core/mace_relax.py::relax_neb` already
   uses, just with a pre-computed-result calculator instead of a live MACE one; endpoints
   are never touched, same guarantee `relax_neb` documents.
4. Takes exactly one **`FIRE(neb, restart=<state-file>).step()`** -- deliberately `step()`
   directly, NOT `opt.run(fmax=..., steps=1)`. Verified live this distinction is load
   -bearing: `ase.optimize.optimize.Dynamics.irun()` (what `run()` drives) re-evaluates
   forces on the JUST-MOVED geometry right after `step()` to decide its own convergence
   status for the next iteration -- but a `SinglePointCalculator` cannot recompute
   anything for a geometry it wasn't constructed with, and this crashed with
   `ase.calculators.calculator.PropertyNotImplementedError` the moment `run()` was tried.
   Calling `step()` directly only ever evaluates forces on the CURRENT, already-known
   geometry (needed for the gradient) and never re-queries afterward. Convergence is
   instead checked explicitly beforehand, against the same pre-step forces already
   computed for the report.
5. Below `--fmax`: writes a `NEB_CONVERGED` sentinel at the root and stops (no new cycle;
   re-running again is a clean no-op once this file exists). Otherwise: writes
   `cycle_{N+1}/image_NN/` with the updated positions (reusing `neb.py::write_image_folder`
   directly -- import, not duplicate -- and each new image's own already-copied
   pseudopotentials as the `pp_path` source for the next cycle, same convention
   `adsorb_bsse.py` already uses for its own ghost-fragment folders).

**Restart mechanism, verified directly against the installed ASE source before writing
any of this** (not assumed): `Optimizer.__init__(..., restart=<path>)` calls `self.read()`
if the file exists, else `self.initialize()`; `FIRE.step()` ends every single call with
`self.dump((self.vel, self.dt))` (JSON via `ase.io.jsonio.write_json`, not pickle) -- so
even one `step()` correctly persists FIRE's velocity/timestep state for the next, wholly
separate process invocation to pick back up. Proved live, not just by reading the source:
two sequential `stb-nebCycle` calls with the exact SAME fabricated forces reapplied each
time gave a step-2 displacement of a tracked interior image that was **exactly 2.000x**
step-1's displacement (`FIRE`'s momentum genuinely accelerating while consecutive
gradients agree) -- and a wholly independent, freshly-started run from the same starting
geometry (no restart file) reproduced step-1's exact displacement again, byte-for-byte, to
`< 1e-9` Ang, confirming the acceleration seen in the first run really did come from the
restart file, not from anything else differing between the two runs.

**A real, live-caught bug, now guarded against explicitly**: `ase.mep.neb.NEB`'s
`"improvedtangent"` tangent estimate normalizes by its own norm
(`tangent /= np.linalg.norm(tangent)`) -- a 0/0 whenever two neighboring images become
(near-)degenerate (identical position and/or energy), silently producing NaN forces
instead of raising. Caught live while stress-testing the restart mechanism with repeated,
disconnected-from-geometry random forces: NaN quietly entered `FIRE`'s own persisted
velocity state, then surfaced as a literal `nan   nan   nan` written into the NEXT cycle's
`structure.fdf` -- a corrupted geometry with no error anywhere until something downstream
inevitably broke. `stb-nebCycle` now evaluates `neb.get_forces()` inside
`np.errstate(invalid="raise")` and refuses to take a step (clear `[ERROR]`, exit 1, no
cycle folder written) the moment a non-finite force appears, rather than silently
propagating it. Verified deterministically (not relying on reproducing the original
organic trigger) with a purpose-built degenerate fixture: `stb-neb -i x.fdf -f x.fdf`
(identical initial/final) puts every image at the exact same position, and any nonzero
fabricated force reliably reproduces the exact same crash-turned-clean-error.

**`core/siesta_log.py::read_fa_forces`** -- promoted from `mlff_analysis.py` (previously
its own only consumer) once `neb_cycle.py` became a second, following this suite's usual
extract-on-second-use policy; `mlff_analysis.py` now calls `siesta_log.read_fa_forces`
directly instead of keeping a local copy, no behavior change.

**`stb-nebAnalysis` (`neb_analysis.py`) becomes a mode-aware router**, reading a new
`# MODE: N` marker `stb-neb` now writes into `neb_setup.txt` (alongside the pre-existing
`# ML_NEB_USED:`/`# IMAGE_TABLE`) to decide how to read a study, with each mode's own
dedicated code path:
- **Mode 1**: `run_mode1_analysis()` reads `neb_mace_result.json` directly -- no
  `calc.out` anywhere, since SCF/residual-force diagnostics simply don't apply to a
  pure-MACE result. `read_mode1_json()` converts the JSON's `images` list into the exact
  same `ImageRow` shape the rest of the module already works with, so
  `fit_spline_barrier()`/`write_curve_plot()` are reused unmodified as a genuine
  cross-check against the JSON's own already-reported barrier (both independently derived
  from the same underlying MACE energies, one from ASE's `NEBTools.get_barrier(fit=True)`
  Hermite-style fit at prep time, one from this module's own cubic-spline fit at analysis
  time). `--apply` is refused with a clear, non-fatal note (mode 1 never writes a
  `structure.fdf` per image at all -- JSON only).
- **Mode 2**: unchanged -- exactly today's pre-existing `image_NN/` code path, no new
  code touches it at all.
- **Modes 3/4**: `find_analysis_cycle()` (thin wrapper reusing `neb_cycle.py`'s own
  `find_latest_cycle`/`CONVERGED_SENTINEL` -- import, not duplicate) picks the
  highest-numbered `cycle_NN/` under `--dir`, and every downstream path that used to read
  `os.path.join(args.dir, label)` now reads `os.path.join(images_root, label)` instead
  (`images_root` defaulting to `args.dir` for modes 2 and the pre-existing fallback path,
  becoming the discovered `cycle_dir` only for 3/4) -- the reaction-coordinate table
  itself still comes from the ORIGINAL `neb_setup.txt` at the root (written once, at
  `cycle_00` prep time; positions shift slightly cycle to cycle but the path's own
  topology/ordering doesn't), matching `adsorb_analysis.py::read_site_table`'s own
  documented convention of never trusting a table's own recorded directory column. A
  `[0]`-section `[NOTE]` reports plainly whether the analyzed cycle is the genuinely
  `NEB_CONVERGED` one or just "the latest refinement so far". A bonus, non-blocking
  `[2b] BARRIER VS. CYCLE` section (`collect_cycle_barrier_history`/
  `write_cycle_convergence_plot`) scans every COMPLETE `cycle_*` (every image's energy
  readable; an in-progress/partial cycle is silently skipped, not an error) and plots the
  forward barrier's trend across the real-DFT refinement history once >= 2 complete
  cycles exist -- a genuinely new capability with no analogue in the old single-path tool,
  made possible only because modes 3/4 keep every cycle's folder instead of overwriting
  in place.

Verified end-to-end for all 4 modes against the same 9-atom graphene+H toy fixture
`test/4-workflow/9-neb/prep/{initial,final}.fdf` already uses: mode 1 produced a valid,
schema-complete JSON with a real (if tiny, toy-fixture-typical) barrier; mode 2 reproduced
the exact old `--ml-neb` folder layout and report wording; mode 3/4 both wrote
`cycle_00/image_NN/` with the correct climb-after suggestion in the printed loop; and
`stb-nebAnalysis` correctly routed all four (including a live two-cycle mode-3 scenario
where a second, refined `cycle_01/` was correctly preferred over `cycle_00/`, the
convergence-history plot correctly picked up both, and the `[NOTE]` wording flipped
correctly once a `NEB_CONVERGED` sentinel was added).

**Two real bugs found and fixed from live user reports on the `stb-nebSites` -> `stb-neb`
handoff, both surfacing as a scrambled/wrong reaction path with no crash**:

1. **Distance-based atom matching (`--autosort-tol`) had no working fallback for a small/
   densely-packed cell.** A real user's `site_A`/`site_B` pair (independently SIESTA-relaxed)
   raised pymatgen's `Unable to reliably match structures with autosort_tol=X` -- expected,
   since real substrate relaxation near an adsorption site routinely moves atoms further than
   any single distance tolerance can unambiguously resolve. `neb.py`'s `main()` now catches this
   `ValueError` and automatically retries once with `autosort_tol=0` (index-based: trust
   `--initial`/`--final`'s own on-disk atom order directly, no distance matching), printing a
   `[WARNING]` instead of a hard `[ERROR]` -- this is the right correspondence whenever both
   endpoints share a guaranteed matching atom order (any `stb-nebSites` pair). Verified live
   that raising `--autosort-tol` instead does NOT reliably fix this for a small cell (tested
   0.3/1.0/2.0 Ang against the real repro, all failing differently, 2.0 Ang making every atom
   ambiguous) -- `autosort_tol=0` was the only fix that actually worked.

2. **`--autosort-tol 0` was itself never PROOF of a correct correspondence, only an assumption**
   -- reported live as the path still coming out "completely scrambled" even after fix #1.
   Traced to a real, confirmed gap: `core/structure_io.py::read_relaxed_or_input` (the function
   `stb-neb` uses to prefer a folder's finished `siesta.XV` over its `structure.fdf` guess) took
   `.XV`'s positions **purely by list index** and re-tagged them with `structure.fdf`'s symbols
   -- the `.XV`'s own per-atom species column was never read or cross-checked, only an atom
   -*count* check existed. Confirmed (by reading pymatgen's own source) that
   `Structure.interpolate(..., autosort_tol=0)` only raises on a species mismatch AT THE SAME
   INDEX -- a same-species atom swap (the common case for any multi-atom-of-one-element
   substrate) passes completely silently, exactly matching the reported symptom (no crash,
   wrong physics). Fixed two ways, mirroring the user's own suggested design (an explicit
   per-atom identity record, not a distance-based guess):
   - `read_relaxed_or_input` now ALSO cross-checks the `.XV`'s own per-atom species (via sisl,
     same `[a.symbol for a in geom.atoms]` idiom already used by `adsorb_analysis.py`/
     `mlff_analysis.py`) against `structure.fdf`'s expected per-index species, unconditionally
     (no new parameter) -- strictly additive, no plausible false positive for a correctly
     -formed folder, and benefits every existing caller (`stb-neb`, `adsorb_bsse.py::
     read_relaxed_site_fdf`) for free.
   - **New `core/neb_manifest.py`** (mirrors `core/dftu_data.py`'s `run_manifest.json`/
     `load_manifest` pattern used by `stb-hubbardu`/`stb-hubbardUAlphas`): `stb-nebSites` now
     writes an identical `neb_manifest.json` into BOTH `site_A/` and `site_B/` (not a shared
     root file -- `stb-neb` only ever receives one folder at a time via `--initial`/`--final`),
     recording the folder's actual on-disk `species_sequence` (re-read back from the
     just-written `structure.fdf`, since `write_fdf` regroups atoms by species before writing
     -- the physical file order is not simply "substrate then adsorbate"), `pair_id` (a
     `uuid.uuid4().hex[:12]` shared by both folders of one `stb-nebSites` run), and atom-count
     breakdown. `stb-neb`'s new `resolve_manifest_pair()` loads+cross-validates both manifests
     (when BOTH `--initial`/`--final` are directories carrying one) two ways: against each
     OTHER (`neb_manifest.validate_manifest_pair`, `species_sequence` list-equality is the
     actual proof of a correct index correspondence -- `pair_id` mismatch alone is only an
     informational `[WARNING]`, not fatal) and against the structure ACTUALLY read back by
     `read_relaxed_or_input` (catches a stale manifest left over after a `structure.fdf`
     regeneration/hand-edit). When both checks pass, `stb-neb` forces `--autosort-tol 0`
     automatically (now PROVEN, not guessed) and prints a clean `[OK] Atom order PROVEN via
     neb_manifest.json` confirmation instead of the fallback `[WARNING]` from fix #1. When
     EITHER side lacks a manifest (a hand-built pair, or an older `stb-nebSites` run predating
     this feature): falls back to exactly the pre-existing distance-based behavior, zero
     regression. Any manifest inconsistency (tampered `species_sequence`, a required field
     missing, a stale/regenerated folder) is a clean `[ERROR]`, never a raw traceback.

   Also added, motivated by a SEPARATE real finding from the same live investigation: running
   the user's actual (manifest-proven-correct) `site_A`/`site_B` pair through `stb-neb --mode 1`
   gave a physically implausible ~47 eV barrier (16/33 atoms moved > 0.3 Ang between endpoints,
   no `--idpp`) -- correct atom identity alone does not prevent plain linear interpolation from
   producing atomic clashes in intermediate images across a large rearrangement. New
   `check_endpoint_displacement()` reports max/mean per-atom displacement between
   `--initial`/`--final` and, only when a broad fraction of the structure moved (>= 3 atoms AND
   > 10% of the structure -- deliberately gated so a normal single-atom reaction hop, e.g. this
   suite's own `initial.fdf`/`final.fdf` fixture at 1/9 atoms, never triggers it) and `--idpp`
   wasn't used, prints a `[WARNING]` recommending `--idpp`. Advisory-only, same convention as
   `check_path_quality` -- never auto-enables `--idpp` itself, since that would silently change
   the actual interpolated path.

   Verified end-to-end with a real `stb-nebSites`-generated `site_A`/`site_B` pair (synthetic
   `siesta.XV` fabricated via sisl, same recipe as `test/4-workflow/8-adsorption/bsse/test.sh`):
   manifest-proven pair skips distance matching entirely and never prints the fallback text;
   removing the manifest from either or both sides reproduces today's exact behavior; a
   manifest edited to disagree with its own folder, one edited to disagree with its pair
   partner (while staying internally self-consistent -- requires editing the raw `.fdf` TEXT
   directly rather than round-tripping through `read_fdf`/`write_fdf`, since `write_fdf`'s own
   species-regrouping silently undoes an in-memory atoms-list reorder between atoms of
   different species), and one missing a required field each produce a distinct, clean
   `[ERROR]`.

**Real bug found and fixed (user-reported): `core/structure_io.py::write_fdf` could write a
non-sequential `%block ChemicalSpeciesLabel`**, which SIESTA requires to be gap-free
(1, 2, 3, ...). `write_fdf` used to write each surviving species' id verbatim from
`structure.species_meta[symbol]['id']` -- but `species_with_atoms` (its own zero-atom-count
filter, e.g. a species declared in the file but with no atoms left after some upstream
transformation) drops species from the OUTPUT without ever renumbering the ids of the ones
that remain, so a `species_meta` like `{H:1, He:2, O:3, C:4}` with no He atoms left produced
`1, 3, 4` on disk -- id `2` never reclaimed. `from_pymatgen` had the same root cause one level
up: reusing a caller-supplied `species_meta` (the universal `species_dict(...)` pattern used
by `stb-defect`/`stb-passivate`/`stb-nebSites`/every structure-transforming tool) after a
species was removed elsewhere carried the same gap forward even before `write_fdf` ever saw
it. Fixed at the one shared chokepoint every `.fdf` write in the suite funnels through:
`write_fdf` now always renumbers `species_with_atoms` fresh to 1..N at write time (ignoring
whatever `species_meta['id']` said), and `from_pymatgen` renumbers its returned
`species_meta` the same way before returning (defense-in-depth, keeps the in-memory
`FdfStructure` consistent with what gets written). A codebase-wide audit (two parallel
research passes) confirmed every OTHER `.fdf`-species-writing path -- `translate.py`
(delegates to `write_fdf`), `adsorb_bsse.py`'s ghost species, `core/passivation.py`/
`passivate.py`/`slab.py --passivate`, `defect.py`, `crystalbuilder.py`, `crystalcast.py`
(all via `ensure_species_id`/`from_pymatgen`/`write_fdf`), and `cohesive_energy.py`'s two
hand-rolled writers (structurally incapable of a gap: one always single-species id=1, the
other a from-empty incrementing counter with a dedup guard) -- either already produced
correct numbering or inherits the fix automatically through the shared functions, so no
other file needed a change. Verified directly (a deliberately gapped `species_meta` fixture
now writes `1, 2, 3`, not `1, 3, 4`) and end-to-end through a real CLI call
(`stb-defect --type vacancy` removing the one atom of a 4-species structure's rarest
species correctly renumbers the remaining 3 species to `1, 2, 3`).

**`stb-nebSites` (`neb_sites.py`) gained `center_slab_in_vacuum()`, a real bug fix reported
live**: a freshly-cut slab conventionally starts near frac `z=0` with all vacuum stacked
above it (pymatgen's own `SlabGenerator` convention, `stb-slab`'s output). During a real
SIESTA relaxation, an atom that starts at frac `z` close to 0 can drift slightly negative
and get wrapped by PBC to frac `z` close to 1 (the TOP of the cell) instead -- which then
reads as an enormous, spurious per-atom displacement to anything comparing this relaxed
site against another snapshot of the same atom (exactly `stb-neb`'s own endpoint-matching
step, see the manifest/`--autosort-tol` fixes above -- this is a THIRD, independent
contributor to the same "structure completamente bagunçada" class of symptom, this time a
genuine coordinate artifact rather than an atom-identity or interpolation-quality issue).
`center_slab_in_vacuum()` runs right after `resolve_slab_orientation` (so vacuum is
already relabeled to c) and before any site is enumerated: computes the slab's frac-z span
`[min, max]`, and if its midpoint is more than 1 Ang (in real Angstrom, via the c-lattice
-vector length) away from the cell's own midpoint (`0.5`), rigidly translates every atom
along c (`Structure.translate_sites(..., to_unit_cell=True)`, wrapping into the cell) so
the slab sits centered in the vacuum gap instead -- a pure coordinate-origin shift, no
change to any bond length/angle. Silent no-op (structure returned unchanged, not even
copied) when already within that 1 Ang tolerance, so a well-behaved input isn't
gratuitously rewritten. Reported in `[1] REFERENCE SLAB` either way (`[INFO] Slab was near
a cell boundary...` with the exact shift applied, or `Slab already centered...`). Verified
live: a synthetic near-boundary fixture (frac `z=0.02`, 20 Ang cell) correctly shifted by
`+9.6 Ang` to frac `z=0.5`, with the adsorbate still placed the requested height above the
now-recentered surface (frac `z=0.60` for `--height 2.0`); the suite's own already-centered
graphene fixture (frac `z=0.5` from the start) correctly triggered no shift at all.

**`stb-adsorb` (item 4.8.1) gained the same `center_slab_in_vacuum()` recentering,
user-requested** -- `neb_sites.py`'s own definition moved into `core/adsorption_sites.py`
(extract-on-second-use, `neb_sites.py` now imports it instead of defining it locally, no
behavior change there) once `adsorb.py` became a second consumer. Called right after
`pmg_structure = structure_io.to_pymatgen(structure)` (so every downstream step --
site-finding, `clean_slab/`, every `sites/site_*/` -- already sees the recentered
geometry), with the `[INFO]`/`Slab already centered` report line moved into `[1]
REFERENCE FOLDERS` (via `print_dual`, so it's captured by `--save-report`) rather than
`adsorb.py`'s own pre-existing `relabel_note` placement (a plain `print()` before the
report file even opens, left exactly as-is -- not a bug this change was asked to fix,
so deliberately not touched). Verified live with the identical repro `stb-nebSites`
itself was verified against: a synthetic frac `z=0.02` structure in a 20 Ang cell
shifted by exactly `+9.600 Ang` to frac `z=0.5`, confirmed in the written
`clean_slab/structure.fdf`; the suite's own already-centered fixture correctly
triggers no shift.

**`stb-adsorb` (4.8.1): DFT-D3 van der Waals correction forced unconditionally on every
folder it writes, and systematic multi-orientation MACE screening for molecular
adsorbates, both user-requested.**

**Van der Waals**: new `core/adsorption_sites.py::VDW_CORRECTION_BLOCK` (`DFTD3
.true.`, same `config_extra.fdf`-include mechanism as `FIXED_CELL_BLOCK`/
`SPIN_POLARIZED_BLOCK`/`DIPOLE_CORRECTION_BLOCK`) and a new `force_vdw=False` parameter
on `write_reference_folder`. `stb-adsorb`'s 3 call sites (`clean_slab/`, the isolated
-adsorbate reference, `sites/site_*/`) all pass `force_vdw=True` unconditionally --
same "hardcoded `True`, no opt-out flag" style `force_dipole=True` already used, not
`force_spin`'s opt-out-able style, since (unlike the dipole correction) dispersion is
physically meaningful for every folder alike, boxed molecule included. **Real-evidence
syntax**: confirmed against an actual, successfully-completed SIESTA run in this exact
workflow (`MESSAGES: "DFT-D3: loading default parameters for functional PBE"`) --
SIESTA's OTHER vdW mechanism (a native `XC.Functional VDW` non-local functional) has no
evidence of use anywhere in this repository and is NOT what this forces. Propagated
everywhere the same "same level of theory as `E_site`" rule already applies to Spin/
dipole (the exact bug class fixed earlier this session for BSSE): `adsorb_bsse.py`'s
and `adsorb_analysis.py`'s own `read_site_theory_flags` now each return a 3-tuple
`(spin, dipole, vdw)` (checking for `"DFTD3"` in the site's `config_extra.fdf`
alongside `"Spin"`/`"Slab.DipoleCorrection"`), propagated into the BSSE ghost fragments
and the Gibbs-prep Hessian displacement folders respectively, same as spin/dipole
already were. Deliberately NOT extended to `stb-nebSites`/`neb_sites.py` (same
`write_reference_folder`, but out of scope this round, same scoping precedent already
set earlier this session for `wrap_into_cell`).

**A second real, latent bug caught while wiring vdW into the Gibbs-prep isolated
-adsorbate reference** (`adsorb_analysis.py`'s `[6] GIBBS FREE ENERGY (DG) PREP`): the
isolated reference's own Hessian displacement folders were built from `open(ads_dir/
calc.fdf).read()` directly -- unlike the site side, which already used
`read_original_calc_text` to STRIP the leading `%include config_extra.fdf` line before
reuse. Since `write_gibbs_folder` never copies a `config_extra.fdf` into `disp_NNN/`
(only `calc.fdf`+`structure.fdf`), every isolated-reference Gibbs folder was carrying a
literal, DANGLING `%include config_extra.fdf` that SIESTA would fail to parse --
present since the ΔG feature first shipped, only surfaced now because adding vdW
required re-deriving what flags to re-apply and forced a closer look at this exact
code path. Fixed by switching to `read_original_calc_text` there too (matching the
site side exactly), then explicitly re-applying whatever `read_site_theory_flags(
ads_dir)` finds (`Spin`, `DFTD3` -- never `Slab.DipoleCorrection`, a boxed molecule
isn't a slab, matching `write_reference_folder`'s own `adsorbate_dir` call having no
`force_dipole`) via new `_force_vdw_correction_text` (mirroring the existing
`_force_spin_polarized_text`/`_force_dipole_correction_text` flat-calc.fdf helpers).
Verified live: `gibbs/<adsorbate>_isolated/disp_*/calc.fdf` no longer contains
`%include config_extra.fdf` at all, and correctly carries `Spin polarized`+`DFTD3
.true.` (no `Slab.DipoleCorrection`) inherited from the isolated reference's own
folder.

**Orientation sampling**: new `core/adsorption_sites.py::generate_systematic_
orientations(mol, n_polar, n_azimuthal)` and `deduplicate_orientations(entries,
n_substrate, rmsd_tol, energy_tol)`, wired into `stb-adsorb --ml-rank` via 4 new flags
(`--n-orientations-polar`/`--n-orientations-azimuthal`, both default `1` = today's
exact single-orientation behavior; `--orientation-top-k`; `--orientation-rmsd-tol`,
default 0.3 Ang) -- deliberately NOT reusing `mladsorb.py::pick_best_orientation`
(fully random `Rotation.random`, single-point energy only, picks ONE winner before a
single relax): the user explicitly asked for a SYSTEMATIC (reproducible, no seed) set,
with a FULL relax of every sampled orientation (not just single-point), keeping the
Top-N unique survivors -- different enough from the existing random/single-point
-then-relax-once pattern that reusing it wasn't practical.
- `generate_systematic_orientations`: `n_polar` reference directions (in the
  molecule's OWN unrotated frame) spread evenly over a sphere via a Fibonacci-sphere
  lattice (deterministic, no pole-clustering the way a naive lat/long grid has, no RNG
  needed) -- exploits that `AdsorbateSiteFinder.add_adsorbate` always translates
  whichever atom(s) sit at the molecule's OWN most-negative-z to touch the surface, so
  pre-rotating a different reference direction to that "most negative z" position is
  exactly how to control which part of the molecule (for H2O: O, either H, or
  in-between) ends up facing the surface. Each direction combined with `n_azimuthal`
  evenly-spaced spins (0-360 deg, exclusive) about that same new z-axis, via
  `scipy.spatial.transform.Rotation.align_vectors`+`from_rotvec` (already a suite
  dependency, same `scipy.spatial.transform` `mladsorb.py` already imports). Verified:
  a rigid rotation (bond lengths and center-of-mass exactly preserved to `1e-6`) with
  `n_polar=6, n_azimuthal=4` -> 24 pairwise-distinct orientations, none identical.
- `deduplicate_orientations`: greedy, RMSD on adsorbate atoms only (no Kabsch alignment
  needed -- same substrate, same reference frame across every candidate at one site)
  computed directly (no periodic minimum-image correction -- documented simplification,
  not expected to matter for independently-relaxed candidates starting from the same
  site coordinate). Drops a candidate only if it's within BOTH `rmsd_tol` AND
  `energy_tol` (default 0.01 eV) of an already-kept one -- requiring both avoids
  merging two genuinely different minima that coincidentally share one metric.
  Verified: near-geometry+near-energy pairs collapse to 1; near-geometry-but-far-energy
  and far-geometry-but-near-energy pairs both correctly stay as 2.

Integration into `adsorb.py`'s existing `[3] ML PRE-SCREENING` `--ml-rank` loop: for
each `(site, height)` candidate, generates the orientation set, runs a full MACE relax
(substrate `FixAtoms`, same `mace_relax.relax` call already used) per orientation,
sorts by energy, deduplicates, optionally trims to `--orientation-top-k`, and emits one
`scored` entry per surviving orientation (tuple grew from 6 to 7 elements, `orient`
appended last: `0` when sampling wasn't requested -- preserving every existing table/
plot/folder-naming output byte-for-byte -- else the 1-based post-dedup energy rank).
Flows into the SAME `[4] WRITING SITE FOLDERS` step unmodified, just with more
`site_..._orientN` labels; the existing `--top-k` (trims across ALL candidates,
existing flag) composes cleanly on top, unchanged. `plot_ml_ranking`'s bar labels grow
a `-oN` suffix when orientations were sampled. Verified live end-to-end on the suite's
own graphene fixture with a real H2O adsorbate (`--n-orientations-polar 4
--n-orientations-azimuthal 2 --orientation-top-k 3`): 8 orientations sampled, 3 unique
survivors, 3 real `site_1_ontop_orientN/` SIESTA folders written, each with genuinely
different adsorbate-atom z-positions (confirmed by direct inspection) and `DFTD3
.true.` correctly present in each.

**A follow-up pass on orientation sampling added 4 more usability features, all
user-requested from real workflow friction**:
- **Single-site orientation sampling**: `--ml-rank` used to hard-require `--all-sites`
  (`parser.error`, "screening a single --site-index site has nothing to rank against") --
  too strict once orientation sampling exists, since a single chosen site's own sampled
  orientations ARE something to rank, even with `--all-sites` never given. The check now
  reads `if args.ml_rank and not args.all_sites and not orientation_flags_used:`, so
  `--ml-rank --site-index N --n-orientations-polar/--azimuthal ...` (no `--all-sites`)
  runs all requested orientations on that ONE chosen site instead of forcing every other
  symmetric site too -- `variants` (the `(slot, height, site_type, coord)` list the
  `--ml-rank` loop already iterates) needed zero changes for this, since a single-site
  selection was already just `n_selected=1` in that same list; only the validation gate
  was too strict. `stb_suite.py`'s `run_adsorb_setup()` menu wrapper was updated to
  match: the orientation-sampling prompt used to live only inside the `--all-sites`
  branch (`if all_sites: ... if ml_rank: ...`); it's now asked whenever `all_sites` is
  True AND `ml_rank` was chosen, OR whenever `all_sites` is False (single site) --
  answering "yes" in the single-site case is now what turns `ml_rank` on (no separate
  "pre-screen with MACE?" question there, since orientation sampling is the only reason
  `--ml-rank` does anything for one site). This also surfaced a real, latent args-
  construction bug: the pre-existing code only ever appended `--ml-rank`/orientation
  flags to `args` inside the `if all_sites:` branch, so a single-site `ml_rank=True`
  (now reachable) would have silently built a command line with NONE of `--ml-rank`/
  `--n-orientations-*` even though the user had just said yes to all of it -- fixed by
  moving the `if ml_rank:` block out to run unconditionally (after the `all_sites`/
  `site_index` args are appended), gating only `--top-k` (site-ranking, meaningless for
  one site) on `all_sites` specifically. Verified live both via the CLI directly
  (`--site-index 0 --ml-rank --n-orientations-polar 2 --n-orientations-azimuthal 2`,
  no `--all-sites`, succeeds and writes `site_1_ontop_orientN/` folders) and through the
  interactive menu (`4.8.1` -> single site -> orientation sampling -> yes), including
  that the plain `--ml-rank --site-index N` case with NO orientation sampling is still
  correctly rejected (`parser.error`, unchanged message).
- **Succinct per-orientation progress counter**: the per-orientation MACE relax loop
  previously printed nothing at all while it ran (only the post-dedup/top-k survivors
  got a `print_dual` line each) -- fine for a small 2x2 grid, invisible/silent-feeling
  for a large one (e.g. the docstring's own suggested 100-orientation, keep-the-best-10
  workflow). New `core/cli.py::print_progress_line(line, step, total)`/
  `finish_progress_line()` (extracted from `amorphize.py`'s pre-existing, tty-aware
  melt/quench progress helpers of the same design -- `amorphize.py`'s own
  `print_progress`/`finish_progress_line` are now thin wrappers calling the shared
  versions, zero behavior change there, same extract-on-second-use policy as the rest of
  `core/`): a self-overwriting `\r`-updated line on an interactive terminal, or a
  periodic full line (every ~10% of the total, plus the final step) when stderr is
  redirected (e.g. into a log file, exactly how this suite's own test scripts capture
  output) -- avoiding both "spam scrollback with hundreds of lines" and "produce nothing
  visible when piped." Called once per orientation inside `adsorb.py`'s relax loop
  (`site {i}/{n_sites} ({type}, h=...): orientation {j}/{n_orientations} relaxed (E =
  ... eV, converged/hit step cap)`), `finish_progress_line()` once the per-site
  orientation loop ends. The post-dedup "N orientation(s) sampled -> M unique kept"
  summary line (unchanged) still prints normally via `print_dual`.
- **`.xyz` trajectory of every generated orientation**: `--n-orientations-polar`/
  `--azimuthal` sampling now also writes `sites/orientation_trajectory<_adsorbate>.xyz`
  (multi-frame extended XYZ, `ase.io.write(..., format="extxyz")`) containing EVERY
  relaxed orientation candidate actually generated across every sampled site for that
  adsorbate -- not just the unique/`--orientation-top-k`-surviving ones written as real
  SIESTA folders -- so the full exploration (including near-duplicates/discarded, higher
  -energy attempts) can be visually reviewed as a trajectory/video in VMD/OVITO,
  something no static ranking table/plot can show. Each frame's `ase.Atoms.info` records
  `site`/`site_type`/`height`/`orientation`/`energy_eV`/`converged` (visible as
  extended-XYZ header key=value pairs, same convention `stb-mlmd`'s own `--out-format
  xyz` already uses for per-frame metadata). **Real bug found and fixed while
  implementing this**: writing straight from the just-relaxed `ase.Atoms` (still
  carrying its live MACE `SinglePointCalculator`/`MACECalculator` results) raised
  `KeyError: 'key from calculator already exists in atoms.info'` -- `ase`'s own extxyz
  writer tries to fold the ATTACHED CALCULATOR's own results (`energy`, etc.) into
  `atoms.info` too, colliding with the `energy_eV`-etc. keys this tool sets manually one
  line earlier. Fixed by detaching the calculator (`ase_atoms.calc = None`) right after
  reading `get_potential_energy()` and before appending to the frame list -- harmless,
  since nothing downstream (`AseAtomsAdaptor.get_structure`, `min_adsorbate_slab_
  distance`) needs the attached calculator, only the already-extracted `energy`/
  `converged` values. Verified live: a 4-orientation single-site run wrote exactly 4
  frames (matching `orient_results`' full count, before any dedup/top-k trimming), each
  with the correct per-frame `site=`/`site_type=`/`energy_eV=` header fields.
- **`--orientation-top-k` already covers "choose how many folders per site, ranked by
  energy"**: no new code needed -- its existing docstring/behavior (`kept =
  kept[:args.orientation_top_k]`, applied AFTER `deduplicate_orientations`) already
  means exactly "of the N unique orientations found for a site (e.g. 100 sampled -> M
  unique), keep only the `--orientation-top-k` lowest-energy ones as real SIESTA
  folders" -- verified live it matches the requested "100 generated, pick the best 10"
  workflow shape exactly (checked with a small 2x2=4 grid and `--orientation-top-k 1`:
  exactly 1 folder written, while `orientation_trajectory.xyz` still carries all 4).

**Real, latent gap found (user-requested check, not a crash report): a MACE relax under
`--ml-rank` could write a fractional coordinate outside `[0, 1)`.** The write-time
`wrap_into_cell` fix earlier this session only ever covered the plain (non-`--ml-rank`)
`add_adsorbate` path -- the `--ml-rank`/orientation-sampling path builds
`relaxed_struct = AseAtomsAdaptor.get_structure(ase_atoms)` straight from the
just-relaxed `ase.Atoms` and writes it unwrapped. The substrate is `FixAtoms`-frozen
during this relax, but nothing constrains the free adsorbate atom(s)' motion to stay
inside the cell -- an in-plane (periodic, non-vacuum) direction is exactly where a
FIRE-relax trajectory starting near a cell edge (e.g. an `ontop` site sitting on a
corner atom, frac x~0/1) can legitimately slide across the boundary while converging,
same underlying risk `wrap_into_cell` was already written to guard against, just
triggered by relaxation motion instead of the raw `add_adsorbate` placement. Fixed by
wrapping `relaxed_struct` (via the same, already-verified `wrap_into_cell` from
`neb.py`) immediately after building it, before `min_adsorbate_slab_distance` is
computed on it (unaffected either way -- pymatgen's `distance_matrix` is periodic-aware,
so wrapping before or after doesn't change that value) -- covers `--ml-rank` with or
without orientation sampling, since both paths funnel through the same `kept` loop.
Verified two ways: (1) live on the suite's own graphene+H fixture with a deliberately
large `--height 19.5` (a genuine relaxation, no fabrication) -- the H atom converged to
frac (0.9999999998, 0.9999999997), a hair below 1.0 from real floating-point noise at
the exact 0/1 boundary (the identical `0.9999999997712383`-class artifact already
documented for the write-time fix, printed as `1.00000000` by `write_fdf`'s own `.8f`
rounding -- cosmetic, not a bug); (2) deterministically, by monkeypatching
`mace_relax.relax` to shift the relaxed adsorbate atom by exactly 1.5 lattice vectors
along both in-plane axes after a real relax completes -- the written `structure.fdf`
correctly shows frac (0.5, 0.5), the same physical position with only the whole-lattice
-vector part removed (not, say, silently clamped to 0.0/1.0, which would have discarded
the genuine +0.5 displacement and misrepresented the physics).

**Follow-up, user-reported: the wrap above was still too late** -- placed on
`relaxed_struct` right before it fed into the site folder write, it never touched
`orientation_frames` (the `orientation_trajectory.xyz` export) or the `ase_atoms`
objects `deduplicate_orientations` compares by RMSD, both built earlier in the same
loop from the SAME `ase.Atoms` object the wrap only ran on later. The user opened
`orientation_trajectory.xyz` directly (VMD) and the written structures via `ase view`
and saw molecules sitting outside the cell box in both. Fixed by moving the wrap to
the earliest possible point instead: `ase_atoms.wrap()` (ASE's own in-place, all-3-axes
Cartesian wrap -- same canonicalization `wrap_into_cell` does via pymatgen, just
operating directly on the `Atoms` object every downstream consumer in this loop already
holds a reference to) is now called immediately after each orientation's relax
completes and its energy is read off, before `orient_results.append(...)` -- so
`orientation_frames.append(ase_atoms)`, `deduplicate_orientations`'s RMSD (now
comparing consistently-wrapped candidates instead of an inconsistent mix, a latent
correctness improvement to the dedup itself, not just a display fix), and the later
`relaxed_struct = AseAtomsAdaptor.get_structure(ase_atoms)` used for the site folder
all see the same, already-wrapped positions -- the redundant `wrap_into_cell()` call at
that last step was removed as now-dead code. Verified live: monkeypatching
`mace_relax.relax` to shift every free (adsorbate) atom by 1.5 lattice vectors along
both in-plane axes, with H2O + orientation sampling on, confirmed BOTH
`orientation_trajectory.xyz`'s frames (every sampled orientation, not just the kept
ones) and the written `structure.fdf` land inside `[0, 1)` at the correct wrapped
fractional position.

`core/adsorption_sites.py::write_reference_folder` gained an opt-in `force_spin=False`
parameter (new `SPIN_POLARIZED_BLOCK` constant, same file as `FIXED_CELL_BLOCK`), and
`stb-nebSites` now passes `force_spin=True` by default for `site_A`/`site_B` (new
`--force-spin`/`--no-force-spin` flags, `--force-spin` default ON, same
`store_true`/`store_false`-sharing-a-`dest` pattern as `stb-neb`'s own
`--ml-freeze-substrate`/`--no-ml-freeze-substrate`) -- a single adsorbate atom (or a
molecule containing one) bonded to a slab commonly leaves the combined system with a net
magnetic moment (most simply, an odd total valence-electron count, which a spin-restricted
SCF cannot represent correctly at all), the same reasoning `adsorb.py`'s own
`force_spin_polarized()` already applies to its isolated-adsorbate reference -- costs
nothing for a genuinely closed-shell site (converges to zero moment). Implemented via
`config_extra.fdf` (not by editing `--calc`'s own `Spin` tag directly, unlike
`adsorb.py`'s version): `Spin polarized` is appended there alongside the existing
`MD.VariableCell false`, and since `config_extra.fdf` is `%include`'d at the very top of
the written `calc.fdf` (`structure_io.prepend_include`) and SIESTA's fdf reader is
first-occurrence-wins for a duplicate label, this correctly overrides any `Spin
non-polarized` (or absent tag) in the user's own `--calc` template regardless of where it
appears. `force_spin` defaults to `False` in `write_reference_folder` itself, so
`adsorb.py`'s existing 3 call sites (`clean_slab_dir`, the isolated-adsorbate `ads_dir`,
and `site_dir`) are unaffected -- deliberately scoped to `stb-nebSites` only, matching
what was asked, though the same reasoning would apply equally to `stb-adsorb`'s own
`site_dir` if a future need arises. Verified live: default run's `config_extra.fdf` gets
`Spin                polarized` and the written `calc.fdf`'s first non-empty line is the
`%include config_extra.fdf`, confirmed ahead of the template's own (contradictory) `Spin
non-polarized` line; `--no-force-spin` leaves `config_extra.fdf` with no `Spin` line at
all, template untouched.

**`stb-adsorb --ml-rank` (item 4.8.1) fixed a real transparency gap, user-reported**: the
relaxed geometry it computes (MACE-MP-0, substrate fixed, adsorbate free -- so the
adsorbate's height DOES already move away from the `--height` initial guess as part of
this relax) was always what actually got written to the `sites/site_*/` SIESTA folders,
but neither fact was visible anywhere in the report -- the `[3] ML PRE-SCREENING` table
only ever showed the pre-relaxation `Height` guess, and nothing stated that `[4] WRITING
SITE FOLDERS` uses the relaxed result rather than that guess. Fixed by computing each
candidate's actual post-relax adsorbate-slab distance (`min_adsorbate_slab_distance` on
`relaxed_struct`, already imported/defined in `adsorb.py`) and adding it as a new
`Relaxed dist (Ang)` table column alongside the renamed `Init. h (Ang)` column, plus an
explicit note in `[4]`: "every site folder below is written from the MACE-MP-0-relaxed
geometry ... not the pre-relaxation 'Init. h' guess". The `scored` tuple grew one element
(`relaxed_dist` appended last); every unpacking site (`plot_ml_ranking`'s two list
comprehensions, `[4]`'s folder-writing loop) updated to match. `site_table.txt` (the
machine-readable file `stb-adsorbAnalysis` parses) deliberately left untouched --
its own comment already warns the column format/order is depended upon downstream, and
this fix is about console/report visibility, not that file's contract. Verified live on
the suite's own graphene+H fixture: 4 candidates that all started from the same `h=2.00`
guess relaxed to visibly different distances (1.211-3.097 Ang), each one correctly
reflected in the corresponding `structure.fdf`'s actual written adsorbate position.

**`stb-adsorb` (item 4.8.1) closed two real physics-config gaps, user-reported**: spin
polarization was inconsistent across the folders it writes (`adsorbate*/` was already
forced `Spin polarized` via `force_spin_polarized`, but `clean_slab/` and every
`sites/site_*/` silently used `--calc`'s own `Spin` setting unmodified -- the suite's
own reference template declares `Spin non-polarized`, so in practice the isolated
reference and the combined site -- the quantity `E_ads` actually depends on -- were
computed at different levels of theory), and SIESTA's slab dipole correction
(`Slab.DipoleCorrection`) was never set anywhere in this tool, despite `stb-her`/
`stb-oer`/`stb-oerIntermediates` already treating the identical physics ("adsorbing on
only one face breaks the slab's inversion/mirror symmetry along the vacuum axis, giving
the cell a net dipole along a PERIODIC direction") as "structurally required," forced
unconditionally on every site folder they write.

- **Spin**: new `--force-spin`/`--no-force-spin` flags (default ON, same
  `store_true`/`store_false`-sharing-a-`dest` pattern as `stb-nebSites`'s own flags of
  the same name) extend `force_spin` (the `core/adsorption_sites.py::
  write_reference_folder` parameter added for `stb-nebSites`, `config_extra.fdf`
  -based, `%include`d at the top so it overrides `--calc`'s own `Spin` tag) to every
  `sites/site_*/` folder too. `clean_slab/` deliberately does NOT get it (no adsorbate
  present, no clear universal reason to force it -- same reasoning already documented
  on `SPIN_POLARIZED_BLOCK` itself). The pre-existing `[NOTE]` (previously gated on
  `len(mol) == 1`, claiming the combined calc.fdf is "used as given, unmodified") is
  now unconditional per adsorbate (the underlying `force_spin_polarized` on the
  isolated reference already was) and its wording branches on `args.force_spin` to
  stay accurate either way.
- **Dipole correction**: new `DIPOLE_CORRECTION_BLOCK` constant in
  `core/adsorption_sites.py` (`Slab.DipoleCorrection      .true.`) and a new
  `force_dipole=False` parameter on `write_reference_folder` itself -- same
  `config_extra.fdf`-based mechanism as `force_spin`, deliberately NOT the calc_text
  -editing `regex.subn`-style `force_dipole_correction` approach `her.py`/`oer.py`/
  `oer_intermediates.py` each already have their own copy of (first implementation
  pass here mirrored theirs via a new `core/calc_directives.py` copy, then rewritten
  to the config_extra.fdf approach instead once asked to match how `force_spin`/
  `MD.VariableCell false` already work in this same function -- the
  `core/calc_directives.py` copy was removed again, unused). Applied unconditionally
  (no flag -- matching HER/OER's own unconditional treatment of the same physics, and
  free when the true dipole is already zero, since the correction itself evaluates to
  zero) via `force_dipole=True` on the `write_reference_folder` calls for
  `clean_slab_dir` and `site_dir`. Deliberately NOT applied to `adsorbate_dir` (the
  isolated-adsorbate reference) -- a molecule in an all-around vacuum box, not a
  2D-periodic-plus-vacuum slab, so `Slab.DipoleCorrection` has no physical meaning
  there.

Verified live: `clean_slab/config_extra.fdf` and every `sites/site_*/config_extra.fdf`
contain `Slab.DipoleCorrection      .true.` (their `calc.fdf` itself carries neither
this nor a `Spin` line directly -- both live only in `config_extra.fdf`);
`sites/site_*/config_extra.fdf` also contains `Spin polarized` by default, with
`%include config_extra.fdf` confirmed as the actual first line of the written
`calc.fdf` (so it correctly overrides the template's own `Spin non-polarized` later in
the same file); `clean_slab/config_extra.fdf` never gets a `Spin` line;
`adsorbate/config_extra.fdf` gets neither override (its own spin-polarization is
forced the older way, directly on `adsorbate_calc_text`, unchanged); `--no-force-spin`
leaves `sites/site_*/config_extra.fdf` with no `Spin` override while the dipole
correction override still applies. Zero regression across the full
adsorption-workflow regression suite (`prep` 132/132, `bsse` 46/46, `analysis`
66/66), plus `stb-nebSites`'s own suites (`sites` 74/74) confirming
`write_reference_folder`'s new `force_dipole` parameter defaulting to `False` is a
no-op for its existing, unrelated `force_spin`-only caller.

**CPU/GPU device selection for every MACE-consuming tool** (`core/mace_relax.py`, all
32 tools that ever load a MACE model). Every one of these tools now exposes a
`--device {cpu,cuda}` flag (named `--ml-device`/`--mace-device` instead, matching each
tool's own pre-existing naming convention, on `adsorb.py`, `neb.py`,
`stackingfault.py`, `oer_intermediates.py`, `phonons_create.py`, `xrdsearch.py`,
`crystalcast.py`), defaulting to `cpu`, plumbed straight through to
`core/mace_relax.py::get_calculator`'s existing `device` parameter -- the single
chokepoint every MACE-consuming tool in the suite already calls through, so GPU
validation is implemented exactly once and every caller inherits it for free just by
threading its own `device` argument to that one function:
- `core/mace_relax.py::gpu_available()` -- returns `(available, detail)`, checking
  `torch.backends.cuda.is_built()` (was PyTorch compiled with CUDA support at all)
  then `torch.cuda.is_available()` (is a CUDA-capable GPU actually visible right now),
  `detail` being the failure reason or (on success) `torch.cuda.get_device_name(0)`.
- `core/mace_relax.py::resolve_device(device)` -- a no-op for `"cpu"`; for `"cuda"`,
  calls `gpu_available()` and raises a clear `ValueError` (pointing at
  https://pytorch.org/get-started/locally/) if no usable GPU is found, rather than
  letting PyTorch fail with its own less-actionable internal error deeper in model
  loading. `get_calculator` calls this as its very first line, before any model
  loading happens.
- Every interactive `stb_suite.py` wrapper for these 32 tools gained a matching
  `Device [cpu/cuda, default: cpu]:` prompt inside its existing advanced-settings /
  model-selection block (right after the model-size prompt, or gated behind the same
  `if <ml_option>:` conditional as the rest of that function's MACE-specific prompts
  when MACE usage is itself opt-in, e.g. `run_defect_generator`'s two separate
  `--ml-rank`/`--ml-relax` gates) -- proactively runs `gpu_available()` itself and
  prints `[OK] GPU detected: <name>` or a `[WARNING]` before ever invoking the tool, so
  a bad `cuda` choice is caught in the menu rather than surfacing as a subprocess
  failure.

Verified live: `gpu_available()` correctly reports a real GPU's name (RTX 3060) on a
machine that has one; mocking `torch.cuda.is_available` to `False` confirms
`resolve_device("cuda")` raises the expected `ValueError`.

**`stb-adsorb` (item 4.8.1) is a deliberate exception to `resolve_device`'s raise
-on-missing-GPU default, per explicit user request**: `--ml-rank`/`--ml-prerelax` are
only ever a fast pre-screening/pre-relax step ahead of the real SIESTA jobs the rest of
the tool writes, so a missing GPU aborting the whole run (as every other MACE-consuming
tool still does, unchanged) was judged more disruptive than useful here. New
`adsorb.py::resolve_ml_device(device, mace_relax_module)` checks `gpu_available()`
itself and returns `("cpu", <warning text>)` instead of letting `get_calculator`'s own
`resolve_device` raise, so both `--ml-prerelax` and `--ml-rank` degrade to CPU with a
`[WARNING]` in the report rather than erroring out. Returns the warning text instead of
printing it directly: both call sites resolve this from inside a
`capture_library_noise()` block, whose stdout redirect would otherwise silently swallow
a `print_dual` call made in there, so the print happens just after the `with` block
exits. Verified live on a GPU-less machine: `stb-adsorb --ml-rank --ml-device cuda`
exits 0, prints `[WARNING] --ml-device cuda requested, but no CUDA-capable GPU was
detected ... -- falling back to CPU.`, and still writes the site folder normally.

**`stb-adsorb` (item 4.8.1), user-reported: the adsorbate could be written to
`sites/site_*/structure.fdf` with a fractional coordinate outside `[0, 1)`** (negative,
or >= 1). Root cause confirmed against pymatgen's own installed source:
`AdsorbateSiteFinder.add_adsorbate`/`adsorb_both_surfaces` place the adsorbate atom(s)
via `struct.append(specie, ads_coord + site.coords, coords_are_cartesian=True)` with no
`to_unit_cell=True`, so the resulting fractional coordinate is never normalized -- and
`core/structure_io.py::from_pymatgen` writes `site.frac_coords` straight through with no
wrap either. All 3 possible origins of a site's structure (plain `add_adsorbate`, the
`--ml-rank` MACE-relaxed structure, `--both-sides`) funnel through the same site-writing
loop in `main()` before `write_reference_folder` is called, so a single fix there covers
all of them: wraps `ads_struct` via `wrap_into_cell` -- already reused (not duplicated,
same policy as everywhere else in `core/`) directly from `neb.py` (`from stb.neb import
wrap_into_cell`), the same function `mlneb.py`/`mladsorb.py`/`mldiffusion.py` already
import from there. `min_adsorbate_slab_distance`'s own `distance_matrix` is already
periodic-aware, so wrapping first (before that overlap check, in the same loop) doesn't
change its result. Scope deliberately limited to `stb-adsorb` alone, per explicit user
request -- `stb-nebSites` has the identical latent issue (same `add_adsorbate` +
`write_reference_folder` pattern) but was left unfixed for now. Verified live with a
deliberately large `--height 15` on the suite's own 20 Ang-cell graphene fixture (slab
top at frac z=0.5, so the unwrapped adsorbate would land at frac z=25/20=1.25): the
written `structure.fdf` correctly shows z=0.25. A second, weaker case surfaced while
building this regression test: a raw frac coordinate that's mathematically just below
1.0 by floating-point noise (e.g. `0.9999999997712383`, from `--ml-rank`'s MACE relax
landing an ontop adsorbate almost exactly back at its starting xy) is legitimately
wrapped and correct, but `write_fdf`'s own `.8f` formatting rounds it to the printed
string `"1.00000000"` -- not a bug (physically indistinguishable from 0.0, ~1e-8 Ang),
but the regression test's own tolerance (`test/4-workflow/8-adsorption/prep/test.sh`)
had to widen from a strict `< 1.0` to `<= 1.0 + 1e-6` to avoid flagging this harmless
display-rounding artifact as a false failure.

**Follow-up, user-reported: `sites/adsorption_sites.png` still showed candidate markers
outside the highlighted unit cell**, even after the write-time fix above -- a different
code path entirely (`core/adsorption_sites.py::write_site_plot`, the generic candidate
-overview plot, was never touched by that fix, and only ever draws pymatgen's OWN
`find_adsorption_sites()` candidates at pymatgen's own defaults, not this run's actual
selected/written sites). The real constraint, already documented on
`cluster_candidate_coords`'s own docstring from an earlier (2026-08-06) investigation:
a group of symmetrically-equivalent candidates that are mutually close to each other can
provably straddle the cell boundary in a way that makes "mutually close" and "fully
inside one fixed `[0,1)` box" mutually exclusive for a SINGLE chosen representative
point per candidate -- an attempt to force it anyway was tried and reverted back then.
A first attempt this round (`periodic_marker_copies`, since replaced) sidestepped that
by drawing EVERY periodic image of every candidate visible in the plotted window
(matching `plot_slab`'s own `repeat=5` tiling of the underlying atoms) -- correct, but
visually noisy (many duplicate "x" marks scattered well outside the cell too) and not
what was actually asked for. Per explicit follow-up request ("mark only the ones inside
the unit cell"), replaced with **`core/adsorption_sites.py::wrap_markers_into_cell
(cart_sites, lattice_matrix, symm_op)`**: exactly one marker per candidate, its unique
periodic image whose in-plane fractional coordinate (solved via `inv([a2d, b2d]) @ xy`
in the plot's own 2D rotated frame, then `% 1.0` on both components -- the 2D-projected
analog of `neb.py::wrap_into_cell`'s `frac_coords % 1.0`, needed separately since
candidates here are bare Cartesian points, not full pymatgen Structures) lands in
`[0,1) x [0,1)`. This deliberately reintroduces `cluster_candidate_coords`'s own
documented tradeoff in the OTHER direction: two genuinely-adjacent-across-a-cell-seam
candidates (e.g. frac x=0.001 and x=0.999) each wrap correctly but independently, so
they can still end up near opposite edges of the same drawn cell -- accepted as the
right tradeoff for this specific plot now that "only inside the cell, no clutter" is
the stated priority. `write_site_plot` calls this in place of `cluster_candidate_coords`
for its own single scatter call; `cluster_candidate_coords` itself is unchanged and
still used everywhere else it always was (site-selection distance logic in
`adsorb.py::main()`, `stb-nebSites`'s "closest candidate pair" suggestion), since those
need one real, distance-comparable coordinate, not a plotting convenience. Benefits
`stb-mladsorb` for free too (`mladsorb.py` reuses the same `write_site_plot`);
`stb-nebSites` is unaffected, it has its own dedicated, unrelated plot functions.
Verified live on the suite's own graphene fixture (`--site-type all --all-sites`, 4
candidates): all 4 markers now render exactly once each, all strictly inside the drawn
cell parallelogram, no duplicates outside it. Unit test in
`test/4-workflow/8-adsorption/prep/test.sh` (6h): a candidate deliberately placed
outside `[0,1)` (frac x=-0.3 in the rotated in-plane basis) comes back wrapped to
alpha~0.7, one that's already inside comes back unchanged, and exactly one marker is
returned per input candidate.

**`stb-adsorbAnalysis` (Stage 3) gained physical diagnostics, a "what to run next"
pointer, and an opt-in Gibbs-free-energy prep step; new `stb-adsorbGibbs` (Stage 4,
item 4.8.4) computes the actual DG.** User-requested three-part addition to the
Adsorption workflow (item 8):

- **`[2b] PHYSICAL DIAGNOSTICS`** -- per site, dipole magnitude (reusing the
  already-existing `core/siesta_log.py::get_electric_dipole`, no new parser needed)
  and net spin/magnetic moment (new `core/siesta_log.py::get_spin_moment`, same
  "last matching line wins, never raises" convention as every other parser there --
  **flagged as unverified**: no spin-polarized `.out` fixture exists anywhere in this
  repository to confirm SIESTA's exact wording against, tries two plausible patterns
  `"Total spin polarization"`/`"Total spin moment"`, degrades to `None`/`"--"` if
  neither matches a real run -- needs live confirmation), and a bond-length-change
  column (`min_adsorbate_slab_distance`, imported directly from `adsorb.py` -- this
  workflow's own established cross-import convention, e.g. `adsorb_bsse.py` already
  imports `read_site_table` from `adsorb_analysis.py` the same way -- comparing the
  pre-relaxation guess `structure.fdf` against the relaxed `.XV`/`read_site_geometry_
  atoms` result).
- **`[5] SUGGESTED NEXT ANALYSES`** -- always-printed, references `apply_source_label`
  (the same winning-site variable `[4] APPLY` already resolves) with conditional
  ("if you also generated `<file>`...") pointers at `stb-bader` (needs `.RHO`),
  `stb-dos` (needs `.PDOS.xml`), `stb-workfunction` (needs `.VT` -- also a natural
  dipole-correction consistency check, since `[2b]` already reports the dipole
  `Slab.DipoleCorrection` is meant to fix), and `stb-coop` (needs a full-BZ `.WFSX` +
  `.HSX`) -- none of these SIESTA outputs are things `stb-adsorb` itself guarantees
  exist, so every suggestion is conditional rather than assumed.
- **`[6] GIBBS FREE ENERGY (DG) PREP`**, opt-in via `--compute-gibbs` (off by
  default) -- the one deliberate exception in this tool to "analysis never generates
  new folders" (explicit user design decision, recorded here so it doesn't get
  "fixed" back by accident later): writes the vibrational-Hessian displacement
  folders the new Stage 4 needs, for the winning site (`apply_source_label`) and its
  isolated-adsorbate reference, once BOTH are confirmed relaxed (`core/structure_io.py
  ::read_relaxed_or_input`, a hard `[ERROR]` -- not silently falling back to the
  unrelaxed guess -- if either lacks a finished `siesta.XV`, since a Hessian at a
  non-stationary point isn't physically meaningful). `--zpe-mode {local,full}`
  (default `local`, confirmed with the user) mirrors `stb-herRefs`/`stb-oerRefs`'s
  own dual-mode design:
  - `local` (default): a partial-Hessian, decoupled-oscillator finite-difference
    calculation -- new `write_local_hessian_folders` (duplicated/adapted from
    `oer_refs.py::write_local_zpe_folders`'s own N-atom generalization of
    `her_refs.py`'s single-atom version -- each workflow stays self-contained,
    re-reading/duplicating rather than cross-importing a SIBLING workflow's module,
    exactly the same convention HER/OER already follow between themselves) --
    displaces only the adsorbate's own atoms (substrate frozen) for the site side, 6
    folders per atom (`disp_001..disp_{6N}`), plus a `gibbs_local_meta.json` sidecar
    (atom indices/symbols/displacement/order/SystemLabel). **No Phonopy involved at
    all** in this mode.
  - `full`: a real Phonopy phonon calculation of the ENTIRE site structure AND the
    clean slab (needed so Stage 4 can subtract the clean slab's own phonon ZPE/
    entropy, same physical reasoning `her_analysis.py`'s own module docstring already
    documents) -- reuses `core/phonon_workflow.py::build_phonon_displacements`/
    `write_displacement_folders` directly, near-verbatim adaptation of
    `her_refs.py`'s own full-mode block.
  - **The isolated-adsorbate reference ALWAYS uses the `local` method, regardless of
    `--zpe-mode`** -- confirmed as the right simplification with the user: a molecule
    in an all-around vacuum box has no periodic substrate phonon to subtract, so
    Phonopy's supercell machinery buys nothing there; a plain full-molecule Hessian
    (every atom of the isolated reference, not just "the adsorbate's atoms" since
    there's no substrate to freeze) is the correct, and cheaper, choice either way.
  - Site/reference `config.fdf` correctness: the site's Hessian folders inherit the
    site's own `Spin polarized`/`Slab.DipoleCorrection` (new `read_site_theory_flags`,
    duplicated -- not cross-imported, same self-contained-workflow reasoning as
    above -- from the identical fix already made for `stb-adsorbBsse`'s own ghost
    fragments earlier this session) rather than silently reverting to the raw
    `--calc` template's own (often non-polarized) `Spin` setting -- the exact class of
    bug already found and fixed once for BSSE, deliberately not reintroduced here.
    The isolated reference's own `calc.fdf` already has `Spin polarized` baked
    directly into its text (from `adsorb.py`'s `force_spin_polarized`, applied
    unconditionally to every isolated-adsorbate reference regardless of `--force-
    spin`) and no `Slab.DipoleCorrection` at all (correctly -- not a slab) -- verified
    live that its Gibbs folders inherit exactly that, unchanged.

New **`stb-adsorbGibbs`** (`adsorb_gibbs.py`, Stage 4, item 4.8.4, deliberately IS
registered in `stb_suite.py`'s `WORKFLOW_TOOLS` -- unlike `stb-nebCycle`, this one has
a normal interactive menu entry, no cluster-submission-loop reason to omit it) reads
whichever single `<site>/`+`<adsorbate>_isolated/`(+`clean_slab_full/` for `full`
mode) folder set exists under `<dir>/gibbs/` (auto-detected by name pattern, no
site label needs to be passed explicitly -- `--compute-gibbs` only ever preps ONE
site's worth of folders per run), rebuilds the mass-weighted `3N x 3N` dynamical
matrix from `.FA` force files via `core/siesta_log.py::read_fa_forces` (duplicated
`read_fa_force` single-atom helper, same "each workflow re-reads its own files"
policy as `her_analysis.py`/`oer_analysis.py`), diagonalizes it (reusing the exact
same `_FREQ_CONVERSION_THZ`/`_EV_PER_THZ`/`_BOLTZMANN_EV_K` constants and
harmonic-oscillator ZPE/entropy formulas as HER/OER), and combines:
```
DG(T) = E_ads (BSSE-corrected if complete, else raw) + [ZPE(site) - ZPE(isolated ref)]
        - T * [S(site) - S(isolated ref)]
```
swept over `--tmin`/`--tmax`/`--tstep`, plotted as `adsorption_gibbs.png`. Reuses (not
duplicates) `adsorb_analysis.py::read_site_table`/`read_bsse_energy`/
`GIBBS_LOCAL_META_FILE` directly -- a legitimate cross-import since `adsorb_gibbs.py`
is a brand-new file that only ever imports FROM `adsorb_analysis.py`, never the
reverse, so no cycle is created (unlike why `adsorb_analysis.py` itself had to
duplicate rather than import `adsorb_bsse.py`'s `read_site_theory_flags`: THAT
direction already has `adsorb_bsse.py` importing `read_site_table` FROM
`adsorb_analysis.py`, so importing back would cycle).

**Documented, permanent limitation** (printed in every Stage 4 report, not just this
file): the isolated-adsorbate reference's own ZPE/entropy term is VIBRATIONAL/
HARMONIC ONLY -- it does NOT include the translational/rotational entropy a real free
gas molecule has (an ideal-gas partition-function term this tool does not compute).
For a single-atom adsorbate this is exact (verified live: a genuinely free atom with
zero force at every displacement gives EXACTLY `zpe_ev = 0.0`/`ts_ev = 0.0` once its 3
pure-translation modes are correctly excluded by the existing near-zero-frequency
filter -- the right physical answer, a free atom has no vibrational modes at all); for
a polyatomic adsorbate it's an approximation, in the same spirit as (though computed
differently from) `stb-herRefs`/`stb-oerRefs`'s own use of a FIXED literature entropy
for their H2/H2O gas-phase references -- generalized here since there is no literature
table for an arbitrary adsorbate.

Verified end-to-end live (not just unit-level): built a real site + isolated O
reference via `stb-adsorb`, fabricated relaxed `siesta.XV` for both (same sisl-write
recipe `../bsse/test.sh` already uses), ran `stb-adsorbAnalysis --compute-gibbs` for
real, then fed the resulting `gibbs/` folders an IDENTICAL synthetic isotropic
harmonic force constant (`k=20 eV/Ang^2`) on both the site's O atom and the isolated
reference's O atom -- same mass, same force constant, so DZPE/DTS must cancel to
EXACTLY zero by construction, making DG land EXACTLY on the electronic-only `E_ads`
(`-0.100000 eV`) -- confirmed to the printed 4th decimal, a strong, exact-value
cross-check of the whole Hessian-to-DG pipeline, not just "did it run without
crashing". New tests: `test/4-workflow/8-adsorption/analysis/test.sh` (`2g`, the
`--compute-gibbs` flow + config-inheritance checks) and a new
`test/4-workflow/8-adsorption/gibbs/test.sh` (33 checks: the exact-cancellation DG
cross-check, the free-atom ZPE=0 physics check, and error handling for a missing/
ambiguous `gibbs/` layout or a missing `.FA` file).

**A real regression, self-discovered via the regression sweep this rollout required**
(not user-reported): `core/structure_io.py::write_fdf`'s species-id renumbering
(added earlier the same session to close a "non-sequential `%block
ChemicalSpeciesLabel`" gap -- SIESTA requires ids to be a gap-free `1, 2, 3, ...`)
renumbered species in `structure.species`' own iteration order (first-occurrence-in
-atoms order), not the RELATIVE order implied by each species' pre-existing
`species_meta['id']` -- silently reversible for any caller that assigns ids in an
intentional order independent of physical atom-list order. Caught running
`test/2-structures/12-crystalcast/test.sh`: `stb-crystalcast --species O Ni` (meaning
"O should be id 1") could come out with Ni as id 1 instead, purely depending on which
species happened to appear first in the generated structure's own atom list. Fixed in
both `write_fdf` and the identical bug in `from_pymatgen` by sorting species by their
OLD id before assigning fresh sequential ones (`sorted(species_with_atoms, key=lambda
s: int(structure.species_meta[s]['id']))`) -- this both closes any id gap AND
preserves whatever ordering intent the caller already encoded, instead of the two
being in tension. Verified directly: `--species O Ni` now keeps O as id 1 while still
renumbering sequentially; the original gap-fixing scenario (H:1, He:2 removed, O:3,
C:4 -> H:1, O:2, C:3) still works.

**A second real, unrelated bug found via the same regression sweep**:
`mladsorb.py`'s `main()` called `structure = resolve_slab_orientation(structure,
args.vacuum_gap)` without unpacking its return value -- `resolve_slab_orientation`
(in `core/adsorption_sites.py`) was changed, in an earlier session, from a
print-and-`sys.exit` function into one that returns `(structure, relabel_note)` (or
raises `ValueError`) so each of its two then-existing callers (`adsorb.py`,
`neb_sites.py`) could report the info/error line in their own report's exact style --
`mladsorb.py` was a pre-existing third consumer that was never updated for the new
signature, so `structure` ended up bound to a `(structure, note)` tuple, crashing the
very next line (`structure_io.to_pymatgen(structure)`) with `AttributeError: 'tuple'
object has no attribute 'atoms'` on literally every run. Fixed to match the other two
callers' pattern (`try: structure, relabel_note = resolve_slab_orientation(...) except
ValueError: fail(...)`, then print the note if present). Caught via
`test/5-mlsimulations/10-mladsorb/test.sh` going from 26/57 to 55/57 once fixed (the 2
remaining failures are unrelated -- `Rotation.random(rng=rng)` in `pick_best_orientation`
needs scipy >= 1.15's `rng` keyword; this environment has 1.13.1, which only accepts
`random_state`; confirmed via `git log -p` that this line predates any of this
session's changes).

## Domain conventions worth knowing

- Structure file formats handled throughout: POSCAR (VASP), CIF, FDF (SIESTA), XYZ
  (+ separate `.lattice` file), XSF, FHI-aims `geometry.in`, DFTB+ `gen`.
- SIESTA output files consumed by analysis tools: `.bands`, `PDOS.xml`, `.RHO`/`.VT`
  grid files, `.XV`, etc.
- Core dependencies: `numpy`, `ase`, `matplotlib`, `pymatgen`, `spglib`, `sisl`,
  `pybader`, `scipy`, `phonopy`. `ase` in particular is often imported inside a
  `try/except ImportError` since it's only required for CIF I/O in some tools.
