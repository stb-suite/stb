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

## Domain conventions worth knowing

- Structure file formats handled throughout: POSCAR (VASP), CIF, FDF (SIESTA), XYZ
  (+ separate `.lattice` file), XSF, FHI-aims `geometry.in`, DFTB+ `gen`.
- SIESTA output files consumed by analysis tools: `.bands`, `PDOS.xml`, `.RHO`/`.VT`
  grid files, `.XV`, etc.
- Core dependencies: `numpy`, `ase`, `matplotlib`, `pymatgen`, `spglib`, `sisl`,
  `pybader`, `scipy`, `phonopy`. `ase` in particular is often imported inside a
  `try/except ImportError` since it's only required for CIF I/O in some tools.
