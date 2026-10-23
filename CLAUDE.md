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
  full rebuilds. Never edit here and never treat it as source of truth.
- `stb-suite/pyproject.toml` — package metadata and the `[project.scripts]` table that
  maps each `stb-*` shell command to a module's `main()`.
- `stb-suite/meta.yaml` — conda-forge recipe (mirrors `pyproject.toml` dependencies).
- `test/` — example input files and manual smoke scripts organized by category,
  mirroring the `stb-suite` main menu 1:1 (`1-inputs`, `2-structures`,
  `3-analysis`, `4-workflow`, `5-mlsimulations`, `6-utils`). `4-workflow/<property>/`
  has `prep/` + `analysis/` (and sometimes further stage-named) subfolders for the
  multi-stage prep+analysis workflow properties. This is **not** an automated test
  suite (no pytest/unittest, no CI config exists). `test/6-utils/1-translator/test.sh`
  is representative: it generates sample structure files in every supported format,
  runs the built CLI against them, and greps the output for expected content.

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

## Versioning

The suite has ONE version, `1.9.<number of commits>` (`1.9` marks the development
series), computed by `stb-suite/src/stb/_version.py` and exposed as `stb.__version__`.
Every tool takes its `VERSION` from it (`from stb import __version__ as VERSION`), for its
`--version` and its banner -- **never give a tool a version literal of its own**. In a git
checkout it is recomputed at each import (`git rev-list --count HEAD`), so an editable
install is always current; an installed copy reports the number baked into its package
metadata when it was built (`pyproject.toml` reads `__version__`); `STB_VERSION`
overrides it, for builds without git history (the conda recipe uses this). Nothing is
edited on a push. The per-tool numbers used before (1.0.0-3.0.0) are retired; the notes
that sat next to them are kept in `stb-suite/CHANGELOG.md`. `bash test/version/test.sh`
checks all of this (needs the package and its `ml` extra).

## Extensions (plugins)

Tools that must not live in this public repository (restricted to a team, for example) ship
as a separate package that plugs into the menu, with no code copied here. `stb_suite.py`
calls `_load_plugins()` before it builds `TOOL_CODES`: every module registered under the
`stb.plugins` entry-point group may define any of `INPUT_TOOLS`, `STRUCTURE_TOOLS`,
`ANALYSIS_TOOLS`, `WORKFLOW_TOOLS`, `MLSIM_TOOLS`, `UTILITY_TOOLS` (same shape as here);
its entries are merged in by menu number and the tables are re-sorted, since the menu lists
in insertion order. A taken number or a plugin that fails to load only warns -- the menu
never fails because of an extension. An extension's entries run their tools like any other,
through `run_tool` and the console commands its own package installs; `stb.core.menu`
(`run_tool`, `prompt_pseudo_source`) is the API it imports.

The documentation tooling reads extension repositories too: `STB_DOCS_EXTRA_ROOTS` (paths
separated by `os.pathsep`), each with a `docs_plugin.toml` naming its menu module,
`pyproject.toml`, help snapshot and examples folder (see `docs/stbdocs/catalog.py`);
`docs/dump_help.py --pyproject ... --out ...` writes such a snapshot. Extension pages get no
"edit this page" link. The public site is built without the variable, so it never shows them.
`bash test/menu/test.sh` and `test/docs/test.sh` (with `test/docs/fixture_plugin/`) cover both.
**Do not name or describe an extension's implementation in this repository.**

## Architecture

**Each module in `stb-suite/src/stb/` is a thin CLI script** (argparse `main()`,
registered as a console script in `pyproject.toml` under `[project.scripts]`, e.g.
`stb-kgrid = "stb.kgrid:main"`) that delegates shared logic to `stb-suite/src/stb/core/`.

When adding a new tool, mirror this pattern: a new standalone `stb/<name>.py` that
imports what it needs from `stb.core`, plus a new entry in `[project.scripts]` in
`pyproject.toml` (and a fixture under `test/` if practical). New format-specific
structure readers/writers (POSCAR, CIF, XYZ, XSF, FHI, DFTB) still live only in
`translate.py` — it's the sole consumer of those formats. If a second tool needs one,
move it into `core/structure_io.py` first rather than copying `translate.py`'s
implementation (this **extract-on-second-use** policy is the general rule for
`core/`: don't pre-extract a helper into `core/` until a second consumer actually
needs it).

### Core shared modules (`stb-suite/src/stb/core/`)

- `structure_io.py` — `read_fdf` honors the `LatticeConstant` unit (Ang, Bohr/au, nm) and
  stores `lattice_constant` in Ang, so `raw_lattice_vectors()` still returns the block as
  written; before, `LatticeConstant 1.0 Bohr` (phonopy's SIESTA writer) was read as Ang and
  such a cell came out 1.89x too large.
  The only `.fdf` reader/writer (`read_fdf`, `write_fdf`,
  `to_pymatgen`, `from_pymatgen`, `rewrite_fdf_lattice`, `rewrite_fdf_positions`,
  `read_relaxed_or_input`, plus thin accessors like `lattice_only`,
  `species_list`/`species_dict`). Use `raw_lattice_vectors()` instead of
  `lattice_only()`/`.lattice` when writing back into an existing file's `%block
  LatticeVectors` while leaving its `LatticeConstant` line untouched — using the
  wrong one double-applies the lattice constant. `rewrite_fdf_positions(source_path,
  new_positions, out_path)` replaces only `%block AtomicCoordinatesAndAtomicSpecies`,
  preserving everything else (basis, SCF, pseudopotential blocks) verbatim — used
  where many rattled/derived copies of one reference calculation are needed without
  hand-reconstructing the whole `.fdf` (`write_fdf` only writes a bare-minimum file).
  `replicate_fdf_cell(source, dims, out)` writes an NA x NB x NC supercell of a full
  calculation `.fdf` (same atom order as `ase.Atoms.repeat`, `NumberOfAtoms` updated).
  `read_relaxed_or_input(path)` resolves a bare `.fdf` (unchanged `read_fdf`
  behavior) or a directory into `(FdfStructure, used_relaxed)`, preferring a `*.XV`
  inside that directory over its `structure.fdf` guess when one exists, and
  cross-checks the `.XV`'s own per-atom species (not just atom count) against the
  expected species at each index — a same-species atom swap at write time is
  otherwise a silent, undetected atom-identity bug. `write_fdf`/`from_pymatgen`
  always renumber `%block ChemicalSpeciesLabel` ids sequentially (1..N), sorted by
  each species' prior id (so relative ordering intent, e.g. "species A must be id
  1", survives even when some species were removed) — SIESTA requires this block to
  be gap-free.
- `siesta_log.py` — parsers for SIESTA `.out` logs: `get_fermi_energy`,
  `get_cell_height`, `get_stress_tensor` (matrix block + Voigt fallback, eV/Å³),
  `get_stress_voigt_kbar` (Voigt-only, raw kBar — kept separate from
  `get_stress_tensor` because their callers do incompatible downstream math),
  `get_free_energy`, `get_outcell` (actual cell SIESTA used, from the last
  `outcell:` block), `get_dynamics_type`/`categorize_dynamics`, `find_out_file`
  (locates a generically-named log, e.g. `calc.out`), `get_md_trajectory`
  (per-MD-step cell/energy/temperature, correct even for variable-cell/NPT runs),
  `read_fa_forces` (all-atom `.FA` force reader), `get_electric_dipole`,
  `get_spin_moment` (tries two plausible SIESTA wordings; not yet confirmed against
  a real spin-polarized fixture). All return `None` on not-found/parse-error rather
  than raising — these are used in loops that scan many folders and must tolerate
  incomplete ones.
- `kspace.py` — `compute_monkhorts` (Monkhorst-Pack grid from lattice + target
  density; raises `ValueError` on zero cell volume), `detect_vacuum_axes`.
- `cli.py` — `COLORS`, `color_text()`, `show_intro(lines, delay=0.2)`, `get_input`/
  `get_float_input`/`get_int_input`, `print_progress_line(line, step, total)`/
  `finish_progress_line()` (tty-aware: self-overwriting `\r` line on an interactive
  terminal, periodic full lines when stdout/stderr is redirected — use this instead
  of a bespoke progress printer for any loop over many expensive candidates).
- `deps.py` — `require_sisl()`, `require_mace()`: the shared heavy-optional-dependency
  guards (print a consistent install hint and exit if missing).
- `symmetry.py` — `reduce_to_unitcell(structure, mode, symprec, angle_tolerance)`,
  wrapping pymatgen's `SpacegroupAnalyzer` (primitive/conventional/refined). The
  output's atom order and coordinate origin are never guaranteed to match the input
  in any mode — spglib rebuilds the cell from scratch and can pick any
  symmetry-equivalent origin. Not a bug; the crystal is the same, just relabeled.
- `passivation.py` — `passivate_dangling_bonds(structure, passivant, cutoff,
  bond_length)`: caps undercoordinated atoms along the missing-bond direction from
  local coordination geometry alone (vector-sum of existing-neighbor directions,
  negated). Only auto-caps single-missing-bond atoms; 2+ missing bonds are
  geometrically underdetermined and are reported, not guessed.
- `mace_relax.py` — `build_cell_mask(vacuum_axes)`, `get_calculator(model, device)`,
  `relax(atoms, calc, cell_mask, ...)`, `relax_neb(...)`: the shared MACE-MP-0
  load/relax/NEB logic. Callers must call `core.deps.require_mace()` themselves
  first — this module only imports `mace`/`ase.optimize`/`ase.filters` lazily inside
  each function, so merely importing the module doesn't force the heavy optional
  dependency to load. `get_calculator`'s `model` accepts a custom `.model` file path
  (autodetected via `os.path.isfile`) in addition to `"small"/"medium"/"large"`.
  `gpu_available()` (checks `torch.backends.cuda.is_built()` +
  `torch.cuda.is_available()`) and `resolve_device(device)` (no-op for `"cpu"`;
  raises a clear `ValueError` for `"cuda"` with no usable GPU) are the single
  chokepoint every MACE-consuming tool's `--ml-device`/`--mace-device {cpu,cuda}`
  flag (default `cpu`) goes through — implement device validation once here, not
  per-tool. `stb-adsorb`'s `--ml-rank`/`--ml-prerelax` are a deliberate exception:
  they degrade to CPU with a `[WARNING]` instead of raising, since MACE there is
  only a pre-screening step ahead of the real SIESTA jobs.
- `md_traj.py` — `read_static_lattice`/`read_frame_lattices` (per-MD-step cell,
  falling back `.out` → `.XV` → `.fdf`), `read_md_timestep_fs`, `unwrap_trajectory`
  (minimum-image PBC unwrap).
- `adsorption_sites.py` — shared by the adsorption/NEB-site tools:
  `resolve_slab_orientation`/`reorient_vacuum_to_c` (raises `ValueError` or returns
  `(structure, relabel_note)` — callers report the note/error in their own style,
  never prints-and-exits itself), `center_slab_in_vacuum` (recenters a slab that
  starts near a cell boundary so PBC wrap during relaxation can't teleport an atom
  to the opposite face and read as a spurious huge displacement),
  `write_reference_folder` (with opt-in `force_spin`/`force_dipole`/`force_vdw`
  parameters), `write_site_plot`, `wrap_markers_into_cell`,
  `generate_systematic_orientations`/`deduplicate_orientations`,
  `min_adsorbate_image_distance`. Also defines the `config_extra.fdf` directive
  blocks (see "config_extra.fdf" convention below).
- `eos_fit.py` — `fit_eos`/`normalize_eos_string`/`invert_pressure`: the
  `ase.eos.EquationOfState` wrapper (GPa conversion, R², out-of-range detection)
  shared by `stb-mleos` and `stb-eosAnalysis`. Kept out of any module that calls
  `core.deps.require_mace()` at import time, so a plain DFT-workflow tool never
  forces the heavy `ml` extra to load just to fit a curve.
- `mace_phonons.py` — `generate_ml_displacements`/`compute_force_constants`: the
  MACE-driven phonon physics shared by `stb-mlphonons` (its sole consumer today).
  `compute_force_constants` subtracts the model's forces on the UNDISPLACED supercell.
- `phonon_workflow.py` — **residual forces**: `write_perfect_supercell` writes
  `<disp dir>/perfect_supercell/` (the undisplaced supercell; not `disp-000`, every
  consumer globs `disp-*`), written by stb-phononsCreate, stb-ir and stb-raman by default
  (`--no-perfect-supercell`), and `load_phonon_with_force_constants` subtracts its `.FA`
  from every displacement's -- centred first, since phonopy removes each force set's drift (mean force per atom) when it builds FORCE_SETS -- (FORCE_SETS rewritten, original kept as FORCE_SETS.raw; warns
  when the folder is absent/unfinished) -- so every SIESTA phonon consumer (phonons, IR,
  Raman, adsorption Gibbs) benefits. A residual that is tiny per atom adds up over a
  sublattice: graphene relaxed to 0.004 eV/Ang, 5x5 supercell, gave Gamma E2g 31.9 instead
  of 41.0 THz. Load with `symprec=0.01` (phonopy's 1e-5 made graphene Cmmm and split LO/TO);
  stb-mlffValidate `--phonon-ref` now does both (the model's own residual subtracted too).
- `xrd.py` — pattern computation/experimental-pattern reading shared by `stb-xrd`,
  `stb-xrdsearch`/`stb-xrdrank`. `looks_like_peak_list` (coefficient-of-variation
  heuristic on point spacing) auto-detects a sparse peak-list-shaped
  `--compare-to` file and Gaussian-broadens it (`broaden_peak_list`, same
  FWHM/resolution as the simulated side) before scoring similarity — comparing a
  raw peak list directly against a Gaussian-broadened simulated profile is
  apples-to-oranges and silently deflates the score. `--raw-experimental` opts out.
- `neb_manifest.py` — `neb_manifest.json` read/write/`validate_manifest_pair`: an
  explicit atom-correspondence proof (`species_sequence`, shared `pair_id`) meant to
  be written by any tool that builds two structures guaranteed to share atom-index
  correspondence (no in-repo tool currently produces one — the former `stb-nebSites`
  did, before its symmetry-site-enumeration approach was removed from the NEB
  workflow — but the format stays as forward-looking infrastructure for a future
  producer). A consumer (`stb-neb`) that finds a valid manifest on both sides skips
  distance-based atom matching entirely instead of guessing; missing/inconsistent
  manifests fall back to the pre-existing distance-based (`--autosort-tol`) behavior.
- `mlff_workflow.py` — load/save of `mlff_workflow.json`, the state file the four MLFF
  stages (4.17.x) communicate through (schema in its docstring; paths relative to the file).
- `atomic_env.py` — `environment_fingerprints`/`EnvironmentSet`: model-independent
  atomic-environment fingerprints (sorted distances to the 16 nearest neighbours), novelty in
  units of a set's typical spacing; stb-mlffActiveLearning's novelty filter and stb-mlff's
  near-duplicate check.
- `model_card.py` — the **model card**: `mlff_model_card.json`, written by `stb-mlffAnalysis` next to every
  member model of a round (`write_round_cards` in `mlff_analysis.py`; `--make-cards [ROUND]` rewrites them for an
  already trained round without training; `stb-mlffValidate` refreshes each card's `validation` section and, on every run, rewrites `<run>/best_model/` -- best round by committee-mean test F RMSE, then its best member, model files + card with a `selection` section, recorded as `best` in the workflow JSON; the folder is overwritten each time, `export_best_model` in `mlff_validate.py`). It holds
  the model's SHA-256, level of theory (incl. whether the labels include DFT-D3), E0s/energy scale, training
  domain (elements, area/volume per atom of the periodic axes, shortest distance, force/energy/pressure ranges),
  recipe, metrics, stability/speed, usage notes and provenance. **A fine-tuned model is always handed over with
  its card**: every tool with `--custom-model` also has `--custom-model-card` (`--ml-custom-model-card`,
  `--mace-custom-model-card`, `--custom-model-cards` for stb-mlconvergence), added by
  `model_card.add_card_argument` and checked right after `parse_args()` by `model_card.setup_from_args`
  (flag missing, unreadable card, or SHA-256 mismatch -> clean error). `setup_from_args` registers the pair;
  `mace_relax.get_calculator` refuses a custom model file with no registered card (`allow_uncarded=True` only
  for the mlff stages themselves, which work from `mlff_workflow.json`), skips `dispersion=True` for a D3
  model, and returns a calculator that raises for an element outside the model's table and warns once per
  message when the cell size / shortest bond / predicted forces leave the trained range. Menu wrappers ask
  for the card with `core.menu.prompt_model_card`. A new tool that adds `--custom-model` must follow the same
  three calls (import, `add_card_argument`, `setup_from_args`).
- `ase_view.py` — `view_structure_interactive`: shared, graceful (no-display-safe)
  wrapper around ASE's interactive 3D viewer, used by every tool's `--view` flag.

### `stb_suite.py` (interactive dispatcher)

Shows a menu with 6 categories, in this order: **1 Inputs, 2 Structures, 3 Analysis,
4 Workflow, 5 ML Simulations, 6 Utils** — backed by six dicts, `INPUT_TOOLS`,
`STRUCTURE_TOOLS`, `ANALYSIS_TOOLS`, `WORKFLOW_TOOLS`, `MLSIM_TOOLS`,
`UTILITY_TOOLS`, each keyed by menu number to `{'title', 'description', 'func'}`.
`INPUT_TOOLS` holds only the 3 tools that configure an actual SIESTA run;
`STRUCTURE_TOOLS` holds everything that builds/generates/transforms a structure
file; `MLSIM_TOOLS` holds tools that RUN a simulation using a MACE potential instead
of driving a real SIESTA calculation (distinct from `WORKFLOW_TOOLS`, which
generates/analyzes DFT calculations, and from the one-shot ML structure-preprocessing
tools in `STRUCTURE_TOOLS` like `stb-mlrelax`/`stb-amorphize`). `WORKFLOW_TOOLS` is
one level deeper: each property has a `'stages'` dict of 2+ entries instead of a
`'func'`; `run_sub_menu()` recurses into `entry['stages']` automatically whenever it
finds one, so a multi-stage workflow needs no separate menu function.
`_flatten_tool_codes()` builds a flat `{"1.1": func, ..., "4.1.2": func, ...}`
lookup (`TOOL_CODES`) from these same dicts at import time — regenerate nothing by
hand, it derives entirely from the 6 dicts. Every `run_*` function builds an `args:
List[str]` from interactive prompts and dispatches via `run_tool(tool_name, args)`,
which shells out to the **installed** console command by name (`subprocess.run`,
must be on `PATH`) and centralizes error handling plus the "Press Enter to
continue" pause. **Don't reintroduce a second dispatch path** (resolving a sibling
script file via `__file__` and invoking it with `sys.executable`) — `run_tool()` is
the only dispatch path.

### Recurring conventions

New tools, and changes to existing ones, should follow these established patterns
rather than inventing a new one:

- **`config_extra.fdf`**: the standard mechanism for forcing calculation directives
  (`Spin polarized`, `Slab.DipoleCorrection`, `DFTD3 .true.`, `MD.VariableCell
  false`, ...) onto a written folder without editing the user's own `--calc`
  template text. A `config_extra.fdf` is written alongside `calc.fdf`, and
  `%include config_extra.fdf` is prepended at the very top of `calc.fdf`
  (`structure_io.prepend_include`) — SIESTA's fdf reader is first-occurrence-wins
  for a duplicate label, so this correctly overrides any conflicting setting later
  in the same file. See `core/adsorption_sites.py`'s block constants for examples.
- **Level-of-theory propagation**: any stage that reads a previous stage's already
  -relaxed geometry and derives a new calculation from it (a BSSE ghost fragment, a
  Hessian/Gibbs displacement folder) must read that stage's *actual*
  `config_extra.fdf` (`read_site_theory_flags`-style helper) to know whether
  Spin/dipole/vdW were applied there, rather than assuming a default — inheriting
  the wrong level of theory here has repeatedly produced physically-impossible
  results (e.g. a "more negative than uncorrected" BSSE correction).
  `read_site_theory_flags` is deliberately duplicated per self-contained workflow
  (adsorption, HER, OER) rather than cross-imported between sibling workflows.
- **Fractional-coordinate wrapping**: any structure written after a step that can
  move atoms (relaxation, MD, ML) must be wrapped into `[0, 1)` before writing —
  use `wrap_into_cell` (defined in `neb.py`, imported by `mlneb.py`/`mladsorb.py`/
  `mldiffusion.py`/`adsorb.py`) or `ase.Atoms.wrap()` directly on the `Atoms` object
  as early as possible (before it's used to build any comparison/trajectory/output
  downstream, not just at the final write). An unwrapped coordinate near a cell
  face is a genuine, silent atom-correspondence/PBC-jump hazard for anything that
  compares two snapshots of "the same" atom later (NEB endpoint matching,
  displacement diagnostics).
- **Trajectory output formats**: xsf (multi-frame AXSF, default), pdb (multi-model
  with `CRYST1`), xyz (extended XYZ, carries a `Lattice=` tag and optional per-frame
  metadata) — the shared 3-format convention across `stb-ani2traj`, `stb-mlmd`,
  `stb-mlphonons`, `stb-mlsearch`, all via ASE's native multi-frame writers
  (`ase.io.write` given a list of `Atoms`), never hand-rolled.
- **Plot conventions differ by category**: `WORKFLOW_TOOLS` (real-DFT prep/analysis
  pairs) write gnuplot `.dat`+`.gplot` pairs; `MLSIM_TOOLS` write matplotlib PNGs
  directly. Match the convention of the category/sibling tool you're extending, not
  personal preference.
- **Foundation-model comparison**: `--custom-model <path>` (+
  `--skip-foundation-comparison` to opt out) is the repeated pattern across
  `stb-mlphonons`/`stb-mlelastic`/`stb-mlff*`/`stb-mladsorb`/`stb-mleos` for
  overlaying a fine-tuned model's result against the raw MACE-MP-0 foundation model
  on the same structure/plot.
- **Consuming an already-finished external DFT result as if it were live**: when a
  tool needs to feed pre-computed energy/forces (read from a finished `calc.out`/
  `.FA`) into an ASE optimizer, wrap the `Atoms` in
  `ase.calculators.singlepoint.SinglePointCalculator(atoms, energy=E, forces=F)` and
  call the optimizer's `.step()` directly rather than `.run(fmax=..., steps=1)` —
  `run()`'s convergence check re-evaluates forces on the just-moved geometry, which
  a `SinglePointCalculator` cannot do and raises `PropertyNotImplementedError`.
  Convergence must instead be checked explicitly against the pre-step forces. See
  `neb_cycle.py` for the reference implementation (also guards against NaN
  tangents/forces from near-degenerate neighboring images via
  `np.errstate(invalid="raise")` before taking a step).
- **CLI-only tools**: a tool meant to be invoked from outside any live interactive
  session (e.g. from inside a cluster job-submission script's loop) should NOT be
  registered in `stb_suite.py`'s menu dicts — see `stb-nebCycle`.

### Tool inventory by category (non-obvious details only)

**2-Structures**: `stb-unitcell` (reduce to primitive/conventional/refined cell via
`core/symmetry.py`); `stb-crystalbuilder` (inverse — builds from space group +
Wyckoff sites via pymatgen); `stb-defect --all-inequivalent-sites`/`--ml-rank`
(auto-enumerate symmetrically distinct sites, optionally MACE-relax and rank each);
`stb-fetch` (network-dependent — COD, Materials Project, or any OPTIMADE-compliant
database; collapses same-element disorder automatically, rejects genuine
multi-element disorder); `stb-passivate`/`stb-slab --passivate` (H-terminate
dangling bonds via `core/passivation.py`); `stb-molecule` (ASE's G2 database,
case-sensitive names); `stb-mlrelax` (`--relax-cell` auto-masks vacuum-padded axes
via `build_cell_mask`); `stb-amorphize` (melt-quench via `NPTBerendsen`, bulk-only,
float32 for MD/float64 for the final static relax).

**3-Analysis**: `stb-aimdAnalysis` (RDF/MSD-diffusion/VDOS from a SIESTA `.ANI` or
any ASE-readable trajectory via `--trajectory`; velocities are finite-differenced
from positions, so VDOS is qualitative only); `stb-xrd` (pyxtal diffraction engine;
see `core/xrd.py` peak-list-vs-broadened-profile note above — this is a permanent
comparator rule, not a one-off fix).

**4-Workflow item 8, Adsorption** — 4 stages, in order:
`stb-adsorb` (Stage 1: enumerates sites, writes `sites/site_*/`; forces
`Slab.DipoleCorrection`+`DFTD3` unconditionally, `Spin` via `--force-spin` default
on; `--ml-rank` optionally MACE-screens candidates/orientations first via
`--n-orientations-polar`/`--n-orientations-azimuthal`) → `stb-adsorbBsse` (Stage 2:
reads each site's *relaxed* `siesta.XV` and writes ghost-fragment BSSE folders at
that geometry — evaluating BSSE at the pre-relaxation guess systematically
under-corrects) → `stb-adsorbAnalysis` (Stage 3: aggregates energies/BSSE/physical
diagnostics; `--compute-gibbs` is the one deliberate, intentional exception to "an
Analysis stage never generates new folders" — it writes vibrational-Hessian
displacement folders for Stage 4) → `stb-adsorbGibbs` (Stage 4: reads those
folders, diagonalizes the mass-weighted Hessian, reports ΔG(T); the isolated
-adsorbate reference's entropy is vibrational-only, no translational/rotational gas
-phase term). An "Analysis stage must never generate new geometry/folders" is the
suite's general workflow-placement rule (see also `stb-hubbardu`'s 3-stage split) —
`--compute-gibbs` is the single documented exception, not a precedent to extend
casually.

**4-Workflow item 9, NEB / Reaction Path** — `stb-neb` (Stage 1: takes
`--initial`/`--final`, two already-relaxed endpoint structures the user provides
directly — no symmetry-derived site enumeration; that was the previous Stage 1,
`stb-nebSites`, removed because restricting NEB endpoints to symmetry-derived
candidate sites defeated the point of NEB. Hard-validates the two endpoints share
the same composition and the same lattice — within `1e-3` Å, `require_lattice_match`
— exiting with a clean `[ERROR]` on either mismatch instead of guessing;
`--force-spin`/`--force-vdw`/`--force-dipole` (all default ON, `--no-force-*` to
opt out) control a `config_extra.fdf` written into every image folder, reusing the
same block constants as `core/adsorption_sites.py`'s `write_reference_folder`. Then
interpolates the band, `--mode 1-3`, default 2). A pure MACE-MP-0 path with no
SIESTA at all does not exist as a mode here — that capability lives exclusively in
`stb-mlneb` (5-ML Simulations), which already existed as "the pure-MACE extraction
of that same capability" (its own module docstring) before this workflow's own
former `--mode 1` (JSON-only, no SIESTA) was removed as a straight duplicate of it;
the numbering below was then closed up (old modes 2/3/4 → new 1/2/3) rather than
leaving a permanent gap. Every generated file/folder (the `neb_setup.txt` report,
`neb_path.xyz`, and the mode-specific output below) is written under a `neb_run/`
subfolder of `--output-dir`, not `--output-dir` directly — matching `stb-adsorb`'s
own `sites/`-nesting convention — so `stb-nebCycle`/`stb-nebAnalysis --dir` must
point at `<output-dir>/neb_run`, not `<output-dir>` itself:

| Mode | Path engine | Output (under `neb_run/`) | SIESTA's role |
|---|---|---|---|
| 1 | 100% MACE-MP-0, then one read | `image_NN/` (single-point) | one round, energy only |
| 2 (default) | MACE-MP-0 + real-DFT refinement cycles | `cycle_00/image_NN/` + loop | several rounds, real forces |
| 3 | 100% real-DFT NEB | `cycle_00/image_NN/` + loop | every round |

`--mode` replaces the older `--ml-neb` boolean outright. If `--initial`/`--final`
each carry a matching `neb_manifest.json` (no in-repo tool currently produces one,
but any future tool that builds two atom-index-corresponding structures may write
one — see `core/neb_manifest.py` below), atom correspondence is proven rather than
guessed and `--autosort-tol 0` is forced automatically; otherwise falls back to
distance-based matching, retrying once at `autosort_tol=0` if the initial distance
match fails (common for a small/densely-relaxed cell).
`check_endpoint_displacement()` warns (advisory only, never auto-enables anything)
when a broad fraction of the structure moved between endpoints and `--idpp` wasn't
used. → `stb-nebCycle` (CLI-only, not in the menu — see "CLI-only tools" above; one
call = one real-DFT NEB step for modes 2/3, using the `SinglePointCalculator`
+`.step()` pattern) → `stb-nebAnalysis` (Stage 2, mode-aware router reading a `#
MODE: N` marker in `neb_setup.txt`; for modes 2/3 always analyzes the
highest-numbered `cycle_NN/`, plots barrier-vs-cycle once 2+ complete cycles exist).

**4-Workflow item 15, GQCA** — `stb-gqca` (Stage 1: pair-cluster (N=2) reference
structures AA/AB/BB of a substitutional alloy on one sublattice, via icet's
Hart-Forcade enumeration; the AB structure is the 50:50 cell whose nearest-neighbour
pair correlation is closest to -1, which a frustrated (FCC/triangular) sublattice cannot
reach -- reported, never hidden; all 3 folders are written for a real CG relaxation,
`MD.VariableCell T` only for a 3D input) -> `stb-gqcaAnalysis` (Stage 2: per-site energies,
`Delta_j`, mass-action solve, `Omega = 2*Delta_1` and `T_c` estimate, level-of-theory and
convergence checks read back from each folder, tables/interpretation, optional dashboard
`--png/--show`, `--view`, `--export-relaxed`). Conventions worth knowing: the forced
directives (`MD.TypeOfRun CG`, `MD.NumCGsteps` > 0 -- SIESTA's default 0 would skip the
relaxation --, `MD.VariableCell T`, the per-cell `kgrid.MonkhorstPack`) go into each folder's
`config_extra.fdf`, `%include`d on top of an untouched `calc.fdf` template (the old append-at-the-end
form lost to a template's own `kgrid` because the fdf reader is first-occurrence-wins);
`--calc` defaults to `calc.fdf`; `--kdensity` defaults to 0.2 1/Ang like the other
workflows. **Normalization in `core/gqca_solver.py`**: the reference energies are per
sublattice site, the cluster gas is made of N-site clusters -- Boltzmann weights use the
cluster energy `N*Delta_j`, the cluster entropy is divided by `N`; the random limit must
reproduce the exact ideal entropy per site (a per-cluster entropy once made it exactly 2x
too large and shifted every result to 2T; `test/4-workflow/15-gqca/analysis/test.sh` guards
this with an ideal-entropy check and a brute-force minimization of `H - T*S`). The entropy's
closed form beyond the random limit is still flagged `[UNVERIFIED]`. The example
(`examples/4.15-gqca/`) ships REAL finished SIESTA runs in `reference_results/` (InGaAs/AlGaAs,
DZP/TZP, one mesh test) so Stage 2 runs genuine data; its README documents an unresolved
AlGaAs `Delta_1` anomaly (~29 meV left after ruling out k-points, cell shape, mesh, strain and
partly the basis) -- do not present AlGaAs numbers as validated.

**4-Workflow item 16, Optical Properties** — `stb-optical` (Stage 1: one folder per polarization
direction of an ALREADY-RELAXED structure -- `xx/yy/zz` diagonal, `xy/xz/yz` biaxial bisectors whose
reconstruction `eps_ij = eps_bisector - (eps_ii + eps_jj)/2` Stage 2 does automatically) ->
`stb-opticalAnalysis` (Stage 2: reads each folder's `.EPSIMG` + the `Optical:` block of `calc.out`,
Kramers-Kronig `eps1`, n/k/alpha/R/L/sigma1, and the diagnostics around them: SIESTA's f-sum check
and `N_eff(E)` recomputed from `eps2` (the two agree), static constants, onset/Tauc/Penn gap, Kohn-Sham
gap read from `<label>.EIG`, peaks and plasmon, anisotropy checked against the symmetry-allowed
tensor, `--scissor` what-if (rigid `eps2` shift), `--dimensionality-correction`, `--show/--png`
9-panel dashboard). Conventions worth knowing: every mandatory directive goes into each folder's
`config_extra.fdf` -- single-point guards (`MD.Steps 0`, `MD.NumCGsteps 0`, `MD.VariableCell false`), the
optical block, and the SCF k-grid, which is forced ONLY when the template has none (SIESTA would sample
Gamma) or on `--kgrid/--kgrid-density`, never over a template's own grid; **`Optical.PolarizationType
polarized` is mandatory** (the default `polycrystal` silently ignores `Optical.Vector`). The tensor
symmetry analysis (`core/optical_setup.symmetry_analysis`) projects onto the subspace invariant under the
space group in the cell's own Cartesian frame, so it works for any orientation; Stage 1 reports the minimal
direction set and Stage 2 checks finished spectra against it. Shared helpers: `core/fdf_directives.py`
(SIESTA's first-occurrence-wins reading of directives/k-grids, also used by `stb-gqcaAnalysis`),
`core/optical_setup.py`, the analysis functions in `core/optical_properties.py`, and
`siesta_log.get_optical_info/read_eig_gap`. Hard-won facts (documented in `examples/4.16-optical/`, whose
README is the page of this function on the site): the optical mesh (`Optical.Mesh`, independent of the SCF
grid) sets peak positions -- converge it first; broadening sets heights; once all bands are used the f-sum
is limited by the BASIS, not the band count; a metal's low-energy "onset" is a k-mesh artifact, not a gap;
the 2D out-of-plane correction is ill-conditioned and depends on the user's `--thickness`. The example
ships REAL finished SIESTA runs (Si with an 11-run convergence study, graphene). Pseudopotential caution:
a bank's quality is not guaranteed -- check the relaxed geometry against a known value first (the
PseudoDojo S gave a wrong Mo-S bond in one test; graphene uses the `virtual_vault` PSF carbon).

**4-Workflow item 17, MLFF** — shared across the 4 stages: one definition of an unphysical config,
`mlff_workflow.UNPHYSICAL_FMAX` (15 eV/Ang) / `UNPHYSICAL_DE` (1000 meV/atom) -- Stage 1's ML-screen flag,
Stage 2's SIESTA-label warning, Stage 4's default `--max-force` and `--auto-sampling` limit (before, 10 / 15 /
none: an auto-sampled 0.15 Ang graphene round gave 33-45 eV/Ang labels and a REGRESSED round); every report
is numbered with ALL warnings in a last section and third-party output (prints, `warnings`, and the root
`logging` logger -- mace's 'No dtype selected' -- via `core.cli.capture_library_noise`) only counted; every
stage has `--show` (Stage 1 adds coverage/mlff_coverage.png next to its gnuplot files) and `--save-report`
writes into the run folder (Stage 2: the round dir); a volume strain changes only the periodic axes in both
Stage 1 `--volume-strains` and Stage 4 `--sampling-method strain`; Stage 3's E(V) needs 5+ distinct test
volumes (Birch-Murnaghan has 4 parameters) and Stage 1 says whether its plan gives them. **Dispersion** (`--d3`, Stage 1; menu question too): `DFTD3 .true.` goes into every folder's `config_extra.fdf` (`mlff_workflow.DFTD3_BLOCK`; SIESTA's own DFTD3 defaults) and `reference.dftd3 = {on, forced, source}` into the JSON; every later folder follows it (`dftd3_block(state)`: Stage 4 picks, Stage 3 `--md-check` frames; Stage 4 also has `--d3` for runs without a JSON); a template's own `DFTD3 .true.` travels with the copied calc.fdf. Stages 2 and 3 print the labels' dispersion (`dispersion_check`, from the folders' own settings) and warn when it disagrees with the JSON; the round records `settings.dftd3`. A D3-labelled model is PBE+D3 already -- never add `dispersion=True`/torch-dftd on top (Stage 2's next steps say so). Device defaults
differ on purpose: Stage 2 auto-detects cuda (training), Stages 3-4 default to cpu like every MACE tool.
4 stages that talk through one `mlff_workflow.json` inside the
run folder `mlff_runs/` (Stage 1's default `-o`, like the other workflows' `*_runs`; it holds the
`mlff_config_*`, `mlff_test/` and `round_NN/` folders; the later stages default to
`mlff_runs/mlff_workflow.json`, or `./mlff_workflow.json` when run inside it — `default_path`) (`core/mlff_workflow.py`: reference .fdf + pseudo source, every
config with `set` train/test, `kind`, origin, proposing round; one entry per training round
with its committee — seed, model path, validation metrics —, settings, test-set scores and
active-learning picks; history of every call; paths relative to the file, atomic writes).
Every stage still accepts explicit paths and works without the file.
`stb-mlff` (Stage 1, Prep: rattled/AIMD-sampled training configs from the reference — the
suite's usual `-s structure.fdf -c calc.fdf` pair (defaults; calc.fdf must `%include` the
structure file, which is the one rewritten per folder, calc.fdf copied verbatim) or one
self-contained `--file`; the JSON records `reference.file` (structure) + `reference.calc`, and
Stage 4 takes the same inputs, adding `NumberOfAtoms` to a supercell's `config_extra.fdf` since
the calc.fdf one is not rewritten —, plus `--volume-strains` cells — isotropic, rattled — without which the training set
has one volume and the model cannot learn E(V)/B0/pressure; `--n-test N` held-out configs in
`mlff_test/test_NNN/`, the same kinds cycled with an independent seed; creates the JSON, an
old one is kept as `.bak`; every folder, Stage 4's too, gets a single-point `config_extra.fdf`,
so a relaxation/variable-cell template can't move the sampled geometry; the report reads every
written folder back and tabulates volume, displacement from the reference and shortest distance
per config, flags close contacts, near-duplicates (`--dup-threshold`, core/atomic_env) and test configs
outside the training range; `--ml-screen` runs a MACE-MP foundation model on every config
before SIESTA (max force, energy above the reference) to catch unphysical geometries; writes
gnuplot coverage plots to `<run folder>/coverage/` (`--no-plots`), stores each config's numbers
in the JSON, `--export-traj` writes all structures as one .xsf/.pdb/.xyz, and `--view` opens
them in ASE's viewer) → `stb-mlffAnalysis`
(Stage 2, Train: configs from the JSON (else `--path`/`mlff_config_*`), a `--committee N`
(default 3) of models with seeds `--seed`, `--seed+1`, ... — each its own train/valid split —
in `round_NN/m_s<seed>/` (`--committee 1` and no JSON: the old single model in `.`), recorded
as a new round. Splits train/valid in Python — not via `mace_run_train --valid_fraction`, which
is ignored once `--valid_file` is given — then fine-tunes without `--multiheads_finetuning`/
`--E0s foundation`, because SIESTA's absolute total energy is not on the foundation model's
energy scale (`average` and `multiheads_finetuning`/`foundation` are mutually exclusive in
`mace_run_train` itself). The recipe defaults come from a 3-seed bulk-Si benchmark against a
held-out SIESTA test set (see `DEFAULT_LR`/`DEFAULT_STRESS_WEIGHT`): E0s default to
`foundation_residual_e0s` (`--e0s foundation-residual`: lstsq of E_SIESTA minus the
foundation's own interaction energy, passed as an explicit `{Z: E0}`) because `average`
assumes zero interaction energy and starts the model eV/atom off; `--lr` 0.001 because
0.01 wrecks the pretrained model in epochs 0-1; stress is fitted only through `--loss
stress` — `--stress_key`/`--compute_stress` alone leave the energy+forces loss, which once
silently ignored it — and `--stress-weight` defaults to 1e4, since stress in eV/Ang^3 is a
~2 % term of the loss even at 100. Energy metrics vs the foundation model are offset-removed:
its absolute E RMSE is the ~1e5 meV/atom scale offset, not an error. Report: dataset per kind
(energy spread, forces, pressure, volume; SCF convergence read from each log), the recipe, a
per-member training table (learning history from mace_run_train's `results/*_train.txt`, kept
epoch from its log), every model evaluated on every config (train/valid errors, an out-of-sample
pool -- each config scored by the members that did not train on it -- per kind vs the
foundation), committee spread vs held-out error; ALL warnings are collected into a last
`[10] WARNINGS` section, library noise captured; `round_NN/mlff_training.png` (6 panels),
`--show` displays it, `--save-data` writes per-config errors. Also: a level-of-theory check
across folders (settings in SIESTA's reading order, k-sampling as k x a so a divided supercell
grid still matches, PAO.Basis block and pseudopotential md5); `--cv K` (K-fold: member k
validates on fold k, every config out-of-sample once; overrides --committee); a comparison table
against earlier rounds; the worst-predicted configs (`--view-worst N` writes them with per-atom
`F_error`, also as initial charges for ASE's viewer coloring, and opens them); an MD stability
test per model (NVT then NVE in an auto supercell of the reference, blow-up criteria, NVE drift;
`--no-stability`) that also times the model (ms/step, ns/day); `--init-from previous|PATH`
continues from the last round's plain `.model` -- mace_run_train accepts a model file as
`--foundation_model` -- while the comparison baseline stays the MACE-MP checkpoint) → `stb-mlffValidate`
(Stage 3, Validate: scores every round's committee mean and members on the test set — E
absolute and offset-removed, F, stress, per `kind` —, the committee calibration (true error /
spread, per-config and per-atom correlation), optionally phonons vs a finished
stb-phononsCreate folder (`--phonon-ref`: the committee-mean model on the SAME displacements,
in an Angstrom/eV phonopy object built from the SIESTA one, frequencies compared q-point by
q-point, imaginary modes counted); verdict STOP (`--target-force` met) / REGRESSED (latest round
worse than the best earlier one by more than `--plateau` -- it names the round to use) / PLATEAU (F RMSE
gain < `--plateau` %, default 5) / CONTINUE; scores stored per round, per config too; writes no folders except with --md-check (below).
Report in the Stage-2 style: numbered sections, ALL warnings in a last `[10] WARNINGS`, library noise
captured; every round's errors with a bootstrap +/- over test configs (how small a difference the test set
can resolve) and F relative to the RMS SIESTA force; the test set compared with the detailed round's
training configs (level of theory via Stage 2's `check_theory`, which test configs fall outside the
training V/F/E range); one round in detail (`--detail-round`, default latest): per kind vs foundation, per
config (worst atom, spread, novelty, range), error vs force/novelty correlations, test vs Stage-2
out-of-sample error (generalization gap), member table with outliers and the averaging gain; test folders
missing from the JSON are flagged; without a JSON the kind comes from the folder name (`large_01` ->
`large`); equilibrium configs (spread < 1 meV/Ang) are left out of error/spread ratios. Also: force
errors per element; a Birch-Murnaghan E(V) of every test kind spanning >= 4 volumes (V0/B0 of each
round and the foundation vs SIESTA -- on the Si benchmark B0 94.4 SIESTA, +3.4 % fine-tuned, -21 %
foundation); a power-law learning curve F = a N^-b across rounds (N = `n_usable`, needs 3+ sizes;
configs needed for `--target-force`); `--compare A B` config by config and per kind. `--md-check` is
the ONE option that writes folders: an NVT run with the best round's first model (optionally in a
`--md-check-supercell`, k-grid divided via stb-mlffActiveLearning's `supercell_kgrid_block`), frames
written with `stb.mlff.write_config_folder` to `<run>/md_check/rNN_<T>K/frame_NNN/` and recorded in the
JSON as their own set `md_check` (never the test set); after SIESTA, every later run reports the
error along the trajectory (growth first->last third), `--md-check-to-train` moves the finished frames
into training. Sections [0]-[13], warnings last. `--show`, `--save-data` (mlff_validation_errors.dat),
`--view-worst N` (mlff_test_worst.xyz), mlff_validation.png: 6 panels plus a third row with whichever
of learning curve / E(V) / MD check / comparison is available) → `stb-mlffActiveLearning` (Stage 4:
committee, reference, pseudo source, training set and output dir from the JSON (latest round,
or `--round K`) unless given; picks recorded as training configs of the proposing round.
Screens candidates — `--sampling-method rattle|md`, `--md-temperature` one NVT trajectory per
value, optional `--supercell NA NB NC` of the reference, written with
`structure_io.replicate_fdf_cell` (ase `repeat()` atom order) and the reference k-grid
divided by the supercell in `config_extra.fdf`, since a small training cell only constrains
what its own periodicity allows (8-atom Si: Gamma and X phonons, not L; a 2x2x2 MD round cut
the Si phonon RMSE 0.49 -> 0.31 THz). Score: with 2+ models the max per-atom force spread
across the committee; with one model the largest predicted force. Benchmark caveat: seeds
fine-tuned from the same foundation model are nearly the same model — on Si the spread tracked
the error per config (corr 0.81) but not per atom (~0) and was 3-15x smaller, so it ranked
candidates about like max|F| does. Greedy pick by score, skipping candidates over
`--max-force` and those whose atomic environments are all within `--min-novelty` of the
training set or of an earlier pick; novelty uses a model-independent fingerprint (sorted
distances to the 16 nearest neighbours) in units of the training set's own typical spacing,
because `MACECalculator.get_descriptors` fails on the compiled model; max-over-atoms novelty is
weak in a supercell — consecutive MD frames of a 64-atom cell all pass it. Report in the Stage-2/3 style ([0]-[8], ALL warnings last, library noise captured, progress while screening): the committee's quality from the JSON (4.17.2 out-of-sample, 4.17.3 test F and verdict) and Validate's error/spread calibration, which turns each spread into an estimated error -- only when the models ARE that round's committee, never for an explicit --model; the labelled training range (SIESTA max|F|, shortest distance, wall time); the pool per sampling setting (max|F|, spread, novelty vs the training set alone, d_min), MD runs that stopped early, spread-vs-max|F| correlation (> 0.9: the spread only ranks by distortion), the share of the pool that is new territory; the selection funnel; each pick (estimated error, novelty, d_min, pool rank, worst atom, beyond-training flags, close contacts), the picks' mutual diversity; a SIESTA cost estimate (median training wall time scaled as N^3). `--dry-run` writes no folder and leaves the JSON alone; `--show`, `--save-data` (mlff_al_<tag>_pool.dat), `--view-picks` (mlff_al_<tag>_picks.xyz, per-atom spread as initial charges), mlff_al_<tag>.png (6 panels), tag rNN_<run>; the JSON's active_learning entry also stores the pool statistics and skip counts; picks are wrapped into the cell before writing. Sampling `strain` (random symmetric strain up to `--strain-max`, periodic axes only, rattled; the pick's own cell goes into `%block LatticeVectors` via `write_config_folder(raw_lattice=cell / template LatticeConstant)`) and `--md-ensemble npt` (bulk only; refused with a vacuum axis). `--score force|energy|both`: the energy spread is taken after removing each member's constant offset over the pool (`center_energy_spread` -- members have their own E0s, the raw spread was ~100 meV/atom of pure offset). `--selection diverse` (default): score x min(1, novelty/`--diversity`) against training + picks so far, reported against the greedy pick (overlap, score ratio, closest pair). `--auto-sampling` scans `--scan-values` (default stdev 0.03-0.3 Ang / 300-1500 K / strain 1-8 %) with `--scan-n` candidates each; physical settings (<= 10 % close contacts, stable MD, --max-force, median energy <= --max-energy) are then limited to those whose median predicted energy above the reference stays within `--frontier` (default 1.5) x the training set's 90th-percentile energy, and the highest median score among them wins. **Energy cap**: every pick, auto-sampled or not, is at most min(`--max-energy` (default UNPHYSICAL_DE), `--frontier` x that 90th percentile) above the reference (committee-mean energies, the reference cell evaluated once); the report gives energies also as a temperature (3/2 kT per atom). Why: the committee spread grows with distortion, so 'highest score that is still physical' always chose the edge of the physical range -- on 2-atom graphene a 0.15 Ang rattle, picks 1.2-2.4 eV/atom above equilibrium (~10^4 K), and two such rounds made the test F RMSE 21 -> 37 meV/Ang and the phonons 1.6 -> 2.4 THz; with the cap the same round picks 0.28-0.55 eV/atom. Stage 3 lists training configs more than one test-set energy spread above the highest test energy ('Far from the test', `far_training_configs`), and a learning curve whose error GROWS with N is reported as such (`F = a N^+b`), pointing there. Level-of-theory check: k x a is measured in reference cells (`_kxa_in_reference_cell`: per axis round(l / l_min) reference cells x k x l_min), so strained cells with the reference's own grid, as Stage 1/4 write them on purpose, no longer warn. `--target-force` gives a verdict STOP / SMALL (only a tail above target) / CONTINUE / UNSURE (no calibration: the spread is only a lower bound). Seeds: the effective seed is `--seed + 1000 x` the round the picks will train (recorded as `seed` in the active_learning entry), so every round screens a fresh pool, and rattle candidates add their own +30000 (`RATTLE_SEED_OFFSET`) -- Stage 1 rattles config i with seed + i through the same generator and ASE draws stdev x one fixed normal vector per seed, so with the shared default 42 candidate i used to be training config i along the very same direction, 3x larger (cos = 1.00; now ~0.1). An ACTIVE-LEARNING HISTORY table lists every earlier run of the workflow (pool score, calibrated max error, new territory). Sections [0]-[10]; the figure gains a third row (scan, energy-vs-force spread, history) when available).

**4-Workflow item 18, EOS** — `stb-eosInputs`/`stb-eosAnalysis`: real-SIESTA analog
of `stb-mleos`, same isotropic volume-scan convention (`factor = (1 +
strain_pct/100)**(1/3)`); both share `core/eos_fit.py`.

**5-ML Simulations** (category order): `stb-mlmd` (NVE/NVT/NPT; all physical units
explicitly converted to ASE's internal unit system before use — `* ase.units.fs`,
`* ase.units.GPa`, `/ ase.units.bar` for `--compressibility`, given in 1/bar; tracks
`E_total`, not just `E_pot`, for NVE drift). The compressibility was once passed raw and
labelled eV/Ang^3 -- ASE wants 1/(eV/Ang^3), so water's 4.57e-5 1/bar was ~1.6e6 x too small
and every NPT cell (stb-mlmd, stb-amorphize, stb-mlffActiveLearning --md-ensemble npt) stayed
put (1e-6 volume change at 10 GPa); converted, an 8-atom Si cell shrinks ~10 % at 10 GPa, as
B0 ~95 GPa implies. `build_dynamics(..., compressibility_bar=)` does the conversion; stb-amorphize
calls NPTBerendsen itself and converts there. `stb-mlphonons` (standalone phonon calc; ASR correction via
`symmetrize_force_constants`; imaginary-mode check uses a small negative tolerance,
`_IMAGINARY_MODE_TOL_THZ`, not a strict `< 0`; `--freeze-unstable-mode`/
`--animate-mode` displacement amplitude must be calibrated via a minimum-image
-convention helper, not raw Cartesian differencing, or PBC wraparound corrupts the
measured amplitude). `stb-mlelastic` (stress-strain fit via MACE, reuses
`elastic_analysis.py`'s post-fit numerics). `stb-mlsearch` (`BasinHopping` or
simulated-annealing MD cooling; positions only, cell shape untouched during the
search — not crystal structure prediction). `stb-mlmelting` (Lindemann-index
bracketing; a small defect-free periodic cell has no free surface to nucleate
melting from, so this method reliably SUPERHEATS past the true melting point —
document this limitation with any estimate). `stb-mlconvergence` (same reference
calc run at multiple MACE-MP-0 sizes; checks whether the *answer* has converged
with model size, not against real DFT). `stb-mlneb`/`stb-mldiffusion` (pure-MACE
NEB and vacancy-migration-barrier screening; `find_neighbor_shell` uses each site's
own local nearest-neighbor distance, not a global minimum). `stb-mlgcmc` (canonical/
grand-canonical MC, single-atom adsorbate only, rigid host; `mu` is extremely
sensitive — `beta` is ~O(40) eV⁻¹ at room temperature, so print the ideal-gas
reference loading up front to catch a bad `--mu` before an expensive run).
`stb-mladsorb` (pure-MACE adsorption-site screening; `--n-orientations` samples
rotations before a full relax, since `AdsorbateSiteFinder.add_adsorbate` always
re-aligns the molecule's z-axis to the surface normal regardless of input
orientation). `stb-mleos` (E(V) curve fit, distinct method from `stb-mlelastic`'s
stress-strain bulk modulus — comparing the two is a genuine cross-check, not a
duplicate feature).

**6-Utils**: `stb-status` (per-folder run-type/SCF/force/energy summary, `--path`
glob support, every field degrades to `None` rather than raising); `stb-archive`
(packages a finished calc into `.tar.gz` + `MANIFEST.txt`; complements `stb-clean`).

## Domain conventions worth knowing

- Structure file formats handled throughout: POSCAR (VASP), CIF, FDF (SIESTA), XYZ
  (+ separate `.lattice` file), XSF, FHI-aims `geometry.in`, DFTB+ `gen`.
- SIESTA output files consumed by analysis tools: `.bands`, `PDOS.xml`, `.RHO`/`.VT`
  grid files, `.XV`, etc.
- Core dependencies: `numpy`, `ase`, `matplotlib`, `pymatgen`, `spglib`, `sisl`,
  `pybader`, `scipy`, `phonopy`. `ase` in particular is often imported inside a
  `try/except ImportError` since it's only required for CIF I/O in some tools.
