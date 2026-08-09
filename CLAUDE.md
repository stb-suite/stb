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

- `structure_io.py` — the only `.fdf` reader/writer (`read_fdf`, `write_fdf`,
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
interpolates the band, `--mode 1-4`, default 3):

| Mode | Path engine | Output | SIESTA's role |
|---|---|---|---|
| 1 | 100% MACE-MP-0 (climbing-image NEB) | `neb_mace_result.json` | none |
| 2 | 100% MACE-MP-0, then one read | `image_NN/` (single-point) | one round, energy only |
| 3 (default) | MACE-MP-0 + real-DFT refinement cycles | `cycle_00/image_NN/` + loop | several rounds, real forces |
| 4 | 100% real-DFT NEB | `cycle_00/image_NN/` + loop | every round |

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
call = one real-DFT NEB step for modes 3/4, using the `SinglePointCalculator`
+`.step()` pattern) → `stb-nebAnalysis` (Stage 2, mode-aware router reading a `#
MODE: N` marker in `neb_setup.txt`; for modes 3/4 always analyzes the
highest-numbered `cycle_NN/`, plots barrier-vs-cycle once 2+ complete cycles exist).

**4-Workflow item 17, MLFF** — `stb-mlff` (Stage 1: generates rattled/AIMD-sampled
training configs from a reference `calc.fdf`) → `stb-mlffAnalysis` (Stage 2: splits
train/valid in Python — not via `mace_run_train --valid_fraction`, which is ignored
once `--valid_file` is given — then fine-tunes with `--E0s average`, not
`--multiheads_finetuning`/`foundation`, because SIESTA's absolute total energy is
not on the foundation model's energy scale; `average` and
`multiheads_finetuning`/`foundation` are mutually exclusive in `mace_run_train`
itself) → `stb-mlffActiveLearning` (Stage 3: screens new candidates with the
fine-tuned model, `--sampling-method rattle|md`, keeps the `--top-k` highest
predicted-force configs for real labeling).

**4-Workflow item 18, EOS** — `stb-eosInputs`/`stb-eosAnalysis`: real-SIESTA analog
of `stb-mleos`, same isotropic volume-scan convention (`factor = (1 +
strain_pct/100)**(1/3)`); both share `core/eos_fit.py`.

**5-ML Simulations** (category order): `stb-mlmd` (NVE/NVT/NPT; all physical units
explicitly converted to ASE's internal unit system before use — `* ase.units.fs`,
`* ase.units.GPa`, not passed through raw; tracks `E_total`, not just `E_pot`, for
NVE drift). `stb-mlphonons` (standalone phonon calc; ASR correction via
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
