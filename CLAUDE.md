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
- `test/` — example input files and manual smoke scripts organized by category
  (`1-inputs`, `2-analysis`, `3-workflow`, `4-utils`), each mirroring one or more
  tools. `3-workflow/<property>/` has `prep/` + `analysis/` subfolders for the 4
  paired prep+analysis tools (strain, elastic, cohesive, phonons). This is **not**
  an automated test suite (no pytest/unittest, no CI config exists).
  `test/test_translate/test.sh` is representative: it generates sample structure
  files in every supported format, runs the built CLI against them, and greps the
  output for expected content.

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
`bash test/test_translate/test.sh`) and inspect the produced output files, following
the same "generate sample input → run the CLI → check the output file" pattern used
in `test.sh`.

## Architecture

**Each module in `stb-suite/src/stb/` is a thin CLI script** (argparse `main()`,
registered as a console script in `stb-suite/pyproject.toml` under
`[project.scripts]`, e.g. `stb-kgrid = "stb.kgrid:main"`) that delegates shared logic
to `stb-suite/src/stb/core/`:

- `core/structure_io.py` — the only `.fdf` reader/writer (`read_fdf`, `write_fdf`,
  `to_pymatgen`, `rewrite_fdf_lattice`, plus thin accessors like `lattice_only`,
  `species_list`/`species_dict`). Use `raw_lattice_vectors()` instead of
  `lattice_only()`/`.lattice` when you're about to write back into an existing file's
  `%block LatticeVectors` while leaving its `LatticeConstant` line untouched (see the
  module docstring) — using the wrong one double-applies the lattice constant.
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

New format-specific structure readers/writers (POSCAR, CIF, XYZ, XSF, FHI, DFTB) still
live only in `translate.py` — it's the sole consumer of those formats, so there's
nothing to extract yet. If a second tool needs one of them, move it into
`core/structure_io.py` first rather than copying `translate.py`'s implementation.

When adding a new tool, mirror this pattern: a new standalone `stb/<name>.py` that
imports what it needs from `stb.core`, plus a new entry in `[project.scripts]` in
`pyproject.toml` (and update `test/` with a fixture if practical).

Newer `1-inputs` tools worth knowing about specifically:
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
shows a menu with 5 categories, in this order: **1 Inputs, 2 Structures, 3 Analysis,
4 Workflow, 5 Utils** — backed by five dicts, `INPUT_TOOLS`, `STRUCTURE_TOOLS`,
`ANALYSIS_TOOLS`, `WORKFLOW_TOOLS`, `UTILITY_TOOLS`, each keyed by menu number to
`{'title', 'description', 'func'}`. `INPUT_TOOLS` (category 1) holds only the 3
tools that configure an actual SIESTA run (input file, k-grid, k-path);
`STRUCTURE_TOOLS` (category 2, e.g. `stb-slab`, `stb-supercell`, `stb-defect`,
`stb-crystalcast`) holds everything that builds/generates/transforms a structure
file, split out once `INPUT_TOOLS` grew to 17 entries mixing the two concerns.
`WORKFLOW_TOOLS` is one level deeper (the 4 paired prep+analysis properties —
strain, elastic, cohesive, phonons — each with a `'stages'` dict of 2 entries
instead of a `'func'`); `run_sub_menu()` recurses into `entry['stages']`
automatically whenever it finds one instead of a `'func'`, so no separate menu
function was needed for the extra level. `_flatten_tool_codes()` builds a flat
`{"1.1": func, ..., "2.1": func, ..., "4.1.2": func, ..., "5.4": func}` lookup
from these same dicts at import time (`TOOL_CODES`), letting the main menu's
prompt accept a dotted code (e.g. `4.1.2`) to jump straight to a tool instead of
navigating level by level — regenerate nothing by hand here, it derives entirely
from the 5 dicts above. Every `run_*` function builds an `args: List[str]`
from interactive prompts (`get_input`/`get_float_input`/`get_int_input`) and dispatches
via `run_tool(tool_name, args)`, which shells out to the **installed** console command
by name (`subprocess.run`, must be on `PATH`) and centralizes error handling plus the
"Press Enter to continue" pause. Don't reintroduce a second dispatch path (resolving a
sibling script file via `__file__` and invoking it with `sys.executable`) — that used to
coexist for about a third of the tools for no documented reason and was consolidated
onto `run_tool()` alone.

## Domain conventions worth knowing

- Structure file formats handled throughout: POSCAR (VASP), CIF, FDF (SIESTA), XYZ
  (+ separate `.lattice` file), XSF, FHI-aims `geometry.in`, DFTB+ `gen`.
- SIESTA output files consumed by analysis tools: `.bands`, `PDOS.xml`, `.RHO`/`.VT`
  grid files, `.XV`, etc.
- Core dependencies: `numpy`, `ase`, `matplotlib`, `pymatgen`, `spglib`, `sisl`,
  `pybader`, `scipy`, `phonopy`. `ase` in particular is often imported inside a
  `try/except ImportError` since it's only required for CIF I/O in some tools.
