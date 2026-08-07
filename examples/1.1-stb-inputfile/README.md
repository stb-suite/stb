# 1.1 — Input File Generator (`stb-inputfile`)

## What this tool does

You almost always start a SIESTA project with a *structure* — a set of atoms
in a box, with no calculation settings attached (that's what `structure.fdf`
in this folder is). `stb-inputfile` takes that bare structure and turns it
into a complete, ready-to-run `calc.fdf`: it works out a basis set, computes
a Monkhorst-Pack k-point grid from the cell's own density, and drops in the
right SIESTA blocks for whatever kind of run you asked for. Instead of
hand-assembling a `.fdf` file from scratch (and re-deriving the k-grid by
hand every time your cell changes), you answer one question — "what am I
calculating?" — and get a working input file back.

Before writing `calc.fdf`, it also always runs a quick structure validation
pass — symmetry (space group / crystal system) plus a handful of
malformation checks (atoms too close together, a left-handed cell, an
implausible density, a stale atom-count header) — and reports whatever it
finds as `[WARNING]` lines. This never blocks generation; it's a heads-up,
not a gate.

## When you'd reach for it

Any time you have a structure file (freshly built with `stb-slab`,
`stb-supercell`, `stb-crystalbuilder`, fetched with `stb-fetch`, converted
with `stb-translate`, ...) and need to turn it into an actual SIESTA
calculation. It's normally the very last step of structure preparation and
the very first step of running SIESTA.

## Two ways to run it

**A — direct CLI**, the fastest once you know the flags:

```bash
stb-inputfile structure.fdf -t relax -p dojo
```

**B — interactive `stb-suite` menu**, which asks the same questions as
guided prompts instead of flags:

```bash
stb-suite
# at the main prompt, type: 1.1
```

`1.1` jumps straight past the category menus (`1` = Inputs, item `1` =
Input File Generator) to this tool. It then asks you, one at a time: which
structure file, which calculation mode (numbered list), whether to enable
DFT-D3, whether to enable spin polarization, an optional pseudopotential
source, whether to save a text report, and whether to open the interactive
3D viewer — building the exact same `stb-inputfile ... -t ... -d3 -s -p
... --view` command underneath and running it for you. Good for exploring
what a tool can do before you've memorized its flags; the CLI is faster
once you have.

Both paths run the identical underlying tool and produce byte-identical
output — `example_1.1.sh` demonstrates this directly, running the same
`relax` calculation both ways and diffing the results.

## Files in this folder

- `structure.fdf` — bulk-silicon example structure (conventional cubic
  diamond cell, 8 atoms, 3D/fully periodic).
- `structure_molecule.fdf` / `structure_chain.fdf` / `structure_graphene.fdf` —
  an isolated CH4 molecule (0D), an infinite carbon chain (1D), and a
  graphene monolayer (2D).
- `example_1.1.sh` — the guided walkthrough (**not** an automated test —
  see `test/1-inputs/1-input_file/test.sh` for that). Pauses between
  sections so you can read before moving on; safe to re-run.
- `output/` — created by `example_1.1.sh` when you run it (git-ignored, not
  checked in). See below.

## Running the walkthrough

```bash
./example_1.1.sh
```

`stb-inputfile` always needs pseudopotentials to be handed off to SIESTA
for real, so every case below is generated with `-p dojo` and written into
its **own folder** under `output/` — each one already self-contained
(structure file + `calc.fdf` + pseudopotentials, and `kpath_bs.fdf` for
`bands/`) and ready to run in SIESTA as-is, no extra assembly needed:

| Folder                | Command (conceptually)                              |
|------------------------|------------------------------------------------------|
| `output/relax/`        | `stb-inputfile structure.fdf -t relax -p dojo`       |
| `output/relax_d3/`     | `stb-inputfile structure.fdf -t relax -d3 -p dojo` — same mode, DFT-D3 switched on |
| `output/relax_spin/`   | `stb-inputfile structure.fdf -t relax -s -p dojo` — same mode, spin-polarized |
| `output/bands/`        | `stb-inputfile structure.fdf -t bands -p dojo`, then `stb-kpath -f structure.fdf` to generate the `kpath_bs.fdf` that `bands` mode references |
| `output/molecule/`     | `stb-inputfile structure_molecule.fdf -t total_energy -p dojo` |
| `output/chain/`        | `stb-inputfile structure_chain.fdf -t total_energy -p dojo`    |
| `output/graphene/`     | `stb-inputfile structure_graphene.fdf -t total_energy -p dojo` |
| `output/silicon/`      | `stb-inputfile structure.fdf -t total_energy -p dojo`          |

The script also drives the exact same `relax` case through the interactive
`stb-suite` menu (non-interactively, via a piped `printf`, no folder kept
for it) and diffs the result against `output/relax/calc.fdf` — proving the
CLI and the menu produce identical output.

### Structure validation (always on, every case above)

Every run prints a `[2] STRUCTURE VALIDATION` section, right after parsing
the structure and before the k-grid is even computed. It always runs --
there's no flag to turn it on, and it never blocks generation. Every check
is listed **explicitly**, one row per check, each marked
`[OK]`/`[WARNING]`/`[SKIPPED]` (via the shared `core/structure_checks.py`
module, also used by `stb-fetch`/`stb-mlrelax`) -- not just a summary line
-- so it's always clear exactly what was verified:

```
Check                  | Result                                            | Status
--------------------------------------------------------------------------------------
Atom proximity         | min. distance = 2.211 Ang (>= 0.5 Ang threshold) | [OK]
Lattice handedness     | right-handed (positive determinant)               | [OK]
Atomic density         | 0.0500 atoms/Ang^3 (within [0.01, 0.15])          | [OK]
NumberOfAtoms header   | declared 8, matches actual atom count             | [OK]
NumberOfSpecies header | declared 1, matches actual species count          | [OK]
```

Five specific things get checked:

1. **Atoms too close together** -- any two atoms closer than 0.5 Ang
   (minimum-image distance). Usually means two overlapping/duplicate atoms
   from a bad edit or a bad supercell/merge operation. `[SKIPPED]` if the
   structure has fewer than 2 atoms.
2. **Left-handed cell** -- the lattice vectors' determinant is negative
   (negative cell volume). A silent geometry bug: SIESTA will still run,
   but stress/chirality-sensitive downstream analysis can pick up sign
   errors from it.
3. **Implausible atomic density** -- atoms-per-Ang^3 outside roughly
   `[0.01, 0.15]`, a deliberately generous band covering everything from
   light/loosely-packed solids to dense covalent ones. The classic trigger
   is a `LatticeConstant` given in the wrong unit (e.g. Bohr instead of
   Angstrom) -- the cell ends up ~6.7x too large or too small in volume.
   **Only checked for a genuine 3D bulk structure** (see the vacuum note
   below) -- a vacuum-padded cell's volume is dominated by the artificial
   vacuum by design, so this metric would be meaningless (and would
   false-positive) there; `[SKIPPED]` instead, explicitly, rather than
   silently absent.
4. **Stale atom-count header** -- the `.fdf` file's own declared
   `NumberOfAtoms`/`NumberOfSpecies` line doesn't match what's actually in
   `%block AtomicCoordinatesAndAtomicSpecies`/`%block
   ChemicalSpeciesLabel`. Catches a header left un-updated after hand-editing
   the structure. `[SKIPPED]` if a given header line isn't present at all
   (they're not required) -- two independent rows, since either header can
   be present/absent on its own.
5. **Symmetry** -- space group and crystal system, always printed
   (informational, not a warning by itself). **This is where dimensionality
   matters**: an ordinary space group only means what it says for a fully
   3D-periodic structure. `molecule/`, `chain/`, and `graphene/` below are
   each vacuum-padded along one or more axes, so `stb-inputfile` still
   reports *a* space group (spglib will happily analyze the padded cell as
   if the vacuum were just another periodic direction), but appends this
   extra `[WARNING]`:

   ```
   [WARNING] Structure has a vacuum-padded axis -- the space group above
   treats it as an ordinary periodic direction and may not reflect the
   true symmetry. Use stb-symmetry (code 3.5) for a dimension-aware
   layer-group/point-group analysis.
   ```

   In other words: for `silicon/` (fully 3D) the printed space group is the
   real, physically meaningful answer. For `molecule/`/`chain/`/`graphene/`
   (0D/1D/2D) it's only a rough approximation of a real 3D crystal that
   doesn't actually exist -- the vacuum thickness itself can even change
   which "space group" comes out, since spglib has nothing else to base a
   3D classification on. If you need the *correct* symmetry for a
   non-3D structure (a proper 2D layer group, or a 0D point group),
   run `stb-symmetry` (code 3.5) directly instead of trusting this
   quick check.

Run `stb-inputfile structure.fdf -t total_energy -p dojo` yourself and look
for `[2] STRUCTURE VALIDATION` in the output to see all of this live.

### Visualizing the structure (`--view`)

`--view` opens the structure in ASE's interactive 3D viewer right before
the tool exits (needs a local display, or `ssh -X`/`-Y`). Off by default --
not used in this script (an unattended walkthrough shouldn't pop a GUI
window), but try it yourself:

```bash
stb-inputfile structure.fdf -t relax -p dojo --view
```

### Why `bands/` needs an extra command

`stb-inputfile`'s `bands` mode writes a `calc.fdf` that references
`%include kpath_bs.fdf`, but doesn't generate that file itself — that's
`stb-kpath`'s (code 1.4) job. The script calls it too, so `output/bands/`
ends up genuinely complete instead of missing a file SIESTA would need.

### Why dimensionality matters

A k-point grid only makes physical sense along a genuinely periodic
direction. `stb-inputfile` measures the largest empty gap along each
lattice vector and, above a threshold (10 Ang), treats that axis as vacuum
instead of periodic when computing the k-grid — automatically, from the
geometry alone, no flag needed:

| Folder       | Structure                             | Dimensionality | Periodic axes | Typical k-grid shape  |
|---------------|-----------------------------------------|:---:|:---:|-------------------------|
| `molecule/`   | `structure_molecule.fdf` (CH4)          | 0D | 0 | `1  1  1`   (isolated)   |
| `chain/`      | `structure_chain.fdf` (carbon chain)    | 1D | 1 | `N  1  1`   (a wire)     |
| `graphene/`   | `structure_graphene.fdf` (graphene)     | 2D | 2 | `N  N  1`   (a slab)     |
| `silicon/`    | `structure.fdf` (bulk silicon)          | 3D | 3 | `N  N  N`   (a crystal)  |

This is why the same `stb-inputfile <structure> -t total_energy` command
works unchanged whether the structure is a molecule, a wire, a slab, or a
bulk crystal — the vacuum gap in the structure file itself is the signal,
so you never need a separate "molecule mode" or "slab mode" flag.

## Try it yourself

Once you've read through the walkthrough, experiment directly:

```bash
stb-inputfile structure.fdf -t total_energy -p dojo
stb-inputfile structure.fdf -t aimd -p dojo
stb-inputfile structure.fdf -t bands -d3 -s -p dojo   # -d3 and -s combine with any mode
```

Open each generated `calc.fdf` and compare it against `structure.fdf` and
against the other modes — that's the fastest way to see exactly what each
mode adds.

## Flag reference

```
stb-inputfile <structure_file> -t <mode> [-d3] [-s] [-p <pseudo-source>] [--view]
```

- `-t/--type` (required): `total_energy`, `relax`, `aimd`, or `bands`.
- `-d3/--d3` (optional): enable the DFT-D3 (Grimme) dispersion correction,
  independent of and combinable with any mode above.
- `-s/--spin-polarized` (optional): enable spin polarization, also
  independent of and combinable with any mode above. Defaults to
  non-polarized when omitted.
- `-p/--pp-path` (optional): a bundled bank name (`dojo`, ...) or a path to a
  folder containing `.psml`/`.psf` pseudopotential files.
- `--view` (optional): open the structure in ASE's interactive 3D viewer
  after generation. Off by default; needs a display.

Structure validation (symmetry + malformation warnings, see below) always
runs before generation -- there's no flag for it, it's always on.

Run `stb-inputfile --help` for the full list of options.

## What's next

Any `output/<case>/` folder is ready to hand to SIESTA as-is. After the run
finishes, `stb-status` gives you a quick summary of what happened
(converged? final energy? which output files were produced), and the
`Analysis` category tools (`stb-bands`, `stb-dos`, ...) pick up from there.
