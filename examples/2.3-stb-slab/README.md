# 2.3 — Slab Builder (`stb-slab`)

## What this tool does

Cuts a Miller-index surface out of a bulk structure and adds vacuum along
the surface normal, producing a 2D-periodic slab ready for a surface-science
SIESTA calculation. Same **Structures** category (2) as `2.1-stb-2Dstacking/`
and `2.2-stb-supercell/` — everything here builds, generates, or transforms
a structure file, as opposed to **Inputs** (1), which configures an actual
SIESTA run.

## Why this matters (a bit of theory)

### Miller indices, terminations, and why more than one slab can come out

A Miller index `(h,k,l)` defines a family of lattice planes, but **cutting**
along that family is not unique: depending on exactly where between two
planes you cut, you can expose different atoms at the surface — different
**terminations**. `stb-slab` (via `pymatgen`'s `SlabGenerator`) enumerates
every symmetrically distinct termination for the requested `(h,k,l)`, sorts
non-polar/symmetric ones first (the usual physically preferred choice), and
by default keeps only the top of that sorted list (`--all` keeps every one,
`-i` lets you pick, `--termination N` picks by index directly).

- **Symmetric** — the same termination on both faces of the slab (related by
  a mirror/inversion). Preferred: an asymmetric slab has two physically
  different surfaces, which usually also means a spurious residual dipole
  across the artificial vacuum gap in a periodic calculation.
- **Polar** — see the dedicated caveat below; `stb-slab`'s own `Polar`
  column needs a caveat before you trust it.

### A real gotcha, verified live: the `Polar` column needs oxidation states it never gets

`pymatgen`'s `Slab.is_polar()` computes a dipole from each site's formal
**oxidation state** — its own docstring says so directly: *"the Slab must be
oxidation state decorated for this to work properly, otherwise the Slab
will always have a dipole moment of 0."* A structure read from a SIESTA
`.fdf` file never carries oxidation states (`.fdf` has no such concept) —
so `stb-slab`'s `Polar` column reads **`No` unconditionally**, regardless of
the real electrostatics of the cut.

Verified directly: cutting the same bulk NaCl both along `(1,0,0)`
(genuinely non-polar — each atomic layer already has equal Na/Cl) and along
`(1,1,1)` (the textbook **polar** rocksalt surface — alternating pure-Na and
pure-Cl layers, Tasker's "Type III" classification) both report `Polar: No`
through `stb-slab`. Decorating the exact same two structures with formal
`Na+1`/`Cl-1` oxidation states in a Python REPL and calling `is_polar()`
directly recovers the real physics: dipole per unit cell `0.0` for `(100)`
vs. `26.05` (arbitrary units, e/Ang) for `(111)` — see `polar-caveat/`
below for the live, side-by-side proof.

**Takeaway**: use `stb-slab`'s `Polar` column only as "pymatgen found no
reason to flag this" — for an ionic/mixed-valence material, cross-check a
surface's true polarity yourself (electronegativity difference + which
species terminates each side) before trusting it blindly.

### Dangling bonds and `--passivate`

Cutting a covalent structure (e.g. tetrahedral Si) necessarily leaves atoms
at the new surfaces with fewer bonds than in the bulk — dangling bonds,
each roughly a half-filled orbital, a strong artificial perturbation on the
electronic structure. `--passivate` caps single-missing-bond sites with a
passivating atom (H by default) placed along the geometrically missing-bond
direction (`core.passivation.py`, shared with `stb-passivate`) — see
`passivate/` below for a real, fully-resolved 8-dangling-bond case.

### Surface relaxation and `--ml-relax`

A freshly cut slab's surface atoms sit at their *bulk* positions, but a
real surface almost always relaxes away from them once the missing
neighbors on one side change the local force balance — "surface
relaxation," a real, physically expected effect (distinct from
*reconstruction*, which changes the surface's periodicity/composition, not
just atomic positions). `--ml-relax` lets a MACE potential estimate this
relaxation in seconds, useful groundwork before committing to an expensive
real SIESTA relaxation — see `ml-relax/` below for a measured example.

## When you'd reach for it

- Preparing a clean surface slab for adsorption (`stb-adsorb`), STM/work
  -function (`stb-workfunction`), or surface-band-structure calculations.
- Comparing multiple terminations of the same Miller index before deciding
  which one is worth real DFT time.
- Passivating a covalent semiconductor's dangling bonds before a
  calculation that isn't actually studying the bare surface.
- Getting a MACE-pre-relaxed starting geometry for a defect/adsorption
  study built on top of a cut slab.

## Two ways to run it

A — direct CLI:
```bash
stb-slab -f si_bulk.fdf --hkl 1 0 0
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 2.3
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_2.3.sh` proves this directly at the end.

## What every run does (always on)

- **A numbered report** (`[0] RUN METADATA` … `[7] SUMMARY & FILES` — `[4]`
  only appears with `--ml-relax`) printed to the console and, with
  `--save-report`, also saved to `stb_slab_report.txt`.
- **Structure validation** — atom proximity, lattice handedness, atomic
  density, each an explicit `[OK]`/`[WARNING]`/`[SKIPPED]` row — run once on
  the bulk input and once per written slab.
- **A before/after symmetry comparison table** per written slab — crystal
  system, 3D space group, layer group, point group, Hall symbol, for the
  bulk structure vs. that slab. The bulk always reports `Layer Group: N/A
  (not 2D-periodic)`; the slab (always genuinely 2D-periodic by
  construction) always gets a real one — see `basic/` below.
- **`references.bib`** — SIESTA always, plus the MACE papers if
  `--ml-relax` was used. Merges with whatever `references.bib` is already
  in the working directory.
- **A provenance header** written into every output `.fdf`: the input file,
  Miller index, termination index (polar/non-polar, symmetric/asymmetric),
  passivation detail if used, and MACE convergence/energy if `--ml-relax`
  was used.

## Optional (off by default)

- **`-i`/`--interactive`** — show every termination found (formula, atoms,
  thickness, vacuum, polar, symmetric) and pick one by hand, instead of the
  default (index 0 after sorting non-polar/symmetric first).
- **`--termination N`** — pick a specific termination index directly.
- **`--all`** — write every termination found, one file per termination.
- **`--symmetrize`** — ask `pymatgen` to symmetrize otherwise-asymmetric
  terminations (usually changes both the atom count and how many
  terminations are found — see `symmetrize/` below).
- **`--primitive`** — reduce the bulk to its primitive cell before cutting.
- **`--lll-reduce`** / **`--center-slab`** — pass straight through to
  `pymatgen`'s own `SlabGenerator` options.
- **`--passivate`** (+ `--passivant`/`--cutoff`/`--bond-length`) — cap
  dangling bonds on each written slab. See `passivate/` below.
- **`--ml-relax`** (needs the optional `ml` extra: `pip install
  stb_suite[ml]`) — pre-relaxes each written slab with a MACE potential.
  Positions only by default; add **`--ml-relax-cell`** to also relax the
  in-plane cell — the vacuum axis always stays exactly fixed.
  `--model small/medium/large` (default `small`) or `--custom-model PATH`
  picks which MACE potential to use.
- **`--save-report`** — also persists the full report to
  `stb_slab_report.txt`.
- **`--view`** — opens the bulk structure and every written slab in ASE's
  interactive 3D viewer (`ase-gui`), as pageable frames. Needs a local
  display. Never exercised by `example_2.3.sh` itself.

## Files in this folder

- `si_bulk.fdf` — bulk silicon, diamond cubic, the real conventional 8-atom
  cell (`a` = 5.431 Å, space group Fd-3m No. 227) — covalent, so every cut
  surface leaves real dangling bonds.
- `nacl_bulk.fdf` — bulk NaCl, the real rocksalt structure (FCC Bravais
  lattice, 2-atom basis), conventional 8-atom cubic cell (`a` = 5.64 Å,
  space group Fm-3m No. 225) — ionic, for the polar-surface caveat.
- `example_2.3.sh` — the guided walkthrough (see below).
- `.gitignore` — excludes `output/`, `references.bib`, and
  `stb_slab_report.txt`.

## Running the walkthrough

```bash
cd examples/2.3-stb-slab
./example_2.3.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Six self-contained cases are generated (the
last one is skipped, with a hint, if the optional `ml` extra isn't
installed):

| Folder              | What it shows                                                        |
|----------------------|-----------------------------------------------------------------------|
| `basic/`             | Si(100), default termination — and the bulk-has-no-layer-group/slab-does symmetry contrast |
| `symmetrize/`        | NaCl(111) — one asymmetric termination vs. two symmetric ones with `--symmetrize` |
| `polar-caveat/`      | NaCl(100)/(111) both report `Polar: No` through stb-slab — the oxidation-state gotcha, proven live |
| `passivate/`         | Si(111), `--passivate` — 8/8 dangling bonds genuinely resolved with H  |
| `ml-relax/`          | Si(111), `--ml-relax` — measured surface relaxation                    |
| `full-report/`       | `--save-report` + the validation checklist + `references.bib`         |

## Try it yourself

```bash
# Show every termination and pick by hand
stb-slab -f si_bulk.fdf --hkl 1 1 1 -i

# Open the bulk and the cut slab in ASE's interactive viewer
stb-slab -f si_bulk.fdf --hkl 1 0 0 --view

# Pre-relax with MACE, cell included, with your own fine-tuned model
stb-slab -f si_bulk.fdf --hkl 1 1 1 --ml-relax --ml-relax-cell --custom-model my_finetuned.model

# A thicker slab, symmetrized, every termination
stb-slab -f nacl_bulk.fdf --hkl 1 1 1 --min-slab-size 20 --symmetrize --all
```

## Flag reference

| Flag                | Meaning                                                                |
|----------------------|-------------------------------------------------------------------------|
| `-f/--file`          | Input bulk structure file, `.fdf` (required)                           |
| `--hkl H K L`        | Miller index of the surface to cut (required)                          |
| `--min-slab-size`    | Minimum slab thickness, Å (default `10.0`)                             |
| `--min-vacuum-size`  | Minimum vacuum thickness, Å (default `15.0`)                           |
| `--lll-reduce`       | LLL-reduce the slab lattice                                            |
| `--center-slab`      | Center the slab in the vacuum                                          |
| `--primitive`        | Reduce the bulk to its primitive cell before cutting                   |
| `--symmetrize`       | Ask pymatgen to symmetrize polar/asymmetric terminations               |
| `--symprec`          | Symmetry tolerance, Å — polar/symmetric diagnostics AND the before/after table (default `0.1`) |
| `-i/--interactive`   | Show all terminations and pick one by hand                             |
| `--termination N`    | Pick a specific termination index directly                             |
| `--all`              | Write every termination found                                          |
| `--passivate`        | Cap dangling bonds on each written slab                                |
| `--passivant`        | Element to passivate with, only with `--passivate` (default `H`)       |
| `--cutoff`           | Neighbor-search radius, Å, only with `--passivate` (default: auto)     |
| `--bond-length`      | Passivant bond length, Å, only with `--passivate` (default: auto)      |
| `-o/--output`        | Output `.fdf` base name (default `slab.fdf`)                           |
| `--ml-relax`         | Pre-relax each written slab with a MACE potential before writing       |
| `--ml-relax-cell`    | With `--ml-relax`, also relax the in-plane cell (vacuum axis stays fixed) |
| `--model`            | MACE-MP-0 size for `--ml-relax`: `small`/`medium`/`large` (default `small`) |
| `--custom-model`     | Path to a custom fine-tuned `.model` file for `--ml-relax`             |
| `--save-report`      | Persist the full report (incl. symmetry table) to disk                 |
| `--view`             | Open the bulk structure and every written slab in ASE's interactive viewer |

## What's next

See `2.2-stb-supercell/` for a different structure-transform tool (an
integer tiling matrix instead of a Miller-index cut), or
`1.7-stb-mlrelax/` for a closer look at the MACE pre-relaxation
`--ml-relax` reuses here.
