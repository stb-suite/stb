# 1.5 — DFT+U / Hubbard Block Generator (`stb-dftu`)

## What this tool does

For materials with localized d or f electrons (transition-metal oxides,
lanthanides, actinides), plain DFT often gets the electronic structure
wrong. `stb-dftu` generates a ready-to-use `%block LDAU.proj` snippet for
SIESTA's DFT+U (Hubbard) correction — one stanza per species, with the
correlated shell auto-detected and the Hubbard `U` (and exchange `J`)
supplied by you. It never estimates `U` on its own; it only offers a
literature table as an optional, clearly-labeled starting point.

## Why this matters (a bit of theory)

Ordinary DFT (LDA/GGA) has a **self-interaction error**: an electron partly
interacts with itself. For most materials this barely matters, but for
*localized* d or f electrons it artificially favors delocalization — the
textbook failure is **NiO**, a real antiferromagnetic insulator that plain
GGA predicts as a metal.

**DFT+U** (the Hubbard correction) patches this with a mean-field term that
penalizes *fractional* occupation of the correlated orbital, pushing it
toward an integer occupation (0 or 1) — closer to how a real localized
electron behaves. It's a cheap, physically-motivated patch, not a
first-principles fix. Two parameters describe it:

- **`U`** — on-site Coulomb repulsion: the energy cost of putting two
  electrons in the same correlated orbital.
- **`J`** — Hund's-rule exchange energy (usually smaller than `U`).

The SIESTA block also needs `n`/`l` — the correlated shell's principal and
angular quantum numbers (`l=2` for a d shell, `l=3` for f) — auto-detected
by periodic-table block: transition metals → `nd`, lanthanides → `4f`,
actinides → `5f` (see `core/dftu_data.py::SHELL_NAMES`/`DEFAULT_SHELL`).

### Why `stb-dftu` never guesses `U` for you

Unlike `n`/`l` (uncontroversial chemistry), `U` is *not* a universal
constant — it depends on the material, the functional, the
pseudopotential, and the basis set (screening, hybridization with
neighboring ligands all shift it). `stb-dftu` bundles one literature table
— **Wang, Maxisch & Ceder, Phys. Rev. B 73, 195107 (2006)**, the Materials
Project's own GGA+U oxide calibration, and the exact reference this tool
writes to `references.bib` — purely as a starting-point *sanity check*,
never a validated value. You always have to pass `--u` explicitly, or opt
in with `--use-reference`.

For a real, first-principles `U` computed from *your own* system (not a
generic oxide table), use `stb-hubbardu`/`stb-hubbarduAnalysis` (Workflow
menu) instead — the Cococcioni & de Gironcoli linear-response method.
`stb-dftu` itself already points there in its own `--help` text.

## When you'd reach for it

Any time you're setting up a SIESTA calculation for a transition-metal,
lanthanide, or actinide compound and need a `DFT+U` correction — pasting
the generated block into `calc.fdf`'s DFT+U section (the same file
`stb-inputfile`, example `1.1`, produces).

## Two ways to run it

**A — direct CLI**:

```bash
stb-dftu --species Mn --u 3.9
```

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 1.5
```

`1.5` itself offers two sub-modes: **[1]** read species off a structure
file and auto-fill `U` from the reference table, or **[2]** enter species
one at a time by hand. Both end with the same save-block/save-report
prompts. `example_1.5.sh` proves the CLI and the menu agree.

## Files in this folder

- `sc_fe_o.fdf` — a synthetic 3-species fixture (Sc, Fe, O) used to
  demonstrate `--fdf --use-reference`'s 3-way split (see below).
- `example_1.5.sh` — the guided walkthrough (**not** an automated test —
  see `test/1-inputs/5-dftu/test.sh` for that). Pauses between sections so
  you can read before moving on; safe to re-run.
- `output/` — created by `example_1.5.sh` when you run it (git-ignored, not
  checked in). See below.

## Running the walkthrough

```bash
./example_1.5.sh
```

Every case is generated with `--save-report`, each into its own folder
under `output/` — `stb_dftu_report.txt` + `references.bib`:

| Folder                    | Command (conceptually)                                          |
|----------------------------|-------------------------------------------------------------------|
| `output/single/`           | `stb-dftu --species Mn --u 3.9` (auto-detected 3d shell)          |
| `output/multi/`            | `stb-dftu --species Fe Co --u 5.3 3.32 --j 0.5 0.5 --shell 3d 3d` |
| `output/from-structure/`   | `stb-dftu --fdf sc_fe_o.fdf --use-reference`                      |
| `output/list-reference/`   | `stb-dftu --list-reference` (no block, table only)                |
| `output/suggest/`          | `stb-dftu --suggest Ni` (no block, one value only)                |

### `from-structure/` — the 3-way split, live

`sc_fe_o.fdf` has 3 species: **Sc** (a metal with *no* tabulated
reference), **Fe** (a metal *with* one), and **O** (a non-metal).
`--use-reference` treats all three differently:

- **Fe** → tabulated, included in the block, with a loud `[WARNING]` that
  a literature (not validated) value is being used.
- **Sc** → a metal, but untabulated → skipped, with its own `[WARNING]` so
  it's never silently dropped.
- **O** → not a metal at all → silently skipped (DFT+U doesn't apply to a
  non-correlated species, no warning needed).

### `list-reference/`/`suggest/` — lookups only

Neither generates a block, so neither's `references.bib` gets a SIESTA
citation (nothing SIESTA-bound was produced) — just the Wang/Maxisch/Ceder
reference-table citation.

### "Never guesses" in practice

The script also runs `stb-dftu --species Si --u 2.0` — Si has no default
correlated shell and no tabulated reference, so this is a real error (exit
code 1, no block, no silent fallback), not a guess.

## Try it yourself

```bash
stb-dftu --species Ni --use-reference     # single species, auto-filled U
stb-dftu --species Cr --u 3.7 --shell 3d
```

## Flag reference

```
stb-dftu (--species <label>... --u <eV>... | --fdf <path> [--use-reference])
         [--j <eV>...] [--shell <3d|4d|5d|4f|5f>...] [-o/--output <file>]
         [--list-reference] [--suggest <element>] [--save-report]
```

- `--species`/`--fdf` — one is required (unless `--list-reference`/
  `--suggest`); `--fdf` reads species off a structure file, `--species`
  filters/overrides.
- `--u` (eV, one per species) — required unless `--use-reference`.
- `--use-reference` — auto-fill `U` from the bundled literature table for
  any species without an explicit `--u`.
- `--j` (eV, default `0.0`) / `--shell` (default: auto-detected).
- `-o/--output` — also save the block to a file.
- `--list-reference`/`--suggest <element>` — print the reference table (or
  one value) and exit; no block generated.
- `--save-report` (optional): also persist the report to
  `stb_dftu_report.txt`. Off by default.

`references.bib` is always written — there's no flag for it.

Run `stb-dftu --help` for the full list of options.

## What's next

Paste the generated block into `calc.fdf`'s DFT+U section (the file
`stb-inputfile`, example `1.1`, produces). If you need a real,
system-specific `U` instead of the bundled literature table, see
`stb-hubbardu`/`stb-hubbarduAnalysis` in the Workflow menu.
