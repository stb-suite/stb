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
structure file, which calculation mode (numbered list), an optional
pseudopotential source, and whether to save a text report — building the
exact same `stb-inputfile ... -t ... -p ...` command underneath and running
it for you. Good for exploring what a tool can do before you've memorized
its flags; the CLI is faster once you have.

Both paths run the identical underlying tool and produce byte-identical
output — `example_1.1.sh` demonstrates this directly, running the same
`relax` calculation both ways and diffing the results.

## Files in this folder

- `structure.fdf` — bulk-silicon example structure (conventional cubic
  diamond cell, 8 atoms, 3D/fully periodic).
- `structure_0D.fdf` / `structure_1D.fdf` / `structure_2D.fdf` — an isolated
  CH4 molecule, an infinite carbon chain, and a graphene monolayer.
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

| Folder           | Command (conceptually)                              |
|-------------------|------------------------------------------------------|
| `output/relax/`   | `stb-inputfile structure.fdf -t relax -p dojo`       |
| `output/bands/`   | `stb-inputfile structure.fdf -t bands -p dojo`, then `stb-kpath -f structure.fdf` to generate the `kpath_bs.fdf` that `bands` mode references |
| `output/0D/`      | `stb-inputfile structure_0D.fdf -t total_energy -p dojo` |
| `output/1D/`      | `stb-inputfile structure_1D.fdf -t total_energy -p dojo` |
| `output/2D/`      | `stb-inputfile structure_2D.fdf -t total_energy -p dojo` |
| `output/3D/`      | `stb-inputfile structure.fdf -t total_energy -p dojo`    |

The script also drives the exact same `relax` case through the interactive
`stb-suite` menu (non-interactively, via a piped `printf`, no folder kept
for it) and diffs the result against `output/relax/calc.fdf` — proving the
CLI and the menu produce identical output.

### Why `bands/` needs an extra command

`stb-inputfile`'s `bands` mode writes a `calc.fdf` that references
`%include kpath_bs.fdf`, but doesn't generate that file itself — that's
`stb-kpath`'s (code 1.3) job. The script calls it too, so `output/bands/`
ends up genuinely complete instead of missing a file SIESTA would need.

### Why dimensionality matters

A k-point grid only makes physical sense along a genuinely periodic
direction. `stb-inputfile` measures the largest empty gap along each
lattice vector and, above a threshold (10 Ang), treats that axis as vacuum
instead of periodic when computing the k-grid — automatically, from the
geometry alone, no flag needed:

| Dimensionality | Structure                          | Periodic axes | Typical k-grid shape   |
|-----------------|-------------------------------------|:---:|-------------------------|
| 0D              | `structure_0D.fdf` (CH4 molecule)   | 0 | `1  1  1`   (isolated)   |
| 1D              | `structure_1D.fdf` (carbon chain)   | 1 | `N  1  1`   (a wire)     |
| 2D              | `structure_2D.fdf` (graphene)       | 2 | `N  N  1`   (a slab)     |
| 3D              | `structure.fdf` (bulk silicon)      | 3 | `N  N  N`   (a crystal)  |

This is why the same `stb-inputfile <structure> -t total_energy` command
works unchanged whether the structure is a molecule, a wire, a slab, or a
bulk crystal — the vacuum gap in the structure file itself is the signal,
so you never need a separate "molecule mode" or "slab mode" flag.

## Try it yourself

Once you've read through the walkthrough, experiment directly:

```bash
stb-inputfile structure.fdf -t total_energy -p dojo
stb-inputfile structure.fdf -t aimd -p dojo
stb-inputfile structure.fdf -t relax+d3 -p dojo   # any mode + "+d3" adds a DFT-D3 correction
```

Open each generated `calc.fdf` and compare it against `structure.fdf` and
against the other modes — that's the fastest way to see exactly what each
mode adds.

## Flag reference

```
stb-inputfile <structure_file> -t <mode> [-p <pseudo-source>]
```

- `-t/--type` (required): `total_energy`, `relax`, `aimd`, or `bands` — each
  also accepts a `+d3` suffix (e.g. `relax+d3`) to enable the DFT-D3
  dispersion correction.
- `-p/--pp-path` (optional): a bundled bank name (`dojo`, ...) or a path to a
  folder containing `.psml`/`.psf` pseudopotential files.

Run `stb-inputfile --help` for the full list of options.

## What's next

Any `output/<case>/` folder is ready to hand to SIESTA as-is. After the run
finishes, `stb-status` gives you a quick summary of what happened
(converged? final energy? which output files were produced), and the
`Analysis` category tools (`stb-bands`, `stb-dos`, ...) pick up from there.
