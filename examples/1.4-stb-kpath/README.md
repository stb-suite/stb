# 1.4 — K-Path Generator (`stb-kpath`)

## What this tool does

For a band-structure calculation, SIESTA needs a path of k-points through
the Brillouin zone (BZ) — the `%block BandLines` in `kpath_bs.fdf`.
`stb-kpath` builds one automatically: given a structure, it detects its
dimensionality and Bravais lattice, then generates the correct
high-symmetry path for that lattice type, dimension-aware (a 1D wire, a 2D
slab, and a 3D bulk crystal each have a genuinely different Brillouin zone
shape, hence a different path).

## Why this matters (a bit of theory)

`stb-kgrid` (example `1.3`) builds a Monkhorst-Pack **grid** — k-points
spread uniformly over the *whole* BZ, used to approximate an *integral*
(total energy, charge density). A band structure plot is a different job:
it shows how the energy `E(k)` varies *along* a route through the BZ, so
it's traced along a **path** connecting the BZ's high-symmetry points
instead — `E(k)` is smooth, and the physically interesting features (band
extrema, gaps, crossings) are captured well by a handful of well-chosen
points and the segments between them.

Which points count as "high-symmetry", and in what order to connect them,
isn't arbitrary — it depends on the shape of the BZ, which in turn depends
on the crystal's **Bravais lattice**. **Setyawan & Curtarolo (2010)** — the
exact reference `stb-kpath` writes to `references.bib` — standardized this
into one path convention per Bravais lattice type, now used across
essentially every DFT code for reproducible band-structure plots. ASE
implements that exact convention (`Cell.bandpath`) and extends it to 1D/2D
systems via a periodic-axes mask; `stb-kpath` is a thin, dimension-aware
wrapper around it — it does not invent or choose paths itself.

### Bravais lattice vs. space group — two different questions

`stb-kpath` reports both, and they answer different things:

- **Bravais lattice** — from the raw *lattice vectors* alone (ignoring the
  atoms). Determines the BZ shape and which path template applies. This is
  what the k-path is actually built from.
- **Space group** — the *full* crystallographic symmetry, lattice + atomic
  basis together (via the same `core/symmetry.py` accessors
  `stb-inputfile` also uses).

These can legitimately disagree: a structure file can list its atoms in a
non-primitive (conventional) cell, so the raw lattice vectors alone look
lower-symmetry than the real crystal actually is once the atomic
arrangement is taken into account. `example_1.4.sh` shows exactly this on
bulk silicon: a plain cubic Bravais lattice, but an `Fd-3m` (diamond) space
group.

### Why a 0D structure is a hard error

Bloch's theorem, and the whole idea of a k-point, only applies along a
genuinely periodic direction (see example `1.3`'s own theory section). An
isolated molecule (0D — vacuum-padded on every axis) has no periodic
direction at all, so a k-path isn't just low quality there, it's not
physically meaningful. `stb-kpath` treats this as a real error (exit code
1, nothing written) rather than a soft warning.

## When you'd reach for it

Whenever you need `kpath_bs.fdf` for a band-structure calculation. In
practice you rarely run this standalone — `stb-inputfile`'s `bands` mode
(example `1.1`) already calls `stb-kpath` itself to produce it. `stb-kpath`
is for inspecting or double-checking the path it chose, or for generating
one on its own.

## Two ways to run it

**A — direct CLI**:

```bash
stb-kpath --file structure.fdf
```

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 1.4
```

`1.4` asks for the structure file, then — the same "advanced settings"
pattern used elsewhere in the suite — offers to configure the vacuum-gap,
the Bravais-lattice tolerance, and the symmetry tolerances, all gated
behind a single `[y/N]` question so the essential flow stays short.
`example_1.4.sh` proves both paths agree.

## Files in this folder

- `structure.fdf` — bulk-silicon example structure (same 8-atom
  conventional cubic cell used throughout `examples/1.1-stb-inputfile/`
  and `examples/1.3-stb-kgrid/`).
- `structure_chain.fdf` / `structure_graphene.fdf` — a carbon chain (1D)
  and a graphene monolayer (2D).
- `structure_molecule.fdf` — an isolated CH4 molecule (0D), used to
  demonstrate the hard error above.
- `example_1.4.sh` — the guided walkthrough (**not** an automated test —
  see `test/1-inputs/4-k-path/test.sh` for that). Pauses between sections
  so you can read before moving on; safe to re-run.
- `output/` — created by `example_1.4.sh` when you run it (git-ignored, not
  checked in). See below.

## Running the walkthrough

```bash
./example_1.4.sh
```

Every successful case is generated with `--save-report`, each into its own
folder under `output/` — `stb_kpath_report.txt`, `references.bib` (SIESTA +
Setyawan-Curtarolo citations), and `kpath_bs.fdf`:

| Folder             | Structure                           | Dimensionality | Bravais lattice | Path |
|---------------------|---------------------------------------|:---:|-----------------|------|
| `output/silicon/`   | `structure.fdf` (bulk silicon)        | 3D | primitive cubic (CUB) | `Γ-X-M-Γ-R-X \| M-R` |
| `output/chain/`     | `structure_chain.fdf` (carbon chain)  | 1D | primitive line (LINE) | `Γ-X` |
| `output/graphene/`  | `structure_graphene.fdf` (graphene)   | 2D | primitive hexagonal (HEX2D) | `Γ-M-K-Γ` |
| `output/molecule/`  | `structure_molecule.fdf` (CH4)        | 0D | *(hard error — no path, no file)* | — |

For the 1D/2D cases, the script also prints the same vacuum-padded-axis
caveat `stb-inputfile` uses: the space group shown treats the vacuum as an
ordinary periodic direction and may not be exact — use `stb-symmetry`
(code `3.5`) for a proper layer-group/point-group analysis of those.

### The `--vacuum-gap` threshold changes the physics, not just a number

Same idea as example `1.3`: a structure with a 12 Ang gap along `c` sits
just above the default 10 Ang vacuum threshold (2D, square BZ, path
`M-Γ-X-M`) but below a stricter 15 Ang one (3D, tetragonal BZ, a
completely different, longer path) — same atoms, different threshold,
different Brillouin zone entirely.

### Proof: CLI and the interactive menu agree

The script also drives the same silicon case through the interactive
`stb-suite` menu's advanced-settings gate (vacuum-gap/eps/symprec/angle all
entered explicitly at their defaults) and diffs the resulting path against
`output/silicon/` — proving the CLI and the menu produce identical output.

## Try it yourself

```bash
stb-kpath --file structure.fdf --vacuum-gap 5    # a tighter vacuum threshold
stb-kpath --file structure_chain.fdf -o my_path.fdf
```

## Flag reference

```
stb-kpath --file <structure.fdf> [-p/--prec <eps>] [--vacuum-gap <Ang>]
          [--symprec <tol>] [--angle-tolerance <deg>] [-o/--output <file>]
          [--save-report]
```

- `-f/--file` (required): path to a SIESTA structure file (`.fdf` only).
- `-p/--prec` (optional, default `0.0002`): tolerance for ASE's
  Bravais-lattice detection.
- `--vacuum-gap` (optional, default `10.0`): minimum empty span (Ang) along
  an axis, wrapped periodically, to treat it as vacuum-padded.
- `--symprec`/`--angle-tolerance` (optional, defaults `0.001`/`5.0`):
  tolerances for the crystallographic space-group analysis.
- `-o/--output` (optional, default `kpath_bs.fdf`): output filename.
- `--save-report` (optional): also persist the report to
  `stb_kpath_report.txt`. Off by default.

`references.bib` (SIESTA + Setyawan-Curtarolo citations) is always written
whenever a path is actually generated — there's no flag for it.

Run `stb-kpath --help` for the full list of options.

## What's next

`kpath_bs.fdf` is exactly what `stb-inputfile`'s `bands` mode (example
`1.1`) references via `%include kpath_bs.fdf` — but you rarely need to
generate it separately: that example's own `bands/` step already calls
`stb-kpath` for you, right after `stb-inputfile` itself.
