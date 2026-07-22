# 1.5 — Structure Fetcher (`stb-fetch`)

## What this tool does

`stb-fetch` looks up a structure in an online database — the
Crystallography Open Database (COD, no account needed), Materials Project
(needs a free `PMG_MAPI_KEY`), or any OPTIMADE-compliant database (AFLOW,
JARVIS, 2D-materials databases like `twodmatpedia`, ...) — and writes it
straight out as a SIESTA `.fdf`. It's the suite's first network-dependent
tool: everything else processes files you already have locally.

Real database entries are frequently disordered even for common materials
(e.g. magnetite's octahedral Fe site is often given as a Fe2+/Fe3+ mixed
-occupancy split rather than plain Fe). Same-element-different-oxidation
-state disorder is collapsed automatically to a single element at full
occupancy; genuine multi-element disorder (two different elements sharing a
site) is refused with a clear error instead of silently guessing.

## When you'd reach for it

Any time the material you want already exists in the literature and you'd
rather fetch it than hand-build a `POSCAR`/CIF/`.fdf` yourself — as a
starting point for `stb-inputfile` (example `1.1`), or for `stb-symmetry`/
`stb-kpath` once you have the structure in hand.

## Two ways to run it

**A — direct CLI**:

```bash
stb-fetch --source cod --cod-id 1010369
```

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 1.5
```

`1.5` walks you through picking a source/provider, an exact id or a formula
search, an optional unit-cell reduction, and the save-report/`--view`
prompts — the same choices the CLI flags below control directly.
`example_1.5.sh` proves the CLI and the menu agree.

## Files in this folder

- `example_1.5.sh` — the guided walkthrough. Pauses between sections so you
  can read before moving on; safe to re-run. **Needs internet access** —
  unlike every other example so far, this tool queries COD/OPTIMADE live;
  if a request fails (offline, or a remote endpoint is down), re-run once
  connectivity is back.
- `output/` — created by `example_1.5.sh` when you run it (git-ignored, not
  checked in). See below.

## Running the walkthrough

```bash
./example_1.5.sh
```

| Folder                      | Command (conceptually)                                                          |
|------------------------------|----------------------------------------------------------------------------------|
| `output/cod-by-id/`          | `stb-fetch --source cod --cod-id 1010369` (magnetite, disordered, bulk 3D)      |
| `output/cod-primitive/`      | same, + `--unitcell primitive`                                                  |
| `output/optimade-2d/`        | `stb-fetch --source optimade --provider twodmatpedia --optimade-id 2dm-3150` (monolayer MoS2, 2D) |
| `output/formula-search/`     | `--formula MoS2` with no exact id — multiple candidates, informative error      |
| `output/merged-citations/`   | two fetches (COD + OPTIMADE) into the same folder — `references.bib` merges     |

Every case that generates a structure is run with `--save-report`, each
into its own `stb_fetch_report.txt` + `references.bib`.

### Structure validation (always on, every fetch)

Every run prints a `[3] STRUCTURE VALIDATION` section, right after the
structure is fetched and disorder-resolved. It always runs — there's no
flag to turn it on, and it never blocks writing the `.fdf`. It checks the
same three malformation cases `stb-inputfile` (example `1.1`) checks on a
hand-built structure:

1. **Atoms too close together** — any two atoms closer than 0.5 Å
   (minimum-image distance).
2. **Left-handed cell** — the lattice vectors' determinant is negative.
3. **Implausible atomic density** — atoms-per-Å³ outside roughly
   `[0.01, 0.15]`. **Only checked for a genuine 3D bulk structure** — a
   vacuum-padded cell's volume is dominated by the artificial vacuum by
   design, so this metric would be meaningless there.

Real database entries are usually clean, so these three rarely fire in
practice — but a fetched structure gets exactly the same scrutiny as one
you built by hand.

### Dimension-aware group label

Right after the malformation checks, the same section reports the
structure's symmetry group — **the label only**, not the full list of
operations or a Wyckoff-site table (that's `stb-symmetry`, code `3.5`'s
job):

- **Bulk (3D, no vacuum axis)**: the ordinary space group, e.g.
  `Fd-3m (No. 227)`.
- **Slab (exactly one vacuum-padded axis, e.g. a 2D material)**: the
  **layer group** instead — the physically correct classification for a
  genuinely 2D-periodic structure (an ordinary 3D space group would treat
  the vacuum gap as just an unusually tall periodic cell). Same detection
  `stb-symmetry` uses (`spglib.get_layergroup`), needs `spglib >= 2.1.0`.
- **Wire (two vacuum-padded axes)**: neither spglib nor pymatgen has a
  rod-group equivalent, so this case only gets a note — same documented
  limitation as `stb-symmetry`.
- **Isolated molecule (all three axes vacuum-padded)**: the **point group**
  instead (pymatgen's `PointGroupAnalyzer`). In practice, none of COD/
  Materials Project/OPTIMADE hand back an isolated 0D structure — they're
  crystal-structure databases — so this case is implemented for
  completeness/robustness but isn't reachable through the walkthrough.

### `references.bib` — always written, and it merges

`references.bib` is always written — there's no flag for it. Citations
always include SIESTA, spglib (symmetry detection now always runs), and
one entry for whichever database the structure came from (COD/Materials
Project/OPTIMADE). Fetching from two different sources into the same
output folder merges citations by BibTeX key instead of the second run
erasing the first's — see `output/merged-citations/`.

### `--view` (optional, off by default)

Opens the final structure in ASE's interactive 3D viewer after writing the
`.fdf` — needs a display. Not run by the walkthrough script (no display in
a CI/headless environment); try it yourself:

```bash
stb-fetch --source cod --cod-id 1010369 --view
```

## Try it yourself

```bash
stb-fetch --list-providers                    # every curated OPTIMADE alias
stb-fetch --source cod --formula Si --limit 5 # a formula search
stb-fetch --source materials-project --material-id mp-19306 --api-key <your PMG_MAPI_KEY>
```

## Flag reference

```
stb-fetch --source <cod|materials-project|optimade>
          (--material-id <id> | --cod-id <id> | --optimade-id <id> | --formula <formula>)
          [--most-stable] [--api-key <key>] [--provider <alias|url>]
          [--unitcell <primitive|conventional|refined>] [--symprec <val>] [--angle-tolerance <deg>]
          [-o/--output <file>] [--save-report] [--view]
```

- `--source` (required unless `--list-providers`): which database to query.
- One of `--material-id`/`--cod-id`/`--optimade-id`/`--formula` is
  required — exact id, or a formula search that may need a follow-up pick.
- `--most-stable` (materials-project + `--formula` only): auto-pick the
  lowest energy-above-hull candidate.
- `--unitcell` (optional): also reduce to the primitive/conventional/
  refined cell (same logic as `stb-unitcell`) before writing.
- `-o/--output` (default `fetched.fdf`).
- `--save-report` (optional): also persist the report to
  `stb_fetch_report.txt`. Off by default.
- `--view` (optional): open the structure in ASE's interactive 3D viewer
  after writing. Off by default; needs a display.

Structure validation and the dimension-aware group label (see above) always
run before writing — there's no flag for either, they're always on.
`references.bib` is always written too.

Run `stb-fetch --help` for the full list of options.

## What's next

The `.fdf` this tool writes is a plain structure file — feed it to
`stb-inputfile` (example `1.1`) to build a full SIESTA calculation, or to
`stb-symmetry` (code `3.5`) for the full operations/Wyckoff-site symmetry
analysis this tool only ever summarizes as a single label.
