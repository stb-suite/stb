# 1.2 — Pseudopotential Resolver (`stb-pseudo`)

## What this tool does

Before SIESTA can run a calculation, every chemical species in your
structure needs a matching pseudopotential file (`.psf` or `.psml`) sitting
in the run folder, named after the species *label* declared in `%block
ChemicalSpeciesLabel`. `stb-pseudo` automates that step: give it a
structure (or a plain list of elements) and a pseudopotential source, and
it tells you exactly which elements are covered, which aren't, and copies
the ones that are into a folder ready for `calc.fdf`. `--fallback-dir`
lets you patch gaps from a second source automatically, instead of
hunting down a missing file by hand.

## Why this matters (a bit of theory)

A real atom's electrons split into two very different roles: a handful of
**core** electrons, tightly bound near the nucleus, essentially inert to
chemistry, and a few **valence** electrons, which is where all the actual
bonding physics happens. Solving the all-electron problem exactly means
resolving both — including the core electrons' deep, rapidly-oscillating
wavefunctions near the nucleus, which is numerically expensive and, for
chemistry, mostly wasted effort.

A **pseudopotential** replaces the true nucleus-plus-core-electrons
Coulomb potential with a smoother *effective* potential, tuned so that
outside a chosen core radius `r_c` the resulting pseudo-wavefunctions
reproduce the real (all-electron) valence wavefunctions almost exactly.
Inside `r_c` the pseudo-wavefunction is deliberately smooth (no more sharp
oscillations to resolve), which is what actually makes the calculation
cheap.

**Norm-conserving** pseudopotentials — the kind both banks `stb-pseudo`
ships (see below) provide — add one more constraint: the pseudo- and
all-electron wavefunctions must enclose the *same charge* inside `r_c`.
This single constraint is what makes a pseudopotential **transferable**:
it keeps the scattering properties (logarithmic derivatives) correct not
just at one reference energy, but over the whole range of energies the
atom actually experiences in different chemical environments. Without it,
a pseudopotential fitted for one bonding situation could silently give
wrong results in another — transferability is the reason a pseudopotential
library can be trusted as a general-purpose building block, rather than
needing a fresh all-electron calculation for every new compound.

Two format notes that matter in practice:

- **`.psf`** is SIESTA's older native pseudopotential format (the
  `virtual_vault` bank ships these). **`.psml`** is a newer,
  code-agnostic exchange format that SIESTA also reads directly (the
  `dojo` bank ships these). Numerically they describe the same kind of
  object — `stb-pseudo` copies whichever is present, preferring `.psml`
  when both exist for a species.
- A pseudopotential is generated for one specific **exchange-correlation
  functional** (both bundled banks here are PBE). Using it in a
  calculation set up for a *different* functional is a real, silent
  correctness trap — `stb-pseudo` does not check this (it has no way to
  know what functional your `calc.fdf` will use), so it's on you to keep
  the two consistent.

## When you'd reach for it

Right before (or right after) `stb-inputfile` (example `1.1`), to make
sure every element your structure needs actually has a pseudopotential
available — and to get them all copied into the run folder in one step,
instead of resolving each one by hand. It's also the tool to reach for
when a bundled bank is missing one particular element: `--fallback-dir`
lets a second source fill just that gap, without abandoning the first
source for everything else.

## Two ways to run it

**A — direct CLI**:

```bash
stb-pseudo --file sic.fdf --pp-path dojo --output out/
```

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 1.2
```

`1.2` first asks what you want to do — resolve from a structure file
**[1]**, resolve from an explicit element list **[2]**, or just browse a
bundled bank's contents **[3]** — then walks through the same primary
source / fallback source / output directory / dry-run / save-report
questions the CLI flags cover. `example_1.2.sh` proves the CLI and the
menu agree.

## Files in this folder

- `sic.fdf` — a 2-atom SiC primitive cell (Si + C, both covered by the
  bundled `dojo` bank).
- `example_1.2.sh` — the guided walkthrough (**not** an automated test —
  see `test/1-inputs/2-pseudo/test.sh` for that). Pauses between sections
  so you can read before moving on; safe to re-run.
- `output/` — created by `example_1.2.sh` when you run it (git-ignored,
  not checked in). See below.

## Running the walkthrough

```bash
./example_1.2.sh
```

| Folder                | Command                                                    | What it shows |
|------------------------|-------------------------------------------------------------|----------------|
| *(none — just prints)* | `stb-pseudo --list-elements dojo`                          | browsing a bank without any structure |
| `output/sic/`          | `stb-pseudo -f sic.fdf -p dojo -o output/sic`               | resolve + copy from a real structure |
| `output/fallback/`     | `stb-pseudo --species Si At -p dojo --fallback-dir virtual_vault -o output/fallback` | filling a gap `dojo` can't cover |
| `output/dry_run/`      | `stb-pseudo -f sic.fdf -p dojo -o output/dry_run --dry-run` | report only, nothing copied |

### Bundled banks: what's actually inside `dojo` vs `virtual_vault`

`stb-pseudo --list-elements` (used above with no structure at all) is the
quickest way to check whether a bank covers what you need *before*
committing to it:

- **`dojo`** — PseudoDojo v0.5, PBE, `.psml` format, 72 elements.
- **`virtual_vault`** — the SIESTA Pseudopotentials Virtual Vault, PBE,
  `.psf` format, a different (overlapping but not identical) element
  coverage.

Neither bank covers every element you might need — which is exactly why
`--fallback-dir` exists, rather than forcing you to pick one bank and
live with its gaps.

### Filling the gaps: `--fallback-dir` made concrete

Astatine (`At`) is a real, if obscure, case: it's absent from `dojo` but
present in `virtual_vault`. Resolving `Si At` against `dojo` alone reports
`At` as `MISSING`; adding `--fallback-dir virtual_vault` resolves it from
the second source instead, and the report says exactly which source
supplied each element (`primary` vs `fallback`) — nothing is silently
mixed up.

### Proof: CLI and the interactive menu agree

The script also drives the same `sic.fdf` / `dojo` case through the
interactive `stb-suite` menu (non-interactively, via a piped `printf`) and
diffs the resulting "Resolved : N/N" line against the direct CLI run —
proving the two paths produce identical output.

## Try it yourself

```bash
stb-pseudo --species Fe O -p dojo -o out/          # a common oxide pair
stb-pseudo --list-elements virtual_vault           # what's in the other bank?
```

## Flag reference

```
stb-pseudo (-f/--file <structure.fdf> | --species EL [EL ...]) -p/--pp-path <BANK_OR_PATH>
           [--fallback-dir <BANK_OR_PATH>] [-o/--output <dir>] [--dry-run] [--save-report]
stb-pseudo --list-elements <BANK>
```

- `-f/--file`: a SIESTA structure file (`.fdf`) to read the required
  species from. Mutually exclusive with `--species`.
- `--species`: element symbol(s) to resolve directly, instead of reading
  a structure file.
- `-p/--pp-path` (required unless `--list-elements` is given): a bundled
  bank name (`dojo`, `virtual_vault`) or a directory of your own
  `.psf`/`.psml` files.
- `--fallback-dir` (optional): a second source, consulted only for
  elements missing from `-p/--pp-path`.
- `-o/--output` (optional, default `.`): directory the resolved
  pseudopotentials are copied into.
- `--dry-run` (optional): only run the availability check — copy
  nothing.
- `--list-elements <BANK>`: print the elements available in a bundled
  bank and exit — no structure or `-p/--pp-path` needed.
- `--save-report` (optional): also persist the report to
  `stb_pseudo_report.txt`. Off by default.

`references.bib` (SIESTA plus a citation for every bank actually used) is
always written — there's no flag for it.

Run `stb-pseudo --help` for the full list of options.

## What's next

Point `stb-inputfile`'s (example `1.1`) `-p/--pp-path`, or your
`calc.fdf`'s pseudopotential setup directly, at the output directory
`stb-pseudo` just filled — every species it resolved is already sitting
there with the right filename.
