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
  code-agnostic exchange format that SIESTA also reads directly (every
  `dojo` table ships these). Numerically they describe the same kind of
  object — `stb-pseudo` copies whichever is present, preferring `.psml`
  when both exist for a species.
- A pseudopotential is generated for one specific **exchange-correlation
  functional** (PseudoDojo has PBE, PBEsol and LDA tables; the Virtual
  Vault is PBE). Using it in a calculation set up for a *different*
  functional is a silent correctness trap: SIESTA runs anyway. Every tool
  that copies pseudopotentials into a run folder therefore reads the
  functional (and the relativistic treatment) written in each file and
  compares it with the folder's own `.fdf` (`XC.Functional`/`XC.Authors`,
  `Spin spin-orbit`), printing a `[WARNING]` on a mismatch. It never
  stops the tool: the decision stays with you.
- **Spin-orbit coupling** needs *full-relativistic* pseudopotentials
  (with `lj` projectors). With scalar-relativistic files SIESTA stops
  with `The PSML file does not contain lj projs`; the `dojo-fr-*` tables
  are the ones to use.

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
**[1]**, resolve from an explicit element list **[2]**, browse a bundled
bank's contents **[3]**, or list every bundled bank **[4]** — then walks
through the same primary source (PseudoDojo: relativity, functional,
accuracy) / fallback source / output directory / dry-run
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

### Bundled banks: the PseudoDojo tables and `virtual_vault`

The PseudoDojo bank holds eleven tables (norm-conserving ONCVPSP, `.psml`).
You choose one by relativity, functional and accuracy, with a single name
given to `-p` (here and in every other tool that takes pseudopotentials):

```
dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent]
```

The words go in any order, and each missing word takes its default
(`sr`, `pbe`, `standard`). `fr` can also be written `soc`, and `pw` can be
written `lda`. `stb-pseudo --list-banks` prints the whole list:

| `-p` name | PseudoDojo table | Relativity | Functional | Elements |
|---|---|---|---|---|
| `dojo` = `dojo-sr-pbe-standard` | nc-sr-05 pbe standard | scalar | PBE | 72 |
| `dojo-stringent` | nc-sr-05 pbe stringent | scalar | PBE | 72 |
| `dojo-pbesol`, `dojo-pbesol-stringent` | nc-sr-04 pbesol | scalar | PBEsol | 72 |
| `dojo-pw`, `dojo-pw-stringent` | nc-sr-04 pw | scalar | LDA (PW92) | 70 |
| `dojo-fr` | nc-fr-04 pbe standard | full (SOC) | PBE | 70 (no La, Lu) |
| `dojo-fr-stringent` | nc-fr-04 pbe stringent | full (SOC) | PBE | 72 |
| `dojo-fr-pbesol`, `dojo-fr-pbesol-stringent` | nc-fr-04 pbesol | full (SOC) | PBEsol | 71 (no Lu) |
| `dojo-sr3plus` | nc-sr-04-3plus pbe standard | scalar | PBE | 14 (Ce-Lu) |
| `virtual_vault` | SIESTA Virtual Vault (`.psf`) | — | PBE | 69 |

- **Standard or stringent.** Standard is the usual choice. For groups
  13–16 (Ga–Se, In–Te, Tl–Po), standard *already* has the d shell in the
  valence; stringent adds the s and p shell below it (Ge and Sn go from
  z = 14 to z = 22). For Hf and Ta, stringent adds the 4f shell (z 12 → 26
  and 13 → 27). More electrons and more orbitals make the calculation more
  expensive. The report shows each file's valence (`z=`).
- **`fr`** tables are required for `Spin spin-orbit`. They also run
  without SOC, but do not compare absolute energies from `fr` and `sr`
  files. PseudoDojo has no full-relativistic LDA table.
- **`sr3plus`** holds the lanthanides Ce–Lu with the 4f electrons frozen in
  the core. Use it only for Ln3+ compounds: it does not describe f physics
  (magnetism, 4f states). Give the other elements a fallback, for example
  `-p dojo-sr3plus --fallback-dir dojo`.
- `dojo` names exactly the files this bank held before it was split into
  tables (PseudoDojo v0.5, PBE, standard), so old commands, scripts and
  `mlff_workflow.json` files that say `dojo` give the same result.

In the interactive menu, every question about the pseudopotential source
offers *PseudoDojo* and then asks for the relativity, the functional and the
accuracy (Enter keeps the default). Only combinations that exist are offered.

`stb-pseudo --list-elements <name>` is the quickest way to check whether a
table covers what you need *before* committing to it.

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
stb-pseudo --list-banks                            # every PseudoDojo table
stb-pseudo --species Bi Se -p dojo-fr -o soc/      # full-relativistic, for spin-orbit
stb-pseudo -f structure.fdf -p dojo-sr3plus --fallback-dir dojo   # Ln3+ with 4f in the core
```

## Flag reference

```
stb-pseudo (-f/--file <structure.fdf> | --species EL [EL ...]) -p/--pp-path <BANK_OR_PATH>
           [--fallback-dir <BANK_OR_PATH>] [-o/--output <dir>] [--dry-run]
stb-pseudo --list-elements <BANK>
stb-pseudo --list-banks
```

- `-f/--file`: a SIESTA structure file (`.fdf`) to read the required
  species from. Mutually exclusive with `--species`.
- `--species`: element symbol(s) to resolve directly, instead of reading
  a structure file.
- `-p/--pp-path` (required unless `--list-elements` is given): a bundled
  bank name (`dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent]`,
  `virtual_vault`) or a directory of your own `.psf`/`.psml` files.
- `--fallback-dir` (optional): a second source, consulted only for
  elements missing from `-p/--pp-path`.
- `-o/--output` (optional, default `.`): directory the resolved
  pseudopotentials are copied into.
- `--dry-run` (optional): only run the availability check — copy
  nothing.
- `--list-elements <BANK>`: print the elements available in a bundled
  bank (any name `-p` accepts) and exit — no structure or `-p/--pp-path`
  needed.
- `--list-banks`: print every bundled table (relativity, functional,
  accuracy, element count) and exit.
- **Report** — always saved to `stb_pseudo_report.txt` in the current folder (`_1`, `_2`, ... when
  that name is taken); there is no option for it.

`references.bib` (SIESTA plus a citation for every bank actually used) is
always written — there's no flag for it.

Run `stb-pseudo --help` for the full list of options.

## What's next

Point `stb-inputfile`'s (example `1.1`) `-p/--pp-path`, or your
`calc.fdf`'s pseudopotential setup directly, at the output directory
`stb-pseudo` just filled — every species it resolved is already sitting
there with the right filename.
