# 3.10 — Fatbands Analyzer (`stb-fatbands`)

## What this tool does

`stb-fatbands` combines a SIESTA `.bands` file (k-path eigenvalues) with a
`.WFSX` file (the wavefunction coefficients for that same k-path) to color
and size every `(k, band)` point on a band-structure plot by its **orbital
character** — which atomic orbital(s), atom(s), or chemical species build
up a given state. Optionally cross-checked/weighted with a `.HSX` file
(Hamiltonian + overlap) for physically exact weights.

Five projection modes (`--projection`):

- `l` — angular momentum only (`s`, `p`, `d`, `f`)
- `ml` — full orbital detail (`s`, `px`, `py`, `pz`, `dxy`, ...)
- `atom` — per atom index
- `species` — per chemical species
- `species_l` (**default**) — species AND angular momentum combined
  (`Sn-s`, `Sn-p`, `O-s`, `O-p`, ...)

For a spin-polarized calculation, every category is automatically split
into `<category>_up`/`<category>_down` — see the theory section below for
why that matters.

## Why this matters (a bit of theory)

### LCAO wavefunctions: what a `.WFSX` file actually holds

SIESTA is a **linear-combination-of-atomic-orbitals (LCAO)** code: every
Bloch state is expanded as

```
psi_n(k, r) = sum_mu  c_mu,n(k) * phi_mu(r)
```

where `phi_mu` are the localized, numerical atomic orbitals of the basis
set (one or more per atom, depending on `PAO.BasisSize`: `s`/`p`/`d`/...,
possibly several per angular momentum for a multiple-zeta basis like DZP),
and `c_mu,n(k)` are complex expansion coefficients — one per `(orbital,
band, k-point)` triple. **These coefficients are exactly what a `.WFSX`
file stores.** A plain `.bands` file only has the eigenvalues (WHERE a
state sits in energy); the `.WFSX` additionally has WHAT that state is
made of.

### The projection weight: a Mulliken-like population analysis

Because the numerical atomic orbitals overlap each other in space, the
basis is **non-orthogonal**: `S_mu,nu = <phi_mu|phi_nu>` is not the
identity matrix (unlike, e.g., a plane-wave basis). The physically correct
weight of orbital `mu` in state `n` at k is not simply `|c_mu,n(k)|^2` but

```
w_mu,n(k) = Re[ c*_mu,n(k) * sum_nu S_mu,nu * c_nu,n(k) ]
```

— the same orbital-resolved decomposition of `<psi|S|psi>` a **Mulliken
population analysis** uses (`stb-fatbands` gets this via sisl's own
`norm2(projection="hadamard")`, applied to `psi` and the real SIESTA
overlap matrix `S`). Summed over every orbital of a category (a species,
an `l` shell, ...), it tells you what *fraction of that state's own norm*
lives in that category — and, up to the `l > 3` exclusion noted below,
summing across **every** category at a fixed `(n, k)` recovers the full
norm (≈1): the categories are a complete partition of the state, not an
arbitrary rescaling.

Computing `S` needs a real Hamiltonian object, read from `<label>.HSX`
(written when `SaveHS T` is set) — **without it**, sisl falls back to an
*implicit orthogonal basis* (`S = identity`), giving the cruder `|c_mu|^2`
instead. Same order of magnitude, not exact — `stb-fatbands` always prints
(and reports) which of the two it actually used, so you know exactly how
much to trust the numbers.

### Why `WFS.Write.For.Bands T` is required in the SIESTA calculation

**This is the one setting this tool cannot work without.** SIESTA computes
the same diagonalization for a `.bands` k-path regardless, but it only
*writes out* the wavefunction coefficients if explicitly told to — the
file can be large (one complex number per orbital, per band, per k-point),
so it is never saved by default. Concretely, in the `.fdf` (see
`spin/Ospin.fdf` in this folder for the complete, minimal, working
example):

```
WFS.Write.For.Bands T
BandLinesScale ReciprocalLatticeVectors
%block BandLines
  1   0.000  0.000  0.000  GAMMA
  5   0.500  0.000  0.000  X
%endblock BandLines
```

Add `SaveHS T` too, so `stb-fatbands` can use the accurate, overlap-aware
weights instead of falling back to the approximation above. Without
`WFS.Write.For.Bands T`, SIESTA still writes `<label>.bands` (the plain
eigenvalues) but **no** `<label>.bands.WFSX` at all — `stb-fatbands`
refuses to run outright in that case (a clear `No .WFSX file found` error)
rather than silently producing nothing or a wrong plot. If you already
have a finished calculation without this flag, you normally don't need a
fresh SCF cycle: add the flag, set `DM.UseSaveDM T`, and re-run — SIESTA
restarts from the converged density matrix and only needs to redo the
band-path diagonalization.

### The five projection modes, and why `species_l` is the default

- `l` and `ml` tell you *which orbital character* dominates a band, but
  not *which atom/species* it belongs to in a multi-species structure.
- `species` (and `atom`) tell you *which chemical element* (or specific
  site) dominates, but not *which orbital character* within it.
- `species_l` gives you both at once (`Sn-s`, `O-p`, ...) — the single
  most informative view for understanding chemical bonding character
  (e.g. "the valence band top is mostly O-p, the conduction band bottom is
  mostly Sn-s") without needing two separate runs. This is why it is the
  default as of v2.1.0 (it used to be plain `l`).
- `atom` is still useful on its own for comparing symmetry-inequivalent
  sites of the *same* species (e.g. surface vs. bulk atoms in a slab) —
  `species`/`species_l` cannot distinguish those, since both group by
  species/orbital, never by individual atom index.

### Spin-polarized calculations: every category is split automatically

For a spin-polarized (`nspin=2`) calculation, the spin-up (majority) and
spin-down (minority) channels can have genuinely different dispersion and
orbital character — not just a rigid energy shift, since that difference
*is* what magnetism means at the electronic-structure level. Merging both
channels into one category would silently combine two physically distinct
band sets into one indistinguishable series. `stb-fatbands` splits every
category into `<category>_up`/`<category>_down` for `nspin=2` (same
`_up`/`_down` naming `stb-dos` already uses for spin-resolved PDOS) —
demonstrated live in this walkthrough on a real, converged spin-polarized
calculation (see `spin/` below), where verifying this exact behavior
caught a genuine bug: the original implementation merged both spin
channels into one series, and two band sets whose conduction-band minima
differ by ~29 eV were being plotted (and averaged into the report's
per-category weight statistics) as if they were one.

### The band gap analysis "for free"

Resolving `--shift vbm`/`--shift cbm` already needs the exact same
VBM/CBM search `stb-bands` performs on the same `.bands` file — so
`stb-fatbands` reports it too (`[2] BAND GAP ANALYSIS`, reusing
`core/siesta_bands.py`'s `cbm_vbm`, spin-resolved for `nspin=2`), at no
extra computational cost and guaranteed to agree with what `stb-bands`
itself would report for the identical file.

### The `.bands`/`.WFSX` correspondence check

There is no file metadata tying a specific `.WFSX` to a specific `.bands`
file — nothing stops you from accidentally pointing `--wfsx` at the wrong
run. `stb-fatbands` guards this two ways: a hard k-point-count and
orbital-count match, and an eigenvalue cross-check (`--k-tol`, default
0.001 eV) at three sampled k-points (first/middle/last) between what the
`.WFSX` itself stored and what the `.bands` file reports. Always prefer
`--label` (auto-pairs same-run files) over hand-picking `--file`/`--wfsx`
from different calculations.

## The report: console output, `--save-report`, `--save-gnuplot`, `--view`

Every run prints a numbered report, the same `[0]...[6]` style every newer
tool in this suite uses:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | input files, shift mode, projection, category filter, output dir, active options |
| `[1] INPUT DATA` | Fermi energy, spin channels, bands x k-points, orbital count, weight accuracy (accurate/approximate), `.bands`/`.WFSX` correspondence check |
| `[2] BAND GAP ANALYSIS (k-path)` | VBM/CBM, direct/indirect gap, gap type — spin-resolved (+ half-metallic flag) for `nspin=2` |
| `[3] ORBITAL PROJECTION` | categories found/used, `l > 3` exclusions (if any), spin-split note, and a per-category table (orbital count, mean/max weight) |
| `[4] WRITING OUTPUT FILES` | whether `fatbands_<category>.dat`/`.gplot` were written |
| `[5] REFERENCES` | writes `references.bib` (SIESTA) |
| `[6] SUMMARY & FILES` | status and a recap of every file written |

- **`--save-report`** — persist the (text) report to
  `stb_fatbands_report.txt`. Off by default.
- **`--save-gnuplot`** — write `fatbands_<category>.dat` (one file per
  category — `k_position`/`energy(eV)`/`value` columns) + `fatbands.gplot`
  (all categories overlaid) together. Off by default — this tool used to
  write both **unconditionally** on every run; that's no longer the case.
- **`--view`** — show an interactive matplotlib preview (single-category
  runs get a colorbar-encoded scatter; multi-category runs get one
  discrete-colored series per category, shared legend). Off by default —
  this tool used to always show it, with no way to skip a blocking window.
  Marker sizes are deliberately smaller than `stb-ipr`/`stb-stm`'s own
  shared defaults (`core/band_scatter.py`): a fatbands plot scatters every
  `(k, band[, category])` point at once, several times denser than those
  tools' single-series plots.

## When you'd reach for it

- Identifying which orbital/species/atom character dominates a valence or
  conduction band, e.g. before discussing chemical bonding or assigning a
  band's character in a paper.
- Distinguishing majority/minority spin channels visually in a magnetic
  material's band structure.
- Comparing symmetry-inequivalent atoms of the same species (`--projection
  atom`), e.g. surface vs. bulk sites in a slab.
- A quick sanity check that a `.WFSX` and `.bands` pair truly belong to the
  same calculation, before trusting a plot built from hand-picked files.

## Two ways to run it

A — direct CLI:
```bash
stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.10
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.10.sh` proves this directly at the end. The menu
additionally defaults `--shift` to Fermi level and `--projection` to
`species_l` (just press Enter through both), plus separate `y/N` prompts
for `--save-report`/`--save-gnuplot`/`--view`.

## Files in this folder

- `Sn3O4.bands` / `Sn3O4.bands.WFSX` / `Sn3O4.HSX` / `calc.fdf` /
  `structure.fdf` / `Sn.ion(.xml)` / `O.ion(.xml)` — a real, converged bulk
  Sn3O4 calculation (copied from `test/3-analysis/10-fatbands/`, itself a
  deliberately SHORT band path — Gamma-X-Y-Gamma, 16 k-points — re-run from
  the same converged density matrix used elsewhere in this suite's tests,
  purely so `.bands.WFSX` stays a few MB instead of the 100+ MB a full,
  finely-sampled path would produce). Real dispersion, a genuine indirect
  gap, and two chemical species — used for most of this walkthrough.
  `calc.fdf`'s own `### BANDS` section shows `WFS.Write.For.Bands T` in a
  real, full production input (alongside DOS/relaxation settings you don't
  need for this tool specifically).
- `spin/Ospin.fdf` / `spin/Ospin.bands` / `spin/Ospin.bands.WFSX` /
  `spin/Ospin.HSX` / `spin/O.ion(.xml)` — a real SIESTA calculation run
  specifically for this walkthrough: a single, spin-polarized oxygen atom
  in a large vacuum box, seeded (`%block DM.InitSpin`) to converge to its
  physical 2 Bohr-magneton triplet ground state (confirmed in the
  calculation's own log: `spin moment ... |S| = 2.00000`). Deliberately
  the smallest possible real system that is genuinely spin-polarized —
  `Ospin.fdf` is short enough to read top-to-bottom as a complete,
  self-contained example of every fdf setting this tool needs
  (`WFS.Write.For.Bands T`, `SaveHS T`, `Spin polarized`).
- `example_3.10.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/10-fatbands/test.sh` for that).
- `.gitignore` — excludes `output/` and every file this tool itself
  generates on a run.

## Running the walkthrough

```bash
cd examples/3.10-stb-fatbands
./example_3.10.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder                    | What it shows                                                              |
|----------------------------|-----------------------------------------------------------------------------|
| `species-l-default/`       | The default `species_l` projection on a real, dispersive, multi-species structure; accurate (`.HSX`-based) overlap weights |
| `projection-modes/`        | The same structure/bands, contrasting `l` vs. `species` vs. `species_l` categories |
| `category-filter/`         | `--category` restricting the output to specific species/orbital combinations |
| `accuracy-fallback/`       | No `.HSX` available: the approximate (`|c|^2`) fallback, via `--geometry-file` |
| `spin-polarized/`          | The real spin-polarized fixture: spin-resolved gap analysis AND `<category>_up`/`<category>_down` splitting |
| `full-report/`             | Default (no report/data files) vs. `--save-report --save-gnuplot`, `references.bib` |

## Try it yourself

```bash
# Your own finished SIESTA calculation (needs WFS.Write.For.Bands T -- see above)
stb-fatbands --label my_calc --shift fermi --save-report --save-gnuplot

# Which orbital character dominates each band, split by species AND l
stb-fatbands --label my_calc --shift fermi --projection species_l --category Fe-d O-p

# A spin-polarized calculation -- categories split into _up/_down automatically
stb-fatbands --label my_magnetic_calc --shift fermi --view
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `--label`           | SIESTA label; shorthand for `--file <label>.bands` + auto-detected `.bands.WFSX`/`.HSX`/`.fdf`. |
| `--file`/`--wfsx`   | Explicit `.bands`/`.WFSX` paths (alternative to `--label`).             |
| `--hsx-file`        | Explicit `.HSX` for accurate, overlap-weighted norms.                  |
| `--geometry-file`   | Explicit `.fdf` (+ `.ion`/`.ion.xml`) fallback parent when no `.HSX` is available. |
| `--shift`           | Energy reference: `vbm`/`cbm`/`fermi`/`manual` (same vocabulary as `stb-bands`). |
| `--manual-value`    | Custom shift value (required with `--shift manual`).                   |
| `--gap-tol`         | Metallic-classification threshold, eV (default `0.01`).                |
| `--projection`      | `l` / `ml` / `atom` / `species` / `species_l` (default).               |
| `--category`        | Restrict output to specific category values (default: all found).      |
| `--k-tol`           | `.bands`/`.WFSX` eigenvalue cross-check tolerance, eV (default `0.001`).|
| `-o/--output-dir`   | Where all generated files (and `references.bib`) land.                 |
| `--save-report`     | Persist the full report to `stb_fatbands_report.txt`.                  |
| `--save-gnuplot`    | Also write `fatbands_<category>.dat`/`fatbands.gplot`.                 |
| `--view`            | Show an interactive matplotlib preview.                                |

## What's next

`core/band_scatter.py`'s scatter-plot machinery is shared with `stb-ipr`
(a single continuous inverse-participation-ratio value per `(k, band)`,
no category) and `stb-spintexture` (signed spin-expectation components).
`core/siesta_wfsx.py`'s `.WFSX`/`.HSX` loading is shared with
`stb-wfdensity`, `stb-sts`, and `stb-coop` — tools that need the same
wavefunction coefficients for a full k-mesh instead of a single k-path.
`stb-bands` covers the plain (non-orbital-resolved) band gap question this
tool's own `[2] BAND GAP ANALYSIS` section reuses.
