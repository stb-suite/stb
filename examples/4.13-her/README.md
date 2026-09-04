# 4.13 — Workflow: HER, Hydrogen Evolution Reaction (`stb-her` / `stb-herRefs` / `stb-herAnalysis`)

Hydrogen production by water electrolysis needs a catalyst that binds
adsorbed hydrogen atoms (H\*) neither too weakly nor too strongly. The
computational hydrogen electrode (CHE, Nørskov et al., J. Phys. Chem. B
108, 17886, 2004) turns that question into a single first-principles
number, `Delta-G_H*`, computable from ground-state DFT total energies
alone — no explicit electrochemical interface, no applied potential, no
reaction-barrier search. This workflow automates the whole pipeline:
finding every symmetrically distinct adsorption site, building the
reference calculations the formula needs (gas-phase H2, a Boys-Bernardi
BSSE counterpoise triad, a vibrational zero-point/entropy correction), and
combining them into `Delta-G_H*` with the numerical safeguards a
by-hand calculation is easy to get subtly wrong on.

Three stages, like `4.7-hubbardu`: site search (Stage 1), reference/ZPE
prep (Stage 2), analysis (Stage 3). None of the three run SIESTA for you —
each is a real, no-shortcuts folder-generation/analysis tool; you run
SIESTA yourself in between. This walkthrough proves the whole chain end to
end with a synthetic dataset whose exact `Delta-G_H*` is known in advance,
so you can see precisely what "correct" looks like before trusting it on
your own material.

## 1. Theory

### 1.1 The CHE descriptor and the Sabatier principle

```
Delta-G_H* = Delta-E_H (BSSE-corrected) + Delta-ZPE - T*Delta-S
```

`Delta-E_H` is the electronic energy of moving one H atom from the
gas-phase H2 reservoir onto the surface:

```
Delta-E_H = E(slab+H) - E(slab) - 1/2 * E(H2)
```

(the `1/2` because H2 supplies two H atoms per molecule). `Delta-ZPE` and
`Delta-TS` add the vibrational zero-point energy and entropy difference
between the adsorbed and gas-phase states, at temperature `T`. The
qualitative reading (`stb-herAnalysis`'s own `[4] FINAL RESULT` verdict):

- `Delta-G_H* ~ 0` (within ~0.2 eV): near-optimal — the **Sabatier
  principle**, the same volcano-plot logic behind essentially every
  heterogeneous-catalysis screening study. H adsorbs and desorbs at
  comparable, low barriers.
- `Delta-G_H* >> 0`: too weak — H\* adsorption itself is the bottleneck;
  the surface barely wants to bind H at all.
- `Delta-G_H* << 0`: too strong — H\* binds so tightly that recombination
  into H2 and desorption becomes the bottleneck instead ("H poisoning").

### 1.2 Why three stages

`Delta-G_H*` needs **seven** separate SIESTA energies, computed
consistently:

1. `E(slab+H)` at the best (lowest-energy) adsorption site.
2. `E(slab)`, the bare substrate, same numerical settings.
3. `E(H2)`, the isolated gas-phase molecule, relaxed to its own bond length.
4. Three more for the Boys-Bernardi BSSE counterpoise correction (Section 1.3).
5. A partial-Hessian force/energy set for `Delta-ZPE`/`Delta-TS` (Section 1.4).

None of these exist until you know **which** site is the winning one —
Stage 1 (`stb-her`) writes one relaxation folder per symmetrically distinct
candidate site; Stage 2 (`stb-herRefs`) reads back whichever one SIESTA
actually converged to the lowest energy and builds every downstream
reference **from that site's own relaxed geometry**; Stage 3
(`stb-herAnalysis`) combines everything.

### 1.3 Basis Set Superposition Error (BSSE)

Same issue `4.8-adsorption`'s own README explains in more depth: because
SIESTA is a localized-basis (LCAO) code, the naive `Delta-E_H` above is
contaminated by BSSE — in the COMBINED `slab+H` calculation, each fragment
"borrows" extra basis functions that belong to the other fragment,
artificially lowering the combined energy and making adsorption look more
strongly bound than it really is. The Boys & Bernardi (1970) fix
re-evaluates each fragment WITH the other fragment's basis functions
present as chargeless "ghost" atoms, at the site's own relaxed geometry:

```
BSSE(slab)  = E(slab, real)      - E(slab, real + H-ghost)     = E_deformed - E_ghost
BSSE(H)     = E(H, real)         - E(H, real + slab-ghost)     = E_H_iso - E_H_ghost
BSSE(total) = BSSE(slab) + BSSE(H)
Delta-E_H (corrected) = Delta-E_H (raw) + BSSE(total)
```

`03_slab_deformed` (the winning site minus H, real atoms only) and
`04_slab_ghost` (same geometry, H turned into a ghost) isolate the slab
side; `06_h_ghost_slab` (H real, slab ghosted) and `07_h_isolated` (H
alone) isolate the H side — all four evaluated at the **exact same,
already-relaxed geometry**, single-point (no further relaxation).

### 1.4 Zero-point energy and entropy: local vs. full mode

H is by far the lightest atom in the system, so its vibrational
zero-point energy and entropy are the only ones that matter for a useful
correction — but *how much* of the substrate's own vibrational response
to get right is a real accuracy/cost trade-off, controlled by
`--zpe-mode`:

- **`local`** (default): a **partial-Hessian, decoupled-oscillator**
  approximation. Only the adsorbed H is displaced (+/-x, +/-y, +/-z,
  6 single-point folders, `05_zpe_calc/disp_00{1..6}/`), every substrate
  atom held fixed. Diagonalizing the resulting 3x3 mass-weighted Hessian
  gives 3 harmonic frequencies; each contributes
  `ZPE = 1/2 * h*f` and a Bose-Einstein entropy term. Fast, but ignores
  H-substrate vibrational coupling.
- **`full`**: complete Phonopy displacement sets for BOTH the winning site
  AND the clean slab (`05_zpe_calc_site/`, `05_zpe_calc_clean/`), reusing
  the exact same FORCE_SETS-building code the Raman/IR/Phonons tools
  already use. Because a full-structure phonon calculation's ZPE includes
  the ENTIRE slab's own zero-point energy (perturbed by H adsorption, not
  just H\*'s own contribution), the physically correct expression
  subtracts the clean slab's own ZPE/entropy first:
  `Delta-ZPE = [ZPE(site) - ZPE(clean_slab)] - 1/2 * ZPE(H2)` (and
  likewise for `Delta-S`) — roughly double the cost of one full phonon run.
- **`standard`**: no folders generated at all — `stb-herAnalysis` uses a
  fixed literature offset, `Delta-ZPE - T*Delta-S ~= +0.24 eV` (the
  Nørskov et al. value), when you don't need (or can't afford) either mode
  above.

H2's own ZPE/entropy are **fixed literature values** in every mode
(`H2_ZPE_EV = 0.270`, `H2_TS_298K_EV = 0.400`, both at 298.15 K regardless
of `--temp`) — H2's vibrational structure (one stretch mode, well-known
rotational structure) is experimentally well characterized, and
re-deriving it from a DFT harmonic calculation would introduce more
xc-functional noise than it removes.

### 1.5 `config_extra.fdf`: forcing directives without touching your template

Every folder this workflow writes gets a `config_extra.fdf` sidecar,
`%include`-d at the very top of `calc.fdf` (SIESTA's fdf reader is
first-occurrence-wins, so this always overrides whatever the same tag
might say later in your own template) — the same mechanism `4.8-adsorption`,
`4.11-raman`, and `4.12-ir` use, never editing your `--calc` text in
place. Three directives are **mandatory, not opt-in**, everywhere a slab
is present:

- `MD.VariableCell false` — a local adsorption-site relaxation has no
  business changing the bulk-derived in-plane lattice constant.
- `Slab.DipoleCorrection .true.` — adsorbing H on only ONE face breaks
  whatever mirror/inversion symmetry the clean slab had along the surface
  normal, giving the cell a net dipole along a PERIODIC direction; without
  the correction, that spurious periodic-image field contaminates every
  site's energy (and therefore the site ranking Stage 2 does).
- `DFTD3 .true.` — standard GGA misses dispersion forces entirely, which
  can matter for a physisorption-like distant site; mandatory rather than
  merely advised, since a dispersion correction applied to one side of an
  energy difference but not the other would bias it.

`Spin polarized` is *also* mandatory almost everywhere (H adsorption
intermediates can have real radical/open-shell character), with exactly
**one** deliberate exception (Section 1.6). Every single-point derived
folder Stage 2 writes (`00_clean_slab`, `03_slab_deformed`, the BSSE
triad, every ZPE displacement folder) inherits the **same** four
directives from the winning site — the "level-of-theory propagation"
convention documented across this whole suite: comparing energies
computed at inconsistent levels of theory has repeatedly produced
physically-impossible results elsewhere in this codebase (a BSSE
correction "more negative than uncorrected", for one).

### 1.6 Why H2 alone is forced spin-UNpolarized

H2's ground state is an unambiguous **closed-shell singlet** — two
electrons filling the bonding sigma orbital, no radical character at all,
unlike the adsorbed site (which genuinely can be open-shell). Forcing
`Spin polarized` there anyway isn't just unnecessary, it's actively
risky: a spin-polarized SCF can converge to a spurious nonzero magnetic
moment depending on the initial guess, silently biasing `1/2 * E(H2)` in
every `Delta-E_H` computed against it. `02_h2_molecule/config_extra.fdf`
is the one folder in this entire workflow with `Spin non-polarized`
instead of `Spin polarized` — forced explicitly (not just left to
whatever your template happens to say), so the choice is
self-documenting and immune to your own `--calc` template's own setting.

### 1.7 Fragment labeling: never confusing the adsorbate with the slab

Every atom Stage 1 writes is relabeled by which fragment it came from —
`<symbol>_slab` for the substrate, `<symbol>_ads` for the adsorbed H (a
`fragment_manifest.json` sidecar records the split). This matters for a
concrete, easy-to-hit failure mode: `write_fdf` always groups atoms by
species, not physical origin, so if your OWN slab already carries H (a
passivated dangling-bond termination, say), a bare `H` label for the
freshly-adsorbed atom would be silently indistinguishable from the slab's
own H after a `structure.fdf` round trip. Stage 2's BSSE ghost-triad
construction (`make_ghost_variant`) and the isolated-H reference
(`isolate_atom`) both rely on this label surviving intact.

### 1.8 What this workflow simplifies — explicitly

By hand: writing 8+ nearly-identical `calc.fdf` variants (single-point vs.
relaxed, Gamma-only vs. your own k-grid, spin-polarized vs. not);
remembering BSSE's Boys-Bernardi ghost-atom convention and getting the
`+`/`-` signs right; deriving harmonic ZPE/entropy from finite-difference
forces (Section 5.2); tracking which site actually won *before* building
any reference; and finally combining seven numbers into one formula
without a sign error. This workflow does all of that mechanically and
consistently — and, just as importantly, tells you **when a residual force
is a real problem vs. an expected artifact of the method itself** (Section
5.3), rather than either staying silent or crying wolf on every folder
uniformly.

## 2. Libraries and external dependencies used

- **`pymatgen` / `spglib`** — `AdsorbateSiteFinder` (site search, Stage 1)
  and `SpacegroupAnalyzer` (the `[1] SLAB SYMMETRY` section, space/point
  group, Wyckoff sites, Section 3.2).
- **`phonopy`** — `--zpe-mode full` only (Section 1.4), the same
  FORCE_SETS-building code the Raman/IR/Phonons tools already use.
- **`numpy`** — the local-mode Hessian diagonalization and Bose-Einstein
  entropy formulas (Section 5).
- **`matplotlib`** — `stb-herAnalysis --show`, an optional on-screen
  preview of the same energy-breakdown chart the gnuplot files draw.
- **`gnuplot`** (optional, `stb-herAnalysis --plot`) — writes
  `plot/her_energy_breakdown.{dat,gplot}`, this workflow category's
  bar-chart convention (`4.7-hubbardu`'s own `--save-gnuplot`, adapted
  here to a categorical breakdown chart instead of a response curve).
- No network access, no ML/MACE dependency anywhere in this workflow —
  every stage is a pure, local, deterministic file operation.

## 3. Stage 1: adsorption site search (`stb-her`, code `4.13.1`)

### 3.1 What it does

Given a slab/2D `structure.fdf` (vacuum along exactly one axis) and a
`calc.fdf` template, `stb-her` finds every symmetrically distinct
H-adsorption site of the requested type(s) (`ontop`/`bridge`/`hollow`/
`all`) and writes one `sites/site_N_<type>/` folder per site — **always
every distinct site**, unlike `stb-adsorb`'s own `--site-index`/
`--all-sites` choice: HER's whole point is finding the *global* most
stable site, so there's no reason to write only some of them.

### 3.2 Symmetry, live

```
$ stb-her -s structure.fdf -c calc.fdf --no-intro
```

```
[1] SLAB SYMMETRY
------------------------------------------------------------
Property       | Value
--------------------------------
Space Group    | P6/mmm (191)
Point Group    | 6/mmm
Crystal System | Hexagonal
Hall Symbol    | -P 6 2
Layer Group    | p6/mmm (No. 80)

  1 symmetrically distinct substrate atom(s) (space group P6/mmm (No. 191)):
Atom idx (0-based) | Element | Wyckoff | Multiplicity
-----------------------------------------------------
0                  | C       | d       | 2
```

This is the **substrate's own** crystallographic symmetry (space/point
group, layer group since a slab/2D structure genuinely has one, Wyckoff
sites) — a real space-group/Wyckoff analysis via `core/symmetry.py`, going
further than `stb-adsorb`'s own Stage 1, which never reports a group label
at all (only a raw-vs-reduced candidate *count*, see Section 3.3). It
answers "how symmetric is my starting material", not directly "how many H
sites are there" — that's Section 3.3's job.

### 3.3 Site search: raw-vs-symmetry-reduced counts

```
[3] ADSORPTION SITES: FINDING & COUNT
------------------------------------------------------------
Site type | Raw candidates | After symm. reduction | Reduced away (equivalent)
------------------------------------------------------------------------------
ontop     | 2              | 1                      | 1
bridge    | 6              | 3                      | 3
hollow    | 0              | 0                      | 0
TOTAL     | 8              | 4                      | 4
```

`pymatgen.AdsorbateSiteFinder` proposes every geometric candidate first
(`Raw candidates`), then collapses symmetry-equivalent ones (`--symprec`,
default `0.01`) down to the physically distinct set actually written —
4 sites total for pristine graphene (1 ontop, 3 bridge, 0 hollow; this
2-atom honeycomb basis never exposes a hollow candidate to pymatgen's
finder, at any cell size). A numbered `sites/adsorption_sites.png` and a
reusable `sites/site_positions.dat` (Section 3.5) are written alongside.

### 3.4 Fragment labels and `config_extra.fdf`, live

```
%block ChemicalSpeciesLabel
 1   6   C_slab
 2   1   H_ads
%endblock ChemicalSpeciesLabel
```

```
$ cat sites/site_1_ontop/config_extra.fdf
```
```
# Auto-generated -- keeps the cell fixed ...
MD.VariableCell false
# Auto-generated -- forces the slab dipole correction ...
Slab.DipoleCorrection      .true.
# Auto-generated -- forces spin-polarized SCF ...
Spin                polarized
# Auto-generated -- forces the Grimme DFT-D3 dispersion correction ...
DFTD3                   .true.
```

See Sections 1.5–1.7 for why each of these is mandatory.

### 3.5 Manual overrides: `--position` and `--positions-file`

For a site pymatgen's automatic finder doesn't propose (or to reproduce
one exact point), `--position X Y` bypasses site-finding entirely — the
height is measured along the **true surface normal** (`finder.mvec`), not
a naive Cartesian z-offset, which matters the moment your cell's in-plane
lattice vectors have any out-of-plane Cartesian component (a tilted or
database-fetched structure that wasn't pre-aligned).

Every run — automatic or manual — always writes `sites/site_positions.dat`
(in-plane fractional coordinates, one line per site, human-editable), which
`--positions-file` reads straight back: hand-delete the sites you don't
want, or reuse an exact site set across multiple `--height` sweeps,
without re-running the site finder.

### 3.6 `--both-sides`

For a genuinely free-standing 2D material (vacuum on both faces, like
this walkthrough's graphene), `--both-sides` mirrors each selected site
onto the bottom face too — needs a symmetry operation actually mapping top
to bottom (a centrosymmetric or free-standing structure; a one-sided slab
fails with a clear error instead of guessing).

### 3.7 The too-close-atoms safety check

```
[WARNING] site_1_ontop: closest slab-H distance is only 0.100 Ang --
likely overlapping atoms (--height too small, or a bad
--position/--positions-file entry). Check before running SIESTA.
```

A sanity check on `--height` (or a hand-picked `--position`), not a
physical statement — no real X-H bond is anywhere near 0.1 Ang.

### 3.8 Report structure

```
[0] RUN METADATA
[1] SLAB SYMMETRY
[2] CLEAN SLAB REFERENCE
[3] ADSORPTION SITES: FINDING & COUNT
[4] WRITING SITE FOLDERS
[5] SUMMARY & NEXT STEPS
[6] LIBRARY WARNINGS
```

`[2]` also recenters the slab in its vacuum gap if it starts near a cell
boundary (avoids an atom drifting across the periodic image during
relaxation and reading as a spurious huge displacement later) — reported,
not silent, either way. `[6]` traps pymatgen/spglib's own deprecation
noise instead of letting it interleave mid-report.

### 3.9 Running it both ways

**A — direct CLI**: every command above.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.13.1
```

`4.13.1` asks for the structure file, the calc.fdf template, a
pseudopotential source (bundled bank / custom path / skip), the site
type, height, both-sides, output directory, then optional
symprec/vacuum-gap overrides, and finally whether to save a report — the
same flags as the CLI, in the same order. `example_4.13.sh` proves both
paths agree.

## 4. Stage 2: references & ZPE prep (`stb-herRefs`, code `4.13.2`)

### 4.1 What it does

Scans Stage 1's `sites/site_*/` for the lowest-`FreeEng` **finished**
site, reads back its relaxed geometry (from the SIESTA `.out` file's own
`outcoor: Relaxed atomic coordinates` block — not a separate `.XV` file),
then writes (Section 1.2–1.4): `00_clean_slab/`, `02_h2_molecule/`, the
diagnostic `03_slab_deformed/`, the BSSE triad
(`04_slab_ghost/`/`06_h_ghost_slab/`/`07_h_isolated/`), and the
`--zpe-mode`-dependent ZPE folder(s).

### 4.2 The reference-folder table, live

```
Folder           | Atoms | Formula        | Run type
-------------------------------------------------------------------
00_clean_slab    | 2     | C2             | single-point
02_h2_molecule   | 2     | H2             | CG relax (Gamma-only)
03_slab_deformed | 2     | C2             | single-point (diagnostic)
04_slab_ghost    | 3     | C2 +H(ghost)   | single-point (BSSE)
06_h_ghost_slab  | 3     | H +C2(ghost)   | single-point (BSSE)
07_h_isolated    | 1     | H              | single-point (BSSE)
```

Ghost atoms are called out explicitly in the formula (`+H(ghost)`) —
real element, real pseudopotential file, zero valence charge, contributing
basis functions only. `02_h2_molecule` is the one folder that **relaxes**
(`MD.TypeOfRun CG`, its own bond length must reach its own equilibrium,
unlike every other folder here which is evaluated at a fixed, already-known
geometry) and the one folder forced Gamma-only + spin-**unpolarized**
(Section 1.6) instead of inheriting the winning site's own `kgrid`/`Spin`.

### 4.3 `--zpe-mode local|full|standard`

See Section 1.4 for the physics. `local` (default) writes
`05_zpe_calc/disp_00{1..6}/` — 6 single-point folders, H displaced
+/-x/y/z from its relaxed position, everything else fixed. `full` writes
complete Phonopy displacement sets under `05_zpe_calc_site/` AND
`05_zpe_calc_clean/`, with the raw-candidate-vs-symmetry-reduced
displacement count reported for each (same spirit as Stage 1's own
site-count table):

```
Site supercell  : 1x1x1 (19 atom(s)) -- 114 raw candidate displacement(s)
[OK] 114 displacement folder(s) under 05_zpe_calc_site/ (symmetry-reduced
     from the 114 raw candidates above)
```

`standard` writes nothing here at all — Stage 3 uses the fixed Nørskov
offset instead.

### 4.4 Report structure

```
[0] RUN METADATA
[1] WINNING SITE
[2] REFERENCE FOLDERS
[3] ZPE PREPARATION
[4] SUMMARY & NEXT STEPS
```

`[1]` lists every scanned site's `FreeEng` (marking the winner), the
energy spread across all readable sites (a quick consistency check — a
huge spread across symmetry-equivalent-looking sites usually means one of
them didn't actually converge), and the same residual-force quality
diagnostic Stage 3 uses (Section 5.3) for the winning site itself.

### 4.5 Running it both ways

**A — direct CLI**: `stb-herRefs --directory her_study --zpe-mode local`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.13.2
```

`4.13.2` asks for the Stage-1 directory, the SIESTA output filename, a
pseudopotential source, the ZPE mode, then (if `local`/`full`) the
displacement size and (if `full`) the supercell, and finally whether to
save a report.

## 5. Stage 3: analysis (`stb-herAnalysis`, code `4.13.3`)

### 5.1 Report structure

```
[0] RUN METADATA
[1] ELECTRONIC ENERGIES (FreeEng)
[2] BSSE CORRECTION
[3] THERMAL CORRECTION
[4] FINAL RESULT
[5] SUMMARY & FILES
```

### 5.2 `[1]`/`[2]`: every energy and every BSSE term, in a table

Every raw `FreeEng` is printed inline (with its own quality diagnostic,
Section 5.3) AND recapped in a clean table at the end of `[1]`; `[2]`
prints the full BSSE derivation — deformation cost, both BSSE terms, the
total, and `Delta-E_H` raw vs. corrected — as one table instead of a wall
of separate lines, so every contribution is visually separated before it
ever reaches the final formula:

```
Term                                          | Value (eV)
----------------------------------------------------------
Deformation cost (diagnostic, not used below) | +0.1000
BSSE (slab side)                              | +0.0500
BSSE (H side)                                 | +0.0500
BSSE (total)                                  | +0.1000
Delta-E_H (raw)                               | -0.5000
Delta-E_H (corrected)                         | -0.4000
```

### 5.3 Quality diagnostics: a real warning vs. an EXPECTED artifact

Every energy is checked for SCF convergence and residual force against
`--force-tolerance` (default `0.05` eV/Ang) — **except** `03_slab_deformed`
and `04_slab_ghost`. Both deliberately evaluate the winning site's own
relaxed geometry with H removed/ghosted; the atoms nearest the former H
position are, by construction, never at their own equilibrium there (a
real, localized effect — only the 2-3 atoms closest to where H used to be
carry a large force, typically of order 1 eV/Ang, while every other atom
and the *net* force stay small). A generic "this geometry may not be
relaxed" warning would fire on **every single run**, for a reason that has
nothing to do with a bad calculation — so the check is skipped entirely
for exactly these two folders, not softened into a note that would still
always fire. Every other folder still gets the ordinary check:

```
[WARNING] Residual force on E_clean (2.5000 eV/Ang) exceeds
--force-tolerance (0.05 eV/Ang) -- this geometry may not be relaxed.
```

**Take this seriously everywhere it DOES appear** — it means exactly what
it says for every folder except the two named above.

### 5.4 `[4]`: the final result, itemized

```
Delta-G_H* = -0.1209 eV
Qualitative HER assessment: near-optimal (Sabatier principle, |Delta-G_H*| small)

Contribution          | Value (eV)
----------------------------------
Delta-E_H (raw)       | -0.5000
BSSE correction       | +0.1000
Delta-E_H (corrected) | -0.4000
Delta-ZPE             | +0.0810
-Delta-TS             | +0.1981
Delta-G_H* (TOTAL)    | -0.1209
```

Every term that sums to the final answer, top to bottom, colored by role
(a plain contribution, the `Delta-E_H (corrected)` running checkpoint, and
the final total) — the same three-way distinction the optional chart below
draws with color too.

### 5.5 `--plot`/`--no-plot` and `--show`/`--no-show`

Two independent yes/no questions, asked interactively unless answered by
a flag (safe for scripts/CI: an unattended run with no TTY on stdin
defaults to "no" for both rather than hanging):

- **`--plot`**: saves `plot/her_energy_breakdown.{dat,gplot}` — this
  workflow category's gnuplot convention (`4.7-hubbardu`'s own
  `--save-gnuplot`, WORKFLOW_TOOLS write gnuplot pairs, not matplotlib
  PNGs directly). `gnuplot her_energy_breakdown.gplot` renders a PDF bar
  chart of the exact same breakdown table above, colored the same way (a
  plain contribution, the running subtotal, the final total), with an
  explicit y-range padded to the data's own span so the rotated category
  labels always land in clear space below the zero line rather than
  running back up into the bars themselves.
- **`--show`**: opens the same chart directly on screen via matplotlib
  (needs a display) — a quick preview alongside, not instead of, the
  gnuplot files.

### 5.6 Running it both ways

**A — direct CLI**: `stb-herAnalysis --directory her_study --save-report`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.13.3
```

`4.13.3` asks for the Stage-1/2 directory, the SIESTA output filename,
temperature, force tolerance, an output-file basename, then whether to
save a report — the plot/show questions are asked by `stb-herAnalysis`
itself afterward, exactly as they would be from the direct CLI.

## 6. Worked example: a known analytic ground truth, end to end

Real SIESTA output isn't available inside this walkthrough (there's no
SIESTA binary to invoke), so — exactly like `4.7-hubbardu`'s own
linear-response verification — Section "output/workflow/" of the script
fabricates `calc.out`/`.FA` files with a **hand-chosen** set of numbers, so
the correct `Delta-G_H*` is known *before* running Stage 3:

```
E_clean (bare slab)        = -400.000000 eV
E_slab+H (winning site)    = -416.250000 eV   ->  Delta-E_H (raw) = -0.5000 eV
E_H2 (gas-phase reference) =  -31.500000 eV
BSSE (slab side + H side)  =   +0.100000 eV   ->  Delta-E_H (corrected) = -0.4000 eV
H* local Hessian           = isotropic harmonic spring, k = 5.0 eV/Ang^2
```

`k = 5.0 eV/Ang^2` gives an exact, hand-verifiable ZPE/entropy via the
standard harmonic-oscillator formulas — 3 degenerate modes at
`f = (1/2*pi) * sqrt(k/m_H)`, each contributing `1/2 * h*f` to the ZPE and
a Bose-Einstein occupation term to `T*S` at 298.15 K:

```
Delta-ZPE  = +0.0810 eV
Delta-TS   = -0.1981 eV
Delta-G_H* = -0.4000 + 0.0810 - (-0.1981) = -0.1209 eV  (near-optimal)
```

The script independently computes this in plain Python (not by calling
`stb-herAnalysis`), then runs the real 3-stage pipeline on the fabricated
data and asserts the two agree:

```
Delta-G_H* = -0.1209 eV
Qualitative HER assessment: near-optimal (Sabatier principle, |Delta-G_H*| small)
```

— confirming the entire pipeline (energy reading → BSSE arithmetic →
harmonic ZPE/entropy → final formula) end to end, the same way
`4.7-hubbardu`'s own synthetic linear-response dataset confirms *its*
fitting pipeline.

## 7. Known, deliberate limitations

- **No SIESTA run is ever performed by this suite.** All three stages
  generate/read input and output files; you run SIESTA in every `sites/
  site_*/` folder and every reference/ZPE folder Stage 2 writes, yourself,
  in between stages.
- **Rerunning SIESTA in the same folder does NOT continue from where it
  left off.** SIESTA reads atomic coordinates from `structure.fdf` itself
  (the ORIGINAL, pre-relaxation input), not automatically from a
  previously-written `.XV` file — verified live, the hard way, in this
  suite's own development: doubling `MD.Steps` and rerunning restarts the
  full relaxation from scratch (the SCF density matrix does get reused via
  SIESTA's own `DM.UseSaveDM`-style continuation, which is why a rerun can
  still converge faster/differently than the first attempt, but the
  GEOMETRY search itself restarts). If a site doesn't converge within
  `MD.Steps`, raise `MD.Steps` in your `--calc` template (or loosen
  `MD.MaxForceTol`) BEFORE the first run, not after.
- **Not every symmetrically distinct site converges equally easily.** In
  real testing, `ontop` sites on a hexagonal 2D material typically
  converge in well under 100 CG steps; `bridge`/`hollow` sites can sit
  near a much flatter part of the potential-energy surface and need
  several times that (or a looser `MD.MaxForceTol`, or a different
  optimizer such as `MD.TypeOfRun Broyden`) before the residual force
  actually drops below tolerance.
- **A single isolated atom cannot use many MPI ranks.** `07_h_isolated`
  (and any true single-atom folder) will make SIESTA abort
  ("You have too many processors for the system size") if you launch it
  with the same `mpirun -np N` you used for the multi-atom folders — run
  it with `-np 1` (or however few ranks SIESTA's own orbital-distribution
  check accepts for a 1-atom system).
- **A large `Delta-G_H*` is not automatically a bug.** A chemically inert,
  large-gap pristine 2D material can genuinely have a strongly unfavorable
  `Delta-G_H*` — that IS the physically correct answer for that surface
  (defect engineering or doping is the usual next step for such a
  material, not a re-run). Check the individual energies and BSSE/thermal
  terms (Sections 5.2–5.4) for internal consistency before assuming a
  surprising number is wrong.
- **`--zpe-mode local` ignores adsorbate-substrate vibrational coupling**
  by construction (Section 1.4) — a real quantitative limitation, not a
  bug, and the reason `full` mode exists for when that coupling matters.
- **The BSSE triad and `03_slab_deformed` inherit the winning site's own
  theory settings unconditionally** (Section 1.5) — there is currently no
  way to evaluate them at a *different* level of theory than the site
  itself, by design (the whole point is numerical consistency).

## 8. Step-by-step: running this workflow on your own structure

1. Have a relaxed slab/2D `structure.fdf` (vacuum along exactly one axis)
   and a working `calc.fdf` template.
2. Run Stage 1: `stb-her -s <structure.fdf> -c <calc.fdf> --site-type all
   [--both-sides] [-p <pseudo-bank-or-path>]`.
3. Run SIESTA yourself in every `her_study/sites/site_*/` folder. If one
   doesn't converge within `MD.Steps`, see Section 7's second bullet
   before just re-launching it.
4. Run Stage 2: `stb-herRefs --directory her_study --zpe-mode local`
   (or `full` if adsorbate-substrate coupling matters for your system, or
   `standard` for the fast literature-offset estimate).
5. Run SIESTA yourself in every folder Stage 2 just wrote — remember
   `07_h_isolated` needs far fewer MPI ranks than the rest (Section 7).
6. Run Stage 3: `stb-herAnalysis --directory her_study --save-report`.
   Read `[2]`/`[3]` for any `[WARNING]` before trusting the final number
   — Section 5.3 explains exactly which ones are expected and which
   aren't. Answer `y` to the plot prompt (or pass `--plot`) for a saved
   energy-breakdown chart.

## 9. Files in this folder

| File | Purpose |
|---|---|
| `structure.fdf` | Free-standing 2-atom graphene primitive cell, vacuum along c — small and fast, purely to exercise the tools; also enables the `--both-sides` demo (Section 3.6). |
| `calc.fdf` | Shared template (non-polarized, `MD.VariableCell true`) — proves `config_extra.fdf`'s forced overrides (Spin polarized, fixed cell) actually win over the template's own settings. |
| `example_4.13.sh` | The guided walkthrough (**not** an automated test — see `test/4-workflow/13-her/{prep,refs,analysis}/test.sh` for that). Pauses between sections so you can read before moving on; safe to re-run. |
| `output/` | Created by `example_4.13.sh` when you run it (git-ignored, not checked in). See below. |

## 10. Running the script

```bash
./example_4.13.sh
```

| Case | Command(s) | What it shows |
|---|---|---|
| `output/stage1/` | `stb-her -s structure.fdf -c calc.fdf` | Symmetry, raw-vs-reduced site table, fragment labels, `config_extra.fdf` |
| `output/stage1_both/` | `stb-her --site-type ontop --both-sides` | Both faces of a free-standing 2D material |
| `output/stage1_position/` | `stb-her --position 0.0 0.0 --height 1.6` | Manual site override, true-surface-normal height |
| `output/stage1_positionsfile/` | `stb-her --positions-file mypositions.dat` | Round-tripping `site_positions.dat` |
| `output/stage1_overlap/` | `stb-her --site-type ontop --height 0.1` | The too-close-atoms warning |
| `output/workflow/` | Stage 1 (single site) -> (fabricated relaxation/energies) -> Stage 2 -> Stage 3 | The full 3-stage chain, known-answer worked example (Section 6), the `[NOTE]`-skip fix (Section 5.3), and `--plot` |
| *(no folder — a diff only)* | Stage 1 via `printf … \| stb-suite` | Proof the interactive menu (`4.13.1`) agrees with the CLI |

## What's next

- **`4.8-adsorption`** — the general-purpose version of Stage 1/2's own
  BSSE and site-search machinery (any adsorbate, multi-adsorbate, ML
  pre-screening, height sweeps) — read this one first if you haven't, for
  a deeper BSSE derivation.
- **`4.9-neb`** — once you know `Delta-G_H*` is favorable, this is how
  you'd estimate the actual reaction *barrier* for H migrating between
  two adsorption sites, rather than just the thermodynamic endpoint this
  workflow computes.
- **`3.5-stb-symmetry`** — the full per-site/operations symmetry analysis
  this workflow's own `[1] SLAB SYMMETRY` section is a summarized slice of.
- Every other `4.x` workflow example in this suite generates real
  candidate/perturbed geometries and expects a real SIESTA run in
  between stages, the same two-way (direct CLI / interactive menu) split,
  and the same numbered-report convention — if this is your first
  workflow example, `4.1-strain` is the shortest two-stage one to start
  with.
