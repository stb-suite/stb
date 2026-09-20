# 4.14 — Workflow: OER, Oxygen Evolution Reaction (`stb-oer` / `stb-oerIntermediates` / `stb-oerRefs` / `stb-oerAnalysis`)

Water electrolysis has two half-reactions: hydrogen evolution (HER,
`4.13-her`) at the cathode, and oxygen evolution (OER, this workflow) at
the anode. OER is almost always the efficiency bottleneck — it moves 4
electrons instead of HER's 2, through three distinct adsorbed
intermediates (OH\*, O\*, OOH\*), and real catalysts routinely need
several hundred extra millivolts beyond the thermodynamic minimum to
drive it at a useful rate. The 4-electron computational hydrogen
electrode (CHE, Rossmeisl et al., *Chem. Phys. Lett.* 2007; Man et al.,
*ChemCatChem* 2011) turns that into a single first-principles number,
`eta` (the theoretical overpotential), plus which of the four steps is
responsible for it (the potential-determining step, PDS) — computable
from ground-state DFT total energies alone, no explicit electrochemical
interface, no applied potential, no reaction-barrier search (exactly
HER's own simplification, extended to 4 electrons and 3 intermediates
instead of 1).

This workflow automates the whole pipeline: finding every symmetrically
distinct OH-adsorption site, deriving O\* and OOH\* from that **same**
site (Section 1.2 explains why this is mandatory, not a convenience),
building every reference calculation the formula needs (H2/H2O gas-phase
references, THREE separate Boys-Bernardi BSSE counterpoise triads, a
vibrational zero-point/entropy correction for each of the three
intermediates and H2O), and combining them into `eta`/PDS with the
numerical safeguards a by-hand calculation is easy to get subtly wrong
on. Four stages, one more than HER's three: site search (Stage 1) →
derive O\*/OOH\* (Stage 2) → reference/BSSE/ZPE prep (Stage 3) → analysis
(Stage 4). None of the four run SIESTA for you — each is a real,
no-shortcuts folder-generation/analysis tool; you run SIESTA yourself in
between. This walkthrough proves the whole chain end to end with a
synthetic dataset whose exact `eta`/PDS are known in advance, so you can
see precisely what "correct" looks like before trusting it on your own
material.

## 1. Theory

### 1.1 The 4-electron CHE descriptor

```
* + H2O  -> OH*  + (H+ + e-)      dG1 = G(OH*)  + 1/2*G(H2) - G(*)    - G(H2O)
OH*      -> O*   + (H+ + e-)      dG2 = G(O*)   + 1/2*G(H2) - G(OH*)
O* + H2O -> OOH* + (H+ + e-)      dG3 = G(OOH*) + 1/2*G(H2) - G(O*)   - G(H2O)
OOH*     -> *    + O2 + (H+ + e-) dG4 = G(*) + G(O2) - G(OOH*) + 1/2*G(H2)

eta = max(dG1, dG2, dG3, dG4) - 1.23 V
PDS = whichever step has the largest dG
```

`1.23 V` is the experimental equilibrium potential of
`2H2O -> O2 + 4H+ + 4e-` (from the standard Gibbs free energy of water
formation, `-nFE` with `n=4`). At exactly this applied potential the
OVERALL 4-electron reaction is thermodynamically flat (net `dG = 0`), but
individual steps can still be uphill — `eta` is precisely the extra
voltage needed to make the WORST of the four steps favorable too. Read
`eta`:

- **`eta` close to 0** (roughly < 0.4 V): a genuinely good OER catalyst —
  the best known ones (IrO2, RuO2) sit around 0.3–0.4 V, close to the
  fundamental floor Section 1.6 explains.
- **`eta` 0.4–0.8 V**: moderate — most 3d-metal oxides/oxyhydroxides
  land somewhere in this range.
- **`eta` > 0.8 V**: poor — often a sign the material lacks a proper
  active site altogether (this is exactly what this walkthrough's own
  graphene toy example shows, Section 8, and it is the ordinary,
  expected answer for a material with no transition-metal center, not a
  broken calculation).

### 1.2 Why O\*/OOH\* MUST be derived from the SAME site as OH\*

This is the single most important methodological rule in this workflow,
and it did **not** used to be enforced: an earlier version of this suite
let O\* and OOH\* be searched *independently* on the clean slab
(`--o-strategy search`/`--ooh-strategy search`), on the theory that the
true lowest-energy O\*/OOH\* site need not coincide with OH\*'s own
winning site. In practice this was found to silently break the
descriptor's own physics and was removed entirely — `stb-oerIntermediates`
now *always* derives O\* and OOH\* from OH\*'s own relaxed geometry
(Section 3), with no flag to opt out. Two independent reasons converged
on the same conclusion:

1. **The reaction pathway itself only makes sense at one site.** `eta`
   models a single active center walking through
   `OH* -> O* -> OOH*` via sequential proton-coupled electron transfers.
   A "pathway" stitched together from the lowest-energy OH\* at one
   location, the lowest-energy O\* at a completely different one, and the
   lowest-energy OOH\* at a third does not correspond to anything a real
   catalytic site could physically do.
2. **The universal OOH\*/OH\* scaling relation assumes it.** Across a
   wide range of oxide catalysts, `dG(OOH*) - dG(OH*) ~= 3.2 +/- 0.2 eV`
   (Rossmeisl et al. 2007; Man et al. 2011) — a correlation that exists
   *because* OH\* and OOH\* both bond to the surface through the same
   metal-oxygen bond. Verified live on a real system during this
   workflow's own development: computing OH\*/O\*/OOH\* at independently
   -searched sites (different site TYPES for each) gave a gap of only
   `1.72 eV`, a `-1.48 eV` deviation from the literature band. Re-running
   the identical system with every intermediate constrained to the SAME
   site changed the gap only slightly (to `1.43 eV`) — proving the large
   deviation was a genuine property of that material (no transition
   -metal center to share the bond the way oxides do), NOT a
   site-consistency artifact — but only *after* removing the possibility
   of mixing sites could that distinction be drawn at all. With mixed
   sites, a suspiciously small/large gap is uninterpretable: you cannot
   tell whether it reflects real chemistry or an accounting error.

`stb-oerIntermediates` still supports **orientation sampling** of OOH\*'s
new O-H group AT the winning site (Section 4.3) — that is a genuinely
different thing (varying how the flexible part of the SAME molecule
points, never which surface location it's anchored to) and is fully
compatible with this rule; see Section 1.7 for why it needs its own,
separate safeguard.

### 1.3 Basis Set Superposition Error (BSSE) — always 3 separate triads

Same underlying issue as `4.8-adsorption`'s and `4.13-her`'s own BSSE
sections: because SIESTA is a localized-basis (LCAO) code, the naive
electronic energy of each intermediate is contaminated by BSSE — in the
combined `slab+adsorbate` calculation, each fragment "borrows" extra
basis functions that belong to the other fragment, artificially lowering
the combined energy. The Boys & Bernardi (1970) fix re-evaluates each
fragment WITH the other fragment's basis functions present as chargeless
"ghost" atoms, at the intermediate's own relaxed geometry:

```
BSSE(slab)      = E(slab, real)      - E(slab, real + adsorbate-ghost)
BSSE(adsorbate) = E(adsorbate, real) - E(adsorbate, real + slab-ghost)
BSSE(total)     = BSSE(slab) + BSSE(adsorbate)
E(corrected)    = E(raw) + BSSE(total)
```

Unlike HER (one adsorbate, one triad), OER writes **three independent
triads** — `05_bsse_OH_*/`, `05_bsse_O_*/`, `05_bsse_OOH_*/`, 12 folders
total, each evaluated at THAT intermediate's own geometry — never one
triad shared/reused across all three (that mode existed once, under the
name `--bsse-mode shared`, and was removed for the same reason as
Section 1.2's site-search removal: verified live, it systematically
UNDER-corrected OH\*/O\* relative to their true per-species BSSE (by
~0.30 eV for OH\* on a real system), inflating `eta` by a comparable
amount once fixed). A perhaps-counterintuitive real finding worth
knowing in advance: **the isolated single O\* atom's own BSSE component
can be the LARGEST of the three**, even larger than OOH\*'s (a bigger
molecule) — a bare atom has the FEWEST of its own basis functions to
fall back on, so the relative impact of "borrowing" the ghost slab's
functions is larger, not smaller. This is a real, physically expected
consequence of using a small, localized (DZP-class) basis for BSSE
counterpoise on isolated atoms/small fragments, not a bug — but it does
mean you should not be alarmed by a surprisingly large O\*-side BSSE
number in your own results.

### 1.4 Zero-point energy and entropy: local mode, generalized to N atoms

`--zpe-mode local` (the only mode this walkthrough exercises; `full`
exists too, same trade-off as HER's own Section 1.4) is a
**partial-Hessian, decoupled-oscillator** approximation — but where
HER's own H\* is always exactly 1 atom (a 3x3 Hessian), OER's local atoms
range from 1 (O\*) to 3 (OOH\*, and H2O's own gas-phase reference), a
genuine 3N x 3N mass-weighted Hessian that captures inter-atom coupling
within the adsorbate (e.g. the O-H stretch inside OH\*/OOH\*/H2O), not
just independent per-atom oscillators. Each local atom is displaced
`+/-x, +/-y, +/-z` (6 single-point folders per atom — `08_zpe_calc_OH/`
has 12 folders for 2 local atoms, `08_zpe_calc_O/` has 6 for 1,
`08_zpe_calc_OOH/`/`08_zpe_calc_H2O/` have 18 for 3), diagonalizing the
resulting Hessian gives 3N harmonic frequencies (with any near-zero/
negative eigenvalue — a translational/rotational residual, not a real
imaginary mode — reported and excluded from the ZPE/entropy sum, never
silently included). H2's own ZPE/entropy remain **fixed literature
values** in every mode (`H2_ZPE_EV = 0.270`, `H2_TS_298K_EV = 0.400`,
both at 298.15 K regardless of `--temp`), the same simplification HER
makes for the same reason — H2's vibrational structure is experimentally
extremely well characterized.

### 1.5 `config_extra.fdf`: forcing directives without touching your template

Every folder this workflow writes gets a `config_extra.fdf` sidecar,
`%include`-d at the very top of `calc.fdf` (SIESTA's fdf reader is
first-occurrence-wins, so this always overrides whatever the same tag
might say later in your own template) — the exact same mechanism HER,
`4.8-adsorption`, `4.11-raman` and `4.12-ir` use. Four directives are
**mandatory, not opt-in**, everywhere a slab is present: `MD.VariableCell
false`, `Slab.DipoleCorrection .true.` (adsorbing OH/O/OOH on only ONE
face breaks whatever mirror/inversion symmetry the clean slab had,
giving the cell a net dipole along a PERIODIC direction), `Spin
polarized` (every intermediate here can genuinely be open-shell — O\* in
particular, a bare oxygen atom, has a real triplet ground state), and
`DFTD3 .true.` (dispersion correction, mandatory rather than merely
advised — see Section 1.7 for why this also has to match on the ML side
whenever `--ml-prerelax`/`--ml-rank` is used). H2/H2O are the two
deliberate exceptions to `Spin polarized` — Section 1.8.

### 1.6 Why no O2 SIESTA reference is ever computed

`G(O2)` never comes from a direct DFT total-energy calculation anywhere
in this workflow — it is **derived**:

```
G(O2) = 2*G(H2O) - 2*G(H2) + 4*1.23 eV   (= + 4.92 eV)
```

Standard GGA functionals systematically misdescribe the O2 triplet
ground state by at least ~0.3 eV (Sargeant et al., *J. Electroanal.
Chem.* 2021) — large enough to distort `eta` by a similar amount if used
directly. Deriving `G(O2)` from the experimental 4-electron reaction free
energy instead sidesteps that error entirely, at the cost of one
built-in identity that must ALWAYS hold, independent of your material's
actual electronic energies:

```
dG1 + dG2 + dG3 + dG4 == 4.92 eV   (exactly, by construction)
```

`stb-oerAnalysis` prints this sum as a live internal sanity check every
run (Section 6) — if it doesn't come out to `4.92` to numerical
precision, something in the arithmetic pipeline itself is broken, not
your input data.

### 1.7 Orientation sampling: a hemisphere, never the full sphere

Both `stb-oer` (OH\* orientation, Section 3.3) and
`stb-oerIntermediates` (OOH\*'s new O-H group, Section 4.3) can
systematically sample multiple starting orientations per site
(`--n-orientations-polar`/`--n-orientations-azimuthal`,
`--ooh-n-orientations-polar`/`--ooh-n-orientations-azimuthal`), then
either write every one unscreened or MACE-MP-0 pre-rank them
(`--ml-rank`/`--ml-prerelax`) and keep only the best/`--orientation
-top-k`. The two tools sample fundamentally differently, though, and
this is a real, previously-discovered pitfall worth understanding:

- `stb-oer`'s OH\* sampling uses `core.adsorption_sites.
  generate_systematic_orientations` (a Fibonacci-sphere over the FULL
  sphere), which is correct there because `AdsorbateSiteFinder.
  add_adsorbate` always TRANSLATES the molecule afterward so whichever
  end ends up most-negative-z touches the surface, regardless of which
  way the un-rotated molecule initially pointed.
- `stb-oerIntermediates`'s OOH\* sampling has no such translation — O1
  (the winning OH\* site's own already-fixed relaxed oxygen, Section 1.2)
  never moves, and the new O2/H atoms are placed DIRECTLY as an offset
  from it. Verified live: reusing the full-sphere sampler here drove a
  downward-pointing sample straight into the substrate (0.32 Ang from a
  slab atom), producing a MACE "energy" around **-134,800 eV** instead of
  the correct ~-172 eV — a textbook force-explosion artifact, not a
  real result. The fix (now the only code path) restricts OOH\*
  orientation sampling to a HEMISPHERE pointing away from the surface
  (polar angle 0 up to 70 degrees from the local surface normal), which
  by construction can never place O2/H below O1's own height.

The practical takeaway: if you ever see MACE pre-relax energies of
extreme, physically nonsensical magnitude while using orientation
sampling anywhere in this suite, suspect exactly this class of bug
before suspecting your SIESTA settings.

### 1.8 Why H2 and H2O alone are forced spin-UNpolarized

Both are unambiguous **closed-shell singlets** — no radical character at
all, unlike every adsorbed intermediate here (which genuinely can be
open-shell, Section 1.5). Forcing `Spin polarized` on them anyway isn't
just unnecessary, it's actively risky: a spin-polarized SCF can converge
to a spurious nonzero magnetic moment depending on the initial guess,
silently biasing `1/2*G(H2)`/`G(H2O)` in every reaction step computed
against them. `02_h2_molecule/` and `03_h2o_molecule/`'s own
`config_extra.fdf` are the only two folders in this entire workflow with
`Spin non-polarized` instead of `Spin polarized` — forced explicitly, so
the choice is self-documenting and immune to your own `--calc`
template's own setting.

### 1.9 Fragment labeling: never confusing an intermediate with the slab

Every atom Stage 1 writes is relabeled by which fragment it came from —
`<symbol>_slab` for the substrate, `O_ads`/`H_ads` for the adsorbed OH —
exactly HER's own convention (its README's Section 1.7), just with two
adsorbate species sharing one suffix instead of one. This survives every
downstream derivation: O\* (Stage 2) reuses O1's own label verbatim after
removing H; OOH\* (Stage 2) gives its two brand-new atoms (O2, new H)
that SAME `_ads` suffix; every BSSE ghost-triad construction (Stage 3)
relies on the label surviving intact to know which atoms to ghost.

### 1.10 What this workflow simplifies — explicitly

By hand: writing 20+ nearly-identical `calc.fdf` variants across 4
stages; deriving O\*/OOH\* from the correct (same) site by hand every
time; getting three independent BSSE triads' sign conventions right;
deriving a coupled 3N x 3N harmonic ZPE/entropy Hessian per intermediate
from finite-difference forces; deriving `G(O2)` from the experimental
4.92 eV identity instead of a doomed direct DFT calculation; and finally
combining a dozen-plus numbers into 4 reaction steps and an overpotential
without a sign error anywhere. This workflow does all of that
mechanically and consistently — and, just as importantly, tells you
**which residual-force warnings are expected artifacts of the method
itself and which are real** (Section 7), and generates three
publication-ready charts (Section 6.4) automatically.

## 2. Libraries and external dependencies used

- **`pymatgen` / `spglib`** — `AdsorbateSiteFinder` (site search,
  Stage 1) and `SpacegroupAnalyzer` (the `[1] SLAB SYMMETRY` section).
- **`numpy`** — every local-mode Hessian diagonalization and
  Bose-Einstein entropy formula (Section 6.3).
- **`matplotlib`** — `stb-oerAnalysis --show`, an optional on-screen
  preview of the same three charts the gnuplot files draw.
- **`gnuplot`** (optional, `stb-oerAnalysis --plot`) — writes all three
  charts under `plot/` (Section 6.4).
- **`mace-torch`/PyTorch** — OPTIONAL, only if you pass `--ml-rank`
  (Stage 1) or `--ml-prerelax` (Stage 2); this walkthrough deliberately
  never uses either, so it needs none of it (no network access, no ML
  dependency, matching HER's own example).

## 3. Stage 1: adsorption site search (`stb-oer`, code `4.14.1`)

### 3.1 What it does

Given a slab/2D `structure.fdf` (vacuum along exactly one axis) and a
`calc.fdf` template, `stb-oer` finds every symmetrically distinct
OH-adsorption site of the requested type(s) and writes one
`sites/site_N_<type>/` folder per site — always every distinct site,
exactly HER's own "OER/HER wants the global minimum, not a sample"
philosophy.

### 3.2 Site search, live

```
$ stb-oer -s structure.fdf -c calc.fdf --no-intro
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

[3] ADSORPTION SITES: FINDING & COUNT
------------------------------------------------------------
Site type | Raw candidates | After symm. reduction | Reduced away (equivalent)
------------------------------------------------------------------------------
ontop     | 2              | 1                     | 1
bridge    | 6              | 3                     | 3
hollow    | 0              | 0                     | 0
TOTAL     | 8              | 4                     | 4
```

The same 4 sites (1 ontop, 3 bridge, 0 hollow) HER's own example finds on
this identical graphene fixture — expected, since Stage 1's site-finding
machinery is the same `AdsorbateSiteFinder`-based approach either way,
just adsorbing OH instead of H.

### 3.3 Orientation sampling — new relative to HER

OH\* is a 2-atom rod, not a single atom, so it has a real orientation
degree of freedom HER's own H\* never had:

```
$ stb-oer --site-type ontop --n-orientations-polar 2 --n-orientations-azimuthal 2
```

writes `site_1_ontop_orient1` .. `site_1_ontop_orient4` — 4 systematically
sampled starting orientations of the SAME site, unscreened (every one
written directly) since `--ml-rank` wasn't passed. See Section 1.7 for
the full-sphere-vs-hemisphere distinction (this tool, unlike Stage 2's
OOH\* sampling, correctly uses the full sphere, because of the touchdown
-translation mechanics explained there).

### 3.4 Fragment labels and `config_extra.fdf`, live

```
%block ChemicalSpeciesLabel
 1   6   C_slab
 2   8   O_ads
 3   1   H_ads
%endblock ChemicalSpeciesLabel
```

```
$ cat sites/site_1_ontop/config_extra.fdf
```
```
MD.VariableCell false
Slab.DipoleCorrection      .true.
Spin                polarized
DFTD3                   .true.
```

See Sections 1.5/1.9 for why each of these is mandatory.

### 3.5 `--both-sides`, `--position`, `--positions-file`

Identical mechanics to HER's own Section 3.5/3.6 (manual override along
the TRUE surface normal, not a naive Cartesian z-offset; a free-standing
2D material's second face; `sites/site_positions.dat` round-tripping) —
`--both-sides` on this same free-standing graphene fixture gives a
6-atom `site_1_ontop_bothsides/` folder (2 C + 1 OH per face).

### 3.6 The too-close-atoms safety check

```
$ stb-oer --site-type ontop --height 0.1
[WARNING] site_1_ontop: closest slab-OH distance is only 0.100 Ang --
likely overlapping atoms (--height too small, or a bad
--position/--positions-file entry). Check before running SIESTA.
```

A sanity check on `--height` (or a hand-picked `--position`), not a
physical statement — no real X-OH bond is anywhere near 0.1 Ang.

### 3.7 Report structure

```
[0] RUN METADATA
[1] SLAB SYMMETRY
[2] CLEAN SLAB REFERENCE
[3] ADSORPTION SITES: FINDING & COUNT
[4] WRITING SITE FOLDERS
[5] SUMMARY & NEXT STEPS
[6] LIBRARY WARNINGS
```

Identical shape to HER's own Stage 1 report — see that README's Section
3.8 for what `[2]`'s recentering and `[6]`'s library-noise trapping do.
`[5]` ends with three `[NOTE]`s worth reading once: later stages always
read sites straight from these exact folders (never re-derive them), no
O2 gas-phase reference is EVER generated (Section 1.6), and this whole
descriptor is AEM-only (Section 7).

### 3.8 Running it both ways

**A — direct CLI**: every command above.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.14.1
```

`4.14.1` asks for the structure file, the calc.fdf template, a
pseudopotential source, the site type, height, both-sides, an
orientation-sampling grid (blank to skip), output directory, then
optional symprec/vacuum-gap overrides, and finally whether to save a
report — the same flags as the CLI, in the same order.
`example_4.14.sh` proves both paths agree.

## 4. Stage 2: derive O\*/OOH\* from the SAME site (`stb-oerIntermediates`, code `4.14.2`)

### 4.1 What it does

Scans Stage 1's `sites/site_*/` for the lowest-`FreeEng` finished site,
reads back its relaxed geometry, then derives (Section 1.2 — **always**
from that same site, no exceptions): `O* = OH*` minus H (3 atoms for this
fixture: 2 C + O), and `OOH* = OH*` plus a second O-H group in an
illustrative starting orientation (5 atoms: 2 C + O + O + H). Neither is
a search — both write exactly ONE folder each by default
(`intermediates/o_star/`, `intermediates/ooh_star/`), unless OOH\*
orientation sampling is requested (Section 4.3).

### 4.2 Live output

```
[2] O* GEOMETRY
------------------------------------------------------------
  [OK] intermediates/o_star (derived from the winning OH* site: H removed)
  [Saved] intermediates/o_trajectory.xyz (1 frame, OVITO/VMD-viewable)

[3] OOH* GEOMETRY
------------------------------------------------------------
  [OK] intermediates/ooh_star (derived from the winning OH* site: +O at
       1.450 Ang, +H bent 100.0 deg -- illustrative starting geometry,
       refined by CG relaxation)
  [Saved] intermediates/ooh_trajectory.xyz (1 frame, OVITO/VMD-viewable)
```

The O-O bond length (`--oo-bond-length`, default 1.45 Ang) and O-O-H bend
angle (`--ooh-bend-deg`, default 100.0, H2O2-like) are explicitly
illustrative STARTING geometries, refined by the real CG relaxation you
run afterward — not literature-fitted equilibrium values, the same
"CG relaxation refines it" spirit as HER's own `--height` default.

### 4.3 OOH\* orientation sampling — same idea as Stage 1, different mechanics

```
--ooh-n-orientations-polar 4 --ooh-n-orientations-azimuthal 4 --ml-prerelax --orientation-top-k 2
```

samples 16 starting orientations of OOH\*'s new O-H group AT the SAME
site (never a different one — Section 1.2), optionally MACE-MP-0 ranks
them (`--ml-prerelax`, needs the optional `ml` extra — skipped in this
walkthrough), and keeps only the `--orientation-top-k` best, writing
`intermediates/ooh_star_orient1/`, `ooh_star_orient2/`, etc. Section 1.7
explains the hemisphere-restricted sampling this uses internally, and
why that restriction exists.

### 4.4 Report structure

```
[0] RUN METADATA
[1] WINNING OH* SITE
[2] O* GEOMETRY
[3] OOH* GEOMETRY
[4] SUMMARY & NEXT STEPS
```

### 4.5 Running it both ways

**A — direct CLI**: `stb-oerIntermediates --directory oer_study`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.14.2
```

`4.14.2` asks for the Stage-1 directory, the SIESTA output filename, a
pseudopotential source, whether to MACE pre-relax, an OOH\*
orientation-sampling grid (blank to skip), then advanced settings
(bond lengths/bend angle, ML device), and finally whether to save a
report.

## 5. Stage 3: references, BSSE (3 triads) & ZPE prep (`stb-oerRefs`, code `4.14.3`)

### 5.1 What it does

Locates OH\*'s winning site and O\*/OOH\*'s derived, relaxed geometries
(Stage 2), then writes: `00_clean_slab/` (pristine slab, single-point),
`02_h2_molecule/`/`03_h2o_molecule/` (RELAX, gas-phase CHE references,
Gamma-only, forced spin-unpolarized — Section 1.8), `04_slab_deformed/`
(OH\*'s geometry with both adsorbate atoms removed — diagnostic only,
Section 7), the 3 separate BSSE triads (Section 1.3), and the
`--zpe-mode`-dependent ZPE folder(s) for OH\*, O\*, OOH\* AND H2O
(Section 1.4).

### 5.2 The BSSE triad table, live

```
[5] BSSE CORRECTION
------------------------------------------------------------
  [OK] 3 separate triads written (05_bsse_OH_*, 05_bsse_O_*, 05_bsse_OOH_*),
       each at its own intermediate's geometry

Folder                           | Atoms | Formula        | Run type
-------------------------------------------------------------------------------
05_bsse_OH_slab_only             | 2     | C2             | single-point (BSSE)
05_bsse_OH_slab_ghost            | 4     | C2 +HO(ghost)  | single-point (BSSE)
05_bsse_OH_adsorbate_ghost_slab  | 4     | HO +C2(ghost)  | single-point (BSSE)
05_bsse_OH_isolated              | 2     | HO             | single-point (BSSE)
05_bsse_O_slab_only              | 2     | C2             | single-point (BSSE)
05_bsse_O_slab_ghost             | 3     | C2 +O(ghost)   | single-point (BSSE)
05_bsse_O_adsorbate_ghost_slab   | 3     | O +C2(ghost)   | single-point (BSSE)
05_bsse_O_isolated               | 1     | O              | single-point (BSSE)
05_bsse_OOH_slab_only            | 2     | C2             | single-point (BSSE)
05_bsse_OOH_slab_ghost           | 5     | C2 +HO2(ghost) | single-point (BSSE)
05_bsse_OOH_adsorbate_ghost_slab | 5     | HO2 +C2(ghost) | single-point (BSSE)
05_bsse_OOH_isolated             | 3     | HO2            | single-point (BSSE)
```

Ghost atoms are called out explicitly in the formula (`+HO(ghost)`) —
real element, real pseudopotential file, zero valence charge,
contributing basis functions only. Notice `05_bsse_O_isolated` is a
genuine 1-atom folder — Section 7 has an important MPI-ranks caveat for
running SIESTA there.

### 5.3 The ZPE folder-count table, live

```
[6] ZPE PREPARATION
------------------------------------------------------------
Intermediate | Local atom(s) | Folder range       | Total folders | Displacement
--------------------------------------------------------------------------------
OH           | 2             | disp_001..disp_012 | 12            | 0.015 Ang
O            | 1             | disp_001..disp_006 | 6             | 0.015 Ang
OOH          | 3             | disp_001..disp_018 | 18            | 0.015 Ang
H2O          | 3             | disp_001..disp_018 | 18            | 0.015 Ang
```

54 folders total (local mode) — O\*'s single local atom needs only 6
(one atom x 3 axes x 2 signs); OOH\*/H2O's three each need 18.

### 5.4 Report structure

```
[0] RUN METADATA
[1] WINNING OH* SITE
[2] O* FINAL GEOMETRY
[3] OOH* FINAL GEOMETRY
[4] REFERENCE FOLDERS
[5] BSSE CORRECTION
[6] ZPE PREPARATION
[7] SUMMARY & NEXT STEPS
```

### 5.5 Running it both ways

**A — direct CLI**: `stb-oerRefs --directory oer_study --zpe-mode local`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.14.3
```

`4.14.3` asks for the Stage-1/2 directory, the SIESTA output filename, a
pseudopotential source, the ZPE mode, then (if `local`/`full`) the
displacement size and (if `full`) the supercell, and finally whether to
save a report. There is no BSSE-mode question at all — it is always the
3-triad behavior, unconditionally (Section 1.3).

## 6. Stage 4: analysis (`stb-oerAnalysis`, code `4.14.4`)

### 6.1 Report structure

```
[0] RUN METADATA
[1] ELECTRONIC ENERGIES (FreeEng)
[2] BSSE CORRECTION
[3] THERMAL CORRECTION
[4] GAS-PHASE O2 REFERENCE (DERIVED)
[5] REACTION STEPS & POTENTIAL-DETERMINING STEP
[6] FINAL RESULT
[7] SUMMARY & FILES
```

### 6.2 `[1]`/`[2]`: every energy, every BSSE term, all three triads, in tables

Every raw `FreeEng` is printed inline (with its own SCF/force quality
diagnostic, Section 7) AND recapped in a table at the end of `[1]`;
`[2]` prints the full derivation for **all three** BSSE triads — slab
term, adsorbate term, total, for OH\*, then O\*, then OOH\* — followed by
a clean recap table of raw vs. corrected energies:

```
Intermediate | E (raw, eV) | BSSE (eV) | E (corrected, eV)
----------------------------------------------------------
OH*          | -851.8709   | +0.1000   | -851.7709
O*           | -834.7918   | +0.1500   | -834.6418
OOH*         | -1287.9627  | +0.0500   | -1287.9127
```

Three genuinely different BSSE values — proof the 3-separate-triad
requirement (Section 1.3) is actually being honored, not silently
collapsed back to one shared number.

### 6.3 `[3]`: thermal correction, all four species

Same local-mode Hessian diagonalization as HER's own Section 5, just run
four times (OH\*, O\*, OOH\*, H2O) instead of once, with the coupled 3N x
3N generalization Section 1.4 describes. Each intermediate's own table of
harmonic frequencies (with any near-zero eigenvalue reported as
`(IMAGINARY) -- excluded from ZPE/entropy`, never silently summed in) is
printed before its `ZPE(...)`/`TS(...)` totals, then a compact recap
table across all four species.

### 6.4 `[6]`/charts: the final result, and THREE saved charts

```
Theoretical overpotential eta = +0.9700 V
Potential-determining step (PDS) = Step 1
```

`--plot`/`--show` (tri-state, asked interactively if not passed, exactly
HER's own convention) now generate/display **three** charts together
(all always embedded, as PNGs, in the always-written `OER_report.md` too
— Section 6.6 — regardless of `--plot`/`--show`):

- **`oer_free_energy_diagram.{dat,gplot,pdf,png}`** — the CHE
  free-energy diagram, drawn in the standard literature STEP/PLATEAU
  form (a short horizontal "shelf" per state, dashed lines only guiding
  the eye between them — never a slanted line-through-a-point, which
  would misleadingly suggest something exists partway between two
  discrete states) at both `U=0 V` and `U=1.23 V` (equilibrium), with
  the PDS step shaded.
- **`oer_bsse_correction.{dat,gplot,pdf,png}`** — a categorical bar
  chart of all 9 BSSE terms (slab/adsorbate/total x 3 intermediates),
  the same numbers as `[2]`'s own table, visual proof the three triads
  differ.
- **`oer_bsse_raw_corrected.{dat,gplot,pdf,png}`** — raw vs. BSSE
  -corrected energy, one SMALL-MULTIPLE PANEL per intermediate, each
  with its OWN tightly-zoomed y-axis. This is deliberate, not a
  simplification: BSSE (~0.05-0.5 eV) is 3-4 orders of magnitude smaller
  than the raw electronic energies (~hundreds to thousands of eV), so a
  single shared/zero-based axis would render the correction as a
  literally invisible sliver on an enormous bar. Each panel draws Raw
  and Corrected as two bars on a window scaled to just those two nearby
  values, with the BSSE shift labeled directly between them.

### 6.5 Quality diagnostics — and a real gap vs. HER worth knowing

Every energy is checked for SCF convergence and residual force against
`--force-tolerance` (default `0.05` eV/Ang). HER's own Stage 3 skips this
check entirely for its two analogous non-equilibrium diagnostic folders
(`03_slab_deformed`/`04_slab_ghost`, its own README Section 5.3) — a real
residual force there is EXPECTED by construction (the atoms nearest a
just-removed/ghosted atom are never at their own equilibrium), so
warning about it every single run would be pure noise. **`stb-oerAnalysis`
does not carry that same skip-list**: `04_slab_deformed` and every
`05_bsse_*/` folder here are equally non-equilibrium single-point
snapshots by construction, but none of them are exempted — expect
`[WARNING] Residual force on ...` on these specific folders in ordinary
use, and read it with that context rather than assuming the calculation
is broken. This is an honest, currently-real gap between the two
sibling workflows, not a design choice being defended — Section 8's
worked example demonstrates it directly with a deliberately large force
on `04_slab_deformed`.

### 6.6 `OER_report.md`: a consolidated Markdown deliverable

Always written (unconditionally, like `OER_report.txt`), consolidating
every number from every stage into one Markdown file with all 9 report
sections as proper tables, PLUS all three charts as embedded PNGs — built
from THIS run's own live-computed values, never by re-parsing/
concatenating the earlier stages' own `.txt` reports (which may not even
exist if you ran an earlier stage without `--save-report`).

### 6.7 Running it both ways

**A — direct CLI**: `stb-oerAnalysis --directory oer_study --save-report`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.14.4
```

`4.14.4` asks for the Stage-1/2/3 directory, the SIESTA output filename,
temperature, force tolerance, an output-file basename, then whether to
save a report — the plot/show questions are asked by `stb-oerAnalysis`
itself afterward, exactly as from the direct CLI.

## 7. Known, deliberate limitations

- **No SIESTA run is ever performed by this suite.** All four stages
  generate/read input and output files; you run SIESTA in every folder
  yourself, in between stages.
- **Rerunning SIESTA in the same folder does NOT continue from where it
  left off** — same behavior as every other workflow in this suite
  (HER's own README documents this in detail); SIESTA restarts geometry
  from `structure.fdf` itself, not a saved `.XV` file.
- **`stb-oerAnalysis` does not skip the force-quality check on
  `04_slab_deformed`/`05_bsse_*/`** the way HER skips its own analogous
  folders (Section 6.5) — expect (and correctly ignore) a residual-force
  warning there in ordinary use.
- **A single isolated atom cannot use many MPI ranks** — `05_bsse_O_
  isolated` (a genuine 1-atom folder) will make SIESTA abort with the
  same `mpirun -np N` you used elsewhere; run it with far fewer ranks
  (same caveat HER's own `07_h_isolated` has).
- **O\*'s own BSSE correction can be surprisingly large** (Section 1.3)
  — a real, expected consequence of counterpoise-correcting a bare atom
  with a small localized basis, not a red flag on its own.
- **A large `eta` is not automatically a bug.** A chemically inert,
  large-gap pristine 2D material (this walkthrough's own graphene
  fixture, for instance) can genuinely have a strongly unfavorable
  `eta` — that IS the physically correct answer for a material with no
  transition-metal active site, not a sign something needs re-running.
- **`--zpe-mode local` ignores adsorbate-substrate vibrational coupling**
  by construction, same limitation as HER's own local mode.
- **AEM descriptor only — this workflow does NOT model the Lattice
  Oxygen Evolution Mechanism (LOER)**, where lattice oxygen atoms
  themselves participate directly (common on some oxide/oxyhydroxide
  catalysts, and linked to catalyst instability under anodic potential;
  see Exner, *ChemCatChem* 2021, "On the Lattice Oxygen Evolution
  Mechanism: Avoiding Pitfalls"). If your material is expected to favor
  LOER (a perovskite, or a Ru/Ir oxide operating at high anodic
  potential, for instance), treat this workflow's `eta` as an AEM-only
  estimate, not the full mechanistic picture — `stb-oerAnalysis` prints
  this same caveat every run.
- **Only ONE site's full pathway is ever evaluated.** Stage 1 finds
  every symmetrically distinct OH\* site, but only the single
  lowest-energy one is carried forward through Stages 2-4 (Section
  1.2) — a genuinely different, higher-energy OH\* site might turn out
  to support a *better* overall pathway (a smaller `eta`) than the
  globally lowest-energy OH\* site does; this workflow does not check.
- **The OOH\*/OH\* scaling-relation gap is a useful sanity check, not a
  hard requirement.** A real material's gap deviating substantially from
  the literature's ~3.2 +/- 0.2 eV band (Section 1.2) is not itself proof
  of an error — it may simply mean the material lacks the shared
  metal-oxygen bonding motif the correlation is built on (verified for a
  real system during this suite's own development, Section 1.2) — but a
  SUSPICIOUSLY far-off gap is still worth double-checking your BSSE/ZPE
  inputs over before trusting the result.

## 8. Worked example: a known analytic ground truth, end to end

Real SIESTA output isn't available inside this walkthrough, so — exactly
like HER's own Section 6 — the script fabricates `calc.out`/`.FA` files
with a **hand-chosen** set of numbers, worked BACKWARD from clean,
verifiable `dG1..dG4` values so the correct `eta`/PDS are known before
running Stage 4:

```
Target (chosen by hand):
  dG1 = +2.20 eV   dG2 = +1.10 eV   dG3 = +0.70 eV   dG4 = +0.92 eV
  sum = 4.92 eV exactly (Section 1.6's identity, by construction)
  eta = max(dG1..dG4) - 1.23 = 2.20 - 1.23 = +0.97 V
  PDS = Step 1 (* + H2O -> OH*)
```

Three DISTINCT BSSE corrections (Section 1.3):

```
BSSE(OH*)  = +0.10 eV (slab +0.06, adsorbate +0.04)
BSSE(O*)   = +0.15 eV (slab +0.09, adsorbate +0.06)
BSSE(OOH*) = +0.05 eV (slab +0.02, adsorbate +0.03)
```

ZPE/TS via isotropic harmonic springs (Section 1.4), `k_O = 8.0 eV/Ang^2`,
`k_H = 5.0 eV/Ang^2` — an exact, hand-verifiable closed form
(3 degenerate modes per atom, `f = 1/(2*pi) * sqrt(k/m)`, each
contributing `1/2*h*f` to the ZPE and a Bose-Einstein term to `T*S` at
298.15 K), computed independently in plain Python (not by calling
`stb-oerAnalysis`) for all three intermediates AND H2O simultaneously:

```
ZPE(OH*)  = 0.2846 eV   TS(OH*)  = 0.0440 eV
ZPE(O*)   = 0.0686 eV   TS(O*)   = 0.0421 eV
ZPE(OOH*) = 0.3532 eV   TS(OOH*) = 0.0860 eV
ZPE(H2O)  = 0.5006 eV   TS(H2O)  = 0.0458 eV
```

The electronic energies needed to hit those exact `dG` targets (given
`e_clean=-400`, `e_h2=-31.5`, `e_h2o=-470` and the BSSE/ZPE numbers
above) were then solved for algebraically — not round numbers themselves
(`e_oh = -851.870886 eV`, etc.), but exactly reproducible, and the whole
point is that the FOUR REACTION STEPS come out clean, not the raw inputs.
The script runs the real 4-stage pipeline on this fabricated data and
asserts the result:

```
eta = +0.9700 V
Potential-determining step (PDS) = Step 1
```

— confirming the entire pipeline (energy reading → three independent
BSSE corrections → four coupled harmonic ZPE/entropy Hessians → the
derived-O2 formula → eta/PDS) end to end, the same way HER's own
synthetic dataset confirms *its* simpler 1-intermediate pipeline.

## 9. Step-by-step: running this workflow on your own structure

1. Have a relaxed slab/2D `structure.fdf` (vacuum along exactly one
   axis) and a working `calc.fdf` template.
2. Run Stage 1: `stb-oer -s <structure.fdf> -c <calc.fdf> --site-type
   all [--both-sides] [-p <pseudo-bank-or-path>]`.
3. Run SIESTA yourself in every `oer_study/sites/site_*/` folder.
4. Run Stage 2: `stb-oerIntermediates --directory oer_study` (add
   `--ooh-n-orientations-polar`/`-azimuthal` + `--ml-prerelax` if OOH\*
   orientation matters for your system, Section 4.3).
5. Run SIESTA yourself in `intermediates/o_star/` and
   `intermediates/ooh_star/` (or every `ooh_star_orientN/`, if sampled).
6. Run Stage 3: `stb-oerRefs --directory oer_study --zpe-mode local`
   (or `full` if adsorbate-substrate coupling matters).
7. Run SIESTA yourself in every folder Stage 3 just wrote — remember
   `05_bsse_O_isolated` needs far fewer MPI ranks than the rest
   (Section 7).
8. Run Stage 4: `stb-oerAnalysis --directory oer_study --save-report`.
   Read `[1]`-`[3]` for any `[WARNING]` before trusting `eta` — Section
   6.5 explains exactly which ones are expected and which aren't.
   Answer `y` to the chart prompts (or pass `--plot`) for all three
   saved charts.

## 10. Files in this folder

| File | Purpose |
|---|---|
| `structure.fdf` | The SAME free-standing 2-atom graphene primitive cell `4.13-her`'s own example uses — small and fast, purely to exercise the tools; keeps the HER/OER example pair directly comparable. |
| `calc.fdf` | Shared template (non-polarized, `MD.VariableCell true`) — proves `config_extra.fdf`'s forced overrides actually win over the template's own settings. |
| `example_4.14.sh` | The guided walkthrough (**not** an automated test — see `test/4-workflow/14-oer/{prep,intermediates,refs,analysis}/test.sh` for that). Pauses between sections so you can read before moving on; safe to re-run. |
| `output/` | Created by `example_4.14.sh` when you run it (git-ignored, not checked in). See below. |

## 11. Running the script

```bash
./example_4.14.sh
```

| Case | Command(s) | What it shows |
|---|---|---|
| `output/stage1/` | `stb-oer -s structure.fdf -c calc.fdf` | Symmetry, raw-vs-reduced site table, fragment labels, `config_extra.fdf` |
| `output/stage1_orient/` | `stb-oer --site-type ontop --n-orientations-polar 2 --n-orientations-azimuthal 2` | OH\* orientation sampling (unscreened) |
| `output/stage1_both/` | `stb-oer --site-type ontop --both-sides` | Both faces of a free-standing 2D material |
| `output/stage1_position/` | `stb-oer --position 0.0 0.0 --height 1.6` | Manual site override, true-surface-normal height |
| `output/stage1_positionsfile/` | `stb-oer --positions-file mypositions.dat` | Round-tripping `site_positions.dat` |
| `output/stage1_overlap/` | `stb-oer --site-type ontop --height 0.1` | The too-close-atoms warning |
| `output/workflow/` | Stage 1 (single site) -> (fabricated relaxation/energies) -> Stage 2 -> Stage 3 (3 BSSE triads + 4-species local ZPE) -> Stage 4 | The full 4-stage chain, known-answer worked example (Section 8), the non-skip force-warning nuance (Section 6.5), and all three `--plot` charts |
| *(no folder — a diff only)* | Stage 1 via `printf … \| stb-suite` | Proof the interactive menu (`4.14.1`) agrees with the CLI |

## What's next

- **`4.13-her`** — the sibling 1-electron/1-intermediate descriptor this
  workflow generalizes to 4 electrons/3 intermediates; read it first if
  you haven't, for the simplest version of every idea here (a single
  BSSE triad, a single-atom local Hessian, one reaction step instead of
  four).
- **`4.8-adsorption`** — the general-purpose version of Stage 1's own
  BSSE and site-search machinery (any adsorbate, multi-adsorbate, ML
  pre-screening, height sweeps) — a deeper BSSE derivation than either
  workflow's own README repeats in full.
- **`4.9-neb`** — once you know a particular step's `dG` is favorable,
  this is how you'd estimate the actual reaction *barrier* between two
  adsorption geometries, rather than just the thermodynamic endpoint
  this workflow computes.
- **`3.5-stb-symmetry`** — the full per-site/operations symmetry
  analysis this workflow's own `[1] SLAB SYMMETRY` section is a
  summarized slice of.
- Every other `4.x` workflow example in this suite generates real
  candidate/perturbed geometries and expects a real SIESTA run in
  between stages, the same two-way (direct CLI / interactive menu)
  split, and the same numbered-report convention — if this is your
  first workflow example, `4.1-strain` is the shortest two-stage one to
  start with.
