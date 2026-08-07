# 4.7 — Workflow: Hubbard U (Linear Response) (`stb-hubbardu` / `stb-hubbarduAlphas` / `stb-hubbarduAnalysis`)

Plain DFT (LDA/GGA) systematically gets transition-metal and rare-earth
compounds wrong — wrong band gaps, wrong magnetic moments, sometimes even
the wrong ground state. The standard fix is a DFT+U correction, but DFT+U
needs a number, `U`, and picking it by hand (or copying one from a paper
that used a different code, functional, and basis set) is exactly the kind
of unverifiable shortcut this suite tries to avoid. `stb-dftu` (example
`1.5`) already lets you *use* a literature `U`; this workflow computes a
**first-principles** one for *your own* structure, functional, pseudopotential,
and basis set, via Cococcioni & de Gironcoli's linear-response method
(Phys. Rev. B 71, 035105, 2005).

Unlike every other workflow example in this suite so far, this one has
**three** stages, not two: a reference calculation (Stage 1), a set of small
perturbations around it (Stage 2), and a fit that turns the perturbed
occupations into a number (Stage 3). None of the three run SIESTA for
you — each is a real, no-shortcuts folder-generation/analysis tool; you run
SIESTA yourself in between. This walkthrough proves the whole chain end to
end with a synthetic dataset whose exact answer is known in advance, so you
can see precisely what "correct" looks like before trusting it on your own
material.

## 1. Theory

### 1.1 Why plain DFT gets correlated compounds wrong

LDA and GGA approximate the exchange-correlation energy as a (semi-)local
functional of the density. For simple metals and covalent semiconductors
this is a good approximation. It breaks down for compounds with partially
filled, spatially localized `d` or `f` shells (transition-metal oxides,
lanthanides, actinides) because of the **self-interaction error**: an
electron in LDA/GGA feels part of its own charge through the mean-field
Hartree term, an error that isn't cancelled cleanly by the approximate
exchange-correlation piece. For delocalized states this error is small; for
a tightly localized `d`/`f` orbital it is not, and it shows up as a
spurious *curvature* of the total energy as a function of that orbital's
occupation — real physics says removing/adding a fraction of an electron
from a localized orbital should cost energy in a specific, sometimes
discontinuous way (the derivative discontinuity of the exact functional);
LDA/GGA instead gives a smooth, convex curve that favours **fractional
occupation**, artificially delocalizing electrons that should stay
localized. The textbook symptom is NiO/MnO/FeO predicted as metals or
small-gap semiconductors with plain GGA, when they are real, wide-gap
antiferromagnetic insulators.

### 1.2 DFT+U: penalizing fractional occupation

DFT+U (Anisimov, Zaanen & Andersen, 1991) adds a Hubbard-model-like term
that penalizes fractional occupation of the chosen correlated shell and
rewards integer occupation (0 or 1 per orbital) — pushing the functional's
energy-vs-occupation curve back toward the physically-correct shape. In
SIESTA this is the `%block LDAU.proj` mechanism `stb-dftu` (example `1.5`)
already builds a snippet for, given a `U` (and `J`) you supply. The
correction only ever does its job if `U` is *right* for your exact setup —
functional, pseudopotential, and basis set all affect it — which is exactly
why a copied literature value is a starting-point sanity check
(`stb-dftu --use-reference`/`--suggest`), never a substitute for computing
one from your own data.

### 1.3 The Cococcioni & de Gironcoli linear-response method

Cococcioni & de Gironcoli's insight: `U` is, by definition, the *curvature*
of the total energy with respect to the correlated shell's occupation —
`U = d²E/dn²`. Instead of computing that second derivative directly (hard),
they compute it indirectly, from how the shell's **occupation itself**
responds to a small external perturbation:

1. Apply a rigid potential shift `alpha` to the correlated shell only (in
   SIESTA: `%block LDAU.proj` with `LDAU.PotentialShift T`, setting the
   would-be `U` slot in that block to the perturbation strength `alpha`
   instead of a real Hubbard correction — the same block, reused as a
   *diagnostic probe* rather than a correction).
2. Measure how the shell's SCF occupation `n` changes as `alpha` is swept
   through a handful of small positive and negative values.
3. The slope `dn/dalpha` is a **response function**, `chi`. Two versions of
   it — self-consistent and frozen-density (next subsection) — combine into
   `U` via a formula (Section 1.4) derived from exactly this occupation
   response, without ever computing `d²E/dn²` directly.

Stage 1 (`stb-hubbardu`) prepares the `alpha=0` reference point; Stage 2
(`stb-hubbarduAlphas`) prepares every other `alpha`; Stage 3
(`stb-hubbarduAnalysis`) fits the slopes and applies the formula.

### 1.4 Screened (`chi`) vs. bare (`chi0`) response, and the U formula

A perturbation on one shell doesn't just change that shell's own
occupation — the rest of the system reorganizes (screens) in response too.
Two different response slopes matter, and the difference between them is
exactly what isolates the correlation-driven curvature from an ordinary,
uninteresting electrostatic response:

- **`chi` — the self-consistent (screened) response**: each perturbed
  calculation is run to full self-consistency (`scf_alpha_*`/`reference`
  folders). This includes the "real" physical response — the rest of the
  electronic density relaxing around the perturbed shell.
- **`chi0` — the frozen-density (bare) response**: each perturbed
  calculation is capped at a handful of SCF iterations
  (`frozen_alpha_*` folders, `MaxSCFIterations` deliberately small) — just
  enough for the perturbed shell's own occupation to react to `alpha`,
  before the rest of the density has a chance to relax around it.

`chi` and `chi0` are each fit as the slope of a straight line, occupation
`n` vs. perturbation `alpha` (Stage 3, `numpy.polyfit`, degree 1). The
final formula:

```
U = 1/chi0 - 1/chi
```

Subtracting the *bare* response's reciprocal from the *screened* response's
reciprocal removes exactly the part of the response that isn't a genuine
electron-correlation effect, leaving the physically meaningful on-site `U`.
Section 6 below fits both branches against a known, hand-constructed
occupation-vs-alpha dataset and recovers the exact expected `U`, which is
the cleanest way to see this formula actually work before trusting it on
noisy, real SCF numbers.

### 1.5 Why every perturbed run needs the exact same starting point

The measured slope is only meaningful if the *only* thing that changes
between runs is `alpha` — not the geometry, not the SCF starting point.
Two safeguards make this true, both enforced automatically by Stage 1/2,
never left to you to remember:

- **Fixed geometry**: every generated folder (reference and every
  perturbed `alpha`) forces a pure single-point SCF (`MD.TypeOfRun CG`,
  `MD.Steps 0`, `MD.NumCGsteps 0`, `MD.VariableCell false`) regardless of
  what your own `calc.fdf` template configures. Letting the ions (or the
  cell) relax under a rigid, unphysical `alpha` shift would change the
  geometry between perturbation points — contaminating the measured
  response with a geometric effect that has nothing to do with the
  electronic occupation response the method needs.
- **Fixed starting density**: Stage 2 physically copies Stage 1's
  converged `.DM` file into every perturbed folder and sets
  `DM.UseSaveDM T`. Without this, each perturbed SCF cycle would start from
  its own default initial guess and could converge to a slightly different
  point on the potential-energy surface (especially likely for a
  correlated, often multi-minima system) — again contaminating the slope
  with numerical noise unrelated to the real electronic response.

### 1.6 Isolating a single atom: species aliasing

`%block LDAU.proj` perturbs an entire species *label*, not one atom. If
the target species appears more than once in the structure (e.g. two
equivalent Mn atoms in a supercell), perturbing the species would perturb
*both* atoms at the exact same time — measuring a **collective** response
of two coupled perturbed sites, not the **single isolated-atom** response
the Cococcioni-de Gironcoli supercell method actually requires (the whole
point of using a big-enough supercell is that the *other* atoms of the
same species stay unperturbed, acting as the "bulk" the one perturbed atom
responds against). `stb-hubbardu` handles this by aliasing the one chosen
atom (`--atom-index`) onto a new, same-`Z` species label
(`<species>_pert`) behind the scenes — only that aliased label goes into
`%block LDAU.proj`; every other atom of the original species is left
completely alone. Section 3.3 demonstrates this live, including the
`%block ChemicalSpeciesLabel`/coordinates-block edits it makes.

### 1.7 What this workflow simplifies — explicitly

By hand, a linear-response `U` calculation means: hand-writing 9+ nearly
identical `calc.fdf` files that differ only in one perturbation-strength
number; manually copying a converged `.DM` file into each one and adding
`DM.UseSaveDM T`; remembering to cap `MaxSCFIterations` for exactly the
"frozen" half of them and not the other half; parsing an `Occupations:`
line out of 9+ `.out` files by hand (including the non-obvious factor-of-2
correction for a non-polarized calculation, Section 5.2); fitting two
slopes; and finally re-deriving `U = 1/chi0 - 1/chi` correctly. This
workflow does all of that mechanically and consistently, and — just as
importantly — reports **when the underlying assumption (a clean, linear
response) doesn't hold** (Section 5.3), rather than silently handing you a
number.

## 2. Libraries and external dependencies used

- **`numpy`** — the linear fit itself (`numpy.polyfit`, degree 1) for both
  the `chi` and `chi0` branches, and the R² computed from its residuals.
- **`pymatgen` / `spglib`** (via `core/symmetry.py`) — only when
  `--atom-index` is ambiguous (Section 3.3): finds symmetry-inequivalent
  atoms of the target species and reports which Wyckoff site each
  candidate `--atom-index` belongs to, so you're not guessing blindly.
- **`matplotlib`** (optional, `--view`) — a live preview of both response
  branches (points + fitted lines) at Stage 3, so you can eyeball the fit
  instead of trusting R² alone.
- **`gnuplot`** (optional, `--save-gnuplot`) — writes the same response
  curve as a portable `.dat`/`.gplot` pair under `hubbardu_runs/plot/`,
  this workflow category's own older plotting convention (matching
  `stb-strainAnalysis`/`stb-elasticAnalysis`/`stb-convergenceAnalysis`).
- No network access, no ML/MACE dependency anywhere in this workflow —
  every stage is a pure, local, deterministic file operation or numerical
  fit.

## 3. Stage 1: reference prep (`stb-hubbardu`, code `4.7.1`)

### 3.1 What it does

Given a structure, a target species, and (optionally) a specific atom
index, `stb-hubbardu` writes one folder, `reference/`, containing:

- `structure.fdf` — your structure, copied in as-is (or with one atom's
  species aliased, Section 3.3).
- `calc.fdf` — your own template, untouched, except for one new line at
  the very top: `%include config_extra.fdf`.
- `config_extra.fdf` — the DFT+U *perturbation* block (`alpha=0.0` here,
  the reference point) plus the single-point-SCF-enforcement lines
  (Section 1.5) — kept in a sidecar file rather than inlined, so nothing
  in your own template is ever rewritten or duplicated.

It also snapshots the (possibly aliased) structure and the raw `calc.fdf`
template as `_template_structure.fdf`/`_template_calc.fdf`, and writes
`run_manifest.json` — the single source of truth Stages 2 and 3 both read
(species, shell quantum numbers, `J`, the SIESTA `SystemLabel`, and a
`runs` dict that grows with every stage).

### 3.2 The single-point-SCF constraint, live

```
$ stb-hubbardu --species Mn --no-intro
```

```
[2] DFT+U PERTURBATION SETUP
------------------------------------------------------------
SystemLabel            : siesta
Species / shell        : Mn (3d: n=3, l=2)
Perturbation           : alpha=0.0 (reference run)
Forced to single-point SCF via '%include config_extra.fdf'.

config_extra.fdf (written into every generated folder):
  # Auto-generated by stb-hubbardu -- DFT+U perturbation (Cococcioni-de Gironcoli) + single-point SCF enforcement.
  LDAU.PotentialShift T
  %block LDAU.proj
  Mn   1
  n=3    2
  0.000    0.000
  0.000    0.000
  %endblock LDAU.proj

  # Forces a pure single-point SCF (no ionic or cell relaxation) -- ...
  MD.TypeOfRun       CG
  MD.Steps           0
  MD.NumCGsteps      0
  MD.VariableCell    false
```

Note `n=3, l=2` — this is Mn's default correlated shell (`3d`), the exact
same `core/dftu_data.py::DEFAULT_SHELL`/`SHELL_NAMES` tables `stb-dftu`
(example `1.5`) already uses, so a shell picked automatically here means
the same thing it would mean there. `reference/calc.fdf` itself is just:

```
%include config_extra.fdf

SystemLabel siesta
%include structure.fdf
...
```

— your template, byte-for-byte, with one new line prepended.

### 3.3 Multi-atom species and `--atom-index`, live

`structure_2mn.fdf` in this folder is the conventional (2-atom) bcc-Mn
cell — space group `Im-3m` (No. 229), both atoms symmetry-equivalent.
Asking to perturb `Mn` here without saying *which* one is a hard error,
not a guess:

```
$ stb-hubbardu -s structure_2mn.fdf --species Mn --no-intro
```

```
[1] ATOM SELECTION
------------------------------------------------------------
Total atoms in structure : 2
Atom(s) of 'Mn' : 2 found (index(es) 1, 2)
[FAIL] --species 'Mn' appears 2 times in 'structure_2mn.fdf' -- the
linear-response perturbation must be isolated to a single atom (perturbing
every atom of the species at once measures a collective response, not the
isolated single-site response the method requires). Pass --atom-index to
pick one:
Space group: Im-3m (No. 229)
--atom-index | Wyckoff | Note
1            | a       | 2 symmetry-equivalent atom(s) -- any one of them gives the same result
```

Passing `--atom-index 1` resolves it — and aliases just that one atom:

```
$ stb-hubbardu -s structure_2mn.fdf --species Mn --atom-index 1 --no-intro
```

```
Perturbed atom : atom #1 of 'Mn', isolated via alias species 'Mn_pert'
(same Z, own pseudopotential) -- the other atom(s) of 'Mn' are left
unperturbed.
```

and the generated `reference/structure.fdf` shows exactly what changed —
one new species label, same atomic number, and only atom #1 repointed to
it:

```
NumberOfSpecies 2
%block ChemicalSpeciesLabel
1 25 Mn
 2   25   Mn_pert
%endblock ChemicalSpeciesLabel
...
%block AtomicCoordinatesAndAtomicSpecies
  0.00   0.00   0.00   2
0.50 0.50 0.50 1
%endblock AtomicCoordinatesAndAtomicSpecies
```

Atom #1 (originally species index `1`, i.e. `Mn`) now points at species
index `2` (`Mn_pert`); atom #2 is untouched, still plain `Mn`. `%block
LDAU.proj` in `config_extra.fdf` now perturbs `Mn_pert` only — exactly the
isolated single-site perturbation Section 1.6 explained, made concrete.

### 3.4 Pseudopotentials, live

`-p/--pseudo-dir` accepts a bundled bank (same `core/pseudopotentials.py`
banks every other tool in this suite shares) or a custom folder, and only
copies the species your structure actually needs:

```
$ stb-hubbardu --species Mn --pseudo-dir dojo --no-intro
```

```
[4] PSEUDOPOTENTIALS
------------------------------------------------------------
[OK] Copied 1 pseudopotential file(s) from '.../pseudopotentials/dojo'
into 'hubbardu_runs/reference' (and saved for stb-hubbarduAlphas to reuse)
```

The snapshot it saves is what lets Stage 2 copy the same pseudopotential
into every `scf_alpha_*`/`frozen_alpha_*` folder too, without you having to
pass `--pseudo-dir` a second time.

### 3.5 Report structure

```
[0] RUN METADATA
[1] ATOM SELECTION
[2] DFT+U PERTURBATION SETUP
[3] REFERENCE FOLDER & TEMPLATE SNAPSHOT
[4] PSEUDOPOTENTIALS
[5] SUMMARY & NEXT STEPS
```

`--save-report` persists this exact text to
`hubbardu_runs/hubbardu_stage1.txt`; no `references.bib` is written by any
of the three tools in this workflow (unlike `stb-inputfile`/`stb-kgrid`/
`stb-kpath`) — the Cococcioni & de Gironcoli citation is instead printed
directly in Stage 3's `[3] COMPUTED U` section, right next to the number
it justifies.

### 3.6 Running it both ways

**A — direct CLI**: every command above.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.7.1
```

`4.7.1` asks for the structure file, the calc.fdf template, the species,
an optional atom index, an optional shell override, `J`, a pseudopotential
source (bundled bank / custom path / skip), the output directory, and
finally whether to save a report — the exact same flags as the CLI, in
the exact same order. `example_4.7.sh` proves both paths agree.

## 4. Stage 2: perturbation prep (`stb-hubbarduAlphas`, code `4.7.2`)

### 4.1 What it does

Reads `reference/` and `run_manifest.json` from Stage 1, and writes one
folder per perturbation strength in `--alphas` (default
`-0.10 -0.05 0.05 0.10`, eV):

- `scf_alpha_<v>/` — full self-consistency, for the **screened** (`chi`)
  branch.
- `frozen_alpha_<v>/` — `MaxSCFIterations` capped (`--frozen-iterations`,
  default `1`), for the **bare** (`chi0`) branch — one per `alpha`
  **plus** `alpha=0.0`, evaluated by the exact same frozen recipe as every
  other point (not just re-using `reference`'s own fully-converged value),
  so every point on the frozen-response line is mutually consistent with
  every other.

### 4.2 Auto-copying the reference DM, live

```
$ stb-hubbarduAlphas --dir hubbardu_runs --no-intro
```

```
[1] VALIDATION & REFERENCE CHECK
------------------------------------------------------------
Species / shell        : Mn (n=3, l=2), J=0.000 eV
Reference DM           : hubbardu_runs/reference/siesta.DM (found, will be copied into every folder)
Template snapshot      : OK (no pre-existing DFT+U block)
```

If `reference/siesta.DM` doesn't exist yet, Stage 2 refuses to proceed —
by design: without it, every perturbed folder would start from its own
independent SCF guess, defeating the entire point of Section 1.5's
fixed-starting-density requirement. (This walkthrough hasn't actually run
SIESTA, so a placeholder file stands in for it below — see Section 6 for
why that's fine for proving the mechanics and the analysis math, but not a
substitute for a real converged density.)

### 4.3 `scf_*` vs. `frozen_*` folders

```
$ cat hubbardu_runs/scf_alpha_0.0500/config_extra.fdf
```
```
LDAU.PotentialShift T
%block LDAU.proj
Mn   1
n=3    2
0.050    0.000
0.000    0.000
%endblock LDAU.proj
DM.UseSaveDM T
MD.TypeOfRun       CG
...
```

```
$ cat hubbardu_runs/frozen_alpha_0.0500/config_extra.fdf
```
```
LDAU.PotentialShift T
%block LDAU.proj
Mn   1
n=3    2
0.050    0.000
0.000    0.000
%endblock LDAU.proj
MaxSCFIterations 1
DM.UseSaveDM T
MD.TypeOfRun       CG
...
```

Identical perturbation, identical starting density (`DM.UseSaveDM T`) —
the *only* difference is `MaxSCFIterations 1` in the frozen folder. That
one line is the entire mechanism separating the bare response from the
screened one: one SCF iteration is barely enough for the perturbed
orbital's own occupation to react, nowhere near enough for the rest of the
density to relax around it.

### 4.4 Report structure

```
[0] RUN METADATA
[1] VALIDATION & REFERENCE CHECK
[2] PSEUDOPOTENTIALS
[3] PERTURBATION FOLDERS
[4] SUMMARY & NEXT STEPS
```

### 4.5 Running it both ways

**A — direct CLI**: `stb-hubbarduAlphas --dir hubbardu_runs`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.7.2
```

`4.7.2` asks for the stage-1 directory, the perturbation strengths
(space-separated, blank = default), the frozen-iteration cap, and whether
to save a report. `example_4.7.sh` proves both paths agree.

## 5. Stage 3: analysis (`stb-hubbarduAnalysis`, code `4.7.3`)

### 5.1 Report structure

```
[0] RUN METADATA
[1] READING RUNS
[2] LINEAR RESPONSE FIT
[3] COMPUTED U
[4] OUTPUT FILES
[5] SUMMARY & NEXT STEPS
[6] LIBRARY WARNINGS
```

`[1]` prints one table row per folder (kind, alpha, parsed occupation, SCF
status); `[2]` the two fitted slopes/intercepts/R²; `[3]` the formula and
the final `U`, plus a same-element literature-value comparison (from
`core/dftu_data.py::REFERENCE_U`, the same table `stb-dftu
--list-reference`/`--suggest` expose) as an *informational* sanity check
only — never a validation, since that table is itself just a
Materials-Project-tabulated GGA+U convenience value for oxides, not a
ground truth for your exact functional/pseudopotential/basis.

### 5.2 The non-polarized occupation-doubling fix

SIESTA's own `Occupations:` line means something different depending on
whether the calculation is spin-polarized. For a **non-polarized**
(`nspin=1`) run, SIESTA prints only the *per-spin-channel* occupation —
literally half the real total, because the code internally divides by
`(3 - nspin)`. `stb-hubbarduAnalysis` detects this from the *shape* of the
parsed line (exactly 2 numbers, vs. 3+ for a polarized run) and doubles it
back:

```python
total = float(nums[-1])
if len(nums) == 2:
    total *= 2.0
```

This was verified against a real non-polarized Mn run: a raw printed value
of `~2.47` is the correct `~4.94` total once doubled — silently trusting
the raw number would have measured a response curve at roughly *half* the
real occupation, propagating straight into a wrong `chi`/`chi0` and hence a
wrong `U`. The worked example in Section 6 deliberately uses a
3-number (already-total) format specifically so its analytic ground truth
stays simple to verify by hand; a real non-polarized SIESTA run would hit
this doubling path automatically and correctly.

### 5.3 Sanity checks: R², intercept deviation, and sign checks

A single number, `U`, can look confident and still be meaningless if the
occupation-vs-`alpha` data isn't actually linear. Three independent checks
guard against that, all demonstrated live in Section 6.3:

- **R² per branch** (default tolerance `0.98`, `--r2-tolerance`): warns
  when either fit doesn't explain the data well — usually a sign the
  `--alphas` range in Stage 2 was too wide for the linear regime.
- **Intercept deviation**: each fit's *intercept* (extrapolated occupation
  at `alpha=0`) is compared against the occupation *actually measured* at
  the real `alpha=0` point; a difference above `0.01` warns that the
  response may not be as linear as R² alone suggests — R² can look
  deceptively good even when the intercept has drifted.
- **Sign / negativity checks**: `U` is a repulsive on-site term and must
  be positive; a negative result, or `chi`/`chi0` with opposite signs, is
  flagged as a strong indicator of a bad fit rather than a real (if
  surprising) physical result.

### 5.4 `[6] LIBRARY WARNINGS`

`numpy.polyfit` silently emits a `RankWarning` when the fit is badly
conditioned (near-duplicate `alpha` values, or too few points for the
requested degree). Rather than letting that warning interleave with the
report mid-print — confusing at best — `stb-hubbarduAnalysis` traps it via
`core/cli.py::capture_library_noise` and reports it in its own dedicated,
always-present final section, `No library warnings.` when the fit was
well-conditioned (as in every case in this walkthrough).

### 5.5 Running it both ways

**A — direct CLI**: `stb-hubbarduAnalysis --dir hubbardu_runs`.

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 4.7.3
```

`4.7.3` asks for the stage-1/2 directory, the SIESTA output filename
inside each folder, the R² tolerance, an optional output filename, then
three yes/no questions (`--save-gnuplot`, `--view`, `--save-report`).
`example_4.7.sh` proves both paths agree.

## 6. Worked example: a known analytic ground truth, end-to-end

Real SIESTA output isn't available inside this walkthrough (there's no
SIESTA binary to invoke), so — exactly like `stb-eosAnalysis`'s own
verification against a synthetic exact Birch-Murnaghan curve — Section 6
fabricates `calc.out` files whose occupation-vs-`alpha` response follows
an **exact, hand-chosen linear law**, so the correct answer is known
*before* running Stage 3, not just plausible after the fact:

- Reference occupation `N0 = 5.0`.
- Self-consistent (screened) slope `chi = 0.30` 1/eV: `n(alpha) = 5.0 + 0.30 * alpha`.
- Frozen-density (bare) slope `chi0 = 0.20` 1/eV: `n(alpha) = 5.0 + 0.20 * alpha`.
- Expected `U = 1/chi0 - 1/chi = 1/0.20 - 1/0.30 = 5.0 - 3.333333 = 1.666667 eV`.

The script writes these nine occupation values (reference + 4 `scf_alpha_*`
+ 5 `frozen_alpha_*`, matching Stage 2's own default `--alphas`) as real
`Occupations:   x   x   TOTAL` lines into the exact folders Stage 1/2 just
generated, then runs Stage 3 for real:

```
[2] LINEAR RESPONSE FIT
------------------------------------------------------------
Branch                     | chi (1/eV) | Intercept | R^2
Self-consistent (screened) | 0.300000   | 5.000000  | 1.0000
Frozen-density (bare)      | 0.200000   | 5.000000  | 1.0000

[3] COMPUTED U
------------------------------------------------------------
Formula                : U = 1/chi0 - 1/chi (Cococcioni & de Gironcoli, PRB 71, 035105, 2005)
Computed U             : 1.6667 eV
```

Both fits recover the exact input slopes (R² = 1.0000, since the data is
exactly linear by construction) and `U` matches the hand-computed
`1.666667 eV` to the 4th decimal — confirming the entire pipeline (folder
scan → occupation parsing → linear fit → `U` formula) end to end, the same
way `stb-eosAnalysis`'s own synthetic-EOS verification confirms *its*
fitting pipeline.

### 6.1 What a bad fit looks like: negative `U`

Swapping which branch is steeper (`chi = 0.20`, `chi0 = 0.30` — the bare
response now *larger* than the screened one) reproduces the exact same
formula honestly, and it comes out negative:

```
[3] COMPUTED U
------------------------------------------------------------
[WARNING] Computed U is NEGATIVE -- the Hubbard U is a repulsive on-site
term and should always be positive. This strongly suggests a bad fit or
bad data rather than a real result.
Computed U             : -1.6667 eV
```

Nothing here is a bug: the arithmetic is correct, and the tool tells you,
plainly, that a negative on-site repulsion isn't physical — a real
material would never converge to this sign combination in practice.

### 6.2 What a bad fit looks like: low R² and intercept drift

Replacing just one of the five `scf_alpha_*` occupations with an outlier
(everything else left exactly on the same line as Section 6's clean
dataset) triggers both remaining safety nets at once:

```
[2] LINEAR RESPONSE FIT
------------------------------------------------------------
Self-consistent (screened) | 0.870000   | 4.943000  | 0.2497
[WARNING] The self-consistent response fit has R^2=0.2497, below the 0.98
tolerance -- the perturbation range may be too wide for the linear regime.
[WARNING] The self-consistent fit's intercept (4.943000) differs from the
actually-measured alpha=0 occupation (5.000000) by more than 0.01 -- the
response may not be as linear as the R^2 alone suggests.
```

One bad SCF point (a common real-world failure: a run that converged to
the wrong local minimum, or simply didn't converge) is enough to wreck the
fit — exactly why Stage 3 checks this instead of reporting `chi`/`U`
unconditionally.

## 7. Known, deliberate limitations

- **No SIESTA run is ever performed by this suite.** All three stages
  generate/read input and output files; you run SIESTA in `reference/`
  and every `scf_alpha_*`/`frozen_alpha_*` folder yourself, in between
  stages, exactly as the real workflow requires.
- **`J` is fixed for the whole workflow, not fit.** Only `U` comes out of
  the linear-response fit; `J` (default `0.0`) is whatever you passed to
  Stage 1 throughout.
- **One species/one atom per run.** Computing `U` for two different
  species (e.g. a spinel with two transition-metal sites) means running
  the whole three-stage workflow twice, independently, into two separate
  `--output-dir`/`--dir` folders.
- **The `REFERENCE_U` comparison in `[3] COMPUTED U` is informational
  only** — a literature GGA+U oxide value from a different code, functional,
  pseudopotential, and basis set is a sanity check on the right order of
  magnitude, never a validation of your own computed value.
- **The linear-response approximation itself assumes small `alpha`.**
  Section 6.2 shows what happens when a single point is inconsistent with
  linearity; a genuinely too-wide `--alphas` range (Section 1.4) would
  show the same symptom (low R²) even with otherwise-perfect SCF
  convergence in every folder.

## 8. Step-by-step: running this workflow on your own structure

1. Have a structure (`.fdf`) and a working `calc.fdf` template — one that
   does **not** already configure `%block LDAU.proj`/`DFTU.*` (Stage 1
   refuses to proceed on a template that does, to avoid silently stacking
   perturbations on top of an existing DFT+U setup).
2. Run Stage 1: `stb-hubbardu --species <El> [--atom-index N] [--shell
   3d|4d|5d|4f|5f] [--j <val>] [--pseudo-dir <bank-or-path>]`. If your
   species appears more than once, follow the printed Wyckoff-site table
   to choose `--atom-index`.
3. Run SIESTA yourself inside `hubbardu_runs/reference/`. Confirm it
   converges (writes a `.DM` file) — if it doesn't, `stb-hubbardu --help`
   itself lists concrete `calc.fdf` settings to try (a lower
   `DM.MixingWeight`, higher `DM.NumberPulay`, looser `SCF.H.Tolerance`).
4. Run Stage 2: `stb-hubbarduAlphas --dir hubbardu_runs [--alphas ...]
   [--frozen-iterations N]`.
5. Run SIESTA yourself inside every generated `scf_alpha_*/` and
   `frozen_alpha_*/` folder.
6. Run Stage 3: `stb-hubbarduAnalysis --dir hubbardu_runs [--save-gnuplot]
   [--view] [--save-report]`. Read `[2]`/`[3]` for any `[WARNING]` before
   trusting the computed `U` — don't just read the final number.
7. Copy the `%block LDAU.proj` from the written `<species>_LDAU.fdf` into
   your **production** `calc.fdf` (this file has no
   `LDAU.PotentialShift` line — it's the real correction, not another
   perturbation probe).

## 9. Files in this folder

| File | Purpose |
|---|---|
| `structure.fdf` | Single-atom simplified bcc-Mn cell — the main Stage 1/2 mechanics demo. Not a physically converged structure, purely for exercising the tools. |
| `structure_2mn.fdf` | Conventional 2-atom bcc-Mn cell (space group `Im-3m`), both atoms symmetry-equivalent — demonstrates the `--atom-index`/species-aliasing requirement (Section 3.3). |
| `calc.fdf` | Shared `calc.fdf` template for both structures above (non-polarized, modest basis/mesh — mechanics only). |
| `example_4.7.sh` | The guided walkthrough (**not** an automated test — see `test/4-workflow/7-hubbardu/{prep,alphas,analysis}/test.sh` for that). Pauses between sections so you can read before moving on; safe to re-run. |
| `output/` | Created by `example_4.7.sh` when you run it (git-ignored, not checked in). See below. |

## 10. Running the script

```bash
./example_4.7.sh
```

| Case | Command(s) | What it shows |
|---|---|---|
| `output/stage1/` | `stb-hubbardu --species Mn -o output/stage1` | Stage 1 mechanics, single-atom structure |
| `output/stage1_2mn/` | `stb-hubbardu -s structure_2mn.fdf --species Mn [--atom-index 1] -o output/stage1_2mn` | `--atom-index`/aliasing, both the ambiguous-error and resolved paths |
| `output/stage1_pp/` | `stb-hubbardu --species Mn --pseudo-dir dojo -o output/stage1_pp` | pseudopotential bank copy |
| `output/workflow/` | Stage 1 → Stage 2 → (synthetic `calc.out` injection) → Stage 3 | the full 3-stage chain, known-answer worked example (Section 6) |
| `output/workflow_badU/`, `output/workflow_noisy/` | same chain, deliberately inverted/noisy occupations | the warning paths (Section 6.1/6.2) |
| *(no folder — a diff only)* | Stages 1/2/3 via `printf … \| stb-suite` | proof the interactive menu (`4.7.1`/`4.7.2`/`4.7.3`) agrees with the CLI |

## What's next

- **`1.5-stb-dftu`** — once you have a computed `U`, this is what you'd
  reach for again to regenerate the exact same `%block LDAU.proj` snippet
  from that number for a *different* structure/species combination, or to
  compare it against the literature `REFERENCE_U` table this workflow's
  own Stage 3 already checks against informally.
- **`1.1-stb-inputfile`** — the production `calc.fdf` this `U` ultimately
  belongs in; paste the `%block LDAU.proj` from `<species>_LDAU.fdf`
  straight into it (no `LDAU.PotentialShift` line to remove — Stage 3
  already omits it).
- Every other `4.x` workflow example in this suite generates real
  candidate/perturbed geometries and expects a real SIESTA run in between
  stages, the same two-way (direct CLI / interactive menu) split, and the
  same numbered-report convention — if this is your first workflow
  example, `4.1-strain` is the shortest two-stage one to start with.
