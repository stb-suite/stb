# 4.6 — Workflow: Structure Solution (XRD) (`stb-xrdsearch` / `stb-xrdrank`)

This workflow has 2 stages: **Stage 1** (`stb-xrdsearch`, code `4.6.1`) takes
a *composition* (not a structure — this is the one workflow pair in this
suite that doesn't start from an already-built `structure.fdf`) and a list
of candidate space groups, and casts one or more random atomic arrangements
per group, ready to relax with real SIESTA. **Stage 2** (`stb-xrdrank`, code
`4.6.2`) simulates each candidate's own powder XRD pattern and ranks them by
similarity to a real experimental pattern — a fast pre-screen for *which*
candidate space group/arrangement is most likely the true structure, before
committing to expensive DFT on all of them.

This answers a genuinely different question from `4.1`-`4.5`: those all
start from a structure you already believe is right and ask how it responds
to load (`4.1`/`4.2`), how strongly it's bound (`4.3`), how it vibrates
(`4.4`), or whether its total energy/relaxed geometry has converged with
some numerical parameter (`4.5`) — this workflow instead asks **which
structure is right in the first place**, using a real experimental
measurement as the arbiter.

Both stages live in this one folder and this one tutorial, same reasoning
as `4.1`-`4.5`: Stage 1's output only exists to feed Stage 2.

## 1. Theory

### 1.1 Powder XRD and Bragg's law

A powder sample contains crystallites in every orientation at once, so a
monochromatic X-ray beam (this workflow defaults to Cu-K-alpha,
`1.54184 Ang`) diffracts off every family of lattice planes `(hkl)`
simultaneously, each at its own angle given by **Bragg's law**:

```
n * lambda = 2 * d_hkl * sin(theta)
```

`d_hkl` is the spacing between `(hkl)` planes (for a cubic cell,
`d_hkl = a / sqrt(h^2 + k^2 + l^2)`). The experiment records intensity vs.
`2*theta` — a 1D "fingerprint" of the 3D crystal structure: peak
*positions* come from the lattice geometry (Bragg's law above), peak
*intensities* come from the structure factor (1.2). This is exactly why a
1D pattern can distinguish between different atomic arrangements sharing
the same composition: different space groups give different lattice
geometries and different structure factors, hence different fingerprints.

### 1.2 The structure factor and systematic absences

Not every `(hkl)` allowed by Bragg's law actually produces a peak — the
**structure factor** `F_hkl = sum_j f_j * exp(2*pi*i*(h*x_j + k*y_j +
l*z_j))` (sum over every atom `j` in the basis, `f_j` its atomic scattering
factor) can be exactly zero for some `(hkl)` combinations, depending on the
space group's own symmetry — **systematic absences**. This is the concrete
mechanism this workflow's decisive example (Section 5) exploits: diamond-
cubic silicon (`Fd-3m`, space group 227, 8 atoms/cell at Wyckoff `8a`) and a
naively-relaxed simple-cubic arrangement of the same 8 Si atoms have
*different* systematic absences, so their simulated patterns look nothing
alike even though both are "just silicon."

### 1.3 From a stick pattern to a comparable profile

A calculated pattern (via `pyxtal.XRD.XRD`, Section 2) is naturally a set
of discrete sticks (`2theta, d, h, k, l, intensity`) — zero width, since
Bragg's law gives an exact angle for a perfect infinite crystal. Real
instruments broaden every peak (finite crystallite size, instrumental
resolution, ...), so before comparing to a real pattern, both this
workflow's own simulated stick pattern (via `pattern.get_profile()`) and
any experimental input that looks like a similarly sparse stick/peak list
(Section 1.5/4.4) are Gaussian-broadened onto the same footing:
`FWHM = 0.1 deg`, `res = 0.01 deg` — pyxtal's own defaults, used
consistently everywhere this suite touches XRD (`stb-xrd`,
`core/xrd.py::broaden_peak_list`).

### 1.4 The similarity metric

`stb-xrdrank` scores every candidate with `pyxtal.XRD.Similarity`, a
cosine-weighted cross-correlation between the two broadened profiles `f`
(simulated) and `g` (experimental):

```
value = | xCorr_fg(w) / sqrt( xCorr_ff(w) * xCorr_gg(w) ) |
```

(`xCorr_fg(w)` is a windowed cross-correlation, weighted by a cosine kernel
over a shift range `l` — this is what makes the metric tolerant of a small
rigid shift between the two patterns, e.g. from a slightly-off lattice
constant, rather than demanding a perfect point-by-point match.) Read
literally: **this is a normalized dot product** — `1.0` means the two
profiles are identical up to scale, `~0` means uncorrelated. Section 4.5
traces the exact place this formula breaks (a real bug this session found
and fixed) when one of the two profiles is degenerate (all zero).

### 1.5 The experimental data file: the model this workflow expects

**This is the canonical model** — `experimental.dat` in this folder is a
real, working example of it:

```
# 2theta   intensity
   28.330    100.000
   47.110     47.434
   55.890     24.408
   ...
```

Two columns. Column 1 is `2theta` in degrees; column 2 is intensity, on
**any** scale (`0-100`, `0-1`, raw counts, ...) — both sides get
max-normalized before comparing, so the absolute scale never matters.
Whitespace- or comma-separated; blank lines and `#`-prefixed comments
(the whole header above) are skipped. `--wavelength` must match whatever
source actually produced this data (default `CuKa`).

**How the code actually picks the intensity column, and why it isn't just
"always the last one":** `core/xrd.py::read_experimental_pattern` reads
*every* column on each line (not just the first and last), then
`_select_intensity_column()` picks the **last column that is not constant**
among the columns after 2-theta. For a plain 2-column file that's the only
candidate anyway. This matters for real instrument export formats with a
3rd (or more) column — Section 4.5 works through a real one, live, and
shows exactly what used to go wrong before this session fixed it.

### 1.6 Why this needs a *search* at all

Given only a composition and a space group, there are usually many
symmetry-equivalent-but-not-identical ways to place atoms on that group's
Wyckoff positions (and, for a low-symmetry group, real continuous
positional freedom within them) — `stb-crystalcast` (driven internally by
`stb-xrdsearch`, Section 3) samples this space randomly, `--count-per-group`
times per group, which is why Stage 1 is a *search* (possibly multiple
candidates per group) feeding a *ranking* (Stage 2), not a single
deterministic answer.

## 2. Libraries and external dependencies used

- **[PyXtal](https://pyxtal.readthedocs.io/)** (`1.1.4` in this
  environment) does essentially all of the domain-specific work: random
  symmetry-constrained structure casting from a space group + composition
  (`stb-crystalcast`, driven by Stage 1), the powder-pattern simulation
  itself (`pyxtal.XRD.XRD`), and the similarity metric (`pyxtal.XRD.
  Similarity`, Section 1.4). `core/xrd.py` wraps the pieces shared between
  `stb-xrd` (Analysis, code `3.9`) and this workflow's Stage 2.
- **[pymatgen](https://pymatgen.org/)** (`2025.10.7`) supplies the
  `Structure` object every other piece of code here operates on, plus
  `SpacegroupAnalyzer` (via **spglib**, `2.7.0`) for the space-group label
  shown in Stage 2's ranking table (`core/symmetry.py::space_group_label`).
- **[ASE](https://wiki.fysik.dtu.dk/ase/)** (`3.27.0`) reads a relaxed
  SIESTA `.STRUCT_OUT` file (`core/structure_io.py::read_siesta_structure`)
  — there is no `.fdf`-style reader for that format, it's ASE's own.
- **[MACE](https://github.com/ACEsuit/mace)** (via `core/mace_relax.py`,
  the same shared module `stb-mlrelax`/`stb-defect --ml-rank`/`stb-
  amorphize` use) — **optional**, only touched by Stage 1's `--mace-relax`
  (Section 3.4) and `--ml-rank` (a *different*, transient, ranking-only
  relax passed through to `stb-crystalcast` — Section 3.4 explains the
  distinction). Needs `pip install stb_suite[ml]`.
- **numpy** / **scipy** — general array math; `scipy.signal.find_peaks`
  is used only by *this tutorial's own* fixture-generation, not by the
  suite itself.
- **matplotlib** (`3.10.8`) — Stage 2's `--view` (Section 4.6), a two-panel
  figure (ranking bar chart + best-match overlay), always via a blocking
  `plt.show()` after the report is fully printed (same convention as
  `stb-xrd`/`stb-density`/`stb-workfunction`'s own `--view`).
- **gnuplot** (`6.0`, external, not a Python dependency) — only touched
  when `--save-gnuplot` is given (Section 4.6); the suite writes `.dat` +
  `.gplot` file pairs, gnuplot itself is never invoked by the suite.

## 3. Stage 1: casting candidates (`stb-xrdsearch`, code `4.6.1`)

### 3.1 What it does, and how

For each space group in `--groups`, shells out to `stb-crystalcast` once
(`--count-per-group` random candidates per call), then writes each surviving
candidate into its own `<output-dir>/group_<G>[_<i>]/structure.fdf` — a bare
structure file, no calc setup yet. From there, **by default**, every
candidate folder is made ready to run as-is:

1. **calc.fdf** — auto-generated (the exact same standard relaxation
   template `stb-inputfile -t relax` produces, reused directly as a library
   call — Section 3.2), or copied verbatim from `--calc <path>` if you
   already have your own template.
2. **Pseudopotentials** — copied in from `-p/--pseudo-dir` (a bundled bank
   name, e.g. `dojo`, or your own folder — Section 3.3), for the real
   species found in the *cast* structure (not `--species` directly, which
   `--molecular` can make a molecule name rather than an element symbol).
3. **MACE pre-relaxation** (opt-in, `--mace-relax`) — positions + cell,
   before calc.fdf is generated, so the k-grid calc.fdf computes already
   reflects the pre-relaxed cell (Section 3.4).

### 3.2 calc.fdf: auto-generate vs. custom template, live

```
$ stb-xrdsearch --species Si --num-ions 8 --groups 227 --no-intro

[1] CALC.FDF SETUP
------------------------------------------------------------
Mode                   : auto-generate (standard relaxation calc.fdf, same as stb-inputfile -t relax)
DFT-D3 correction      : disabled
Spin polarization      : non-polarized
```

`--calc my_template.fdf` switches to the second mode instead — copied
verbatim into every candidate folder as `calc.fdf` (always this exact name
in *either* mode, even though a custom template's own basename might
differ, so every folder in a batch is runnable the identical way
regardless of which mode produced it).

### 3.3 Pseudopotentials, live

```
[4] PSEUDOPOTENTIALS
------------------------------------------------------------
Source                 : .../stb-suite/src/stb/pseudopotentials/dojo
Species | Status
----------------
Si      | found
[OK] All required pseudopotentials found -- will be copied into every generated folder.
```

Copied (`shutil.copy2`), not symlinked — a symlink into a bundled bank or
an external `--pseudo-dir` breaks the moment a candidate folder is
archived/rsynced/scp'd elsewhere (e.g. to an HPC cluster) without also
bringing the link target; a real copy always travels with the folder
(`core/pseudopotentials.py`, shared by every prep tool in this suite).

### 3.4 MACE pre-relaxation, live

```
$ stb-xrdsearch --species Si --num-ions 8 --groups 227,225,229 --mace-relax -p dojo --no-intro

[3] MACE PRE-RELAXATION
------------------------------------------------------------
Model                  : MACE-MP-0 (small)
Folder    | E before (eV) | E after (eV) | Vol change | Max disp (Ang) | Steps
--------------------------------------------------------------------------------------
group_227 | -42.9479      | -42.9666     | +1.91%     | 0.0444         | 7 (converged)
group_225 | -38.7782      | -41.4365     | -24.48%    | 0.7294         | 9 (converged)
group_229 | -38.3101      | -41.4366     | +27.56%    | 0.5867         | 7 (converged)
```

Real, physically meaningful result, not illustrative: `227` (`Fd-3m`,
diamond) relaxes to `-42.9666 eV`; `225` and `229` (both physically
*wrong* space groups for elemental Si — a single-site FCC/BCC packing,
not the tetrahedrally-bonded diamond motif) both collapse under their own
imposed symmetry constraint into the same higher-energy `-41.44 eV`
minimum — **already, before any XRD comparison at all, MACE's own energy
shows a ~1.5 eV gap favoring the true structure.** This is *not* the same
thing as `--ml-rank`: that flag (passed straight through to
`stb-crystalcast`) relaxes candidates positions-only, transiently, purely
to *select/rank* which of several random attempts to keep — the relaxed
geometry it computes is discarded, never written to disk. `--mace-relax`
relaxes (positions + cell) the *final* candidates that survive casting,
and **keeps** the result, overwriting `structure.fdf` — the two are
independent and can be combined.

### 3.5 Report structure

```
[0] RUN METADATA
[1] CALC.FDF SETUP
[2] CANDIDATE CASTING
[3] MACE PRE-RELAXATION
[4] PSEUDOPOTENTIALS
[5] CALC.FDF & PSEUDOPOTENTIAL GENERATION
[6] SUMMARY & NEXT STEPS
[7] LIBRARY WARNINGS
```

`[3]`-`[5]` always print, even when their feature wasn't requested (e.g.
`[3]` prints "Not requested" rather than being skipped) — every section
number means the same thing on every run, whether or not you used that
particular flag.

### 3.6 Running it both ways

Both CLI and `stb-suite` -> `4.6.1` produce identical `xrd_search/` output
for identical answers (Section 9's script exercises the CLI path; the
interactive menu asks the same questions in the same order, including the
MACE pre-relax model choice, Section 10's script for `4.4` shows the
general pattern this suite uses for every workflow's interactive twin).

## 4. Stage 2: ranking candidates (`stb-xrdrank`, code `4.6.2`)

### 4.1 Report structure

```
[0] RUN METADATA
[1] CANDIDATES FOUND
[2] EXPERIMENTAL PATTERN
[3] RANKING
[4] OUTPUT DATA & PLOTS
[5] SUMMARY & NEXT STEPS
[6] LIBRARY WARNINGS
```

### 4.2 `--input-dir` defaults to `xrd_search`

`--input-dir` is **not required** — it defaults to `xrd_search`, Stage 1's
own default `--output-dir`, so the common case (Stage 2 right after Stage
1, same folder) needs zero extra flags:

```
$ stb-xrdrank --experimental experimental.dat --no-intro
...
Input directory        : xrd_search
```

If a candidate subfolder also contains a `*.STRUCT_OUT` file (you've run
real SIESTA there since Stage 1), that relaxed structure is used instead
of the raw `structure.fdf` — no flag needed, `find_candidates()` checks
for it automatically per subfolder (this is exactly what happens for the
real production data in Section 5).

### 4.3 The experimental pattern, live — a clean structure solution

Using this folder's own `experimental.dat` (real Si, Cu-K-alpha, the exact
2-column model Section 1.5 describes) against the 3 candidates cast in
Section 3.4 (`227` correct, `225`/`229` wrong):

```
$ stb-xrdrank --experimental experimental.dat --no-intro

[2] EXPERIMENTAL PATTERN
------------------------------------------------------------
Points                 : 11900
2-theta span (file)    : 1.000 - 119.990 deg
Peak-list detected     : no -- compared as a continuous scan, unbroadened.

[3] RANKING
------------------------------------------------------------
  -> group_225: similarity 0.0326 (Si, Pm-3m (No. 221), raw)
  -> group_227: similarity 0.9972 (Si, Fd-3m (No. 227), raw)
  -> group_229: similarity 0.0376 (Si, Pm-3m (No. 221), raw)
Rank | Name      | Similarity | Formula | Space group     | Source
------------------------------------------------------------------
1    | group_227 | 0.9972     | Si      | Fd-3m (No. 227) | raw
2    | group_229 | 0.0376     | Si      | Pm-3m (No. 221) | raw
3    | group_225 | 0.0326     | Si      | Pm-3m (No. 221) | raw
Score gap (#1 vs #2)   : 0.9596
```

A textbook-clean decisive result — `0.9972` for the true structure vs.
`~0.03-0.04` for two physically wrong ones, exactly the kind of
unambiguous separation Section 1.2's systematic-absences argument predicts.
Compare this against Section 5's real production data, where all 3
candidates are the *same* correct structure at slightly different
relaxation quality — a genuinely different, subtler question the score gap
also answers correctly.

### 4.4 Peak-list auto-detection and broadening — the fix, demonstrated live

`experimental_peaks.dat` (this folder) is the *same* real Si pattern, but
as a **sparse 9-point peak list** — the shape a literature/indexed peak
table or `stb-xrd`'s own stick-pattern output has, not a dense continuous
scan. `core/xrd.py::looks_like_peak_list()` auto-detects this (highly
uneven spacing between sorted 2-theta points) and Gaussian-broadens it
(Section 1.3) before comparing:

```
$ stb-xrdrank --experimental experimental_peaks.dat --no-intro
[2] EXPERIMENTAL PATTERN
Points                 : 9000
Peak-list detected     : yes -- auto-broadened (Gaussian, FWHM=0.1 deg) before comparing...
[3] RANKING
  -> group_227: similarity 0.9765   (correct structure -- still wins decisively)
  -> group_229: similarity 0.0159
  -> group_225: similarity 0.0138
```

**Now the same file, `--raw-experimental` (skip broadening):**

```
$ stb-xrdrank --experimental experimental_peaks.dat --raw-experimental --no-intro
[2] EXPERIMENTAL PATTERN
Points                 : 9
Peak-list detected     : yes -- but --raw-experimental was given, so it is compared unbroadened...
[3] RANKING
  -> group_225: similarity 0.3485
  -> group_229: similarity 0.3476
  -> group_227: similarity 0.2897   (correct structure -- now LOSES)
Score gap (#1 vs #2)   : 0.0009  [close call -- consider relaxing both before deciding]
```

**This is a live, real demonstration of why the fix matters**, not a
synthetic illustration: comparing a zero-width stick list directly against
a broadened simulated profile via cubic interpolation (what
`pyxtal.XRD.Similarity` does internally with no broadening) draws a
spurious curve through the empty gaps between peaks — the *wrong*
structures happen to score higher, a complete inversion of the correct
ranking. Broadening both sides the same way before comparing is not a
cosmetic default; without it, this workflow gives the wrong answer.

### 4.5 The 3-column bug — reproduced and shown fixed, live

Real diffractometer export formats aren't always 2 columns. A `.int`-style
export this session encountered in real production use
(`/home/carlos/test/Si/xrd.int`, Section 5) has **3** columns: `2theta`,
intensity, and a 3rd field (an unpopulated uncertainty/esd column) that was
**exactly zero on all 11900 lines**. `read_experimental_pattern` used to
always take the *last* column as intensity (a rule that happened to be
correct for the only 2 formats seen before: a plain 2-column file, and
`stb-xrd`'s own 6-column stick-pattern output) — for this file, that meant
reading the all-zero 3rd column as "intensity."

Traced to the exact line: `pyxtal.XRD.similarity_calculate()` computes
`xCorrfg_w / sqrt(aCorrff_w * aCorrgg_w)` (Section 1.4) — with an all-zero
`gy` (the misread "experimental" intensity), `aCorrgg_w = 0` and
`xCorrfg_w = 0` too, so the division is a literal `0/0` -> **`nan`**, with
no exception and no visible cause beyond a bare `RuntimeWarning: invalid
value encountered in scalar divide` — which, before this session's own
`[6] LIBRARY WARNINGS` section existed, had nowhere organized to surface at
all (Section 4.7).

**The fix, reproduced live in this folder** — same real data
(`experimental.dat`) with a synthetic bogus 3rd column appended, mimicking
the exact real-world mistake:

```
$ awk '{print $1, $2, "0.00000"}' experimental.dat > experimental_3col.dat
$ stb-xrdrank --experimental experimental_3col.dat --no-intro
[3] RANKING
  -> group_227: similarity 0.9972   (identical to the clean 2-column result, Section 4.3)
```

`_select_intensity_column()` now reads every column, excludes any column
that's constant across the whole file (a constant column can never be a
real diffraction intensity signal), and picks the **last non-constant**
one — correct for this 3-column case, and still correct (unchanged
behavior) for the 2-column and 6-column stick-pattern formats, since none
of *their* columns are ever constant for a real multi-peak pattern. If
*every* column after 2-theta is degenerate (nothing usable at all), the
tool now raises a clear error instead of silently proceeding:

```
$ printf "1 0\n2 0\n3 0\n4 0\n5 0\n" > broken.dat
$ stb-xrdrank --experimental broken.dat --no-intro
[FAIL] 'broken.dat': every column after 2-theta is constant (no usable intensity signal) -- check the file's column order/format.
```

This same fix also benefits `stb-xrd --compare-to` (Analysis, code `3.9`),
which shares the exact same `core/xrd.py::read_experimental_pattern`.

### 4.6 Plotting: gnuplot and matplotlib, both opt-in

`--save-gnuplot` writes **two** `.dat`+`.gplot` pairs, both under their own
`plot/` subfolder (kept separate from the plain-text ranking/report files
next to it):

```
$ stb-xrdrank --experimental experimental.dat --save-gnuplot --no-intro
[4] OUTPUT DATA & PLOTS
[OK] Ranking bar chart written to './plot/xrd_rank.dat' / './plot/xrd_rank.gplot'.
[OK] Best-match overlay written to './plot/xrd_rank_top_sim.dat' / './plot/xrd_rank_top_exp.dat' / './plot/xrd_rank_overlay.gplot'.
```

- **`xrd_rank.gplot`** — a bar chart, one bar per candidate, `similarity`
  on the y-axis.
- **`xrd_rank_overlay.gplot`** — the best candidate's own simulated pattern
  overlaid on the experimental one (two separate small `.dat` files, since
  they're generally on different grids — pyxtal's own broadened-profile
  resolution for one, the experimental file's own sampling for the other).

`--view` shows the **same two plots**, as one matplotlib figure (bar chart
on top, overlay below) instead of two `.gplot` scripts — a live preview,
not a file output, always shown **last**, after the report is fully
printed/closed, so a blocking `plt.show()` window never delays it
(`MPLBACKEND=Agg` makes this a safe no-op for a non-interactive script,
same convention `example_4.4.sh` uses).

### 4.7 `[6] LIBRARY WARNINGS` — what it captures, and how

Every third-party call this tool makes (`pyxtal`'s pattern computation and
`Similarity`, `pymatgen`/`spglib`'s space-group detection, the peak-list
broadening) is wrapped in `capture_library_noise()` (`core/cli.py`, shared
with `stb-phononsCreate`/`stb-xrdsearch`): both `stdout` prints **and**
`warnings.warn()` calls made *inside* the wrapped block are captured into
one collector list, instead of leaking straight to the terminal interleaved
with the report. This tool's *own* `print_dual` report lines are never
inside a wrapped block, so the report itself is never delayed or hidden.

A real, subtle gap this session found and fixed (verified against
`/home/carlos/test/Si`'s own real data, Section 5): the space-group lookup
(`core/symmetry.py::space_group_label`, called once per candidate to fill
the ranking table's own column) used to sit in its *own* `try/except`,
**outside** the wrapped block that computes the pattern/similarity for that
same candidate. Python's own warning system shows a given warning only
*once per call site* by default — depending on whether an identical
warning had already fired earlier inside some *other* wrapped block in the
same run, this specific one would sometimes vanish untracked and sometimes
leak straight to the terminal, unpredictably. Moving the
`space_group_label()` call *inside* the same `capture_library_noise` block
as the rest of that candidate's comparison closes the gap regardless of
ordering — reproduced live:

```
$ stb-xrdrank --experimental experimental.dat --no-intro
...
[6] LIBRARY WARNINGS
------------------------------------------------------------
Messages emitted by external libraries (pyxtal/scipy/numpy/pymatgen/spglib) during this run...
[pyxtal XRD ranking] DeprecationWarning: Set OLD_ERROR_HANDLING to false and catch the errors directly.
```

That one line is the *only* place this message appears anywhere in the
tool's output — never before the report banner, never interleaved with
`[3] RANKING`'s own live progress lines.

### 4.8 Running it both ways

Both CLI and `stb-suite` -> `4.6.2` produce the identical report for
identical answers — the interactive menu prompts for the same choices
(input folder, defaulting to `xrd_search`; experimental file; auto-
broadening; `--top`; save-gnuplot/view/save-report) in the same order as
the CLI flags above.

## 5. Worked example: real production data, end-to-end

This is not illustrative data — this is a real calculation the user ran
this session, at `/home/carlos/test/Si`: `stb-xrdsearch` cast 3 candidates
for `Si8`, space group `227` (`Fd-3m`), each pre-relaxed with MACE, then
**actually relaxed with real SIESTA** (GGA-PBE, DZP, `320 Ry` mesh cutoff,
`6x6x6` Monkhorst-Pack — `132` irreducible k-points, force tolerance
`0.01 eV/Ang`) — a genuinely different question from Section 4.3's clean
right-vs-wrong-space-group case: all 3 candidates here **are** the same
correct structure; the question is which SIESTA run relaxed it best.

### 5.1 First run: `nan` for every candidate

The user's real experimental file, `xrd.int`, is the exact 3-column format
Section 4.5 works through (2theta, intensity, an always-zero 3rd column).
Before this session's fix, every candidate scored `similarity nan`:

```
[3] RANKING
  -> group_227_1: similarity nan (Si, Fd-3m (No. 227), relaxed)
  -> group_227_2: similarity nan (Si, Fd-3m (No. 227), relaxed)
  -> group_227_3: similarity nan (Si, Fd-3m (No. 227), relaxed)
[6] LIBRARY WARNINGS
[pyxtal XRD ranking] RuntimeWarning: invalid value encountered in scalar divide
```

Exactly the mechanism Section 4.5 traces, on real data, in real production
use — not a contrived test case.

### 5.2 After the fix: real numbers, and a real relaxation-quality finding

```
Rank | Name        | Similarity | Formula | Space group     | Source
1    | group_227_2 | 0.9922     | Si      | Fd-3m (No. 227) | relaxed
2    | group_227_3 | 0.9539     | Si      | Fd-3m (No. 227) | relaxed
3    | group_227_1 | 0.9537     | Si      | Fd-3m (No. 227) | relaxed
Score gap (#1 vs #2)   : 0.0383
```

Investigating **why** `group_227_2` scored differently (SIESTA's own
`calc.out` files, not this workflow's report) found a real relaxation
-quality difference: `group_227_1`/`group_227_3` fully relaxed over 4 CG
steps to `a = 5.4800 Ang` with a near-zero residual stress
(`~2e-4 eV/Ang^3`); `group_227_2` stopped after just **1** CG step at
`a = 5.4608 Ang`, with a residual stress (`~0.0058 eV/Ang^3`, `~0.92 GPa`)
just under SIESTA's default stress tolerance (`MD.MaxStressTol` was left
commented out in the auto-generated `calc.fdf`, Section 3.2) — **not**
fully converged. Reconstructing the experimental pattern's own implied
lattice constant from its 6 strongest peaks (Bragg's law, Section 1.1) gives
`a_exp = 5.4561 +/- 0.0002 Ang` — closer to the *real* literature value for
Si (`5.4310 Ang`) than either PBE result, but still measurably different
from both, and notably close to typical published PBE-GGA silicon values
(`~5.463-5.470 Ang`) — suggesting `xrd.int`, despite its name, is more
likely a theoretical/simulated reference pattern than a genuine lab
measurement. **`group_227_2` "won" this ranking because its under-converged
cell happened to land closer to that reference value than the fully
-relaxed cells did — not because it is the more correct DFT result.**
A real, physically subtle finding this workflow's own numbers surface, but
cannot resolve on their own (Section 6's limitations, Section 7's
recommended next step: re-relax `group_227_2` with an explicit, tighter
`MD.MaxStressTol` before trusting this ranking).

## 6. Known, deliberate limitations

- **Not a Rietveld refinement.** This is a fast similarity pre-screen
  (Section 1.4) between simulated and experimental *patterns* — it never
  fits atomic positions/thermal factors/preferred orientation against the
  data the way a real Rietveld refinement does. Use it to prioritize which
  candidate(s) deserve real DFT, not as final proof (every report says so
  explicitly, Section 4.1's own `[5]`).
- **The similarity score can't distinguish "wrong structure" from
  "right structure, poorly relaxed."** Section 5.2's real
  `group_227_2` case is the concrete demonstration — a single number
  conflates both effects; only reading the underlying SIESTA relaxation
  quality (residual stress, number of CG steps) separates them.
- **`_select_intensity_column`'s "last non-constant column" heuristic
  (Section 4.5) is not a general-purpose column-format detector** — it
  correctly handles every real format this suite has encountered so far
  (2-column, `stb-xrd`'s 6-column stick pattern, this `.int` 3-column
  export), but a hypothetical file where intensity is genuinely *not* the
  last varying column would still be misread. No `--intensity-column`
  override exists yet (would follow this suite's own "expose on first
  genuine need" convention if one ever comes up).
- **`--mace-relax` (Section 3.4) is a heuristic pre-screen, not a DFT
  replacement** — same documented limitation as every other MACE-MP-0
  consumer in this suite (`stb-mlrelax` etc.).
- **Casting is random, not exhaustive** — `--count-per-group` samples the
  Wyckoff configuration space randomly; a low `--count-per-group` on a
  low-symmetry space group with real continuous positional freedom can
  miss the true minimum-energy arrangement entirely.
- **Wavelength/2-theta-range mismatches between the simulated and
  experimental sides are the user's responsibility** — `--wavelength`
  must match whatever source produced the experimental file; the
  comparison is only ever computed over the *overlap* of both ranges, with
  no explicit warning if that overlap is small.

## 7. Step-by-step: running this workflow on your own composition

1. **Decide the composition and candidate space group(s)** (`--species`/
   `--num-ions`/`--groups`) — if you have prior knowledge (e.g. from a
   related compound, or a short list from indexing software), list only
   those; otherwise cast a broader net and let Stage 2's ranking narrow it
   down (Section 4.3's clean 3-way example is exactly this).

2. **Decide `--count-per-group`** and whether to add `--ml-rank` (a cheap,
   transient MACE-based pre-screen *during* casting, Section 3.4) if
   you're casting many random attempts per group and want only the most
   promising ones kept (`--keep-top-per-group`).

3. **Run Stage 1** (`stb-xrdsearch`, code `4.6.1`):
   - CLI: `stb-xrdsearch --species <S...> --num-ions <N...> --groups
     <G,...> --mace-relax -p <bank>`.
   - Interactive: `stb-suite` -> `4.6.1`.
   - Writes `xrd_search/group_<G>[_<i>]/` (Section 3.1), each folder
     already containing `structure.fdf` + `calc.fdf` + pseudopotentials,
     ready to run.

4. **(Optional) run SIESTA yourself** in each `xrd_search/group_*/`
   folder — Stage 2 works equally well on the raw, pre-SIESTA candidates
   (Section 4.3) or on relaxed `*.STRUCT_OUT` results (Section 5) if you
   already have them; it auto-detects which is present per folder.

5. **Get (or prepare) your experimental pattern as the 2-column model**
   (Section 1.5/this folder's `experimental.dat`) — export it directly if
   your instrument software supports it, or extract just the `2theta` and
   intensity columns from whatever multi-column format it actually
   produces (Section 4.5's exact real-world lesson).

6. **Run Stage 2** (`stb-xrdrank`, code `4.6.2`) once every candidate is
   ready:
   - CLI: `stb-xrdrank --experimental <file> [--wavelength <src>]
     [--save-gnuplot] [--view]` (`--input-dir` needs no flag if it's still
     the default `xrd_search`, Section 4.2).
   - Interactive: `stb-suite` -> `4.6.2`.

7. **Read `[2]` before trusting `[3]`** — confirm `Peak-list detected`
   correctly matches what your file actually is (Section 4.4); a wrongly
   -skipped auto-broadening can invert the ranking (Section 4.4's own
   live before/after).

8. **Check `[6] LIBRARY WARNINGS`** — usually harmless, but a genuine
   `RuntimeWarning: invalid value encountered` there (Section 4.5/5.1) is
   a strong signal something upstream (most likely the experimental
   file's own column format) needs a second look, even if `[3]` itself
   doesn't obviously look wrong yet.

9. **Relax the top-ranked candidate(s) with real SIESTA**, then rerun
   Stage 2 on the relaxed `*.STRUCT_OUT` results to confirm (Section 5's
   own real example is exactly this second pass) — and read the *SIESTA*
   convergence details (residual stress, CG steps), not just the
   similarity score, before picking a final winner (Section 5.2/6).

## 8. Files in this folder

| File | What it is |
|---|---|
| `experimental.dat` | Real Si (diamond cubic) powder pattern, Cu-K-alpha, the canonical 2-column model (Section 1.5) -- derived from the real `/home/carlos/test/Si/xrd.int` production data with its bogus 3rd column removed |
| `experimental_peaks.dat` | The same real pattern, as a sparse 9-point peak list (Section 4.4) |
| `example_4.6.sh` | This walkthrough's runnable script (both stages) |

No `structure.fdf`/`calc.fdf` fixture is checked in here, unlike `4.1`-
`4.4` — this workflow is the one pair in the suite that starts from a bare
composition, not an existing structure (Section intro).

## 9. Running the script

```bash
bash example_4.6.sh
```

## 10. What's next

- **`4.1-strain`** / **`4.2-elastic`** / **`4.3-cohesive`** / **`4.4
  -phonons`** / **`4.5-convergence`** — once this workflow has told you
  *which* structure is right, these are the natural next questions about
  it (mechanical response, stiffness, binding energy, vibrational
  stability, numerical convergence).
- **`stb-xrd`** (Analysis, code `3.9`) — the single-structure XRD tool
  this workflow's own `core/xrd.py` pattern-computation code is shared
  with; useful on its own for simulating/inspecting one structure's
  pattern without a full candidate search.
- **`stb-mlrelax`** (Structures, code `1.7`) — Section 5's own
  recommendation for a fast MACE pre-relax pass before committing to a
  second real-SIESTA relaxation of a borderline candidate like
  `group_227_2`.
