# 4.9 — Workflow: NEB / Reaction Path (`stb-neb` / `stb-nebCycle` / `stb-nebAnalysis`)

A reaction or diffusion barrier is the energy cost of the **highest point**
along the path connecting two already-relaxed structures — an adsorbate
hopping to a neighboring site, a vacancy migrating, a bond breaking and
reforming. The Nudged Elastic Band (NEB) method finds that path, and its
saddle point, without you having to guess where it is by hand. This item
has **two** stages, plus one CLI-only tool that sits between them:

- **Stage 1 (`stb-neb`, code `4.9.1`)** takes two already-relaxed endpoint
  structures you provide directly, interpolates a chain of images between
  them, optionally shapes that chain with a real climbing-image NEB on the
  MACE-MP-0 foundation model, and writes single-point SIESTA folders for
  every image.
- **`stb-nebCycle`** (CLI-only, **not** in the interactive `stb-suite`
  menu — see Section 6) takes one real-DFT NEB step at a time, meant to be
  called from your own cluster submission script's loop, alternating with
  real SIESTA runs.
- **Stage 2 (`stb-nebAnalysis`, code `4.9.2`)** reads every image's
  SIESTA energy back, estimates the barrier (both from the highest
  discrete image and from a spline fit), and re-validates the band as it
  actually exists on disk right now — not just what Stage 1 reported at
  prep time.

`example_4.9.sh` is a guided, runnable walkthrough of all of this —
**including two live demonstrations of exactly the things this README was
asked to make clear: what stb-neb requires of your input files (Section
2), and what happens when those requirements aren't met yet (the same
section, live).**

## 1. The theory, briefly

Plain linear interpolation between two endpoints is not a reaction path —
it is one **specific, usually wrong** guess at one. Two well-known
failure modes motivate everything the rest of this document explains:

- **Corner-cutting.** A straight-line interpolation between two
  Cartesian endpoints can walk an atom directly through where another
  atom already sits, if the true reaction path has to bend around it.
- **Uneven spacing / kinks.** Nothing keeps successive images evenly
  spread out along a curved path — they can bunch up in easy regions and
  leave gaps in the interesting (high-energy) one, wasting SIESTA time
  and under-resolving the actual barrier.

NEB fixes the second problem directly: it connects neighboring images
with virtual springs (spring constant `k`, `--ml-k`/`--k` below) and
relaxes every image **perpendicular** to the path (the springs' own force
along the path keeps the spacing even; nothing but the springs acts along
it). The **climbing-image** refinement (Henkelman & Jónsson, 2000) then
turns this off for exactly one image — the current highest-energy one —
and lets it climb **uphill** along the path direction too, so it
converges onto the true saddle point instead of staying at whichever
quantized image happened to sample it best. `stb-neb`/`stb-nebCycle` run
this in two stages, always: stage 1 with climbing disabled (shape the
whole band first), stage 2 (or, for the cluster loop, cycles past
`--climb-after`) with it enabled — climbing too early risks locking onto
the wrong image as "the" saddle before the band has even found its
overall shape.

`--idpp` (Image Dependent Pair Potential — Smidstrup, Pedersen, Vegge,
Jónsson, 2014) is a cheaper, purely classical/geometric alternative that
runs **before** any of the above, directly on the plain linear
interpolation: it smooths the path by minimizing pairwise interatomic
distance changes, with no notion of chemistry or energy at all. It fixes
corner-cutting far more reliably than linear interpolation alone, at
essentially no cost — reach for it whenever the endpoints differ by more
than just the reacting atom(s) themselves (Section 4 explains the
advisory that tells you when).

## 2. Input file requirements — read this before running your own study

This is the part of the workflow the tool can only partially check for
you. Two of the four requirements below are **hard, enforced,
[ERROR]-and-exit checks**; the other two are **entirely your own
responsibility** — `stb-neb` has no way to verify them and does not try
to.

### 2.1 Enforced: `--initial`/`--final` must already be relaxed geometries

`--initial`/`--final` each accept either a bare `.fdf` file (used exactly
as given — no relaxation-status concept applies to a bare file) or a
**directory**, in which case `stb-neb` looks for that directory's own
finished `siesta.XV` and prefers it over the pre-relaxation
`structure.fdf` guess sitting next to it
(`structure_io.read_relaxed_or_input`, the same helper
`stb-adsorbBsse` uses for the analogous "read the RELAXED geometry, not
the guess" need). If no `siesta.XV` exists yet, `stb-neb` still runs —
useful to preview the path — but prints a `[WARNING]`:

```
[WARNING] init_dir, final_dir: no finished 'siesta.XV' found in this directory yet --
using the pre-relaxation structure.fdf guess instead. stb-neb officially expects
already-relaxed endpoints ...; the band below is still generated, useful to preview
the path, but re-run once SIESTA has actually relaxed this endpoint.
```

`example_4.9.sh`'s `output/endpoints_unrelaxed/` section triggers this
warning live, then fabricates a `siesta.XV` and re-runs to show it go
away (`output/endpoints_relaxed/`).

**Why this matters, not just as a formality:** every image folder Stage 1
writes is a **single-point** evaluation (Section 3) — nothing ever
relaxes an endpoint further once it's inside this workflow. If
`--initial`/`--final` are not actually at their true relaxed geometry,
every energy along the "path" is measured relative to the wrong starting
and ending points, and the reported barrier is not a real reaction
barrier at all, just a number computed on an arbitrary pair of
structures.

### 2.2 Your own responsibility: `calc.fdf` must match the level of theory `--initial`/`--final` were actually relaxed at

`stb-neb` has no way to know what basis set, k-point density, mesh
cutoff, or XC functional you used to relax `--initial`/`--final` in the
first place — it only reads their **geometry**. Whatever `--calc`
template you pass here is copied into every `image_NN/`, with
`MD.TypeOfRun`/`MD.Steps` forced to a single-point evaluation
(unconditionally, Section 3) — but **everything else in that template is
used exactly as you wrote it**, including basis size, k-grid, and XC
functional.

If those settings **don't** match what actually produced `--initial`/
`--final`'s relaxed geometry, the energies `stb-nebAnalysis` compares are
not on a consistent footing with each other — you'd be comparing a
single-point energy computed with (say) a DZP basis against a geometry
that was only ever actually a minimum at SZP. This is exactly the same
"level-of-theory propagation" concern `stb-adsorbBsse`'s BSSE folders and
`stb-hubbardu`'s perturbation probes already have to get right (CLAUDE.md
calls this out as a repeated, previously-buggy pattern across this
suite) — here, it's on you to get it right by hand, since `--initial`/
`--final` can come from **any** prior calculation, not necessarily one
this suite produced (so there's no `config_extra.fdf`/report file for
`stb-neb` to read the prior level of theory back out of, unlike those
other workflows).

**A practical checklist**, before trusting any barrier this workflow
reports:

- Same `PAO.BasisType`/`PAO.BasisSize`/`PAO.EnergyShift` as the
  calculation that relaxed `--initial`/`--final`.
- Same `kgrid.MonkhorstPack` density (or at least dense enough that the
  energy is converged with respect to it — Section 4.1's own
  `stb-kgrid` walkthrough covers how to check this).
- Same `Mesh.CutOff`/`FilterCutoff`.
- Same `XC.Functional`/`XC.Authors`.
- Same `DFTD3`/`Spin` settings, if either mattered for the original
  relaxation.

`calc.fdf` in this folder documents this directly at the top: its
`MD.TypeOfRun CG` / `MD.Steps 100` represent whatever relaxation you'd
have used to converge `initial.fdf`/`final.fdf` in the first place — and
are the one part of the file `stb-neb` silently ignores and overrides,
since only the geometry-determining settings above actually matter to a
single-point evaluation.

### 2.3 Enforced: same composition, same lattice

```
[ERROR] Initial and final structures have different composition (H: 1 vs 0) --
NEB requires the exact same atoms on both endpoints (a reaction that creates or
destroys atoms is not a reaction path in this sense).
```

```
[ERROR] Initial and final structures have different lattices (largest component
difference: 0.1200 Ang, tolerance: 0.001 Ang) -- ase.mep.neb.NEB and pymatgen's
interpolation both require every image to share one exact cell (no variable-cell
NEB support in ASE). stb-neb requires --initial/--final to already share the same
lattice ... re-relax --final with --initial's own cell ... before retrying.
```

Both are hard `[ERROR]`s, exit 1, nothing written — `example_4.9.sh`
triggers both live (`output/composition_mismatch/`,
`output/lattice_mismatch/`). Composition must match exactly (a reaction
that creates or destroys atoms is a different kind of problem —
thermodynamics, not a reaction path). The lattice check has a **1e-3 Ang**
tolerance for ordinary floating-point/rounding noise between two
independently-converged calculations (below it, `--final` is silently
rebuilt bit-for-bit onto `--initial`'s own lattice matrix, no warning
needed); above it, `stb-neb` refuses outright rather than guessing which
endpoint's cell is "right" — **this used to silently override one with
the other**, which is worse than refusing, since a real physical
difference in relaxed lattice constant between two calculations usually
means something is actually wrong (not fully converged, or a genuinely
different structure), not that one should just be discarded.

This check exists because no downstream library involved supports a
per-image cell change at all: `ase.mep.neb.NEB`/`idpp_interpolate` raise
`NotImplementedError` on any per-image cell mismatch, and pymatgen's
`Structure.interpolate` raises `ValueError` on unequal lattices unless
told to interpolate them too (which `stb-neb` deliberately never does,
for the same reason Section 2.1 of `examples/4.8-adsorption/README.md`
keeps a slab's cell fixed during adsorption: a single reacting/adsorbing
atom has no physical business changing a bulk lattice constant that's
common to every image on the band).

### 2.4 Atom correspondence between the two endpoints

A subtler, related requirement: `--initial`/`--final` must list their
atoms in an order pymatgen can match up correctly — `--autosort-tol`
(default 0.5 Ang, pymatgen's own suggested value) matches atom `i` in one
endpoint to whichever atom in the other sits within that distance. If
both files were derived from the same base structure (the common case —
e.g. one is the other with a single atom moved), they already share the
exact same atom order; pass `--autosort-tol 0` to use it directly instead
of re-matching by distance (more robust for a small/densely-packed cell,
where several atoms of the same species can sit within any reasonable
distance of each other). If matching fails at the requested tolerance,
`stb-neb` automatically retries once at `--autosort-tol 0` (a
`[WARNING]`, not a dead end) before giving up.

## 3. Why every image is forced to a single-point evaluation

Every `image_NN/config_extra.fdf` Stage 1 writes contains, unconditionally
(no opt-out flag — this is the one thing about an image folder that isn't
configurable):

```
MD.TypeOfRun        CG
MD.Steps             0
```

`%include`-d ahead of your own `--calc` template (same `config_extra.fdf`
+ `%include` convention used everywhere else in this suite — SIESTA's fdf
reader is first-occurrence-wins for a duplicate directive, so the override
has to physically precede anything your own template already sets). Each
image is one independent, uncoupled sample point **on** the path — if
SIESTA were allowed to relax it, it would walk straight downhill off the
path and back toward whichever endpoint is closer, destroying the whole
point of having a chain of images in the first place. The NEB coupling
between images (the spring force, and the climbing-image correction) is
applied by `stb-neb`'s own MACE relax (modes 1/2) or by `stb-nebCycle`
(the real-DFT refinement loop) — never by SIESTA's own internal
relaxation.

`--force-spin`/`--force-vdw`/`--force-dipole` (all default **ON**, same
`config_extra.fdf` block constants `core/adsorption_sites.py` already
defines) additionally force `Spin polarized`/`DFTD3 .true.`/
`Slab.DipoleCorrection .true.` on every image, for the same reasons
`stb-adsorb`'s own reference folders do (a reacting/adsorbing system
commonly has a non-zero net spin; dispersion matters for any
physisorption-like point on the path; a one-sided slab reaction breaks
inversion symmetry along the vacuum axis). `--no-force-*` opts out of any
of the three individually if your system genuinely doesn't need it (e.g.
`--no-force-dipole` for a bulk/fully-periodic reaction path with no
vacuum gap at all).

## 4. The three modes

```
$ stb-neb -i initial.fdf -f final.fdf -c calc.fdf --mode {1,2,3}
```

| Mode | Path engine | Output | SIESTA's role |
|---|---|---|---|
| 1 | 100% MACE-MP-0 (climbing-image NEB), then ONE single-point read | `image_NN/` | one round, energy only |
| 2 (default) | MACE-MP-0 shapes the path, then a FEW real-DFT refinement cycles | `cycle_00/image_NN/` + printed loop | several rounds, real forces |
| 3 | 100% real-DFT NEB from a plain interpolated path, no MACE at all | `cycle_00/image_NN/` + printed loop | every round |

A pure MACE-MP-0 path with **no** SIESTA at all — useful for a fast,
purely-ML barrier estimate, e.g. screening many candidate diffusion
events before committing to any real DFT — is **not** a mode here at all;
that capability lives exclusively in `stb-mlneb` (ML Simulations menu,
code `5.x`), which already existed as a dedicated, pure-MACE extraction
of the same underlying interpolation/NEB mechanics `stb-neb` itself
reuses (`neb.py`'s own `check_composition_match`, `wrap_into_cell`,
`linear_interpolate_images`, `idpp_refine_images`,
`cumulative_reaction_coordinates`, `check_path_quality`,
`write_path_trajectory`, `write_ml_preview_plot` are all shared, not
duplicated, between the two tools).

**Choosing between 1/2/3:**

- **Mode 1** is the cheapest way to get a real reaction-path estimate:
  MACE-MP-0 does all the path-shaping work, SIESTA is only asked to
  confirm each image's energy once. Good for a quick look, or for a
  system MACE-MP-0's foundation training already covers reasonably well.
- **Mode 2** (the default) is the right choice for a production barrier:
  MACE-MP-0 gives the real-DFT refinement loop a physically sensible
  starting shape (fewer `stb-nebCycle` rounds needed to converge), and
  `--climb-after 0` is safe from the very first cycle, since the band
  already has a sensible overall shape by the time SIESTA sees it.
- **Mode 3** skips MACE-MP-0 entirely — appropriate when you don't trust
  (or don't have) a MACE-MP-0 potential for your chemistry, at the cost
  of needing more real-DFT refinement cycles, and a later
  `--climb-after` (the printed cluster snippet suggests 5, vs. mode 2's
  0) — climbing too early on a raw interpolated path risks locking onto
  the wrong image as the saddle before the band has found its shape.

`--ml-freeze-substrate` (default ON, modes 1/2 only) freezes any atom
whose position barely differs between `--initial`/`--final`
(`--ml-freeze-threshold`, default 0.3 Ang) during the MACE relax — fewer
degrees of freedom for the optimizer, and it stops a spectator atom from
spuriously wandering during a fast/loose ML relax. `--ml-prerelax-endpoints`
(independent of `--mode`) additionally relaxes both endpoints' own
positions with MACE-MP-0 before interpolating — a cheap safety net, not a
substitute for Section 2.1's real requirement.

Interactive path: `stb-suite`, type `4.9.1` — same questions, same
defaults, as the flags above. `example_4.9.sh` proves the two agree.

## 5. Stage 1 mechanics, live (`stb-neb`, code `4.9.1`)

```
$ stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 3
```

writes everything under a self-contained `neb_run/` subfolder of
`--output-dir` (matching `stb-adsorb`'s own `sites/`-nesting convention),
never `--output-dir` directly — `stb-nebCycle`/`stb-nebAnalysis --dir`
must point at `<output-dir>/neb_run`, not `<output-dir>` itself (both
default to exactly that, so most of the time you never need to think
about it).

Report structure (`[0]`–`[9]`, always saved to `neb_run/neb_setup.txt`):

```
[0] RUN METADATA
[1] ML PRE-RELAX ENDPOINTS      (only with --ml-prerelax-endpoints)
[2] ENDPOINT CHECKS
[3] PATH INTERPOLATION
[4] MACE-MP-0 PATH SHAPING      (modes 1/2 only)
[5] PATH QUALITY
[6] IMAGE FOLDERS
[7] CLUSTER SUBMISSION          (modes 2/3 only)
[8] SUMMARY
[9] LIBRARY WARNINGS
```

`[5]` always writes `neb_run/neb_path.xyz` — every image, one multi-frame
extended-XYZ file, viewable directly in VESTA/OVITO/ASE-GUI — **before**
any real-DFT refinement has happened (it's the pre-refinement, MACE
-shaped-or-not path). `[5]` also runs `check_path_quality`: a `[WARNING]`
if two neighboring images are nearly identical (step below 0.05 Ang —
consider fewer `-n`) or if one step is unusually large relative to the
mean (usually a sign of a bad `--autosort-tol` match, not a genuinely
non-uniform path, unless atom correspondence is already proven — Section
7). `check_endpoint_displacement` (in `[2]`) separately recommends
`--idpp` whenever a broad fraction of the structure (not just the
reacting atom itself) moved between the two endpoints.

Interactive path: `stb-suite`, type `4.9.1`.

## 6. `stb-nebCycle` mechanics, live (CLI-only — modes 2/3 only)

```
$ stb-nebCycle --dir neb_run --fmax 0.05 --climb-after 0
```

Deliberately **not** registered in the interactive `stb-suite` menu (see
"CLI-only tools" in `CLAUDE.md`) — it's meant to be called once per real
-DFT NEB step, from inside your own cluster submission script's loop
(Stage 1's own printed snippet, Section 4 above, is exactly that loop),
not from an interactive session. Each call:

1. Reads the **latest** `cycle_NN/image_*/` with finished SIESTA results
   (`calc.out` energy + `<SystemLabel>.FA` forces).
2. Takes exactly **one** climbing-image NEB step, using ASE's `FIRE`
   optimizer — its velocity/timestep state is persisted to
   `neb_cycle_state.json` between calls (`ase.optimize.FIRE`'s own native
   restart mechanism), since each call is a separate, short-lived process
   with no live optimizer to keep running in memory between real SIESTA
   jobs that may finish hours apart.
3. Either writes `cycle_{N+1}/` (the next geometry — `structure.fdf`/
   `calc.fdf` only, no `calc.out` yet, waiting for you to run SIESTA
   there) or, once the max residual force drops below `--fmax`, writes a
   `NEB_CONVERGED` sentinel file and stops — your submission-script loop
   should check for this file and break out (the printed snippet already
   does).

`--climb-after N` only enables the climbing-image correction once the
**current** cycle number is `>= N` — `0` (mode 2's suggestion) is safe
once MACE-MP-0 has already shaped the band; a later value (mode 3's
suggestion: `5`) gives a raw interpolated path a few cycles to find its
overall shape first, same reasoning `core/mace_relax.py::relax_neb`'s own
two-stage `climb=False`-then-`True` approach uses, just spread across
separate cluster-queued cycles instead of one live process.

`example_4.9.sh`'s `output/mode2_cycle/` section fabricates two rounds of
`.FA`/`calc.out` (a smooth harmonic "pull toward a fixed target" force
field — **not** random noise, which was found live to trigger a real ASE
`improvedtangent` 0/0 tangent edge case for two near-degenerate images)
and shows both outcomes: a normal step (`cycle_01/` written) and, with a
looser `--fmax`, convergence (`NEB_CONVERGED` written).

## 7. Stage 2 mechanics, live (`stb-nebAnalysis`, code `4.9.2`)

```
$ stb-nebAnalysis --dir neb_run --save-report --save-gnuplot --save-path-xyz
```

Report structure (`[0]`–`[7]`, `--save-report` persists it to
`neb_run/neb_report.txt`):

```
[0] RUN METADATA
[1] IMAGE ENERGIES
[2] CONSISTENCY CHECKS
[3] BARRIER ANALYSIS
[3b] BARRIER VS. CYCLE          (modes 2/3 only, once 2+ complete cycles exist)
[4] ENERGY PROFILE PLOT
[5] SUMMARY
[6] APPLY                       (only with --apply)
[6b] PATH EXPORT                (only with --save-path-xyz)
[7] LIBRARY WARNINGS
```

`--dir` defaults to `neb_run` (auto-detected as `.` instead if you're
already standing inside it — e.g. right after Stage 1's own printed
cluster-submission snippet's `cd "$run_root"`). For modes 2/3, it
auto-discovers the **right** cycle to analyze: not necessarily the
highest-numbered `cycle_NN/` — while not yet converged, that folder is
always `stb-nebCycle`'s freshly-written **next** geometry guess, with
zero SIESTA results yet (a real, previously-live bug reproduced against a
real 30+-cycle run: analyzing the highest-numbered folder between two
`stb-nebCycle` calls used to hard-fail with "No valid image results
found" instead of reporting the previous, real cycle's numbers). Once
`NEB_CONVERGED` exists, the highest-numbered cycle correctly **is** the
converged one.

`[2] CONSISTENCY CHECKS` re-validates the band **as it exists on disk
right now** — not just what Stage 1 reported at prep time: declared vs.
found image count, reaction-coordinate monotonicity, cross-image
composition, a recomputed path length cross-checked against Stage 1's own
`neb_setup.txt` record, and a plausibility check on the spline-fitted
transition state (warns if it sits within 5% of either endpoint). `[3]`
reports both the discrete highest-energy image's barrier and a
cubic-spline-fitted one (needs >= 4 images with distinct reaction
coordinates) — a smoother, less grid-dependent estimate, the DFT-side
analogue of `stb-neb`'s own `NEBTools.get_barrier(fit=True)` for the MACE
band. `--apply <path>` copies the highest-energy image's `structure.fdf`
(the approximate transition-state guess) forward; `--save-path-xyz`
writes the **currently**-analyzed band (whichever cycle `[0]` reports,
converged or not) as `neb_run/neb_path_current.xyz` — distinct from
Stage 1's own `neb_path.xyz`, which is the pre-refinement path and is
never updated afterwards.

Interactive path: `stb-suite`, type `4.9.2`.

## 8. Worked example: a symmetric hop, live numbers

`example_4.9.sh`'s `output/mode3/` section fabricates `calc.out` for all
5 images of a mode-3 (no MACE, fastest) band, hand-chosen (same "no real
SIESTA binary available" convention every example in this suite uses) so
that `image_00` and `image_04` — physically **equivalent** carbons, by
the lattice's own symmetry — share the exact same energy:

| Image | E (eV) | Max force (eV/Ang) |
|---|---|---|
| `image_00` | `-300.000000` | `0.010000` |
| `image_01` | `-299.700000` | `0.015000` |
| `image_02` | `-299.500000` | `0.400000` |
| `image_03` | `-299.700000` | `0.015000` |
| `image_04` | `-300.000000` | `0.010000` |

`stb-nebAnalysis`'s report, verbatim:

```
Highest-energy image (approx. TS) : image_02  (E = -299.500000 eV)
Forward barrier  (TS - initial)   : 0.500000 eV
Backward barrier (TS - final)     : 0.500000 eV
Reaction energy  (final - initial) : 0.000000 eV, endothermic (unfavorable)
Spline-fitted barrier (smoothed)  : 0.500000 eV at reaction coordinate 0.7101 Ang (cubic spline through the 5 energies)
```

Forward and backward barriers coming out **identical**, and the
spline-fitted transition state sitting at **exactly 50%** of the path
(`[2]`'s plausibility check confirms this explicitly), is exactly what a
genuinely symmetric hop between two equivalent sites should produce —
were `example_4.9.sh` to have injected an asymmetric pair of energies
instead, this cross-check is what would have caught it. (The `[3]`
wording "endothermic (unfavorable)" at a reaction energy of exactly `0`
is a labeling edge case worth knowing about, not a bug: the check is a
plain `< 0`/`>= 0` split, so an exactly thermoneutral reaction — the
`0.000000 eV` here — prints on the "unfavorable" side of that split by
convention, rather than a separate "thermoneutral" label.)

The deliberately large residual force on `image_02` (`0.400000` eV/Ang,
well above `--force-tolerance`'s default `0.05`) also fires the expected
advisory:

```
[WARNING] 1 image(s) have residual force above --force-tolerance (0.05 eV/Ang), possibly
not single-point-converged: image_02.
```

— never fatal to the barrier estimate, just a flag that this particular
image's own SCF/relaxation quality is worth a second look before trusting
the number too far.

## 9. Known, deliberate limitations

- **No SIESTA run is ever performed by this suite.** `stb-neb` and
  `stb-nebCycle` generate/read input and output files; you run SIESTA
  yourself in every `image_NN/`, alternating with `stb-nebCycle` calls
  for modes 2/3.
- **`stb-neb` cannot verify Section 2.2's level-of-theory requirement.**
  Unlike `stb-adsorbBsse`/`stb-hubbardu`, which read a PRIOR stage's own
  `config_extra.fdf` to inherit the right settings automatically, `--initial`/
  `--final` can come from any prior calculation at all — there is no
  report file for `stb-neb` to read a "what level of theory was this
  relaxed at" answer back out of. Getting `calc.fdf` right is entirely on
  you (the checklist in Section 2.2).
- **The lattice-mismatch tolerance (`1e-3` Ang) is a hard-coded
  constant**, not user-configurable — appropriate for ordinary
  floating-point noise between two independently-converged calculations,
  not for a genuinely different (even if small) relaxed lattice constant,
  which is refused outright rather than silently rebuilt.
- **A single fixed lattice for the whole band is a hard requirement, not
  a simplification this tool chose** — no variable-cell NEB exists in
  either `ase.mep.neb` or pymatgen's interpolation, so there is no
  "relax the cell too" option to offer even if it were physically
  desirable for your reaction (e.g. a reconstruction that genuinely
  changes the lattice constant along the way is out of scope entirely).
- **The reported barrier is only as good as the path.** Mode 3's plain
  interpolated-then-real-DFT-refined path, without a MACE-shaped starting
  guess, may need many more `stb-nebCycle` rounds to reach a physically
  meaningful saddle point than modes 1/2 — `[3b] BARRIER VS. CYCLE`'s
  convergence-history plot (once 2+ complete cycles exist) is the honest
  way to check whether the barrier has actually settled down, rather than
  trusting a single early cycle's number.
- **`--ml-freeze-substrate`'s threshold (0.3 Ang default) is a blunt,
  global cutoff** — every atom below it is frozen during the MACE stage,
  regardless of whether it's truly a spectator or just happens to move a
  little; for a reaction with genuine, broad substrate relaxation, lower
  the threshold or use `--no-ml-freeze-substrate`.
- **Mode 1's single-point SIESTA read has no coupling between images at
  the DFT level at all** — Section 4's caveat table calls this out
  directly (`stb-nebAnalysis`'s own `[3]` prints a matching `[NOTE]` for
  this exact case): a real chemical effect MACE-MP-0's foundation
  training doesn't capture well can still produce a wrong path shape that
  SIESTA never gets a chance to correct, since nothing relaxes the images
  further at the DFT level in this mode.
- **A pure MACE-MP-0 path with no SIESTA at all is not available here**
  (Section 4) — use `stb-mlneb` (ML Simulations menu) if that's genuinely
  what you want; forcing it through this workflow's modes isn't possible
  by design, not an oversight.

## 10. Step-by-step: running this workflow on your own structures

1. Have two **already-relaxed** endpoint structures (Section 2.1) sharing
   the exact same composition and lattice (Section 2.3), and a `calc.fdf`
   template using the **same level of theory** that relaxed them (Section
   2.2).
2. Run Stage 1: `stb-neb -i <initial> -f <final> -c <calc.fdf> --mode 2`
   (or `1`/`3` — Section 4). Read `[2]`/`[5]`'s advisories before
   proceeding — consider `--idpp` if recommended.
3. For mode 1: run SIESTA yourself in every `neb_run/image_NN/`, then
   skip to step 5.
4. For modes 2/3: run SIESTA yourself in every
   `neb_run/cycle_00/image_NN/`, then loop `stb-nebCycle --dir neb_run
   --fmax <target> --climb-after <N>` (Stage 1's own printed snippet is
   ready to paste into a real submission script), running SIESTA in each
   newly-written `cycle_{N+1}/image_NN/` in between, until
   `neb_run/NEB_CONVERGED` appears.
5. Run Stage 2: `stb-nebAnalysis --dir neb_run --save-report
   --save-gnuplot --save-path-xyz`. Read `[2]` for the re-validated
   consistency checks and `[3]` for the barrier — both the discrete and
   spline-fitted estimates. For modes 2/3, check `[3b]`'s convergence
   history once available.
6. `--apply <path>` if you want the approximate transition-state
   structure carried forward for a separate, dedicated TS-refinement
   calculation.

## 11. Files in this folder

| File | Purpose |
|---|---|
| `structure.fdf` | Bare 2-atom graphene primitive cell (same as `examples/4.8-adsorption/structure.fdf`) — not itself fed to `stb-neb`; kept as the documented starting point `initial.fdf`/`final.fdf` were built from (Section 10 of this file explains how, via `stb-supercell`/`stb-adsorb`). |
| `initial.fdf` / `final.fdf` | A 3x3 graphene supercell (18 C) with a single H atom adsorbed ontop of two nearest-neighbor carbons — a small, genuine surface-diffusion hop, not an arbitrary displacement. Treated by this walkthrough as already-relaxed (Section 2.1). |
| `calc.fdf` | Shared `calc.fdf` template (non-polarized, Gamma-only, modest basis/mesh — mechanics only, not physically converged; documents Section 2.2's level-of-theory point directly at the top). |
| `example_4.9.sh` | The guided walkthrough (**not** an automated test — see `test/4-workflow/9-neb/{prep,cycle,analysis}/test.sh` for that). Pauses between sections so you can read before moving on; safe to re-run. |
| `output/` | Created by `example_4.9.sh` when you run it (git-ignored, not checked in). See below. |

## 12. Running the script

```bash
./example_4.9.sh
```

| Case | Command(s) | What it shows |
|---|---|---|
| `output/endpoints_unrelaxed/` | `stb-neb -i init_dir -f final_dir --mode 3` | Section 2.1's `[WARNING]`, live, when no `siesta.XV` exists yet |
| `output/endpoints_relaxed/` | same, after fabricating `siesta.XV` in both | the warning gone, same command otherwise |
| `output/composition_mismatch/` | `stb-neb` against a deliberately-mismatched `final.fdf` | Section 2.3's hard `[ERROR]`, composition |
| `output/lattice_mismatch/` | `stb-neb` against a deliberately-rescaled `final.fdf` | Section 2.3's hard `[ERROR]`, lattice |
| `output/mode3/` | `stb-neb --mode 3` then `stb-nebAnalysis` | Stage 1/2 mechanics, single-point folders, **the worked example (Section 8)** |
| `output/mode1/` (needs the `ml` extra) | `stb-neb --mode 1 --ml-max-steps 30` | a real MACE-MP-0 climbing-image relax shaping the path |
| `output/mode2_cycle/` | `stb-neb --mode 2` then `stb-nebCycle` x2 | the cluster-submission snippet, one real-DFT NEB step, then convergence |
| *(no folder — a diff only)* | Stage 1 via `printf … \| stb-suite` | proof the interactive menu (`4.9.1`) agrees with the CLI |

## What's next

- **`stb-mlneb`** (ML Simulations menu) — the pure-MACE, no-SIESTA
  extraction of this same interpolation/climbing-image-NEB machinery, for
  a fast pre-DFT barrier estimate with no `--calc`/pseudopotentials
  needed at all.
- **`stb-mldiffusion`** (ML Simulations menu) — a pure-MACE
  vacancy-migration-barrier screening tool, a related but distinct use
  case (many candidate hops at once, not one specific reaction path).
- **`stb-kgrid`/`stb-mlrelax`** — useful before this workflow even starts,
  for getting `--initial`/`--final` themselves converged and relaxed in
  the first place (Section 2.1/2.2).
