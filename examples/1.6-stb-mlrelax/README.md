# 1.6 — ML Pre-Relaxation (`stb-mlrelax`)

## What this tool does

`stb-mlrelax` runs a fast, heuristic geometry pre-relaxation using
**MACE-MP-0** — a "foundation" machine-learned interatomic potential (MLIP)
trained on ~1.5 million DFT relaxation snapshots from the Materials Project
— or a custom model you fine-tuned yourself with `stb-mlffAnalysis`, before
you commit to an expensive real SIESTA relaxation. It cleans up obviously
-wrong geometry (a bad lattice guess, a freshly-added defect/passivant atom
at the wrong bond length) in seconds, so the follow-up DFT relaxation starts
much closer to equilibrium and needs far fewer, cheaper ionic steps.

**This is a heuristic, not a substitute for DFT.** Universal MLIPs are
trained mostly on Materials-Project-like inorganic crystals and are less
reliable far from that distribution. Always relax the result with SIESTA
afterward.

## Why this matters (a bit of theory)

A real SIESTA relaxation is expensive because every ionic step needs a full
DFT self-consistent-field (SCF) cycle — solving the many-electron problem
from (nearly) scratch just to take one small step downhill in energy. An
MLIP sidesteps that entirely: it's a function, trained in advance to
*approximate* the DFT potential energy surface (PES), so evaluating an
energy and its forces on a new structure is just a forward pass through a
neural network — microseconds to milliseconds, not minutes.

**MACE** (the architecture) is an *equivariant* message-passing graph
neural network: atoms are nodes, atoms within a cutoff radius are connected
by edges, and each message-passing "round" lets an atom build up a
description of its local environment from its neighbors — two rounds in
practice, since each one implicitly reaches one hop further (MACE's own
efficiency trick over earlier, more rounds-hungry equivariant GNNs).
*Equivariant* is the physically important part: if you rotate the whole
structure, the network's predicted forces rotate the exact same way — a
hard architectural constraint, not something the model has to learn
approximately, so a MACE potential can't accidentally break rotational
symmetry the way an ordinary (non-equivariant) network could.

**MACE-MP-0** is one specific *trained* MACE model — the "foundation"
part means it was fit once, broadly, across the Materials Project's public
inorganic-crystal database, rather than for one particular material. That's
exactly what makes it directly usable out of the box on almost any
structure you hand it, at the cost of being a generalist, not a specialist
(a model you fine-tune yourself on your own SIESTA data, via
`stb-mlffAnalysis`, trades that generality for accuracy on your specific
system — usable here via `--custom-model`).

### Convergence methods: FIRE, BFGS, LBFGS, and the `fmax` criterion

Once MACE can predict a force on every atom (`F = -dE/dR`, the same physics
as DFT's Hellmann-Feynman forces), an ordinary geometry optimizer walks
"downhill" in energy, one step at a time, until the structure stops moving
appreciably. `--optimizer` picks *which* algorithm does that walk — all
three converge to the same kind of answer (a local energy minimum), they
just take different-shaped paths to get there, with different costs per
step:

- **`FIRE`** (Fast Inertial Relaxation Engine, the default) — treats the
  atoms like a damped mechanical system: move along the force direction as
  inertia would carry it, accelerating while the motion keeps pointing
  downhill, and adaptively resetting/damping the "velocity" whenever it
  doesn't. Each step only needs the current forces (no extra bookkeeping,
  no matrix), so it's cheap and very robust even from a badly-guessed
  starting geometry — but because it only ever looks at the local force
  direction, it can need more steps to fully settle once it's already close
  to the minimum.
- **`BFGS`** (Broyden-Fletcher-Goldfarb-Shanno) — a *quasi-Newton* method:
  instead of only following the current force, it builds up an approximate
  curvature (Hessian) of the energy surface from the history of positions
  and forces seen so far, and uses that curvature to jump more directly
  toward the minimum. This typically needs noticeably fewer steps once the
  structure is already reasonably close to equilibrium (superlinear
  convergence), at the cost of maintaining a full curvature matrix --
  memory that grows with the *square* of the number of degrees of freedom
  (3 x atoms, +6 more if `--relax-cell`), and a per-step cost that grows
  with it too.
- **`LBFGS`** (Limited-memory BFGS) — the same curvature-based idea as
  BFGS, but keeps only a short rolling history of recent position/force
  changes instead of the full curvature matrix. Nearly all of BFGS's
  fast, few-steps-to-converge behavior, at a memory/step cost that stays
  small regardless of system size -- the practical choice once a system is
  too large for full BFGS to be cheap.

Convergence itself is judged by **`fmax`**: the largest force component on
any single atom, in eV/Ang (default target `0.05`, tunable via `--fmax`) —
the same style of criterion SIESTA's own relaxation uses, so "relaxed
enough for `stb-mlrelax`" and "relaxed enough for SIESTA" mean the same
physical thing, just possibly at different thresholds. All three optimizers
above are judged by the exact same `fmax` target -- `--optimizer` only
changes *how fast* that target is reached, never what "converged" means.

Measured live on the same 8-atom `si_defect.fdf` (see
`output/optimizer-comparison/` below):

| `--optimizer` | Steps to converge | Wall time | Final energy |
|----------------|-------------------:|----------:|--------------:|
| `FIRE` (default) | 24 | 1.2 s | -42.952733 eV |
| `BFGS`         | 11 | 0.7 s | -42.953736 eV |
| `LBFGS`        | 11 | 0.7 s | -42.953736 eV |

BFGS/LBFGS needed under half the steps here -- the expected quasi-Newton
advantage once already reasonably close to a minimum. BFGS and LBFGS land
on identical numbers for this tiny 8-atom case: with so few degrees of
freedom, "limited memory" isn't actually limiting anything yet -- the gap
between them only opens up on much larger structures, where full BFGS's
curvature matrix becomes expensive to build and store.

## When you'd reach for it

- After `stb-defect`/`stb-passivate`/`stb-slab` adds a new atom at an
  approximate bond length.
- After hand-editing or converting a structure from another source (unit
  slips, wrong lattice constants).
- Before a real SIESTA relaxation, to cut down the number of ionic steps
  SIESTA itself needs.

## Two ways to run it

**A — direct CLI**:

```bash
stb-mlrelax -f defect.fdf
```

**B — interactive `stb-suite` menu**:

```bash
stb-suite
# at the main prompt, type: 1.6
```

`example_1.6.sh` proves the CLI and the menu agree.

## What every run does (always on)

- **`[2]`/`[6]` structure validation**, before *and* after relaxing — the
  same shared checklist `stb-inputfile`/`stb-fetch` use (see
  `core/structure_checks.py`): atom proximity, lattice handedness, and
  atomic density, each printed as its own row with an explicit
  `[OK]`/`[WARNING]`/`[SKIPPED]` status, e.g.:
  ```
  Check              | Result                                            | Status
  -------------------------------------------------------------------------------
  Atom proximity     | min. distance = 2.211 Ang (>= 0.5 Ang threshold)  | [OK]
  Lattice handedness | right-handed (positive determinant)               | [OK]
  Atomic density     | 0.0500 atoms/Ang^3 (within [0.01, 0.15])          | [OK]
  ```
  Every check that ran is listed explicitly -- not just a summary line --
  so it's always clear exactly what was verified, in both the console
  output and the report file.
- **`references.bib`** — always written, no flag: SIESTA (the output is a
  `.fdf`), the MACE architecture paper, and (unless `--custom-model`) the
  MACE-MP-0 foundation-model paper. Merges with an existing
  `references.bib` in the same folder instead of overwriting it.
- **A detailed before/after comparison table** (`[5]`, `Quantity | Before |
  After | Change`) — energy (total and per-atom), max force, full cell
  parameters + volume + density (if `--relax-cell`), mean/RMS/max atomic
  displacement (and which atom moved most), the nearest-neighbor distance
  before vs. after, the **space group** before vs. after, and the max
  **stress**-tensor component in GPa. See "Two more physical diagnostics"
  below for what these last two actually mean.
- **A convergence plot** (`<stem>_relax_convergence.png`, always) — energy
  and max force vs. optimizer step, so you can see *how* the relaxation
  converged, not just the final numbers.
- **MACE simulation detail** (`[3]`) — model type/size, device, precision,
  and, best-effort, the loaded model's approximate parameter count and
  cutoff radius (see the model-size table below).

## Optional (off by default)

- `--save-report` — also persist the full numbered report to
  `stb_mlrelax_report.txt`.
- `--save-data` — also dump the raw per-step convergence trace to
  `<stem>_relax_convergence.dat`.
- `--view` — open an interactive ASE 3D viewer with **both** the structure
  before and after relaxation as two pageable frames (needs a display).

## Two more physical diagnostics: space group and stress

**Space group before vs. after** is a genuine diagnostic, not decoration:
if a structure is supposed to be a specific symmetric phase but relaxes
into a *lower* symmetry, that's often a sign the starting guess (or the
potential itself, this far from its training distribution) is misbehaving
— whereas recovering or keeping *high* symmetry through relaxation is a
good physical sanity check (see `si_wrong_a.fdf` below: a badly-guessed
cubic Si cell relaxes back into a properly cubic one, not something
lower-symmetry).

**Stress** is force's generalization to the cell itself — how hard the
structure is "pushing" outward (or pulling inward) on its own periodic
boundary, per unit area, reported here in GPa (best-effort: a custom model
trained without stress targets just reports "not available" instead of
failing the run). A converged *positions-only* relaxation at a fixed cell
will generally still show nonzero stress — expected, since only
`--relax-cell` explicitly relaxes the cell toward near-zero stress; it's
extra information here, not a second convergence criterion.

## Model size: small vs. medium vs. large

MACE-MP-0 ships in 3 sizes — more parameters/message-passing capacity costs
more compute but (usually, across many structures — not guaranteed for any
one structure) predicts closer to the underlying DFT reference. Numbers
below were measured live in this session on the same 8-atom `si_defect.fdf`
(see `example_1.6.sh`'s `output/model-comparison/`, which reproduces this
exact table from the tool's own `[3]`/`[4]` report sections, not a separate
claim):

| `--model` | Parameters (approx.) | Cutoff radius | Wall time (8 atoms, CPU) | Final energy |
|-----------|----------------------:|:---:|---:|---:|
| `small` (default) | ~3.85 M | 6.00 Ang | 1.2 s | -42.9527 eV |
| `medium`           | ~4.69 M | 6.00 Ang | 2.3 s | -42.7278 eV |
| `large`            | ~15.85 M | 4.50 Ang | 5.0 s | -42.8940 eV |

Use `small` (the default) for routine pre-relaxation; reach for
`medium`/`large` when a small-model result looks suspicious (e.g. the
>15% lattice-change warning fired) or as an extra sanity check before
committing to expensive real DFT. `--custom-model` bypasses this table
entirely with your own fine-tuned model.

## Files in this folder

- `si_wrong_a.fdf` — bulk Si, correct chemistry, but a deliberately wrong
  lattice constant (5.00 Ang instead of ~5.43-5.46 Ang) — for `--relax-cell`.
- `si_defect.fdf` — bulk Si, correct cell, one atom nudged ~0.16 Ang off its
  ideal site — for the default positions-only mode.
- `si_too_close.fdf` — bulk Si with a mis-typed atom sitting ~0.09 Ang from
  its neighbor — deliberately malformed, to demonstrate structure
  validation actually catching a real problem.
- `example_1.6.sh` — the guided walkthrough (**not** an automated test —
  see `test/2-structures/12-mlrelax/test.sh` for that). Pauses between
  sections so you can read before moving on; safe to re-run.
- `output/` — created by `example_1.6.sh` when you run it (git-ignored, not
  checked in). See below.

## Running the walkthrough

```bash
./example_1.6.sh
```

| Folder                      | Command (conceptually)                                        | What it shows                                    |
|------------------------------|----------------------------------------------------------------|---------------------------------------------------|
| `output/positions-only/`     | `stb-mlrelax -f si_defect.fdf`                                  | default mode; cell is bit-identical before/after  |
| `output/relax-cell/`         | `stb-mlrelax -f si_wrong_a.fdf --relax-cell`                    | a=5.00 Ang recovers to ~5.46 Ang                  |
| `output/bad-structure/`      | `stb-mlrelax -f si_too_close.fdf`                               | pre-relax `[WARNING]`, clean post-relax check     |
| `output/model-comparison/`   | `stb-mlrelax -f si_defect.fdf --model {small,medium,large}`     | the speed/accuracy table above, measured live     |
| `output/optimizer-comparison/` | `stb-mlrelax -f si_defect.fdf --optimizer {FIRE,BFGS,LBFGS}`  | the FIRE/BFGS/LBFGS table above, measured live    |
| `output/full-report/`        | `stb-mlrelax -f si_defect.fdf --save-report --save-data`        | full report, `references.bib`, convergence files  |

### `bad-structure/` — the structure check catching a real problem

`si_too_close.fdf` has one atom ~0.09 Ang from its neighbor (well under the
0.5 Ang threshold). `[2] STRUCTURE VALIDATION (pre-relaxation)` flags it
with a `[WARNING]` before anything is relaxed.

*Why the starting energy is enormous*: at ~0.09 Ang apart, the two atoms'
electron clouds are forced to overlap almost completely. Physically this is
**Pauli repulsion** — electrons can't occupy the same state, so pushing two
atoms' electron clouds together costs a steeply rising energy penalty (this
is also, fundamentally, why atoms have a "size" at all). MACE learned this
repulsive wall from its DFT training data just like everything else, so the
predicted force here is huge and points straight apart — resolving the
overlap within a handful of optimizer steps.

After relaxation, `[6] STRUCTURE VALIDATION (post-relaxation)` reports "No
malformation issues detected." — so you can see the problem *and* confirm
it's actually resolved, not just relaxed into a different problem.

## Try it yourself

```bash
stb-mlrelax -f your_structure.fdf --view
stb-mlrelax -f your_slab.fdf --relax-cell   # vacuum axis stays exactly fixed
stb-mlrelax -f your_structure.fdf --custom-model my_finetuned.model
stb-mlrelax -f your_structure.fdf --model large   # slower, usually more accurate
stb-mlrelax -f your_structure.fdf --optimizer BFGS   # fewer steps once already close
```

## Flag reference

```
stb-mlrelax -f <structure.fdf> [--relax-cell] [--vacuum-gap <Ang>]
            [--model small|medium|large | --custom-model <path>]
            [--fmax <eV/Ang>] [--max-steps <N>] [--optimizer FIRE|BFGS|LBFGS]
            [--device cpu|cuda] [-o/--output <file>]
            [--save-report] [--save-data] [--view]
```

- `-f/--file` — required, the input structure (`.fdf`).
- `--relax-cell` — also relax the lattice (auto-adapts to periodicity: a
  vacuum-padded axis, e.g. a slab, is held exactly fixed; a fully isolated
  structure, e.g. a molecule, has nothing to relax and the flag is a no-op).
- `--model`/`--custom-model` — MACE-MP-0 foundation size, or your own
  fine-tuned model file.
- `--fmax`/`--max-steps`/`--optimizer`/`--device` — relaxation controls.
- `-o/--output` — output `.fdf` name (default `relaxed.fdf`).
- `--save-report`/`--save-data`/`--view` — all optional, off by default.

`references.bib` is always written — there's no flag for it. Structure
validation (before and after relaxation) always runs too — there's no flag
for it either.

Run `stb-mlrelax --help` for the full list of options.

## What's next

Take the relaxed `.fdf` and run a real SIESTA relaxation on it (e.g. via
`stb-inputfile`, example `1.1`) — it should converge in far fewer ionic
steps than starting from the original, un-relaxed guess. If you have your
own SIESTA training data, `stb-mlffAnalysis` (Workflow menu) can fine-tune a
custom MACE model for `--custom-model` that's more accurate for your
specific material than the generic foundation checkpoint.
