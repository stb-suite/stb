# 3.13 — `stb-sts`: Simulated STS Spectroscopy

`stb-sts` simulates a scanning tunneling spectroscopy (STS) curve — dI/dV
as a function of energy — at **one fixed point in real space**, from a
SIESTA full-Brillouin-zone `.WFSX` file. It answers the question "if I
parked a tip at this exact (x, y, z) and swept the bias voltage, what
would the tunneling spectrum look like?" — the energy-resolved twin of
`stb-stm`'s spatial STM images.

This example walks through the theory, the real-world file conventions
you need to get right (a full-BZ `.WFSX`, not a band-path one), and two
real bugs found and fixed while reviewing this tool this session — using
a real SIESTA calculation of monolayer graphene.

## 1. Theory: the Tersoff-Hamann proxy, extended to an energy axis

### 1.1 What a real STS measures

In a scanning tunneling microscope, the tunneling current at bias voltage
`V` and tip position **r** is, in the simplest (Tersoff-Hamann, s-wave
tip) approximation, proportional to the **integrated local density of
states**:

```
I(V, r) ~  integral_0^{eV}  LDOS(r, E) dE
```

so the differential conductance dI/dV is directly proportional to the
LDOS itself at that one bias energy:

```
dI/dV(r, E) ~ LDOS(r, E) = sum_{n,k} |psi_{n,k}(r)|^2 * delta(E - eps_{n,k})
```

`stb-stm` already uses this exact same LDOS proxy for a **spatial** STM
image: fix the energy window, scan `r` over a grid. `stb-sts` does the
opposite reduction: fix `r` at one point, scan `E`. Conceptually it is
"PDOS projected onto a single real-space point instead of onto atomic
orbitals" — the delta function is Gaussian-broadened (same
`--sigma`/`--fwhm` convention as `stb-dos`/`stb-coop`) and every
`(k, band, spin)` eigenstate is weighted by its own Brillouin-zone weight,
exactly the k-mesh sum a real DOS/LDOS needs.

### 1.2 Where this tool sits among the suite's other real-space/energy tools

This suite has four tools that each keep a different pair of "real space"
and "energy" fixed vs. integrated. It's worth seeing them side by side —
each is genuinely useful for a different question:

| Tool             | Real space (r)        | Energy (E)             | Answers                                  |
|------------------|------------------------|-------------------------|-------------------------------------------|
| `stb-wfdensity`  | full 3D grid           | ONE state (k, band)     | "What does this one Bloch state look like?" |
| `stb-dos`/`stb-coop` | integrated over all r | full curve, BZ-summed | "How many states (or bonds) exist at each E?" |
| `stb-stm`        | full 2D/3D grid        | ONE window (near E_F)  | "What does the surface look like to a tip at fixed bias?" |
| `stb-sts` (here) | ONE fixed point        | full curve, BZ-summed  | "What spectrum would a tip see, parked right here?" |

`stb-sts` and `stb-stm` are the two members of this family that actually
use the Tersoff-Hamann tunneling proxy (they simulate a tip); `stb-dos`/
`stb-coop` are pure electronic-structure quantities with no tip model at
all. `stb-wfdensity` is the only one of the four that looks at a single
eigenstate rather than a BZ-summed quantity.

### 1.3 Why a full-BZ `.WFSX`, not a band-path one

`stb-fatbands` reads a `.WFSX` written along a **band path**
(`WFS.Write.For.Bands T` + `%block BandLines`) — perfect for a band
structure, useless for a DOS-like quantity, since it isn't a real,
properly-weighted Brillouin-zone sample. `stb-sts` needs the OTHER
mechanism: `WriteWaveFunctions T` + an explicit `%block WaveFuncKPoints`
listing a genuine k-mesh, which SIESTA names `<label>.selected.WFSX` by
its own convention. This example's fixture (`calc.fdf`) uses exactly
that — a modest 9-point mesh (spanning one quadrant of the graphene
Brillouin zone, `kx, ky in {0, 0.25, 0.5}` in reciprocal-lattice units) is
enough for a didactic, fast-to-run demo; a real study would want a much
denser mesh for a converged spectrum, the same way a real DOS calculation
needs a converged k-grid.

Confirmed live while building this fixture: `WriteWaveFunctions T` alone,
with no `WaveFuncKPoints` block, silently writes nothing at all — SIESTA
needs the explicit k-point list to know which states to save.

### 1.4 The real per-orbital basis requirement

Just like `stb-wfdensity`, evaluating `psi(r)` at an arbitrary point needs
the orbitals' actual numerical radial shape — only available from a real
`.fdf` (with its `.ion`/`.ion.xml` files alongside it), never from a bare
`.XV`/`.HSX`. This example's fixture confirms the same asymmetry already
verified for `stb-wfdensity`: `calc.fdf` gives 26 real orbitals (13 per
carbon atom in this DZP basis: 2 radial functions each for 2s/2p, plus a
polarization d-shell), while `Graphene.XV` alone carries none.

### 1.5 The confined-basis cutoff radius, verified live

SIESTA's numerical atomic orbitals (PAOs) are strictly zero beyond a
finite cutoff radius — unlike a plane-wave or Gaussian basis, this is not
a numerically-small tail, it is *exactly* zero. For this fixture's carbon
DZP basis:

```python
>>> max(a.maxR() for a in geometry.atoms)
2.576278646253102   # Ang
```

Sweeping `--height` above the graphene sheet (topmost/only atomic plane
at Z = 10.000 Ang, half the 20 Ang cell) confirms this exactly:

| `--height` (Ang) | Absolute Z (Ang) | Peak dI/dV     |
|-------------------|------------------|----------------|
| 1.5                | 11.5             | 3.25e-02 (real signal) |
| 3.0                | 13.0             | 0.0 (beyond the 2.576 Ang cutoff) |
| 5.0 - 15.0         | 15.0 - 25.0      | 0.0 |

`--height 3.0` already lands 0.42 Ang beyond the cutoff radius and gives
an exactly-flat, zero curve with no error — `stb-sts` warns about this
(`[WARNING] |psi(r_tip)|^2 is exactly zero for every state ...`), but it's
easy to miss if you're not watching for it. Always keep `--height` a
couple of Angstrom at most for a typical DZP basis.

### 1.6 A real, verified physics finding: the vacuum axis is still periodic

This one is worth knowing before trusting a large `--height`. SIESTA's
"vacuum" axis in a slab calculation is a real periodic lattice vector —
just padded with empty space, not actually non-periodic. `stb-sts`
evaluates `psi(r_tip)` by translating the geometry so `r_tip` lands at the
origin of a periodic grid; if `r_tip` sits far enough along the vacuum
axis, it silently wraps around to the **next periodic image of the
surface**, not to "farther into empty space."

Verified live on this exact fixture (cell length along Z = 20 Ang, atoms
at Z = 10 Ang, so the periodic image of the atomic plane recurs at
Z = 10 + 20 = 30 Ang, i.e. `--height` = 20):

| `--height` (Ang) | Absolute Z (Ang) | Peak dI/dV | What's really being sampled |
|-------------------|------------------|------------|------------------------------|
| 8.0 - 15.0         | 18.0 - 25.0      | 0.0        | genuine vacuum (beyond the PAO cutoff either way) |
| 19.0               | 29.0             | 2.46e-01   | wrapping back toward the periodic image at Z=30 |
| 20.0               | 30.0             | 3.90e-01   | **exactly** the periodic image of the atomic plane |

Both the near-tip region (Z = 11.5-13) and the near-periodic-image region
(Z = 29-30) give strong, non-zero signal for the same underlying reason
(both are within ~2.6 Ang of a real carbon atom — one is the original,
the other its periodic copy). A large `--height` does NOT probe some
Platonic "far from any surface" limit; it eventually samples the *other*
side of the same periodic slab. In practice this only matters if your
vacuum padding is small relative to `--height` — this suite's own
`stb-slab`/`stb-2Dstacking` typically pad with 15-20+ Ang of vacuum
specifically so this doesn't bite a normal STM/STS height (a few Ang), but
it's worth keeping in mind if you ever go looking for the true vacuum
plateau by cranking `--height` up.

## 2. What changed this session — the same fixes as `stb-wfdensity`

`stb-sts` had accumulated the same class of gaps `stb-wfdensity` did, and
was rewritten (v1.0.0 -> v2.0.0) with the identical fixes applied:

- **Numbered `[0]`...`[6]` report** (`RUN METADATA`, `INPUT DATA`,
  `TIP POSITION`, `STS CURVE`, `OUTPUT DATA & PLOTS`, `REFERENCES`,
  `SUMMARY & FILES`), matching `stb-bands`/`stb-dos`/`stb-wfdensity`.
- **`--save-report`** persists the full report to `stb_sts_report.txt`
  (off by default).
- **`--save-gnuplot`** — a real bug fix: this tool used to write `sts.dat`
  unconditionally with no way to opt out, and never wrote a `.gplot`
  script at all, despite being "gnuplot output" in spirit. Now both the
  `.dat` write and a real, working `.gplot` script are together behind
  one off-by-default flag.
- **`--view`** replaces the old `--no-plot`: the matplotlib preview is now
  off by default and opted INTO, instead of on by default and opted out
  of (the same convention flip every rewritten Analysis tool this session
  has gotten).
- **Fermi-energy source decoupled from `--label`**: `--shift fermi` used
  to accept only an explicit `--fermi` value. It now has the same
  priority-ordered hierarchy as `stb-wfdensity`'s `--band vbm/cbm`:
  `--fermi` (explicit) > `--bands-file` > `--fermi-file` > an
  auto-detected `.out` log in the current directory (via
  `core.siesta_log.find_out_file`, NOT assumed to be named `<label>.out`
  -- many real SIESTA jobs redirect stdout to a generic name instead).
  This priority-ordered resolution was extracted into
  `core.siesta_bands.resolve_fermi_energy_hierarchy` once `stb-sts`
  became a second consumer of the exact logic `stb-wfdensity` already had
  — both tools now share one implementation instead of two.
- **`--label` + `--geometry-file` together used to be rejected outright**
  — the identical overly-strict validation bug found and fixed in
  `stb-wfdensity`. This example's own fixture demonstrates exactly why it
  matters: `SystemLabel` is `Graphene`, but the real input file is
  `calc.fdf` — there is no `Graphene.fdf` anywhere.
- **A silently-dropped warning, fixed**: the surface-normal axis-alignment
  check (shared with `stb-stm`, for a sheared cell where "height above the
  surface" along one axis isn't simply Cartesian) was being computed but
  its result was never printed. Now it prints, exactly like `stb-stm`
  already does.
- **The same "naive topmost atom" bug `stb-stm` already had fixed,
  inherited here too**: the `--xy`/`--height` tip-height calculation used
  a plain `xyz[:, axis].max()`, which silently picks the wrong bounding
  atom for a structure whose atoms straddle the periodic cell boundary
  (verified previously on a real CrS monolayer with atoms at fractional
  Z = 0, 0, 0.066, 0.934 — see the `3.11-stm` example for that live
  verification). `stb-sts` now reuses the same
  `core.kspace.find_surface_reference` fix `stb-stm` uses. This example's
  own graphene fixture happens to have its atoms exactly centered
  (fractional Z = 0.5 for both), so the naive and gap-aware methods agree
  here — the fix is inherited defense against structures that aren't this
  well-behaved, not something this particular fixture can show diverging.
- The interactive `stb-suite` menu (item `3.13`) now asks for the label
  and the `.fdf` path SEPARATELY, and gained an energy-shift submenu
  offering the same Fermi-source options — see the "Two ways to run it"
  section below.

## 3. Known, deliberate limitations (unchanged this session)

- **Not a full STS simulation** — the tip is a structureless point probe
  (same philosophy as `stb-stm`'s constant-height/current modes), not a
  real tip-orbital/decay model.
- **No `--shift vbm/cbm`** (unlike `stb-fatbands`/`stb-wfdensity`) — would
  need a second, expensive full pass over the whole `.WFSX` just to find
  the extremum; use `--shift fermi` directly.
- **Only the first spinor component is evaluated** — non-collinear/SOC
  wavefunctions aren't fully represented (this fixture is a plain
  non-polarized calculation, so this doesn't matter here).
- **Cost scales as O(n_k * n_bands)** point-wavefunction evaluations —
  234 contributions for this tiny 9 k-point / 26-orbital fixture already
  takes a couple of seconds; a real, converged k-mesh will take
  considerably longer.
- **No WFSX-mesh-quality metadata** — `stb-sts` can only print an
  advisory if the k-point count/weights look suspiciously like a
  band-path file (few points, all equal weight), never block outright.

## 4. Files in this folder

| File | What it is |
|------|------------|
| `calc.fdf` | The real SIESTA input (`SystemLabel Graphene` -- note the mismatched filename, see section 2 above) |
| `structure.fdf` | `%include`d geometry: 2-atom primitive hexagonal graphene cell, 20 Ang vacuum along Z |
| `C.ion`/`C.ion.xml` | The DZP carbon basis set (needed alongside `calc.fdf` for real orbital shapes) |
| `Graphene.selected.WFSX` | The full-BZ wavefunction file (9 explicit k-points via `WaveFuncKPoints`) |
| `Graphene.XV` | Geometry-only file (kept to demonstrate why it's NOT sufficient on its own -- no orbitals) |
| `example_3.13.sh` | This walkthrough's runnable script |

## 5. Running it

```bash
bash example_3.13.sh
```

Or try the tool directly:

```bash
stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 \
    --shift fermi --fermi-file calc.out \
    --save-report --save-gnuplot --view
```

## 6. What's next

- `stb-stm` (`3.11-stm/`) — the spatial-map twin of this tool: fix the
  energy window, scan real space instead of fixing real space and
  scanning energy.
- `stb-wfdensity` (`3.12-wfdensity/`) — a single eigenstate's full 3D
  `|psi|^2`, the tool this one shares its Fermi-source hierarchy and
  `--label`/`--geometry-file` fix with.
- `stb-dos` (`3.2-stb-dos/`) — the same BZ-summed, energy-resolved idea,
  but integrated over all of real space instead of evaluated at one
  point.
