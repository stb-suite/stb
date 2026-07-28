# 3.16 — `stb-effmass`: Effective Mass / Band Velocity

Welcome! `stb-effmass` computes two of the most fundamental single-state
electronic-structure quantities: the **band velocity** (the group velocity
of a wavepacket at that state) and the **effective mass tensor** (how the
band curves in energy around that k-point) for one chosen (k, band)
eigenstate — most usefully the valence band maximum (VBM) or conduction
band minimum (CBM), where these quantities govern carrier transport.

This example walks through the theory, a real crash a user hit on their
own spin-polarized calculation (now fixed), and a second, deeper physics
subtlety this session's review uncovered along the way.

## 1. Theory

### 1.1 Band velocity and effective mass, from first principles

Near a band extremum, the energy dispersion can be Taylor-expanded:

```
E_n(k) ~ E_n(k0) + hbar * v . (k - k0) + (hbar^2/2) (k-k0)^T [M^-1] (k-k0) + ...
```

- **v** (band velocity) is the *first* derivative, `v = (1/hbar) dE/dk` —
  how fast a wavepacket built from this state moves through the crystal.
- **M** (the effective mass tensor) comes from the *second* derivative,
  `[M^-1]_ij = (1/hbar^2) d^2E/dk_i dk_j` — how strongly the band curves,
  which determines how heavy or light a carrier behaves under an applied
  field (`F = M^-1 . p`, Newton's second law with the vacuum mass replaced
  by this tensor).

`stb-effmass` gets both from sisl's `EigenstateElectron.velocity()`/
`.effective_mass()`, which compute these derivatives **analytically**
from the Hamiltonian's own k-derivatives (`dH/dk`, `d²H/dk²`) — not by
finite-differencing a band structure, so there's no step-size to tune and
no numerical-noise tradeoff.

**A real, verified subtlety**: the 2nd-order curvature correction needs
energy differences between ALL bands at that k-point (a Berry-curvature-
like term), so it must always be computed on the FULL multi-band
eigenstate — verified live that calling it after `.sub(band_index)` gives
wildly wrong values (differing by orders of magnitude on a real test
system). This tool always computes on the full state and only indexes
the band of interest afterward.

### 1.2 A real bug found and fixed: spin-polarized calculations used to crash

A user ran this tool on their own real, spin-polarized SIESTA calculation
and hit:

```
TypeError: SparseOrbitalBZ._ddPk() got an unexpected keyword argument 'spin'
```

Verified directly against sisl's own source: `Hamiltonian.dPk()` (the
1st-order derivative `velocity()` uses) accepts a `spin` argument, but
`Hamiltonian.ddPk()` (the 2nd-order derivative `effective_mass()` needs)
**does not have a `spin` parameter at all** in this sisl version. This
tool used to assume only non-collinear/spin-orbit (nspin=4/8)
calculations were affected — verified live that it's actually **any**
spin-polarized Hamiltonian (nspin=2 too). Fixed: `stb-effmass` now
detects any spin-resolved Hamiltonian up front and reports velocity only,
with a clear warning instead of a raw traceback (plus a `try/except` as a
second line of defense). See section 3 below for the live demonstration.

### 1.3 A second, deeper finding: per-axis values can be misleading

sisl's own `effective_mass()` returns each Voigt component (`m*_xx`,
`m*_yy`, `m*_zz`, `m*_yz`, `m*_xz`, `m*_xy`) as an **independent,
element-wise reciprocal** of the corresponding curvature term — NOT a
proper matrix inversion of the full 3x3 curvature tensor. This means
`m*_xx`/`yy`/`zz` are only rigorously "the effective mass along x/y/z"
when the off-diagonal terms are small.

Verified live on this example's own Sn3O4 VBM — assembling the FULL
curvature tensor from all 6 Voigt values, inverting THAT matrix, and
diagonalizing it gives the properly basis-independent **principal
effective masses**:

| | m*_xx | m*_yy | m*_zz | Principal masses |
|---|-------|-------|-------|-------------------|
| Sn3O4 VBM | -0.629 | -0.757 | -0.515 | **-1.320, -0.844, -0.343** |

Substantially different (up to 2-4x) from the naive per-axis reading of
the exact same state! This is common, not a rare corner case: away from
a crystal's principal axes (Sn3O4 has low symmetry), the curvature tensor
is generically NOT diagonal in the Cartesian basis. `stb-effmass` now
reports BOTH values, with an automatic `[OK]`/`[WARNING]` flag comparing
how much the off-diagonal curvature contributes.

## 2. Worked examples: what a real oxide's states actually look like

### 2.1 A deep, core-like state at Gamma (band 1)

```
$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --k-index 0 --band 1
```

| Component | Value (m0) |
|-----------|------------|
| m*_xx | 0.444 |
| m*_yy | 0.121 |
| m*_zz | 0.085 |
| m*_yz | 0.118 |
| m*_xz | 0.427 |
| m*_xy | 0.493 |

Off-diagonal/diagonal curvature ratio: **1.144** — off-diagonal terms
actually EXCEED the diagonal ones here, flagging `[WARNING]`. The
properly diagonalized principal masses (0.052, 0.585, 0.762 m0) are all
positive and physically ordinary — but along directions that are NOT
simply x, y, or z.

### 2.2 The valence band maximum: negative (hole-like) effective mass

```
$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --fermi -3.200055
```

Every principal effective mass here is **negative** (-1.320, -0.844,
-0.343 m0) — exactly the expected sign for a valence band maximum: the
band curves DOWNWARD in energy away from the VBM, so a hole (missing
electron) there behaves as if it has positive mass while the electron
itself has negative curvature-derived mass. A physically satisfying,
textbook-correct sign, now backed by a properly diagonalized value
instead of 6 independent per-axis numbers.

### 2.3 The conduction band minimum: a caution about coarse k-meshes

```
$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band cbm --fermi -3.200055
```

Principal masses here come out **mixed-sign** (-0.820, -0.522, +0.734
m0) — not the all-positive signature a genuine local energy minimum
should have in every direction. This fixture's `.WFSX` only samples a
coarse 2x2x2-equivalent k-mesh (8 explicit k-points); the point this
mesh identifies as the "global CBM" may not be an exact true minimum
along all three directions of the real, continuous band structure. A
real, honest limitation of ANY effective-mass calculation done at a
mesh-identified extremum, not a bug in this tool — worth knowing before
trusting a single coarse-mesh CBM/VBM effective mass at face value.

## 3. The real bug, demonstrated live

```
$ stb-effmass --wfsx Ospin.bands.WFSX --hsx-file Ospin.HSX --k-index 0 --band 3
```

This fixture (`spin/Ospin.*`) is a real, spin-polarized isolated oxygen
atom (converged to its physical 2 Bohr-magneton triplet ground state).
Before the fix, this exact command would crash with the raw `TypeError`
shown in section 1.2. Now:

```
[WARNING] Effective mass (2nd-order k-derivative) is not supported by
sisl for spin-resolved Hamiltonians (nspin=2) in this environment --
reporting velocity only.
...
[3] EFFECTIVE MASS (per-axis Voigt, ...): N/A (unsupported for nspin=2, ...)
[4] EFFECTIVE MASS (principal, full tensor): N/A (unsupported for nspin=2, ...)
[5] BAND VELOCITY: (still computed normally)
```

A non-collinear/spin-orbit fixture (`nc/IsolatedO.*`, nspin=4) shows the
exact same graceful degradation — this was already handled before the
fix; the new fix simply extends the SAME treatment to nspin=2, which
used to be (incorrectly) assumed safe.

## 4. What changed this session

- **Numbered `[0]`...`[8]` report** (`RUN METADATA`, `INPUT DATA`,
  `STATE SELECTION`, `EFFECTIVE MASS (per-axis Voigt)`,
  `EFFECTIVE MASS (principal, full tensor)`, `BAND VELOCITY`,
  `OUTPUT DATA & PLOTS`, `REFERENCES`, `SUMMARY & FILES`).
- **`--save-report`**, **`--save-gnuplot`** (writes `effmass.dat`/
  `velocity.dat` + a real 2-panel bar-chart `.gplot` script, off by
  default), and **`--view`** (matplotlib bar charts, off by default) —
  this tool previously only ever wrote a single `effmass.txt`, with no
  plot at all.
- **The real nspin=2 crash, fixed** (section 1.2/3 above).
- **The principal-effective-mass addition** (section 1.3 above) — new
  information this tool never reported before.
- **A real bug fixed**: `--label` + `--hsx-file` together used to be
  rejected outright — the same overly strict validation bug already
  fixed in `stb-wfdensity`/`stb-sts`/`stb-coop`.
- **`--band vbm/cbm`'s Fermi-energy source decoupled from `--label`**:
  now the same priority-ordered hierarchy (`--fermi` > `--bands-file` >
  `--fermi-file` > an auto-detected `.out` log) those three tools have,
  via the shared `core.siesta_bands.resolve_fermi_energy_hierarchy`.

## 5. Known, deliberate limitations (unchanged this session)

- Only a single (k, band) point, not swept along a path — most
  meaningful at a true band extremum (VBM/CBM); a coarse k-mesh's
  "extremum" may not be an exact local minimum in every direction
  (section 2.3 above).
- A component reading exactly 0 along a vacuum-padded (non-periodic)
  axis is expected, not an error — and the principal-mass diagonalization
  is skipped there (the curvature tensor is singular).
- Berry curvature correction to the velocity itself is NOT included
  (sisl's own documented limitation).
- Effective mass (not velocity) is unavailable for ANY spin-resolved
  Hamiltonian (section 1.2/3).

## 6. Files in this folder

| File | What it is |
|------|------------|
| `Sn3O4.selected.WFSX`/`Sn3O4.HSX` | The real full-BZ Sn3O4 calculation (same fixture `3.14-coop/` uses) |
| `spin/Ospin.bands.WFSX`/`spin/Ospin.HSX` | A real spin-polarized isolated-O-atom calculation (the exact bug-fix fixture) |
| `nc/IsolatedO.selected.WFSX`/`nc/IsolatedO.HSX` | A real non-collinear/SOC calculation (pre-existing graceful-degradation case) |
| `example_3.16.sh` | This walkthrough's runnable script |

## 7. Running it

```bash
bash example_3.16.sh
```

Or try the tool directly:

```bash
stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm \
    --fermi -3.200055 --save-report --save-gnuplot --view
```

## 8. What's next

- `stb-coop` (`3.14-coop/`) — shares this exact Sn3O4 fixture, a
  different (bonding-character) question about the same electronic
  structure.
- `stb-fatbands` (`3.10-stb-fatbands/`) — the same spin-polarized
  `Ospin` fixture, used there to find/fix the analogous nspin=2
  category-merging bug in orbital-projected weights.
- `stb-bands`/`stb-dos` — the band structure / density of states this
  single-state quantity sits within.
