# 3.14 — `stb-coop`: COOP/COHP Bonding Analysis

`stb-coop` computes energy-resolved **Crystal Orbital Overlap Population
(COOP)** and **Crystal Orbital Hamilton Population (COHP)** curves between
user-chosen atom pairs, from a SIESTA full-Brillouin-zone `.WFSX` + `.HSX`.
Where `stb-dos` tells you *how many* electronic states exist at each
energy, COOP/COHP tell you *whether those states are bonding or
antibonding* between two specific atoms — the quantum-chemical bridge
between a band structure and a picture of chemical bonds.

This example walks through the theory, and a real, serious bug found and
fixed while reviewing this tool this session — using a real SIESTA
calculation of the tin oxide Sn3O4.

## 1. Theory

### 1.1 What COOP/COHP actually measure

For a crystal orbital (Bloch state) `psi_{n,k} = sum_mu c_mu,n(k) phi_mu`,
the overlap population between two atomic orbitals `mu` (on atom I) and
`nu` (on atom J) is `c_mu,n(k) * c_nu,n(k) * S_{mu,nu}(k)` — a term that is
**positive when the orbital combination is bonding** (in-phase, builds up
density between the atoms) and **negative when it's antibonding**
(out-of-phase, depletes density between the atoms). Summing this over
every orbital pair between atoms I and J, over every band `n`, and
Gaussian-broadening around each band's eigenvalue, gives a continuous
curve:

```
COOP_IJ(E) = sum_(n,k) w_k * [sum_(mu in I, nu in J) c_mu,n(k) c_nu,n(k) S_mu,nu(k)] * delta(E - eps_n(k))
```

COHP replaces the overlap matrix `S` with the Hamiltonian matrix `H` in
the same formula — a Hamilton-weighted version of the same idea (some
literature plots `-COHP` so that positive also means bonding; this suite
keeps sisl's own raw sign, i.e. **COOP/COHP > 0 = bonding**, verified
below rather than just asserted).

This suite reuses sisl's own `EigenstateElectron.COOP()`/`.COHP()`
(explicitly marked "experimental" in sisl's own docs), which computes one
k-point at a time; `stb-coop` itself loops every k-point, weights by the
`.WFSX`'s own Brillouin-zone weights, and reduces the resulting
orbital-pair matrix down to the atom pairs you actually asked for.

### 1.2 Why a full-BZ `.WFSX`, and why a real `.HSX`

Same requirement as `stb-sts`: COOP/COHP is a Brillouin-zone-integrated,
DOS-like quantity, so it needs a genuine k-mesh `.WFSX`
(`WriteWaveFunctions T` + an explicit `%block WaveFuncKPoints`, named
`<label>.selected.WFSX` by SIESTA's own convention) -- NOT the band-path
`.WFSX` `stb-fatbands` uses. Unlike `stb-wfdensity`/`stb-sts`, there is
also no approximate fallback here: COOP needs the real overlap matrix
`Sk()` and COHP needs the real Hamiltonian `Hk()`, neither of which exists
on a bare `.fdf`-derived `Geometry` — only a real `.HSX` carries them.

This example's fixture (`Sn3O4.selected.WFSX` + `Sn3O4.HSX`) is a real,
non-polarized GGA-PBE (+DFT-D3) calculation of the tin oxide Sn3O4 (14
atoms, 182 orbitals, DZP basis), sampled at 8 explicit k-points (a
2x2x2-equivalent mesh).

### 1.3 A worked example: a real bond vs. a non-bond

This fixture's atom 0 (Sn) has a real Sn-O bond to atom 6 (O) at 2.077
Ang, and sits 3.800 Ang from atom 1 (another Sn) -- too far to be a real
bond. Watch how differently these two pairs behave.

**With a window that includes UNoccupied states** (`--erange -10 5`, well
above the real Fermi energy of -3.204 eV), the picture is misleading:

| Pair | Largest &#124;value&#124; | Integrated over `--erange` |
|------|-----------------|------------------------------|
| 0-6 (Sn-O, bonded) | -1.101e-01 at E=3.75 eV | -1.490e-01 (**antibonding**) |
| 0-1 (Sn-Sn, not bonded) | -1.508e-02 at E=5.00 eV | +7.507e-03 (bonding) |

That looks backwards -- the real bond reads "antibonding" and the
non-bond reads "bonding"! The reason: this window captures a lot of
UNoccupied antibonding character for the real Sn-O bond (states above
E_F that would be antibonding if they were ever filled), which dominates
the integral. **Restricting the window to the actual occupied states**
(`--erange -25 -3.204`, i.e. up to the real Fermi energy) gives the
physically meaningful picture instead:

| Pair | Largest &#124;value&#124; | Integrated over `--erange` |
|------|-----------------|------------------------------|
| 0-6 (Sn-O, bonded) | +2.861e-02 at E=-21.65 eV | **+1.330e-01 (bonding)** |
| 0-1 (Sn-Sn, not bonded) | -3.501e-03 at E=-5.998 eV | **-7.660e-04 (essentially zero)** |

Now it makes physical sense: the real bond shows net bonding character
over the occupied states, and the non-bonded pair is essentially flat at
zero. **The energy window you integrate over matters** -- `stb-coop`'s
own `[3] COOP/COHP CURVE` report section prints this integrated value
with an explicit reminder that it's over `--erange`, not automatically
"the occupied states".

### 1.4 A real, serious bug found and fixed: `--bond-order`

`--bond-order` is meant to be a cheap, energy-INTEGRATED complementary
sanity check on the COOP curve above: the **Mulliken bond order**,
`B_IJ = 2 * sum M_ij * S_ij` (summed over orbitals), where `M` is the
**density matrix**. This tool used to call `H.bond_order(...)` directly
on the Hamiltonian object `H` -- but sisl's `bond_order()` just uses
whatever sparse matrix `self` holds as `M` in that formula. Passing the
Hamiltonian silently substituted `H_ij` (Hamiltonian matrix elements, in
eV, tens of eV for on-site terms) for a density-matrix population
(dimensionless, typically O(0.1-2)).

Verified live on this exact fixture, before the fix:

| Pair | Old (broken) "bond order" |
|------|------------------------------|
| 0-0 (Sn on-site) | -92.38 |
| 0-6 (Sn-O, REAL bond, 2.08 Ang) | -38.32 |
| 0-1 (Sn-Sn, NOT bonded, 3.80 Ang) | -29.04 |

A genuinely bonded pair and a genuinely non-bonded pair gave numbers of
the same order of magnitude -- the smoking gun that this measured nothing
about bonding at all, just Hamiltonian-matrix-element scale.

**The fix**: read a REAL density matrix (`--dm-file`, or an
auto-detected `<label>.DM`) and call `bond_order()` on THAT instead. One
more wrinkle needed fixing too: SIESTA's `.DM` file format doesn't store
the overlap matrix at all (confirmed live: `DM.Sk()` gave an exactly-zero
overlap trace right after a bare read) -- only `.HSX`/`.TSHS` do. So the
already-loaded Hamiltonian's own overlap column is spliced into the
freshly-read density matrix before calling `bond_order()`.

After the fix, on the exact same pairs:

| Pair | Fixed bond order | Physical read |
|------|-------------------|----------------|
| 0-6 (Sn-O, bonded) | **+0.527** | a real, single-bond-scale value |
| 0-1 (Sn-Sn, not bonded) | **-0.010** | essentially zero, as expected |

...and, as section 1.3 already showed, the energy-integrated COOP curve
(a COMPLETELY different computation -- from eigenstate wavefunction
coefficients, not the density matrix) **agrees in sign** with the fixed
bond order for both pairs (+0.133 vs +0.527 for the bond, -0.0008 vs
-0.010 for the non-bond). Two independent methods agreeing is real
evidence, not just an assertion, that this suite's COOP sign convention
(positive = bonding) is correct.

## 2. What changed this session

- **Numbered `[0]`...`[7]` report** (`RUN METADATA`, `INPUT DATA`,
  `PAIR SELECTION`, `COOP/COHP CURVE`, `BOND ORDER`, `OUTPUT DATA & PLOTS`,
  `REFERENCES`, `SUMMARY & FILES`), matching `stb-wfdensity`/`stb-sts`.
- **`--save-report`** persists the report to `stb_coop_report.txt` (off by
  default).
- **`--save-gnuplot`** -- a real gap closed: this tool used to write
  `coop.dat`/`cohp.dat` unconditionally with no way to opt out, and never
  wrote a `.gplot` script at all. Now both are together behind one
  off-by-default flag, with a real multi-pair gnuplot script (one curve
  per selected pair, via `columnheader()`).
- **`--view`** replaces the old `--no-plot`: off by default, opted INTO
  (same convention flip every rewritten Analysis tool has gotten).
- **`--shift fermi`'s Fermi-energy source decoupled from `--label`**: now
  the same priority-ordered hierarchy as `stb-wfdensity`/`stb-sts`
  (`--fermi` > `--bands-file` > `--fermi-file` > an auto-detected `.out`
  log), via the shared `core.siesta_bands.resolve_fermi_energy_hierarchy`.
- **`--label` + `--hsx-file` together used to be rejected outright** --
  the same overly strict validation bug already fixed in
  `stb-wfdensity`/`stb-sts`. `load_parent()` already preferred an
  explicit `--hsx-file` over `<label>.HSX` -- the CLI validation was
  simply stricter than it needed to be.
- **`--bond-order` fixed** to use a real density matrix (`--dm-file`/
  auto-detected `<label>.DM`) instead of the Hamiltonian -- see section
  1.4 above. Errors clearly if no usable `.DM` is found, instead of
  silently returning a wrong number.
- The interactive `stb-suite` menu (item `3.14`) now asks for the label
  and `.HSX` path separately, gained a Fermi-source submenu, and a
  bond-order/`.DM` prompt.

## 3. Known, deliberate limitations (unchanged this session)

- sisl's COOP/COHP API is explicitly "experimental" and memory-heavy for
  many energy points/large systems -- start modest with `--npoints`.
- No metadata distinguishes a full-BZ-mesh `.WFSX` from a band-path one;
  `stb-coop` only prints an advisory (few k-points, all equal weight),
  never blocks.
- No `--shift vbm/cbm` (would need a second, expensive full pass over the
  whole `.WFSX`) -- use `--shift fermi` instead.
- Not tested for non-collinear states (sisl's own docstring caveat,
  passed through unchanged).

## 4. Files in this folder

| File | What it is |
|------|------------|
| `Sn3O4.HSX` | The real Hamiltonian + overlap matrix (mandatory -- no approximate fallback) |
| `Sn3O4.selected.WFSX` | The full-BZ wavefunction file (8 explicit k-points via `WaveFuncKPoints`) |
| `Sn3O4.DM` | The real density matrix, needed for `--bond-order` (see section 1.4) |
| `Sn3O4.XV` | Geometry-only file (not actually needed by `stb-coop`, kept for reference) |
| `example_3.14.sh` | This walkthrough's runnable script |

## 5. Running it

```bash
bash example_3.14.sh
```

Or try the tool directly:

```bash
stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 \
    --erange -25 -3.204 --sigma 300 --bond-order \
    --save-report --save-gnuplot --view
```

## 6. What's next

- `stb-sts` (`3.13-sts/`) -- the other Brillouin-zone-integrated,
  energy-resolved tool this session gave the same Fermi-source hierarchy
  to.
- `stb-dos`/`stb-ipr` -- related energy-resolved electronic-structure
  quantities, without the atom-pair bonding character COOP/COHP add.
- `stb-fatbands` (`3.10-stb-fatbands/`) -- orbital-projected band
  structure, the band-path analog (not BZ-integrated) of the same
  per-orbital coefficients COOP/COHP consume here.
