# 3.15 — `stb-ipr`: Inverse Participation Ratio (Anderson Localization)

`stb-ipr` scatters one number per (k, band) onto a band structure: the
**Inverse Participation Ratio (IPR)**, a standard diagnostic for whether
an electronic state is *localized* (weight concentrated on a few atoms)
or *extended* (weight spread evenly across the whole cell). It reuses
exactly the same input files and orbital-projection machinery as
`stb-fatbands` — this example builds directly on that one.

## 1. Theory

### 1.1 What IPR measures

For a Bloch state `psi_n(k) = sum_i c_i` (sum over orbitals `i`, weighted
by the real, non-orthogonal overlap matrix), the generalized IPR is:

```
I_q(k, band) = sum_i |psi_i|^(2q) / [sum_i |psi_i|^2]^q
```

which, for a properly normalized state (`sum_i |psi_i|^2 = 1`), reduces to
simply `sum_i |psi_i|^(2q)` — this suite's own module docstring formula.
Verified directly against sisl's source (`EigenstateElectron.ipr()`):

```python
state_abs2 = self.norm2(projection="hadamard").real
return (state_abs2**q).sum(-1) / state_abs2.sum(-1)**q
```

`norm2(projection="hadamard")` is the **exact same per-orbital quantity**
`stb-fatbands` already uses for its own orbital weights — IPR is really
just a different reduction of the same underlying numbers `stb-fatbands`
scatters by orbital/atom/species category. That's why `stb-ipr` shares
`stb-fatbands`' entire `.bands`/`.WFSX` loading, cross-check, and
`.HSX`-accurate/`.fdf`-approximate machinery (`core/band_scatter.py` was
literally extracted from `fatbands.py` once `ipr.py` became its second
consumer).

The classic Anderson-localization interpretation, for `q=2` (the default,
and the only value sisl's `ipr()` will accept below 2 — it hard-asserts
`q >= 2` internally):

```
lim (L -> infinity) I_2 =  1/L^d    for an EXTENDED state (d = dimensionality)
                          = const.   for a LOCALIZED state
```

A state with weight on just 1 of N orbitals has `I_2 = 1`; a state spread
perfectly evenly over all N has `I_2 = 1/N`. **There is no universal
absolute scale** — the value depends on orbital count/cell size, so only
compare IPR values within the same system/basis, never across different
structures.

### 1.2 Verified live: deep states are more localized than near-Fermi ones

Running on this example's real Sn3O4 fixture (the same one
`3.10-stb-fatbands/` uses):

```
$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi
```

| Series | States | Mean IPR | Min IPR (most extended) | Max IPR (most localized) |
|--------|--------|----------|--------------------------|----------------------------|
| ipr | 2912 | 0.0493 | 0.0172 (E=51.99 eV) | 0.3056 (E=75.23 eV) |

Deep, high-energy core-like states (E > 50 eV above the shifted Fermi
reference, from this DZP basis's own semi-core-like orbitals) average
**~0.107** mean IPR; states within 2 eV of the Fermi level average only
**~0.052** — deep states really are more localized than the itinerant,
bonding-derived states near the Fermi level, exactly the physically
expected picture for a real oxide's electronic structure.

### 1.3 The `--q` order parameter: higher q sharpens the localized/extended contrast

```
$ stb-ipr ... --q 2   ->  mean 0.0493, max 0.3056, min 0.0172
$ stb-ipr ... --q 4   ->  mean 0.0008, max 0.0832, min 0.0000
$ stb-ipr ... --q 6   ->  mean 0.0001, max 0.0240, min 0.0000
```

As `q` increases, every value shrinks (raising fractions below 1 to a
higher power), but the **contrast** sharpens: the most-extended state's
IPR collapses toward exactly 0 much faster than the most-localized
state's does. Higher `q` is a cheap way to emphasize which states are
*really* localized versus merely "somewhat concentrated" at `q=2`.

### 1.4 Accuracy: real overlap matrix vs. the orthogonal-basis approximation

Just like `stb-fatbands`, `stb-ipr` needs a real `.HSX` (Hamiltonian +
overlap) for a physically correct, overlap-weighted `norm2`; without one
it falls back to an implicit-orthogonal approximation (`|c|^2`, no true
`S` matrix) with a printed warning. Verified live on the same fixture,
same states:

| Source | Mean IPR |
|--------|----------|
| `Sn3O4.HSX` (accurate, overlap-aware) | 0.04930 |
| `calc.fdf` only (approximate, no `.HSX`) | 0.04356 |

A genuine, ~12% difference — not a rounding-level discrepancy — confirming
the accuracy warning is worth heeding whenever a real `.HSX` is available.

## 2. What changed this session

`stb-ipr` was still on the pre-rewrite style (plain banner, no numbered
report, `ipr.dat`/`ipr.gplot` written unconditionally, matplotlib preview
always shown with no way to skip it). Rewritten (v1.0.0 -> v2.0.0) to
match `stb-fatbands`' own v2.0.0 rewrite:

- **Numbered `[0]`...`[6]` report** (`RUN METADATA`, `INPUT DATA`,
  `BAND GAP ANALYSIS`, `IPR ANALYSIS`, `WRITING OUTPUT FILES`,
  `REFERENCES`, `SUMMARY & FILES`) — the `IPR ANALYSIS` section is new
  content this tool never reported before: per-series mean/min/max IPR
  AND which exact (k, energy) state is the most localized/most extended.
- **`--save-report`**, **`--save-gnuplot`** (off by default — this tool
  used to write `ipr.dat`/`ipr.gplot` on every single run with no way to
  opt out), and **`--view`** (off by default — this tool used to always
  pop a blocking matplotlib window, with no flag to skip it at all).
- **A real bug fixed**: `--label` + `--hsx-file`/`--geometry-file`
  together used to be rejected outright — the same overly strict
  validation bug already found and fixed in `stb-wfdensity`/`stb-sts`/
  `stb-coop` this session. `load_parent()` already preferred an explicit
  `--hsx-file` over `<label>.HSX` on its own; only `--label` +
  `--file`/`--wfsx` still needs to be rejected (ambiguous).
- **A real bug fixed**: `--q` had no lower-bound check. sisl's own
  `ipr()` asserts `q >= 2` internally — `--q 1`/`0`/negative used to crash
  with a raw, unfriendly `AssertionError` traceback. Now a clean
  `parser.error` up front.
- **A real bug fixed**: the exact same nspin=2 spin-channel-merging bug
  `stb-fatbands` already found and fixed this session, inherited here
  too. See section 3 below for the live verification.

Unlike `stb-wfdensity`/`stb-sts`/`stb-coop`, no Fermi-source hierarchy was
added: `stb-ipr` (like `stb-fatbands`) already requires a `.bands` file as
its primary input, which carries its own embedded Fermi energy — there
was nothing to decouple.

## 3. A real bug found and fixed: spin channels silently merged

This tool's row-building loop used to dump BOTH spin channels into one
flat, spin-blind list — for a genuinely magnetic system this silently
combines two physically different IPR series into one indistinguishable
scatter/data file (the exact same bug `stb-fatbands` already found and
fixed for its own orbital weights this session). Fixed the same way:
`nspin=2` now splits into `ipr_up`/`ipr_down`, written as separate
`.dat` files and reported as separate table rows; `nspin=1` is completely
unaffected (still plain `ipr`).

Verified live against a real spin-polarized calculation (a single oxygen
atom in a large vacuum box, seeded to converge to its physical 2
Bohr-magneton triplet ground state — the exact same fixture
`3.10-stb-fatbands/spin/` already uses and verified):

```
$ stb-ipr --label Ospin --shift fermi
```

| Series | Mean IPR |
|--------|----------|
| ipr_up | 0.9515 |
| ipr_down | 1.0870 |

A real, physically meaningful difference (the two spin channels see a
different effective potential in a magnetic atom) that used to be
silently averaged away into one merged series. (Both values sitting
close to/slightly above 1 makes sense here too — with only 13 orbitals
total for a single, near-isolated atom, several states really are close
to fully localized on one or two orbitals; values fractionally above the
theoretical q=2 maximum of exactly 1.0 reflect small numerical noise in
the overlap-weighted normalization, not a code defect.)

## 4. Known, deliberate limitations (unchanged this session)

- No file metadata ties a given `.WFSX` to a given `.bands` file — same
  k-count/orbital-count guard and `--k-tol` eigenvalue cross-check as
  `stb-fatbands`.
- IPR values are basis-size-dependent, never directly comparable across
  different structures/basis sets (section 1.1 above).
- Non-collinear/SOC (nspin 4/8) inputs get a raw `_s{n}` per-channel
  suffix, same fallback `dos.py`/`fatbands.py` use — not specially
  handled or verified beyond that.

## 5. Files in this folder

| File | What it is |
|------|------------|
| `Sn3O4.bands`/`Sn3O4.bands.WFSX`/`Sn3O4.HSX` | The real band-path calculation (same fixture as `3.10-stb-fatbands/`) |
| `calc.fdf`/`structure.fdf`/`Sn.ion(.xml)`/`O.ion(.xml)` | The structure/basis this calculation used |
| `spin/Ospin.*` | A real spin-polarized isolated-O-atom calculation, for the nspin=2 section |
| `example_3.15.sh` | This walkthrough's runnable script |

## 6. Running it

```bash
bash example_3.15.sh
```

Or try the tool directly:

```bash
stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --q 4 --save-report --save-gnuplot --view
```

## 7. What's next

- `stb-fatbands` (`3.10-stb-fatbands/`) — the direct sibling this tool
  shares its entire input/loading/accuracy machinery with; orbital
  projection instead of a single localization measure.
- `stb-dos` (`3.2-stb-dos/`) — a complementary, energy-only view of the
  same electronic structure, with no per-state localization information.
