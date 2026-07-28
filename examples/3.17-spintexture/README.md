# 3.17 — `stb-spintexture`: Spin Texture Analyzer

`stb-spintexture` scatters `<Sx>`, `<Sy>`, `<Sz>` — the quantum-mechanical
expectation values of the Pauli spin operators — onto a band structure,
one point per (k, band). It's the tool you reach for after a
non-collinear or spin-orbit-coupled (SOC) SIESTA calculation, to see how
each electronic state's spin is oriented: Rashba splitting, topological
surface states, and canted magnetic textures all show up here.

## 1. Theory

### 1.1 What a spin texture actually is

A collinear (non-polarized or simple spin-up/down) calculation only ever
needs a single number per state: "how much spin-up vs. spin-down". A
non-collinear/SOC calculation instead represents each electron as a full
2-component spinor, and its spin direction can point ANYWHERE in 3D —
described by the expectation values of the three Pauli matrices:

```
<Sx> = <psi| sigma_x |psi>,  <Sy> = <psi| sigma_y |psi>,  <Sz> = <psi| sigma_z |psi>
```

For a properly normalized spinor, these three numbers always satisfy:

```
|S| = sqrt(<Sx>^2 + <Sy>^2 + <Sz>^2)  <=  1
```

`stb-spintexture` gets these from sisl's `EigenstateElectron.
spin_moment()`, which needs the real overlap matrix `Sk()` — the same
"hamiltonian_preferred" accuracy tradeoff as `stb-fatbands`/`stb-ipr`: a
real `.HSX` gives the true overlap-weighted value; a bare `.fdf` falls
back to an implicit-orthogonal approximation with a printed warning.

### 1.2 A real, live-verified reason the accuracy warning matters

This isn't a theoretical concern. Running this example's own fixture
(a real non-collinear isolated oxygen atom) both ways:

```
$ stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX      # accurate
$ stb-spintexture --wfsx IsolatedO.selected.WFSX --geometry-file calc.fdf      # approximate
```

| Source | Sz values | max &#124;S&#124; |
|--------|-----------|------------|
| `IsolatedO.HSX` (accurate) | all exactly ±1 (to numerical noise) | **1.000001** |
| `calc.fdf` only (approximate) | ranges from -18.16 to +18.45 | **18.451336** |

The approximate fallback doesn't just give a slightly-off number here —
it gives values that are **physically impossible** (a Pauli expectation
value can never exceed 1 in magnitude for a normalized spinor). This
tool's new `[3] SPIN TEXTURE ANALYSIS` report section catches this
automatically:

```
Sz        | -0.000000 | -18.163734 | 18.451336
|S| (Pauli expectation magnitude): mean=6.216040, max=18.451336
[WARNING] |S| exceeds 1 by more than 2% for at least one state --
this should never happen for a properly normalized spinor; check the
.WFSX/.HSX correspondence.
```

**Always pass a real `--hsx-file` for spin texture work** — the
approximate fallback exists for tools like `stb-fatbands`/`stb-ipr`
where it's merely less accurate; for spin texture, without the true
overlap matrix the result can be outright unphysical.

### 1.3 No `.bands` file needed — and why the x-axis is different

Unlike `stb-fatbands`/`stb-ipr`, this tool does NOT need a companion
`.bands` file: `spin_moment()` works directly from whatever k-points the
`.WFSX` itself contains (a band path, a full-BZ mesh, or even — as in
this example — a single Gamma-only point). Without a `.bands` file's own
k-path geometry, the x-axis is simply the 0-based k-INDEX, not a
physical k-path arc length. For a real SOC band-structure spin-texture
plot (Rashba splitting, topological surface states), pass a genuine
band-path `.WFSX` and expect evenly-spaced k-index ticks, not a
`stb-bands`-style scaled k-path.

## 2. This example's fixture: what it does and doesn't show

`IsolatedO.HSX`/`IsolatedO.selected.WFSX` is a real, non-collinear SIESTA
calculation (a single oxygen atom in a large vacuum box, Gamma-point
only). Without an explicit initial spin canting (`%block DM.InitSpin`),
it converges with:

```
$ stb-spintexture --label IsolatedO
Sz        | 0.000000  | -1.000000 | 1.000001
Most spin-polarized state: |S|=1.000001 at k-index=0, E=23.908838 eV
  (Sx=0.0000, Sy=-0.0000, Sz=1.0000)
```

`Sz` is strongly non-zero (±1.0, the physically expected value for this
open-shell atom's spin), while `Sx`/`Sy` are ~0 (numerical noise,
~1e-15). This **validates the numerical pipeline is correct** (a real,
non-trivial |S| = 1 result), but does NOT demonstrate a genuinely
"textured" spin (canted, k-dependent `Sx`/`Sy`) — that needs real
spin-orbit coupling on a heavier element with a genuine k-path, which
this suite doesn't bundle a fixture for. The `[3]` report section's
"Most spin-polarized state" line is exactly the piece of information
you'd use to find where in a real k-path the most interesting spin
texture occurs.

## 3. What changed this session

- **Numbered `[0]`...`[6]` report** (`RUN METADATA`, `INPUT DATA`,
  `ENERGY REFERENCE`, `SPIN TEXTURE ANALYSIS`, `OUTPUT DATA & PLOTS`,
  `REFERENCES`, `SUMMARY & FILES`) — `SPIN TEXTURE ANALYSIS` is new
  content this tool never reported before: per-component mean/min/max,
  the `|S| <= 1` normalization check (section 1.2), and the single most
  spin-polarized state found.
- **`--save-report`**, **`--save-gnuplot`** (previously wrote
  `spintexture_S{x,y,z}.dat`/`.gplot` unconditionally on every run), and
  **`--view`** (previously the matplotlib preview was always shown,
  blocking, with no way to skip it) — all off by default now.
- **A real bug fixed**: `--label` + `--hsx-file`/`--geometry-file`
  together used to be rejected outright — the same overly strict
  validation bug already fixed in `stb-wfdensity`/`stb-sts`/`stb-coop`/
  `stb-effmass`.
- **`--shift fermi/vbm/cbm`'s Fermi-energy source decoupled from
  `--label`**: the same priority-ordered hierarchy `stb-effmass` has
  (`--fermi` > `--bands-file` > `--fermi-file` > an auto-detected `.out`
  log) — this tool previously only ever accepted an explicit `--fermi`
  value.

## 4. Known, deliberate limitations (unchanged this session)

- Needs `nspin=4` (non-collinear) or `nspin=8` (spin-orbit); a clean
  error for any other `nspin`.
- This example's own fixture doesn't demonstrate a genuinely canted
  (k-dependent) spin texture (section 2 above) — a real SOC calculation
  on a heavier element with a genuine k-path would be needed for that.
- Only "diagonal" projection (`spin_moment` per band), not the
  orbital-resolved "hadamard" projection sisl also supports.

## 5. Files in this folder

| File | What it is |
|------|------------|
| `IsolatedO.HSX`/`IsolatedO.selected.WFSX` | The real non-collinear calculation (mandatory `.WFSX`, optional `.HSX` for accuracy) |
| `calc.fdf` | Geometry-only fallback input (used to demonstrate the real accuracy gap in section 1.2) |
| `example_3.17.sh` | This walkthrough's runnable script |

## 6. Running it

```bash
bash example_3.17.sh
```

Or try the tool directly:

```bash
stb-spintexture --label IsolatedO --shift fermi --fermi -12.0 \
    --save-report --save-gnuplot --view
```

## 7. What's next

- `stb-fatbands` (`3.10-stb-fatbands/`)/`stb-ipr` (`3.15-ipr/`) — share
  this tool's exact `.HSX`-accurate/`.fdf`-approximate accuracy
  tradeoff, and its `core/band_scatter.py` scatter-plot machinery.
- `stb-effmass` (`3.16-effmass/`) — shares this tool's exact Fermi-energy
  source hierarchy for `--shift fermi/vbm/cbm`.
