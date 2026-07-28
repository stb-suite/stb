# 3.8 — Density Plotter (`stb-density`)

## What this tool does

`stb-density` converts a SIESTA charge-density grid (`.RHO`) into 2D
slices, 3D point clouds, or a planar-averaged 1D profile — with an
interactive matplotlib preview and/or gnuplot data+script export:

- **2D slice** (default) — a planar cut through the density, at a given
  position along a chosen axis.
- **3D point cloud** (`--3d`) — the full volumetric grid as an (X, Y, Z,
  value) point cloud, optionally filtered by `--iso-min`.
- **1D profile** (`--profile`) — the planar average along one axis,
  useful for slabs/interfaces/superlattices.
- If the `.RHO` is **spin-polarized**, the net spin (magnetization)
  density is now **detected and processed automatically**, alongside the
  total charge — no flag needed to notice it.
- `--rho2` subtracts a second `.RHO` (charge transfer / bonding analysis).
- `--cube` additionally writes a Gaussian `.cube` file for real 3D
  isosurface rendering in VESTA/VMD/Avogadro.

## Why this matters (a bit of theory)

### The charge density itself

`ρ(r)` is one of DFT's most fundamental outputs — literally the electron
probability density that determines every ground-state property, per the
Hohenberg-Kohn theorems. Slicing/averaging/clouding it is the most direct
possible way to *see* a calculation's result: bonding character (covalent
bridges of density between atoms vs. ionic, spherical-atom-like density),
where charge accumulates or depletes upon adsorption/defect formation
(`--rho2`), or how density decays into a slab's vacuum region
(`--profile`).

### Spin-polarized calculations: a real, verified pitfall

For a **collinear spin-polarized** SIESTA calculation, `.RHO` stores
**two raw components** — but critically, they are the **up-spin** and
**down-spin** densities directly, *not* "total" and "spin" as separately
meaningful quantities. Reading component 0 alone is only half the
electrons; reading component 1 alone is only the *other* half — neither
is the total charge, and neither is the net magnetization on its own.
The physically meaningful quantities are:

```
Total charge density   = up + down
Net spin (moment) density = up - down
```

**This was verified against a real calculation**, not assumed: a
collinear spin-polarized SIESTA run on an O₂ molecule (`o2.RHO` in this
folder — a textbook case, since O₂'s real ground state is a **triplet**,
2 unpaired electrons, SIESTA's own log reports `|S| = 2.0`). Naively
reading raw component 0 as "charge" integrates to **7.0 e** (the up
channel alone) instead of the correct **12.0 e** (2 oxygen atoms × 6
valence electrons each); reading raw component 1 as "spin" integrates to
**5.0 e** (the down channel alone) instead of the correct **2.0 e** net
moment — which *does* match SIESTA's own reported value exactly. An
earlier version of this tool made exactly this mistake. It's fixed now
by using sisl's own combination convention (`index='total'` for `up+down`,
`index='z'` for `up-down`) — the walkthrough below reproduces both the
wrong and the corrected numbers directly, so you can see the difference.

### Planar averaging (`--profile`)

Averaging `ρ(x,y,z)` over the two in-plane directions at each point along
one axis washes out atomic-scale in-plane wiggles while keeping the
physically meaningful variation *along* that axis — e.g. how charge
density decays from a slab's surface out into vacuum. Same technique
`stb-workfunction` (3.7) uses for the electrostatic potential.

### Colorbar convention

A **non-negative** quantity (total charge) is anchored at `cbrange
[0:max]` — white always means "no charge" on a consistent baseline. A
**signed** quantity (net spin density, or a `--rho2` difference) instead
gets a range *symmetric around zero*, so the diverging blue-white-red
palette's white midpoint actually lands on zero — otherwise autorange
would center white at `(min+max)/2`, which is only 0 by coincidence.

## Limitations

- **Non-collinear / spin-orbit calculations aren't handled** — those have
  a genuinely different 4-component spin texture (charge, mₓ, mᵧ, m_z),
  not the simple 2-component up/down split this tool (and the `total`/`z`
  combination) assumes. `stb-cube`'s own `detect_spin_configuration`
  (reading `.HSX`) is the authoritative way to tell collinear-polarized
  apart from non-collinear/spin-orbit, if you need to check first.
- **A skewed cut plane** (common for hexagonal cells) can render visibly
  distorted in both the gnuplot `pm3d map` and the matplotlib `imshow`
  preview, even though the underlying data file's coordinates are
  correct — flagged with a warning, not auto-corrected.
- **`--vmin`/`--vmax` apply only to the primary (charge, or spin if
  `--spin`) quantity** — an auto-detected spin section always uses its
  own zero-symmetric range, since a manually-fixed range for one quantity
  rarely makes sense for the other.
- **`--cube` applies to one quantity at a time** (charge, or spin if
  `--spin` was given) — not both automatically in the same run.
- **A large `--3d` matplotlib preview can be slow/cluttered** — `--iso-min`
  is the same filter used for the exported point cloud; a tighter value
  gives a clearer preview too.

## The report: console output, `--save-report`, `--save-gnuplot`, `--view`

Every run prints a numbered report to the console:

| Section | Content |
|---|---|
| `[0] RUN METADATA` | label, mode, axis, spin-polarization status, output dir, active options |
| `[1] CHARGE DENSITY` | integrated total charge (sanity check), mode-specific details, files written |
| `[2] SPIN DENSITY` | (conditional — auto-detected, or `--spin`) integrated net moment, mode-specific details, files written |
| `[3] CUBE FILE` | (conditional, `--cube`) |
| `[4] REFERENCES` | writes `references.bib` (SIESTA) |
| `[5] SUMMARY & FILES` | status and a recap of every file written |

- **`--save-report`** — persist the full report to `stb_density_report.txt`.
- **`--save-gnuplot`** — also write a `.gplot` script next to each `.dat`
  file (off by default — this tool used to write the `.gplot` script
  **unconditionally** on every run; that's no longer the case, matching
  `--save-gnuplot` elsewhere in this suite).
- **`--view`** — show an interactive matplotlib preview of each quantity
  (heatmap for a slice, line plot for a profile, 3D scatter for a point
  cloud) — the same blue-white-red / white-yellow-red palette convention
  as the gnuplot script, so the two views read as the same plot.

## When you'd reach for it

- A quick visual sanity check on a finished SCF calculation's charge
  distribution.
- Visualizing charge transfer upon adsorption/defect formation (`--rho2`).
- Checking a slab's charge-density decay into vacuum (`--profile`).
- Seeing where a magnetic moment actually sits, spatially, for a
  spin-polarized system (the auto-detected spin section).
- A `.cube` file for a proper 3D isosurface in VESTA/VMD/Avogadro.

## Two ways to run it

A — direct CLI:
```bash
stb-density --label Sn3O4
```

B — interactive `stb-suite` menu:
```
$ stb-suite
Select an option (0-6, or a tool code like 4.1.2): 3.8
```
Both paths call the exact same underlying tool and produce the exact same
output — `example_3.8.sh` proves this directly at the end.

## Files in this folder

- `Sn3O4.RHO` — a real, non-spin-polarized SIESTA calculation (copied from
  `test/3-analysis/8-density/`), a 14-atom Sn₃O₄ structure.
- `o2.RHO`/`o2.XV`/`o2.fdf` — a real, **spin-polarized** SIESTA calculation
  on an isolated O₂ molecule (copied from `test/6-utils/3-cube/`) — O₂'s
  real triplet ground state, the exact fixture used above to verify the
  up/down-vs-total/spin fix.
- `example_3.8.sh` — the guided walkthrough (**not** an automated test —
  see `test/3-analysis/8-density/test.sh` for that).
- `.gitignore` — excludes generated `.dat`/`.gplot`/`.cube` files,
  `references.bib`, the report, and `output/`.

## Running the walkthrough

```bash
cd examples/3.8-stb-density
./example_3.8.sh
```

It pauses between sections (`[Press Enter to continue]`) and always starts
by wiping its own `output/`. Self-contained cases are generated:

| Folder               | What it shows                                                              |
|-----------------------|-----------------------------------------------------------------------------|
| `basic-slice/`        | A 2D slice of Sn₃O₄'s total charge density |
| `spin-auto/`          | The O₂ spin-polarization auto-detection: correct total (12 e) and net moment (2 e) |
| `rho2-diff/`          | `--rho2`: an (admittedly trivial, self-subtracted) charge-transfer difference, symmetric colorbar |
| `full-report/`        | Default (no report/gnuplot files) vs. `--save-report --save-gnuplot`, `references.bib` |

## Try it yourself

```bash
# A quick look at a finished calculation
stb-density --label my_calc --save-gnuplot --view

# Charge transfer upon adsorption
stb-density --label adsorbed --rho2 isolated.RHO --view

# A slab's charge decay into vacuum
stb-density --label slab --profile --axis 2 --view

# Real 3D isosurfaces in VESTA/VMD
stb-density --label my_calc --cube
```

## Flag reference

| Flag                | Meaning                                                                |
|-----------------------|-------------------------------------------------------------------------|
| `-l/--label`        | SystemLabel (looks for `<label>.RHO`).                                  |
| `-o/--output-dir`   | Where all generated files (and `references.bib`) land.                  |
| `--3d`/`--profile`  | 3D point cloud / 1D planar-averaged profile (default: 2D slice).       |
| `-a/--axis`         | Axis normal to the cut plane, or the profile axis (0/1/2).             |
| `-p/--pos`          | Slice position in Å (default: cell center).                            |
| `--spin`            | Process ONLY the net spin density, skipping charge.                    |
| `--rho2`            | A second `.RHO` to subtract (charge/spin transfer).                    |
| `--iso-min`         | For `--3d`: filter out points below this `|density|`.                  |
| `--cube`            | Also write a Gaussian `.cube` file.                                     |
| `--vmin/--vmax`     | Fix the colorbar range manually (primary quantity only).               |
| `--contour`         | Overlay contour lines (2D slice only).                                 |
| `--save-report`     | Persist the full report to `stb_density_report.txt`.                   |
| `--save-gnuplot`    | Also write `.gplot` script(s).                                          |
| `--view`            | Show an interactive matplotlib preview.                                 |

## What's next

`stb-workfunction` (3.7) reads the same class of SIESTA real-space grid
file (`.VT` instead of `.RHO`) for a planar-averaged electrostatic
potential; `stb-bader` (3.6) partitions the charge density into atomic
basins instead of slicing/averaging it directly.
