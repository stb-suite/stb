# mldiffusion — ML vacancy-migration screening (former `stb-mldiffusion`, menu 5.8)

**Status: removed from STB-SUITE on 2026-10-29.** No console command, no menu entry, no page on the site or in the
manual; its test was deleted. The code is kept here unchanged (`mldiffusion.py`). The other menu codes were not
renumbered then (5.8 stayed empty).

> **Renumbering (2026-10-31).** The codes on this page are the ones these tools had *before* the menu was
> renumbered; the menu now has no gaps, and ML Simulations is section 7 (old 5.1 `stb-mlmd` = 7.2, 5.2 `stb-mlphonons`
> = 7.3, 5.3 `stb-mlelastic` = 7.4, 5.7 `stb-mlneb` = 7.5, 5.10 `stb-mladsorb` = 7.6, the MLFF workflow 4.17 = 7.1,
> the EOS workflow 4.18 = 5.15, GQCA 4.15 = 5.13). The full table is in `stb-suite/CHANGELOG.md`. A tool brought back
> gets the next free number at the end of its section.

## What it does

`mldiffusion.py` (291 lines) screens the hops of one vacancy with a MACE potential (the MACE-MP-0 foundation model, or
a fine-tuned model with its card):

1. reads a PERFECT bulk structure (`--file`) and the index of the atom to remove (`--vacancy-index`, 0-based);
2. takes every atom within `--shell-tolerance` (15 %) of that site's nearest-neighbour distance as a hop target;
3. for each target builds the two endpoints — the vacancy at the site, and the target atom moved into it — exactly as
   `examples/7.5-mlneb/make_endpoints.py` does for one hop;
4. runs a climbing-image NEB on each pair through `stb.mlneb.run_single_neb` (the engine of `stb-mlneb`), and lists
   the barriers with the lowest one marked; `migration_barriers.png`, `lowest_barrier_preview.png`,
   `lowest_barrier_path.xyz`, `--save-report`.

## Why it was removed

- **It is `stb-mlneb` in a loop.** The NEB is the same code; the only addition is building the endpoints of each hop,
  a few lines that any user script (or `make_endpoints.py`) does.
- **It repeats identical calculations.** There is no symmetry reduction: in a monoatomic crystal every first-shell hop
  is equivalent, so silicon runs 4 identical NEBs and an fcc metal 12, all with the same barrier.
- **It is wrong for compounds by construction.** The neighbour search ignores the species: in GaAs the first
  neighbours of a Ga vacancy are As, and "hopping" one into it makes an antisite, not a migration. The real hop is to
  the second geometric shell (the same sublattice), which it never considers. Its own docstring admits it.
- **No check of the cell.** It removes an atom from whatever cell it is given; in a small cell the vacancy interacts
  with its periodic images, and nothing warns.
- **Below the standard of `stb-mlneb`.** Old report (no model card in full, no trained range, no committee, warnings
  not collected), and none of what makes an ML barrier trustworthy: `--siesta-check` and `--to-mlff`. The `5.7`
  example shows why that matters: a model validated on the perfect crystal gave 0.72 eV for the silicon vacancy hop,
  SIESTA 0.25 eV on the same path, with every other indicator silent.

## What it would take to bring it back

Better as an option of `stb-mlneb` than as a separate tool, so that every hop gets the card, the committee, the SIESTA
check and the export to training:

1. `stb-mlneb --from-bulk FILE --vacancy-index I [--hop-to J | --all-hops]`, building the endpoints itself;
2. hop targets of the **same species** only, nearest shell of that sublattice;
3. **one NEB per symmetry-inequivalent hop** (spglib site symmetry of the vacancy), with the multiplicity reported;
4. a warning when the vacancy is closer than ~10 Å to its periodic image;
5. a summary table of the hops; then the README/manual/test of `stb-mlneb` (now 7.5) extended; no menu code of its own is needed.

## Running it as it is

It still works with `stb_suite[ml]` installed (it imports `stb.core`, `stb.neb` and `stb.mlneb`):

```bash
python3 experimental/mldiffusion/mldiffusion.py --file si64.fdf --vacancy-index 0 --n-images 5 --save-report
```
