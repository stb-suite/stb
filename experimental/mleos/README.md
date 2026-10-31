# mleos — ML equation of state (former `stb-mleos`, menu 5.11)

**Status: removed from STB-SUITE on 2026-10-30.** No console command, no menu entry, no page on the site or in the
manual; its test was deleted. The code is kept here unchanged (`mleos.py`). The other menu codes were not renumbered then
(5.11 stayed empty).

> **Renumbering (2026-10-31).** The codes on this page are the ones these tools had *before* the menu was
> renumbered; the menu now has no gaps, and ML Simulations is section 7 (old 5.1 `stb-mlmd` = 7.2, 5.2 `stb-mlphonons`
> = 7.3, 5.3 `stb-mlelastic` = 7.4, 5.7 `stb-mlneb` = 7.5, 5.10 `stb-mladsorb` = 7.6, the MLFF workflow 4.17 = 7.1,
> the EOS workflow 4.18 = 5.15, GQCA 4.15 = 5.13). The full table is in `stb-suite/CHANGELOG.md`. A tool brought back
> gets the next free number at the end of its section.

## What it does

`mleos.py` (355 lines) computes a static equation of state with a MACE potential (the MACE-MP-0 foundation model, or a
fine-tuned model with its card):

1. reads a bulk structure (`--file`; a structure with vacuum is refused) and, by default, relaxes it with the model
   (`--no-pre-relax` skips it);
2. scales the cell isotropically to `--n-volumes` volumes over ±`--strain-range` % (factor `(1 + s/100)^(1/3)`),
   relaxing the positions at each volume (`--no-relax` skips it);
3. fits `E(V)` with `core/eos_fit.py` (`--eos birchmurnaghan|vinet|murnaghan|sj`): `V0`, `E0`, `B0`, `B0'` and R²;
4. with a custom model, repeats the scan with the foundation model and overlays both; a PNG of the curve,
   `--save-data`, `--save-report`.

## Why it was removed

- **Its numbers are already given elsewhere.**
  - `B0` of a model: `stb-mlelastic` (`5.3`), from the full stress–strain tensor, with the model card, the trained
    range and the comparison against SIESTA.
  - The same isotropic volume scan: `stb-mlphonons --qha` (`5.2`), which adds the vibrational free energy.
  - The SIESTA equation of state: the `4.18` workflow (`stb-eosInputs` → `stb-eosAnalysis`), with uncertainties,
    stress cross-check and its own example.
- **It was never brought to the standard of the other 5.x tools**: no model card in the report, no trained-range check
  along the scan, no numbered report with warnings last, no guided example.

## Bringing it back

As an option of an existing tool rather than a tool of its own: the static E(V) fit is a few lines on top of
`stb-mlphonons --qha`'s volume scan (`core/eos_fit.py` is still in the package, used by `stb-eosAnalysis`). As a tool:
a report in the 5.1–5.3 style (model card, trained range per volume, warnings last), a test under
`test/5-mlsimulations/`, the `[project.scripts]` line, the menu entry, the help snapshot and a manual section.

## Running it as it is

```bash
python3 experimental/mleos/mleos.py --file structure.fdf --help
```

It imports `stb.core` (`eos_fit`, `mace_relax`, `model_card`, `structure_io`), so it runs with `stb_suite[ml]`
installed.
