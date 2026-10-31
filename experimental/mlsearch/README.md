# mlsearch — ML structure search (former `stb-mlsearch`, menu 5.4)

**Status: removed from STB-SUITE on 2026-10-29.** No console command, no menu entry, no page on the site or in the
manual; its test was deleted. The code is kept here unchanged (`mlsearch.py`). The other menu codes were not
renumbered then (5.4 stayed empty).

> **Renumbering (2026-10-31).** The codes on this page are the ones these tools had *before* the menu was
> renumbered; the menu now has no gaps, and ML Simulations is section 7 (old 5.1 `stb-mlmd` = 7.2, 5.2 `stb-mlphonons`
> = 7.3, 5.3 `stb-mlelastic` = 7.4, 5.7 `stb-mlneb` = 7.5, 5.10 `stb-mladsorb` = 7.6, the MLFF workflow 4.17 = 7.1,
> the EOS workflow 4.18 = 5.15, GQCA 4.15 = 5.13). The full table is in `stb-suite/CHANGELOG.md`. A tool brought back
> gets the next free number at the end of its section.

## What it does

`mlsearch.py` (547 lines) searches for a lower-energy arrangement of the atoms of a structure with a MACE potential
(the MACE-MP-0 foundation model, or a fine-tuned model with its card), **at a fixed cell**:

1. reads the `.fdf` and pre-relaxes cell and atoms with MACE (`--no-pre-relax` to skip);
2. searches, by one of two algorithms (`--algorithm`):
   - **basin hopping** (default): ASE's `BasinHopping`. Each step displaces every atom by a uniform random amount in
     ±`--dr` (0.3 Å) per Cartesian component, relaxes (FIRE, `--fmax`), and accepts or rejects with the Metropolis
     criterion at `--temperature` (1000 K);
   - **simulated annealing**: one NVT Langevin MD run (`stb.mlmd.build_dynamics`) cooled from `--temp-start` (1500 K)
     to `--temp-end` (50 K), exponential or linear; every `--snapshot-interval` steps (20) the configuration is relaxed
     and the lowest one kept;
3. relaxes the best configuration once more, cell included (`--no-final-cell-relax` to skip), and writes
   `best_structure.fdf`, `search_history.png`, the trajectory of the local minima (xsf/pdb/xyz) and a short report
   (`--save-report`, `--save-data`).

## Why it was removed

- **The use cases it advertised do not work.** Its documentation offered defect-site occupation, alloy
  configurations and vacancy arrangements. The only move is a continuous random displacement of at most 0.3 Å per
  atom; atoms of different species are never swapped, and a relaxation after 0.3 Å usually falls back into the same
  basin. Those problems need an enumeration of configurations, which the suite already has:
  `stb-defect --all-inequivalent-sites --ml-rank` (defect sites), `stb-sqs` and the GQCA workflow (`4.15`, alloys).
- **The annealing duplicates other tools.** It is a weaker version of `stb-amorphize` (melt-quench with a variable
  cell) and of `stb-mlmd`.
- **Where basin hopping does make sense — clusters in vacuum, local surface reconstructions, its original use — it
  was never validated**: no example, no comparison with SIESTA.
- **A global search leaves the model's training domain on purpose.** With a fine-tuned model that is exactly where it
  is least reliable; with the foundation model the minima belong to another level of theory (25–60 % too soft in the
  `5.2` and `5.3` examples). Every result would need a SIESTA relaxation, and the tool neither does nor suggests one.
- **Below the standard of `5.1`–`5.3`**: the report does not print the model card or its trained range; no
  deduplication of minima and no symmetry analysis of the result; "new minima found" only counts the times the running
  best energy went down; the relaxation inside ASE's `BasinHopping` ignores `--max-steps`; only a smoke test.

## What it would take to bring it back

1. **Moves for the real problems**: species swaps and atom-to-vacancy jumps, besides the displacement.
2. **Deduplication of minima** by atomic-environment fingerprint (`stb.core.atomic_env`), and the space group of each
   distinct minimum.
3. **The model card** in the report, with every visited configuration checked against the trained range.
4. **SIESTA verification**: write the N best distinct minima as folders (a `calc.fdf` template + a single-point or CG
   `config_extra.fdf`, as the workflows do) and compare energies and ordering.
5. **A guided example** (a cluster, or a small alloy cell with a known ordering) with real runs, in the style of
   `examples/5.1`–`5.3`.
6. Then: a menu entry (the next free number of ML Simulations), the `[project.scripts]` line `stb-mlsearch="stb.mlsearch:main"`, a test
   under `test/7-mlsimulations/`, the help snapshot (`python docs/dump_help.py`) and a manual section.

## Running it as it is

It still works with `stb_suite[ml]` installed (it imports `stb.core` and `stb.mlmd`):

```bash
python3 experimental/mlsearch/mlsearch.py --file structure.fdf --steps 100 --seed 1 --save-report
python3 experimental/mlsearch/mlsearch.py --file structure.fdf --algorithm simulated-annealing --steps 2000
```
