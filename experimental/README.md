# experimental/

Tools taken **out** of STB-SUITE. They are not installed (no `stb-*` command), not in the `stb-suite` menu and not
in the documentation site or the manual; their tests were deleted. Each folder keeps the code as it was. When they
were removed, the codes of the remaining tools were not renumbered (5.4, 5.5, 5.6, 5.8, 5.9 and 5.11 stayed empty).

> **Renumbering (2026-10-31).** The codes on this page are the ones these tools had *before* the menu was
> renumbered; the menu now has no gaps, and ML Simulations is section 7 (old 5.1 `stb-mlmd` = 7.2, 5.2 `stb-mlphonons`
> = 7.3, 5.3 `stb-mlelastic` = 7.4, 5.7 `stb-mlneb` = 7.5, 5.10 `stb-mladsorb` = 7.6, the MLFF workflow 4.17 = 7.1,
> the EOS workflow 4.18 = 5.15, GQCA 4.15 = 5.13). The full table is in `stb-suite/CHANGELOG.md`. A tool brought back
> gets the next free number at the end of its section.

| Folder | Former menu code | Removed | Why |
|---|---|---|---|
| [`mlsearch/`](mlsearch/) | 5.4 (`stb-mlsearch`) | 2026-10-29 | basin hopping / annealing at fixed cell: its moves never swap species, so the defect/alloy searches it advertised cannot work; the annealing duplicates `stb-amorphize`/`stb-mlmd`. Details in its README |
| `mlmelting/` | 5.5 (`stb-mlmelting`) | 2026-10-29 | melting point from the Lindemann index of short NVT runs: a perfect periodic cell at fixed volume superheats, so the answer is biased high by construction; the proper method (two-phase coexistence) is not there, and the rest is a loop over `stb-mlmd` |
| `mlconvergence/` | 5.6 (`stb-mlconvergence`) | 2026-10-29 | compares the MACE-MP-0 sizes with each other, never with DFT: in the `5.2`/`5.3` examples the largest model was the worst, and agreement between sizes says nothing about being right; fine-tuned models are compared against SIESTA by `stb-mlffValidate` (4.17.3) |
| [`mldiffusion/`](mldiffusion/) | 5.8 (`stb-mldiffusion`) | 2026-10-29 | vacancy-hop screening: `stb-mlneb` in a loop over the first neighbour shell, with no symmetry reduction (4 identical NEBs in Si), species-blind (an antisite instead of a hop in a compound), and without the model checks and the SIESTA check of `stb-mlneb`. Details in its README |
| `mlgcmc/` | 5.9 (`stb-mlgcmc`) | 2026-10-29 | GCMC of a single-atom adsorbate in a rigid host: narrow scope, and mu is ill-defined -- the insertion energy includes the model's absolute energy of the isolated atom, meaningless for a fine-tuned model on SIESTA's energy scale, with no warning |
| [`mleos/`](mleos/) | 5.11 (`stb-mleos`) | 2026-10-30 | static equation of state E(V) with MACE: the bulk modulus is `stb-mlelastic`'s (5.3), the volume scan `stb-mlphonons --qha`'s (5.2), and the SIESTA E(V) is the 4.18 workflow; never brought to the report standard of the other 5.x tools (no model card, no example). Details in its README |

## Running one as it is

The code still imports `stb.core` (and `stb.mlmd` / `stb.aimd_analysis` where it did), so with `stb_suite[ml]`
installed it runs as a script, e.g.:

```bash
python3 experimental/mlmelting/mlmelting.py --help
```

Bringing one back means: fixing what the table says, a test under `test/5-mlsimulations/`, the
`[project.scripts]` line in `stb-suite/pyproject.toml`, the menu entry in `stb_suite.py`, the help snapshot
(`python docs/dump_help.py`) and a manual section.
