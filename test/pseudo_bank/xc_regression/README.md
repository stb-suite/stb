# SIESTA regression of the bundled PseudoDojo bank

Needs SIESTA (5.4.2 was used to choose the bank) and MPI, so it is run by hand on a cluster whenever the
bank or the blocks the suite writes for it (`core/pseudopotentials.py`) change. `../test.sh` covers what
needs no SIESTA.

```bash
python make_cases.py -o cases              # ~870 cases, every element of every table
bash run_cases.sh cases 1 600              # serial
bash run_cases.sh cases 8 600              # 8 MPI ranks
python summarize.py cases                  # table per bank and rank count, cases/results.tsv
```

Each case is an X-C molecule (C-C for carbon) in a 12 Ang box, DZP, `PAO.EnergyShift 0.02 Ry`,
`Mesh.Cutoff 420 Ry`, 3 SCF iterations, `SCF.MustConverge false`; it passes when SIESTA writes
`0_NORMAL_EXIT`. Full-relativistic tables run with and without spin-orbit (`fr-pbesol` without it only
Ge, Sb and Te, as X-X). The lanthanides of `sr3plus` use C from `dojo`.

Expected: every case OK with 1 and 8 ranks. Checked locally on 2026-11-20 (SIESTA 5.4.2, 1 rank) for
Ce, Gd, Lu (sr3plus), Be, Li, Pb, S, I (dojo), Xe (dojo-stringent), Ge dimer and Fe with spin-orbit
(dojo-fr-pbesol), Pb (dojo-fr): 12/12; Ce without the synthetic-atoms block stops as expected. Failure categories in `summarize.py`: `polarization` (Li/Be
without the non-perturbative scheme), `ekb0` (`propor: IMAX = 0`), `scalar`, `basis-charge` (sr3plus
without synthetic-atoms or at a wrong index), `ranks` (more ranks than the molecule can use), `nan`,
`timeout`, `other`.
