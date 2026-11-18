# Changelog

What changed in each version of STB-SUITE, newest first. Under each version: **Fixed** (a wrong
result, a crash, a misleading message), **Changed** (a default or a behaviour that differs from the
version before) and **Added** (new options and tools).

## 2.0.1 (unreleased)

### Fixed

- `stb-bader` did not start on a fresh install ("pybader library not found", also for `--version`):
  pybader 0.3.12 imports `pkg_resources`, which setuptools 81 and later no longer ship. The suite
  now gives pybader what it needs from the standard library, so any setuptools version works. When
  an installed dependency fails to import, the message shows the original error instead of telling
  you to install it again.
- Elastic constants (5.2), `--symmetry-method full` and `--method energy`: a cubic crystal with
  numerical noise in its cell (a relaxed structure, 3e-6 Ang) got one strain direction, `xx`, instead
  of two, and `stb-elasticAnalysis` then printed a shear constant with no meaning (Si: C44 = -6288 GPa,
  E = 941 GPa, Poisson -1.26) with zero quality warnings.
  - `stb-elasticInputs` now builds the symmetry operations from the cell with the symmetry made exact,
    and counts a constant as determined by a direction only when that direction reaches it through a
    response no more than 100 times weaker than the direct ones. The choice of directions no longer
    depends on the noise (tested up to 1e-3 Ang on cubic, hexagonal and triclinic cells).
  - `stb-elasticAnalysis` stops with an `[ERROR]` when the folders that were run leave open a constant
    it reports, names those constants and the directions to add. Before, this was a warning in 3D,
    silent in 2D, and the tensor was printed anyway. Constants the report never shows (out of plane in
    2D) may stay open, as before.
  - Folders made by 2.0.0 for a cubic crystal with only `xx`: run `stb-elasticInputs --dirs xy` (or
    `yz`) for the missing shear direction; the `xx` folders stay valid.
- `stb-slab`: the cell's c vector was in general not perpendicular to the surface (fcc primitive
  cell, (111): 35.3 degrees off), so the surface was tilted with respect to the xy plane, while
  everything that follows takes z as the normal: `stb-adsorb` places the adsorbate along z, SIESTA's
  slab dipole correction acts along the vacuum axis, and the thickness and vacuum in the report were
  measured along z (Ag(111), 5 layers: 11.70 Ang printed, 9.55 Ang real).
  - The slab is now written with c along the surface normal, a along x and b in the xy plane. Only c
    changes, to its component along the normal: the atoms, the in-plane cell and the vacuum are the
    same (checked on fcc, diamond, quartz and a triclinic cell: every interatomic distance is kept).
    The report lists the three cell vectors and the angle the generated c made with the normal.
  - Thickness and vacuum are measured along the normal, also in the termination table of `-i`.
  - A slab whose generated cell already had c along the normal is written exactly as before. For any
    other, the file differs from the one 2.0.0 wrote: same slab, different cell.

- `stb-adsorb` (also `stb-mladsorb` and every tool that detects vacuum): a slab with exactly the
  threshold of vacuum (10 Ang by default) was rejected as "no vacuum axis", because the gap read back
  from the file is 9.99999 Ang. Gaps within 0.01 Ang of the threshold now count, and the error of
  `stb-adsorb` gives the gap measured along each axis and the `--vacuum-gap` value to pass.
- `stb-adsorb` with a `--calc` that has no `kgrid.MonkhorstPack [nx ny nz]` line ended in a Python
  traceback; it is a clean error now.

- DFT-D3 and the stress tensor: SIESTA 5.4.2 writes the D3 contribution to the stress with the wrong
  sign (Ag and Cu, 13 volumes: P(D3) - P(PBE) = +4.75 GPa where -dE_D3/dV = -4.75 GPa; energies and
  forces are right). The defect is SIESTA's, but the suite used that stress with no notice. It now
  warns wherever the stress of a D3 run is used:
  - `stb-mlffAnalysis` trained on the stress by default (weight 10000) with wrong stress labels: it
    warns and recommends `--no-stress-training` when the labels include D3 and the logs are from an
    affected version (or the version cannot be read).
  - `stb-eosAnalysis`: the pressure cross-check `[5]` says that P_SIESTA and the B0 from the stress are
    not reliable for these runs, and that V0 and B0 from the energies are.
  - `stb-elasticAnalysis --method stress`: warns and points to `--method energy`.
  - `stb-convergence` (every folder relaxes the cell) and `stb-inputfile -t relax -d3`: warn that the
    cell relaxes to the wrong volume with that version.

- MACE tools on a CPU-only node with an old NVIDIA driver: torch's `UserWarning: CUDA initialization:
  The NVIDIA driver on your system is too old` no longer reaches the output or the report's library
  warnings (it made the `stb-mlmd` test fail). A run that asks for the GPU still gets its own error.

- `stb-strainAnalysis` said "slope under a clamped transverse cell" and "transverse stress
  (diagnostic, clamped cell)" also for runs made with `stb-strain --relax-mode stress-constrained`,
  where the transverse cell relaxes and the slope is Young's modulus. It now reads the mode from each
  folder's `config_extra.fdf`, words the report accordingly, and for the relaxed mode gives Poisson's
  ratio from the transverse strain of the final cells (`*.STRUCT_OUT`).
- `stb-cohesive`: the log called the whole `--vacuum` a "buffer beyond the outermost atom" and, lines
  below, half of it; the first line now says what the value is (the vacuum gap between images, box
  side = 2 x cutoff + vacuum) and the buffer per side.
- Example 5.3 (cohesive energy): the orbital radius `rc = 7.12` read from `calc.out` is in Bohr (3.77
  Ang), not Angstrom, so the old box did not let periodic images overlap; the text is corrected. The
  README now shows the ghost-cutoff scan up to 12 Ang: the corrected GaAs value at 6 Ang (-3.1741
  eV/atom) is not converged (-3.0514 at 12 Ang).

- `stb-mlelastic` on a 2D material (or a 1D one, or 3D with missing directions): the clamped-ion scan
  was run (twice the evaluations) but its table, `INTERNAL RELAXATION (clamped vs relaxed ions)`, was
  printed only for a complete 3D tensor. It is printed for every dimensionality now.

- `stb-opticalAnalysis` called graphene "dielectric (non-metallic)" with an absorption edge at 0.25 eV
  when the run folder had no `<label>.EIG` to read the Kohn-Sham gap from: the edge came from the
  finite k-mesh and the 0.1 eV broadening. An onset below 3 x the broadening (at least 0.3 eV), or a
  Kohn-Sham gap below the broadening, is now reported as "gap not resolved": no dielectric
  classification, no Tauc gap, and a warning that the spectrum below about twice the broadening is
  not reliable. A missing `.EIG` is said.

- `stb-mlphonons --qha --n-volumes 4` was accepted and failed at the very end, after every volume had
  been computed (`At least 5 volume points are needed to fit S(V)`). The minimum is 5 and the option
  is checked before anything runs.
- `stb-eosAnalysis` warned `B0' within 4.5-7.0` for ordinary fits (Si: 4.25-4.47 with R^2 0.999998), so
  the warning appeared in almost every run. The accepted range is 3.0-7.0, and the message says what a
  value outside it means (a scan too narrow or noisy to fix B0'), not that the fit is wrong.
- Test of `stb-crystalcast --ml-rank` without the `ml` extra: the tool does exit with code 1 and the
  install hint; the test read the exit code after an `if` and an `echo`, so it always saw 0.

- `stb-nebCycle` on a band whose images coincide (identical initial and final states) wrote the next
  cycle and ended with `Success` and exit code 0: its guard only caught NaN forces, which ASE does not
  produce in that case. The geometries are now checked before any step: neighbouring images closer
  than 0.001 Ang stop the cycle with an `[ERROR]` and exit code 1, and no folder is written.
- `stb-mlneb` with identical endpoints (or two inputs that relax into the same minimum) ended in
  `numpy.linalg.LinAlgError: Singular matrix`. It now stops with a message: nearly identical states,
  no barrier along this band.
- `stb-mlneb --siesta-check/--to-mlff --mlff-workflow`: a workflow whose reference structure file is
  missing gave a Python `FileNotFoundError`; it is a clean error naming the file. (In the test suite
  the file was missing because three tests still pointed to fixture folders by their numbers from
  before the 2026-10-31 renumbering; fixed.)

- `stb-mlmd`: the initial velocities are rescaled to exactly the temperature asked for (a
  Maxwell-Boltzmann draw of 72 atoms gave 263.5 K for 300 K). Trajectories with the same seed
  therefore differ from 2.0.0's. The lines `Initial state` / `Final state` printed a potential energy
  as `E = ...`, which read as a 1.27 eV change of the total energy in an NVE run: they now give E_pot,
  E_kin and E_tot. An NVE run from an equilibrium geometry whose mean temperature ends below 75 % of
  `--temperature` (equipartition halves it) gets a warning in the final section, with what to do.
- `stb-phononsPos` (and every loader of SIESTA phonons: IR, Raman, adsorption Gibbs, `stb-mlffValidate
  --phonon-ref`): the space group phonopy used to build and check the force constants was not always the
  one at `--symprec`. On the cluster, a relaxed 5x5 graphene (positions 0.66667844 instead of 2/3) was
  analysed as Cmmm, the tolerance written in `phonopy_disp.yaml` (1e-5); the force constants built from
  the single P6/mmm displacement split the E2g pair (45.55 and 45.96 THz, with or without
  `--symmetrize-fc`) and the force-constant symmetry check reported 0.000 THz.
  - The space group is now computed with spglib at `--symprec` and compared with phonopy's; on a
    mismatch the phonopy object is rebuilt with the requested tolerance and the force constants are
    recomputed (with a `[WARNING]`; the NAC parameters read from `BORN` are kept). The report gives the
    space group used, the tolerance written in `phonopy_disp.yaml`, and the number of displacements in the
    folders against the number that group needs, with a warning when there are fewer (the set was reduced
    with a higher symmetry than the one now used). Same data: E2g 38.575 / 38.576 THz.
  - The force-constant symmetry check of `[2]` was skipped with phonopy 4.5 ('index 2 is out of
    bounds'): it needs the full force constants, which are now built for it.
- `stb-convergence` (relaxation mode, the default): `config_extra.fdf` forced a variable-cell relaxation
  but no stress tolerance, so SIESTA's default of 1 GPa applied -- about 1 % of the volume of Si (B0 ~
  88 GPa), twice the 0.5 % the analysis judges by. On Si (DZP, 350 Ry) the relaxations from 0.0005 to
  0.01 Ry stopped at their first step (0.96 GPa left at 0.0005 Ry) and the analysis tabulated the input
  volume as the relaxed one.
  - Every folder now gets `MD.MaxStressTol` (new `--max-stress-tol`, default 0.02 GPa) and
    `MD.MaxForceTol` (new `--max-force-tol`, default 0.01 eV/Ang); `--single-point` is unchanged.
  - `stb-convergenceAnalysis` reads the tolerance back and the final stress from each log, and flags
    `P>tol` (stress left above the tolerance), `Ptol=1GPa` (no tolerance below SIESTA's default: folders
    from earlier versions) and `V=input` (the cell never moved under a loose tolerance). Those runs are
    left out of the volume verdict, with a `[WARNING]`; when every run is one of them there is no verdict.
    Regenerate such folders and rerun them.
- `stb-elasticAnalysis`, folders prepared with `stb-elasticInputs --symmetry-method full` and analysed
  with the default `--symmetry-method basic`: the constants of the directions 'full' had left out came
  out 0 with no error (Si with `xx` only: C44 = 0; graphene and hBN with `xx` only: C22 = C66 = 0), and
  the report ended in "VERDICT: UNSTABLE".
  - `stb-elasticInputs` writes `elastic_scan.json` (symmetry method, directions) into every series
    folder, always. Without `--symmetry-method`, `stb-elasticAnalysis` uses the method found there (or
    in the `Symmetry method` line of `elastic_stage1.txt`, for folders of earlier versions); the report
    says where it was read. With no record the default is still `basic`.
  - With `basic`, a constant the point group allows that comes out 0 because neither of its directions
    was run nor reconstructed is now an `[ERROR]` (no tensor is reported). It says whether the
    directions run determine the tensor with `--symmetry-method full`, or which directions to add.
  - The `[ERROR]` of `full` for an undetermined constant proposes the direction `stb-elasticInputs`
    itself would choose (`xy` for cubic, was `yz`).
  - The menu (5.2.2) no longer answers `basic` when Enter is pressed at the symmetry question: Enter keeps
    the method the folders were prepared with; `1` / `2` still choose `basic` / `full`.
- `stb-cohesive --calc`: a template that sets its k-points with `kgrid.Cutoff` (or `kgrid.File`) was
  taken as having no k-grid, so the full structure's `config_extra.fdf` got a Monkhorst-Pack block from
  the default density, which SIESTA reads first (E05 template, `kgrid.Cutoff 15.0 Ang`: a 10x10x10 grid,
  560 k-points). These count as the template's own k-grid now, as the help says; `-k` still replaces it.
- `stb-convergenceAnalysis --parameter basis`: the orbital filter was not detected in SIESTA 5.4.2 logs
  (Cu with `FilterCutoff 150 Ry`: the prep warned, the analysis did not). The log check now also knows
  5.4.2's `Filter: Requested cutoff` and `PAO Filter: Cutoff used` lines.
- `stb-mlffValidate --phonon-ref` on a 2D material: the default q-mesh was 8x8x8, so the frequencies
  were also compared at q along the vacuum, where the dispersion is an artefact of the periodic images
  (graphene 3x3: the largest deviation was at q = (0.3125, 0.3125, 0.0625)). The default `--phonon-mesh`
  now has 1 along every vacuum axis of the reference cell (8x8x1 for a slab or a sheet); the report says
  so. An explicit `--phonon-mesh` is used as given.
- `stb-convergenceAnalysis --parameter basis` on a 2D material (or a wire) tabulated V0 in Ang^3 and B0
  in GPa of the cell with its vacuum, with a note that only their relative change meant something. It now
  gives what `stb-eosAnalysis` gives for the same kind of scan: the equilibrium area A0 (Ang^2), a0 (the
  first periodic vector at A0) and the in-plane stiffness A d2E/dA2 in N/m -- (C11 + C12)/2 of an
  isotropic sheet; for a wire the length and the axial stiffness in nN. The criteria are unchanged.
- `stb-cohesiveAnalysis`: the warning on the electronic entropy (-T*S above 1 meV) explained every term
  as an isolated atom's fractional occupation, also when the only such term was the solid's (Ag: -1.0 meV,
  the smearing at the Fermi level of a metal). Atoms and solid now get separate messages: a `[WARNING]`
  for the isolated-atom/BSSE runs as before, a `[NOTE]` for the full structure that calls it the usual
  smearing correction.
- `stb-cohesive --calc` showed `FilterCutoff` in its level-of-theory line but gave no warning, unlike
  `stb-eosInputs`, `stb-elasticInputs` and `stb-convergence`; it now gives the same `[WARNING]`.
- `stb-convergence` (and every tool that puts `%include config_extra.fdf` on top of a template): a
  template that already had the line got a second identical one. Harmless, but it read as a mistake; the
  file now keeps one, at the top.
- `stb-convergenceAnalysis --parameter basis`: the command it prints for the next stage had no `-p`, so
  the folders of that stage would have been written without pseudopotentials. `stb-convergence` now
  records the pseudopotential source in `basis_scan_<stage>.json` (as typed: a bank name stays a bank
  name) and the command carries it; for a scan written by an earlier version the command shows
  `-p <pseudopotential source>` to fill in.
- `stb-convergenceAnalysis --parameter basis`, mesh and range stages: a setting closed when its V0 agreed
  with the next TWO tighter ones, which accepts a false plateau. On Ge (500 Ry) V0 was 24.517 / 24.426 /
  24.464 / 25.180 / 25.215 / 25.106 Ang^3/atom from 0.02 to 0.0005 Ry: 0.02 Ry agreed with 0.01 and 0.005
  within 0.4 % and was chosen, though V0 jumps +3 % at 0.002 Ry. A setting now has to agree with EVERY
  tighter setting of the scan (Ge: 0.002 Ry). On the nine test systems: Si and graphene 0.005 Ry; Ge,
  NaCl, Al, Cu and hBN 0.002 Ry; GaAs and Ag 0.001 Ry. In 2D the verdict speaks of the area.
- `stb-mlmd` (and every tool that prints the model card): the settings line said `DFTD3 default` for a
  model trained without D3, which reads as D3 on. It now says `DFT-D3 on` or `DFT-D3 off`, from the
  card's own `dftd3.on`.

- `stb-coop --quantity cohp` called every bond "antibonding": the report applied the COOP rule
  (positive = bonding) to the COHP too, while sisl's COHP uses the raw H(k), Dronskowski's convention,
  in which a bond is NEGATIVE (-ICOHP > 0). On the Sn3O4 example a Sn-O bond (ICOHP -1.44 eV) was
  labelled antibonding; on the X-C bonds of XC8, ICOHP -1.87 eV. The curve and `cohp.dat` were right;
  the label, the sign-convention sentence and the plot title were not. The bond-order cross-check now
  runs for COHP too, comparing -ICOHP with the bond order.

- `stb-bader` stopped with "[ERROR] Analysis failed: object of type 'NoneType' has no len()" when
  the whole density forms a single basin (AlP with a valence-only pseudo-density: the Al has no
  density maximum of its own). PyBader then finds no basin edge and returns no surface distances.
  The report is now written, Al 0 and P 8 e- as the Henkelman program gives, with the
  zero-population warning and no surface-distance column.

### Changed

- `stb-convergence --parameter energyshift`: the default sweep is the geometric set 0.02, 0.01, 0.005,
  0.002, 0.001, 0.0005 Ry (`--energyshift-values` for another). The old linear default, 0.001 to 0.05 in
  steps of 0.01, had no point between 0.001 and 0.011 Ry, the range in use, while the parameter acts
  logarithmically. `--energyshift-range MIN MAX STEP` still gives a linear sweep.
- `stb-convergenceAnalysis` no longer judges a PAO.EnergyShift sweep by the total energy, which has no
  practical plateau in this parameter (0.01 to 1.5 eV/atom between 0.011 and 0.001 Ry, not monotonic for
  graphene): the verdict was always "not converged" or the smallest value. It uses the relaxed volume,
  says so, and tabulates the smallest and largest orbital cutoff radius of each run.

- **The generated `calc.fdf` no longer turns the orbital filter on.** The templates of `stb-inputfile`
  (all four types) and `stb-cohesive` wrote `FilterCutoff 150 Ry`, which makes SIESTA filter the basis
  orbitals (effective cutoff 73.5 Ry): the basis changes with it, and E(V) gets noisy. Cu and Al (PBE,
  DZP, 9 volumes): fit residual 2.13 and 1.43 meV/atom with the line, 0.11 and 0.07 without; Cu V0
  11.715 -> 11.904 Ang^3 (all-electron 11.951) and B0 128 -> 141 GPa (141.3); Al B0 71 -> 80 GPa (78.1).
  Doubling `Mesh.Cutoff` or the k-grid did not help. The line is now written commented out.
  - **Numbers change** for every input generated from now on (large for Al and Cu, medium for Ag, none
    for Si, graphene, hBN).
  - Without the filter `Mesh.Cutoff` has to be converged per system: 350 Ry is enough for Si, NaCl, Al,
    Cu, Ag, graphene and hBN, but Ge and GaAs need about 500 Ry (V0 at 350 vs 1000 Ry: GaAs +0.56 %).
  - A template that still has `FilterCutoff`/`FilterTol` gets a warning from `stb-eosInputs`,
    `stb-eosAnalysis` (which also names the filter among the causes of a large fit residual),
    `stb-elasticInputs` and `stb-convergence`.

- **`stb-cohesiveAnalysis` uses the Kohn-Sham energy E_KS, not SIESTA's `FreeEng`.** `FreeEng` includes
  -T*S of the electronic smearing, and an open-shell isolated atom with degenerate orbitals has
  fractional occupations at `ElectronicTemperature 300 K`: -49.4 meV for Si (p2) and Ga (p1), 0 for As
  and Ag, 0 for a gapped solid. The cohesive energy then depended on a smearing temperature that is not
  physical (Si: -4.5138 with FreeEng, -4.5632 eV/atom with E_KS; GaAs 25 meV/atom).
  - **Numbers change**: cohesive energies of systems with such atoms come out more negative than in
    2.0.0. `--energy free` gives the 2.0.0 values; `--energy zero` the T -> 0 extrapolation
    (E_KS + FreeEng)/2.
  - New section `[2b] ENERGIES USED`: the convention, and E_KS, FreeEng and -T*S of every run, with a
    warning when a reference has |T*S| above 1 meV.

- `stb-stackingfault --scan surface`: the default grid is 6 x 6 on a hexagonal lattice (7 x 7 on any
  other, as before). The AB/BA registries of a hexagonal bilayer are at shifts (1/3, 2/3) and
  (2/3, 1/3), which a 7 x 7 grid does not contain: the reported minimum was 0.8 meV above the exact AB
  and the corrugation 3.62 instead of 3.78 meV/Ang^2. A grid that is not a multiple of 3 on such a
  lattice now warns, in the prep and in `stb-stackingfaultAnalysis`. `-ny` defaults to `-nx`.

- `stb-phononsPos` and `stb-mlphonons` use a Gamma-centred q-mesh. phonopy shifts an even mesh by half
  a step, and on a hexagonal lattice that mesh loses the point group (`-m 30 30 1`:
  `MeshSymmetryFallbackWarning`, only time reversal used). Thermal properties change in the last
  digits (graphene, 30x30x1: zero-point energy 170.64 -> 170.63 meV/atom). `stb-phononsPos
  --shifted-mesh` gives the old mesh.

- Supported Python versions are 3.10 to 3.12 (`requires-python` now says so). 3.9 was listed but never
  worked: phonopy does not build there, so `pip install` failed late with an unrelated-looking build
  error; pip now refuses at once with the reason.

- **Elastic constants (5.2) are relaxed-ion by default.** `stb-elasticInputs` used to force a single
  point in every folder, so the workflow could only give clamped-ion constants, while `stb-mlelastic`
  relaxes the atoms: the same strain gave two different quantities with no notice. In a crystal with
  an atom off an inversion centre (diamond, zincblende, graphene, hBN, wurtzite, perovskites) the
  strain leaves forces on the atoms (Si `xy` at 2 %: 0.39 eV/Ang) and the measured constant is the
  relaxed one (Si C44: 102 GPa clamped, 75 GPa relaxed; graphene C12 +23 %).
  - `stb-elasticInputs --ions relaxed|clamped|both` (default `relaxed`): the cell stays fixed at the
    imposed strain; `relaxed` writes a CG relaxation into `config_extra.fdf` (`--relax-steps`, 200;
    `--fmax`, 0.01 eV/Ang), `clamped` the single point 2.0.0 wrote, `both` the two series in
    `<output-dir>/relaxed` and `<output-dir>/clamped`. It holds for `--method stress` and
    `--method energy`. Also asked by the menu (5.2.1).
  - **To get the inputs 2.0.0 generated, pass `--ions clamped`.**
  - `stb-elasticAnalysis` has a new section `[0b]`: whether the folders are clamped- or relaxed-ion
    (read from each log), the CG moves and the last maximum force per direction. It warns when a
    relaxation stopped before its force tolerance, when the folders mix the two kinds, when a folder
    ran with a variable cell, and, for clamped-ion runs, when a strained cell shows a force above
    0.05 eV/Ang (the constants are then not the measured ones). The summary states the kind.
  - It writes the fitted tensor to `elastic_constants.json` in the run directory and, for the two
    series of `--ions both` (or `--compare-with DIR`), prints relaxed minus clamped for every
    reported constant: the internal-relaxation contribution.
  - The stress and the energy were already read from the last geometry of each log, so folders
    relaxed by hand (a replaced `config_extra.fdf`) give the same numbers as before.

- `stb-elasticAnalysis --symmetry-method`: the default is now the method the folders were prepared with
  (from `elastic_scan.json` or `elastic_stage1.txt`), `basic` only when neither exists. The menu (5.2.2)
  does the same on Enter.
- `stb-mlffValidate --phonon-mesh`: the default is 8 per periodic axis and 1 per vacuum axis (8x8x1 for a
  slab or a 2D material; was 8x8x8). Phonon RMSEs on the mesh of 2D systems change.
- `stb-convergenceAnalysis --parameter basis`: a setting of the mesh or range stage closes only when its V0
  agrees with every tighter setting of the scan (was: the next two); verdicts can move to tighter values.

### Added

- **PseudoDojo tables by relativity, functional and accuracy.** The bundled `dojo` bank now holds
  eleven PseudoDojo tables (PSML): scalar-relativistic PBE, PBEsol and LDA (PW92), full-relativistic
  PBE and PBEsol (for spin-orbit), each standard and stringent, plus the lanthanides with 4f in the core
  (Ce-Lu, PBE standard). Every tool's `-p`/`--pseudo-dir`/`--pp-path` picks one by name,
  `dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent]`, words in any order (`fr` may be
  written `soc`, `pw` `lda`), a missing word taking its default (sr, pbe, standard): `-p dojo-pbesol`,
  `-p dojo-fr`, `-p dojo-stringent`. A table that does not exist (`dojo-fr-pw`: there is no
  full-relativistic LDA table) is a clean error that lists the ones that do. `dojo` alone is
  `dojo-sr-pbe-standard`, the same files as before, so existing commands, scripts and
  `mlff_workflow.json` files give the same result. The menu asks for the relativity, the functional and
  the accuracy after "PseudoDojo", offering only the combinations that exist (Enter keeps the default).
  `stb-pseudo --list-banks` lists every table, `--list-elements` accepts any table name, and the report
  gives each file's functional, relativistic treatment and valence (`z=`, e.g. Ge 14 standard, 22
  stringent); `references.bib` cites PseudoDojo for every table.
- **Pseudopotentials vs the calculation: warnings.** Every tool that copies pseudopotentials into run
  folders now reads the functional and the relativistic treatment written in each `.psml`/`.psf` and
  compares them, when it finishes, with the folder's own `.fdf` (`%include` files followed,
  `config_extra.fdf` included): a `[WARNING]` when the input's functional (`XC.Functional`/`XC.Authors`,
  or SIESTA's LDA default when none is set) is not the files' one, when `Spin spin-orbit` meets files
  with no lj projectors (SIESTA then stops with 'The PSML file does not contain lj projs'), and when one
  folder mixes functionals or full- and scalar-relativistic files. Each warning is printed once, naming
  the first folder and how many more share it. Nothing is blocked: the files are copied as asked.

- `stb-convergence --parameter basis` and its report in `stb-convergenceAnalysis`: basis-set convergence
  by equations of state, in three stages (`--basis-stage mesh | range | size`). Each setting is judged by
  the equilibrium volume (2D: area) of an E(V) scan of 9 single points, not by the total energy:
  `mesh` scans Mesh.Cutoff 350-1000 Ry at PAO.EnergyShift 0.002 Ry, `range` scans PAO.EnergyShift
  0.02-0.0005 Ry at `--mesh-cutoff`, `size` compares the template's PAO.BasisSize with the next one at
  `--mesh-cutoff` and `--energyshift`. The analysis fits every scan, gives the verdict of each stage (the
  first setting within 0.2 % / 0.5 % / 1 % of the next two), lists the orbital radii and the fit
  residual, warns about an active orbital filter and a range that does not close, and prints the command
  of the next stage. `--reference-v0` compares with an all-electron volume.
- `stb-convergence --single-point`: every folder is a single point at the input geometry, so that only
  the swept parameter changes (by default each value relaxes positions and cell, which mixes the
  parameter's effect with the relaxation's noise).
- `stb-eosInputs` / `stb-eosAnalysis` accept 2D materials and wires. A structure with vacuum was
  refused ("bulk only"), and area scans had to be built by hand. Now only the periodic lattice vectors
  are scaled (`--strain-range` is an area strain for a 2D material, a length strain for a wire; the
  vacuum vector is not touched), and Stage 2 reports the equilibrium area and in-plane lattice vectors
  with the in-plane stiffness `A d2E/dA2` in N/m (wire: length, axial stiffness in nN). Its pressure
  cross-check uses the stress conjugate to the scan. Bulk structures are handled exactly as before.
- `stb-adsorb --optimizer cg|broyden|fire`, `--force-tol`, `--max-displ`: relaxation settings written
  into every folder's `config_extra.fdf` (default: whatever `--calc` sets). A weakly bound molecule on a
  metal is ill-conditioned for SIESTA's CG: H2O on Ag(111) 3x3 was still at 0.03-1.3 eV/Ang after
  130-146 steps (11 h on 32 cores). With a molecular adsorbate and a CG template the tool now says so.
- `stb-adsorbAnalysis` warns about every run (clean slab or site) that ended at the step limit without
  reaching the force tolerance. SIESTA ends such a job normally (`0_NORMAL_EXIT`), and the energy was
  used as if it were a minimum.
- `stb-cohesive -c/--calc FILE`: your own `calc.fdf` template for the full structure and every atomic
  reference, so the whole study runs at one level of theory. Before, the tool wrote a fixed template
  (`PAO.EnergyShift 0.02 Ry`, `Mesh.Cutoff 320 Ry`) and any other choice meant editing every generated
  `calc.fdf` by hand (127 files for 13 solids in one study), with nothing recorded in the report. The
  template is copied untouched; what the workflow must impose goes into a `config_extra.fdf` included on
  top of it: Spin polarized, a Gamma-only k-grid and a single point for the atomic references; `--spin`,
  `--dispersion` and the k-grid (only when the template has none, or `-k` is given) for the full
  structure. The report lists the template's level of theory. Without `--calc` the generated files are
  identical to 2.0.0's. Also asked by the menu (5.3.1).
- `stb-phononsPos`: section `[2]` reports how far the force constants are from the crystal's symmetry
  (largest change of a Gamma frequency when they are symmetrized by the space group) and warns above
  0.05 THz, the sign of modes split that symmetry makes degenerate. `--symmetrize-fc` applies the
  symmetrization.
- `stb-mlffValidate --phonon-ref`: the phonon table now has a row for the foundation model next to
  every round (it was given for the fine-tuned model only), says where the largest deviation is
  (q-point and band), and gives the error twice, labelled: (a) on the q-mesh and (b) along the band
  path a band plot shows. The two differ for the same model, and the report did not say which it was.
- `stb-adsorb --fix-layers N` (`--layer-tol`): holds the N lowest atomic layers of the slab fixed, with
  the same `%block Geometry.Constraints` in `clean_slab/` and in every site folder. Before, the whole
  slab always relaxed. Also asked by the menu (5.8.1).
- `stb-adsorb --anchor-atom ATOM` (element or atom number): the atom of a molecular adsorbate that
  faces the surface; `--height` is then the height of that atom. Without it H2O is placed with the
  hydrogens towards the surface, and the orientation could only be changed through `--ml-rank`.
- `stb-slab --layers N`: the slab's size as a number of atomic layers instead of a minimum thickness
  in Angstrom, for a convergence series in layers (3, 4, 5, ...). The slab is cut thick enough and
  layers are removed from its top until N are left; the bottom face is the chosen termination and the
  vacuum is the same for every N. `--layer-tol` (default 0.1 Ang) sets which atoms count as one
  layer: a buckled plane, like the two halves of a Si(111) bilayer, is two layers unless the tolerance
  is raised above its buckling. Also asked by the menu (4.3).
- `stb-slab` reports the number of atomic layers of every slab it writes, the atoms in each and the
  spacing between them.
- `stb-adsorbAnalysis`: the orientation of the adsorbate at the start and at the end of each site's
  relaxation, in section `[2b]` (tilt of the molecular plane, or axis, from the surface plane; lying /
  tilted / upright; the atom closest to the surface; vertical extent). Sites whose orientation CHANGED are
  marked, the final orientations reached are listed with the lowest E_ads of each, and a warning says when
  more than one was reached. In a 12-point H2O/Ag(111) series, 4 relaxations had ended lying flat, ~0.1 eV
  more bound, and the ranking mixed them with upright ones without saying so.
- `stb-convergence --max-stress-tol GPA` (default 0.02) and `--max-force-tol EV_PER_ANG` (default 0.01):
  `MD.MaxStressTol` and `MD.MaxForceTol` of the relaxation mode (see Fixed).
- `stb-elasticInputs` writes `elastic_scan.json` (symmetry method, directions) into every series folder;
  `stb-elasticAnalysis` reads the method from it (see Fixed).

## 2.0.0 (2026-10-31) -- first stable release

- Every tool is version 2.0.0. The `main` branch holds the stable releases (2.0.1, ... for bug
  fixes); development of the next version happens on `dev`, which reports `2.1.0.dev<N>` (N =
  commits since the `v2.0.0` tag). See `src/stb/_version.py`.
- The menu has seven sections with no gaps: 1 Inputs, 2 Analysis, 3 Utils, 4 Structures,
  5 Workflows, 6 Catalysis (filled by an extension package), 7 ML Simulations. The old -> new
  code table is in "Menu renumbering" at the end of this file.
- `stb-suite` answers `--help` and `--version` like every other command.
- Raman workflow (5.11), found while building its guided example (`examples/5.11-raman`, real SIESTA runs of
  bulk Si and graphene):
  - `stb-ramanModes` forces an SCF k-grid into the optical folders when the template has none (suggested at
    0.2 1/Ang, or `--kgrid` / `--kgrid-density`, also in the menu); before, those runs sampled Gamma only.
  - The optical block (also `stb-optical`'s) forces `Diag.ParallelOverK false`: with it on, SIESTA 5.4's optical
    calculation aborts on more than one process ("cdiag: Error in Cholesky factorisation").
  - `stb-ramanAnalysis` warns when a diagonal-only run sees a numerically zero diagonal (a T2g mode of a cubic
    crystal has only off-diagonal entries: Si came out "polarized", activity 2e-4 instead of 2891), and checks
    that degenerate partners computed separately have equal activities (they do not for zero-gap graphene,
    whose eps1(0) is not defined). Large numbers no longer run into the next column of the summary table.

# Per-tool history before the versions were unified

Until the suite moved to a single version (`stb.__version__`, `1.9.<number of commits>`, see
`src/stb/_version.py`), every tool carried a version of its own, `VERSION`, in the range 1.0.0 to
3.0.0, and some kept a comment above it describing their last changes. Those numbers are retired.
This file keeps what they said: each tool's last number, and the comments exactly as they were
written. Anything older is in the git history.

## Last per-tool version of every command

| Command | Last per-tool version |
|---|---|
| `stb-2Dstacking` | 2.2.0 |
| `stb-adsorb` | 1.0.0 |
| `stb-adsorbAnalysis` | 1.5.0 |
| `stb-adsorbBsse` | 1.0.0 |
| `stb-adsorbGibbs` | 1.1.0 |
| `stb-aimdAnalysis` | 2.0.0 |
| `stb-amorphize` | 2.1.0 |
| `stb-ani2traj` | 1.0.0 |
| `stb-archive` | 1.0.0 |
| `stb-bader` | 2.0.1 |
| `stb-bands` | 2.0.0 |
| `stb-chargediffAnalysis` | 2.0.0 |
| `stb-chargediffPrep` | 1.1.0 |
| `stb-clean` | 1.9.1 |
| `stb-cohesive` | 2.0.0 |
| `stb-cohesiveAnalysis` | 2.1.0 |
| `stb-convdos` | 2.0.0 |
| `stb-convergence` | 2.2.0 |
| `stb-convergenceAnalysis` | 2.0.0 |
| `stb-coop` | 2.0.0 |
| `stb-crystalbuilder` | 2.0.0 |
| `stb-crystalcast` | 2.0.0 |
| `stb-cube` | 1.9.1 |
| `stb-defect` | 2.0.0 |
| `stb-density` | 2.0.0 |
| `stb-dftu` | 2.0.0 |
| `stb-dos` | 2.0.0 |
| `stb-effmass` | 2.0.0 |
| `stb-elasticAnalysis` | 1.13.2 |
| `stb-elasticInputs` | 2.0.0 |
| `stb-eosAnalysis` | 1.1.0 |
| `stb-eosInputs` | 1.0.0 |
| `stb-fatbands` | 2.1.0 |
| `stb-fetch` | 1.14.1 |
| `stb-gqca` | 1.0.0 |
| `stb-gqcaAnalysis` | 1.0.0 |
| `stb-hirshfeldAnalysis` | 2.1.0 |
| `stb-hirshfeldIons` | 2.1.0 |
| `stb-hirshfeldPrep` | 2.0.0 |
| `stb-hubbardu` | 2.2.0 |
| `stb-hubbarduAlphas` | 1.3.0 |
| `stb-hubbarduAnalysis` | 2.0.0 |
| `stb-inputfile` | 1.9.1 |
| `stb-ipr` | 2.0.0 |
| `stb-ir` | 1.2.1 |
| `stb-irAnalysis` | 1.2.1 |
| `stb-irModes` | 1.4.1 |
| `stb-kgrid` | 2.0.0 |
| `stb-kpath` | 2.0.0 |
| `stb-mladsorb` | 1.1.0 |
| `stb-mlconvergence` | 1.0.0 (removed from the suite; code on the `dev` branch, `experimental/mlconvergence/`) |
| `stb-mldiffusion` | 1.0.0 (removed from the suite; code on the `dev` branch, `experimental/mldiffusion/`) |
| `stb-mlelastic` | 1.1.0 |
| `stb-mleos` | 1.0.2 (removed from the suite; code on the `dev` branch, `experimental/mleos/`) |
| `stb-mlff` | 1.0.0 |
| `stb-mlffActiveLearning` | 1.0.0 |
| `stb-mlffAnalysis` | 1.0.0 |
| `stb-mlgcmc` | 1.0.0 (removed from the suite; code on the `dev` branch, `experimental/mlgcmc/`) |
| `stb-mlmd` | 1.0.0 |
| `stb-mlmelting` | 1.0.0 (removed from the suite; code on the `dev` branch, `experimental/mlmelting/`) |
| `stb-mlneb` | 1.0.0 |
| `stb-mlphonons` | 2.1.0 |
| `stb-mlrelax` | 1.16.0 |
| `stb-mlsearch` | 1.1.0 (removed from the suite; code on the `dev` branch, `experimental/mlsearch/`) |
| `stb-molecule` | 2.0.0 |
| `stb-nanotube` | 2.0.0 |
| `stb-nativecharges` | 1.0.0 |
| `stb-neb` | 1.0.0 |
| `stb-nebAnalysis` | 1.0.0 |
| `stb-nebCycle` | 1.0.0 |
| `stb-optical` | 1.0.0 |
| `stb-opticalAnalysis` | 1.0.0 |
| `stb-passivate` | 2.0.0 |
| `stb-phononsCreate` | 1.13.0 |
| `stb-phononsPos` | 1.14.0 |
| `stb-pseudo` | 1.0.0 |
| `stb-raman` | 1.1.1 |
| `stb-ramanAnalysis` | 1.0.0 |
| `stb-ramanModes` | 1.3.0 |
| `stb-siesta2wtb` | 1.9.1 |
| `stb-slab` | 2.0.0 |
| `stb-spintexture` | 2.0.0 |
| `stb-sqs` | 2.0.0 |
| `stb-stackingfault` | 1.0.0 |
| `stb-stackingfaultAnalysis` | 1.0.0 |
| `stb-status` | 1.0.0 |
| `stb-stm` | 2.1.0 |
| `stb-strain` | 3.0.0 |
| `stb-strainAnalysis` | 2.0.0 |
| `stb-structural` | 2.0.0 |
| `stb-sts` | 2.0.0 |
| `stb-suite` | 1.9.1 |
| `stb-supercell` | 2.0.0 |
| `stb-symmetry` | 2.0.0 |
| `stb-translate` | 1.9.1 |
| `stb-unitcell` | 2.0.0 |
| `stb-wfdensity` | 2.2.0 |
| `stb-workfunction` | 2.0.0 |
| `stb-xrd` | 2.1.0 |
| `stb-xrdrank` | 2.0.1 |
| `stb-xrdsearch` | 1.14.0 |

## Notes recorded next to the version

### `stb-adsorbAnalysis` (last per-tool version 1.5.0)

```text
New: [2b] PHYSICAL DIAGNOSTICS (dipole/spin moment/bond-length change),
[5] SUGGESTED NEXT ANALYSES, and an opt-in [6] GIBBS FREE ENERGY (DG) PREP
step (--compute-gibbs) that writes the displacement folders stb-adsorbGibbs
(the new Stage 4) needs -- see their own docstrings for the design.
```

### `stb-cohesive` (last per-tool version 2.0.0)

```text
stb-standard numbered report ([0]-[7], --save-report) and a
single wrapping --output-dir (default cohesive_runs, matching
stb-strain/stb-elasticInputs' own <property>_runs convention)
instead of writing structure/atoms/atoms_bsse[_check] loose
into the current directory by default.
```

### `stb-cohesiveAnalysis` (last per-tool version 2.1.0)

```text
multi-cutoff BSSE convergence scan: detects and reports/plots
'atoms_bsse_check_<cutoff>/' (2+ points, from stb-cohesive's
--bsse-convergence-increment given multiple values), a new
[3b] section, alongside the existing single-point
'atoms_bsse_check/' comparison (unchanged behavior)
```

### `stb-elasticAnalysis` (last per-tool version 1.13.2)

```text
fix: eggbox cross-check double-applied the Lz/cross_section dilution factor for 2D/1D, inflating energy_value
```

### `stb-eosAnalysis` (last per-tool version 1.1.0)

```text
--eos all, --target-pressure, v0-outside-range warning (core/eos_fit.py)
```

### `stb-fetch` (last per-tool version 1.14.1)

```text
retry transient connection resets against COD (observed live:
crystallography.net intermittently resets the connection on
both /result and /<id>.cif -- not a stb-fetch bug, but worth
a few automatic retries before giving up)
```

### `stb-hirshfeldAnalysis` (last per-tool version 2.1.0)

```text
reads BOTH cation and anion per species and lets each atom
pick its own, per the literal literature formulation
(previously blended toward a single per-species ion sign);
honors stb-hirshfeldIons' "cation_zero_density" flag
(hydrogen's cation is a bare proton -- exactly zero
electron density, not read from a SIESTA folder)
```

### `stb-hirshfeldIons` (last per-tool version 2.1.0)

```text
writes BOTH cation and anion per species unconditionally
(previously picked one sign per species from a pass-0 vote);
skips the cation SIESTA folder (uses an analytical
zero-density reference instead) for any species with
Z_val <= 1 (hydrogen) -- NetCharge +1.0 would otherwise
leave 0 electrons, a bare-nucleus SCF SIESTA cannot converge
```

### `stb-hirshfeldPrep` (last per-tool version 2.0.0)

```text
combined/ folder now written by this stage itself (previously
required an already-computed combined .RHO to exist beforehand)
```

### `stb-ir` (last per-tool version 1.2.1)

```text
New [LIMITATION] note in [1] INPUT STRUCTURE's bulk-path advisory: a
polar bulk crystal's Stage 1 phonon frequencies get no non-analytic
(LO-TO) correction, so they can land between the true TO and LO values
(root-caused, not fixed, on NaCl -- see stb-irModes' VERSION 1.4.1 comment
for the full explanation and literature references). Documented as a known
limitation only; a real fix is planned for a future version.
(previously 1.2.0: --symprec (default 0.01, pymatgen's own default) threaded through to
Phonopy itself, and --kgrid/--kgrid-density auto-suggest + force the
supercell's own k-grid via a new [1] SINGLE-POINT SCF ENFORCEMENT section --
both ported from stb-raman's own fix (VERSION 1.1.1): Phonopy's raw symprec
default (1e-5) is far tighter than any DFT relaxation's real numerical noise
floor and can silently misdetect the true point group (downstream, in
stb-irModes, this drives IR-active/forbidden classification and degenerate-
mode grouping); and the previous "review the k-grid yourself" NOTE left the
supercell's electronic sampling entirely up to the user copying --calc
verbatim, with no enforcement that it was even single-point (a real risk if
--calc is a reused relaxation template) -- see stb-raman's own graphene G-band
case, where k-density (not supercell size) was the dominant convergence lever.
Also: the Stage 1 -> Stage 2 advisory NOTE now distinguishes a 2D slab (one
vacuum axis) from a genuine 0D/1D non-bulk structure, since Stage 2 now
routes a 2D slab through its own new HYBRID path (see ir_modes.py's VERSION
comment) rather than the plain non-bulk dipole-difference path.
```

### `stb-irAnalysis` (last per-tool version 1.2.1)

```text
New [LIMITATION] note at the top of [2] IR-ACTIVE MODES SUMMARY, shown
whenever any mode used the BULK path: those frequencies have no
non-analytic (LO-TO) correction, so for a polar bulk crystal they can land
between the true TO and LO values -- see stb-irModes' VERSION 1.4.1 comment
for the full root-cause story and literature references. Documented as a
known limitation only; a real fix is planned for a future version.
(previously 1.2.0: New HYBRID path (2D slab, exactly one vacuum axis): combines the BULK
path's Born-charge x eigendisplacement formula (in-plane/periodic axes)
with the NONBULK path's dipole-difference formula (vacuum axis) into one
physically valid dmu/dQ per mode -- see stb-irModes' own VERSION comment
for the full root-cause story (a naive dipole moment silently came back
exactly zero for h-BN's in-plane E' mode, a mode the literature puts as
the dominant IR peak). The two contributions are additive by construction
(each is ~0 exactly where the other is valid); axis_mask selects/verifies
this when the vacuum axis's Cartesian alignment can be confirmed, and
falls back to plain addition with a warning when it can't.
Also: [3] SPECTRUM now always lists detected peaks (not just with
--experimental), [2] reports combined intensity for near-degenerate mode
groups (the basis-independent quantity, since individual dmu/dQ directions
within a degenerate subspace are an arbitrary Phonopy eigenvector choice),
the gnuplot .dat/.gplot pair is now opt-in via --save-gnuplot (written
under <directory>/plot/, same convention stb-ramanAnalysis already uses)
instead of always-on, --view opens an interactive matplotlib preview, and
a new [3c] MODE VIBRATIONS section always reloads Stage 1's force constants
to write an extended-XYZ animation per analyzed mode under
<directory>/mode_animations/ (skipped gracefully, never a hard error, if
phonon_disp/ is no longer present) -- --view-modes additionally opens each
one interactively in ASE's viewer. build_mode_animation_frames/
write_mode_animation/phonopy_atoms_to_ase moved to core/phonon_workflow.py
(were duplicated in raman_modes.py/ir_modes.py) now that this module is a
third consumer -- extract-on-second-use, per this suite's own policy.
```

### `stb-irModes` (last per-tool version 1.4.1)

```text
[1b] SYMMETRY ANALYSIS's degenerate-groups listing no longer prints an empty
"modes  (T1u)" line for a group made entirely of ACOUSTIC bands (not in
band_to_k, e.g. NaCl's own T1u translational mode shares its irrep label with
the real IR-active T1u optical mode) -- now built as (ms, mode_ids) pairs and
filtered to drop any group with no displayed (non-acoustic) member.
Also: new [LIMITATION] note after [1]'s per-mode frequency listing, for the
BULK (3D) path only -- these Gamma frequencies have no non-analytic (LO-TO)
correction applied, so for a genuinely POLAR bulk crystal they can land
between the true TO and LO values rather than matching either. Root-caused
(not fixed) on NaCl: T1u computed at 185.6 cm^-1 vs. ~164/~264 cm^-1
(TO/LO, Raunio et al. 1969), landing almost exactly midway -- the textbook
symptom of a finite-supercell frozen-phonon calculation missing the Gonze &
Lee (1997) dipole-dipole subtraction, which needs eps_inf (not currently
computed by any stb tool) alongside the Z* this stage already computes.
Documented as a known limitation, not fixed -- a real fix needs a new
eps_inf-producing SIESTA Optical-module run plus a dynamical-matrix-level
change in core/phonon_workflow.py, planned for a future version.
(previously 1.4.0: --symprec (default 0.01, pymatgen's own default) threaded through to
Phonopy itself, same fix as stb-ramanModes (VERSION 1.2.0): Phonopy's own
raw default (1e-5) is far tighter than any DFT relaxation's real numerical
noise floor and can silently misdetect the true point group -- this is the
actual tolerance driving [1b] SYMMETRY ANALYSIS/--use-symmetry, so a
too-tight symprec here can artificially split a truly degenerate IR-active
mode into separate near-identical frequencies and mislabel a symmetry-silent
mode as IR-active. [1b] SYMMETRY ANALYSIS (space group, point group,
symmetry op count, per-mode Mulliken label/activity) is now always printed,
not just with --use-symmetry, along with an explicit [CHECK] warning to
verify the reported group against the crystal's known symmetry and loosen
--symprec if it looks wrong -- misdetection here is silent, never an error.
Also: every candidate mode being symmetry-forbidden (e.g. a centrosymmetric
crystal with no IR-active Gamma modes at all) now finishes as a [WARNING]
with exit code 0, not a fatal [ERROR]/exit(1) -- it's a valid physical
result, not a misconfiguration. Still a hard [ERROR]/exit(1) when
--modes/--freq-min/--freq-max (not symmetry) are what emptied the selection.
New HYBRID path for exactly-one-vacuum-axis (2D slab) structures. A naive
real-space dipole moment (the old "non-bulk" path, still used as-is for
0D/1D) is only physically valid along a genuinely NON-periodic direction --
confirmed live on monolayer h-BN: SIESTA itself printed an EXACT ZERO
dipole change for the in-plane E' mode's +/-delta displacement (not a
parsing bug), silently reporting zero IR intensity for a mode the
literature puts as the DOMINANT IR peak (~1420 cm^-1). The two in-plane
(periodic) axes of a slab need the same Born-effective-charge/Berry-phase
machinery the "bulk" path already uses; the vacuum axis is still fine via
the old dipole-difference method (verified: the out-of-plane A2'' mode's
dmu_z/dQ came out correct, matching the literature TO/LO range). HYBRID
writes BOTH a shared born_charge_disp/equilibrium/ folder (PolarizationGrids
built vacuum-aware via build_slab_polarization_grid -- the vacuum axis's own
row zeroed, since SIESTA skips the polarization calculation entirely for any
zero grid-count direction) AND per-mode dipole_disp/mode_XX_plus|minus/
folders; stb-irAnalysis combines one component from each. Empirically
validated before writing this code: a manual BornCharge T run on the relaxed
h-BN structure with this exact grid shape gave a sane, sign-correct,
sum-rule-satisfying in-plane Z* (B ~+2.76/+2.78, N ~-2.76/-2.74 |e|, in-plane
sum ~0.01-0.04 |e|, out-of-plane row correctly ~0). Deliberately scoped to 2D
only (exactly 1 vacuum axis) -- 1D (2 vacuum axes) stays on the old non-bulk
path for now, same limitation, no test case yet.
Also: [1] now tags each mode with a Cartesian-polarization character label
(mode_character_label, e.g. "98% z" or "51% x, 47% y", from its
eigendisplacement -- purely descriptive, compare against the Vacuum
axes/Lattice vectors already printed to interpret in-plane vs. out-of-plane
for a given structure), and [1b] always lists degenerate mode groups by
symmetry (band indices sharing one irrep -- Phonopy's own eigenvector basis
choice within such a group is arbitrary, informational regardless of
--use-symmetry/--skip-degenerate).)
```

### `stb-mleos` (last per-tool version 1.0.2)

```text
v0-outside-scanned-range warning surfaced (core/eos_fit.py's fit_eos fix)
```

### `stb-phononsCreate` (last per-tool version 1.13.0)

```text
--symprec (default 0.01, pymatgen's own default) now threaded through
to Phonopy itself -- a real bug: Phonopy's own raw default (1e-5) is far
tighter than any DFT relaxation's real numerical noise floor and was
silently misdetecting the true space group (verified live on a real
relaxed AlP structure: 1e-5 -> wrong R3m, 0.01 -> correct F-43m). Also:
library prints/warnings (MACE/torch/phonopy) captured and moved to a new
final [LIBRARY WARNINGS] section instead of interleaving with the report;
[0b] ML PRE-FLIGHT CHECK gets an explicit Verdict line; symmetry section
rewritten as a print_table with symprec shown and the result highlighted
```

### `stb-raman` (last per-tool version 1.1.1)

```text
--symprec (default 0.01, pymatgen's own default) now threaded through to
Phonopy itself, same fix as stb-phononsCreate: Phonopy's own raw default
(1e-5) is far tighter than any DFT relaxation's real numerical noise floor
and can silently misdetect the true point group -- which, downstream in
stb-ramanModes, artificially splits a truly degenerate Raman-active mode
into separate near-identical frequencies and can mislabel a
symmetry-silent mode as Raman-active. Epilog expanded with a didactic
Notes section spelling out --symprec's real (downstream, stb-ramanModes)
effect and --kgrid-density's role/independence from symmetry -- verified
live on graphene's G-band (Kohn anomaly): >100 cm^-1 shift from doubling
k-density alone, after supercell size (--dim) had stopped mattering.
```

### `stb-ramanModes` (last per-tool version 1.3.0)

```text
--symprec (default 0.01, pymatgen's own default) now threaded through to
Phonopy itself, same fix as stb-phononsCreate: Phonopy's own raw default
(1e-5) is far tighter than any DFT relaxation's real numerical noise floor
and can silently misdetect the true point group -- this is the actual
tolerance driving [1b] SYMMETRY ANALYSIS/--use-symmetry, so a too-tight
symprec here can artificially split a truly degenerate Raman-active mode
into separate near-identical frequencies and mislabel a symmetry-silent
mode as Raman-active. [1b] SYMMETRY ANALYSIS (space group, point group,
symmetry op count, per-mode Mulliken label/activity) is now always printed,
not just with --use-symmetry, along with an explicit [CHECK] warning to
verify the reported group against the crystal's known symmetry and loosen
--symprec if it looks wrong -- misdetection here is silent, never an error.
Also: every candidate mode being symmetry-forbidden (e.g. a centrosymmetric
crystal with no Raman-active Gamma modes at all) now finishes as a
[WARNING] with exit code 0, not a fatal [ERROR]/exit(1) -- it's a valid
physical result, not a misconfiguration. Still a hard [ERROR]/exit(1) when
--modes/--freq-min/--freq-max (not symmetry) are what emptied the selection.
```

## Menu renumbering (2026-10-31)

The menu was reorganized and renumbered so that no section has gaps. New section order: 1 Inputs, 2 Analysis
(was 3), 3 Utils (was 6), 4 Structures (was 2), 5 Workflows (was 4), 6 Catalysis (new: HER and OER, which were
4.13/4.14 and come from an extension package), 7 ML Simulations (was 5). `stb-nativecharges` moved from Utils to
Analysis; the MLFF workflow moved from Workflows to the top of ML Simulations. Workflow stages keep their stage
number (4.1.2 is now 5.1.2, 4.17.3 is 7.1.3). Example folders (`examples/<code>-<name>/`), their scripts
(`example_<code>.sh`) and the test folders (`test/<n>-<section>/<k>-<name>/`) were renamed to match.

| Old code | New code | Item |
|---|---|---|
| 1.1 | 1.1 | Input File Generator (stb-inputfile) |
| 1.2 | 1.2 | Pseudopotential Resolver (stb-pseudo) |
| 1.3 | 1.3 | K-Grid Generator (stb-kgrid) |
| 1.4 | 1.4 | K-Path Generator (stb-kpath) |
| 1.5 | 1.5 | DFT+U / Hubbard Block Generator (stb-dftu) |
| 1.6 | 1.6 | Structure Fetcher (stb-fetch) |
| 1.7 | 1.7 | ML Pre-Relaxation (stb-mlrelax) |
| 3.1 | 2.1 | Bands Analyzer (stb-bands) |
| 3.2 | 2.2 | PDOS XML Parser (stb-dos) |
| 3.3 | 2.3 | DOS Processor (Convolution) (stb-convdos) |
| 3.4 | 2.4 | Structure Analyzer (stb-structural) |
| 3.5 | 2.5 | Symmetry Analyzer (stb-symmetry) |
| 3.6 | 2.6 | Bader Charge Analysis |
| 3.7 | 2.7 | Work Function Calculator |
| 3.8 | 2.8 | Density Plotter (RHO) |
| 3.9 | 2.9 | Powder XRD Simulator (stb-xrd) |
| 3.10 | 2.10 | Fatbands Analyzer (stb-fatbands) |
| 3.11 | 2.11 | STM Simulator (stb-stm) |
| 3.12 | 2.12 | Wavefunction Density (stb-wfdensity) |
| 3.13 | 2.13 | STS Spectroscopy Simulator (stb-sts) |
| 3.14 | 2.14 | COOP/COHP Bonding Analysis (stb-coop) |
| 3.15 | 2.15 | IPR Analyzer (stb-ipr) |
| 3.16 | 2.16 | Effective Mass / Velocity (stb-effmass) |
| 3.17 | 2.17 | Spin Texture Analyzer (stb-spintexture) |
| 3.18 | 2.18 | AIMD Trajectory Analysis (stb-aimdAnalysis) |
| 6.8 | 2.19 | Native Charges (stb-nativecharges) |
| 6.1 | 3.1 | File Translator (stb-translate) |
| 6.2 | 3.2 | Clean File Tools (stb-clean) |
| 6.3 | 3.3 | Grid to Cube Converter (stb-cube) |
| 6.4 | 3.4 | Wantibexos Interface (stb-siesta2wtb) |
| 6.5 | 3.5 | AIMD Trajectory Converter (stb-ani2traj) |
| 6.6 | 3.6 | Status Checker (stb-status) |
| 6.7 | 3.7 | Archive Packager (stb-archive) |
| 2.1 | 4.1 | 2D Monolayer Stacker (stb-2Dstacking) |
| 2.2 | 4.2 | Supercell Builder (stb-supercell) |
| 2.3 | 4.3 | Slab Builder (stb-slab) |
| 2.4 | 4.4 | Nanotube/Nanoribbon Builder (stb-nanotube) |
| 2.5 | 4.5 | Point Defect Generator (stb-defect) |
| 2.6 | 4.6 | SQS Generator (stb-sqs) |
| 2.7 | 4.7 | Unit Cell Finder (stb-unitcell) |
| 2.8 | 4.8 | Crystal Builder (stb-crystalbuilder) |
| 2.9 | 4.9 | Surface Passivator (stb-passivate) |
| 2.10 | 4.10 | Reference Molecule Builder (stb-molecule) |
| 2.11 | 4.11 | Amorphous Structure Generator (stb-amorphize) |
| 2.12 | 4.12 | Random Crystal Generator (stb-crystalcast) |
| 4.1 | 5.1 | Stress-Strain |
| 4.2 | 5.2 | Elastic Constants |
| 4.3 | 5.3 | Cohesive Energy |
| 4.4 | 5.4 | Phonons |
| 4.5 | 5.5 | Convergence Tests |
| 4.6 | 5.6 | Structure Solution (XRD) |
| 4.7 | 5.7 | Hubbard U (Linear Response) |
| 4.8 | 5.8 | Adsorption |
| 4.9 | 5.9 | NEB / Reaction Path |
| 4.10 | 5.10 | 2D Stacking Fault (Gamma-Surface) |
| 4.11 | 5.11 | Raman Spectrum |
| 4.12 | 5.12 | IR Spectrum |
| 4.15 | 5.13 | Generalized Quasi-Chemical Approximation (GQCA) |
| 4.16 | 5.14 | Optical Properties |
| 4.18 | 5.15 | Equation of State |
| 4.19 | 5.16 | Charge Density Difference |
| 4.20 | 5.17 | Hirshfeld-I Charge Partitioning |
| 4.13 | 6.1 | HER (extension) |
| 4.14 | 6.2 | OER (extension) |
| 4.17 | 7.1 | Custom ML Force Field (MACE fine-tuning) |
| 5.1 | 7.2 | ML Molecular Dynamics (stb-mlmd) |
| 5.2 | 7.3 | ML Phonons (stb-mlphonons) |
| 5.3 | 7.4 | ML Elastic Constants (stb-mlelastic) |
| 5.7 | 7.5 | ML NEB (stb-mlneb) |
| 5.10 | 7.6 | ML Adsorption-Site Screening (stb-mladsorb) |
