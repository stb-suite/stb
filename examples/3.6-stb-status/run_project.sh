#!/bin/bash
# Builds and runs the SIESTA calculations of the 3.6 example (project/), all small and real:
#   bulk_sp/        Si single point at the PBE/DZP lattice constant, with bands and the density
#   relax/          Si variable-cell CG relaxation from a = 5.55 Ang -- converges
#   relax_short/    the same relaxation with MD.NumCGsteps 3 -- stops before converging
#   scf_fail/       Si with MaxSCFIterations 3 and a tight tolerance -- the SCF does not converge, SIESTA aborts
#   md/             Si 8 atoms, velocity Verlet from 600 K, 20 steps
#   strain/a_X/     Si single points at a = 5.40, 5.476, 5.55 Ang (a series for relative energies)
#   fe_spin/        bcc Fe, ferromagnetic -- a spin moment and a metal
# SIESTA 5.x + mpirun; about 3 minutes on 4 processes. Pseudopotentials from the suite's bundled bank (not shipped).
#     bash run_project.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
P="$HERE/project"
NP=${1:-4}
BANK="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"

si_fdf() {   # si_fdf <folder> <a (Ang)> <extra fdf lines...>: a 2-atom fcc Si input
    local d="$P/$1" a="$2"; shift 2
    mkdir -p "$d"
    cat > "$d/calc.fdf" <<EOT
SystemLabel      si
SystemName       Si-bulk
NumberOfAtoms    2
NumberOfSpecies  1
%block ChemicalSpeciesLabel
 1 14 Si
%endblock ChemicalSpeciesLabel
LatticeConstant  $a Ang
%block LatticeVectors
 0.0 0.5 0.5
 0.5 0.0 0.5
 0.5 0.5 0.0
%endblock LatticeVectors
AtomicCoordinatesFormat  Fractional
%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
 0.25 0.25 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize       DZP
PAO.EnergyShift     0.01 Ry
Mesh.Cutoff         200 Ry
%block kgrid.MonkhorstPack
 6 0 0 0.5
 0 6 0 0.5
 0 0 6 0.5
%endblock kgrid.MonkhorstPack
XC.Functional       GGA
XC.Authors          PBE
SCF.DM.Tolerance    1.0d-5
EOT
    for line in "$@"; do echo "$line" >> "$d/calc.fdf"; done
}

si_fdf bulk_sp 5.476 "SaveRho true" "BandLinesScale ReciprocalLatticeVectors" "%block BandLines" \
    " 1  0.5 0.5 0.5  L" " 40 0.0 0.0 0.0  \\Gamma" " 45 0.0 0.5 0.5  X" "%endblock BandLines"
si_fdf relax 5.55 "MD.TypeOfRun CG" "MD.VariableCell true" "MD.NumCGsteps 30" "MD.MaxForceTol 0.01 eV/Ang" \
    "MD.MaxStressTol 0.1 GPa"
si_fdf relax_short 5.55 "MD.TypeOfRun CG" "MD.VariableCell true" "MD.NumCGsteps 3" "MD.MaxForceTol 0.01 eV/Ang" \
    "MD.MaxStressTol 0.1 GPa"
si_fdf scf_fail 5.476 "MaxSCFIterations 3" "SCF.DM.Tolerance 1.0d-8" "SCF.Mixer.Weight 0.02"
for a in 5.40 5.476 5.55; do si_fdf "strain/a_$a" "$a"; done

mkdir -p "$P/md"
cat > "$P/md/calc.fdf" <<'EOT'
SystemLabel      si
SystemName       Si8-NVE
NumberOfAtoms    8
NumberOfSpecies  1
%block ChemicalSpeciesLabel
 1 14 Si
%endblock ChemicalSpeciesLabel
LatticeConstant  5.476 Ang
%block LatticeVectors
 1.0 0.0 0.0
 0.0 1.0 0.0
 0.0 0.0 1.0
%endblock LatticeVectors
AtomicCoordinatesFormat  Fractional
%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
 0.00 0.50 0.50 1
 0.50 0.00 0.50 1
 0.50 0.50 0.00 1
 0.25 0.25 0.25 1
 0.25 0.75 0.75 1
 0.75 0.25 0.75 1
 0.75 0.75 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry
Mesh.Cutoff         150 Ry
%block kgrid.MonkhorstPack
 2 0 0 0.5
 0 2 0 0.5
 0 0 2 0.5
%endblock kgrid.MonkhorstPack
XC.Functional       GGA
XC.Authors          PBE
MD.TypeOfRun        Verlet
MD.InitialTimeStep  1
MD.FinalTimeStep    20
MD.LengthTimeStep   1.0 fs
MD.InitialTemperature  600 K
WriteMDXmol         true
WriteMDHistory      true
EOT

mkdir -p "$P/fe_spin"
cat > "$P/fe_spin/calc.fdf" <<'EOT'
SystemLabel      fe
SystemName       Fe-bcc-FM
NumberOfAtoms    1
NumberOfSpecies  1
%block ChemicalSpeciesLabel
 1 26 Fe
%endblock ChemicalSpeciesLabel
LatticeConstant  2.83 Ang
%block LatticeVectors
 -0.5  0.5  0.5
  0.5 -0.5  0.5
  0.5  0.5 -0.5
%endblock LatticeVectors
AtomicCoordinatesFormat  Fractional
%block AtomicCoordinatesAndAtomicSpecies
  0.0 0.0 0.0 1
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry
Mesh.Cutoff         300 Ry
%block kgrid.MonkhorstPack
 10 0 0 0.0
 0 10 0 0.0
 0 0 10 0.0
%endblock kgrid.MonkhorstPack
XC.Functional       GGA
XC.Authors          PBE
Spin                polarized
OccupationFunction  MP
ElectronicTemperature  0.01 Ry
SCF.Mixer.Weight    0.1
MaxSCFIterations    300
EOT

run() {   # run <folder>: pseudopotentials in, SIESTA, keep what the example uses
    local d="$P/$1"
    for el in $(awk '/%block ChemicalSpeciesLabel/{f=1;next} /%endblock/{f=0} f{print $3}' "$d/calc.fdf"); do
        cp "$BANK/$el.psml" "$d/"
    done
    echo "[siesta] $1"
    (cd "$d" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err) || true
    (cd "$d" && for f in *; do case "$f" in
        calc.fdf|calc.out|MESSAGES|CLOCK|0_NORMAL_EXIT|*.XV|*.STRUCT_OUT|*.EIG|*.KP|*.bands|*.ANI|*.MDE|*.RHO) ;;
        *) rm -rf "$f";; esac; done)
}
for d in bulk_sp relax relax_short scf_fail strain/a_5.40 strain/a_5.476 strain/a_5.55 fe_spin md; do run "$d"; done
