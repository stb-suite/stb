#!/bin/bash
# Runs the SIESTA calculations of the 2.19 example (runs/), each with Charge.Hirshfeld, Charge.Voronoi and Mulliken
# populations in the log:
#   runs/h2o_basis/{SZ,DZ,DZP,TZP}/   one water molecule (experimental geometry) in a 10 Ang box, four basis sizes
#   runs/nacl/                        rock-salt NaCl, 2-atom cell, a = 5.64 Ang
#   runs/o2/                          the O2 triplet, spin-polarized
# SIESTA 5.x + mpirun, about a minute on 4 processes. Pseudopotentials from the suite's bundled PseudoDojo bank
# (copied in, then removed: not kept in the repository). Only the input, the log and the .XV are kept.
#     bash run_example.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
R="$HERE/runs"
NP=${1:-4}
BANK="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"
CHARGES='Charge.Hirshfeld    end           # SIESTA'"'"'s own Hirshfeld charges, at the end of the run
Charge.Voronoi      end           # ... and Voronoi
WriteMullikenPop    1             # ... and Mulliken (atomic populations)'

run() {   # run <folder>
    local d="$R/$1"
    for el in $(awk '/%block ChemicalSpeciesLabel/{f=1;next} /%endblock/{f=0} f{print $3}' "$d/calc.fdf"); do
        cp "$BANK/$el.psml" "$d/"
    done
    echo "[siesta] $1"
    (cd "$d" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err)
    (cd "$d" && for f in *; do case "$f" in calc.fdf|calc.out|*.XV) ;; *) rm -rf "$f";; esac; done)
}

for b in SZ DZ DZP TZP; do
    mkdir -p "$R/h2o_basis/$b"
    cat > "$R/h2o_basis/$b/calc.fdf" <<EOT
# Water, experimental geometry (O-H 0.9572 Ang, H-O-H 104.52 deg), 10 Ang box, basis $b
SystemLabel      h2o
SystemName       H2O-$b
NumberOfAtoms    3
NumberOfSpecies  2
%block ChemicalSpeciesLabel
 1 8 O
 2 1 H
%endblock ChemicalSpeciesLabel
LatticeConstant  1.0 Ang
%block LatticeVectors
 10.0  0.0  0.0
  0.0 10.0  0.0
  0.0  0.0 10.0
%endblock LatticeVectors
AtomicCoordinatesFormat  Ang
%block AtomicCoordinatesAndAtomicSpecies
 5.000000  5.000000  5.000000  1
 5.000000  5.756950  5.585882  2
 5.000000  4.243050  5.585882  2
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize       $b
PAO.EnergyShift     0.01 Ry
Mesh.Cutoff         250 Ry
XC.Functional       GGA
XC.Authors          PBE
SCF.DM.Tolerance    1.0d-5
$CHARGES
EOT
    run "h2o_basis/$b"
done

mkdir -p "$R/nacl"
cat > "$R/nacl/calc.fdf" <<EOT
# Rock-salt NaCl, primitive fcc cell, a = 5.64 Ang (experiment)
SystemLabel      nacl
SystemName       NaCl
NumberOfAtoms    2
NumberOfSpecies  2
%block ChemicalSpeciesLabel
 1 11 Na
 2 17 Cl
%endblock ChemicalSpeciesLabel
LatticeConstant  5.64 Ang
%block LatticeVectors
 0.0 0.5 0.5
 0.5 0.0 0.5
 0.5 0.5 0.0
%endblock LatticeVectors
AtomicCoordinatesFormat  Fractional
%block AtomicCoordinatesAndAtomicSpecies
 0.0 0.0 0.0 1
 0.5 0.5 0.5 2
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize       DZP
PAO.EnergyShift     0.01 Ry
Mesh.Cutoff         300 Ry
%block kgrid.MonkhorstPack
 6 0 0 0.5
 0 6 0 0.5
 0 0 6 0.5
%endblock kgrid.MonkhorstPack
XC.Functional       GGA
XC.Authors          PBE
SCF.DM.Tolerance    1.0d-5
$CHARGES
EOT
run nacl

mkdir -p "$R/o2"
cat > "$R/o2/calc.fdf" <<EOT
# O2 triplet in an 8 Ang box, spin-polarized
SystemLabel      o2
SystemName       O2-triplet
NumberOfAtoms    2
NumberOfSpecies  1
%block ChemicalSpeciesLabel
 1 8 O
%endblock ChemicalSpeciesLabel
LatticeConstant  1.0 Ang
%block LatticeVectors
 8.0 0.0 0.0
 0.0 8.0 0.0
 0.0 0.0 8.0
%endblock LatticeVectors
AtomicCoordinatesFormat  Ang
%block AtomicCoordinatesAndAtomicSpecies
 4.0000  4.0000  3.3950  1
 4.0000  4.0000  4.6050  1
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize       DZP
PAO.EnergyShift     0.01 Ry
Mesh.Cutoff         200 Ry
XC.Functional       GGA
XC.Authors          PBE
Spin                polarized
SCF.DM.Tolerance    1.0d-5
SCF.Mixer.Weight    0.1
MaxSCFIterations    300
$CHARGES
EOT
run o2
