#!/bin/bash
# Builds and runs the SIESTA calculations of the 3.7 example (runs/), kept exactly as SIESTA leaves them:
#   runs/si_relax/      bulk Si: variable-cell CG relaxation + bands + DOS/PDOS + density + H/S, with a job script,
#                       an %include'd structure, a BandLines block read from bandlines.txt, and analysis/ (a plot)
#   runs/strain/a_X/    three single points at a = 5.40, 5.476, 5.55 Ang (a series to archive one folder each)
# SIESTA 5.x + mpirun, about 30 s on 4 processes. Pseudopotential from the suite's bundled bank (an input: it belongs
# in the archive, but it is not kept in the repository -- removed at the end, copied back in by the example scripts).
#     bash run_example.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
R="$HERE/runs"
NP=${1:-4}
PSML="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')/Si.psml"
rm -rf "$R"; mkdir -p "$R/si_relax"
cd "$R/si_relax"
cat > structure.fdf <<'EOT'
NumberOfAtoms    2
NumberOfSpecies  1
%block ChemicalSpeciesLabel
 1 14 Si
%endblock ChemicalSpeciesLabel
LatticeConstant  5.55 Ang
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
EOT
cat > bandlines.txt <<'EOT'
  1  0.5 0.5 0.5  L
 40  0.0 0.0 0.0  \Gamma
 45  0.0 0.5 0.5  X
EOT
cat > calc.fdf <<'EOT'
# Bulk Si: variable-cell CG relaxation, then bands, DOS/PDOS, the density and H/S
SystemLabel      si
SystemName       Si-bulk
%include structure.fdf
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
MD.TypeOfRun        CG
MD.VariableCell     true
MD.NumCGsteps       30
MD.MaxForceTol      0.01 eV/Ang
MD.MaxStressTol     0.1 GPa
BandLinesScale      ReciprocalLatticeVectors
%block BandLines < bandlines.txt
%block ProjectedDensityOfStates
 -15.0 5.0 0.05 400 eV
%endblock ProjectedDensityOfStates
SaveRho             true
SaveHS              true
WriteMDXmol         true
EOT
cat > submit.slurm <<'EOT'
#!/bin/bash
#SBATCH --job-name=si_relax
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
mpirun -np 4 siesta calc.fdf > calc.out
EOT
cp "$PSML" .
echo "[siesta] si_relax"
OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err
mkdir -p analysis
python3 - <<'PY'
# a post-processing result next to the run: the band structure from si.bands
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
t = open("si.bands").read().split()
ef = float(t[0]); nb, ns, nk = int(t[5]), int(t[6]), int(t[7])
v = np.array(t[8:8 + nk * (1 + nb * ns)], float).reshape(nk, 1 + nb * ns)
plt.figure(figsize=(4, 4)); plt.plot(v[:, 0], v[:, 1:] - ef, color="#1f6f8b", lw=1)
plt.ylim(-13, 6); plt.ylabel("E - E_F (eV)"); plt.xticks([]); plt.title("Si, PBE/DZP")
plt.tight_layout(); plt.savefig("analysis/si_bands.png", dpi=100)
np.savetxt("analysis/si_bands.dat", np.column_stack([v[:, 0], v[:, 1:] - ef]), fmt="%.5f")
PY

for a in 5.40 5.476 5.55; do
    d="$R/strain/a_$a"; mkdir -p "$d"; cd "$d"
    sed -e "s/^LatticeConstant  5.55 Ang/LatticeConstant  $a Ang/" "$R/si_relax/structure.fdf" > structure.fdf
    sed -e '/^MD\./d; /BandLines/d; /ProjectedDensity/,/endblock ProjectedDensity/d; /^SaveHS/d; /^WriteMDXmol/d' \
        -e 's/^# Bulk Si.*/# Bulk Si single point/' "$R/si_relax/calc.fdf" > calc.fdf
    cp "$PSML" .
    echo "[siesta] strain/a_$a"
    OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err
    rm -f Si.ion Si.ion.xml Si.ion.nc si.DM si.HSX si.ORB_INDX si.RHO ./*.grid.nc   # trimmed for the repository
done
# The pseudopotential is an input and goes into the archives, but it is not kept in the repository: the example
# scripts copy it in from the bundled bank (bank_path("dojo")/Si.psml) before archiving.
rm -f "$R"/si_relax/Si.psml "$R"/strain/*/Si.psml
