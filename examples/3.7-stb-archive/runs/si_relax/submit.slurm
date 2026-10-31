#!/bin/bash
#SBATCH --job-name=si_relax
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
mpirun -np 4 siesta calc.fdf > calc.out
