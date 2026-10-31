#!/bin/bash
#SBATCH --job-name=si_run
#SBATCH --ntasks=4
#SBATCH --output=slurm-%j.out
srun siesta calc.fdf > calc.out
