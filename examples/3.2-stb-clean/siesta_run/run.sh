#!/bin/bash
#PBS -N si_run
#PBS -l nodes=1:ppn=4
#PBS -o job.o
#PBS -e job.err
cd $PBS_O_WORKDIR
mpirun -np 4 siesta calc.fdf > calc.out
