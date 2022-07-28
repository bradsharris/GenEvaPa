#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --requeue
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high
#SBATCH --exclusive
##SBATCH --exclude=
#SBATCH --job-name=Evap
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=emailhere

module load fftw
module load gromacs
module load bio3

export OMP_NUM_THREADS=1

./AutoEvap.py -tol 10 -nd 385
