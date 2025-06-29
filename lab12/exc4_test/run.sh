#!/bin/bash
#SBATCH --job-name=jacobi
#SBATCH --output=jacobi.out
#SBATCH --error=jacobi.err
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

srun --mpi=pmix ./jacobi