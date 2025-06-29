#!/bin/bash
#SBATCH --job-name=jacobi
#SBATCH --output=jacobi.out
#SBATCH --error=jacobi.err
#SBATCH --ntasks=8                     # 8 MPI processes (1 per socket)
#SBATCH --cpus-per-task=6              # 6 threads per MPI process
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=6
export OMP_PROC_BIND=true
export OMP_PLACES=cores

srun --mpi=pmix ./jacobi
