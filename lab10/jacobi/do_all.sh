#!/bin/bash
git pull
rm jacobi.txt
make
sbatch slurm.sh
squeue
