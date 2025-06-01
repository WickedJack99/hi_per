#!/bin/bash
git pull
make
sbatch slurm.sh
squeue
