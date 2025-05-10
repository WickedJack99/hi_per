#!/bin/bash
#SBATCH --job-name=kaiJacobiWave
#SBATCH --output=kaiJacobiWaveSplitWork.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

# Load modules or activate environment if needed
# module load gcc cmake ...

# Run the benchmark binary
./kaiJacobiWave --benchmark_out=kai_results.json --benchmark_out_format=json
