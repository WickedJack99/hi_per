#!/bin/bash
#SBATCH --job-name=jacobiWaveSplitWorkTest
#SBATCH --output=benchmark_jacobiWaveSplitWork.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

# Load modules or activate environment if needed
# module load gcc cmake ...

# Run the benchmark binary
./jacobiWaveSplitWork --benchmark_out=results.json --benchmark_out_format=json
