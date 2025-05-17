#!/bin/bash                                                                                                             
#SBATCH --job-name=excercise2_vectorTriad                                                                                        
#SBATCH --output=excercise2_vectorTriad.out                                                                             
#SBATCH --ntasks=1                                                                                                      
#SBATCH --cpus-per-task=6                                                                                               
export OMP_PLACES=cores                                                                                                 
export OMP_PROC_BIND=close                                                                                                                                                                                                                      
# Load modules or activate environment if needed                                                                        
# module load gcc cmake ...                                                                                                                                                                                                                     
# Run the benchmark binary                                                                                              
g++ vectorTriad.cpp -o vector_triad -fopenmp -O3 -lbenchmark                                                              
./vector_triad --benchmark_out=excercise2_vectorTriad.json --benchmark_out_format=json