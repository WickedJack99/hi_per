#!/bin/bash
#SBATCH --job-name=G5_jacobi         # Name des Jobs
#SBATCH --output=jacobi.txt      # Output-Datei (%j = Job-ID)
#SBATCH --nodes=4              # Nutze z. B. 4 Nodes
#SBATCH --ntasks-per-node=12

# Programm starten
mpirun ./jacobi
