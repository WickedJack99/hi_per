#!/bin/bash
#SBATCH --job-name=G5_jacobi         # Name des Jobs
#SBATCH --output=jacobi.txt      # Output-Datei (%j = Job-ID)
#SBATCH --ntasks=48                  # Anzahl Tasks (z. B. MPI-Prozesse)
#SBATCH --cpus-per-task=1           # Anzahl CPU-Kerne pro Task

# Programm starten
./jacobi
