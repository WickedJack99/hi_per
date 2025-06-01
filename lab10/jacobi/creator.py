import math
import subprocess

core_counts = [1, 2, 4, 6, 8, 12, 16, 18, 24, 32, 36, 48]
cores_per_node = 12
max_nodes = 4

for y in core_counts:
    x = math.ceil(y / cores_per_node)
    if x > max_nodes:
        print(f"Skipping {y} cores: requires more than {max_nodes} nodes.")
        continue

    filename = f"slurm_{y}.sh"
    script = f"""#!/bin/bash
#SBATCH --job-name=mpi_g5_{y}
#SBATCH --output=jacobi_{y}.txt
#SBATCH --nodes={x}
#SBATCH --ntasks={y}

srun --mpi=pmix ./jacobi
"""

    with open(filename, "w") as f:
        f.write(script)

    # Submit the job
    result = subprocess.run(["sbatch", filename], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Submitted {filename}: {result.stdout.strip()}")
    else:
        print(f"Failed to submit {filename}: {result.stderr.strip()}")
