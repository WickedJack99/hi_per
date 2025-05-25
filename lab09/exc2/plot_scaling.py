import glob
import re
import matplotlib.pyplot as plt

# Find all result files like results1.txt, results2.txt, etc.
files = sorted(glob.glob("results*.txt"), key=lambda f: int(re.search(r"\d+", f).group()))

cores = [1,2,4,6,8,12,16,18,24,32,36,48]
core_increase = [core / cores[i - 1] for i, core in enumerate(cores[1:], start=1)]
core_increase.insert(0, 1)

runtimes = {}
for file in files:
    with open(file, "r") as f:
        line = f.readline()
        match = re.search(r"Runtime:\s+([\d.]+)\s+ms, Tasks:\s+(\d+)", line)
        if match:
            runtime = float(match.group(1))
            tasks = int(match.group(2))
            runtimes[tasks] = runtime
        else:
            print(f"Skipping unrecognized format in {file}")

# Sort by number of tasks
tasks_sorted = sorted(runtimes.keys())
T1 = runtimes[1]

speedups = [T1 / runtimes[t] for t in tasks_sorted]
efficiencies = [s / t for s, t in zip(speedups, cores)]

print(core_increase)
print(speedups)

# Plot speedup
plt.figure()
plt.plot(tasks_sorted, speedups, marker='o')
plt.title("Speedup vs. Number of Tasks")
plt.xlabel("Tasks")
plt.ylabel("Speedup")
plt.grid(True)

# Plot efficiency
plt.figure()
plt.plot(tasks_sorted, efficiencies, marker='o')
plt.title("Efficiency vs. Number of Tasks")
plt.xlabel("Tasks")
plt.ylabel("Efficiency")
plt.grid(True)

plt.show()
