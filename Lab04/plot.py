import pandas as pd
import matplotlib.pyplot as plt

count_threads = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
average_flops = []
average_runtimes = []
for i in range(0, 12):
    avg_flops = 0
    avg_runtime = 0
    with open(f"./stats/stats{i}.csv", 'r') as file:
        line = file.read()
        stats = line.split(";")
        for stat in stats:
            splited = stat.split(",")
            if len(splited) == 2:
                avg_flops += float(splited[0])
                avg_runtime += float(splited[1])
    average_flops.append(avg_flops / 1000)
    if i > 0:
        average_runtimes.append(average_runtimes[0]/avg_runtime)
    else:
        average_runtimes.append(1)


efficiencies = []

for j in range(0, 12):
    efficiencies.append(average_runtimes[j] / count_threads[j] * 100) 
    
# First plot: FLOPs vs Threads
plt.figure(figsize=(10, 4))
plt.plot(count_threads, average_flops, marker='o')
plt.title("Average FLOPs vs Number of Threads")
plt.xlabel("Number of Threads")
plt.ylabel("Average FLOPs")
plt.grid(True)
plt.tight_layout()
plt.show()

# Second plot: Runtime vs Threads
plt.figure(figsize=(10, 4))
plt.plot(count_threads, average_runtimes, marker='o', color='orange')
plt.title("Speedup")
plt.xlabel("Number of Threads")
plt.ylabel("Factor in %")
plt.grid(True)
plt.tight_layout()
plt.show()

# Second plot: Runtime vs Threads
plt.figure(figsize=(10, 4))
plt.plot(count_threads, efficiencies, marker='o', color='red')
plt.title("Efficiency")
plt.xlabel("Number of Threads")
plt.ylabel("Speedup per Number of Threads (/count)")
plt.grid(True)
plt.tight_layout()
plt.show()  