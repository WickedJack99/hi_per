import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Configuration
data_dir = "svals"
max_threads = 12
colors = cm.get_cmap("tab10", max_threads)  # or "tab20", "Set3", etc.

# Data: thread_count -> benchmark_iteration -> rows_per_thread list
data = {t: {} for t in range(1, max_threads + 1)}

# Load CSVs
for filename in os.listdir(data_dir):
    if filename.startswith("stats") and filename.endswith(".csv"):
        try:
            name_part = filename.removeprefix("stats").removesuffix(".csv")
            thread_count, iteration = map(int, name_part.split("_"))

            if thread_count <= max_threads:
                path = os.path.join(data_dir, filename)
                with open(path, "r") as f:
                    line = f.read().strip().rstrip(',')
                    values = list(map(int, line.split(",")))
                    data[thread_count][iteration] = values
        except Exception as e:
            print(f"Failed to parse {filename}: {e}")

# Plot
plt.figure(figsize=(10, 6))
linestyles = ['-', '--', '-.', ':']

for thread_count, iterations_data in data.items():
    if not iterations_data:
        continue

    # Get all iterations sorted
    iterations = sorted(iterations_data.keys())

    # Transpose: thread index -> list of values over iterations
    max_threads_in_file = max(len(v) for v in iterations_data.values())
    thread_rows = [[] for _ in range(max_threads_in_file)]

    for iter_num in iterations:
        row_data = iterations_data[iter_num]
        for i, val in enumerate(row_data):
            thread_rows[i].append(val)
    style = linestyles[(thread_count - 1) % len(linestyles)]
    for i, row_counts in enumerate(thread_rows):
        plt.plot(iterations, row_counts, label=f"{thread_count}t-T{i}", color=colors(thread_count - 1), alpha=0.7, linestyle=style)

# Clean up legend: optional to show only one label per thread count
from collections import defaultdict
shown = defaultdict(bool)
handles, labels = plt.gca().get_legend_handles_labels()
new_handles, new_labels = [], []

for h, l in zip(handles, labels):
    key = l.split("t-")[0]
    if not shown[key]:
        new_handles.append(h)
        new_labels.append(f"{key} threads")
        shown[key] = True

plt.legend(new_handles, new_labels, title="Thread Count", loc="upper right")

plt.xlabel("Benchmark Iteration")
plt.ylabel("Rows Computed")
plt.title("Per-Thread Row Distribution by Thread Count")
plt.grid(True)
plt.tight_layout()
plt.show()
