import json
import matplotlib.pyplot as plt

# Path to the benchmark results JSON
results = ""

# Load the data
def load_data(filename):
    with open(filename, "r") as f:
        data = json.load(f)

    # Constants
    total_flops = 18000

    # Data storage
    thread_counts = []
    flops = []

    # Extract relevant data
    for entry in data["benchmarks"]:
        thread_count = int(entry["name"].split("/")[-1])
        time_per_iter_ms = entry["real_time"]
        time_total_s = (time_per_iter_ms) / 1000  # convert ms to s

        flops_per_second = total_flops / time_total_s

        thread_counts.append(thread_count)
        flops.append(flops_per_second)

    return thread_counts, flops


threads_VPRC, flops_VPRC = load_data(results)

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(threads_VPRC, flops_VPRC, marker="o", color="blue", label="VPRC")
plt.xlabel("Thread Count")
plt.ylabel("Performance (FLOPS)")
plt.title("Performance (FLOPS) vs Thread Count")
plt.grid(True)
plt.tight_layout()
plt.show()
