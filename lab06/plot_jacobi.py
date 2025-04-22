import json
import matplotlib.pyplot as plt
import pandas

# Path to the benchmark results JSON
file_static = "benchmark.csv"
file_dyn = "benchmark_dyn.csv"
file_gui = "benchmark_gui.csv"

data_static = pandas.read_csv(file_static)
data_dyn = pandas.read_csv(file_dyn)
data_gui = pandas.read_csv(file_gui)

#Plotting
plt.figure(figsize=(8, 5))
plt.plot(data_static["threads"], data_static["time"], marker="o", color="blue", label="static")
plt.plot(data_dyn["threads"], data_dyn["time"], marker="o", color="red", label="dynamic")
plt.plot(data_gui["threads"], data_gui["time"], marker="o", color="green", label="guided")
plt.xlabel("Thread Count")
plt.ylabel("Execution Time (ms)")
plt.title("Performance of Jacobi Benchmark")
plt.grid(True)
plt.tight_layout()
plt.show()
