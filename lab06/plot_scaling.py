import json
import math
import matplotlib.pyplot as plt
import pandas

# Path to the benchmark results JSON
plot_diff = "./results/diff_theads.csv"
plot_dyn = "./results/dyn_n.csv"


data_diff = pandas.read_csv(plot_diff)
data_dyn = pandas.read_csv(plot_dyn)

plt.figure(figsize=(8, 5))
plt.plot(
    data_diff["threads"],
    data_diff["time"],
    color="blue",
)
plt.xlabel("threads")
plt.ylabel("time")
plt.grid(True)
plt.tight_layout()
plt.savefig("diff.png")

plt.figure(figsize=(8, 5))
plt.plot(
    data_dyn["max_n"],
    data_dyn["time"],
    color="red",
)
plt.xlabel("n")
plt.ylabel("time")
plt.grid(True)
plt.tight_layout()
plt.savefig("dyn.png")
