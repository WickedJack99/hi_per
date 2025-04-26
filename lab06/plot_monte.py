import json
import math
import matplotlib.pyplot as plt
import pandas

# Path to the benchmark results JSON
plot_single = "./results/really_single_monte.csv"
plot_multi = "./results/multi_monte.csv"


data_single = pandas.read_csv(plot_single)
data_multi = pandas.read_csv(plot_multi)


def calculate_error(pi: float) -> float:
    return abs(math.pi - pi)


data_single["error"] = data_single.apply(lambda x: calculate_error(x["pi"]), axis=1)
data_multi["error"] = data_multi.apply(lambda x: calculate_error(x["pi"]), axis=1)

plt.figure(figsize=(8, 5))
plt.plot(
    data_single["error"].index * 1000,
    data_single["error"],
    color="blue",
    label="single",
    linewidth=0.1,
)
plt.plot(
    data_multi["error"].index * 5000,
    data_multi["error"],
    color="red",
    label="multi",
    linewidth=0.1,
)
plt.xlabel("time")
plt.ylabel("error")
plt.title("pin Monte Carlo Benchmark")
plt.grid(True)
plt.tight_layout()
plt.savefig("mein_plot.png")
