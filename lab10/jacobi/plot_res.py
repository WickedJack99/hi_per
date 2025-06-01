import matplotlib.pyplot as plt

# Pfad zum Verzeichnis mit den Dateien
data_dir = "resluts/"  # ggf. anpassen

threads = [1, 2, 4, 6, 8, 12, 16, 18, 24, 32, 36, 48]
times = [
    175246,
    89030.8,
    45182.7,
    30198.1,
    22710.1,
    15256.1,
    13512.1,
    13073.2,
    10036.8,
    9845.38,
    8423.65,
    7676.71,
]
speedups = [
    1,
    1.96838,
    3.90997,
    5.8091,
    7.73292,
    11.4934,
    13.0374,
    13.4461,
    17.4421,
    17.9058,
    20.8097,
    22.8728,
]

efficiency = [speedup / thread for speedup, thread in zip(speedups, threads)]

# plt.figure(figsize=(10, 5))
plt.plot(threads, efficiency, marker="o", label="Speedup")
plt.title("Efficiency")
plt.xlabel("Threads")
plt.ylabel("Efficiency")
plt.grid(True)
# plt.legend()
# plt.tight_layout()
plt.savefig("efficiency_plot.png")

plt.clf()
# --- Plot 2: Performance ---
plt.plot(threads, times, marker="o")
plt.title("Performance")
plt.xlabel("Threads")
plt.ylabel("Performance (ms)")
plt.grid(True)
plt.savefig("performance_plot.png")
