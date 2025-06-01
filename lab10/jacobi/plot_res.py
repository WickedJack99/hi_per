import matplotlib.pyplot as plt
import os
import re

# Pfad zum Verzeichnis mit den Dateien
data_dir = "resluts/"  # ggf. anpassen

threads = [1, 2,4,6,8,12,16,18,24,32,36,48]
times = []
speedups = []


plt.figure(figsize=(10, 5))
plt.plot(threads, speedups, marker='o', label="Speedup")
# plt.plot(threads, threads, linestyle='--', label="Ideal (Linear Speedup)")
plt.title("Speedup")
plt.xlabel("Threads")
plt.ylabel("Speedup")
plt.grid(True)
plt.legend()
# plt.tight_layout()
plt.savefig("speedup_plot.png")

# --- Plot 2: Performance ---
plt.figure(figsize=(10, 5))
plt.plot(threads, times, marker='o')
plt.title("Performance")
plt.xlabel("Threads")
plt.ylabel("Performance (ms)")
plt.grid(True)
# plt.tight_layout()
plt.savefig("performance_plot.png")

print("Plots gespeichert: speedup_plot.png und performance_plot.png")

