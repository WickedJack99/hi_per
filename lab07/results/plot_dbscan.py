import csv
import matplotlib.pyplot as plt

# Read performance data from CSV file
with open('lab07\\results\\dbscna_results.csv', 'r') as f:
    reader = csv.reader(f)
    data = next(reader)  # Read first line
    times = list(map(int, data))  # Convert to integers

# X: thread count (1 to n), Y: performance (1/time)
threads = list(range(1, len(times) + 1))
performance = [1 / t for t in times]  # You could multiply by a constant to scale if needed

speedup = [times[0] / t for t in times]

efficiency = []
for i in range (1, len(speedup) + 1):
    efficiency.append(speedup[i-1] / i)

# Plot
plt.plot(threads, efficiency, marker='o')
plt.xlabel('Thread Count')
plt.ylabel('Efficiency (Speedup / Thread Count)')
plt.title('Thread Count vs Efficiency')
plt.grid(True)
plt.show()
