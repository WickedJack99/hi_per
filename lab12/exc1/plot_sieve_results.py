import os
import re
import matplotlib.pyplot as plt

def parse_result_files(fixed_N, directory="."):
    pattern = re.compile(r"runtime_p_(\d+)_N_({})\.txt".format(fixed_N))
    results = []

    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            p = int(match.group(1))
            with open(os.path.join(directory, filename), "r") as f:
                line = f.readline().strip()
                time_str, N_str = line.split(",")
                time = float(time_str)
                results.append((p, time))

    results.sort(key=lambda x: x[0])  # sort by number of processes
    return results

def plot_results(results, N, output_file="sieve_scaling_plot.png"):
    ps = [p for p, _ in results]
    times = [t for _, t in results]
    
    performance = [N / t for t in times]
    speedup = [times[0] / t for t in times]
    efficiency = [s / p for s, p in zip(speedup, ps)]

    plt.figure(figsize=(12, 8))

    # Performance plot
    plt.subplot(3, 1, 1)
    plt.plot(ps, performance, marker='o')
    plt.title(f"Performance (N = {N})")
    plt.xlabel("Processes")
    plt.ylabel("Performance (N / time)")

    # Speedup plot
    plt.subplot(3, 1, 2)
    plt.plot(ps, speedup, marker='o')
    plt.title("Speedup")
    plt.xlabel("Processes")
    plt.ylabel("Speedup (T1 / Tp)")

    # Efficiency plot
    plt.subplot(3, 1, 3)
    plt.plot(ps, efficiency, marker='o')
    plt.title("Efficiency")
    plt.xlabel("Processes")
    plt.ylabel("Efficiency (Speedup / p)")

    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved to: {output_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int, help="Fixed N to plot results for")
    parser.add_argument("--dir", type=str, default=".", help="Directory with result files")
    parser.add_argument("--out", type=str, default="sieve_scaling_plot.png", help="Output image file")
    args = parser.parse_args()

    data = parse_result_files(args.N, args.dir)
    if not data:
        print(f"No matching files found for N = {args.N}")
    else:
        plot_results(data, args.N, args.out)
