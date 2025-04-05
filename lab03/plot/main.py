import json
import matplotlib.pyplot as plt

def parse_benchmark_data(json_file):
    """
    Parses a JSON file containing benchmark data and extracts Size and real_time.

    Args:
        json_file (str): Path to the JSON file.

    Returns:
        tuple: A tuple containing two lists: sizes and real_times.
               Returns empty lists if the file is not found or data is malformed.
    """
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            benchmarks = data.get('benchmarks', [])
            sizes = []
            real_times = []
            flops = []
            for benchmark in benchmarks:
                iterations = benchmark.get('Iterations')
                size = benchmark.get('Size')
                real_time = benchmark.get('real_time')
                if size is not None and real_time is not None:
                    flops.append(7 * iterations / real_time)
                    sizes.append(size)
                    real_times.append(real_time)
            return sizes, real_times, flops
    except FileNotFoundError:
        print(f"Error: File not found: {json_file}")
        return [], []
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from: {json_file}")
        return [], []
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return [], []

def create_plot(sizes, real_times, flops):
    """
    Creates a scatter plot of Size vs. real_time.

    Args:
        sizes (list): List of 'Size' values.
        real_times (list): List of 'real_time' values.
    """
    if not sizes or not real_times:
        print("No data to plot.")
        return

    plt.figure(figsize=(10, 6))
    plt.scatter(sizes, real_times, marker='o', color='blue', label='Real Time (ms)')
    plt.scatter(sizes, flops, marker='x', color='red', label='FLOPs')
    plt.xlabel("Size")
    plt.ylabel("Values")
    plt.title("Benchmark: Real Time and FLOPs vs. Size")
    plt.grid(True)
    plt.legend()  # Show the legend to distinguish the plots
    plt.twinx()  # Create a second y-axis to avoid overlapping labels if scales differ significantly
    plt.ylabel("FLOPs", color='red')
    plt.tick_params(axis='y', labelcolor='red')
    plt.gca().yaxis.grid(False) # Turn off grid for the second y-axis

    plt.show()

if __name__ == "__main__":
    json_file = 'results.json'
    sizes, real_times, flops = parse_benchmark_data(json_file)

    if sizes and real_times:
        create_plot(sizes, real_times, flops)
