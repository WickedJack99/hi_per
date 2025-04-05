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
            for benchmark in benchmarks:
                size = benchmark.get('Size')
                real_time = benchmark.get('real_time')
                if size is not None and real_time is not None:
                    sizes.append(size)
                    real_times.append(real_time)
            return sizes, real_times
    except FileNotFoundError:
        print(f"Error: File not found: {json_file}")
        return [], []
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from: {json_file}")
        return [], []
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return [], []

def create_plot(sizes, real_times):
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
    plt.scatter(sizes, real_times, marker='o')
    plt.xlabel("Size")
    plt.ylabel("Real Time (ms)")
    plt.title("Benchmark: Real Time vs. Size")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    json_file = 'result.json'
    sizes, real_times = parse_benchmark_data(json_file)

    if sizes and real_times:
        create_plot(sizes, real_times)
