# Makefile for DBSCAN program

# ----------------------------------------------------
# Parameters
# Change these parameters according to your needs.

# SOURCE_FILES: The source files of the algorithm, used for each build.
# You can add more source files here if needed.
SOURCE_FILES = dbscan.cpp point.cpp

# Main rogram, used to cluster the data and save the result.
# PROGRAM_NAME: The name of the program that will be generated after compilation.
PROGRAM_NAME = dbscan
RUN_MAIN = run.cpp

# Benchmark program: This program is used to benchmark the performance of the algorithm.
# It is not used for the actual clustering process.
BENCHMARK_PROGRAM_NAME = dbscan_bench
BENCHMARK_MAIN = benchmark.cpp

COMPILER_FLAGS = -fopenmp -std=c++17 -lpthread

# ----------------------------------------------------
# The actual makefile rules, only change these if you really need to.

# Default target
# The default target is the one that will be executed when you run 'make' without any arguments.
default: release

release: $(RUN_MAIN) $(SOURCE_FILES)
	g++ $(RUN_MAIN) $(SOURCE_FILES) $(COMPILER_FLAGS) -o $(PROGRAM_NAME) -O3

debug: $(RUN_MAIN) $(SOURCE_FILES)
	g++ $(RUN_MAIN) $(SOURCE_FILES) $(COMPILER_FLAGS) -o $(PROGRAM_NAME) -O0 -g

benchmark: $(BENCHMARK_MAIN) $(SOURCE_FILES)
	g++ $(BENCHMARK_MAIN) $(SOURCE_FILES) $(COMPILER_FLAGS) -o $(BENCHMARK_PROGRAM_NAME) -O3 -lbenchmark

run_bench: benchmark
	./$(BENCHMARK_PROGRAM_NAME) 

run: release
	./$(PROGRAM_NAME)
