#include <iostream>
#include <fstream>
#include <chrono>
#include <benchmark/benchmark.h>
#include "dbscan.h"

using namespace HPC;

static void BM_DBSCAN(benchmark::State& state) {
  // Load points from file
  std::vector<Point> points = readPointsFromFile("data");

  // Create DBSCAN object with parameters from the benchmark state
  DBSCAN ds(5, 0.01);

  // Measure the time taken to run DBSCAN
  for (auto _ : state) {
    ds.run(points);
  }
}

BENCHMARK(BM_DBSCAN)->Unit(benchmark::kMillisecond)->Iterations(10);
BENCHMARK_MAIN();