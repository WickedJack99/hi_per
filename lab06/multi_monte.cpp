#include "matrix.h"
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>

struct bench_result {
  double time;
  double pi;
};

double monteCarloPi(int64_t N) {
  int localCount = 0;
  std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
#pragma omp parallel for reduction(+ : localCount)
  for (int i = 0; i < N; ++i) {
    double x = dist(generator);
    double y = dist(generator);
    if (x * x + y * y <= 1.0) {
      localCount++;
    }
  }

  return 4.0 * localCount / N;
}

bench_result benchmark(int numThreads, int64_t N, double function(int64_t)) {
  omp_set_num_threads(numThreads);

  // start time measurement
  auto start = std::chrono::system_clock::now();
  auto result = bench_result{0.0, 0.0};

  // perform benchmark
  const int maxNumExec = 10;
  for (int nExec = 0; nExec < maxNumExec; ++nExec) {
    result.pi = function(N);
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double, std::milli> diff = end - start;
  result.time = diff.count() / maxNumExec;

  return result;
}

int main() {
  const int64_t max_n = 10000000;
  const int numThreads = 12;
  std::ofstream fout("multi_monte.csv", std::ios::app);
  fout << "pi" << ",";
  fout << "time" << "\n";

  for (int64_t i = 1000; i <= max_n; i += 1000) {
    const bench_result res = benchmark(numThreads, i, monteCarloPi);

    fout << res.pi << ",";
    fout << res.time << "\n";
  }
}
