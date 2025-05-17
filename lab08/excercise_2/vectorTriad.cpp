#include <benchmark/benchmark.h>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

static void vectorTriad(benchmark::State& state) {
  const int n = 1e9;
#pragma omp parallel num_threads(6)
  {
    std::cout << "CPU: " << sched_getcpu() << " Thread: " << omp_get_thread_num() << std::endl;
    const int numThreads = omp_get_num_threads();
    const int nPerThread = n / numThreads;
    std::vector<double> A(nPerThread);
    std::vector<double> B(nPerThread, 1.);
    std::vector<double> C(nPerThread, 2.);
    std::vector<double> D(nPerThread, 3.);
    for (auto _ : state) {
      for (int i = 0; i < nPerThread; ++i) {
        A[i] = B[i] + C[i] * D[i];
      }
      volatile double dummy = A[0];
    }
  }
}

BENCHMARK(vectorTriad)->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();