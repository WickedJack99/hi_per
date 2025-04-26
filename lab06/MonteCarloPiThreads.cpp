#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>

struct bench_result {
  double time;
  int threads;
};

double monteCarloPi(int64_t N) {
  int localCount = 0;
#pragma omp parallel
  {
    std::mt19937 generator(std::random_device{}() + omp_get_thread_num());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
#pragma omp for reduction(+ : localCount)
    for (int i = 0; i < N; ++i) {
      double x = dist(generator);
      double y = dist(generator);
      if (x * x + y * y <= 1.0) {
        localCount++;
      }
    }
  }

  return 4.0 * localCount / N;
}

bench_result benchmark(int numThreads, int64_t N, double function(int64_t)) {
  omp_set_num_threads(numThreads);

  // start time measurement
  auto start = std::chrono::system_clock::now();
  auto result = bench_result{0.0, 0};

  // perform benchmark
  const int maxNumExec = 10;
  for (int nExec = 0; nExec < maxNumExec; ++nExec) {
    volatile double pi = function(N);
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double, std::milli> diff = end - start;
  result.time = diff.count() / maxNumExec;
  result.threads = omp_get_thread_num();

  return result;
}

int main() {
  const int64_t max_n = 10000000;
  const int max_threads = 12;
  std::ofstream fout("diff_theads.csv", std::ios::app);
  fout << "threads" << ",";
  fout << "max_n" << ",";
  fout << "time" << "\n";

  for (int i = 1; i <= max_threads; i++) {
    const bench_result res = benchmark(i, max_n, monteCarloPi);

    fout << res.threads << ",";
    fout << max_n << ",";
    fout << res.time << "\n";
  }
}
