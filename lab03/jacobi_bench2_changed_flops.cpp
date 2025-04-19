#include "matrix.h"
#include <iostream>
#include <fstream>
#include <benchmark/benchmark.h>

Matrix jacobi(const Matrix& init, double eps, int maxNumIter) {
  Matrix phi(init);
  Matrix tmp(init.dim());

  const int n = phi.dim1();
  const int m = phi.dim2();

  for (int i = 0; i < n; ++i) {
    tmp(i, 0) = phi(i, 0);
    tmp(i, m - 1) = phi(i, m - 1);
  }
  for (int j = 0; j < m; ++j) {
    tmp(0, j) = phi(0, j);
    tmp(n - 1, j) = phi(n - 1, j);
  }

  double dist = 1e10;
  int nIter = 0;
  while (dist > eps && nIter < maxNumIter) {
    dist = 0;
    for (int i = 1; i < n - 1; ++i) {
      for (int j = 1; j < m - 1; ++j) {
        tmp(i, j) = .25 * (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) +
                                 phi(i, j - 1));
        dist = std::max(dist, std::abs(tmp(i, j) - phi(i, j)));
      }
    }

    std::swap(tmp, phi);

    nIter++;
  }

  return phi;
}

Matrix initialCondition(int n) {
  Matrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const double x = 1. * i / n;

    for (int j = 0; j < n; ++j) {
      const double y = 1. * j / n;

      // The lower left corner is hot
      m(i, j) += exp(-0.5 * (x * x + y * y) * 10) * exp(-100 * x * y);

      // The right side is cold
      m(i, j) += -1.0 * exp(-20 * (1 - y));
    }
  }
  return m;
}

void storeMatrix(const Matrix& m, const std::string& filename) {
  std::ofstream fout(filename);
  fout << m << std::endl;
}

// int main() {
//   const Matrix init = initialCondition(200);
//   storeMatrix(init, "init.asc");
//   Matrix phi = jacobi(init, 1e-5, 100000);
//   storeMatrix(phi, "solution.asc");
// }

void benchmarkJacobi(benchmark::State& state) {
  int n = state.range(0);
  const Matrix init = initialCondition(n);
  Matrix phi;

  for (auto _ : state) {
    phi = jacobi(init, 1e-5, 100000);  // Reduce max iterations for quicker tests
  }

  state.counters["Size"] = n; // Optional: include n in the output
}

BENCHMARK(benchmarkJacobi)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->DenseRange(100, 2000, 20);
BENCHMARK_MAIN();

