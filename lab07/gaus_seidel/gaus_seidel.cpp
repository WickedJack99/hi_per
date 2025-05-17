#include "matrix.h"
#include <benchmark/benchmark.h>
#include <iostream>
#include <omp.h>

void gauss_seidel(Matrix &phi, int maxNumIter) {
  const int m = phi.dim1();
  const int n = phi.dim2();

  const double osth = 1. / 4;

  for (int iter = 0; iter < maxNumIter; ++iter) {
    for (int i = 1; i < m - 1; ++i) {
      for (int j = 1; j < n - 1; ++j) {
        phi(i, j) = osth * (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) +
                            phi(i, j - 1));
      }
    }
  }
}

void gauss_seidel_par(Matrix &phi, int maxNumIter) {
  const int m = phi.dim1();
  const int n = phi.dim2();

#pragma omp parallel num_threads(12)
  {
    const double osth = 1. / 4;
    const int num_theads = omp_get_num_threads();
    const int thread_num = omp_get_thread_num();
    const int chunk = (m - 2) / num_theads;
    const int start = 1 + chunk * thread_num;
    const int end = chunk * (thread_num + 1);

    for (int iter = 0; iter < maxNumIter; ++iter) {
      // printf("thread %d, start: %d, end: %d\n", omp_get_thread_num(), start,
      //        end);

      for (int j = 1; j < n + num_theads - 1; ++j) {
        int k = j - thread_num;

        if (k > 0 && k < n - 1) {
          for (int i = start; i <= end; ++i) {
            // printf("thread %d, i: %d, j: %d, k: %d\n", omp_get_thread_num(),
            // i, j,
            //        k);

            phi(i, k) = osth * (phi(i + 1, k) + phi(i - 1, k) + phi(i, k + 1) +
                                phi(i, k - 1));
          }
        }
#pragma omp barrier
      }
    }
  }
}

void fill_matrix(Matrix &matrix, const int filler) {
  const int m = matrix.dim1();
  const int n = matrix.dim2();

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == 0 || i == m - 1 || j == 0 || j == n - 1) {
        matrix(i, j) = filler;
      } else
        matrix(i, j) = 0;
    }
  }
}

void check(const Matrix &a, const Matrix &b) {
  const int m = a.dim1();
  const int n = a.dim2();

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (a(i, j) != b(i, j)) {
        std::cout << "Not equal at (" << i << ", " << j << "), a: " << a(i, j)
                  << " != " << b(i, j) << " :b" << std::endl;
      }
    }
  }
}

void print_matrix(const Matrix &matrix) {
  const int m = matrix.dim1();
  const int n = matrix.dim2();

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      std::cout << matrix(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

void benchmarkParallel(benchmark::State &state) {
  int iterations = state.range(0);
  Matrix matrix = Matrix(30002, 20002);

  fill_matrix(matrix, 10);

  for (auto _ : state) {
    gauss_seidel_par(matrix, iterations);
    benchmark::DoNotOptimize(matrix);
  }
}

void benchmarkSerial(benchmark::State &state) {
  int iterations = state.range(0);
  Matrix matrix = Matrix(30002, 20002);

  fill_matrix(matrix, 10);

  for (auto _ : state) {
    gauss_seidel(matrix, iterations);
    benchmark::DoNotOptimize(matrix);
  }
}

int main(int argc, char **argv) {
  ::benchmark::Initialize(&argc, argv);

  benchmark::RegisterBenchmark("gaus_seidel_parallel", benchmarkParallel)
      ->Arg(10)
      ->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark("gaus_seidel_serial", benchmarkSerial)
      ->Arg(10)
      ->Unit(benchmark::kMillisecond);

  ::benchmark::RunSpecifiedBenchmarks();

  return 0;
}

// int main(int argc, char **argv) {
//   int iterations = 2;
//   Matrix matrix = Matrix(2402, 2402);
//
//   fill_matrix(matrix, 10);
//
//   Matrix matrix2 = Matrix(2402, 2402);
//
//   fill_matrix(matrix2, 10);
//
//   gauss_seidel_par(matrix, iterations);
//   gauss_seidel(matrix2, iterations);
//
//   check(matrix, matrix2);
//
//   return 0;
// }
