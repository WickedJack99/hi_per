#include <benchmark/benchmark.h>
#include <fstream>
#include <iostream>
#include "matrix.h"

Matrix jacobi(const Matrix& init, double eps, int maxNumIter) {
  Matrix *phi = new Matrix(init);
  Matrix *tmp = new Matrix(init);
  const int n = (*phi).dim1();
  const int m = (*phi).dim2();

  double dist = 1e10;
  int nIter = 0;
  while (dist > eps && nIter < maxNumIter) {
    dist = 0;

    for (int i = 1; i < n - 1; ++i) {
      for (int j = 1; j < m - 1; ++j) {
        // if (nIter % 2 == 0) {
          (*tmp)(i, j) = .25 * ((*phi)(i + 1, j) + (*phi)(i - 1, j) +
                               (*phi)(i, j + 1) + (*phi)(i, j - 1));
        // } else {
        //   (phi)(i, j) = .25 * ((tmp)(i + 1, j) + (tmp)(i - 1, j) +
        //                        (tmp)(i, j + 1) + (tmp)(i, j - 1));
        // }
        dist = std::max(dist, std::abs((*tmp)(i, j) - (*phi)(i, j)));
      }
    }
    // tmp = Matrix(phi);
    std::swap(tmp, phi);
    // for (int i = 1; i < n - 1; ++i) {
    //   for (int j = 1; j < n - 1; ++j) {
    //     phi(i, j) = tmp(i,j);
    //   }
    // }

    // std::cout << l2phi << ", " << l2tmp << ", " << dist << std::endl;

    nIter++;
  }

  std::cout << "Finished Jacobi after " << nIter
            << " iterations with an error of " << dist << std::endl;
  return (*phi);
}

Matrix jacobi1(const Matrix& init, double eps, int maxNumIter) {
  Matrix phi(init);
  Matrix tmp(phi.dim());
  const int n = phi.dim1();
  const int m = phi.dim2();

  double dist = 1e10;
  int nIter = 0;
  while (dist > eps && nIter < maxNumIter) {
    dist = 0;
    for (int i = 1; i < n - 1; ++i) {
      for (int j = 1; j < m - 1; ++j) {
        tmp(i, j) = .25 * (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) + phi(i, j - 1));
        dist = std::max(dist, std::abs(tmp(i, j) - phi(i, j)));
      }
    }

    for (int i = 1; i < n - 1; ++i) {
      for (int j = 1; j < n - 1; ++j) {
        phi(i, j) = tmp(i,j);
      }
    }


    // std::cout << l2phi << ", " << l2tmp << ", " << dist << std::endl;

    nIter++;
  }

  std::cout << "Finished Jacobi after " << nIter
            << " iterations with an error of " << dist << std::endl;
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
  const Matrix init = initialCondition(1000);

  Matrix phi;

  for (auto _ : state) {
    phi = jacobi1(init, 1e-5, 100000);
  }
}

BENCHMARK(benchmarkJacobi)->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();
