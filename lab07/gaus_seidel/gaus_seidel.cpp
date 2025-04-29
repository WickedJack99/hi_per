#include "matrix.h"
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

  const double osth = 1. / 4;

#pragma omp parallel
  {
    for (int iter = 0; iter < maxNumIter; ++iter) {
#pragma omp for
      for (int i = 1; i < m - 1; ++i) {
        int setp = i - omp_get_num_threads();
        for (int j = 1; j < n - 1; ++j) {
          phi(i, j) = osth * (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) +
                              phi(i, j - 1));
        }
        // #pragma omp barrier
      }
    }
  }
}

void fill_matrix(Matrix &matrix, const int filler) {
  const int m = matrix.dim1();
  const int n = matrix.dim2();

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      matrix(i, j) = filler;
    }
  }
}

bool check(Matrix &a, Matrix &b) {
  const int m = a.dim1();
  const int n = a.dim2();

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (a(i, j) != b(i, j)) {
        std::cout << "Not equal at (" << i << ", " << j << ")" << std::endl;
        return false;
      }
    }
  }
  return true;
}

int main() {
  Matrix first = Matrix(10, 10);
  Matrix sec = Matrix(10, 10);

  fill_matrix(first, 10);
  fill_matrix(sec, 10);

  gauss_seidel(first, 1000);
  gauss_seidel_par(sec, 1000);

  std::cout << check(first, sec) << std::endl;

  return 0;
}
