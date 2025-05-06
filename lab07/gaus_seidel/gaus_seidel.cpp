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

  for (int iter = 0; iter < maxNumIter; ++iter) {
#pragma omp parallel num_threads(3)
    {
      int num_theads = omp_get_num_threads();
      int chunk = (m - 2) / num_theads;
      int start = 1 + chunk * omp_get_thread_num();
      int end = chunk * (omp_get_thread_num() + 1);

      // printf("thread %d, start: %d, end: %d\n", omp_get_thread_num(), start,
      //        end);

      for (int j = 1; j < n + omp_get_num_threads() - 1; ++j) {
        for (int i = start; i <= end; ++i) {

          int k = j - omp_get_thread_num();

          // printf("thread %d, i: %d, j: %d, k: %d\n", omp_get_thread_num(), i,
          // j,
          //        k);

          if (k > 0 && k < n - 1) {
            phi(i, k) = osth * (phi(i + 1, k) + phi(i - 1, k) + phi(i, k + 1) +
                                phi(i, k - 1));
          }
#pragma omp barrier
        }
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

int main() {
  Matrix first = Matrix(8, 8);
  Matrix sec = Matrix(8, 8);

  fill_matrix(first, 10);
  fill_matrix(sec, 10);

  print_matrix(first);
  print_matrix(sec);

  gauss_seidel(first, 12);
  gauss_seidel_par(sec, 12);

  check(first, sec);

  print_matrix(first);
  print_matrix(sec);

  return 0;
}
