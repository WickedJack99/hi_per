#include <cmath>
#include <iostream>
#include <thread>
#include "matrix.h"
#include "test.h"

// Create the matrix and vector to be multiplied and fill them
// with some sensible initial values.
std::pair<Matrix, std::vector<double>> createMatrixAndVector() {
  const int n = 1e3*9;
  Matrix mat(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      mat(i, j) = pow(-1, i) * (i + j);
    }
  }

  std::vector<double> vec(n);
  for (int i = 0; i < n; ++i) {
    vec[i] = 1. / (i + 1);
  }

  return std::pair(mat, vec);
}

// Verify that the computed result is correct. Rather inefficient,
// since it runs on a single core.
void verifyResult(const std::vector<double> result) {
  auto [mat, vec] = createMatrixAndVector();
  const int n = vec.size();

  for (int i = 0; i < n; ++i) {
    double expected = 0;
    for (int j = 0; j < n; ++j) {
      expected += mat(i, j) * vec[j];
    }
    check(result[i], expected);
  }
}

void computeResult(const Matrix& mat, const std::vector<double>& vec, int start, int end, std::vector<double>& result) {
  int n = vec.size();
  for (int x = start; x < end; x++) {
    for (int y = 0; y < n; y++) {
      result[y] = mat(x,y) * vec[x];
    }
  }
}

int main() {
  auto [mat, vec] = createMatrixAndVector();
  std::vector<double> result(vec.size(), 0);
  result = mat * vec;
  // TODO: compute result = mat * vec with multiple threads

  verifyResult(result);
}