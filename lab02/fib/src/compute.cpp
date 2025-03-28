
#include "matrix.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

int fib(int n) {
  if (n < 2)
    return n;

  int a = 0;
  int b = 1;

  for (int i = 2; i <= n; ++i) {
    int temp = a + b;
    a = b;
    b = temp;
  }

  return b;
}

std::vector<double> get_lookup(const std::vector<int> &v) {
  std::vector<double> lookup = std::vector<double>(256);

  for (int i = 0; i < 256; i++) {
    lookup[i] = fib(i);
  }

  return lookup;
}

Matrix compute(const Matrix &s, const std::vector<int> &v) {
  std::vector lookup = get_lookup(v);

  Matrix m(s.dim());
  const int n = v.size();
  for (int j = 0; j < n; ++j) {
    double val = static_cast<double>(lookup[v[j] % 256]);
    double val2 = (sin(val) * tan(val) / sqrt(cos(val) + 2));
    for (int i = 0; i < n; ++i) {
      m(j, i) = s(j, i) * val2;
    }
  }
  return m;
}

