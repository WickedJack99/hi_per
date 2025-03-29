#include <iostream>
#include "matrix.h"
#include "compute.cpp"

int main() {
  std::vector v = std::vector<int>(3);
  v[0] = 10;
  v[1] = 2;
  v[2] = 20;

  Matrix m = Matrix(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int y = 0; y < 3; y++) {
      m(i, y) = 1;
    }
  }

  auto erg = compute(m, v);

  std::cout << v[0] << std::endl;
  std::cout << v[1] << std::endl;
  std::cout << v[2] << std::endl;

  std::cout << std::endl;

  for (int i = 0; i < 3; i++) {
    for (int y = 0; y < 3; y++) {
      std::cout << erg(i, y) << "\t";
    }
    std::cout << std::endl;
  }

  return 0;
}

