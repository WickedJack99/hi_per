#include <iostream>
#include "matrix.h"
#include "compute.cpp"
#include "computeNew.cpp"
#include <chrono>

void print(Matrix &m) {
  std::cout << std::endl;

  for (int i = 0; i < 3; i++) {
    for (int y = 0; y < 3; y++) {
      std::cout << m(i, y) << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

Matrix create_matrix(int x, int n) {
  Matrix m = Matrix(3, 3);

  for (int i = 0; i < 3; i++) {
    for (int y = 0; y < 3; y++) {
      m(i, y) = 1;
    }
  }

  return m;
}

std::vector<int> create_vector(int x, int n) {
  std::vector v = std::vector<int>(n);
  for (int i = 0; i < n; i++) {
    v[i] = x;
  }
  return v;
}

int main() {

  std::vector<int> v = create_vector(1, 3);
  Matrix m = create_matrix(1, 1000);

  std::chrono::duration<double, std::milli> sumTimeOld;
  std::chrono::duration<double, std::milli> sumTimeNew;
  auto startTimeOld = std::chrono::system_clock::now();
  volatile auto erg1 = compute(m, v);
  auto endTimeOld = std::chrono::system_clock::now();
  sumTimeOld = endTimeOld - startTimeOld;

  auto startTimeNew = std::chrono::system_clock::now();
  volatile auto erg2 = computeNew(m, v);
  auto endTimeNew = std::chrono::system_clock::now();
  sumTimeNew = endTimeNew - startTimeNew;

  std::cout << "Time old: " << sumTimeOld.count() << "\nTime new: " << sumTimeNew.count() << std::endl;


  return 0;
}

