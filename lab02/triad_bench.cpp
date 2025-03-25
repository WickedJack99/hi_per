#include <chrono>
#include <vector>
#include <iostream>
#include <fstream> 

std::vector<double> B(1000000, 1);
std::vector<double> C(1000000, 2);
std::vector<double> D(1000000, 3);

static double ownMethod(int n) {
  std::vector<double> A(n, 0);
  
  std::chrono::duration<double, std::milli> sumTime;
  for (int i = 0; i < 20; ++i) {
    auto startTime = std::chrono::system_clock::now();
    for (int i = 0; i < n; ++i) {
      A[i] = B[i] + C[i] * D[i];
    }
    // prevent the compiler from optimizing everything away
    volatile double dummy = A[0];
    auto endTime = std::chrono::system_clock::now();
    sumTime += endTime - startTime;
  }
  double averageTime = sumTime.count() / 20;
  return averageTime;
}

int main() {
  std::vector<double> outcomes;
  for (int i = 100; i < 1000000; i+=10) {
    outcomes.push_back(ownMethod(i));
  }
  std::ofstream MyFile("outcomes.csv");
  for (int j = 0; j < outcomes.size(); j++) {
    MyFile << outcomes[j] << ",";
  }
  MyFile.close();
}