#include <chrono>
#include <vector>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <thread>

std::vector<double> B(1000000, 1);
std::vector<double> C(1000000, 2);
std::vector<double> D(1000000, 3);

static double ownMethod(int n) {
  std::vector<double> A(n, 0);
  std::chrono::duration<double, std::milli> sumTime;
  for (int i = 0; i < 20; ++i) {
    auto startTime = std::chrono::system_clock::now();
    for (int j = 0; j < n; ++j) {
      A[j] = B[j] + C[j] * D[j];
    }
    // prevent the compiler from optimizing everything away
    volatile double dummy = A[0];
    auto endTime = std::chrono::system_clock::now();
    sumTime += endTime - startTime;
  }
  double averageTime = sumTime.count() / 20;
  return averageTime;
}

static void workerThread(int start, int end) {
  std::vector<double> outcomes;
  for (int i = start; i < end; i+=10) {
    outcomes.push_back(ownMethod(i));
  }
  std::stringstream fileName;
  fileName << "outcomes-" << start << "-" << end << ".csv";
  std::ofstream MyFile(fileName.str());
  for (int j = 0; j < outcomes.size(); j++) {
    MyFile << outcomes[j] << ",";
  }
  MyFile.close();
}

int main() {
  int range_start = 100;
  int range_end = 1'000'000;
  workerThread(range_start, range_end);

  // int num_threads = std::thread::hardware_concurrency(); // use system's core count

  // int total = range_end - range_start;
  // int chunk_size = total / num_threads;

  // std::vector<std::thread> threads;

  // for (int i = 0; i < num_threads; ++i) {
  //     int start = range_start + i * chunk_size;
  //     int end = (i == num_threads - 1) ? range_end : start + chunk_size;

  //     threads.emplace_back(workerThread, start, end);
  // }

  // for (auto& t : threads) t.join();

  // std::cout << "All threads finished.\n";
  return 0;
}