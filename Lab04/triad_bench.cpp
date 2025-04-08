#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
#include <string>

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

static std::vector<std::string> ownMethod2(int n) {
  std::vector<double> A(n, 0);
  std::vector<std::string> results(2);
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
  double flops = (2.0 * n) / (averageTime / 1000);
  results[0] = std::to_string(averageTime);
  results[1] = std::to_string(flops);
  return results;
}

static void workerThread(int start, int end, int threadCount) {
  std::vector<std::vector<std::string>> outcomes;
  for (int i = start; i < end; i += 10) {
    outcomes.push_back(ownMethod2(i));
  }
  std::stringstream fileName;
  fileName << "outcomes-" << start << "-" << end << ".csv";
  std::ofstream MyFile(fileName.str());
  for (int i = 0; i < outcomes.size(); i++) {
    MyFile << "(" << outcomes[i][0] << "," << outcomes[i][1] << "," << threadCount << ")";
  }
  MyFile.close();
}

int main() {
  int range_start = 100;
  int range_end = 1'000'000;
  int countThreads[] = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64};
  // int num_threads = std::thread::hardware_concurrency(); // use system's core
  // count
  for (int j = 0; j < (sizeof(countThreads) / sizeof(int)); j++) {
    int total = range_end - range_start;
    int chunk_size = total / countThreads[j];

    std::vector<std::thread> threads;

    for (int i = 1; i < countThreads[j]; ++i) {
      int start = range_start + i * chunk_size;
      int end = (i == countThreads[j] - 1) ? range_end : start + chunk_size;

      threads.emplace_back(workerThread, start, end, countThreads[j]);
    }
    workerThread(range_start, range_end, countThreads[j]);

    for (auto& t : threads)
      t.join();

    std::cout << "All threads finished.\n";
  }

  return 0;
}