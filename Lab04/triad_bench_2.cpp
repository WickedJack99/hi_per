#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
#include <string>
#include <cmath>

const int N = pow(2,17) * 9;

std::vector<double> B(N, 1);
std::vector<double> C(N, 2);
std::vector<double> D(N, 3);

struct Statistics {
  double avg_time;
  double flops;
};

void ownMethod(std::vector<double>& A, int start, int end) {
    for (int j = start; j < end; ++j) {
      A[j] = B[j] + C[j] * D[j];
    }
}

void workerThread(std::vector<double>& A, int start, int end, int threadCount, Statistics& stat) {
  
  std::chrono::duration<double, std::nano> sumTime;
  for (int i = 0; i < 20; i++) {
    auto startTime = std::chrono::system_clock::now();
    ownMethod(A, start, end);
    auto endTime = std::chrono::system_clock::now();
    sumTime += endTime - startTime;
  }
  
  stat.avg_time = sumTime.count() / 20;
  stat.flops = (2.0 * N) / (stat.avg_time  / 1000000);
}

int main() {
  
  std::vector<double> A(N, 0);
  int range_start = 0;
  int range_end =  N;
  int countThreads[] = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64};
  
  // int num_threads = std::thread::hardware_concurrency(); // use system's core
  // count
  for (int j = 0; j < (sizeof(countThreads) / sizeof(int)); j++) {
    std::vector<Statistics> stats(countThreads[j]);

    std::vector<std::thread> threads;

    for (int i = 0; i < countThreads[j]; ++i) {
      int start = N / countThreads[j] * i;
      int end = N / countThreads[j] * (i + 1);

      threads.emplace_back(workerThread, std::ref(A), start, end, countThreads[j], std::ref(stats[i]));
    }

    for (auto& t : threads) {
      t.join();
    }

    //std::cout << A[1179647];
    volatile std::vector<double> dummy(A);

    std::stringstream fileName;
    fileName << "stats" << j << ".csv";
    std::ofstream MyFile(fileName.str());
    for (auto entry : stats) {
      MyFile << entry.avg_time << "," << entry.flops << ";";
    }
    MyFile.close();

    std::cout << "All threads finished.\n";
  }

  return 0;
}