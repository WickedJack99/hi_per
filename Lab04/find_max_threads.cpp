#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <thread>

/**
 * Function to find the maximum element in a vector using multiple threads.
 * 
 * @param maxElement Reference to the variable where the maximum element will be stored.
 * @param v The vector of integers to search.
 */
void findMaxElement(const std::vector<int>& v, const int start, const int end, int& result) {
  int maxElement = std::numeric_limits<int>::min();
  // Find max element in the range
  for (int i = start; i < end; ++i) {
    if (v[i] > maxElement) {
      maxElement = v[i];
    }
  }
  result = maxElement;
}

int main() {
  // Create a vector with some random values
  const int n = 10240;
  const int threadCounter = 4;
  std::vector<std::thread> threads;
  std::vector<int> res(4, 0);
  std::vector<int> v(n);
  v[0] = 1;
  for (int i = 1; i < n; ++i) {
    v[i] = (75 * v[i - 1]) % 65537;
  }

  for(int i = 0; i < threadCounter; i++) {
    int start = n / threadCounter * i;
    int end = n / threadCounter * (i + 1);

    threads.emplace_back(findMaxElement, std::ref(v), start, end, std::ref(res[i]));
  }
  
  for(auto&t : threads) {
    t.join();
  }

  int max = res[0];
  for(int i = 1; i < res.size(); i++) {
    if (max < res[i]){
      max = res[i];
    }
  }

  // Find the maximum element in the vector
  // This should be done in parallel using multiple threads
  //int maxElement = findMaxElement(v);

  std::cout << "Maximum element: " << max << std::endl;

  return 0;
}
