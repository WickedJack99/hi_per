#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

/**
 * Function to find the maximum element in a vector using multiple threads.
 * 
 * @param maxElement Reference to the variable where the maximum element will be stored.
 * @param v The vector of integers to search.
 */
int findMaxElement(const std::vector<int>& v) {
  int maxElement = std::numeric_limits<int>::min();
  // Find max element in the range
  for (int i = 0; i < v.size(); ++i) {
    if (v[i] > maxElement) {
      maxElement = v[i];
    }
  }
  return maxElement;
}

int main() {
  // Create a vector with some random values
  const int n = 10240;
  std::vector<int> v(n);
  v[0] = 1;
  for (int i = 1; i < n; ++i) {
    v[i] = (75 * v[i - 1]) % 65537;
  }

  // Find the maximum element in the vector
  // This should be done in parallel using multiple threads
  int maxElement = findMaxElement(v);

  std::cout << "Maximum element: " << maxElement << std::endl;

  return 0;
}