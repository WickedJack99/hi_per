#include <chrono>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <fstream>

double computePartialSum(int64_t i0, int64_t i1, double step)
{
  // compute part of the integral
  double sum = 0.0;
  for (int64_t i = i0; i < i1; i++)
  {
    double x = (i + 0.5) * step;
    sum += 4.0 / (1.0 + x * x);
  }
  return sum;
}

int main(int argc, char *argv[])
{
  auto start = std::chrono::system_clock::now();

  MPI_Init(&argc, &argv);

  int rank, numTasks;
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // roughly 1e11, but divisible through 2, 4, 6, 8, 12, 16, 18, ...
  const int64_t numSteps = 1024LL * 1024 * 1024 * 9 * 9;

  // compute boundaries
  int64_t chunkSize = numSteps / numTasks;
  int64_t i0 = rank * chunkSize;
  int64_t i1 = (rank == numTasks - 1) ? numSteps : i0 + chunkSize;

  double step = 1.0 / (double)numSteps;

  // compute partial sum
  double sum = computePartialSum(i0, i1, step);

  double globalSum = 0.0;
  MPI_Reduce(&sum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    double pi = step * globalSum;
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double, std::milli> diff = end - start;

    std::string filename = "results_reduce" + std::to_string(numTasks) + ".txt";
    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(16);
    outFile << "π ≈ " << pi << ", Runtime: " << diff.count() << " ms, Tasks:  " << numTasks << ", Step: " << step << ", Sum: " << sum << std::endl;
    outFile.close();
  }

  MPI_Finalize();
  return 0;
}
