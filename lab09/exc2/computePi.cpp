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

  if (rank != 0)
  {
    // if not rank 0, send msg to 0 with partial sum
    MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  else
  {
    // if rank 0 receive all partial sums and sum them up
    double recvValue;
    for (int i = 1; i < numTasks; i++)
    {
      MPI_Recv(&recvValue, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      sum += recvValue;
    }
  }

  if (rank == 0)
  {
    double pi = step * sum;
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double, std::milli> diff = end - start;

    std::string filename = "results" + std::to_string(numTasks) + ".txt";
    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(16);
    outFile << "π ≈ " << pi << ", Runtime: " << diff.count() << " ms, Tasks:  " << numTasks << ", Step: " << step << ", Sum: " << sum << std::endl;
    outFile.close();
  }

  MPI_Finalize();
  return 0;
}
