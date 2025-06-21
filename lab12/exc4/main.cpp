#include <iostream>
#include <stdexcept>
#include <chrono>
#include <mpi.h>
#include "jacobi.h"
#include "matrix.h"
#include "matrix_io.h"


Matrix getInitialCondition(int n, int rank, int numProc) {
  if (n % numProc != 0) {
    throw std::runtime_error(
        "Matrix dimension not divisible by number of processes.");
  }

  const int numRowsLocal = n / numProc;
  const int numRowsTotal = n;
  const int numCols = n;

  // const int numCols = n;
  const int i0 = rank * numRowsLocal;
  Matrix m = Matrix::zeros(numRowsLocal, numCols);
  for (int i = 0; i < numRowsLocal; ++i) {
    const double x = 1. * (i + i0) / numRowsTotal;
    for (int j = 0; j < numCols; ++j) {
      const double y = 1. * j / numCols;

      // The lower left corner is hot
      m(i, j) += exp(-0.5 * (x * x + y * y) * 10) * exp(-100 * x * y);

      // The right side is cold
      m(i, j) += -1.0 * exp(-20 * (1 - y));
    }
  }

  return m;
}

Matrix getDistributedInitialCondition(int n) {
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  return getInitialCondition(n, rank, numProc);
}

Matrix getCompleteInitialCondition(int n) {
  return getInitialCondition(n, 0, 1);
}

double benchmark(const Matrix& init,
                 double eps,
                 int maxNumIter,
                 int numBenchmarkRuns) {
  double err = 0.0;
  int numIter = 0;

  auto start = std::chrono::system_clock::now();
  for (int i = 0; i < numBenchmarkRuns; ++i) {
    auto result = Jacobi().run(init, eps, maxNumIter);
    err = result.finalDist;
    numIter = result.numIter;
  }
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double, std::milli> diff = end - start;
  double avgTime = diff.count() / numBenchmarkRuns;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Finished running Jacobi, nIter=" << numIter
              << ", error=" << err << ", time=" << avgTime << "ms" << std::endl;
  }

  return avgTime;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  int n = 3072;
  double eps = 3e-6;
  int maxNumIter = 100000;

  Matrix init = getDistributedInitialCondition(n);
  benchmark(init, eps, maxNumIter, 10);

  MPI_Finalize();
}