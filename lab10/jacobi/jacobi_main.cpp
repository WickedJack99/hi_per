#include "jacobi_main.h"

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

Matrix run(Jacobi &jacobi, const Matrix &init, double eps, int maxNumIter) {
  auto [m, dist, nIter] = jacobi.run(init, eps, maxNumIter);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Finished running Jacobi, nIter=" << nIter << ", dist=" << dist
              << std::endl;
  }
  return m;
}

double benchmark(Jacobi &jacobi, const Matrix &init, double eps,
                 int maxNumIter) {
  auto start = std::chrono::system_clock::now();
  jacobi.run(init, eps, maxNumIter);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double, std::milli> diff = end - start;
  return diff.count();
}

double serialBenchmark(const int n, double eps, int maxNumIter) {
  // Benchmark the serial version
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double serialTime = 0;
  if (rank == 0) {
    std::cout << "Starting serial benchmark" << std::endl;
    JacobiSerial jacobiSerial;
    serialTime = benchmark(jacobiSerial, getCompleteInitialCondition(n), eps,
                           maxNumIter);

    std::cout << "Serial benchmark finished, time=" << serialTime << "ms"
              << std::endl;
  }
  MPI_Bcast(&serialTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return serialTime;
}

double parallelBenchmark(const int n, double eps, int maxNumIter) {
  JacobiMPI jacobi;
  Matrix init = getDistributedInitialCondition(n);
  double time = benchmark(jacobi, init, eps, maxNumIter);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::cout << "Parallel benchmark finished, time=" << time << "ms"
              << std::endl;
  }
  return time;
}

void runSerial(int n, double eps, int maxNumIter) {
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  // To run the serial version, we just execute the code on rank 0,
  // all other ranks are idle.
  if (rank == 0) {
    std::cout << "Starting serial Jacobi run" << std::endl;
    Matrix init = getCompleteInitialCondition(n);
    const std::string initFilename = "serial_init.asc";
    MatrixIO().saveSerial(init, initFilename);
    JacobiSerial jacobi;
    Matrix phi = run(jacobi, init, eps, maxNumIter);
    const std::string resultFilename = "result_serial.asc";
    MatrixIO().saveSerial(phi, resultFilename);
    std::cout << "Finished serial Jacobi run" << std::endl << std::endl;
  }
}

void runParallel(int n, double eps, int maxNumIter) {
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  if (rank == 0) {
    std::cout << "Starting parallel Jacobi run" << std::endl;
  }

  Matrix init = getDistributedInitialCondition(n);

  const std::string initFilename = "parallel_init.asc";
  MatrixIO().saveDistributed(init, initFilename);

  JacobiMPI jacobi;
  Matrix phi = run(jacobi, init, eps, maxNumIter);
  const std::string resultFilename = "result_parallel.asc";
  MatrixIO().saveDistributed(phi, resultFilename);

  if (rank == 0) {
    std::cout << "Finished parallel Jacobi run" << std::endl;
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  // Set the Jacobi parameters
  const int n = 9 * 64;
  const double eps = 1e-7;
  const int maxNumIter = 100000;

  // Run the serial and parallel versions
  // Result is saved to file.
  // Use this to graphically verify the correctness of the parallel
  // implementation.
  // runSerial(n, eps, maxNumIter);

  runParallel(n, eps, maxNumIter);

  // Run the benchmark
  // double serialTime = 0;
  // double parallelTime = 0;

  // serialTime = serialBenchmark(n, eps, maxNumIter);
  //  parallelTime = parallelBenchmark(n, eps, maxNumIter);
  //
  //  if (rank == 0) {
  //    //std::cout << "Serial time: " << serialTime << "ms" << std::endl;
  //    std::cout << "Serial time: ms" << std::endl;
  //    std::cout << "Parallel time: " << parallelTime << "ms" << std::endl;
  //    //std::cout << "Speedup: " << serialTime / parallelTime << std::endl;
  //    std::ofstream fout("benchmark.txt", std::ios::app);
  //    fout << numProc << "\t" << parallelTime << "\n";
  //  }

  MPI_Finalize();
}
