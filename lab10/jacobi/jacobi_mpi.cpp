#include "jacobi.h"
#include "jacobi_main.h"

Result JacobiMPI::run(const Matrix& init, double epsilon, int maxNumIter){
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  std::vector<Matrix> phi(2, init);
  const int numRows = phi[0].rows();
  const int numCols = phi[0].cols();
  std::vector<double> halo_up(numCols);
  std::vector<double> halo_low(numCols);

  double dist = std::numeric_limits<double>::max();
  int nIter = 0;

  int t0 = 0;  // array index of current timestep
  int t1 = 1;  // array index of next timestep
  while (dist > epsilon && nIter < maxNumIter) {
    dist = 0;

    if (rank == 0)
    {
      MPI_Send(&phi[t0][numRows - 1]);

      for (int i = 1; i < numRows - 1; ++i) {
        for (int j = 1; j < numCols - 1; ++j) {
          phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + phi[t0](i - 1, j) +
                                 phi[t0](i, j + 1) + phi[t0](i, j - 1));
          const double diff = phi[t1](i, j) - phi[t0](i, j);
          dist = std::max(dist, std::abs(diff));
        }
      }

    }
    

    for (int i = 1; i < numRows - 1; ++i) {
      for (int j = 1; j < numCols - 1; ++j) {
        phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + phi[t0](i - 1, j) +
                               phi[t0](i, j + 1) + phi[t0](i, j - 1));
        const double diff = phi[t1](i, j) - phi[t0](i, j);
        dist = std::max(dist, std::abs(diff));
      }
    }

    // if (nIter % 1000 == 0) {
    //   std::cout << "Iteration " << nIter << ", dist=" << dist << "\n";
    // }

    nIter++;
    std::swap(t0, t1);
  }

  return Result{phi[t0], dist, nIter};
}
