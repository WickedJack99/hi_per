#include "jacobi.h"
#include "jacobi_main.h"
#include <iostream>
#include <vector>

std::vector<double> getFilledBuffer(const Matrix &matrix, int row) {
  int numCols = matrix.cols();
  std::vector<double> rowBuffer(numCols);

  // Extract the row (matrix is column-major)
  for (int j = 0; j < numCols; ++j) {
    rowBuffer[j] = matrix(row, j);
  }
  return rowBuffer;
}

Jacobi::Result JacobiMPI::run(const Matrix &init, double epsilon,
                              int maxNumIter) {
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  std::cout << "rank: " << rank << std::endl;
  // std::cout << "numProc: " << numProc << std::endl;

  std::vector<Matrix> phi(2, init);

  const int numRows = phi[0].rows();
  const int numCols = phi[0].cols();

  const int neighborUpper = rank - 1;
  const int neighborLower = rank + 1;

  std::cout << neighborUpper << ":" << neighborLower << std::endl;

  const int indexRowGlobalStart = (rank * numRows);
  const int indexRowGlobalEnd = (rank * numRows) + (numRows - 1);

  const int rowHaloUpperIndex = indexRowGlobalStart - 1;
  const int rowHaloLowerIndex = indexRowGlobalEnd + 1;

  std::vector<double> haloUpper(numCols);
  std::vector<double> haloLower(numCols);

  double dist = std::numeric_limits<double>::max();
  int nIter = 0;

  int t0 = 0; // array index of current timestep
  int t1 = 1; // array index of next timestep
  while (dist > epsilon && nIter < maxNumIter) {
    dist = 0;

    if (rank == 0) {
      MPI_Request send_lower;
      std::vector<double> send_vec_low = phi[t1].get_row(numRows - 1);
      MPI_Isend(&send_vec_low, numCols, MPI_DOUBLE, neighborLower, 0,
                MPI_COMM_WORLD, &send_lower);

      MPI_Request requestLower;
      MPI_Irecv(haloLower.data(), numCols, MPI_DOUBLE, neighborLower, 0,
                MPI_COMM_WORLD, &requestLower);

      MPI_Wait(&send_lower, MPI_STATUS_IGNORE);
      MPI_Wait(&requestLower, MPI_STATUS_IGNORE);

      for (int i = 1; i < numRows; ++i) {
        for (int j = 1; j < numCols - 1; ++j) {
          double valueBottom;
          (i == numRows - 1) ? valueBottom = haloLower[j]
                             : valueBottom = phi[t0](i + 1, j);
          phi[t1](i, j) = .25 * (valueBottom + phi[t0](i - 1, j) +
                                 phi[t0](i, j + 1) + phi[t0](i, j - 1));
          const double diff = phi[t1](i, j) - phi[t0](i, j);
          dist = std::max(dist, std::abs(diff));
        }
      }
    }

    else if (rank == (numProc - 1)) {
      MPI_Request send_upper;
      std::vector<double> send_vec_up = phi[t1].get_row(rank * numRows);
      MPI_Isend(&send_vec_up, numCols, MPI_DOUBLE, neighborUpper, 0,
                MPI_COMM_WORLD, &send_upper);

      MPI_Request request_upper;
      MPI_Irecv(haloUpper.data(), numCols, MPI_DOUBLE, neighborUpper, 0,
                MPI_COMM_WORLD, &request_upper);

      MPI_Wait(&send_upper, MPI_STATUS_IGNORE);
      MPI_Wait(&request_upper, MPI_STATUS_IGNORE);

      for (int i = 0; i < numRows - 1; ++i) {
        for (int j = 1; j < numCols - 1; ++j) {
          double valueTop;
          (i == 0) ? valueTop = haloUpper[j] : valueTop = phi[t0](i - 1, j);
          phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + valueTop +
                                 phi[t0](i, j + 1) + phi[t0](i, j - 1));
          const double diff = phi[t1](i, j) - phi[t0](i, j);
          dist = std::max(dist, std::abs(diff));
        }
      }
    }

    else {
      MPI_Request send_upper;
      MPI_Request send_lower;

      std::vector<double> send_vec_up = phi[t1].get_row(rank * numRows);
      MPI_Isend(&send_vec_up, numCols, MPI_DOUBLE, neighborUpper, 0,
                MPI_COMM_WORLD, &send_upper);

      std::vector<double> send_vec_low =
          phi[t1].get_row(rank * numRows + numRows);
      MPI_Isend(&send_vec_low, numCols, MPI_DOUBLE, neighborLower, 0,
                MPI_COMM_WORLD, &send_lower);

      MPI_Request request_upper;
      MPI_Request request_lower;
      MPI_Irecv(haloUpper.data(), numCols, MPI_DOUBLE, neighborUpper, 0,
                MPI_COMM_WORLD, &request_upper);
      MPI_Irecv(haloLower.data(), numCols, MPI_DOUBLE, neighborLower, 0,
                MPI_COMM_WORLD, &request_lower);

      MPI_Wait(&send_upper, MPI_STATUS_IGNORE);
      MPI_Wait(&send_lower, MPI_STATUS_IGNORE);
      MPI_Wait(&request_upper, MPI_STATUS_IGNORE);
      MPI_Wait(&request_lower, MPI_STATUS_IGNORE);

      for (int i = 0; i < numRows; ++i) {
        for (int j = 1; j < numCols - 1; ++j) {
          if (i == 0) {
            double valueTop = haloUpper[j];
            phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + valueTop +
                                   phi[t0](i, j + 1) + phi[t0](i, j - 1));
            const double diff = phi[t1](i, j) - phi[t0](i, j);
            dist = std::max(dist, std::abs(diff));
          } else if (i == numRows - 1) {
            double valueBottom = haloLower[j];
            phi[t1](i, j) = .25 * (valueBottom + phi[t0](i - 1, j) +
                                   phi[t0](i, j + 1) + phi[t0](i, j - 1));
            const double diff = phi[t1](i, j) - phi[t0](i, j);
            dist = std::max(dist, std::abs(diff));
          } else {
            phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + phi[t0](i - 1, j) +
                                   phi[t0](i, j + 1) + phi[t0](i, j - 1));
            const double diff = phi[t1](i, j) - phi[t0](i, j);
            dist = std::max(dist, std::abs(diff));
          }
        }
      }
    }

    // if (nIter % 1000 == 0) {
    //   std::cout << "Iteration " << nIter << ", dist=" << dist << "\n";
    // }
    double globalDist;
    MPI_Allreduce(&dist, &globalDist, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dist = globalDist;

    nIter++;
    std::swap(t0, t1);
  }

  Matrix result = Matrix::zeros(numProc * numRows, numCols);

  MPI_Gather(&phi[t0], numRows * numCols, MPI_DOUBLE, result.data(),
             numRows * numCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return Jacobi::Result{phi[t0], dist, nIter};
}
