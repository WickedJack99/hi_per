#include "jacobi.h"
#include "jacobi_main.h"
#include <iostream>
#include <iterator>
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

  // Tags für MPI senden/empfangen besser konstant machen
  const int TAG_UP = 0;
  const int TAG_DOWN = 1;

  // Globale Distanz synchronisieren
  while (dist > epsilon && nIter < maxNumIter) {
    dist = 0;

    std::vector<double> send_vec_up, send_vec_low;

    MPI_Request send_upper, send_lower;
    MPI_Request request_upper, request_lower;

    if (rank != 0) {
      send_vec_up = phi[t0].get_row(0); // oberste lokale Zeile
      MPI_Isend(send_vec_up.data(), numCols, MPI_DOUBLE, neighborUpper, TAG_UP,
                MPI_COMM_WORLD, &send_upper);
      MPI_Irecv(haloUpper.data(), numCols, MPI_DOUBLE, neighborUpper, TAG_DOWN,
                MPI_COMM_WORLD, &request_upper);
    }

    if (rank != numProc - 1) {
      send_vec_low = phi[t0].get_row(numRows - 1); // unterste lokale Zeile
      MPI_Isend(send_vec_low.data(), numCols, MPI_DOUBLE, neighborLower,
                TAG_DOWN, MPI_COMM_WORLD, &send_lower);
      MPI_Irecv(haloLower.data(), numCols, MPI_DOUBLE, neighborLower, TAG_UP,
                MPI_COMM_WORLD, &request_lower);
    }

    if (rank != 0) {
      MPI_Wait(&send_upper, MPI_STATUS_IGNORE);
      MPI_Wait(&request_upper, MPI_STATUS_IGNORE);
    }

    if (rank != numProc - 1) {
      MPI_Wait(&send_lower, MPI_STATUS_IGNORE);
      MPI_Wait(&request_lower, MPI_STATUS_IGNORE);
    }

    for (int i = 0; i < numRows; ++i) {
      for (int j = 1; j < numCols - 1; ++j) {
        double top = (i == 0) ? (rank == 0 ? phi[t0](i, j) : haloUpper[j])
                              : phi[t0](i - 1, j);
        double bottom =
            (i == numRows - 1)
                ? (rank == numProc - 1 ? phi[t0](i, j) : haloLower[j])
                : phi[t0](i + 1, j);

        phi[t1](i, j) =
            0.25 * (top + bottom + phi[t0](i, j - 1) + phi[t0](i, j + 1));
        double diff = phi[t1](i, j) - phi[t0](i, j);
        dist = std::max(dist, std::abs(diff));
      }
    }

    // Globale max-Distanz berechnen
    double globalDist = 0;
    MPI_Allreduce(&dist, &globalDist, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dist = globalDist;

    nIter++;
    std::swap(t0, t1);
  }

  // std::cout << phi[t1] << std::endl;
  return Jacobi::Result{phi[t0], dist, nIter};
}
