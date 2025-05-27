#include "jacobi.h"
#include "jacobi_main.h"

std::vector<double> getFilledBuffer(const Matrix &matrix, int row)
{
  int numCols = matrix.cols();
  std::vector<double> rowBuffer(numCols);

  // Extract the row (matrix is column-major)
  for (int j = 0; j < numCols; ++j)
  {
    rowBuffer[j] = matrix(row, j);
  }
  return rowBuffer;
}

/**
 * Takes a matrix and an index indexRowLocal which defines which row of the local matrix has to be sent.
 *
 * Sends data of that row via MPI_Send to process with rank receiverRank.
 *
 * indexRowGlobal identifies row in the view of the whole matrix, used from receiving process
 * to identify if lower or upper halo by comparing against rowHaloUpperIndex and rowHaloLowerIndex.
 */
void sendLocalRow(MPI_Request &request, std::vector<double> &row, const Matrix &matrix, const int indexRowLocal, const int receiverRank, const int indexRowGlobal)
{
  row = getFilledBuffer(matrix, indexRowLocal);
  MPI_Isend(row.data(), row.size(), MPI_DOUBLE, receiverRank, indexRowGlobal, MPI_COMM_WORLD, &request);
}

Jacobi::Result JacobiMPI::run(const Matrix &init, double epsilon, int maxNumIter)
{
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  std::vector<Matrix> phi(2, init);

  const int numRows = phi[0].rows();
  const int numCols = phi[0].cols();

  const int neighborUpper = rank - 1;
  const int neighborLower = rank + 1;

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
  while (dist > epsilon && nIter < maxNumIter)
  {
    dist = 0;

    if (rank == 0)
    {
      std::vector<double> row;
      MPI_Request request;
      sendLocalRow(request, row, phi[t0], numRows - 1, neighborLower, indexRowGlobalEnd);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      MPI_Request requestLower;
      MPI_Irecv(haloLower.data(), haloLower.size(), MPI_DOUBLE,
                neighborLower, rowHaloLowerIndex, MPI_COMM_WORLD, &requestLower);
      MPI_Wait(&requestLower, MPI_STATUS_IGNORE);

      for (int i = 1; i < numRows; ++i)
      {
        for (int j = 1; j < numCols - 1; ++j)
        {
          double valueBottom;
          (i == numRows - 1) ? valueBottom = haloLower[j] : valueBottom = phi[t0](i + 1, j);
          phi[t1](i, j) = .25 * (valueBottom + phi[t0](i - 1, j) +
                                 phi[t0](i, j + 1) + phi[t0](i, j - 1));
          const double diff = phi[t1](i, j) - phi[t0](i, j);
          dist = std::max(dist, std::abs(diff));
        }
      }
    }

    else if (rank == (numProc - 1))
    {
      std::vector<double> row;

      MPI_Request request;
      sendLocalRow(request, row, phi[t0], 0, neighborUpper, indexRowGlobalStart);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      MPI_Request requestUpper;
      MPI_Irecv(haloUpper.data(), haloUpper.size(), MPI_DOUBLE,
                neighborUpper, rowHaloUpperIndex, MPI_COMM_WORLD, &requestUpper);
      MPI_Wait(&requestUpper, MPI_STATUS_IGNORE);

      for (int i = 0; i < numRows - 1; ++i)
      {
        for (int j = 1; j < numCols - 1; ++j)
        {
          double valueTop;
          (i == 0) ? valueTop = haloUpper[j] : valueTop = phi[t0](i - 1, j);
          phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + valueTop +
                                 phi[t0](i, j + 1) + phi[t0](i, j - 1));
          const double diff = phi[t1](i, j) - phi[t0](i, j);
          dist = std::max(dist, std::abs(diff));
        }
      }
    }

    else
    {
      std::vector<double> rowUpper;
      std::vector<double> rowLower;

      MPI_Request requestSendUpper;
      MPI_Request requestSendLower;
      sendLocalRow(requestSendUpper, rowUpper, phi[t0], 0, neighborUpper, indexRowGlobalStart);
      sendLocalRow(requestSendLower, rowLower, phi[t0], numRows - 1, neighborLower, indexRowGlobalEnd);

      MPI_Request sendRequests[2] = {requestSendUpper, requestSendLower};
      MPI_Waitall(2, sendRequests, MPI_STATUSES_IGNORE);

      MPI_Request requestUpper;
      MPI_Irecv(haloUpper.data(), haloUpper.size(), MPI_DOUBLE,
                neighborUpper, rowHaloUpperIndex, MPI_COMM_WORLD, &requestUpper);
      MPI_Request requestLower;
      MPI_Irecv(haloLower.data(), haloLower.size(), MPI_DOUBLE,
                neighborLower, rowHaloLowerIndex, MPI_COMM_WORLD, &requestLower);

      MPI_Request requests[2] = {requestUpper, requestLower};
      MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

      for (int i = 0; i < numRows; ++i)
      {
        for (int j = 1; j < numCols - 1; ++j)
        {
          if (i == 0)
          {
            double valueTop = haloUpper[j];
            phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + valueTop +
                                   phi[t0](i, j + 1) + phi[t0](i, j - 1));
            const double diff = phi[t1](i, j) - phi[t0](i, j);
            dist = std::max(dist, std::abs(diff));
          }
          else if (i == numRows - 1)
          {
            double valueBottom = haloLower[j];
            phi[t1](i, j) = .25 * (valueBottom + phi[t0](i - 1, j) +
                                   phi[t0](i, j + 1) + phi[t0](i, j - 1));
            const double diff = phi[t1](i, j) - phi[t0](i, j);
            dist = std::max(dist, std::abs(diff));
          }
          else
          {
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

  return Jacobi::Result{phi[t0], dist, nIter};
}
