#include <cmath>
#include <mpi.h>
#include <cassert>

#include "matrix.h"
#include "jacobi.h"
#include "matrix_io.h"

Jacobi::Jacobi()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc_);
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm_);
  MPI_Comm_rank(shm_comm_, &shm_rank_);
  is_first_on_node_ = (shm_rank_ == 0);
}

bool Jacobi::isFirstRank() const
{
  return rank_ == 0;
}

bool Jacobi::isLastRank() const
{
  return rank_ == numProc_ - 1;
}

int Jacobi::numRowsWithHalos(int numRowsWithoutHalos) const
{
  if (isFirstRank() || isLastRank())
  {
    return numRowsWithoutHalos + 1;
  }
  else
  {
    return numRowsWithoutHalos + 2;
  }
}

int Jacobi::numRowsWithoutHalos(int numRowsWithHalos) const
{
  if (isFirstRank() || isLastRank())
  {
    return numRowsWithHalos - 1;
  }
  else
  {
    return numRowsWithHalos - 2;
  }
}

int Jacobi::lowerOffset() const
{
  return isFirstRank() ? 0 : 1;
}

int Jacobi::upperOffset() const
{
  return isLastRank() ? 0 : 1;
}

void Jacobi::exchangeHaloLayers(Matrix &phi)
{
  if (is_first_on_node_)
  {
    exchangeHaloLayersNodeMPIProcFirst(phi);
  }
  else
  {
    exchangeHaloLayersNodeMPIProcSecond(phi);
  }
}

void Jacobi::exchangeHaloLayersNodeMPIProcFirst(Matrix &phi)
{
  int n = phi.rows();
  int sendSize = phi.cols();
  int tag = 0;

  std::vector<MPI_Request> req;

  // Communication with upper partner
  if (!isLastRank())
  {
    int upper = rank_ + 1;
    MPI_Isend(&phi(n - 2, 0), sendSize, MPI_DOUBLE, upper, tag, MPI_COMM_WORLD,
              &req.emplace_back());
    MPI_Irecv(&phi(n - 1, 0), sendSize, MPI_DOUBLE, upper, tag, MPI_COMM_WORLD,
              &req.emplace_back());
  }

  // Communication with lower partner
  if (!isFirstRank())
  {
    SharedmemStates *states = reinterpret_cast<SharedmemStates *>(baseptr_);
    double *shm0 = reinterpret_cast<double *>(states + 1); // row 0
    double *shm1 = shm0 + sendSize;                        // row 1

    // communication with second rank on same node via shared memory
    // We write our send row to shared memory row 0

    while (states->shmStates[0] == SharedmemState::Unread)
    {
      if (rank_ == (numProc_ - 1))
        std::cout << "102" << std::endl;
      MPI_Win_sync(win_);
    }

    for (int j = 0; j < sendSize; ++j)
      shm0[j] = phi(1, j);
    states->shmStates[0] = SharedmemState::Unread;
    MPI_Win_sync(win_); // ensure memory visibility

    // Wait for second proc to write its row back to shared memory row 1
    while (states->shmStates[1] == SharedmemState::Read)
    {
      MPI_Win_sync(win_);
    }

    for (int j = 0; j < sendSize; ++j)
      phi(0, j) = shm1[j]; // halo from second proc
    states->shmStates[1] = SharedmemState::Read;
  }

  // Wait for communication to finish
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(req.size(), req.data(), stat.data());
}

void Jacobi::exchangeHaloLayersNodeMPIProcSecond(Matrix &phi)
{
  int n = phi.rows();
  int sendSize = phi.cols();
  int tag = 0;

  std::vector<MPI_Request> req;

  // Communication with upper partner
  if (!isLastRank())
  {
    // All ranks query pointer to the shared memory (rank 0 gets what it allocated)
    void *raw_ptr;
    MPI_Aint query_size;
    int query_disp_unit;

    // We always query rank 0’s memory (shared across node)
    MPI_Win_shared_query(win_, 0, &query_size, &query_disp_unit, &raw_ptr);

    SharedmemStates *states = reinterpret_cast<SharedmemStates *>(raw_ptr);
    double *shm0 = reinterpret_cast<double *>(states + 1); // row 0
    double *shm1 = shm0 + sendSize;                        // row 1

    // communication with second rank on same node via shared memory
    // We write our send row to shared memory row 0

    while (states->shmStates[1] == SharedmemState::Unread)
    {
      MPI_Win_sync(win_);
    }

    for (int j = 0; j < sendSize; ++j)
      shm1[j] = phi(n - 2, j); // our last inner row
    states->shmStates[1] = SharedmemState::Unread;
    MPI_Win_sync(win_); // ensure memory visibility

    // Wait for first proc to write its row back to shared memory row 0
    while (states->shmStates[0] == SharedmemState::Read)
    {
      MPI_Win_sync(win_);
    }

    for (int j = 0; j < sendSize; ++j)
      phi(n - 1, j) = shm0[j]; // halo from first proc
    states->shmStates[0] = SharedmemState::Read;
  }

  // Communication with lower partner
  if (!isFirstRank())
  {
    int lower = rank_ - 1;
    MPI_Isend(&phi(1, 0), sendSize, MPI_DOUBLE, lower, tag, MPI_COMM_WORLD,
              &req.emplace_back());
    MPI_Irecv(&phi(0, 0), sendSize, MPI_DOUBLE, lower, tag, MPI_COMM_WORLD,
              &req.emplace_back());
  }

  // Wait for communication to finish
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(req.size(), req.data(), stat.data());
}

Matrix Jacobi::addHaloLayers(const Matrix &withoutHalos)
{
  const int numRows = withoutHalos.rows();
  const int numCols = withoutHalos.cols();

  Matrix withHalos = Matrix::zeros(numRowsWithHalos(numRows), numCols);

  for (int i = 0; i < numRows; ++i)
  {
    for (int j = 0; j < numCols; ++j)
    {
      withHalos(i + lowerOffset(), j) = withoutHalos(i, j);
    }
  }
  return withHalos;
}

Matrix Jacobi::removeHaloLayers(const Matrix &withHalos)
{
  const int numRows = numRowsWithoutHalos(withHalos.rows());
  const int numCols = withHalos.cols();

  Matrix withoutHalos = Matrix::zeros(numRows, numCols);

  for (int i = 0; i < numRows; ++i)
  {
    for (int j = 0; j < numCols; ++j)
    {
      withoutHalos(i, j) = withHalos(i + lowerOffset(), j);
    }
  }
  return withoutHalos;
}

Jacobi::Result Jacobi::run(const Matrix &init, double eps, int maxNumIter)
{
  std::vector<Matrix> phi(2, addHaloLayers(init));
  const int numRows = phi[0].rows();
  const int numCols = phi[0].cols();

  MPI_Aint size = is_first_on_node_ ? (numCols * 2 * sizeof(double) + sizeof(SharedmemStates)) : 0;

  MPI_Win_allocate_shared(size, sizeof(char), MPI_INFO_NULL, shm_comm_, &baseptr_, &win_);

  int nIter = 0;
  double dist = std::numeric_limits<double>::max();

  int t0 = 0; // array index of current timestep
  int t1 = 1; // array index of next timestep
  while (dist > eps && nIter < maxNumIter)
  {
    dist = 0;

    exchangeHaloLayers(phi[t0]);

    for (int i = 1; i < numRows - 1; ++i)
    {
      for (int j = 1; j < numCols - 1; ++j)
      {
        phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + phi[t0](i - 1, j) +
                               phi[t0](i, j + 1) + phi[t0](i, j - 1));

        const double diff = phi[t1](i, j) - phi[t0](i, j);
        dist = std::max(dist, std::abs(diff));
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &dist, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    nIter++;
    std::swap(t0, t1);
  }

  return {removeHaloLayers(phi[t0]), dist, nIter};
}
