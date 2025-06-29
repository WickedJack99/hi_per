#ifndef JACOBI_H
#define JACOBI_H

#include "matrix.h"

enum SharedmemState
{
  Unread = 0,
  Read = 1
};

struct SharedmemStates
{
  SharedmemState shmStates[2];        // Flags: one for each row
};

/**
 * Base class for Jacobi algorithms.
 * Defines a struct result which is returned by the algorithm.
 */
class Jacobi {
 public:
  /**
   * Result of the computation
   */
  struct Result {
    Matrix phi = Matrix::zeros(0, 0);  // The converged solution matrix
    double finalDist = 0;              // The final distance
    int numIter = 0;                   // The number of iterations made
  };

  Jacobi();

  /**
   * Run the Jacobi algorithm on the initial matrix init.
   * The algorithm terminates if either the maximum number of iterations
   * exceeds maxNumIter or if the difference between two steps is smaller
   * than eps.
   */
  Result run(const Matrix& init, double eps, int maxNumIter);

 private:
  void exchangeHaloLayers(Matrix& phi);

  void exchangeHaloLayersNodeMPIProcFirst(Matrix& phi);
  void exchangeHaloLayersNodeMPIProcSecond(Matrix& phi);

  Matrix addHaloLayers(const Matrix& phi);
  Matrix removeHaloLayers(const Matrix& phi);

  bool isFirstRank() const;
  bool isLastRank() const;

  int lowerOffset() const;
  int upperOffset() const;

  int numRowsWithHalos(int numRowsWithoutHalos) const;
  int numRowsWithoutHalos(int numRowsWithHalos) const;

  // Global rank on cluster
  int rank_ = MPI_PROC_NULL;

  // Local rank on node
  int shm_rank_ = MPI_PROC_NULL;

  // Communicator for local domain of node
  MPI_Comm shm_comm_;

  // True if this MPI proc is first on node
  bool is_first_on_node_;

  // count of MPI procs
  int numProc_ = 1;

  void* baseptr_;

  MPI_Win win_;

  SharedmemStates shm_states_ = {
    {Read, Read}
  };
};

// 4 times horizontal split | with mpi

// 12 times vertical split - with openmp



#endif  // JACOBI_H