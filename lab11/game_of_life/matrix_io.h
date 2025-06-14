#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include <fstream>
#include <mpi.h>
#include <exception>
#include <memory>
#include "matrix.h"
#include "common.h"

/**
 * File i/o for distributed matrices.
 * Matrices is distributed along 2 dimensions (i.e. rows and columns are
 * distributed).
 *
 * File format is compatible with python numpy loadtxt
 * First line: # numRows numCols
 * One line of the matrix per line in the file, separated by whitespace
 */
class MatrixIO {
 public:
  MatrixIO(MPIGridSize mpiProcs);

  /**
   * Store a matrix on disk, that is not distribued. Should be called by a
   * single rank only.
   */
  void saveSerial(const Matrix& m, const std::string& filename);

  /**
   * Load a matrix from disk, that is not distribued. Should be called by a
   * single rank only.
   */
  Matrix loadSerial(const std::string& filename);

  /**
   * Store a distributed matrix on disk. Should be called collectively by
   * all ranks.
   */
  void saveDistributed(const Matrix& m, const std::string& filename);

  /**
   * Load a matrix from disk and distribute it along the first axis.
   * Should be called collectively by all ranks.
   */
  Matrix load(const std::string& filename);

  /**
   * Gather a distributed matrix on the root rank.
   * The matrix is gathered in a single matrix on the root rank.
   * Other ranks return an empty matrix.
   * Should be called collectively by all ranks.
   */
  Matrix gatherMatrixOnRoot(const Matrix& distributedMatrix);

  /**
   * Scatter a matrix from the root rank to all other ranks.
   * The matrix is scattered along the first axis.
   * Should be called collectively by all ranks.
   */
  Matrix scatterMatrixFromRoot(const Matrix& matrixOnRoot);

 private:
  int rank() const { return rank_; };
  int nProc() const {
    int size;
    MPI_Comm_size(comm_, &size);
    return size;
  }

  int rank_ = 0;
  MPIGridSize mpiProcs_ = {0, 0};  // Number of processes in each dimension

  MPI_Comm comm_ = MPI_COMM_NULL;  // The communicator for the matrix operations
};

#endif  // MATRIX_IO_H