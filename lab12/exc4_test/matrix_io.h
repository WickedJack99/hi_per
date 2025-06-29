#include "matrix.h"
#include <fstream>
#include <mpi.h>
#include <exception>

/**
 * File i/o for distributed matrices.
 * Matrices may be distributed along the first dimension (i.e. rows are
 * distributed). Only works if number of rows is divisible by number of
 * processes.
 *
 * File format is compatible with python numpy loadtxt
 * First line: # numRows numCols
 * One line of the matrix per line in the file, separated by whitespace
 */
class MatrixIO {
 public:
  MatrixIO();

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
  Matrix loadDistributed(const std::string& filename);

 private:
  Matrix gatherMatrixOnRoot(const Matrix& distributedMatrix);

  Matrix scatterMatrixFromRoot(const Matrix& matrixOnRoot);

  int rank_ = MPI_PROC_NULL;
  int numProc_ = 0;
};
