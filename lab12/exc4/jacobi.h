#ifndef JACOBI_H
#define JACOBI_H

#include "matrix.h"

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

  Matrix addHaloLayers(const Matrix& phi);
  Matrix removeHaloLayers(const Matrix& phi);

  bool isFirstRank() const;
  bool isLastRank() const;

  int lowerOffset() const;
  int upperOffset() const;

  int numRowsWithHalos(int numRowsWithoutHalos) const;
  int numRowsWithoutHalos(int numRowsWithHalos) const;

  int rank_ = MPI_PROC_NULL;
  int numProc_ = 1;
};

#endif  // JACOBI_H