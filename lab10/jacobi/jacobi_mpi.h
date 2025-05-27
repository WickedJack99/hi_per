#ifndef JACOBI_MPI_H
#define JACOBI_MPI_H

#include "jacobi.h"

class JacobiMPI : public Jacobi {
 public:
  Result run(const Matrix& init, double epsilon, int maxNumIter) {
    return Result{};
  }
};

#endif  // JABOBI_MPI_H
