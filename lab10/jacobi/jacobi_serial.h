#ifndef JACOBI_SERIAL_H
#define JACOBI_SERIAL_H


#include "jacobi.h"

class JacobiSerial : public Jacobi {
 public:
  Result run(const Matrix& phi, double epsilon, int maxNumIter);
};

#endif  // JABOBI_SERIAL_H
