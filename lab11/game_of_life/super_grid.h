#ifndef SUPER_GRID_H
#define SUPER_GRID_H

#include "matrix.h"

class SuperGrid {
public:
  static SuperGrid zeros(int rows, int cols);

  int rows() const;
  int cols() const;
  double &operator()(int i, int j);

private:
  Matrix grid_;
};

inline SuperGrid SuperGrid::zeros(int rows, int cols)
    : grid_(Matrix::zeros(rows + 2, cols + 2)), {}

inline double &SuperGrid::operator()(int i, int j) {
  return this->grid_(i + 1, j + 1);
}

inline int SuperGrid::rows() const { return this->grid_.rows() - 2; }

inline int SuperGrid::cols() const { return this->grid_.cols() - 2; }

#endif // SUPER_GRID_H
