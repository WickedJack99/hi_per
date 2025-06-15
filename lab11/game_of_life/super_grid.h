#ifndef SUPER_GRID_H
#define SUPER_GRID_H

#include "matrix.h"

class SuperGrid {
public:
  static SuperGrid zeros(int rows, int cols);

  SuperGrid(const Matrix &other);

  double &operator()(int i, int j);
  const double &operator()(int i, int j) const;

  const Matrix get_matrix() const;

  int rows() const;
  int cols() const;

private:
  Matrix grid_;
};

inline SuperGrid::SuperGrid(const Matrix &other)
    : grid_(
          Matrix::zeros(other.rows() + 2, other.cols() + 2)) // initialize grid_
{
  for (int i = 0; i < other.rows(); i++) {
    for (int j = 0; j < other.cols(); j++) {
      grid_(i + 1, j + 1) = other(i, j); // copy into grid_ directly
    }
  }
}

inline SuperGrid SuperGrid::zeros(int rows, int cols) {
  return SuperGrid(Matrix::zeros(rows, cols));
}

inline double &SuperGrid::operator()(int i, int j) {
  return this->grid_(i + 1, j + 1);
}

inline const double &SuperGrid::operator()(int i, int j) const {
  return this->grid_(i + 1, j + 1);
}

inline int SuperGrid::rows() const { return this->grid_.rows() - 2; }

inline int SuperGrid::cols() const { return this->grid_.cols() - 2; }

inline const Matrix SuperGrid::get_matrix() const {
  Matrix mat = Matrix::zeros(this->rows(), this->cols());

  for (int i = 0; i < this->rows(); i++) {
    for (int j = 0; j < this->rows(); j++) {
      mat(i, j) = this->grid_(i + 1, j + 1);
    }
  }

  return mat;
}

#endif // SUPER_GRID_H
