#ifndef SUPER_GRID_H
#define SUPER_GRID_H

#include "matrix.h"
#include <cstdio>
#include <mpi.h>

class Neighbors {
  public:
    int top_left, top, top_right;
    int left, right;
    int bottom_left, bottom, bottom_right;

    bool operator==(const Neighbors& other) const {
      return top_left     == other.top_left &&
            top          == other.top &&
            top_right    == other.top_right &&
            left         == other.left &&
            right        == other.right &&
            bottom_left  == other.bottom_left &&
            bottom       == other.bottom &&
            bottom_right == other.bottom_right;
    }

  Neighbors()
      : top_left(-1), top(-1), top_right(-1),
        left(-1), right(-1),
        bottom_left(-1), bottom(-1), bottom_right(-1) {}
};

class SuperGrid {
public:
  static SuperGrid zeros(int rows, int cols, MPI_Comm communicator);

  SuperGrid(const Matrix &other);

  double &operator()(int i, int j);
  const double &operator()(int i, int j) const;

  const Matrix get_matrix() const;

  int rows() const;
  int cols() const;

  Neighbors get_neighbors();
  void find_neighbors();

  MPI_Comm &get_communicator();
  void set_communicator(MPI_Comm communicator);

private:
  Matrix grid_;
  MPI_Comm comm_;
  Neighbors neighbors_;
};

inline SuperGrid::SuperGrid(const Matrix &other)
    : grid_(
          Matrix::zeros(other.rows() + 2, other.cols() + 2)),
          neighbors_(Neighbors()) // initialize grid_
{
  for (int i = 0; i < other.rows(); i++) {
    for (int j = 0; j < other.cols(); j++) {
      grid_(i + 1, j + 1) = other(i, j); // copy into grid_ directly
    }
  }
}

inline SuperGrid SuperGrid::zeros(int rows, int cols, MPI_Comm communicator) {
  SuperGrid grid = SuperGrid(Matrix::zeros(rows, cols));
  grid.set_communicator(communicator);
  return grid;
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

inline MPI_Comm &SuperGrid::get_communicator() { return this->comm_; }

inline void SuperGrid::set_communicator(MPI_Comm communicator) {
  this->comm_ = communicator;
}

inline void SuperGrid::find_neighbors() {
  if (this->comm_ == MPI_COMM_NULL) {
    std::cerr << "Communicator is NULL!\n";
  }
  int rank;
  MPI_Comm_rank(this->comm_, &rank);
  int coords[2];
  MPI_Cart_coords(this->comm_, rank, 2, coords);
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      if (dx == 0 && dy == 0) continue; // skip self

      int neighbor_coords[2] = { coords[0] + dx, coords[1] + dy };
      int neighbor_rank;

      int err = MPI_Cart_rank(this->comm_, neighbor_coords, &neighbor_rank);
      if (err != MPI_SUCCESS) continue; // skip invalid neighbors (if non-periodic)

      if (dx == -1 && dy == -1) neighbors_.top_left = neighbor_rank;
      if (dx == -1 && dy ==  0) neighbors_.top      = neighbor_rank;
      if (dx == -1 && dy ==  1) neighbors_.top_right = neighbor_rank;
      if (dx ==  0 && dy == -1) neighbors_.left     = neighbor_rank;
      if (dx ==  0 && dy ==  1) neighbors_.right    = neighbor_rank;
      if (dx ==  1 && dy == -1) neighbors_.bottom_left = neighbor_rank;
      if (dx ==  1 && dy ==  0) neighbors_.bottom   = neighbor_rank;
      if (dx ==  1 && dy ==  1) neighbors_.bottom_right = neighbor_rank;
    }
  }
}

inline Neighbors SuperGrid::get_neighbors() {
  return this->neighbors_;
}

#endif // SUPER_GRID_H
