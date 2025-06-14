#include "patterns.h"

Pattern::Pattern(int rows, int cols, MPIGridSize mpiProcs)
    : mpiProcs_(mpiProcs),
      grid_(Matrix::zeros(rows / np0(), cols / np1())) {
  if (rows <= 0 || cols <= 0) {
    throw std::invalid_argument("Rows and columns must be positive");
  }
  if (rows % np0() != 0) {
    throw std::invalid_argument(
        "Rows must be divisible by the number of processes in the first "
        "dimension");
  }
  if (cols % np1() != 0) {
    throw std::invalid_argument(
        "Columns must be divisible by the number of processes in the second "
        "dimension");
  }

  std::array<int, 2> periods = {1, 1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, mpiProcs.data(), periods.data(), true,
                  &comm_);
}

int Pattern::np0() const {
  return mpiProcs_[0];
}

int Pattern::np1() const {
  return mpiProcs_[1];
}

int Pattern::numRowsLocal() const {
  return grid_.rows();
}

int Pattern::numRowsTotal() const {
  return grid_.rows() * np0();
}

int Pattern::numColsLocal() const {
  return grid_.cols();
}

int Pattern::numColsTotal() const {
  return grid_.cols() * np1();
}

Pattern Pattern::beeHive(int row, int col) {
  if (row < 0 || col < 0 || row + 3 >= numRowsTotal() ||
      col + 3 >= numColsTotal()) {
    throw std::out_of_range("Bee hive pattern exceeds grid bounds");
  }
  setCell(row, col + 1);
  setCell(row, col + 2);
  setCell(row + 1, col);
  setCell(row + 1, col + 3);
  setCell(row + 2, col);
  setCell(row + 2, col + 3);
  setCell(row + 3, col + 1);
  setCell(row + 3, col + 2);
  return *this;
}

Pattern Pattern::glider(int row, int col) {
  if (row < 0 || col < 0 || row + 3 >= numRowsTotal() ||
      col + 3 >= numColsTotal()) {
    throw std::out_of_range("Glider pattern exceeds grid bounds");
  }
  setCell(row, col + 1);
  setCell(row + 1, col + 2);
  setCell(row + 2, col);
  setCell(row + 2, col + 1);
  setCell(row + 2, col + 2);
  return *this;
}

Pattern Pattern::octagon(int row, int col) {
  if (row < 0 || col < 0 || row + 7 >= numRowsTotal() ||
      col + 7 >= numColsTotal()) {
    throw std::out_of_range("Octagon pattern exceeds grid bounds");
  }
  setCell(row + 0, col + 3);
  setCell(row + 0, col + 4);
  setCell(row + 1, col + 2);
  setCell(row + 1, col + 5);
  setCell(row + 2, col + 1);
  setCell(row + 2, col + 6);
  setCell(row + 3, col + 0);
  setCell(row + 3, col + 7);
  setCell(row + 4, col + 0);
  setCell(row + 4, col + 7);
  setCell(row + 5, col + 1);
  setCell(row + 5, col + 6);
  setCell(row + 6, col + 2);
  setCell(row + 6, col + 5);
  setCell(row + 7, col + 3);
  setCell(row + 7, col + 4);

  return *this;
}

Matrix Pattern::getGrid() const {
  return grid_;
}

Pattern Pattern::setCell(int globalRow, int globalCol) {
  if (globalRow < 0 || globalCol < 0 || globalRow >= numRowsTotal() ||
      globalCol >= numColsTotal()) {
    throw std::out_of_range("Cell indices are out of bounds");
  }

  int rank;
  MPI_Comm_rank(comm_, &rank);

  std::array<int, 2> coords;
  MPI_Cart_coords(comm_, rank, 2, coords.data());
  int localRow = globalRow - coords[0] * numRowsLocal();
  int localCol = globalCol - coords[1] * numColsLocal();
  if (localRow < 0 || localCol < 0 || localRow >= grid_.rows() ||
      localCol >= grid_.cols()) {
    return *this;  // Ignore out-of-bounds indices
  }
  grid_(localRow, localCol) = 1;
  return *this;
}
