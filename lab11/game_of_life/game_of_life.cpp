#include "game_of_life.h"

GameOfLife::GameOfLife(const SuperGrid &grid, MPIGridSize mpiProcs)
    : grid_(grid), mpiProcs_(mpiProcs) {}

void GameOfLife::step() {
  SuperGrid next = SuperGrid::zeros(grid_.rows(), grid_.cols(), grid_.get_communicator());
  const int rows = grid_.rows();
  const int cols = grid_.cols();
  grid_.update();
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      const int numLiveNeighbors = countLiveNeighbors(i, j);
      next(i, j) = updateCell(grid_(i, j), numLiveNeighbors);
    }
  }

  grid_ = next;
}

int GameOfLife::countLiveNeighbors(int row, int col) const {
  int count = 0;
  const int rows = grid_.rows();
  const int cols = grid_.cols();
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      if (i == 0 && j == 0)
        continue;                            // Skip the cell itself
      int nextRow = (row + i + rows) % rows; // Wrap around
      int nextCol = (col + j + cols) % cols; // Wrap around
      if (grid_(nextRow, nextCol) == 1) {
        count++;
      }
    }
  }

  return count;
}

int GameOfLife::updateCell(int currentState, int numLiveNeighbors) const {
  if (numLiveNeighbors == 3) {
    return 1;
  } else if (numLiveNeighbors == 2) {
    return currentState;
  } else {
    return 0;
  }
}

Matrix GameOfLife::getGrid() const { return grid_.get_matrix(); }

MPIGridSize GameOfLife::mpiProcs() const { return mpiProcs_; }
