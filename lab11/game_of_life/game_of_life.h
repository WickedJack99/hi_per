#ifndef GAME_OF_LIFE_H
#define GAME_OF_LIFE_H

#include <iostream>
#include <vector>
#include <fstream>
#include <memory>
#include <array>
#include <mpi.h>

#include "common.h"
#include "matrix.h"

class GameOfLife {
 public:
  GameOfLife(const Matrix& grid, MPIGridSize mpiProcs);

  void step();

  Matrix getGrid() const;

  MPIGridSize mpiProcs() const;

 private:
  int countLiveNeighbors(int row, int col) const;

  int updateCell(int currentState, int numLiveNeighbors) const;

  Matrix grid_;
  MPIGridSize mpiProcs_;
  int myRank_ = 0;

  std::array<std::array<int, 3>, 3> neighborRanks_;

  friend class GameOfLifeTest;
};

#endif  // GAME_OF_LIFE_H