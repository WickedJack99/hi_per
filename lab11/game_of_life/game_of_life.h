#ifndef GAME_OF_LIFE_H
#define GAME_OF_LIFE_H

#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <vector>

#include "common.h"
#include "matrix.h"
#include "super_grid.h"

class GameOfLife {
public:
  GameOfLife(const SuperGrid &grid, MPIGridSize mpiProcs);

  void step();

  Matrix getGrid() const;

  MPIGridSize mpiProcs() const;

private:
  int countLiveNeighbors(int row, int col) const;

  int updateCell(int currentState, int numLiveNeighbors) const;

  SuperGrid grid_;
  MPIGridSize mpiProcs_;
  int myRank_ = 0;

  std::array<std::array<int, 3>, 3> neighborRanks_;

  friend class GameOfLifeTest;
};

#endif // GAME_OF_LIFE_H
