#include "test.h"
#include "common.h"
#include "game_of_life.h"
#include "matrix.h"

Matrix init_step(Matrix start) {
  GameOfLife gol = GameOfLife(start, MPIGridSize());
  gol.step();
  return gol.getGrid();
}

TEST(test_underpopulation) {
  Matrix start = Matrix::zeros(10);
  start(0, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 0);
}

TEST(test_survive) {
  Matrix start = Matrix::zeros(10);
  start(0, 0) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_overpopulation) {
  Matrix start = Matrix::zeros(10);
  start(1, 1) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  start(2, 1) = 1;
  start(1, 2) = 1;
  Matrix end = init_step(start);
  check(end(1, 1), 0);
}

TEST(test_reproduction) {
  Matrix start = Matrix::zeros(10);
  start(0, 1) = 1;
  start(1, 0) = 1;
  start(2, 1) = 1;
  Matrix end = init_step(start);
  check(end(1, 1), 1);
}

int main() { return 0; }
