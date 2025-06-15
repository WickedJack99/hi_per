#include "test.h"
#include "common.h"
#include "game_of_life.h"
#include "matrix.h"
#include "super_grid.h"
#include "utils.h"

Matrix init_step(SuperGrid start) {
  GameOfLife gol = GameOfLife(start, MPIGridSize());
  gol.step();
  return gol.getGrid();
}

TEST(test_underpopulation) {
  SuperGrid start = SuperGrid::zeros(10, 10, nullptr);
  start(0, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 0);
}

TEST(test_survive_2) {
  SuperGrid start = SuperGrid::zeros(10, 10, nullptr);
  start(0, 0) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_survive_3) {
  SuperGrid start = SuperGrid::zeros(10, 10, nullptr);
  start(0, 0) = 1;
  start(1, 1) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_overpopulation) {
  SuperGrid start = SuperGrid::zeros(10, 10, nullptr);
  start(1, 1) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  start(2, 1) = 1;
  start(1, 2) = 1;
  Matrix end = init_step(start);
  check(end(1, 1), 0);
}

TEST(test_reproduction) {
  SuperGrid start = SuperGrid::zeros(10, 10, nullptr);
  start(0, 1) = 1;
  start(1, 0) = 1;
  start(2, 1) = 1;
  Matrix end = init_step(start);
  check(end(1, 1), 1);
}

TEST(test_survive_edge) {
  SuperGrid start = SuperGrid::zeros(10, 10, nullptr);
  start(0, 0) = 1;
  start(9, 0) = 1;
  start(0, 9) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_super_grid) {
  SuperGrid su_grid = SuperGrid::zeros(10, 10, nullptr);
  check(su_grid.rows(), 10);
  check(su_grid.cols(), 10);
}

int main() { return 0; }
