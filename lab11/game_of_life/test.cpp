#include "test.h"
#include "common.h"
#include "game_of_life.h"
#include "matrix.h"
#include "super_grid.h"
#include "utils.h"
#include <iostream>

Matrix init_step(SuperGrid start)
{
  GameOfLife gol = GameOfLife(start, MPIGridSize());
  gol.step();
  return gol.getGrid();
}

SuperGrid init()
{
  MPI_Comm comm_;
  int num_procs = 16;
  MPIGridSize mpiProcs = {4, 4};
  MPI_Dims_create(num_procs, 2, mpiProcs.data());

  std::array<int, 2> periods = {1, 1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, mpiProcs.data(), periods.data(), true,
                  &comm_);

  return SuperGrid::zeros(10, 10, comm_);
}

TEST(initialize)
{
  int a;
  MPI_Init(&a, nullptr);
  check(true, true);
}

TEST(test_underpopulation)
{
  SuperGrid start = init();
  start(0, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 0);
}

TEST(test_survive_2)
{
  SuperGrid start = init();
  start(0, 0) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_survive_3)
{
  SuperGrid start = init();
  start(0, 0) = 1;
  start(1, 1) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_overpopulation)
{
  SuperGrid start = init();
  start(1, 1) = 1;
  start(0, 1) = 1;
  start(1, 0) = 1;
  start(2, 1) = 1;
  start(1, 2) = 1;
  Matrix end = init_step(start);
  check(end(1, 1), 0);
}

TEST(test_reproduction)
{
  SuperGrid start = init();
  start(0, 1) = 1;
  start(1, 0) = 1;
  start(2, 1) = 1;
  Matrix end = init_step(start);
  check(end(1, 1), 1);
}

TEST(test_survive_edge)
{
  SuperGrid start = init();
  start(0, 0) = 1;
  start(9, 0) = 1;
  start(0, 9) = 1;
  Matrix end = init_step(start);
  check(end(0, 0), 1);
}

TEST(test_super_grid)
{
  SuperGrid su_grid = init();
  check(su_grid.rows(), 10);
  check(su_grid.cols(), 10);
}

TEST(test_find_neighbors)
{
  SuperGrid su_grid = init();
  su_grid.find_neighbors();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    Neighbors neighbors = su_grid.get_neighbors();
    Neighbors expected;
    expected.top_left = 15;
    expected.top = 12;
    expected.top_right = 13;
    expected.left = 3;
    expected.right = 1;
    expected.bottom_left = 7;
    expected.bottom = 4;
    expected.bottom_right = 5;
    check(neighbors == expected, true);
  }
  MPI_Finalize();
}

TEST(communication)
{
  SuperGrid su_grid = init();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 1) {
    su_grid(0, 0) = 1;
    su_grid(0, 2) = 1;
    su_grid(0, 4) = 1;
  }
  su_grid.update();
  if (rank == 0) {
    Matrix result = su_grid.get_grid();
    check(result(11, 1), 1);
    check(result(11, 3), 1);
    check(result(11, 5), 1);
  }
}

TEST(get_halo_layers)
{
}

int main(int argc, char *argv[]) { 
  std::cout << "main" << std::endl;
  return 0; }
