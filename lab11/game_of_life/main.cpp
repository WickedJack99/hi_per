#include <iostream>
#include <fstream>
#include <memory>
#include <mpi.h>
#include <string>
#include <filesystem>

#include "game_of_life.h"
#include "patterns.h"
#include "matrix_io.h"
#include "utils.h"
#include "common.h"

/**
 * Main function to run the simulation of the game of life
 */
void gameOfLife(MPIGridSize mpiProcs) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create a grid that results in some interesting patterns
  const Matrix grid = Pattern(20, 40, mpiProcs)
                          .glider(10, 17)
                          .beeHive(7, 10)
                          .octagon(6, 27)
                          .octagon(12, 0)
                          .getGrid();

  GameOfLife game(grid, mpiProcs);

  if (rank == 0)
    std::cout << "Initial State:" << std::endl;

  print(game);

  for (int i = 0; i < 50; ++i) {
    game.step();
  }
  if (rank == 0)
    std::cout << "Final state" << std::endl;
  print(game);

  storeAnimation("output", grid, 150, mpiProcs);
}

/**
 * Main entry point for the MPI program.
 * Initializes MPI, checks command line arguments, and starts the game of life
 * simulation.
 */
int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  if (argc != 3) {
    std::cout << "Specify number of processes in x and y as arguments\n";
    std::cout << "jacobi <np0> <np1>\n";
    return 1;
  }

  const int np0 = std::stoi(argv[1]);
  const int np1 = std::stoi(argv[2]);
  int numProc;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  if (np0 * np1 != numProc) {
    std::cout << "Error: nproc != np0 x np1 (" << numProc << "!= " << np0 << "x"
              << np1 << ")\n";
    return 2;
  }

  try {
    gameOfLife({np0, np1});
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  MPI_Finalize();
}