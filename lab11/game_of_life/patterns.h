#ifndef PATTERNS_H
#define PATTERNS_H

#include "common.h"
#include "matrix.h"
#include "super_grid.h"
#include <array>
#include <memory>
#include <mpi.h>

/**
 * Class to create some patterns for the Game of Life.
 *
 * This class allows you to create predefined patterns like bee hive, glider,
 * octagon, etc., on a grid that is distributed across multiple processes.
 *
 * Pattern creation methods return a reference to the Pattern object itself,
 * allowing for method chaining.
 *
 * Example usage for creating a single glider on a 20x40 grid:
 * Pattern pattern(20, 40, comm).glider(10, 17).getGrid();
 */
class Pattern {
public:
  Pattern(int rows, int cols, MPIGridSize mpiProcs);

  SuperGrid getGrid() const;

  Pattern setCell(int globalRow, int globalCol);

  Pattern beeHive(int row, int col);

  Pattern glider(int row, int col);

  Pattern octagon(int row, int col);

private:
  int numRowsLocal() const;
  int numRowsTotal() const;
  int numColsLocal() const;
  int numColsTotal() const;

  int np0() const;
  int np1() const;

  MPIGridSize mpiProcs_;
  SuperGrid grid_;

  MPI_Comm comm_;
};

#endif // PATTERNS_H
