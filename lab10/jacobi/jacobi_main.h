#ifndef JACOBI_MAIN_H
#define JACOBI_MAIN_H

#include <iostream>
#include <exception>

#include <mpi.h>
#include <chrono>

#include "matrix.h"
#include "jacobi.h"
#include "jacobi_serial.h"
#include "jacobi_mpi.h"
#include "matrix_io.h"

/**
 * Generate initial condition for Jacobi iteration.
 * Data is already appropriately distributed.
 */
Matrix getInitialCondition(int n, int rank, int numProc);

/**
 * Generate initial condition for Jacobi iteration.
 * Data is distriburted along the first axis.
 */
Matrix getDistributedInitialCondition(int n);

/**
 * Generate initial condition for Jacobi iteration.
 * Data is not distributed, i.e. all ranks have the same data.
 */
Matrix getCompleteInitialCondition(int n);

/**
 * Create an initial condition and run Jacobi algorithm until convergence.
 * @return The converged steady state.
 */
Matrix run(Jacobi& jacobi,
           const Matrix& init,
           double eps,
           int maxNumIter);

/**
 * Perform a single benchmark run: create a initial condition and measure
 * the time until Jacobi converges.
 * 
 * Allows to compare the performance of the serial and parallel
 * Jacobi implementations by taking an instance of the Jacobi base class.
 * 
 * @return time for jacobi convergence in ms
 */
double benchmark(Jacobi& jacobi,
                 const Matrix& init,
                 double eps,
                 int maxNumIter);
/**
 * Run the serial Jacobi algorithm.
 * 
 * Create the initial condition, run the Jacobi algorithm until convergence,
 * and save the result to a file.
 */
void runSerial(int n, double eps, int maxNumIter);

/**
 * Run the parallel Jacobi algorithm.
 *
 * Create the initial condition, run the Jacobi algorithm until convergence,
 * and save the result to a file.
 */
void runParallel(int n, double eps, int maxNumIter);

int main(int argc, char* argv[]);

#endif  // JACOBI_MAIN_H