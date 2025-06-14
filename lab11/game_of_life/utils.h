#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <memory>
#include <fstream>
#include <filesystem>

#include "game_of_life.h"
#include "matrix_io.h"
#include "matrix.h"
#include "common.h"

/**
 * Clears the contents of a folder or create it if it does not exist.
 * If the folder exists, the user is asked if he wants to proceed and if yes,
 * all folder contents are removed.
 * If the folder does not exist, it is created.
 */
void clearOrCreateFolder(const std::string& foldername);

/**
 * Stores the animation of the Game of Life in a bunch of text files.
 * Each step of the game is saved in a separate file named "step_<n>.txt",
 * where <n> is the step number.
 * The files are stored in the specified folder.
 */
void storeAnimation(const std::string& foldername,
                    const Matrix& initstate,
                    int numSteps,
                    MPIGridSize mpiProcs);

/** Print the current grid of the game on the console */
void print(const GameOfLife& game);

#endif  // UTILS_H