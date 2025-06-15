#include "utils.h"
#include "matrix.h"

void clearOrCreateFolder(const std::string &foldername) {
  namespace fs = std::filesystem;
  fs::path folder(foldername);
  if (fs::exists(folder)) {
    if (!fs::is_directory(folder)) {
      throw std::runtime_error("Path exists but is not a directory: " +
                               foldername);
    }
    if (!fs::is_empty(folder)) {
      char choice;
      std::cout << "Folder '" << foldername
                << "' already exists and is not empty.\n";
      std::cout << "Do you want to clear it? (y/N): ";
      std::cin >> choice;
      if (choice != 'y' && choice != 'Y') {
        throw std::runtime_error(
            "Folder is not empty and user chose not to clear it.");
      }

      std::cout << "Clearing folder '" << foldername << "'...\n";
      // Clear the folder
      for (const auto &entry : fs::directory_iterator(folder)) {
        if (fs::is_directory(entry)) {
          fs::remove_all(entry);
        } else if (fs::is_regular_file(entry)) {
          fs::remove(entry);
        }
      }
    }
  } else {
    // Folder does not exist, create it
    fs::create_directories(folder);
  }
}

void storeAnimation(const std::string &foldername, const Matrix &initstate,
                    int numSteps, MPIGridSize mpiProcs) {
  GameOfLife game(initstate, mpiProcs);

  MatrixIO io(mpiProcs);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Only rank 0 handles file system operations
  if (rank == 0) {
    clearOrCreateFolder(foldername);
  }

  if (rank == 0) {
    std::cout << "Storing animation in folder: " << foldername << std::endl;
  }
  for (int step = 0; step < numSteps; ++step) {
    if (rank == 0) {
      std::string filename =
          foldername + "/step_" + std::to_string(step) + ".txt";
      io.saveDistributed(game.getGrid(), filename);
    }
    game.step();
  }
  if (rank == 0) {
    std::cout << "Animation finished" << std::endl;
  }
}

void print(const GameOfLife &game) {
  MatrixIO io(game.mpiProcs());
  Matrix grid = io.gatherMatrixOnRoot(game.getGrid());

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    for (int i = 0; i < grid.rows(); ++i) {
      for (int j = 0; j < grid.cols(); ++j) {
        std::cout << ((grid(i, j) == 1) ? "X " : ". ");
      }
      std::cout << std::endl;
    }
  }
}

void print(Matrix &grid) {
  for (int i = 0; i < grid.rows(); ++i) {
    for (int j = 0; j < grid.cols(); ++j) {
      std::cout << ((grid(i, j) == 1) ? "X " : ". ");
    }
    std::cout << std::endl;
  }
}
