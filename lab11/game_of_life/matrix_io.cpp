#include <fstream>
#include <mpi.h>
#include <exception>
#include "matrix.h"
#include "matrix_io.h"

MatrixIO::MatrixIO(MPIGridSize mpiProcs) : mpiProcs_(mpiProcs) {
  std::array<int, 2> periods = {1, 1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, mpiProcs.data(), periods.data(), true,
                  &comm_);
  MPI_Comm_rank(comm_, &rank_);
}

void MatrixIO::saveDistributed(const Matrix& distributedMatrix,
                               const std::string& filename) {
  Matrix mat = gatherMatrixOnRoot(distributedMatrix);

  if (rank() == 0) {
    saveSerial(mat, filename);
  }
}

void MatrixIO::saveSerial(const Matrix& m, const std::string& filename) {
  std::ofstream fout(filename);
  const int numRows = m.rows();
  const int numCols = m.cols();
  fout << "# " << numRows << "\t" << numCols << "\n";
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      fout << m(i, j) << "\t";
    }
    fout << "\n";
  }
}

Matrix MatrixIO::load(const std::string& filename) {
  Matrix matrixOnRoot = Matrix::zeros(0, 0);  // Initialize empty matrix
  if (rank() == 0) {
    matrixOnRoot = loadSerial(filename);
  }

  Matrix distributedMatrix = scatterMatrixFromRoot(matrixOnRoot);
  return distributedMatrix;
}

Matrix MatrixIO::loadSerial(const std::string& filename) {
  std::ifstream fin(filename);

  // Check first character (has to be #)
  std::string s;
  fin >> s;
  if (s != "#") {
    throw std::runtime_error("Error, not reading expecte character #\n");
  }

  // Read in number of rows and cols and create matrix
  int numRows, numCols;
  fin >> numRows >> numCols;
  Matrix mat = Matrix::zeros(numRows, numCols);

  // Read in matrix contents
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      fin >> mat(i, j);
    }
  }

  return mat;
}

Matrix MatrixIO::gatherMatrixOnRoot(const Matrix& distributedMatrix) {
  const int numRowsLocal = distributedMatrix.rows();
  const int numColsLocal = distributedMatrix.cols();
  Matrix subMatrix(distributedMatrix);
  Matrix matrixOnRoot = Matrix::zeros(0, 0);  // Initialize empty matrix
  if (rank() == 0) {
    const int numRowsTotal = numRowsLocal * mpiProcs_[0];
    const int numColsTotal = numColsLocal * mpiProcs_[1];
    matrixOnRoot = Matrix::zeros(numRowsTotal, numColsTotal);

    for (int proc = 0; proc < nProc(); ++proc) {
      if (proc != 0) {
        MPI_Recv(&subMatrix(0, 0), subMatrix.numEntries(), MPI_DOUBLE, proc,
                 proc, comm_, MPI_STATUS_IGNORE);
      }

      int coords[2];
      MPI_Cart_coords(comm_, proc, 2, coords);
      const int i0 = coords[0] * numRowsLocal;
      const int j0 = coords[1] * numColsLocal;
      for (int i = 0; i < numRowsLocal; ++i) {
        for (int j = 0; j < numColsLocal; ++j) {
          matrixOnRoot(i + i0, j + j0) = subMatrix(i, j);
        }
      }
    }
  } else {
    MPI_Send(&subMatrix(0, 0), subMatrix.numEntries(), MPI_DOUBLE, 0, rank(),
             comm_);
  }

  return matrixOnRoot;
}

Matrix MatrixIO::scatterMatrixFromRoot(const Matrix& matrixOnRoot) {
  int numRowsTotal = matrixOnRoot.rows();
  int numColsTotal = matrixOnRoot.cols();

  MPI_Bcast(&numRowsTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&numColsTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int numRowsLocal = numRowsTotal / mpiProcs_[0];
  int numColsLocal = numColsTotal / mpiProcs_[1];

  Matrix distributedMatrix = Matrix::zeros(numRowsLocal, numColsLocal);

  if (rank() == 0) {
    for (int proc = nProc() - 1; proc >= 0; --proc) {  // iterate backwards
      int coords[2];
      MPI_Cart_coords(comm_, proc, 2, coords);
      const int i0 = coords[0] * numRowsLocal;
      const int j0 = coords[1] * numColsLocal;
      for (int i = 0; i < numRowsLocal; ++i) {
        for (int j = 0; j < numColsLocal; ++j) {
          distributedMatrix(i, j) = matrixOnRoot(i + i0, j + j0);
        }
      }

      if (proc != 0) {
        MPI_Send(&distributedMatrix(0, 0), distributedMatrix.numEntries(),
                 MPI_DOUBLE, proc, proc, comm_);
      }
    }

  } else {
    MPI_Recv(&distributedMatrix(0, 0), distributedMatrix.numEntries(),
             MPI_DOUBLE, 0, rank(), comm_, MPI_STATUS_IGNORE);
  }

  return distributedMatrix;
}
