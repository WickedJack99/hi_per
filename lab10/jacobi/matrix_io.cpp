#include "matrix_io.h"
#include "matrix.h"
#include <exception>
#include <fstream>
#include <mpi.h>

MatrixIO::MatrixIO() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc_);
}

void MatrixIO::saveDistributed(const Matrix &distributedMatrix,
                               const std::string &filename) {
  Matrix mat = gatherMatrixOnRoot(distributedMatrix);

  if (rank_ == 0) {
    saveSerial(mat, filename);
  }
}

void MatrixIO::saveSerial(const Matrix &m, const std::string &filename) {
  std::ofstream fout(filename);
  const int numRows = m.rows();
  const int numCols = m.cols();
  // fout << numRows << "," << numCols << "\n";
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      fout << m(i, j) << ",";
    }
    fout << "\n";
  }
}

Matrix MatrixIO::loadDistributed(const std::string &filename) {
  Matrix matrixOnRoot = Matrix::zeros(0, 0);
  if (rank_ == 0) {
    matrixOnRoot = loadSerial(filename);
  }

  Matrix distributedMatrix = scatterMatrixFromRoot(matrixOnRoot);
  return distributedMatrix;
}

Matrix MatrixIO::loadSerial(const std::string &filename) {
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

Matrix MatrixIO::gatherMatrixOnRoot(const Matrix &distributedMatrix) {
  const int numRowsLocal = distributedMatrix.rows();
  const int numCols = distributedMatrix.cols();
  const int numRowsTotal = numRowsLocal * numProc_;

  Matrix matrixOnRoot = Matrix::zeros(0, 0);
  if (rank_ == 0) {
    matrixOnRoot = Matrix::zeros(numRowsTotal, numCols);
  }

  const int sendSize = distributedMatrix.numEntries();

  MPI_Gather(&distributedMatrix(0, 0), sendSize, MPI_DOUBLE,
             &matrixOnRoot(0, 0), sendSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return matrixOnRoot;
}

Matrix MatrixIO::scatterMatrixFromRoot(const Matrix &matrixOnRoot) {
  int numRowsTotal = matrixOnRoot.rows();
  int numCols = matrixOnRoot.cols();

  MPI_Bcast(&numRowsTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&numCols, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (numRowsTotal % numProc_ != 0)
    throw std::runtime_error(
        "Number of rows not divisible by number of processes.");

  const int numRowsLocal = numRowsTotal / numProc_;

  Matrix distributedMatrix = Matrix::zeros(numRowsLocal, numCols);

  const int sendSize = distributedMatrix.numEntries();

  MPI_Scatter(matrixOnRoot.data(), sendSize, MPI_DOUBLE,
              distributedMatrix.data(), sendSize, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  const int i0 = rank_ * numRowsLocal;

  for (int i = 0; i < numRowsLocal; ++i) {
    std::cout << "Scattered Line " << i0 + i << ": "
              << distributedMatrix(i, numCols - 1) << "\n";
  }

  return distributedMatrix;
}
