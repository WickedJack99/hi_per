#include "matrix.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdexcept>

Matrix jacobi(const Matrix &init, double eps, int maxNumIter) {
  std::vector<Matrix> phi(2, init);
  const int n = phi[0].dim1();
  const int m = phi[0].dim2();

  double dist = std::numeric_limits<double>::max();
  int nIter = 0;

  int t0 = 0;
  int t1 = 1;
  while (dist > eps && nIter < maxNumIter) {
    dist = 0;
    for (int i = 1; i < n - 1; ++i) {
      for (int j = 1; j < m - 1; ++j) {
        phi[t1](i, j) = .25 * (phi[t0](i + 1, j) + phi[t0](i - 1, j) +
                               phi[t0](i, j + 1) + phi[t0](i, j - 1));
        const double diff = phi[t1](i, j) - phi[t0](i, j);
        dist = std::max(dist, std::abs(diff));
      }
    }

    nIter++;
    std::swap(t0, t1);
  }

  std::cout << "Finished Jacobi after " << nIter
            << " iterations with an error of " << dist << std::endl;
  return phi[t0];
}

Matrix initialCondition(int n) {
  Matrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const double x = 1. * i / n;

    for (int j = 0; j < n; ++j) {
      const double y = 1. * j / n;

      // The lower left corner is hot
      m(i, j) += exp(-0.5 * (x * x + y * y) * 10) * exp(-100 * x * y);

      // The right side is cold
      m(i, j) += -1.0 * exp(-20 * (1 - y));
    }
  }
  return m;
}

/**
 * Store a matrix to a file.
 * May be loaded with python numpy.loadtxt() or with the function loadMatrix
 * below.
 */
void storeMatrix(const Matrix &mat, const std::string &filename) {
  std::ofstream fout(filename);
  // first line: store keyword "#matrix" and dimensions
  // the hash allows loading the matrix with pythons loadtxt
  fout << "#matrix " << mat.dim1() << " " << mat.dim2() << "\n";

  // Store entries with maximum precisions
  fout << std::setprecision(16);

  // Store matrix entry by entry
  for (int i = 0; i < mat.dim1(); ++i) {
    for (int j = 0; j < mat.dim2(); ++j) {
      fout << mat(i, j) << " ";
    }
    fout << "\n";
  }
}

/**
 * Load a matrix from a file that has been stored with storeMatrix
 */
Matrix loadMatrix(const std::string &filename) {
  std::ifstream fin(filename);

  // Read in matrix dimensions
  std::string s;
  int m, n;
  fin >> s >> m >> n;
  if (s != "#matrix") {
    throw std::runtime_error(std::string("File '") + filename +
                             "' does not contain a matrix");
  }

  Matrix mat(m, n);
  for (int i = 0; i < mat.dim1(); ++i) {
    for (int j = 0; j < mat.dim2(); ++j) {
      fin >> mat(i, j);
    }
  }
  return mat;
}

/**
 * Compare the matrix computed to a matrix loaded from a file.
 * @param computed A matrix to be compared
 * @param referencefilename The file from which the reference matrix is loaded
 * @param eps the accuraccy from which on the comparison fails
 * @return: true if computed and the loaded matrix are identical
 */
bool verifyAgainstStoredReference(const Matrix &computed,
                                  const std::string &referenceFilename,
                                  double eps = 1e-8) {
  Matrix ref = loadMatrix(referenceFilename);
  if (!(ref.dim1() == computed.dim1() && ref.dim2() == computed.dim2())) {
    std::cout << "Dimension error";
    return false;
  }

  bool correct = true;
  for (int i = 0; i < ref.dim1(); ++i) {
    for (int j = 0; j < ref.dim2(); ++j) {
      const double diff = fabs(ref(i, j) - computed(i, j));
      if (diff > eps) {
        std::cout << "Error: entry (" << i << ", " << j << ") differs by "
                  << diff << "\n";
        correct = false;
      }
    }
  }
  return correct;
}

struct JacobiParameters {
  int n;
  int maxNumIter;
  double eps;
};

/**
 * Verify the correctness of the Jacobi algorithm by iteration a matrix
 * generated from the function initialCondition() and comparing it to
 * a reference stored in a file.
 */
bool verify(JacobiParameters paramters, const std::string &referenceFilename) {
  const Matrix init = initialCondition(paramters.n);
  Matrix phi = jacobi(init, paramters.eps, paramters.maxNumIter);
  return verifyAgainstStoredReference(phi, referenceFilename);
}

/**
 * Benchmark the Jacobi algorithm for a given number of threads.
 * Number of threads is set inside the benchmark funciton.
 * @return average exeuction time in milliseconds
 */
double benchmark(int numThreads, JacobiParameters paramters) {
  const Matrix init = initialCondition(paramters.n);
  omp_set_num_threads(numThreads);

  // start time measurement
  auto start = std::chrono::system_clock::now();

  // perform benchmark
  const int maxNumExec = 10;
  Matrix phi;
  for (int nExec = 0; nExec < maxNumExec; ++nExec) {
    phi = jacobi(init, paramters.eps, paramters.maxNumIter);
    volatile double dummy = phi(0, 0);
  }

  // compute elapsed time
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double, std::milli> diff = end - start;
  const double averageTime = diff.count() / maxNumExec;

  return averageTime;
}

int main() {
  const int numThreads = 8;
  JacobiParameters parameters{.n = 1024, .maxNumIter = 1000, .eps = 1e-5};

  // Uncomment the following lines to store a matrix on the harddisk as
  // reference.
  const Matrix init = initialCondition(parameters.n);
  Matrix phi = jacobi(init, parameters.eps, parameters.maxNumIter);
  storeMatrix(phi, "ref.asc");
  std::cout << "Stored reference matrix to 'ref.asc'\n";

  // Uncomment the following lines to verify the correctness of the
  // paralleliztion
  // std::cout << "Verifying correct result: ";
  // bool isCorrect = verify(parameters, "ref.asc");
  // if (isCorrect) {
  //   std::cout << "OK, parallel result is unchanged\n";
  // } else {
  //   std::cout << "Failed, parallelization changed the result!\n";
  //   return 1;
  // }

  // // Perform the benchmark
  // std::cout << "Starting benchmark\n";
  // const double time = benchmark(numThreads, n, numIter);
  // std::cout << "Benchmark results (" << numThreads << "): " << time <<
  // "ms\n";
}
