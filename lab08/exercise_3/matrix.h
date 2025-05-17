/**
 * matrix.h a very simplistic class for m times n matrices.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

// A very simple class for m times n matrices
class Matrix {
 public:
  // constructors
  Matrix() : Matrix(0, 0) {}

  Matrix(int m, int n, bool init = true) : m_(m), n_(n) {
    data_ = new double[m_ * n_];
    if (data_ == nullptr) {
      std::cout << "Error: not enough memory for matrix\n";
      exit(1);
    }
    if (init) {
      for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < n_; ++j) {
          operator()(i, j) = 0;
        }
      }
    }
  }

  Matrix(std::pair<int, int> dim, bool init = true)
      : Matrix(dim.first, dim.second, init) {}

  Matrix(int n) : Matrix(n, n) {}

  Matrix(const Matrix& other) {
    m_ = other.m_;
    n_ = other.n_;
    if (m_ == 0 || n_ == 0) {
      data_ = nullptr;
      return;
    }
    data_ = new double[m_ * n_];
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        operator()(i, j) = other(i, j);
      }
    }
  }
  Matrix(Matrix&& other) {
    m_ = other.m_;
    n_ = other.n_;
    data_ = other.data_;
    other.data_ = nullptr;
    other.m_ = 0;
    other.n_ = 0;
  }

  ~Matrix() {
    if (data_ != nullptr) {
      delete[] data_;
      data_ = nullptr;
    }
  }

  // assignment operators
  Matrix& operator=(const Matrix& other) {
    if (this != &other) {
      if (data_ != nullptr) {
        delete[] data_;
      }
      m_ = other.m_;
      n_ = other.n_;
      data_ = new double[m_ * n_];
      for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < n_; ++j) {
          operator()(i, j) = other(i, j);
        }
      }
    }
    return *this;
  }

  Matrix& operator=(Matrix&& other) {
    if (this != &other) {
      if (data_ != nullptr) {
        delete[] data_;
      }
      m_ = other.m_;
      n_ = other.n_;
      data_ = other.data_;
      other.data_ = nullptr;
      other.m_ = 0;
      other.n_ = 0;
    }
    return *this;
  }

  // element access
  double& operator()(int i, int j) { return data_[i * n_ + j]; }
  const double& operator()(int i, int j) const { return data_[i * n_ + j]; }

  // getter functions for the dimensions
  std::pair<int, int> dim() const { return std::pair<int, int>(m_, n_); }
  int dim1() const { return m_; }
  int dim2() const { return n_; }
  int numEntries() const { return m_ * n_; }

  // comparison operators
  bool operator==(const Matrix& b) {
    const double eps = 1e-12;
    if (m_ != b.m_ || n_ != b.n_) {
      return false;
    }
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        if (fabs(operator()(i, j) - b(i, j)) > eps) {
          return false;
        }
      }
    }
    return true;
  }
  bool operator!=(const Matrix& b) { return !operator==(b); }
  // comparison operators

  // addition
  Matrix& operator+=(const Matrix& b) {
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        operator()(i, j) += b(i, j);
      }
    }
    return *this;
  }

  // subtraction
  Matrix& operator-=(const Matrix& b) {
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        operator()(i, j) -= b(i, j);
      }
    }
    return *this;
  }

  // scalar multiplication
  Matrix& operator*=(double x) {
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        operator()(i, j) *= x;
      }
    }
    return *this;
  }

  // scalar division
  Matrix& operator/=(double x) {
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        operator()(i, j) /= x;
      }
    }
    return *this;
  }

  // matrix product (only for square matrices of equal dimension)
  Matrix& operator*=(const Matrix& b) {
    if (dim1() != dim2()) {
      std::cout << "Error in matrix multiplication: no square matrix\n";
    } else if (dim1() != b.dim1() || dim2() != b.dim2()) {
      std::cout << "Error in matrix multiplication: dimensions do not match\n";
    } else {
      Matrix a = *this;
      Matrix& c = *this;
      const int m = dim1();
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          for (int k = 0; k < m; ++k) {
            c(i, j) += a(i, k) * b(k, j);
          }
        }
      }
    }

    return *this;
  }

 public:
  int m_;         // first dimension
  int n_;         // second dimension
  double* data_;  // the matrix' entries
};

// Print the matrix as a table
inline std::ostream& operator<<(std::ostream& os, const Matrix& a) {
  const int width = 10;
  const int precision = 4;

  const auto originalPrecision = os.precision();
  os << std::setprecision(precision);

  for (int i = 0; i < a.dim1(); ++i) {
    for (int j = 0; j < a.dim2(); ++j) {
      os << std::setw(width) << a(i, j) << " ";
    }
    if (i != a.dim1() - 1)
      os << "\n";
  }

  os << std::setprecision(originalPrecision);
  return os;
}

// matrix product
inline Matrix operator*(const Matrix& a, const Matrix& b) {
  if (a.dim2() == b.dim1()) {
    int m = a.dim1();
    int n = a.dim2();
    int p = b.dim2();
    Matrix c(m, p);
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < n; ++k) {
          c(i, j) += a(i, k) * b(k, j);
        }
      }
    }
    return c;
  } else {
    return Matrix(0, 0);
  }
}

inline bool equalWithinRange(const Matrix& a,
                             const Matrix& b,
                             double eps = 1e-12) {
  if (a.dim1() != b.dim1() || a.dim2() != b.dim2())
    return false;

  int m = a.dim1();
  int n = a.dim2();
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (fabs(a(i, j) - b(i, j)) > eps) {
        return false;
      }
    }
  }

  return true;
}

#endif  // MATRIX_H