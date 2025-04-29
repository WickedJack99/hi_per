/**
 * matrix.h a very simplistic class for m times n matrices.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

// A very simplistic vector class for vectors of size n
class Vector {
 public:
  // constructors
  Vector(int n) : n_(n), data_(n_, 0) {}
  Vector(const Vector& other) = default;
  Vector(Vector&& other) = default;
  ~Vector() = default;

  // assignment operators
  Vector& operator=(const Vector& other) = default;
  Vector& operator=(Vector&& other) = default;

  // element access
  double& operator()(int i) { return data_[i]; }
  const double& operator()(int i) const { return data_[i]; }

  // getter functions for the dimensions
  int dim() const { return n_; }

  // comparison operators
  bool operator==(const Vector& b) { return (data_ == b.data_); }
  bool operator!=(const Vector& b) { return (data_ != b.data_); }

  // addition
  Vector& operator+=(const Vector& b) {
    for (int i = 0; i < n_; ++i) {
      operator()(i) += b(i);
    }
    return *this;
  }

  // subtraction
  Vector& operator-=(const Vector& b) {
    for (int i = 0; i < n_; ++i) {
      operator()(i) -= b(i);
    }
    return *this;
  }

  // scalar multiplication
  Vector& operator*=(double x) {
    for (int i = 0; i < n_; ++i) {
      operator()(i) *= x;
    }
    return *this;
  }

  // dot product between two vectors
  double dot(const Vector& other) const {
    double sum = 0;
    for (int i = 0; i < n_; ++i) {
      sum += operator()(i) * other(i);
    }
    return sum;
  }

 private:
  int n_;                     // vector dimension
  std::vector<double> data_;  // the vectors entries
};

inline double dot(const Vector& v1, const Vector& v2) {
  return v1.dot(v2);
}

// Print the vector as a table
inline std::ostream& operator<<(std::ostream& os, const Vector& a) {
  const int width = 10;
  const int precision = 4;

  const auto originalPrecision = os.precision();
  os << std::setprecision(precision);

  for (int i = 0; i < a.dim(); ++i) {
    os << std::setw(width) << a(i) << " ";
  }

  os << "\n";

  os << std::setprecision(originalPrecision);
  return os;
}

// A very simple class for m times n matrices
class Matrix {
 public:
  // constructors
  Matrix() : Matrix(0, 0) {}
  Matrix(int m, int n) : m_(m), n_(n), data_(m_ * n_, 0) {}
  Matrix(std::pair<int, int> dim) : Matrix(dim.first, dim.second) {}
  Matrix(int n) : Matrix(n, n) {}
  Matrix(const Matrix& other) = default;
  Matrix(Matrix&& other) = default;
  ~Matrix() = default;

  // assignment operators
  Matrix& operator=(const Matrix& other) = default;
  Matrix& operator=(Matrix&& other) = default;

  // element access
  double& operator()(int i, int j) { return data_[i * n_ + j]; }
  const double& operator()(int i, int j) const { return data_[i * n_ + j]; }

  // getter functions for the dimensions
  std::pair<int, int> dim() const { return std::pair<int, int>(m_, n_); }
  int dim1() const { return m_; }
  int dim2() const { return n_; }
  int numEntries() const { return data_.size(); }

  // comparison operators
  bool operator==(const Matrix& b) { return (data_ == b.data_); }
  bool operator!=(const Matrix& b) { return (data_ != b.data_); }

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
  int m_;                     // first dimension
  int n_;                     // second dimension
  std::vector<double> data_;  // the matrix' entries
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

// A very simple class for "3D-Matrices" (tensors) with dimension l x m x n
class Matrix3D {
 public:
  // constructors
  Matrix3D(int l, int m, int n) : l_(l), m_(m), n_(n), data_(l) {
    for (int i = 0; i < l_; ++i) {
      data_[i] = std::vector<std::vector<double>>(m_);
      for (int j = 0; j < m_; ++j) {
        data_[i][j] = std::vector<double>(n_, 0);
      }
    }
  }
  Matrix3D(int n) : Matrix3D(n, n, n) {}
  Matrix3D(const Matrix3D& other) = default;
  Matrix3D(Matrix3D&& other) = default;
  ~Matrix3D() = default;

  // assignment operators
  Matrix3D& operator=(const Matrix3D& other) = default;
  Matrix3D& operator=(Matrix3D&& other) = default;

  // element access
  double& operator()(int i, int j, int k) { return data_[i][j][k]; }
  const double& operator()(int i, int j, int k) const { return data_[i][j][k]; }

  // getter functions for the dimensions
  int dim1() const { return l_; }
  int dim2() const { return m_; }
  int dim3() const { return n_; }

  // comparison operators
  bool operator==(const Matrix3D& b) { return (data_ == b.data_); }
  bool operator!=(const Matrix3D& b) { return (data_ != b.data_); }

  // addition
  Matrix3D& operator+=(const Matrix3D& b) {
    for (int i = 0; i < l_; ++i) {
      for (int j = 0; j < m_; ++j) {
        for (int k = 0; k < n_; ++k) {
          operator()(i, j, k) += b(i, j, k);
        }
      }
    }
    return *this;
  }

  // substraction
  Matrix3D& operator-=(const Matrix3D& b) {
    for (int i = 0; i < l_; ++i) {
      for (int j = 0; j < m_; ++j) {
        for (int k = 0; k < n_; ++k) {
          operator()(i, j, k) -= b(i, j, k);
        }
      }
    }
    return *this;
  }

  // scalar multiplication
  Matrix3D& operator*=(double x) {
    for (int i = 0; i < l_; ++i) {
      for (int j = 0; j < m_; ++j) {
        for (int k = 0; k < n_; ++k) {
          operator()(i, j, k) *= x;
        }
      }
    }
    return *this;
  }

  // scalar division
  Matrix3D& operator/=(double x) {
    for (int i = 0; i < l_; ++i) {
      for (int j = 0; j < m_; ++j) {
        for (int k = 0; k < n_; ++k) {
          operator()(i, j, k) /= x;
        }
      }
    }
    return *this;
  }

 private:
  int l_;                                               // first dimension
  int m_;                                               // second dimension
  int n_;                                               // third dimension
  std::vector<std::vector<std::vector<double>>> data_;  // the tensors' entries
};

#endif  // MATRIX_H