/* matrix.h, an extrmely simple matrix class.
 * Version 2.1
 * Copyright (C) 2022-2025 Tobias Kreilos, Offenburg University of Applied
 * Sciences
 *
 * Licensed under the Apache License, Version 2.0(the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * Matrix class
 *
 * This class implements a matrix of size m x n. The matrix is stored in
 * a one-dimensional array of size m*n. The matrix is stored in column-major
 * order, i.e. the first column is stored first, then the second column, etc.
 *
 * Static functions are provided to create a matrix of zeros or an uninitialized
 * matrix. The uninitialized matrix is not initialized, i.e. the entries are
 * not set to zero. This is useful for performance reasons, e.g. to control
 * the placement of matrix entries in locality-domain memory.
 *
 *
 */
class Matrix {
public:
  /**
   * Create a matrix of size m x n and initialize all entries to zero.
   * @param rows number of rows
   * @param cols number of columns
   */
  static Matrix zeros(int numRows, int numCols);

  /**
   * Create a square matrix of size n x n and initialize all entries to zero.
   * @param n number of rows and columns
   */
  static Matrix zeros(int n);

  /**
   * Create a matrix of size m x n and initialize all entries to zero.
   * @param dim number of rows and columns
   */
  static Matrix zeros(std::pair<int, int> dim);

  /**
   * Create a matrix of size m x n and do not initialize the entries.
   * @param rows number of rows
   * @param cols number of columns
   */
  static Matrix uninit(int m, int n);

  /**
   * Create a square matrix of size n x n and do not initialize the entries.
   * @param n number of rows and columns
   */
  static Matrix uninit(int n);

  /**
   * Create a matrix of size m x n and do not initialize the entries.
   * @param dim number of rows and columns
   */
  static Matrix uninit(std::pair<int, int> dim);

  Matrix(const Matrix &other);
  Matrix(Matrix &&other);
  ~Matrix();
  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &&other);

  // Access element a_ij of the matrix
  double &operator()(int i, int j);
  const double &operator()(int i, int j) const;

  // Obtain a pointer to the underlying data
  double *data();
  const double *data() const;

  // Getter functions for the dimensions
  std::pair<int, int> dim() const;
  int rows() const;
  int cols() const;
  int numEntries() const;

  // Comparison operators
  bool operator==(const Matrix &b);
  bool operator!=(const Matrix &b);

  // addition
  Matrix &operator+=(const Matrix &b);

  // subtraction
  Matrix &operator-=(const Matrix &b);

  // scalar multiplication
  Matrix &operator*=(double x);

  // scalar division
  Matrix &operator/=(double x);

  std::vector<double> get_row(int i) {
    std::vector<double> row(0);
    int start = this->cols() * i;
    int end = this->cols() * i + cols();

    for (int j = start; j < end; j++) {
      row.push_back(this->data()[j]);
    }
    return row;
  }

private:
  // Constructor is private to prevent creating an uninitialized matrix
  // accidentally. Use Matrix::zeros() or Matrix::uninit() instead
  Matrix(int m, int n);

  int numRows_;  // number of rows
  int numCols_;  // number of columns
  double *data_; // the matrix' entries
};

/**
 * Vector class
 *
 * This class implements a vector of size n. The vector is stored in a
 * Matrix of size n x 1.
 *
 * Constructors are provided to create a vector of zeros or an uninitialized
 * vector. The uninitialized vector is not initialized, i.e. the entries are
 * not set to zero. This is useful for performance reasons, e.g. to control
 * the placement of vector entries in locality-domain memory.
 */
class Vector {
public:
  static Vector zeros(int n) {
    Vector vec;
    vec.data_ = Matrix::zeros(n, 1);
    return vec;
  }

  static Vector uninit(int n) {
    Vector vec;
    vec.data_ = Matrix::uninit(n, 1);
    return vec;
  }

  bool operator==(const Vector &b) { return data_ == b.data_; }
  bool operator!=(const Vector &b) { return !operator==(b); }
  double &operator()(int i) { return data_(i, 0); }
  const double &operator()(int i) const { return data_(i, 0); }
  double *data() { return data_.data(); }
  const double *data() const { return data_.data(); }

  Vector operator+=(const Vector &b) {
    data_ += b.data_;
    return *this;
  }

  Vector operator-=(const Vector &b) {
    data_ -= b.data_;
    return *this;
  }
  Vector operator*=(double x) {
    data_ *= x;
    return *this;
  }
  Vector operator/=(double x) {
    data_ /= x;
    return *this;
  }

  int size() const { return data_.rows(); }

private:
  Matrix data_ = Matrix::zeros(0, 0);
};

/********** Implementation below ********************/

inline Matrix Matrix::zeros(int m, int n) {
  Matrix mat(m, n);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      mat(i, j) = 0;
    }
  }
  return mat;
}

inline Matrix Matrix::zeros(int n) { return zeros(n, n); }

inline Matrix Matrix::zeros(std::pair<int, int> dim) {
  return zeros(dim.first, dim.second);
}

inline Matrix Matrix::uninit(int m, int n) { return Matrix(m, n); }

inline Matrix Matrix::uninit(int n) { return uninit(n, n); }

inline Matrix Matrix::uninit(std::pair<int, int> dim) {
  return uninit(dim.first, dim.second);
}

inline Matrix::Matrix(const Matrix &other) {
  numRows_ = other.numRows_;
  numCols_ = other.numCols_;
  if (numRows_ == 0 || numCols_ == 0) {
    data_ = nullptr;
    return;
  }
  data_ = new double[numRows_ * numCols_];
  for (int i = 0; i < numRows_; ++i) {
    for (int j = 0; j < numCols_; ++j) {
      operator()(i, j) = other(i, j);
    }
  }
}

inline Matrix::Matrix(Matrix &&other) {
  numRows_ = other.numRows_;
  numCols_ = other.numCols_;
  data_ = other.data_;
  other.data_ = nullptr;
  other.numRows_ = 0;
  other.numCols_ = 0;
}

inline Matrix::~Matrix() {
  if (data_ != nullptr) {
    delete[] data_;
    data_ = nullptr;
  }
}

inline Matrix &Matrix::operator=(const Matrix &other) {
  if (this != &other) {
    if (data_ != nullptr) {
      delete[] data_;
    }
    numRows_ = other.numRows_;
    numCols_ = other.numCols_;
    data_ = new double[numRows_ * numCols_];
    for (int i = 0; i < numRows_; ++i) {
      for (int j = 0; j < numCols_; ++j) {
        operator()(i, j) = other(i, j);
      }
    }
  }
  return *this;
}

inline Matrix &Matrix::operator=(Matrix &&other) {
  if (this != &other) {
    if (data_ != nullptr) {
      delete[] data_;
    }
    numRows_ = other.numRows_;
    numCols_ = other.numCols_;
    data_ = other.data_;
    other.data_ = nullptr;
    other.numRows_ = 0;
    other.numCols_ = 0;
  }
  return *this;
}

inline double &Matrix::operator()(int i, int j) {
  return data_[i * numCols_ + j];
}

inline const double &Matrix::operator()(int i, int j) const {
  return data_[i * numCols_ + j];
}

inline double *Matrix::data() { return data_; }

inline const double *Matrix::data() const { return data_; }

inline std::pair<int, int> Matrix::dim() const {
  return std::pair<int, int>(numRows_, numCols_);
}

inline int Matrix::rows() const { return numRows_; }

inline int Matrix::cols() const { return numCols_; }

inline int Matrix::numEntries() const { return numRows_ * numCols_; }

inline bool Matrix::operator==(const Matrix &b) {
  const double eps = 1e-12;
  if (numRows_ != b.numRows_ || numCols_ != b.numCols_) {
    return false;
  }
  for (int i = 0; i < numRows_; ++i) {
    for (int j = 0; j < numCols_; ++j) {
      if (fabs(operator()(i, j) - b(i, j)) > eps) {
        return false;
      }
    }
  }
  return true;
}

inline bool Matrix::operator!=(const Matrix &b) { return !operator==(b); }

inline Matrix &Matrix::operator+=(const Matrix &b) {
  for (int i = 0; i < numRows_; ++i) {
    for (int j = 0; j < numCols_; ++j) {
      operator()(i, j) += b(i, j);
    }
  }
  return *this;
}

inline Matrix &Matrix::operator-=(const Matrix &b) {
  for (int i = 0; i < numRows_; ++i) {
    for (int j = 0; j < numCols_; ++j) {
      operator()(i, j) -= b(i, j);
    }
  }
  return *this;
}

inline Matrix &Matrix::operator*=(double x) {
  for (int i = 0; i < numRows_; ++i) {
    for (int j = 0; j < numCols_; ++j) {
      operator()(i, j) *= x;
    }
  }
  return *this;
}

inline Matrix &Matrix::operator/=(double x) {
  for (int i = 0; i < numRows_; ++i) {
    for (int j = 0; j < numCols_; ++j) {
      operator()(i, j) /= x;
    }
  }
  return *this;
}

inline Matrix::Matrix(int m, int n) : numRows_(m), numCols_(n) {
  data_ = new double[numRows_ * numCols_];
  if (data_ == nullptr) {
    throw std::bad_alloc();
  }
}

inline std::ostream &operator<<(std::ostream &os, const Matrix &a) {
  const int width = 10;
  const int precision = 4;

  const auto originalPrecision = os.precision();
  os << std::setprecision(precision);

  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      os << std::setw(width) << a(i, j) << " ";
    }
    if (i != a.rows() - 1)
      os << "\n";
  }

  os << std::setprecision(originalPrecision);
  return os;
}

inline bool equalWithinRange(const Matrix &a, const Matrix &b,
                             double eps = 1e-12) {
  if (a.rows() != b.rows() || a.cols() != b.cols())
    return false;

  int m = a.rows();
  int n = a.cols();
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (fabs(a(i, j) - b(i, j)) > eps) {
        return false;
      }
    }
  }

  return true;
}

#endif // MATRIX_H
