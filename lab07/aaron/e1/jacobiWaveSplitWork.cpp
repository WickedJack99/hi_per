/**
 *
 * Schritt 1
 * -------------
 * |T0|--|--|--|
 * -------------
 * |--|--|--|--|
 * -------------
 * |--|--|--|--|
 * -------------
 * |--|--|--|--|
 * -------------
 *
 * Schritt 2
 * -------------
 * |xx|T0|--|--|
 * -------------
 * |T1|--|--|--|
 * -------------
 * |--|--|--|--|
 * -------------
 * |--|--|--|--|
 * -------------
 *
 * Schritt 3
 * -------------
 * |xx|xx|T0|--|
 * -------------
 * |xx|T1|--|--|
 * -------------
 * |T2|--|--|--|
 * -------------
 * |--|--|--|--|
 * -------------
 *
 * Schritt 4
 * -------------
 * |xx|xx|xx|T0|
 * -------------
 * |xx|xx|T1|--|
 * -------------
 * |xx|T2|--|--|
 * -------------
 * |T3|--|--|--|
 * -------------
 */

#include <benchmark/benchmark.h>
#include "matrix.h"
#include <omp.h>
#include <atomic>
#include <vector>

void gauss_seidel(Matrix &phi, int maxNumIter)
{
  const int m = phi.dim1();
  const int n = phi.dim2();
  const double osth = 1. / 4;
  for (int iIter = 0; iIter < maxNumIter; ++iIter)
  {
    // Shared data per iteration
    std::vector<std::atomic<int>> column(m);
    for (int i = 0; i < m; ++i)
    {
      column[i].store(1);
    }
    std::atomic<int> threadsCount(0);

#pragma omp parallel for schedule(static, 10)
    for (int rowToCalculate = 1; rowToCalculate < (n - 1); rowToCalculate++)
    {
      int row = rowToCalculate;
      while (column[row] < n - 1)
      {
        // All threads beneath T_row1 have to check specific circumstances, where
        // T_row1 has no condition to wait for
        if (row != 1)
        {
          // T_rowx has to wait at least until T_row(x-1) has calculated values from
          // its last row
          // T_rowx will calculate value x at row,column if T_row(x-1) has calculated
          // its value at row-1, column-1
          while (column[row] == column[row - 1])
          {
            // Wait for wave (thread) above
          }
        }

        // Central jacobi calculation
        phi(row, column[row]) = osth * (phi(row + 1, column[row]) + phi(row - 1, column[row]) + phi(row, column[row] + 1) +
                                        phi(row, column[row] - 1));

        // Increment column index
        column[row]++;
      }
    }
  }
}

void gauss_seidel_lin(Matrix &phi, int maxNumIter)
{
  const int m = phi.dim1();
  const int n = phi.dim2();

  const double osth = 1. / 4;

  for (int iter = 0; iter < maxNumIter; ++iter)
  {
    for (int i = 1; i < m - 1; ++i)
    {
      for (int j = 1; j < n - 1; ++j)
      {
        phi(i, j) = osth * (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) +
                            phi(i, j - 1));
      }
    }
  }
}

void fill_matrix(Matrix &matrix, const int filler)
{
  const int m = matrix.dim1();
  const int n = matrix.dim2();

  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      if (i == 0 || i == m - 1 || j == 0 || j == n - 1)
      {
        matrix(i, j) = filler;
      }
      else
        matrix(i, j) = 0;
    }
  }
}

void check(const Matrix &a, const Matrix &b)
{
  const int m = a.dim1();
  const int n = a.dim2();

  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      if (a(i, j) != b(i, j))
      {
        std::cout << "Not equal at (" << i << ", " << j << "), a: " << a(i, j)
                  << " != " << b(i, j) << " :b" << std::endl;
      }
    }
  }
}

void print_matrix(const Matrix &matrix)
{
  const int m = matrix.dim1();
  const int n = matrix.dim2();

  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      std::cout << matrix(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

void benchmarkComputeResult(benchmark::State& state) {
  int iterations = state.range(0);
  Matrix sec = Matrix(1000, 1000);
  fill_matrix(sec, 10);

  for (auto _ : state) {
    gauss_seidel(sec, iterations);
    benchmark::DoNotOptimize(sec);
  }
}

int main(int argc, char** argv) {
  ::benchmark::Initialize(&argc, argv);

  for (int iterations = 0; iterations < 10000; iterations++) {
    benchmark::RegisterBenchmark("idk", benchmarkComputeResult)
        ->Arg(iterations)
        ->Unit(benchmark::kMillisecond);
  }

  ::benchmark::RunSpecifiedBenchmarks();

  return 0;
}

// int main()
// {
//   Matrix first = Matrix(8, 8);
//   Matrix sec = Matrix(8, 8);

//   fill_matrix(first, 10);
//   fill_matrix(sec, 10);

//   print_matrix(first);
//   print_matrix(sec);

//   gauss_seidel_lin(first, 1000);
//   gauss_seidel(sec, 1000);

//   check(first, sec);

//   print_matrix(first);
//   print_matrix(sec);

//   return 0;
// }