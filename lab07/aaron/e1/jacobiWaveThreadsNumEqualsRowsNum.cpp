// WORKS!

/**
 *         Spalte 0 1 2 3 4 5 6
 * Schritt
 *  0
 *
 *  1
 *
 *  2
 *
 *  3
 *
 * Über Schritte iterieren statt ?
 *
 * Beispiel A: jeder Thread bearbeitet eine Zeile, es gibt genau m Threads.
 * D.h. m / t = 1
 *
 * Beispiel B: es gibt eine konstante Anzahl Threads t es gibt m Zeilen.
 * Thread tx nimmt sich irgendeine Zeile und
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
    int threadNum = m - 2;
    std::vector<std::atomic<int>> columnsCalculated(threadNum);
    for (int i = 0; i < threadNum; ++i)
    {
      columnsCalculated[i].store(1);
    }
#pragma omp parallel num_threads(threadNum)
    {
      int threadId = omp_get_thread_num();
      int row = threadId + 1;
      while (columnsCalculated[threadId] < n - 1)
      {
        int currentCol = columnsCalculated[threadId];
        if (threadId != 0)
        {
          while (currentCol == columnsCalculated[threadId - 1])
          {
          }
        }
        phi(row, currentCol) = osth * (phi(row + 1, currentCol) + phi(row - 1, currentCol) + phi(row, currentCol + 1) +
                                       phi(row, currentCol - 1));
        columnsCalculated[threadId]++;
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

int main()
{
  Matrix first = Matrix(8, 8);
  Matrix sec = Matrix(8, 8);

  fill_matrix(first, 10);
  fill_matrix(sec, 10);

  print_matrix(first);
  print_matrix(sec);

  gauss_seidel_lin(first, 1);
  gauss_seidel(sec, 1);

  check(first, sec);

  print_matrix(first);
  print_matrix(sec);

  return 0;
}