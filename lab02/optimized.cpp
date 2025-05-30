#include "matrix.h"
#include <cmath>
int fib(int n)
{
    if (n < 2)
        return n;

    int a = 0;
    int b = 1;

    for (int i = 2; i <= n; ++i)
    {
        int temp = a + b;
        a = b;
        b = temp;
    }

    return b;
}

Matrix compute(const Matrix &s, const std::vector<int> &v)
{
    Matrix m(s.dim());
    const int n = v.size();
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            double val = static_cast<double>(fib(v[j] % 256));
            m(j, i) = s(j, i) * (sin(val) * tan(val) / sqrt(cos(val) + 2));
        }
    }
    return m;
}