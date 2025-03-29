#include "matrix.h"
#include <cmath>
#include <iostream>
#include <vector>

int fib2(int i) {
if (i == 0) {
return 0;
} else if (i == 1) {
return 1;
} else {
return fib2(i - 1) + fib2(i - 2);
}
}
Matrix compute2(const Matrix& s, const std::vector<int>& v) {
Matrix m(s.dim());
const int n = v.size();
for (int j = 0; j < n; ++j) {
for (int i = 0; i < n; ++i) {
double val = static_cast<double>(fib2(v[j] % 256));
m(j, i) = s(j, i) * (sin(val) * tan(val) / sqrt(cos(val) + 2));
}
}
return m;
}
