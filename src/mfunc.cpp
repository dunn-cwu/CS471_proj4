#include "mfunc.h"
#include <cmath>

double mfunc::schwefel(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += (-1 * v[i]) * sin(sqrt(fabs(v[i])));
    }

    return (418.9829 * n) - f;
}