#define _USE_MATH_DEFINES

#include "mfunc.h"
#include <cmath>

inline double nthroot(double x, double n)
{
    return pow(x, 1.0 / n);
}

// Function 1
double mfunc::schwefel(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += (-1.0 * v[i]) * sin(sqrt(fabs(v[i])));
    }

    return (418.9829 * n) - f;
}

// Function 2
double mfunc::dejong(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += v[i] * v[i];
    }

    return f;
}

// Function 3
double mfunc::rosenbrok(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = ((v[i] * v[i]) - v[i+1]);
        double b = (1.0 - v[i]);
        f += 100.0 * a * a;
        f += b * b;
    }

    return f;
}

// Function 4
double mfunc::rastrigin(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += (v[i] * v[i]) - (10.0 * cos(2.0 * M_PI * v[i]));
    }

    return 10.0 * n * f;
}

// Function 5
double mfunc::griewangk(double* v, size_t n)
{
    double sum = 0.0;
    double product = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        sum += (v[i] * v[i]) / 4000.0;
    }

    for (size_t i = 0; i < n; ++i)
    {
        product *= cos(v[i] / sqrt(i + 1.0));
    }

    return 1.0 + sum - product;
}

// Function 6
double mfunc::sineEnvelopeSineWave(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = sin(v[i]*v[i] + v[i+1]*v[i+1] - 0.5);
        a *= a;
        double b = (1 + 0.001*(v[i]*v[i] + v[i+1]*v[i+1]));
        b *= b;
        f += 0.5 + (a / b);
    }

    return -1.0 * f;
}

// Function 7
double mfunc::stretchedVSineWave(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = nthroot(v[i]*v[i] + v[i+1]*v[i+1], 4.0);
        double b = sin(50.0 * nthroot(v[i]*v[i] + v[i+1]*v[i+1], 10.0));
        b *= b;
        f += a * b + 1.0;
    }

    return f;
}

// Function 8
double mfunc::ackleysOne(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = (1.0 / pow(M_E, 0.2)) * sqrt(v[i]*v[i] + v[i+1]*v[i+1]);
        double b = 3.0 * (cos(2.0*v[i]) + sin(2.0*v[i+1]));
        f += a + b;
    }

    return f;
}

// Function 9
double mfunc::ackleysTwo(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = 20.0 / pow(M_E, 0.2 * sqrt((v[i]*v[i] + v[i+1]*v[i+1]) / 2.0));
        double b = pow(M_E, 0.5 * (cos(2.0 * M_PI * v[i]) + cos(2.0 * M_PI * v[i+1])));
        f += 20.0 + M_E - a - b;
    }

    return f;
}

// Function 10
double mfunc::eggHolder(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = -1.0 * v[i] * sin(sqrt(fabs(v[i] - v[i+1] - 47.0)));
        double b = (v[i+1] + 47) * sin(sqrt(fabs(v[i+1] + 47.0 + (v[i]/2.0))));
        f += a - b;
    }

    return f;
}

// Function 11
double mfunc::rana(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = v[i] * sin(sqrt(fabs(v[i+1] - v[i] + 1.0))) * cos(sqrt(fabs(v[i+1] + v[i] + 1.0)));
        double b = (v[i+1] + 1.0) * cos(sqrt(fabs(v[i+1] - v[i] + 1.0))) * sin(sqrt(fabs(v[i+1] + v[i] + 1.0)));
        f += a + b;
    }

    return f;
}

// Function 12
double mfunc::pathological(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = sin(sqrt(100.0*v[i]*v[i] + v[i+1]*v[i+1]));
        a = (a*a) - 0.5;
        double b = (v[i]*v[i] - 2*v[i]*v[i+1] + v[i+1]*v[i+1]);
        b = 1.0 + 0.001*b*b;
        f += 0.5 + (a/b);
    }

    return f;
}

// Function 13
double mfunc::michalewicz(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += sin(v[i]) * pow(sin(((i+1) * v[i] * v[i]) / M_PI), 20);
    }

    return -1.0 * f;
}

// Function 14
double mfunc::mastersCosineWave(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = pow(M_E, (-1.0/8.0)*(v[i]*v[i] + v[i+1]*v[i+1] + 0.5*v[i+1]*v[i]));
        double b = cos(4 * sqrt(v[i]*v[i] + v[i+1]*v[i+1] + 0.5*v[i]*v[i+1]));
        f += a * b;
    }

    return -1.0 * f;
}

// Function 15
double mfunc::quartic(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += (i+1) * v[i] * v[i] * v[i] * v[i];
    }

    return f;
}

// Helper function for mfunc::levy()
inline double w(double x)
{
    return 1.0 + (x - 1.0) / 4.0;
}

// Function 16
double mfunc::levy(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        double a = w(v[i]) - 1.0;
        a *= a;
        double b = sin(M_PI * w(v[i]) + 1.0);
        b *= b;
        double c = w(v[n - 1]) - 1.0;
        c *= c;
        double d = sin(2.0 * M_PI * w(v[n - 1]));
        d *= d;
        f += a * (1.0 + 10.0 * b) + c * (1.0 + d);
    }

    double e = sin(M_PI * w(v[0]));
    return e*e + f;
}

// Function 17
double mfunc::step(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double a = fabs(v[i]) + 0.5;
        f += a * a;
    }

    return f;
}

// Function 18
double mfunc::alpine(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        f += fabs(v[i] * sin(v[i]) + 0.1*v[i]);
    }

    return f;
}
