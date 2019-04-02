/**
 * @file mfunc.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementations for various math functions defined in mfunc.h
 * @version 0.1
 * @date 2019-03-29
 * 
 * @copyright Copyright (c) 2019
 */

#define _USE_MATH_DEFINES

#include "mfunc.h"
#include <cmath>

#define _schwefelDesc "Schwefel’s function"
#define _dejongDesc "1st De Jong’s function"
#define _rosenbrokDesc "Rosenbrock"
#define _rastriginDesc "Rastrigin"
#define _griewangkDesc "Griewangk"
#define _sineEnvelopeSineWaveDesc "Sine Envelope Sine Wave"
#define _stretchedVSineWaveDesc "Stretched V Sine Wave"
#define _ackleysOneDesc "Ackley’s One"
#define _ackleysTwoDesc "Ackley’s Two"
#define _eggHolderDesc "Egg Holder"
#define _ranaDesc "Rana"
#define _pathologicalDesc "Pathological"
#define _michalewiczDesc "Michalewicz"
#define _mastersCosineWaveDesc "Masters Cosine Wave"
#define _quarticDesc "Quartic"
#define _levyDesc "Levy"
#define _stepDesc "Step"
#define _alpineDesc "Alpine"

/**
 * Simple inline helper function that returns the nth-root
 * @param x Value to be taken to the nth power
 * @param n root degree
 * @return The value of the nth-root of x
 */
inline double nthroot(double x, double n)
{
    return pow(x, 1.0 / n);
}

/**
 * Constant value for the total number of math functions contained in this namespace
 */
const unsigned int mfunc::NUM_FUNCTIONS = 18;

// ================================================

/**
 * Returns a string description of the schwefel() function
 * @return C-string description
 */
const char* mfunc::schwefelDesc()
{
    return _schwefelDesc;
}

/**
 * @brief Function 1.
 * Implementation of Schwefel’s mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::schwefel(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += (-1.0 * v[i]) * sin(sqrt(fabs(v[i])));
    }

    return (418.9829 * n) - f;
}

// ================================================

/**
 * Returns a string description of the dejong() function
 * @return C-string description
 */
const char* mfunc::dejongDesc()
{
    return _dejongDesc;
}

/**
 * @brief Function 2.
 * Implementation of 1st De Jong’s mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::dejong(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += v[i] * v[i];
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the rosenbrok() function
 * @return C-string description
 */
const char* mfunc::rosenbrokDesc()
{
    return _rosenbrokDesc;
}

/**
 * @brief Function 3.
 * Implementation of the Rosenbrock mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::rosenbrok(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = ((v[i] * v[i]) - v[i+1]);
        double b = (1.0 - v[i]);
        f += 100.0 * a * a;
        f += b * b;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the rastrigin() function
 * @return C-string description
 */
const char* mfunc::rastriginDesc()
{
    return _rastriginDesc;
}

/**
 * @brief Function 4.
 * Implementation of the Rastrigin mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::rastrigin(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += (v[i] * v[i]) - (10.0 * cos(2.0 * M_PI * v[i]));
    }

    return 10.0 * n * f;
}

// ================================================

/**
 * Returns a string description of the griewangk() function
 * @return C-string description
 */
const char* mfunc::griewangkDesc()
{
    return _griewangkDesc;
}

/**
 * @brief Function 5.
 * Implementation of the Griewangk mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::griewangk(double* v, size_t n)
{
    double sum = 0.0;
    double product = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        sum += (v[i] * v[i]) / 4000.0;
    }

    for (size_t i = 0; i < n; i++)
    {
        product *= cos(v[i] / sqrt(i + 1.0));
    }

    return 1.0 + sum - product;
}

// ================================================

/**
 * Returns a string description of the sineEnvelopeSineWave() function
 * @return C-string description
 */
const char* mfunc::sineEnvelopeSineWaveDesc()
{
    return _sineEnvelopeSineWaveDesc;
}

/**
 * @brief Function 6.
 * Implementation of the Sine Envelope Sine Wave mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::sineEnvelopeSineWave(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = sin(v[i]*v[i] + v[i+1]*v[i+1] - 0.5);
        a *= a;
        double b = (1 + 0.001*(v[i]*v[i] + v[i+1]*v[i+1]));
        b *= b;
        f += 0.5 + (a / b);
    }

    return -1.0 * f;
}

// ================================================

/**
 * Returns a string description of the stretchedVSineWave() function
 * @return C-string description
 */
const char* mfunc::stretchedVSineWaveDesc()
{
    return _stretchedVSineWaveDesc;
}

/**
 * @brief Function 7.
 * Implementation of the Stretched V Sine Wave mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::stretchedVSineWave(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = nthroot(v[i]*v[i] + v[i+1]*v[i+1], 4.0);
        double b = sin(50.0 * nthroot(v[i]*v[i] + v[i+1]*v[i+1], 10.0));
        b *= b;
        f += a * b + 1.0;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the ackleysOne() function
 * @return C-string description
 */
const char* mfunc::ackleysOneDesc()
{
    return _ackleysOneDesc;
}

/**
 * @brief Function 8.
 * Implementation of Ackley’s One mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::ackleysOne(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = (1.0 / pow(M_E, 0.2)) * sqrt(v[i]*v[i] + v[i+1]*v[i+1]);
        double b = 3.0 * (cos(2.0*v[i]) + sin(2.0*v[i+1]));
        f += a + b;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the ackleysTwo() function
 * @return C-string description
 */
const char* mfunc::ackleysTwoDesc()
{
    return _ackleysTwoDesc;
}

/**
 * @brief Function 9.
 * Implementation of Ackley’s Two mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::ackleysTwo(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = 20.0 / pow(M_E, 0.2 * sqrt((v[i]*v[i] + v[i+1]*v[i+1]) / 2.0));
        double b = pow(M_E, 0.5 * (cos(2.0 * M_PI * v[i]) + cos(2.0 * M_PI * v[i+1])));
        f += 20.0 + M_E - a - b;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the eggHolder() function
 * @return C-string description
 */
const char* mfunc::eggHolderDesc()
{
    return _eggHolderDesc;
}

/**
 * @brief Function 10.
 * Implementation of the Egg Holder mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::eggHolder(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = -1.0 * v[i] * sin(sqrt(fabs(v[i] - v[i+1] - 47.0)));
        double b = (v[i+1] + 47) * sin(sqrt(fabs(v[i+1] + 47.0 + (v[i]/2.0))));
        f += a - b;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the rana() function
 * @return C-string description
 */
const char* mfunc::ranaDesc()
{
    return _ranaDesc;
}

/**
 * @brief Function 11.
 * Implementation of the Rana mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::rana(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = v[i] * sin(sqrt(fabs(v[i+1] - v[i] + 1.0))) * cos(sqrt(fabs(v[i+1] + v[i] + 1.0)));
        double b = (v[i+1] + 1.0) * cos(sqrt(fabs(v[i+1] - v[i] + 1.0))) * sin(sqrt(fabs(v[i+1] + v[i] + 1.0)));
        f += a + b;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the pathological() function
 * @return C-string description
 */
const char* mfunc::pathologicalDesc()
{
    return _pathologicalDesc;
}

/**
 * @brief Function 12.
 * Implementation of the Pathological mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::pathological(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = sin(sqrt(100.0*v[i]*v[i] + v[i+1]*v[i+1]));
        a = (a*a) - 0.5;
        double b = (v[i]*v[i] - 2*v[i]*v[i+1] + v[i+1]*v[i+1]);
        b = 1.0 + 0.001 * b*b;
        f += 0.5 + (a/b);
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the michalewicz() function
 * @return C-string description
 */
const char* mfunc::michalewiczDesc()
{
    return _michalewiczDesc;
}

/**
 * @brief Function 13.
 * Implementation of the Michalewicz mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::michalewicz(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += sin(v[i]) * pow(sin(((i+1) * v[i] * v[i]) / M_PI), 20);
    }

    return -1.0 * f;
}

// ================================================

/**
 * Returns a string description of the mastersCosineWave() function
 * @return C-string description
 */
const char* mfunc::mastersCosineWaveDesc()
{
    return _mastersCosineWaveDesc;
}

/**
 * @brief Function 14.
 * Implementation of the Masters Cosine Wave mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::mastersCosineWave(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        double a = pow(M_E, (-1.0/8.0)*(v[i]*v[i] + v[i+1]*v[i+1] + 0.5*v[i+1]*v[i]));
        double b = cos(4 * sqrt(v[i]*v[i] + v[i+1]*v[i+1] + 0.5*v[i]*v[i+1]));
        f += a * b;
    }

    return -1.0 * f;
}

// ================================================

/**
 * Returns a string description of the quartic() function
 * @return C-string description
 */
const char* mfunc::quarticDesc()
{
    return _quarticDesc;
}

/**
 * @brief Function 15.
 * Implementation of the Quartic mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::quartic(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += (i+1) * v[i] * v[i] * v[i] * v[i];
    }

    return f;
}

// ================================================

/**
 * Helper math function used in levy()
 */
inline double w(double x)
{
    return 1.0 + (x - 1.0) / 4.0;
}

/**
 * Returns a string description of the levy() function
 * @return C-string description
 */
const char* mfunc::levyDesc()
{
    return _levyDesc;
}

/**
 * @brief Function 16.
 * Implementation of the Levy mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::levy(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
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

// ================================================

/**
 * Returns a string description of the step() function
 * @return C-string description
 */
const char* mfunc::stepDesc()
{
    return _stepDesc;
}

/**
 * @brief Function 17.
 * Implementation of the Step mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::step(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        double a = fabs(v[i]) + 0.5;
        f += a * a;
    }

    return f;
}

// ================================================

/**
 * Returns a string description of the alpine() function
 * @return C-string description
 */
const char* mfunc::alpineDesc()
{
    return _alpineDesc;
}

/**
 * @brief Function 18.
 * Implementation of the Alpine mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
double mfunc::alpine(double* v, size_t n)
{
    double f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += fabs(v[i] * sin(v[i]) + 0.1*v[i]);
    }

    return f;
}

// ================================================

/**
 * @brief Executes a specific function
 * Executes the function with the given id and returns true on success.
 * Otherwise returns false if id is invalid. 
 * @param f Function id to execute
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @param outResult Output reference variable for the result of the mathematical function 
 * @return true if 'f' is a valid id and the function was ran. Otherwise false.
 */
bool mfunc::fExec(unsigned int f, double* v, size_t n, double& outResult)
{
    switch (f)
    {
        case 1:
            outResult = schwefel(v, n);
            return true;
        case 2:
            outResult = dejong(v, n);
            return true;
        case 3:
            outResult = rosenbrok(v, n);
            return true;
        case 4:
            outResult = rastrigin(v, n);
            return true;
        case 5:
            outResult = griewangk(v, n);
            return true;
        case 6:
            outResult = sineEnvelopeSineWave(v, n);
            return true;
        case 7:
            outResult = stretchedVSineWave(v, n);
            return true;
        case 8:
            outResult = ackleysOne(v, n);
            return true;
        case 9:
            outResult = ackleysTwo(v, n);
            return true;
        case 10:
            outResult = eggHolder(v, n);
            return true;
        case 11:
            outResult = rana(v, n);
            return true;
        case 12:
            outResult = pathological(v, n);
            return true;
        case 13:
            outResult = michalewicz(v, n);
            return true;
        case 14:
            outResult = mastersCosineWave(v, n);
            return true;
        case 15:
            outResult = quartic(v, n);
            return true;
        case 16:
            outResult = levy(v, n);
            return true;
        case 17:
            outResult = step(v, n);
            return true;
        case 18:
            outResult = alpine(v, n);
            return true;
        default:
            return false;
    }
}

// ================================================

/**
 * @brief Returns a function's description
 * Returns a C-string description for the given function id if the id is valid.
 * Otherwise returns null
 * @param f Function id to retrieve the description for
 * @return A C-string containing the function description if id is valid, otherwise null.
 */
const char* mfunc::fDesc(unsigned int f)
{
    switch (f)
    {
        case 1:
            return schwefelDesc();
        case 2:
            return dejongDesc();
        case 3:
            return rosenbrokDesc();
        case 4:
            return rastriginDesc();
        case 5:
            return griewangkDesc();
        case 6:
            return sineEnvelopeSineWaveDesc();
        case 7:
            return stretchedVSineWaveDesc();
        case 8:
            return ackleysOneDesc();
        case 9:
            return ackleysTwoDesc();
        case 10:
            return eggHolderDesc();
        case 11:
            return ranaDesc();
        case 12:
            return pathologicalDesc();
        case 13:
            return michalewiczDesc();
        case 14:
            return mastersCosineWaveDesc();
        case 15:
            return quarticDesc();
        case 16:
            return levyDesc();
        case 17:
            return stepDesc();
        case 18:
            return alpineDesc();
        default:
            return NULL;
    }
}

// =========================
// End of mfunc.cpp
// =========================
