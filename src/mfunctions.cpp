/**
 * @file mfunctions.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementations for various math functions defined in mfunc.h
 * @version 0.1
 * @date 2019-03-29
 * 
 * @copyright Copyright (c) 2019
 */

#define _USE_MATH_DEFINES

#include "mfunctions.h"
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
template <class T>
T mfunc::schwefel(T* v, size_t n)
{
    T f = 0.0;

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
template <class T>
T mfunc::dejong(T* v, size_t n)
{
    T f = 0.0;

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
template <class T>
T mfunc::rosenbrok(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = ((v[i] * v[i]) - v[i+1]);
        T b = (1.0 - v[i]);
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
template <class T>
T mfunc::rastrigin(T* v, size_t n)
{
    T f = 0.0;

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
template <class T>
T mfunc::griewangk(T* v, size_t n)
{
    T sum = 0.0;
    T product = 0.0;

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
template <class T>
T mfunc::sineEnvelopeSineWave(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = sin(v[i]*v[i] + v[i+1]*v[i+1] - 0.5);
        a *= a;
        T b = (1 + 0.001*(v[i]*v[i] + v[i+1]*v[i+1]));
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
template <class T>
T mfunc::stretchedVSineWave(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = nthroot(v[i]*v[i] + v[i+1]*v[i+1], 4.0);
        T b = sin(50.0 * nthroot(v[i]*v[i] + v[i+1]*v[i+1], 10.0));
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
template <class T>
T mfunc::ackleysOne(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = (1.0 / pow(M_E, 0.2)) * sqrt(v[i]*v[i] + v[i+1]*v[i+1]);
        T b = 3.0 * (cos(2.0*v[i]) + sin(2.0*v[i+1]));
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
template <class T>
T mfunc::ackleysTwo(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = 20.0 / pow(M_E, 0.2 * sqrt((v[i]*v[i] + v[i+1]*v[i+1]) / 2.0));
        T b = pow(M_E, 0.5 * (cos(2.0 * M_PI * v[i]) + cos(2.0 * M_PI * v[i+1])));
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
template <class T>
T mfunc::eggHolder(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = -1.0 * v[i] * sin(sqrt(fabs(v[i] - v[i+1] - 47.0)));
        T b = (v[i+1] + 47) * sin(sqrt(fabs(v[i+1] + 47.0 + (v[i]/2.0))));
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
template <class T>
T mfunc::rana(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = v[i] * sin(sqrt(fabs(v[i+1] - v[i] + 1.0))) * cos(sqrt(fabs(v[i+1] + v[i] + 1.0)));
        T b = (v[i+1] + 1.0) * cos(sqrt(fabs(v[i+1] - v[i] + 1.0))) * sin(sqrt(fabs(v[i+1] + v[i] + 1.0)));
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
template <class T>
T mfunc::pathological(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = sin(sqrt(100.0*v[i]*v[i] + v[i+1]*v[i+1]));
        a = (a*a) - 0.5;
        T b = (v[i]*v[i] - 2*v[i]*v[i+1] + v[i+1]*v[i+1]);
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
template <class T>
T mfunc::michalewicz(T* v, size_t n)
{
    T f = 0.0;

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
template <class T>
T mfunc::mastersCosineWave(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = pow(M_E, (-1.0/8.0)*(v[i]*v[i] + v[i+1]*v[i+1] + 0.5*v[i+1]*v[i]));
        T b = cos(4 * sqrt(v[i]*v[i] + v[i+1]*v[i+1] + 0.5*v[i]*v[i+1]));
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
template <class T>
T mfunc::quartic(T* v, size_t n)
{
    T f = 0.0;

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
template <class T>
inline T w(T x)
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
template <class T>
T mfunc::levy(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = w<T>(v[i]) - 1.0;
        a *= a;
        T b = sin(M_PI * w<T>(v[i]) + 1.0);
        b *= b;
        T c = w<T>(v[n - 1]) - 1.0;
        c *= c;
        T d = sin(2.0 * M_PI * w<T>(v[n - 1]));
        d *= d;
        f += a * (1.0 + 10.0 * b) + c * (1.0 + d);
    }

    T e = sin(M_PI * w<T>(v[0]));
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
template <class T>
T mfunc::step(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        T a = fabs(v[i]) + 0.5;
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
template <class T>
T mfunc::alpine(T* v, size_t n)
{
    T f = 0.0;

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
template <class T>
bool mfunc::fExec(unsigned int f, T* v, size_t n, T& outResult)
{
    switch (f)
    {
        case 1:
            outResult = schwefel<T>(v, n);
            return true;
        case 2:
            outResult = dejong<T>(v, n);
            return true;
        case 3:
            outResult = rosenbrok<T>(v, n);
            return true;
        case 4:
            outResult = rastrigin<T>(v, n);
            return true;
        case 5:
            outResult = griewangk<T>(v, n);
            return true;
        case 6:
            outResult = sineEnvelopeSineWave<T>(v, n);
            return true;
        case 7:
            outResult = stretchedVSineWave<T>(v, n);
            return true;
        case 8:
            outResult = ackleysOne<T>(v, n);
            return true;
        case 9:
            outResult = ackleysTwo<T>(v, n);
            return true;
        case 10:
            outResult = eggHolder<T>(v, n);
            return true;
        case 11:
            outResult = rana<T>(v, n);
            return true;
        case 12:
            outResult = pathological<T>(v, n);
            return true;
        case 13:
            outResult = michalewicz<T>(v, n);
            return true;
        case 14:
            outResult = mastersCosineWave<T>(v, n);
            return true;
        case 15:
            outResult = quartic<T>(v, n);
            return true;
        case 16:
            outResult = levy<T>(v, n);
            return true;
        case 17:
            outResult = step<T>(v, n);
            return true;
        case 18:
            outResult = alpine<T>(v, n);
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

template <class T>
mfunc::mfuncPtr<T> mfunc::fGet(unsigned int f)
{
    switch (f)
    {
        case 1:
            return mfunc::schwefel<T>;
        case 2:
            return mfunc::dejong<T>;
        case 3:
            return mfunc::rosenbrok<T>;
        case 4:
            return mfunc::rastrigin<T>;
        case 5:
            return mfunc::griewangk<T>;
        case 6:
            return mfunc::sineEnvelopeSineWave<T>;
        case 7:
            return mfunc::stretchedVSineWave<T>;
        case 8:
            return mfunc::ackleysOne<T>;
        case 9:
            return mfunc::ackleysTwo<T>;
        case 10:
            return mfunc::eggHolder<T>;
        case 11:
            return mfunc::rana<T>;
        case 12:
            return mfunc::pathological<T>;
        case 13:
            return mfunc::michalewicz<T>;
        case 14:
            return mfunc::mastersCosineWave<T>;
        case 15:
            return mfunc::quartic<T>;
        case 16:
            return mfunc::levy<T>;
        case 17:
            return mfunc::step<T>;
        case 18:
            return mfunc::alpine<T>;
        default:
            return nullptr;
    }
}

template bool mfunc::fExec<double>(unsigned int f, double* v, size_t n, double& outResult);
template mfunc::mfuncPtr<double> mfunc::fGet(unsigned int f);

// =========================
// End of mfunctions.cpp
// =========================
