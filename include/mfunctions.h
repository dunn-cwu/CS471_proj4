/**
 * @file mfunctions.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains various math function definitions
 * @version 0.1
 * @date 2019-03-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __MFUNCIONS_H
#define __MFUNCIONS_H

#define _USE_MATH_DEFINES

#include <cmath>
#include "mfuncPtr.h"

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
 * Scope for all math functions
 */
namespace mfunc
{
    /**
     * Constant value for the total number of math functions contained in this namespace
     */
    constexpr const unsigned int NUM_FUNCTIONS = 18;

    /**
    * @brief get() returns a function's description
    * Returns a C-string description for the given function id if the id is valid.
    * Otherwise returns null
    * @param f Function id to retrieve the description for
    * @return A C-string containing the function description if id is valid, otherwise null.
    */
    struct FunctionDesc
    {
        static const char* get(unsigned int f)
        {
            switch (f)
            {
                case 1:
                    return _schwefelDesc;
                case 2:
                    return _dejongDesc;
                case 3:
                    return _rosenbrokDesc;
                case 4:
                    return _rastriginDesc;
                case 5:
                    return _griewangkDesc;
                case 6:
                    return _sineEnvelopeSineWaveDesc;
                case 7:
                    return _stretchedVSineWaveDesc;
                case 8:
                    return _ackleysOneDesc;
                case 9:
                    return _ackleysTwoDesc;
                case 10:
                    return _eggHolderDesc;
                case 11:
                    return _ranaDesc;
                case 12:
                    return _pathologicalDesc;
                case 13:
                    return _michalewiczDesc;
                case 14:
                    return _mastersCosineWaveDesc;
                case 15:
                    return _quarticDesc;
                case 16:
                    return _levyDesc;
                case 17:
                    return _stepDesc;
                case 18:
                    return _alpineDesc;
                default:
                    return NULL;
            }
        }
    };

    template <class T>
    struct Functions
    {
        static T schwefel(T* v, size_t n);
        static T dejong(T* v, size_t n);
        static T rosenbrok(T* v, size_t n);
        static T rastrigin(T* v, size_t n);
        static T griewangk(T* v, size_t n);
        static T sineEnvelopeSineWave(T* v, size_t n);
        static T stretchedVSineWave(T* v, size_t n);
        static T ackleysOne(T* v, size_t n);
        static T ackleysTwo(T* v, size_t n);
        static T eggHolder(T* v, size_t n);
        static T rana(T* v, size_t n);
        static T pathological(T* v, size_t n);
        static T mastersCosineWave(T* v, size_t n);
        static T michalewicz(T* v, size_t n);
        static T quartic(T* v, size_t n);
        static T levy(T* v, size_t n);
        static T step(T* v, size_t n);
        static T alpine(T* v, size_t n);
        static mfuncPtr<T> get(unsigned int f);
        static bool exec(unsigned int f, T* v, size_t n, T& outResult);
        static double nthroot(double x, double n);
        static T w(T x);
    };
}

/**
 * Simple helper function that returns the nth-root
 * @param x Value to be taken to the nth power
 * @param n root degree
 * @return The value of the nth-root of x
 */
template <class T>
double mfunc::Functions<T>::nthroot(double x, double n)
{
    return pow(x, 1.0 / n);
}

// ================================================

/**
 * @brief Function 1.
 * Implementation of Schwefel’s mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::schwefel(T* v, size_t n)
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
 * @brief Function 2.
 * Implementation of 1st De Jong’s mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::dejong(T* v, size_t n)
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
 * @brief Function 3.
 * Implementation of the Rosenbrock mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::rosenbrok(T* v, size_t n)
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
 * @brief Function 4.
 * Implementation of the Rastrigin mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::rastrigin(T* v, size_t n)
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
 * @brief Function 5.
 * Implementation of the Griewangk mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::griewangk(T* v, size_t n)
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
 * @brief Function 6.
 * Implementation of the Sine Envelope Sine Wave mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::sineEnvelopeSineWave(T* v, size_t n)
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
 * @brief Function 7.
 * Implementation of the Stretched V Sine Wave mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::stretchedVSineWave(T* v, size_t n)
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
 * @brief Function 8.
 * Implementation of Ackley’s One mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::ackleysOne(T* v, size_t n)
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
 * @brief Function 9.
 * Implementation of Ackley’s Two mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::ackleysTwo(T* v, size_t n)
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
 * @brief Function 10.
 * Implementation of the Egg Holder mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::eggHolder(T* v, size_t n)
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
 * @brief Function 11.
 * Implementation of the Rana mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::rana(T* v, size_t n)
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
 * @brief Function 12.
 * Implementation of the Pathological mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::pathological(T* v, size_t n)
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
 * @brief Function 13.
 * Implementation of the Michalewicz mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::michalewicz(T* v, size_t n)
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
 * @brief Function 14.
 * Implementation of the Masters Cosine Wave mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::mastersCosineWave(T* v, size_t n)
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
 * @brief Function 15.
 * Implementation of the Quartic mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::quartic(T* v, size_t n)
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
T  mfunc::Functions<T>::w(T x)
{
    return 1.0 + (x - 1.0) / 4.0;
}

/**
 * @brief Function 16.
 * Implementation of the Levy mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::levy(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n - 1; i++)
    {
        T a = w(v[i]) - 1.0;
        a *= a;
        T b = sin(M_PI * w(v[i]) + 1.0);
        b *= b;
        T c = w(v[n - 1]) - 1.0;
        c *= c;
        T d = sin(2.0 * M_PI * w(v[n - 1]));
        d *= d;
        f += a * (1.0 + 10.0 * b) + c * (1.0 + d);
    }

    T e = sin(M_PI * w(v[0]));
    return e*e + f;
}

// ================================================

/**
 * @brief Function 17.
 * Implementation of the Step mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::step(T* v, size_t n)
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
 * @brief Function 18.
 * Implementation of the Alpine mathematical function
 * @param v Vector as a double array
 * @param n Size of the vector 'v'
 * @return The result of the mathematical function 
 */
template <class T>
T mfunc::Functions<T>::alpine(T* v, size_t n)
{
    T f = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        f += fabs(v[i] * sin(v[i]) + 0.1*v[i]);
    }

    return f;
}

// ================================================

template <class T>
mfunc::mfuncPtr<T> mfunc::Functions<T>::get(unsigned int f)
{
    switch (f)
    {
        case 1:
            return Functions<T>::schwefel;
        case 2:
            return Functions<T>::dejong;
        case 3:
            return Functions<T>::rosenbrok;
        case 4:
            return Functions<T>::rastrigin;
        case 5:
            return Functions<T>::griewangk;
        case 6:
            return Functions<T>::sineEnvelopeSineWave;
        case 7:
            return Functions<T>::stretchedVSineWave;
        case 8:
            return Functions<T>::ackleysOne;
        case 9:
            return Functions<T>::ackleysTwo;
        case 10:
            return Functions<T>::eggHolder;
        case 11:
            return Functions<T>::rana;
        case 12:
            return Functions<T>::pathological;
        case 13:
            return Functions<T>::michalewicz;
        case 14:
            return Functions<T>::mastersCosineWave;
        case 15:
            return Functions<T>::quartic;
        case 16:
            return Functions<T>::levy;
        case 17:
            return Functions<T>::step;
        case 18:
            return Functions<T>::alpine;
        default:
            return nullptr;
    }
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
bool mfunc::Functions<T>::exec(unsigned int f, T* v, size_t n, T& outResult)
{
    auto fPtr = get(f);
    if (fPtr == nullptr) return false;

    outResult = fPtr(v, n);
    return true;
}

#endif

// =========================
// End of mfunctions.h
// =========================
