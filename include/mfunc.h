/**
 * @file mfunc.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains various math function definitions
 * @version 0.1
 * @date 2019-03-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __MFUNC_H
#define __MFUNC_H

#include <cstddef> // size_t definition

/**
 * Scope for all math functions
 */
namespace mfunc
{
    extern const unsigned int NUM_FUNCTIONS;

    template <class T>
    using mfuncPtr = T (*)(T*, size_t);

    const char* schwefelDesc();
    const char* dejongDesc();
    const char* rosenbrokDesc();
    const char* rastriginDesc();
    const char* griewangkDesc();
    const char* sineEnvelopeSineWaveDesc();
    const char* stretchedVSineWaveDesc();
    const char* ackleysOneDesc();
    const char* ackleysTwoDesc();
    const char* eggHolderDesc();
    const char* ranaDesc();
    const char* pathologicalDesc();
    const char* michalewiczDesc();
    const char* mastersCosineWaveDesc();
    const char* quarticDesc();
    const char* levyDesc();
    const char* stepDesc();
    const char* alpineDesc();

    template <class T>
    T schwefel(T* v, size_t n);

    template <class T>
    T dejong(T* v, size_t n);

    template <class T>
    T rosenbrok(T* v, size_t n);

    template <class T>
    T rastrigin(T* v, size_t n);

    template <class T>
    T griewangk(T* v, size_t n);
    
    template <class T>
    T sineEnvelopeSineWave(T* v, size_t n);

    template <class T>
    T stretchedVSineWave(T* v, size_t n);

    template <class T>
    T ackleysOne(T* v, size_t n);

    template <class T>
    T ackleysTwo(T* v, size_t n);

    template <class T>
    T eggHolder(T* v, size_t n);

    template <class T>
    T rana(T* v, size_t n);

    template <class T>
    T pathological(T* v, size_t n);

    template <class T>
    T michalewicz(T* v, size_t n);

    template <class T>
    T mastersCosineWave(T* v, size_t n);

    template <class T>
    T quartic(T* v, size_t n);

    template <class T>
    T levy(T* v, size_t n);

    template <class T>
    T step(T* v, size_t n);

    template <class T>
    T alpine(T* v, size_t n);

    template <class T>
    bool fExec(unsigned int f, T* v, size_t n, T& outResult);

    const char* fDesc(unsigned int f);

    template <class T>
    mfuncPtr<T> fGet(unsigned int f);
}

#endif

// =========================
// End of mfunc.h
// =========================
