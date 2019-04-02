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

    const char* schwefelDesc();
    double schwefel(double* v, size_t n);

    const char* dejongDesc();
    double dejong(double* v, size_t n);

    const char* rosenbrokDesc();
    double rosenbrok(double* v, size_t n);

    const char* rastriginDesc();
    double rastrigin(double* v, size_t n);

    const char* griewangkDesc();
    double griewangk(double* v, size_t n);

    const char* sineEnvelopeSineWaveDesc();
    double sineEnvelopeSineWave(double* v, size_t n);

    const char* stretchedVSineWaveDesc();
    double stretchedVSineWave(double* v, size_t n);

    const char* ackleysOneDesc();
    double ackleysOne(double* v, size_t n);

    const char* ackleysTwoDesc();
    double ackleysTwo(double* v, size_t n);

    const char* eggHolderDesc();
    double eggHolder(double* v, size_t n);

    const char* ranaDesc();
    double rana(double* v, size_t n);

    const char* pathologicalDesc();
    double pathological(double* v, size_t n);

    const char* michalewiczDesc();
    double michalewicz(double* v, size_t n);

    const char* mastersCosineWaveDesc();
    double mastersCosineWave(double* v, size_t n);

    const char* quarticDesc();
    double quartic(double* v, size_t n);

    const char* levyDesc();
    double levy(double* v, size_t n);

    const char* stepDesc();
    double step(double* v, size_t n);

    const char* alpineDesc();
    double alpine(double* v, size_t n);

    bool fExec(unsigned int f, double* v, size_t n, double& outResult);
    const char* fDesc(unsigned int f);
}

#endif

// =========================
// End of mfunc.h
// =========================
