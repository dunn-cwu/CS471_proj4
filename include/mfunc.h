#ifndef __MFUNC_H
#define __MFUNC_H


#include <cstddef> // size_t definition

namespace mfunc
{
    double schwefel(double* v, size_t n);
    double dejong(double* v, size_t n);
    double rosenbrok(double* v, size_t n);
    double rastrigin(double* v, size_t n);
    double griewangk(double* v, size_t n);
    double sineEnvelopeSineWave(double* v, size_t n);
    double stretchedVSineWave(double* v, size_t n);
    double ackleysOne(double* v, size_t n);
    double ackleysTwo(double* v, size_t n);
    double eggHolder(double* v, size_t n);
}

#endif