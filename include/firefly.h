#ifndef __FIREFLY_H
#define __FIREFLY_H

#define _USE_MATH_DEFINES

#include <cmath>
#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"
#include "mem.h"

namespace mfunc
{
    template <class T>
    struct FFParams
    {
        mdata::DataTable<T>* fitnessTable;
        size_t fitTableCol;
        mdata::Population<T>* mainPop;
        mfuncPtr<T> fPtr;
        T fMinBound;
        T fMaxBound;
        unsigned int iterations;
        double alpha;
        double betamin;
        double gamma;

        FFParams()
        {
            fitnessTable = nullptr;
            fitTableCol = 0;
            mainPop = nullptr;
            fPtr = nullptr;
            fMinBound = 0;
            fMaxBound = 0;
            iterations = 0;
            alpha = 0;
            betamin = 0;
            gamma = 0;
        }
    };

    template <class T>
    class Firefly
    {
    public:
        Firefly();
        ~Firefly() = default;
        int run(FFParams<T> p);
    private:
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<double> rchance;

        void evaluate(FFParams<T>& p, size_t firefly);
        void move(FFParams<T>& p, size_t firefly1, size_t firefly2);
        T calcDistance(mdata::Population<T>* pop, size_t firefly1, size_t firefly2);
    };
}

template <class T>
mfunc::Firefly<T>::Firefly()
    : seed(), engine(seed()), rchance(0, 1)
{
}

template <class T>
int mfunc::Firefly<T>::run(FFParams<T> p)
{
    if (p.mainPop == nullptr || p.fPtr == nullptr)
        return 1;

    // Get population information
    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();

    if (!p.mainPop->generate(p.fMinBound, p.fMaxBound))
        return 2;

    if (!p.mainPop->calcAllFitness(p.fPtr))
        return 3;


    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        for (size_t firefly = 0; firefly < popSize; firefly++)
        {
            evaluate(p, firefly);
        }

        p.fitnessTable->setEntry(iter, p.fitTableCol, p.mainPop->getBestFitness());
    }

    return 0;
}

template <class T>
void mfunc::Firefly<T>::evaluate(FFParams<T>& p, size_t firefly)
{
    const size_t popSize = p.mainPop->getPopulationSize();
    T lightI = p.mainPop->getFitness(firefly);

    for (size_t otherFirefly = 0; otherFirefly < popSize; otherFirefly++)
    {
        const T lightO = p.mainPop->getFitness(otherFirefly);

        if (lightO < lightI)
        {
            move(p, firefly, otherFirefly);
            p.mainPop->calcFitness(firefly, p.fPtr);
            lightI = p.mainPop->getFitness(firefly);
        }
    }
}

template <class T>
void mfunc::Firefly<T>::move(FFParams<T>& p, size_t firefly1, size_t firefly2)
{
    auto alphaDist = std::normal_distribution<T>(0, p.fMaxBound - p.fMinBound);

    auto pop = p.mainPop;
    const size_t dimSize = pop->getDimensionsSize();
    auto fv1 = pop->getPopulationPtr(firefly1);
    auto fv2 = pop->getPopulationPtr(firefly2);

    T r = calcDistance(p.mainPop, firefly1, firefly2);

    for (size_t d = 0; d < dimSize; d++)
    {
        T beta = p.betamin * std::pow(static_cast<T>(M_E), -1 * p.gamma * r * r);
        fv1[d] += beta * (fv2[d] - fv1[d]) + p.alpha * alphaDist(engine);

        if (fv1[d] < p.fMinBound)
            fv1[d] = p.fMinBound;
        else if (fv1[d] > p.fMaxBound)
            fv1[d] = p.fMaxBound;
    }
}

template <class T>
T mfunc::Firefly<T>::calcDistance(mdata::Population<T>* pop, size_t firefly1, size_t firefly2)
{
    const size_t dimSize = pop->getDimensionsSize();
    auto fv1 = pop->getPopulationPtr(firefly1);
    auto fv2 = pop->getPopulationPtr(firefly2);

    T sum = 0;
    for (size_t d = 0; d < dimSize; d++)
    {
        T diff = fv1[d] - fv2[d];
        sum += diff * diff;
    }

    return std::sqrt(sum);
}

#endif