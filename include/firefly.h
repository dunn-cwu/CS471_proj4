#ifndef __FIREFLY_H
#define __FIREFLY_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"
#include "mem.h"

#define BETA_INIT 1.0

namespace mfunc
{
    template <class T>
    struct FAParams
    {
        std:string popFile;
        mdata::DataTable<T>* fitnessTable;
        size_t fitTableCol;
        mdata::Population<T>* mainPop;
        mdata::Population<T>* nextPop;
        mfuncPtr<T> fPtr;
        T fMinBound;
        T fMaxBound;
        unsigned int iterations;
        double alpha;
        double betamin;
        double gamma;

        FAParams()
        {
            popFile = nullptr;
            fitnessTable = nullptr;
            fitTableCol = 0;
            mainPop = nullptr;
            nextPop = nullptr;
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
        int run(FAParams<T> p);
    private:
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<T> rchance;

        void evaluate(FAParams<T>& p, T* solBuffer, size_t firefly);
        void move(FAParams<T>& p, T* solBuffer, size_t firefly_i, size_t firefly_j);
        T calcDistance(T* fv_i, T* fv_j, size_t dimSize);
    };
}

template <class T>
mfunc::Firefly<T>::Firefly()
    : seed(), engine(seed()), rchance(0, 1)
{
}

template <class T>
int mfunc::Firefly<T>::run(FAParams<T> p)
{
    if (p.mainPop == nullptr || p.nextPop == nullptr || p.fPtr == nullptr)
        return 1;

    // Get population information
    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();

    T* solBuffer = util::allocArray<T>(dimSize);
    if (solBuffer == nullptr)
        return 2;

    if (!p.nextPop->generate(p.fMinBound, p.fMaxBound))
        return 3;

    if (!p.nextPop->calcAllFitness(p.fPtr))
        return 4;

    p.nextPop->sortFitnessDescend();

    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        // p.nextPop->sortFitnessAscend();
        p.mainPop->copyAllFrom(p.nextPop);

        for (size_t firefly_i = 0; firefly_i < popSize; firefly_i++)
        {
            evaluate(p, solBuffer, firefly_i);
        }

        // p.nextPop->generateSingle(popSize - 1, p.fMinBound, p.fMaxBound);
        // p.nextPop->calcFitness(popSize - 1, p.fPtr);

        p.nextPop->sortFitnessDescend();
        p.fitnessTable->setEntry(iter, p.fitTableCol, p.nextPop->getFitness(popSize - 1));
    }

    util::releaseArray(solBuffer);

    return 0;
}

template <class T>
void mfunc::Firefly<T>::evaluate(FAParams<T>& p, T* solBuffer, size_t firefly_i)
{
    const size_t popSize = p.mainPop->getPopulationSize();
    
    for (size_t firefly_j = 0; firefly_j < popSize; firefly_j++)
    {
        const T light_j = p.mainPop->getFitness(firefly_j);

        if (p.nextPop->getFitness(firefly_i) < light_j)
        {
            move(p, solBuffer, firefly_j, firefly_i);
        }
    }
}

template <class T>
void mfunc::Firefly<T>::move(FAParams<T>& p, T* solBuffer, size_t firefly_j, size_t firefly_i)
{
    auto alphaDist = std::normal_distribution<T>(0, p.fMaxBound - p.fMinBound);

    const size_t dimSize = p.mainPop->getDimensionsSize();

    auto fv_j = p.mainPop->getPopulationPtr(firefly_j);
    auto fv_i_next = p.nextPop->getPopulationPtr(firefly_i);

    T r = calcDistance(fv_i_next, fv_j, dimSize);
    T beta = (BETA_INIT - p.betamin) * std::pow(static_cast<T>(M_E), -1 * p.gamma * r * r) + p.betamin;

    for (size_t d = 0; d < dimSize; d++)
    {
        T alpha = p.alpha * (rchance(engine) - 0.5) * (std::abs(p.fMaxBound - p.fMinBound));
        solBuffer[d] = fv_j[d] + (beta * (fv_i_next[d] - fv_j[d])) + alpha;

        if (solBuffer[d] < p.fMinBound)
            solBuffer[d] = p.fMinBound;
        else if (solBuffer[d] > p.fMaxBound)
            solBuffer[d] = p.fMaxBound;
    }

    T newFit = p.fPtr(solBuffer, dimSize);
    T oldFit = p.nextPop->getFitness(firefly_j);

    if (newFit < oldFit);
    {
        p.nextPop->copyPopulation(solBuffer, firefly_j);
        p.nextPop->setFitness(firefly_j, newFit);
    }
}

template <class T>
T mfunc::Firefly<T>::calcDistance(T* fv_i, T* fv_j, size_t dimSize)
{
    T sum = 0;
    for (size_t d = 0; d < dimSize; d++)
    {
        T diff = fv_i[d] - fv_j[d];
        sum += diff * diff;
    }

    return std::sqrt(sum);
}

#endif