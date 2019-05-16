#ifndef __HARMSEARCH_H
#define __HARMSEARCH_H

#include <cmath>
#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"
#include "mem.h"

namespace mfunc
{
    template <class T>
    struct HSParams
    {
        mdata::DataTable<T>* fitnessTable;
        size_t fitTableCol;
        mdata::Population<T>* mainPop;
        mfuncPtr<T> fPtr;
        T fMinBound;
        T fMaxBound;
        unsigned int iterations;
        double hmcr;
        double par;
        double bw;

        HSParams()
        {
            fitnessTable = nullptr;
            fitTableCol = 0;
            mainPop = nullptr;
            fPtr = nullptr;
            fMinBound = 0;
            fMaxBound = 0;
            iterations = 0;
            hmcr = 0;
            par = 0;
            bw = 0;
        }
    };

    template <class T>
    class HarmonySearch
    {
    public:
        HarmonySearch();
        ~HarmonySearch() = default;
        int run(HSParams<T> p);
    private:
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<T> rchance;
        std::uniform_real_distribution<T> rrange;

        void adjustPitch(HSParams<T>& p, T* solBuffer, const size_t numDim);
    };
} // namespace mfunc

template <class T>
mfunc::HarmonySearch<T>::HarmonySearch()
    : seed(), engine(seed()), rchance(0, 1), rrange(-1, 1)
{
}

template <class T>
int mfunc::HarmonySearch<T>::run(HSParams<T> p)
{
    if (p.mainPop == nullptr || p.fPtr == nullptr)
        return 1;

    // Get population information
    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();

    T* solBuffer = util::allocArray<T>(dimSize);
    if (solBuffer == nullptr)
        return 2;

    if (!p.mainPop->generate(p.fMinBound, p.fMaxBound))
        return 3;

    if (!p.mainPop->calcAllFitness(p.fPtr))
        return 4;

    p.mainPop->sortFitnessAscend();

    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        adjustPitch(p, solBuffer, dimSize);

        T newAesthetic = p.fPtr(solBuffer, dimSize);
        T oldAesthetic = p.mainPop->getFitness(popSize - 1);
        if (newAesthetic < oldAesthetic)
        {
            p.mainPop->copyPopulation(solBuffer, popSize - 1);
            p.mainPop->setFitness(popSize - 1, newAesthetic);
        }

        p.mainPop->sortFitnessAscend();
        p.fitnessTable->setEntry(iter, p.fitTableCol, p.mainPop->getFitness(0));
    }

    return 0;
}

template <class T>
void mfunc::HarmonySearch<T>::adjustPitch(HSParams<T>& p, T* solBuffer, const size_t numDim)
{
    auto randDim = std::uniform_int_distribution<size_t>(0, numDim - 1);

    for (size_t dim = 0; dim < numDim; dim++)
    {
        T newPitch = 0;
        if (rchance(engine) <= p.hmcr)
        {
            newPitch = p.mainPop->getPopulationPtr(randDim(engine))[dim];
            if (rchance(engine) <= p.par)
            {
                newPitch += rrange(engine) * p.bw;
            }
        }
        else
        {
            newPitch = (rchance(engine) - 0.5) * std::abs(p.fMaxBound - p.fMinBound);
        }

        solBuffer[dim] = newPitch;
    }
}


#endif