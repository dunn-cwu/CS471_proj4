#ifndef __BLINDSEARCH_H
#define __BLINDSEARCH_H

#include "searchalg.h"

namespace mdata
{
    template<class T>
    class BlindSearch : public SearchAlgorithm<T>
    {
        using SearchAlgorithm<T>::startTimer;
        using SearchAlgorithm<T>::stopTimer;

    public:
        virtual TestResult<T> run(mfunc::mfuncPtr<T> funcPtr, const T fMin, const T fMax, Population<T>* const pop, const T alpha)
        {
            size_t popSize = pop->getPopulationSize();
            size_t dimSize = pop->getDimensionsSize();

            if (funcPtr == nullptr) return TestResult<T>(1, 0, 0.0); // Invalid function id, return with error code 1

            startTimer();

            pop->generate(fMin, fMax);

            for (size_t sol = 0; sol < popSize; sol++)
            {
                // Populate fitness values using given math function pointer
                if (!pop->calcFitness(sol, funcPtr))
                    return TestResult<T>(2, 0, 0.0); // Invalid fitness index, return with error code 2
            }

            return TestResult<T>(0, pop->getBestFitness(), stopTimer());
        }
    };
}

#endif
