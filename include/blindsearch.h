#ifndef __BLINDSEARCH_H
#define __BLINDSEARCH_H

#include "searchalg.h"
#include "mfunc.h"

namespace mdata
{
    template<class T>
    class BlindSearch : public SearchAlgorithm<T>
    {
        using SearchAlgorithm<T>::startTimer;
        using SearchAlgorithm<T>::stopTimer;

    public:
        virtual TestResult run(const unsigned int funcId, const T fMin, const T fMax, Population<T>* const pop, const size_t iterations, const T alpha)
        {
            pop->clearBest();

            size_t popSize = pop->getPopulationSize();
            size_t dimSize = pop->getDimensionsSize();

            // Get function pointer for math function we want to test
            mfunc::mfuncPtr<T> fPtr = mfunc::fGet<T>(funcId);
            if (fPtr == nullptr) return mdata::TestResult(1, 0.0); // Invalid function id, return with error code 1

            startTimer();

            for (size_t i = 0; i < iterations; i++)
            {
                pop->generate(fMin, fMax);

                for (size_t sol = 0; sol < popSize; sol++)
                {
                    // Populate fitness values from given math function pointer
                    if (!pop->setFitness(sol, fPtr))
                        return mdata::TestResult(2, 0.0); // Invalid fitness index, return with error code 2
                }

                pop->storeBest();
            }

            return TestResult(0, stopTimer());
        }
    };
}

#endif
