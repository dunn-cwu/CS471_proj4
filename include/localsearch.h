#ifndef __LOCALSEARCH_H
#define __LOCALSEARCH_H

#include <iostream>
#include "searchalg.h"
#include "mem.h"

namespace mdata
{
    template<class T>
    class LocalSearch : public SearchAlgorithm<T>
    {
        using SearchAlgorithm<T>::startTimer;
        using SearchAlgorithm<T>::stopTimer;

    public:
        virtual TestResult<T> run(mfunc::mfuncPtr<T> funcPtr, const T fMin, const T fMax, Population<T>* const pop, const T alpha)
        {
            const size_t popSize = pop->getPopulationSize();
            const size_t dimSize = pop->getDimensionsSize();

            if (funcPtr == nullptr) return TestResult<T>(1, 0, 0.0); // Invalid function id, return with error code 1
            T result = 0;

            bool stop = false;
            size_t pIndex = 0;
            startTimer();

            pop->generate(fMin, fMax);

            for (size_t sol = 0; sol < popSize; sol++)
            {
                // Populate fitness values using given math function pointer
                if (!pop->calcFitness(sol, funcPtr))
                    return TestResult<T>(2, 0, 0.0); // Invalid fitness index, return with error code 2
            }

            pIndex = pop->getBestFitnessIndex();

            T* x = pop->getPopulationPtr(pIndex);
            T xFit = pop->getFitness(pIndex);

            T* y = util::allocArray<T>(dimSize);
            T yFit = 0;

            T* z = util::allocArray<T>(dimSize);
            T zFit = 0;

            if (x == nullptr || y == nullptr || z == nullptr)
            {
                std::cerr << "Error in Local Search: Memory allocation failed" << std::endl;
                return TestResult<T>(3, 0, 0.0);
            }

            while (!stop)
            {   
                stop = true;
                
                util::copyArray<T>(x, y, dimSize);

                for (size_t a = 0; a < dimSize; a++)
                {
                    y[a] = x[a] + alpha;
                    lockBounds(y[a], fMin, fMax);

                    yFit = funcPtr(y, dimSize);
                    z[a] = x[a] - (alpha * (yFit - xFit));
                    lockBounds(z[a], fMin, fMax);

                    y[a] = x[a]; // Reset y[a]
                }

                zFit = funcPtr(z, dimSize);

                if (zFit < xFit)
                {
                    stop = false;

                    T* tmp = x;
                    x = z;
                    xFit = zFit;
                    z = tmp;
                }
            }

            return TestResult<T>(0, xFit, stopTimer());
        }
    private:
        void lockBounds(T& val, const T& min, const T& max)
        {
            if (val < min) val = min;
            else if (val > max) val = max;
        }
    };
}

#endif
