/**
 * @file localsearch.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief 
 * @version 0.1
 * @date 2019-04-19
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __LOCALSEARCH_H
#define __LOCALSEARCH_H

#include <iostream>
#include <cmath>
#include "searchalg.h"
#include "mem.h"

// Precision when checking neighbor fitness in number of decimal places right side of zero.
#define DEC_PRECISION 12

namespace mdata
{
    /**
     * @brief The LocalSearch class implements the Local Search
     * algorithm, which is ran using the overridden SearchAlgorithm::run() function.
     * 
     * @tparam T Data type used
     */
    template<class T>
    class LocalSearch : public SearchAlgorithm<T>
    {
        using SearchAlgorithm<T>::startTimer;
        using SearchAlgorithm<T>::stopTimer;
        const T MIN_IMPROVEMENT = pow(static_cast<T>(10), static_cast<T>(-1 * DEC_PRECISION));

    public:
        /**
         * @brief Executes Local Search with the given population and parameters
         * 
         * @param funcPtr Function pointer to the math function being used to generate the population
         * @param fMin Minimum bound for the population matrix vector components
         * @param fMax Maximum bound for the population matrix vector components
         * @param pop Pointer to a population object that will be used in the local search
         * @param alpha Alpha value for local search neighbor generation
         * @return TestResult<T> Returns a TestResult struct containing the error code, fitness, and execution time
         */
        virtual TestResult<T> run(mfunc::mfuncPtr<T> funcPtr, const T fMin, const T fMax, Population<T>* const pop, const T alpha)
        {
            // Get population size and dimensions
            const size_t popSize = pop->getPopulationSize();
            const size_t dimSize = pop->getDimensionsSize();

            // Make sure funcPtr is valid;
            if (funcPtr == nullptr) return TestResult<T>(1, 0, 0.0); // Invalid function id, return with error code 1

            // Algorithm related variables
            bool stop = false;
            size_t pIndex = 0;

            startTimer();

            // Start recording execution time
            pop->generate(fMin, fMax);

            for (size_t sol = 0; sol < popSize; sol++)
            {
                // Populate fitness values using given math function pointer
                if (!pop->calcFitness(sol, funcPtr))
                    return TestResult<T>(2, 0, 0.0); // Invalid fitness index, return with error code 2
            }

            // Get the index for the best fitness in the population
            pIndex = pop->getBestFitnessIndex();

            // Get population vector and fitness of best solution
            T* x = pop->getPopulationPtr(pIndex);
            T xFit = pop->getFitness(pIndex);

            // Create empty Y vector
            T* y = util::allocArray<T>(dimSize);
            T yFit = 0;

            // Create empty Z vector
            T* z = util::allocArray<T>(dimSize);
            T zFit = 0;

            if (x == nullptr || y == nullptr || z == nullptr)
            {
                std::cerr << "Error in Local Search: Memory allocation failed" << std::endl;
                return TestResult<T>(3, 0, 0.0);
            }

            // Keep looping until search fails to improve
            while (!stop)
            {   
                stop = true;
                
                // Copy values from X vector into Y vector
                util::copyArray<T>(x, y, dimSize);

                // Loop through each dimension in vector
                for (size_t a = 0; a < dimSize; a++)
                {
                    // Add alpha value to y[a]
                    y[a] = x[a] + alpha;

                    // Make sure y[a] is within the function bounds
                    lockBounds(y[a], fMin, fMax); 

                    // Calculate fitness for y vector
                    yFit = funcPtr(y, dimSize);

                    // Update Z[a] vector value based on the difference between Y and X
                    z[a] = x[a] - (alpha * (yFit - xFit));

                    // Make sure z[a] is within the function bounds
                    lockBounds(z[a], fMin, fMax);

                    y[a] = x[a]; // Reset y[a] to prepare for next loop
                }

                zFit = funcPtr(z, dimSize);

                // The following 'if' statement may cause extreme execution 
                // times for some functions due to floating point precision:
                // if (zFit < xFit)
                //
                // The replacement 'if' statement below places a limit
                // on the minimum acknowledged improvement fitness,
                // hopefully preventing extreme run-times:
                if (xFit - zFit > MIN_IMPROVEMENT)
                {
                    // Z is an improvement on X, so keep searching
                    stop = false;

                    // Swap Z and X for next loop
                    T* tmp = x;
                    x = z;
                    xFit = zFit;
                    z = tmp;
                }
            }

            // Return best result
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

// =========================
// End of localsearch.h
// =========================
