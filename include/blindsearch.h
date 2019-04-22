/**
 * @file blindsearch.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implements the BlindSearch class, which inherits SearchAlgorithm.
 * BlindSearch::run executes the blind search algorithm on a given population.
 * @version 0.1
 * @date 2019-04-21
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __BLINDSEARCH_H
#define __BLINDSEARCH_H

#include "searchalg.h"

namespace mdata
{
    /**
     * @brief The BlindSearch class implements the Blind Search
     * algorithm, which is ran using the overridden SearchAlgorithm::run() function.
     * 
     * @tparam T Data type used
     */
    template<class T>
    class BlindSearch : public SearchAlgorithm<T>
    {
        // Declaration needed due to template base class
        using SearchAlgorithm<T>::startTimer;
        using SearchAlgorithm<T>::stopTimer;

    public:
        /**
         * @brief Executes Blind Search with the given population and parameters
         * 
         * @param funcPtr Function pointer to the math function being used to generate the population
         * @param fMin Minimum bound for the population matrix vector components
         * @param fMax Maximum bound for the population matrix vector components
         * @param pop Pointer to a population object that will be used in the blind search
         * @param alpha Unused in this algorithm
         * @return TestResult<T> Returns a TestResult struct containing the error code, fitness, and execution time
         */
        virtual TestResult<T> run(mfunc::mfuncPtr<T> funcPtr, const T fMin, const T fMax, Population<T>* const pop, const T alpha)
        {
            // Get population size and dimensions
            size_t popSize = pop->getPopulationSize();
            size_t dimSize = pop->getDimensionsSize();

            // Make sure funcPtr is valid;
            if (funcPtr == nullptr) return TestResult<T>(1, 0, 0.0); // Invalid function id, return with error code 1

            // Start recording execution time
            startTimer();

            // Generate values for population vector matrix
            pop->generate(fMin, fMax);

            // For each population vector, calculate the fitness using the funcPtr
            for (size_t sol = 0; sol < popSize; sol++)
            {
                // Populate fitness values using given math function pointer
                if (!pop->calcFitness(sol, funcPtr))
                    return TestResult<T>(2, 0, 0.0); // Invalid fitness index, return with error code 2
            }

            // Return best fitness value in population
            return TestResult<T>(0, *pop->getBestFitnessPtr(), stopTimer());
        }
    };
}

#endif

// =========================
// End of blindsearch.h
// =========================
