/**
 * @file harmsearch.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the HarmonySearch class, which runs the
 * harmony search algorithm using the given parameters
 * @version 0.1
 * @date 2019-05-13
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __HARMSEARCH_H
#define __HARMSEARCH_H

#include <cmath>
#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"
#include "mem.h"
#include "stringutils.h"

#define POPFILE_GEN_PATTERN "%GEN%"

namespace mfunc
{
    /**
     * @brief The HSParams struct contains various parameters that
     * are required to be passed to the HarmonySearch.run() method.
     * 
     * @tparam T Data type used by the search algorithm
     */
    template <class T>
    struct HSParams
    {
        std::string popFile; // String file name for population dump file
        mdata::DataTable<T>* bestFitnessTable; // Data table for best fitness values
        mdata::DataTable<T>* worstFitnessTable; // Data table for worst fitness values
        size_t fitTableCol; // Data table column for best and worst fitness values
        mdata::Population<T>* mainPop; // Pointer to main population object
        mfuncPtr<T> fPtr; // Function pointer to the objective function being tested
        T fMinBound; // Minimum population vector bounds for objective function
        T fMaxBound; // Maximum population vector bounds for objective function
        unsigned int iterations; // Number of iterations to run search algorithm
        double hmcr; // HMCR parameter for harmony search
        double par; // PAR parameter for harmony search
        double bw; // BW parameter for harmony search

        /**
         * @brief Construct a new HSParams object
         */
        HSParams()
        {
            popFile = "";
            bestFitnessTable = nullptr;
            worstFitnessTable = nullptr;
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

    /**
     * @brief The HarmonySearch class runs the harmony search algorithm
     * based on the parameters passed to the run() method.
     * 
     * @tparam T Data type used by the algorithm
     */
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

/**
 * @brief Construct a new HarmonySearch object
 * 
 * @tparam T Data type used by the algorithm
 */
template <class T>
mfunc::HarmonySearch<T>::HarmonySearch()
    : seed(), engine(seed()), rchance(0, 1), rrange(-1, 1)
{
}

/**
 * @brief Runs the harmony search algorithm with the given parameters
 * 
 * @tparam T Data type used by the algorithm
 * @param p Parameters for the search algorithm
 * @return Returns a non-zero error code on failure, or zero on success
 */
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

    // Generate random population vectors
    if (!p.mainPop->generate(p.fMinBound, p.fMaxBound))
        return 3;

    // Calculate fitness values for entire population
    if (!p.mainPop->calcAllFitness(p.fPtr))
        return 4;

    // Sort fitness from best to worst
    p.mainPop->sortFitnessAscend();

    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        // Generate new solution
        adjustPitch(p, solBuffer, dimSize);

        // Calculate the new fitness, and replace worst if new solution is better
        T newAesthetic = p.fPtr(solBuffer, dimSize);
        T oldAesthetic = p.mainPop->getFitness(popSize - 1);
        if (newAesthetic < oldAesthetic)
        {
            p.mainPop->copyPopulation(solBuffer, popSize - 1);
            p.mainPop->setFitness(popSize - 1, newAesthetic);
        }

        // Resort population
        p.mainPop->sortFitnessAscend();

        // Store best fitness value for this iteration
        if (p.bestFitnessTable != nullptr)
            p.bestFitnessTable->setEntry(iter, p.fitTableCol, p.mainPop->getFitness(0));

        // Store worst fitness value for this iteration
        if (p.worstFitnessTable != nullptr)
            p.worstFitnessTable->setEntry(iter, p.fitTableCol, p.mainPop->getFitness(popSize - 1));

        // Dump population vectors to a file
        if (!p.popFile.empty())
            p.mainPop->outputPopulationCsv(util::s_replace(p.popFile, std::string(POPFILE_GEN_PATTERN), std::to_string(iter)));
    }

    util::releaseArray<T>(solBuffer);

    return 0;
}

/**
 * @brief Private helper function that generates a new solution vector and stores
 * the values in solBuffer
 * 
 * @tparam T Data type used by the algorithm
 * @param p Reference to the search algorithm parameters
 * @param solBuffer Pointer to new solution vector array
 * @param numDim Number of dimensions in solution vector
 */
template <class T>
void mfunc::HarmonySearch<T>::adjustPitch(HSParams<T>& p, T* solBuffer, const size_t numDim)
{
    // Set up random number distribution for a random population vector
    const size_t popSize = p.mainPop->getPopulationSize();
    auto randPop = std::uniform_int_distribution<size_t>(0, popSize - 1);

    for (size_t dim = 0; dim < numDim; dim++)
    {
        T newPitch = 0;
        if (rchance(engine) <= p.hmcr)
        {
            // Get random value from existing population
            newPitch = p.mainPop->getPopulationPtr(randPop(engine))[dim];
            if (rchance(engine) <= p.par)
            {
                // Adjust pitch of selected value
                newPitch += rrange(engine) * p.bw;
            }
        }
        else
        {
            // Generate a new completely random value for this dimension
            newPitch = (rchance(engine) - 0.5) * std::abs(p.fMaxBound - p.fMinBound);
        }

        solBuffer[dim] = newPitch;
    }
}


#endif