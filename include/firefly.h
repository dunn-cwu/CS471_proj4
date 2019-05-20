/**
 * @file firefly.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the Firefly class, which runs the firefly algorithm
 * using the given parameters
 * @version 0.1
 * @date 2019-05-12
 * 
 * @copyright Copyright (c) 2019
 * 
 */

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
#include "stringutils.h"

#define BETA_INIT 1.0
#define POPFILE_GEN_PATTERN "%GEN%"

namespace mfunc
{
    /**
     * @brief The FAParams struct contains various parameters that
     * are required to be passed to the Firefly.run() method.
     * 
     * @tparam T Data type used by the search algorithm
     */
    template <class T>
    struct FAParams
    {
        std::string popFile; // String file name for population dump file
        mdata::DataTable<T>* bestFitnessTable; // Data table for best fitness values
        mdata::DataTable<T>* worstFitnessTable; // Data table for worst fitness values
        size_t fitTableCol; // Data table column for best and worst fitness values
        mdata::Population<T>* mainPop; // Pointer to main population object
        mdata::Population<T>* nextPop; // Pointer to next population object
        mfuncPtr<T> fPtr; // Function pointer to the objective function being tested
        T fMinBound; // Minimum population vector bounds for objective function
        T fMaxBound; // Maximum population vector bounds for objective function
        unsigned int iterations; // Number of iterations to run search algorithm
        double alpha; // Alpha parameter for firefly algorithm
        double betamin; // Betamin parameter for firefly algorithm
        double gamma; // Gamma parameter for firefly algorithm

        /**
         * @brief Construct a new FAParams object
         */
        FAParams()
        {
            popFile = "";
            bestFitnessTable = nullptr;
            worstFitnessTable = nullptr;
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

    /**
     * @brief The Firefly class runs the firefly algorithm with the given
     * parameters passed to the run() method.
     * 
     * @tparam T Data type used by the algorithm
     */
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

/**
 * @brief Construct a new Firefly object
 * 
 * @tparam T Data type used by the algorithm
 */
template <class T>
mfunc::Firefly<T>::Firefly()
    : seed(), engine(seed()), rchance(0, 1)
{
}

/**
 * @brief Runs the firefly algorithm with the given parameters
 * 
 * @tparam T Data type used by the algorithm
 * @param p Parameters for the algorithm
 * @return Returns a non-zero error code on failure, or zero on success
 */
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

    // Generate population vectors
    if (!p.nextPop->generate(p.fMinBound, p.fMaxBound))
        return 3;

    // Calculate fitness for all population vectors
    if (!p.nextPop->calcAllFitness(p.fPtr))
        return 4;

    // Sort population from worst to best
    p.nextPop->sortFitnessDescend();

    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        p.mainPop->copyAllFrom(p.nextPop);

        for (size_t firefly_i = 0; firefly_i < popSize; firefly_i++)
        {
            evaluate(p, solBuffer, firefly_i);
        }

        p.nextPop->sortFitnessDescend();

        // Store best fitness for this iteration
        if (p.bestFitnessTable != nullptr)
            p.bestFitnessTable->setEntry(iter, p.fitTableCol, p.nextPop->getFitness(popSize - 1));

        // Store worst fitness for this iteration
        if (p.worstFitnessTable != nullptr)
            p.worstFitnessTable->setEntry(iter, p.fitTableCol, p.nextPop->getFitness(0));

        // Dump population vectors to file
        if (!p.popFile.empty())
            p.nextPop->outputPopulationCsv(util::s_replace(p.popFile, std::string(POPFILE_GEN_PATTERN), std::to_string(iter)));
    }

    util::releaseArray(solBuffer);

    return 0;
}

/**
 * @brief Private helper function that evaluates a single firefly in the population
 * 
 * @tparam T Data type used by the algorithm
 * @param p Reference to algorithm parameters
 * @param solBuffer Buffer used to store the new generated solution
 * @param firefly_i index of the firefly being evaluated
 */
template <class T>
void mfunc::Firefly<T>::evaluate(FAParams<T>& p, T* solBuffer, size_t firefly_i)
{
    const size_t popSize = p.mainPop->getPopulationSize();
    
    // Compare every other firefly with firefly_i, and move it
    // towards firefly_i if fitness is worse
    for (size_t firefly_j = 0; firefly_j < popSize; firefly_j++)
    {
        const T light_j = p.mainPop->getFitness(firefly_j);

        if (p.nextPop->getFitness(firefly_i) < light_j)
        {
            move(p, solBuffer, firefly_j, firefly_i);
        }
    }
}

/**
 * @brief Private helper function that moves one firefly towards another
 * based off it's relative distance
 * 
 * @tparam T Data type used by the algorithm
 * @param p Reference to algorithm parameters
 * @param solBuffer Buffer used to store the new generated solution
 * @param firefly_j Index of firefly being moved
 * @param firefly_i Index of firefly being moved towards
 */
template <class T>
void mfunc::Firefly<T>::move(FAParams<T>& p, T* solBuffer, size_t firefly_j, size_t firefly_i)
{
    const size_t dimSize = p.mainPop->getDimensionsSize();

    auto fv_j = p.mainPop->getPopulationPtr(firefly_j);
    auto fv_i_next = p.nextPop->getPopulationPtr(firefly_i);

    // Calculate distance between the two fireflies and then their beta value
    T r = calcDistance(fv_i_next, fv_j, dimSize);
    T betaDist = std::pow(static_cast<T>(M_E), -1 * p.gamma * r);
    T beta = (BETA_INIT - p.betamin) * betaDist + p.betamin;

    for (size_t d = 0; d < dimSize; d++)
    {
        // Calculate new value for current dimension
        T alpha = p.alpha * (rchance(engine) - 0.5) * (std::abs(p.fMaxBound - p.fMinBound));
        solBuffer[d] = fv_j[d] + (beta * (fv_i_next[d] - fv_j[d])) + alpha;

        if (solBuffer[d] < p.fMinBound)
            solBuffer[d] = p.fMinBound;
        else if (solBuffer[d] > p.fMaxBound)
            solBuffer[d] = p.fMaxBound;
    }

    // Calculate fitness for new firefly
    T newFit = p.fPtr(solBuffer, dimSize);
    T oldFit = p.nextPop->getFitness(firefly_j);

    // Update firefly if new is better than old
    if (newFit < oldFit);
    {
        p.nextPop->copyPopulation(solBuffer, firefly_j);
        p.nextPop->setFitness(firefly_j, newFit);
    }
}

/**
 * @brief Private helper function that calculates the distance between two vectors
 * 
 * @tparam T Data type used by the algorithm
 * @param fv_i Vector 1
 * @param fv_j Vector 2
 * @param dimSize Number of dimensions in the vectors
 * @return T Distance between the two vectors
 */
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