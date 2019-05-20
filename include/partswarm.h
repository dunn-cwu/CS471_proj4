/**
 * @file partswarm.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the ParticleSwarm class, which runs the
 * particle swarm algorithm using the given parameters
 * @version 0.1
 * @date 2019-05-10
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __PARTSWARM_H
#define __PARTSWARM_H

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
     * @brief The Particle struct is a simple data structure used to
     * store the global best particle along with it's fitness
     * 
     * @tparam T 
     */
    template <class T>
    struct Particle
    {
        T* vector;
        T fitness;

        Particle()
            : vector(nullptr), fitness(0)
        {
        }
    };

    /**
     * @brief The PSParams struct contains various parameters that
     * are required to be passed to the ParticleSwarm.run() method.
     * 
     * @tparam T Data type used by the search algorithm
     */
    template <class T>
    struct PSParams
    {
        std::string popFile; // String file name for population dump file
        mdata::DataTable<T>* bestFitnessTable; // Data table for best fitness values
        mdata::DataTable<T>* worstFitnessTable; // Data table for worst fitness values
        size_t fitTableCol; // Data table column for best and worst fitness values
        mdata::Population<T>* mainPop; // Pointer to main population object
        mdata::Population<T>* pbPop; // Pointer to personal best population object
        mfuncPtr<T> fPtr; // Function pointer to the objective function being tested
        T fMinBound; // Minimum population vector bounds for objective function
        T fMaxBound; // Maximum population vector bounds for objective function
        unsigned int iterations; // Number of iterations to run search algorithm
        double c1; // C1 parameter for particle swarm
        double c2; // C2 parameter for particle swarm
        double k;  // k dampening factor parameter for particle swarm

        /**
         * @brief Construct a new PSParams object
         */
        PSParams()
        {
            popFile = "";
            bestFitnessTable = nullptr;
            worstFitnessTable = nullptr;
            fitTableCol = 0;
            mainPop = nullptr;
            pbPop = nullptr;
            fPtr = nullptr;
            fMinBound = 0;
            fMaxBound = 0;
            iterations = 0;
            c1 = 0;
            c2 = 0;
            k = 0;
        }
    };

    /**
     * @brief The ParticleSwarm class runs the particle swarm algorithm
     * with the given parameters passed to the run() method.
     * 
     * @tparam T Data type used by the search algorithm
     */
    template <class T>
    class ParticleSwarm
    {
    public:
        ParticleSwarm();
        ~ParticleSwarm() = default;
        int run(PSParams<T> params);
    private:
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<double> rchance;

        void updateParticle(PSParams<T>& p, const Particle<T>& globalBest, T** velMatrix, size_t pIndex);
        void randomizeVelocity(T** vMatrix, size_t popSize, size_t dimSize, T fMin, T fMax);
    };
}

/**
 * @brief Construct a new ParticleSwarm object
 * @tparam T Data type used by the search algorithm
 */
template <class T>
mfunc::ParticleSwarm<T>::ParticleSwarm()
    : seed(), engine(seed()), rchance(0, 1)
{
}

/**
 * @brief Runs the particle swarm algorithm with the given parameters
 * 
 * @tparam T Data type used by the search algorithm
 * @param p Parameters used by the search algorithm
 * @return Returns a non-zero error code on failure, or zero on success
 */
template <class T>
int mfunc::ParticleSwarm<T>::run(PSParams<T> p)
{
    if (p.mainPop == nullptr || p.pbPop == nullptr || p.fPtr == nullptr)
        return 1;

    // Get population information
    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();

    if (popSize != p.pbPop->getPopulationSize() ||
        dimSize != p.pbPop->getDimensionsSize())
        return 2;


    // Construct global best particle and allocate gBest vector
    Particle<T> globalBest;
    globalBest.vector = util::allocArray<T>(dimSize);

    // Allocate velocity matrix
    T** velocityMatrix = util::allocMatrix<T>(popSize, dimSize);

    if (globalBest.vector == nullptr || velocityMatrix == nullptr)
        return 3;

    if (!p.mainPop->generate(p.fMinBound, p.fMaxBound))
        return 4;
    
    if (!p.mainPop->calcAllFitness(p.fPtr))
        return 5;

    if (!p.pbPop->copyAllFrom(p.mainPop))
        return 6;

    
    // Randomize the velocities for all particles
    randomizeVelocity(velocityMatrix, popSize, dimSize, p.fMinBound, p.fMaxBound);

    auto bestFitIndex = p.mainPop->getBestFitnessIndex();
    util::copyArray<T>(p.mainPop->getPopulationPtr(bestFitIndex), globalBest.vector, dimSize);
    globalBest.fitness = p.mainPop->getFitness(bestFitIndex);

    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        for (size_t pIndex = 0; pIndex < popSize; pIndex++)
        {
            // Update the particles and their velocities
            updateParticle(p, globalBest, velocityMatrix, pIndex);
        }

        // Get the index of current the best solution, and the associated fitness
        bestFitIndex = p.mainPop->getBestFitnessIndex();
        T bestFitVal = p.mainPop->getFitness(bestFitIndex);

        // Update global best if current best is better
        if (bestFitVal < globalBest.fitness)
        {
            util::copyArray<T>(p.mainPop->getPopulationPtr(bestFitIndex), globalBest.vector, dimSize);
            globalBest.fitness = bestFitVal;
        }

        // Store best fitness for this iteration
        if (p.bestFitnessTable != nullptr)
            p.bestFitnessTable->setEntry(iter, p.fitTableCol, globalBest.fitness);

        // Store worst fitness for this iteration
        if (p.worstFitnessTable != nullptr)
            p.worstFitnessTable->setEntry(iter, p.fitTableCol, p.mainPop->getWorstFitness());

        // Dump population vectors to file
        if (!p.popFile.empty())
            p.mainPop->outputPopulationCsv(util::s_replace(p.popFile, std::string(POPFILE_GEN_PATTERN), std::to_string(iter)));
    }

    util::releaseArray<T>(globalBest.vector);
    util::releaseMatrix<T>(velocityMatrix, popSize);

    return 0;
}

/**
 * @brief Private helper function that updates a specific particle and it's velocity
 * 
 * @tparam T Data type used by the search algorithm
 * @param p Reference to the search algorithm parameters
 * @param globalBest Reference to the current global best particle
 * @param velMatrix Pointer to the particle velocity matrix
 * @param pIndex Index of the particle being updated
 */
template <class T>
void mfunc::ParticleSwarm<T>::updateParticle(PSParams<T>& p, const Particle<T>& globalBest, T** velMatrix, size_t pIndex)
{
    const size_t dimSize = p.mainPop->getDimensionsSize();
    auto pBestVector = p.pbPop->getPopulationPtr(pIndex);
    auto curVector = p.mainPop->getPopulationPtr(pIndex);

    // Update particle's velocity and position
    for (size_t d = 0; d < dimSize; d++)
    {
        velMatrix[pIndex][d] += p.c1 * rchance(engine) * (pBestVector[d] - curVector[d]) + 
            p.c2 * rchance(engine) * (globalBest.vector[d] - curVector[d]);
        velMatrix[pIndex][d] *= p.k;

        curVector[d] += velMatrix[pIndex][d];

        if (curVector[d] < p.fMinBound)
            curVector[d] = p.fMinBound;
        else if (curVector[d] > p.fMaxBound)
            curVector[d] = p.fMaxBound;
    }
    
    p.mainPop->calcFitness(pIndex, p.fPtr);
    T newFitness = p.mainPop->getFitness(pIndex);
    T pbFitness = p.pbPop->getFitness(pIndex);

    // Update personal best if current position is better
    if (newFitness < pbFitness)
    {
        p.pbPop->copyFrom(p.mainPop, pIndex, pIndex);
    }
}

/**
 * @brief Helper function that randomizes the velocities in the velocity matrix
 * 
 * @tparam T Data type used by the search algorithm
 * @param vMatrix Pointer to the velocity matrix
 * @param popSize Number of particles in population
 * @param dimSize Number of dimensions in each particle
 * @param fMin Minimum bound for objective function
 * @param fMax Maximum bound for objective function
 */
template <class T>
void mfunc::ParticleSwarm<T>::randomizeVelocity(T** vMatrix, size_t popSize, size_t dimSize, T fMin, T fMax)
{
    std::uniform_real_distribution<T> velDist(0, 0.5 * (fMax - fMin));

    for (size_t s = 0; s < popSize; s++)
    {
        for (size_t d = 0; d < dimSize; d++)
        {
            vMatrix[s][d] = velDist(engine);
        }
    }
}



#endif