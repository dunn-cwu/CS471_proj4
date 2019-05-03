/**
 * @file diffevoalg.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation of the DifferentialEvolution class.
 * Executes the differential evolution algorithm with the specified parameters.
 * @version 0.1
 * @date 2019-04-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef __DIFFEVOALG_H
#define __DIFFEVOALG_H

#include <set>
#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"

#define MAX_RAND_VECTOR_SELECT 4

namespace mfunc
{

    /**
     * @brief Enum used to specify which differential evolution strategy should be used
     */
    enum class DEStrategy
    {
        Best1Exp = 0,
        Rand1Exp = 1,
        RandToBest1Exp = 2,
        Best2Exp = 3,
        Rand2Exp = 4,
        Best1Bin = 5,
        Rand1Bin = 6,
        RandToBest1Bin = 7,
        Best2Bin = 8,
        Rand2Bin = 9,
        Count = 10
    };

    /**
     * @brief Enum used to specify what vector should be perturbed during mutation
     */
    enum class PerturbedVector
    {
        Best,
        Random,
        RandToBest
    };

    /**
     * @brief Enum used to specify the number of random difference vectors used during mutation
     */
    enum class NumberDiffVectors
    {
        One,
        Two
    };

    /**
     * @brief Enum used to specify a crossover strategy
     */
    enum class CrossoverStrat
    {
        Exponential,
        Binomial
    };

    /**
     * @brief Simple structure that holds various parameters to run the
     * differential evolutionary algorithm.
     * 
     * @tparam T Datatype used by the algorithm
     */
    template <class T>
    struct DEParams
    {
        mdata::DataTable<T>* fitnessTable;
        size_t fitTableCol;
        mdata::Population<T>* mainPop;
        mdata::Population<T>* nextPop;
        mfuncPtr<T> fPtr;
        T fMinBound;
        T fMaxBound;
        unsigned int generations;
        double crFactor;
        double scalingFactor1;
        double scalingFactor2;
        DEStrategy strategy;

        DEParams()
        {
            fitnessTable = nullptr;
            fitTableCol = 0;
            mainPop = nullptr;
            nextPop = nullptr;
            fPtr = nullptr;
            fMinBound = 0;
            fMaxBound = 0;
            generations = 0;
            crFactor = 0;
            scalingFactor1 = 0;
            scalingFactor2 = 0;
            strategy = DEStrategy::Best1Exp;
        }
    };

    /**
     * @brief The DifferentialEvolution class executes the differential evolution algorithm
     * on a given population using the given parameters passed via a DEParams structure.
     * To start, call the "run" function.
     * 
     * @tparam T 
     */
    template <class T>
    class DifferentialEvolution
    {
    public:
        DifferentialEvolution();
        ~DifferentialEvolution() = default;
        int run(DEParams<T> params);
    private:
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<double> rchance;

        void mutateAndCrossover(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, const size_t vIndex, const size_t bestIndex,
            const PerturbedVector pertV, const NumberDiffVectors diffV, const CrossoverStrat crossStrat, const double sf1, const double sf2, const double crFactor);
        void select(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, size_t vIndex, mfunc::mfuncPtr<T> fPtr);
        void parseDEStrat(DEStrategy mainStrat, PerturbedVector& oPertv, NumberDiffVectors& oDiffv, CrossoverStrat& oCross);
    };
}

/**
 * @brief Construct a DifferentialEvolution object
 * 
 * @tparam T Datatype used by the algorithm
 */
template <class T>
mfunc::DifferentialEvolution<T>::DifferentialEvolution()
    : seed(), engine(seed()), rchance(0, 1)
{ 
}

/**
 * @brief Runs the algorithm with the given parameters
 * 
 * @tparam T Datatype used by the algorithm
 * @param p Algorithm parameters
 * @return int Non-zero error code on error
 */
template <class T>
int mfunc::DifferentialEvolution<T>::run(DEParams<T> p)
{
    if (p.mainPop == nullptr || p.nextPop == nullptr || p.fPtr == nullptr)
        return 1;

    if (p.mainPop->getPopulationSize() != p.nextPop->getPopulationSize() ||
        p.mainPop->getDimensionsSize() != p.nextPop->getDimensionsSize())
        return 2;

    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();

    // Parse DE strategy
    PerturbedVector pertV;
    NumberDiffVectors diffV;
    CrossoverStrat crStrat;
    parseDEStrat(p.strategy, pertV, diffV, crStrat);

    // Prepare populations
    p.mainPop->setFitnessNormalization(false);
    p.nextPop->setFitnessNormalization(false);
    p.mainPop->generate(p.fMinBound, p.fMaxBound);
    p.mainPop->calcAllFitness(p.fPtr);

    // Calc best fitness pop index
    size_t bestFitIndex = p.mainPop->getBestFitnessIndex();

    for (unsigned int gen = 0; gen < p.generations; gen++)
    {
        for (size_t i = 0; i < popSize; i++)
        {
            // For each population in the next generation, mutate and crossover, then select
            mutateAndCrossover(p.mainPop, p.nextPop, i, bestFitIndex, pertV, diffV, crStrat, p.scalingFactor1, p.scalingFactor2, p.crFactor);
            p.nextPop->boundPopulation(i, p.fMinBound, p.fMaxBound);
            select(p.mainPop, p.nextPop, i, p.fPtr);
        }

        // Swap the two populations
        auto tmp = p.mainPop;
        p.mainPop = p.nextPop;
        p.nextPop = tmp;

        // Recalculate best fitness index and add result to results table
        bestFitIndex = p.mainPop->getBestFitnessIndex();
        p.fitnessTable->setEntry(gen, p.fitTableCol, p.mainPop->getFitness(bestFitIndex));
    }

    return 0;
}

/**
 * @brief Most important part of DE algorithm. Mutates and conducts crossover for a single index of the
 * next population.
 * 
 * @tparam T Datatype used by the algorithm
 * @param mainPop Parent population
 * @param nextPop Next generation population
 * @param vIndex Current vector index
 * @param bestIndex Index of best parent population
 * @param pertV Specify the vector to be perturbed
 * @param diffV Specify the number of random difference vectors
 * @param crossStrat Specify the crossover strategy
 * @param sf1 Scaling factor 1 used in mutation
 * @param sf2 Scaling factor 2 used in mutation
 * @param crFactor Crossover probability factor
 */
template <class T>
void mfunc::DifferentialEvolution<T>::mutateAndCrossover(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, const size_t vIndex, 
    const size_t bestIndex, const PerturbedVector pertV, const NumberDiffVectors diffV, const CrossoverStrat crossStrat, const double sf1, const double sf2, const double crFactor)
{
    const size_t dim = mainPop->getDimensionsSize();
    std::uniform_int_distribution<long> rdim(0, dim - 1);
    size_t pertIndex = 0;

    // Get the index of the perturbed vector
    switch (pertV)
    {
        case PerturbedVector::Best:
            pertIndex = bestIndex;
            break;
        case PerturbedVector::Random:
            do
            {
                pertIndex = rdim(engine);
            } while (pertIndex == vIndex);
            break;
        case PerturbedVector::RandToBest:
            pertIndex = vIndex;
            break;
    }

    // Calculate unique random vector indices
    std::set<size_t> randVectorIndices;

    unsigned int numRandVectors = 2;
    if (diffV == NumberDiffVectors::Two)
        numRandVectors = 4;

    while (randVectorIndices.size() < numRandVectors)
    {
        auto rIndex = rdim(engine);
        if (rIndex != pertIndex && rIndex != bestIndex)
            randVectorIndices.insert(rIndex);
    }

    T* randV[MAX_RAND_VECTOR_SELECT] = { };
    int rvIndex = 0;

    // Convert random vector indices to population vector pointers
    for (auto it = randVectorIndices.begin(); it != randVectorIndices.end(); it++)
    {
        randV[rvIndex] = mainPop->getPopulationPtr(*it);
        rvIndex++;
    }

    // Get population vector points for parent vector, purturbed vector, etc
    T* curVParent = mainPop->getPopulationPtr(vIndex);
    T* curVPert = mainPop->getPopulationPtr(pertIndex);

    nextPop->copyPopulation(vIndex, curVParent);

    T* curVNext = nextPop->getPopulationPtr(vIndex);
    T* curVBest = mainPop->getPopulationPtr(bestIndex);

    size_t count = 0;

    if (crossStrat == CrossoverStrat::Exponential)
    {
        // Mutate with exponential crossover
        // Select random starting dimension
        size_t d = rdim(engine);

        do
        {
            if (pertV == PerturbedVector::RandToBest)
            {
                curVNext[d] = curVParent[d] + (sf2 * (curVBest[d] - curVParent[d])) + (sf1 * (randV[0][d] - randV[1][d]));
            }
            else
            {
                if (numRandVectors == 2)
                    curVNext[d] = curVPert[d] + (sf1 * (randV[0][d] - randV[1][d]));
                else
                    curVNext[d] = curVPert[d] + (sf1 * (randV[0][d] + randV[1][d] - randV[2][d] - randV[3][d]));
            }
            d = (d + 1) % dim;
        } while (rchance(engine) < crFactor);
    }
    else
    {
        // Mutate with binomial crossover

        for (size_t d = 0; d < dim; d++)
        {
            if (rchance(engine) > crFactor)
            continue;

            if (pertV == PerturbedVector::RandToBest)
            {
                curVNext[d] = curVParent[d] + (sf2 * (curVBest[d] - curVParent[d])) + (sf1 * (randV[0][d] - randV[1][d]));
            }
            else
            {
                if (numRandVectors == 2)
                    curVNext[d] = curVPert[d] + (sf1 * (randV[0][d] - randV[1][d]));
                else
                    curVNext[d] = curVPert[d] + (sf1 * (randV[0][d] + randV[1][d] - randV[2][d] - randV[3][d]));
            }
        }
    }
}

/**
 * @brief Selects either the new population vector or the parent population vector
 * depending on which one has the better fitness value
 * 
 * @tparam T Datatype used by the algorithm
 * @param mainPop Parent population pointer
 * @param nextPop Next generation population pointer
 * @param vIndex Index of the population vector being selected
 * @param fPtr Pointer to the objective function for fitness calculation
 */
template <class T>
void mfunc::DifferentialEvolution<T>::select(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, size_t vIndex, mfunc::mfuncPtr<T> fPtr)
{
    nextPop->calcFitness(vIndex, fPtr);
    auto newFit = nextPop->getFitness(vIndex);
    auto oldFit = mainPop->getFitness(vIndex);

    if (newFit > oldFit)
    {
        // Reset (Discard) new pop vector back to last generation
        nextPop->copyPopulation(vIndex, mainPop->getPopulationPtr(vIndex));
        nextPop->setFitness(vIndex, oldFit);
    }
}

/**
 * @brief Splits a DEStrategy into individual components
 * 
 * @tparam T Datatype used by the algorithm
 * @param mainStrat Selected main strategy
 * @param oPertv Vector to be perturbed
 * @param oDiffv Number of random difference vectors during mutation
 * @param oCross Crossover strategy
 */
template <class T>
void mfunc::DifferentialEvolution<T>::parseDEStrat(DEStrategy mainStrat, PerturbedVector& oPertv, NumberDiffVectors& oDiffv, CrossoverStrat& oCross)
{
    switch (mainStrat)
    {
        case DEStrategy::Best1Exp:
            oPertv = PerturbedVector::Best;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Exponential;
            return;
        case DEStrategy::Rand1Exp:
            oPertv = PerturbedVector::Random;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Exponential;
            return;
        case DEStrategy::RandToBest1Exp:
            oPertv = PerturbedVector::RandToBest;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Exponential;
            return;
        case DEStrategy::Best2Exp:
            oPertv = PerturbedVector::Best;
            oDiffv = NumberDiffVectors::Two;
            oCross = CrossoverStrat::Exponential;
            return;
        case DEStrategy::Rand2Exp:
            oPertv = PerturbedVector::Random;
            oDiffv = NumberDiffVectors::Two;
            oCross = CrossoverStrat::Exponential;
            return;
        case DEStrategy::Best1Bin:
            oPertv = PerturbedVector::Best;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Binomial;
            return;
        case DEStrategy::Rand1Bin:
            oPertv = PerturbedVector::Random;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Binomial;
            return;
        case DEStrategy::RandToBest1Bin:
            oPertv = PerturbedVector::RandToBest;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Binomial;
            return;
        case DEStrategy::Best2Bin:
            oPertv = PerturbedVector::Best;
            oDiffv = NumberDiffVectors::Two;
            oCross = CrossoverStrat::Binomial;
            return;
        case DEStrategy::Rand2Bin:
            oPertv = PerturbedVector::Random;
            oDiffv = NumberDiffVectors::Two;
            oCross = CrossoverStrat::Binomial;
            return;
        default:
            oPertv = PerturbedVector::Best;
            oDiffv = NumberDiffVectors::One;
            oCross = CrossoverStrat::Exponential;
            return;
    }
}

#endif