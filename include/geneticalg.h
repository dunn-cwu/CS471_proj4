/**
 * @file geneticalg.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation of the GeneticAlgorithm class.
 * Executes the genetic algorithm with the specified parameters.
 * @version 0.1
 * @date 2019-04-27
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __GENETICALG_H
#define __GENETICALG_H

#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"
#include "mem.h"

namespace mfunc
{
    /**
     * @brief Simple structure that holds various parameters for the
     * genetic algorithm
     * 
     * @tparam T Datatype used by the algorithm
     */
    template <class T>
    struct GAParams
    {
        mdata::DataTable<T>* fitnessTable;
        size_t fitTableCol;
        mdata::Population<T>* mainPop;
        mdata::Population<T>* auxPop;
        mfuncPtr<T> fPtr;
        T fMinBound;
        T fMaxBound;
        unsigned int generations;
        double crProb;
        double mutProb;
        double mutRange;
        double mutPrec;
        double elitismRate;

        GAParams()
        {
            fitnessTable = nullptr;
            fitTableCol = 0;
            mainPop = nullptr;
            auxPop = nullptr;
            fPtr = nullptr;
            fMinBound = 0;
            fMaxBound = 0;
            generations = 0;
            crProb = 0;
            mutProb = 0;
            mutRange = 0;
            mutPrec = 0;
            elitismRate = 0;
        }
    };

    /**
     * @brief The GeneticAlgorithm class executes the genetic algorithm
     * with the specified parameters. To start, execute the run() function.
     * 
     * @tparam T Datatype used by the algorithm
     */
    template <class T>
    class GeneticAlgorithm
    {
    public:
        GeneticAlgorithm();
        ~GeneticAlgorithm() = default;
        int run(GAParams<T> params);
        
    private:
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<double> rchance;

        void select(mdata::Population<T>* pop, size_t& outP1, size_t& outP2);
        size_t selRW(mdata::Population<T>* pop);
        void crossover(size_t dim, T* p1, T* p2, double cr, T* outCh1, T* outCh2);
        void mutate(size_t dim, T* s, double mutProb, double mutRange, double mutPrec, T fMin, T fMax);
        void reduce(mdata::Population<T>* oldPop, mdata::Population<T>* newPop, size_t elitism);
    };
}

/**
 * @brief Construct a new GeneticAlgorithm object
 * 
 * @tparam T Datatype used by the algorithm
 */
template <class T>
mfunc::GeneticAlgorithm<T>::GeneticAlgorithm()
    : seed(), engine(seed()), rchance(0, 1)
{
}

/**
 * @brief Executes the genetic algorithm with the specified parameters
 * 
 * @tparam T Datatype used by the algorithm
 * @param p GA parameters
 * @return int Non-zero error code on failure
 */
template <class T>
int mfunc::GeneticAlgorithm<T>::run(GAParams<T> p)
{
    if (p.mainPop == nullptr || p.auxPop == nullptr || p.fPtr == nullptr)
        return 1;

    if (p.mainPop->getPopulationSize() != p.auxPop->getPopulationSize() ||
        p.mainPop->getDimensionsSize() != p.auxPop->getDimensionsSize())
        return 2;

    // Get population information
    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();
    const size_t elitism = p.elitismRate * (double)popSize;

    // Allocate child buffers
    T* childOne = util::allocArray<T>(dimSize);
    T* childTwo = util::allocArray<T>(dimSize);

    // Prepare populations
    p.mainPop->setFitnessNormalization(true);
    p.auxPop->setFitnessNormalization(true);
    p.mainPop->generate(p.fMinBound, p.fMaxBound);
    p.mainPop->calcAllFitness(p.fPtr);

    // Loop for a number of generations
    for (unsigned int gen = 0; gen < p.generations; gen++)
    {
        for (size_t s = 0; s < popSize;)
        {
            size_t p1Index = 0;
            size_t p2Index = 0;

            // Select parent indices
            select(p.mainPop, p1Index, p2Index);

            // Produce new children by computing the crossover between parents
            crossover(dimSize, p.mainPop->getPopulationPtr(p1Index), p.mainPop->getPopulationPtr(p2Index), 
                p.crProb, childOne, childTwo);

            // Mutate children
            mutate(dimSize, childOne, p.mutProb, p.mutRange, p.mutPrec, p.fMinBound, p.fMaxBound);
            mutate(dimSize, childTwo, p.mutProb, p.mutRange, p.mutPrec, p.fMinBound, p.fMaxBound);

            // Copy new children into next generation
            p.auxPop->copyPopulation(s, childOne);
            s++;
            if (s < popSize) p.auxPop->copyPopulation(s, childTwo);
            s++;
        }

        // Recalculate all fitness values for next generation
        p.auxPop->calcAllFitness(p.fPtr);

        // Select and combine some of the best populations from the previous generation
        // with the next generation
        reduce(p.mainPop, p.auxPop, elitism);

        // Swap the two populations
        auto tmp = p.mainPop;
        p.mainPop = p.auxPop;
        p.auxPop = tmp;

        // Record current best population in results table
        p.fitnessTable->setEntry(gen, p.fitTableCol, p.mainPop->getMinCost());
    }

    // Delete children buffers
    util::releaseArray<T>(childOne);
    util::releaseArray<T>(childTwo);

    return 0;
}

/**
 * @brief Performs a selection using the roulette wheel selection method
 * There were supposed to be multiple selection algorithms here, but
 * due to multiple coding setbacks I ran out of time to implement them
 * correctly.
 * 
 * @tparam T Datatype used by the algorithm
 * @param pop Pointer to the parent population
 * @param outP1 Out reference variable for first parent selection index
 * @param outP2 Out reference variable for second parent selection index
 */
template <class T>
void mfunc::GeneticAlgorithm<T>::select(mdata::Population<T>* pop, size_t& outP1, size_t& outP2)
{
    outP1 = selRW(pop);
    outP2 = selRW(pop);
}

/**
 * @brief Performs a selection using the roulette wheel selection method
 * 
 * @tparam T Datatype used by the algorithm
 * @param pop pop Pointer to the parent population
 * @return size_t Selected parent index
 */
template <class T>
size_t mfunc::GeneticAlgorithm<T>::selRW(mdata::Population<T>* pop)
{
    std::uniform_real_distribution<T> dist(0, pop->getTotalFitness());
    const size_t pSize = pop->getPopulationSize();
    T r = dist(engine);
    size_t s = 0;

    while (s < pSize - 1 && r > 0)
    {
        r -= pop->getFitness(s);
        s += 1;
    }

    return s;
}

/**
 * @brief Performs a crossover between the two parents to produce two children
 * 
 * @tparam T Datatype used by the algorithm
 * @param dim Number of dimensions for the population
 * @param p1 Pointer to the first parent
 * @param p2 Pointer to the second parent
 * @param cr Crossover probability rate
 * @param outCh1 Pointer to the first child that will be produced
 * @param outCh2 Pointer to the second child that will be produced
 */
template <class T>
void mfunc::GeneticAlgorithm<T>::crossover(size_t dim, T* p1, T* p2, double cr, T* outCh1, T* outCh2)
{
    std::uniform_int_distribution<long> rdim(0, dim - 1);

    if (rchance(engine) < cr)
    {
        auto crIndex = rdim(engine);
        for (size_t i = 0; i < dim; i++)
        {
            if (i < crIndex)
            {
                outCh1[i] = p1[i];
                outCh2[i] = p2[i];
            }
            else
            {
                outCh1[i] = p2[i];
                outCh2[i] = p1[i];
            } 
        }
    }
    else
    {
        for (size_t i = 0; i < dim; i++)
        {
            outCh1[i] = p1[i];
            outCh2[i] = p2[i];
        }
    }
}

/**
 * @brief Performs a mutation on the given child population
 * 
 * @tparam T Datatype used by the algorithm
 * @param dim Number of dimensions for the population
 * @param s Pointer to the child population vector
 * @param mutProb Mutation probability factor
 * @param mutRange Mutation range
 * @param mutPrec Mutation precision
 * @param fMin Minimum objective function bounds
 * @param fMax Maximum objective function bounds
 */
template <class T>
void mfunc::GeneticAlgorithm<T>::mutate(size_t dim, T* s, double mutProb, double mutRange, double mutPrec, T fMin, T fMax)
{
    std::uniform_real_distribution<double> rflip(-1, 1);

    for (size_t i = 0; i < dim; i++)
    {
        if (rchance(engine) < mutProb)
        {
            s[i] += rflip(engine) * (fMax - fMin) * mutRange * std::pow(2, (-1 * rchance(engine) * mutPrec));

            if (s[i] < fMin) s[i] = fMin;
            else if (s[i] > fMax) s[i] = fMax;
        }
    }
}

/**
 * @brief Selects a portion of the best last generation population to combine with the next generation population
 * 
 * @tparam T Datatype used by the algorithm
 * @param oldPop Pointer to the last gen population
 * @param newPop Pointer to the next gen population
 * @param elitism The number of populations that will be taken from the last generation and placed in the next generation
 */
template <class T>
void mfunc::GeneticAlgorithm<T>::reduce(mdata::Population<T>* oldPop, mdata::Population<T>* newPop, size_t elitism)
{
    oldPop->sortDescendByFitness();
    newPop->sortDescendByFitness();
    auto const lastIndex = newPop->getPopulationSize() - 1;

    for (size_t i = 0; i < elitism; i++)
    {
        newPop->copyPopulation(lastIndex - i, oldPop->getPopulationPtr(i));
    }
}

#endif