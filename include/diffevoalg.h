#ifndef __DIFFEVOALG_H
#define __DIFFEVOALG_H

#if defined (_MSC_VER)  // Visual studio
    #define thread_local __declspec( thread )
#elif defined (__GCC__) // GCC
    #define thread_local __thread
#endif

#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"

#define MAX_RAND_VECTOR_SELECT 4

namespace mfunc
{

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
        Rand2Bin = 9
    };

    enum class PerturbedVector
    {
        Best,
        Random,
        RandToBest
    };

    enum class NumberDiffVectors
    {
        One,
        Two
    };

    enum class CrossoverStrat
    {
        Exponential,
        Binomial
    };

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
    };

    template <class T>
    struct DifferentialEvolution
    {
        static int run(DEParams<T> params);
        static void mutateAndCrossover(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, size_t vIndex,
            PerturbedVector pertV, NumberDiffVectors diffV, CrossoverStrat crossStrat);
        static void select(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, size_t vIndex);
        static T perturb(T* diffVectorPtrs, NumberDiffVectors numDiffVectors, size_t vIndex, double scalingFactor);
        static void parseDEStrat(DEStrategy mainStrat, PerturbedVector& oPertv, NumberDiffVectors& oDiffv, CrossoverStrat& oCross);
        static void crossover(size_t dim, T* p1, T* p2, double cr, T* outCh1, T* outCh2);
        static void mutate(size_t dim, T* s, double mutProb, double mutRange, double mutPrec, T fMin, T fMax);
        static void reduce(mdata::Population<T>* oldPop, mdata::Population<T>* newPop, size_t elitism);
    };
}

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
    const size_t elitism = p.elitismRate * (double)popSize;

    PerturbedVector pertV;
    NumberDiffVectors diffV;
    CrossoverStrat crStrat;
    parseDEStrat(p.strategy, pertV, diffV, crStrat);

    p.mainPop->setFitnessNormalization(false);
    p.nextPop->setFitnessNormalization(false);
    p.mainPop->generate(p.fMinBound, p.fMaxBound);
    p.mainPop->calcAllFitness(p.fPtr);

    for (unsigned int gen = 0; gen < p.generations; gen++)
    {
        for (size_t i = 0; i < popSize; i++)
        {
            mutateAndCrossover(p.mainPop, p.nextPop, i, pertV, diffV, crStrat);
            select(p.mainPop, p.nextPop, i);
        }

        // Swap the two populations
        auto tmp = p.mainPop;
        p.mainPop = p.auxPop;
        p.auxPop = tmp;

        p.fitnessTable->setEntry(gen, p.fitTableCol, p.mainPop->getBestFitness());
    }

    return 0;
}

template <class T>
void mfunc::DifferentialEvolution<T>::mutateAndCrossover(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, 
    size_t vIndex, PerturbedVector pertV, NumberDiffVectors diffV, CrossoverStrat crossStrat)
{
    static thread_local std::random_device seed;
    static thread_local std::mt19937 engine(seed());
    static thread_local std::uniform_real_distribution<double> rchance(0, 1);

    const dim = mainPop->getDimensionsSize();
    std::uniform_int_distribution<long> rdim(0, dim - 1);
    size_t pertIndex = 0;

    switch (pertV)
    {
        case PerturbedVector::Best:
            pertIndex = mainPop->getBestFitnessIndex();
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

    T* diffVectors[MAX_RAND_VECTOR_SELECT];
    unsigned int numRandIndex = 2;
    if (diffV == NumberDiffVectors::Two)
        numRandIndex = 4;

    for (unsigned int i = 0; i < numRandIndex;)
    {

    }
}

template <class T>
void mfunc::DifferentialEvolution<T>::select(mdata::Population<T>* mainPop, mdata::Population<T>* nextPop, size_t vIndex)
{
    nextPop->calcFitness(vIndex, p.fPtr);
    auto newFit = nextPop->getFitness(vIndex);
    auto oldFit = mainPop->getFitness(vIndex);

    if (newFit > oldFit)
    {
        // Reset (Discard) new pop vector back to last generation
        p.nextPop->copyPopulation(vIndex, p.mainPop->getFitnessPtr(vIndex));
        p.nextPop->setFitness(vIndex, oldFit);
    }
}

template <class T>
T mfunc::DifferentialEvolution<T>::perturb(T* diffVectorPtrs, NumberDiffVectors numDiffVectors, size_t vIndex, double scalingFactor)
{
    if (numDiffVectors == NumberDiffVectors::One)
    {
        return scalingFactor * (diffVectorPtrs[0] - diffVectorPtrs[1]);
    }
    else
    {
        return scalingFactor * (diffVectorPtrs[0] + diffVectorPtrs[1] - diffVectorPtrs[2] - diffVectorPtrs[3]);
    }
}
\
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
        case DEStrategy::Best2Exp
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
        case DEStrategy::Best2Bin
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

template <class T>
void mfunc::DifferentialEvolution<T>::crossover(size_t dim, T* p1, T* p2, double cr, T* outCh1, T* outCh2)
{
    static thread_local std::random_device seed;
    static thread_local std::mt19937 engine(seed());
    static thread_local std::uniform_real_distribution<double> rchance(0, 1);
    std::uniform_int_distribution<long> rdim(0, dim - 1);

    if (rchance(engine) < cr)
    {
        auto crIndex = rdim(engine);
        for (size_t i = 0; i <= dim; i++)
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

template <class T>
void mfunc::DifferentialEvolution<T>::mutate(size_t dim, T* s, double mutProb, double mutRange, double mutPrec, T fMin, T fMax)
{
    static thread_local std::random_device seed;
    static thread_local std::mt19937 engine(seed());
    static thread_local std::uniform_real_distribution<double> rchance(0, 1);
    static thread_local std::uniform_real_distribution<double> rflip(-1, 1);

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

template <class T>
void mfunc::DifferentialEvolution<T>::reduce(mdata::Population<T>* oldPop, mdata::Population<T>* newPop, size_t elitism)
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