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

    template <class T>
    struct PSParams
    {
        std::string popFile;
        mdata::DataTable<T>* bestFitnessTable;
        mdata::DataTable<T>* worstFitnessTable;
        size_t fitTableCol;
        mdata::Population<T>* mainPop;
        mdata::Population<T>* pbPop;
        mfuncPtr<T> fPtr;
        T fMinBound;
        T fMaxBound;
        unsigned int iterations;
        double c1;
        double c2;
        double k;

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

template <class T>
mfunc::ParticleSwarm<T>::ParticleSwarm()
    : seed(), engine(seed()), rchance(0, 1)
{
}

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

    Particle<T> globalBest;
    globalBest.vector = util::allocArray<T>(dimSize);
    T** velocityMatrix = util::allocMatrix<T>(popSize, dimSize);

    if (globalBest.vector == nullptr || velocityMatrix == nullptr)
        return 3;

    if (!p.mainPop->generate(p.fMinBound, p.fMaxBound))
        return 4;
    
    if (!p.mainPop->calcAllFitness(p.fPtr))
        return 5;

    if (!p.pbPop->copyAllFrom(p.mainPop))
        return 6;

    
    randomizeVelocity(velocityMatrix, popSize, dimSize, p.fMinBound, p.fMaxBound);

    auto bestFitIndex = p.mainPop->getBestFitnessIndex();
    util::copyArray<T>(p.mainPop->getPopulationPtr(bestFitIndex), globalBest.vector, dimSize);
    globalBest.fitness = p.mainPop->getFitness(bestFitIndex);

    for (unsigned int iter = 0; iter < p.iterations; iter++)
    {
        for (size_t pIndex = 0; pIndex < popSize; pIndex++)
        {
            updateParticle(p, globalBest, velocityMatrix, pIndex);
        }

        bestFitIndex = p.mainPop->getBestFitnessIndex();
        T bestFitVal = p.mainPop->getFitness(bestFitIndex);

        if (bestFitVal < globalBest.fitness)
        {
            util::copyArray<T>(p.mainPop->getPopulationPtr(bestFitIndex), globalBest.vector, dimSize);
            globalBest.fitness = bestFitVal;
        }

        if (p.bestFitnessTable != nullptr)
            p.bestFitnessTable->setEntry(iter, p.fitTableCol, globalBest.fitness);

        if (p.worstFitnessTable != nullptr)
            p.worstFitnessTable->setEntry(iter, p.fitTableCol, p.mainPop->getWorstFitness());

        if (!p.popFile.empty())
            p.mainPop->outputPopulationCsv(util::s_replace(p.popFile, std::string(POPFILE_GEN_PATTERN), std::to_string(iter)));
    }

    return 0;
}

template <class T>
void mfunc::ParticleSwarm<T>::updateParticle(PSParams<T>& p, const Particle<T>& globalBest, T** velMatrix, size_t pIndex)
{
    const size_t dimSize = p.mainPop->getDimensionsSize();
    auto pBestVector = p.pbPop->getPopulationPtr(pIndex);
    auto curVector = p.mainPop->getPopulationPtr(pIndex);

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

    if (newFitness < pbFitness)
    {
        p.pbPop->copyFrom(p.mainPop, pIndex, pIndex);
    }
}

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