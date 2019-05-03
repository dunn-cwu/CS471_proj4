#ifndef __GENETICALG_H
#define __GENETICALG_H

#include <vector>
#include "population.h"
#include "mfuncptr.h"
#include "datatable.h"
#include "random"
#include "mem.h"

namespace mfunc
{
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
    };

    template <class T>
    class GeneticAlgorithm
    {
    public:
        GeneticAlgorithm();
        ~GeneticAlgorithm() = default;
        int run(GAParams<T> params);
        
    private:
        std::vector<T> childOne;
        std::vector<T> childTwo;
        std::random_device seed;
        std::mt19937 engine;
        std::uniform_real_distribution<double> rchance;

        void select(mdata::Population<T>* pop, size_t& outP1, size_t& outP2);
        size_t selRW(mdata::Population<T>* pop);
        void crossover(size_t dim, T* p1, T* p2, double cr, std::vector<T>& outCh1, std::vector<T>& outCh2);
        void mutate(size_t dim, std::vector<T>& s, double mutProb, double mutRange, double mutPrec, T fMin, T fMax);
        void reduce(mdata::Population<T>* oldPop, mdata::Population<T>* newPop, size_t elitism);
    };
}

template <class T>
mfunc::GeneticAlgorithm<T>::GeneticAlgorithm()
    : seed(), engine(seed()), rchance(0, 1)
{
}

template <class T>
int mfunc::GeneticAlgorithm<T>::run(GAParams<T> p)
{
    if (p.mainPop == nullptr || p.auxPop == nullptr || p.fPtr == nullptr)
        return 1;

    if (p.mainPop->getPopulationSize() != p.auxPop->getPopulationSize() ||
        p.mainPop->getDimensionsSize() != p.auxPop->getDimensionsSize())
        return 2;

    const size_t popSize = p.mainPop->getPopulationSize();
    const size_t dimSize = p.mainPop->getDimensionsSize();
    const size_t elitism = p.elitismRate * (double)popSize;

    childOne.resize(dimSize);
    childTwo.resize(dimSize);

    p.mainPop->setFitnessNormalization(true);
    p.auxPop->setFitnessNormalization(true);
    p.mainPop->generate(p.fMinBound, p.fMaxBound);
    p.mainPop->calcAllFitness(p.fPtr);

    for (unsigned int gen = 0; gen < p.generations; gen++)
    {
        for (size_t s = 0; s < popSize;)
        {
            size_t p1Index = 0;
            size_t p2Index = 0;

            select(p.mainPop, p1Index, p2Index);
            crossover(dimSize, p.mainPop->getPopulationPtr(p1Index), p.mainPop->getPopulationPtr(p2Index), 
                p.crProb, childOne, childTwo);

            mutate(dimSize, childOne, p.mutProb, p.mutRange, p.mutPrec, p.fMinBound, p.fMaxBound);
            mutate(dimSize, childTwo, p.mutProb, p.mutRange, p.mutPrec, p.fMinBound, p.fMaxBound);

            p.auxPop->copyPopulation(s, childOne);
            s++;
            if (s < popSize) p.auxPop->copyPopulation(s, childTwo);
            s++;
        }

        p.auxPop->calcAllFitness(p.fPtr);
        reduce(p.mainPop, p.auxPop, elitism);

        // Swap the two populations
        auto tmp = p.mainPop;
        p.mainPop = p.auxPop;
        p.auxPop = tmp;

        p.fitnessTable->setEntry(gen, p.fitTableCol, p.mainPop->getMinCost());
    }

    return 0;
}

template <class T>
void mfunc::GeneticAlgorithm<T>::select(mdata::Population<T>* pop, size_t& outP1, size_t& outP2)
{
    outP1 = selRW(pop);
    outP2 = selRW(pop);
}

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

template <class T>
void mfunc::GeneticAlgorithm<T>::crossover(size_t dim, T* p1, T* p2, double cr, std::vector<T>& outCh1, std::vector<T>& outCh2)
{
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
void mfunc::GeneticAlgorithm<T>::mutate(size_t dim, std::vector<T>& s, double mutProb, double mutRange, double mutPrec, T fMin, T fMax)
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