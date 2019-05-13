/**
 * @file population.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the Population class.
 * Stores a population and resulting fitness values.
 * @version 0.2
 * @date 2019-04-04
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef __POPULATION_H
#define __POPULATION_H

#include <cstddef> // size_t definition
#include <random>
#include <ostream>
#include "mfuncptr.h"

namespace mdata
{
    /**
     * @brief Data class for storing a multi-dimensional population of data
     * with the associated fitness.
     * 
     * @tparam T Data type of the population.
     */
    template<class T>
    class Population
    {
    public:
        Population(size_t popSize, size_t dimensions);
        ~Population();

        bool isReady();
        size_t getPopulationSize();
        size_t getDimensionsSize();
        T* getPopulationPtr(size_t popIndex);
        T* getBestPopulationPtr();

        bool generate(T minBound, T maxBound);
        bool setFitness(size_t popIndex, T value);
        bool calcFitness(size_t popIndex, mfunc::mfuncPtr<T> funcPtr);
        bool calcAllFitness(mfunc::mfuncPtr<T> funcPtr);

        T getFitness(size_t popIndex);
        T* getFitnessPtr(size_t popIndex);

        T* getBestFitnessPtr();
        size_t getBestFitnessIndex();
        T getBestFitness();

        bool copyFrom(Population<T>* srcPtr, size_t srcIndex, size_t destIndex);
        bool copyAllFrom(Population<T>* srcPtr);

        void outputPopulation(std::ostream& outStream, const char* delim, const char* lineBreak);
        void outputFitness(std::ostream& outStream, const char* delim, const char* lineBreak);
    private:
        const size_t popSize; /** Size of the population, and the number of rows in the popMatrix */
        const size_t popDim;  /** Dimensions of the population, and the number of columns in the popMatrix */
        
        T** popMatrix; /** Matrix of population values */
        T* popFitness; /** Array of fitness values */

        std::random_device rdev; /** Random seed for random number generator */
        std::mt19937 rgen; /** Mersenne twister random number generator engine */

        bool allocPopMatrix();
        void releasePopMatrix();

        bool allocPopFitness();
        void releasePopFitness();
    };
}

#endif

// =========================
// End of population.h
// =========================