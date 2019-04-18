/**
 * @file population.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the Population class.
 * Stores a population and fitness values.
 * Includes functions to analyze the fitness data.
 * @version 0.1
 * @date 2019-04-04
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef __POPULATION_H
#define __POPULATION_H

#include <cstddef> // size_t definition
#include <vector>
#include <random>
#include <ostream>

namespace mdata
{
    template<class T>
    class VectFitPair
    {
    public:
        VectFitPair() : vArr(nullptr), realSize(0) {}

        VectFitPair(size_t vSize) : vArr(nullptr), realSize(0)
        {
            realSize = vSize + 1;
            alloc(realSize);
            init(0);
        }

        VectFitPair(const T* vector, size_t vSize) : VectFitPair(vector, 0, vSize)
        {
        }

        VectFitPair(const T* vector, T fitness, size_t vSize) : vArr(nullptr), realSize(0)
        {
            realSize = vSize + 1;
            alloc(realSize);
            vArr[0] = fitness;
            loadVector(vector);
        }

        ~VectFitPair()
        {
            dealloc();
        }
    
        void init(T val = 0)
        {
            if (vArr == nullptr) return;

            for (size_t i = 0; i < realSize; i++)
                vArr[i] = val;
        }

        size_t vectorSize()
        {
            return realSize - 1;
        }

        T& fitness()
        {
            return vArr[0];
        }
        
        void loadVector(const T*& vOther)
        {
            if (vArr == nullptr) return;

            for (size_t i = 1; i < realSize; i++)
            {
                vArr[i-1] = vOther[i-1];
            }
        }

        const T*& vector()
        {
            return &(vArr[1]);
        }

        T& operator[](const size_t& vectInd)
        {
            return vArr[vectInd + 1];
        }

        // Copy constructor
        VectFitPair(const VectFitPair& other)
        {
            realSize = other.realSize;
            alloc(other.realSize);
            for (size_t i = 0; i < realSize; i++)
                vArr[i] = other.vArr[i];
        }

        // Copy assignment
        VectFitPair& operator=(const VectFitPair& other)
        {
            realSize = other.realSize;
            alloc(other.realSize);
            for (size_t i = 0; i < realSize; i++)
                vArr[i] = other.vArr[i];

            return *this;
        }

        // Move Constructor
        VectFitPair(VectFitPair&& other)
        {
            dealloc();
            realSize = other.realSize;
            vArr = other.vArr;
            other.vArr = nullptr;
            other.realSize = 0;
        }

        // Move assignment
        VectFitPair& operator=(VectFitPair&& other)
        {
            dealloc();
            realSize = other.realSize;
            vArr = other.vArr;
            other.vArr = nullptr;
            other.realSize = 0;

            return *this;
        }
    private:
        size_t realSize;
        T* vArr;

        void alloc(size_t size)
        {
            dealloc();
            vArr = new T[size];
        }

        void dealloc()
        {
            if (vArr == nullptr) return;
            delete[] vArr;
            vArr = nullptr;
        }
    };

    /**
     * @brief Data class for storing a multi-dimensional population of data.
     * Includes fitness analysis functions.
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
        size_t getBestSize();
        T* getPopulation(size_t popIndex);

        bool generate(T minBound, T maxBound);
        bool setFitness(size_t popIndex, T value);
        bool setFitness(size_t popIndex, T (*f)(T*, size_t v));

        void storeBest();
        VectFitPair<T>& getBest(size_t index);
        void clearBest();

        T getFitnessValue(size_t popIndex);
        std::vector<T> getAllFitness();
        std::vector<T> getAllBestFitness();
        T getFitnessAverage();
        T getFitnessStandardDev();
        T getFitnessRange();
        T getFitnessMedian();

        void outputPopulation(std::ostream& outStream, const char* delim, const char* lineBreak);
        void outputFitness(std::ostream& outStream, const char* delim, const char* lineBreak);
    private:
        const size_t popSize; /** Size of the population, and the number of rows in the popMatrix */
        const size_t popDim;  /** Dimensions of the population, and the number of columns in the popMatrix */
        
        T** popMatrix; /** Matrix of population values */
        T* popFitness; /** Array of fitness values */
        std::vector<VectFitPair<T>> bestPop;

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
