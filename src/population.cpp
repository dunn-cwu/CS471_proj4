/**
 * @file population.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the Population class.
 * Stores a population and fitness values.
 * Includes functions to analyze the fitness data.
 * @version 0.1
 * @date 2019-04-04
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "population.h"
#include "mem.h"
#include <new>

using namespace mdata;
using namespace util;

/**
 * @brief Construct a new Population object.
 * 
 * @tparam T Data type of the population.
 * @param pSize Size of the population.
 * @param dimensions Dimensions of the population.
 */
template <class T>
Population<T>::Population(size_t pSize, size_t dimensions) : popMatrix(nullptr), popSize(pSize), popDim(dimensions)
{
    if (!allocPopMatrix() || !allocPopFitness())
        throw std::bad_alloc();
}

/**
 * @brief Destroy Population object.
 * 
 * @tparam T Data type of the population.
 */
template <class T>
Population<T>::~Population()
{
    releasePopMatrix();
    releasePopFitness();
}

/**
 * @brief Returns true if the population instance is allocated
 * and ready to be used.
 * 
 * @tparam T Data type of the population.
 * @return Returns true if the population instance is in a valid state.
 */
template <class T>
bool Population<T>::isReady()
{
    return popMatrix != nullptr && popFitness != nullptr;
}

/**
 * @brief Returns the size of the population.
 * 
 * @tparam T Data type of the population.
 * @return The size of the population.
 */
template <class T>
size_t Population<T>::getPopulationSize()
{
    return popSize;
}

/**
 * @brief Returns the dimensions of the population.
 * 
 * @tparam T Data type of the population.
 * @return The number of dimensions in the population.
 */
template <class T>
size_t Population<T>::getDimensionsSize()
{
    return popDim;
}

/**
 * @brief Returns an array for the population with the given index.
 * 
 * @tparam T Data type of the population.
 * @param popIndex Index of the population vector you wish to retrieve.
 * @return Pointer to population vector array at the given index.
 */
template <class T>
T* Population<T>::getPopulationPtr(size_t popIndex)
{
    if (popFitness == nullptr || popIndex >= popSize) return nullptr;
    
    return popMatrix[popIndex];
}

/**
 * @brief Generates new random values for this population
 * that are within the given bounds. Resets all fitness values
 * to zero.
 * 
 * @tparam T Data type of the population.
 * @param minBound The minimum bound for a population value.
 * @param maxBound The maximum bound for a population value.
 * @return Returns true of the population was succesfully generated, otherwise false.
 */
template <class T>
bool Population<T>::generate(T minBound, T maxBound)
{   
    if (popMatrix == nullptr) return false;

    // Generate a new seed for the mersenne twister engine
    rgen = std::mt19937(rdev());

    // Set up a normal (bell-shaped) distribution for the random number generator with the correct function bounds
    std::uniform_real_distribution<double> dist((double)minBound, (double)maxBound);

    // Generate values for all vectors in popMatrix
    for (size_t s = 0; s < popSize; s++)
    {
        for (size_t d = 0; d < popDim; d++)
        {
            T rand = (T)dist(rgen);
            popMatrix[s][d] = rand;
        }
    }

    // Reset popFitness values to 0
    initArray<T>(popFitness, popSize, (T)0.0);

    return true;
}

/**
 * @brief Sets the fitness value for a specific population vector index.
 * 
 * @tparam T Data type of the population.
 * @param popIndex Index of the population vector you wish to set the fitness for.
 * @param value The value of the fitness.
 * @return Returns true if the fitness was succesfully set, otherwise false.
 */
template<class T>
bool Population<T>::setFitness(size_t popIndex, T value)
{
    if (popFitness == nullptr || popIndex >= popSize) return false;

    popFitness[popIndex] = value;

    return true;
}

template<class T>
bool Population<T>::calcFitness(size_t popIndex, mfunc::mfuncPtr<T> funcPtr)
{
    if (popFitness == nullptr || popIndex >= popSize) return false;

    popFitness[popIndex] = funcPtr(popMatrix[popIndex], popDim);

    return true;
}

/**
 * @brief Returns the fitness value for a specific population vector index.
 * 
 * @tparam T Data type of the population.
 * @param popIndex Index of the population vector you wish to retrieve the fitness from.
 * @return Returns the fitness value if popIndex is valid. Otherwise zero.
 */
template<class T>
T Population<T>::getFitness(size_t popIndex)
{
    if (popFitness == nullptr || popIndex >= popSize) return 0;

    return popFitness[popIndex];
}

/**
 * @brief Returns the fitness value for a specific population vector index.
 * 
 * @tparam T Data type of the population.
 * @param popIndex Index of the population vector you wish to retrieve the fitness from.
 * @return Returns the fitness value if popIndex is valid. Otherwise zero.
 */
template<class T>
T* Population<T>::getFitnessPtr(size_t popIndex)
{
    if (popFitness == nullptr || popIndex >= popSize) return 0;

    return &popFitness[popIndex];
}

template<class T>
std::vector<T> Population<T>::getAllFitness()
{
    return std::vector<T>(popFitness[0], popFitness[popSize]);
}

template<class T>
T* Population<T>::getBestFitnessPtr()
{
    return &popFitness[getBestFitnessIndex()];
}

template<class T>
size_t Population<T>::getBestFitnessIndex()
{
    size_t bestIndex = 0;

    for (size_t i = 1; i < popSize; i++)
    {
        if (popFitness[i] < popFitness[bestIndex])
            bestIndex = i;
    }

    return bestIndex;
}

/**
 * @brief Outputs all population data to the given output stream.
 * 
 * @tparam T Data type of the population.
 * @param outStream Output stream to write the data to.
 * @param delim Delimiter characters to separate columns.
 * @param lineBreak Delimiter characters to separate rows.
 */
template<class T>
void Population<T>::outputPopulation(std::ostream& outStream, const char* delim, const char* lineBreak)
{
    if (popMatrix == nullptr) return;

    for (size_t j = 0; j < popSize; j++)
    {
        for (size_t k = 0; k < popDim; k++)
        {
            outStream << popMatrix[j][k];
            if (k < popDim - 1)
                outStream << delim;
        }

        outStream << lineBreak;
    }
}

/**
 * @brief Outputs all fitness data to the given output stream.
 * 
 * @tparam T Data type of the population.
 * @param outStream Output stream to write the data to.
 * @param delim Delimiter characters to separate columns.
 * @param lineBreak Delimiter characters to separate rows.
 */
template<class T>
void Population<T>::outputFitness(std::ostream& outStream, const char* delim, const char* lineBreak)
{
    if (popFitness == nullptr) return;

    for (size_t j = 0; j < popSize; j++)
    {
        outStream << popFitness[j];
            if (j < popSize - 1) 
                outStream << delim;
    }

    if (lineBreak != nullptr)
        outStream << lineBreak;
}

/**
 * @brief Helper function that allocates the population matrix.
 * 
 * @tparam T Data type of the population.
 * @return Returns true if the population matrix was succesfully allocated, otherwise false.
 */
template <class T>
bool Population<T>::allocPopMatrix()
{
    if (popSize == 0 || popDim == 0) return false;

    popMatrix = allocMatrix<T>(popSize, popDim);
    initMatrix<T>(popMatrix, popSize, popDim, 0);

    return popMatrix != nullptr;
}

/**
 * @brief Helper function that releases the population matrix memory.
 * 
 * @tparam T Data type of the population.
 */
template <class T>
void Population<T>::releasePopMatrix()
{
    releaseMatrix<T>(popMatrix, popSize);
}

/**
 * @brief Helper function that allocates the population fitness array.
 * 
 * @tparam T Data type of the population.
 * @return Returns true if the fitness array was succesfully allocated, otherwise false.
 */
template <class T>
bool Population<T>::allocPopFitness()
{
    if (popSize == 0 || popDim == 0) return false;

    popFitness = allocArray<T>(popSize);
    initArray<T>(popFitness, popSize, 0);

    return popFitness != nullptr;
}

/**
 * @brief Helper function that releases the population fitness array memory.
 * 
 * @tparam T Data type of the population.
 */
template <class T>
void Population<T>::releasePopFitness()
{
    releaseArray<T>(popFitness);
}

template class mdata::Population<float>;
template class mdata::Population<double>;
template class mdata::Population<long double>;

// =========================
// End of population.cpp
// =========================
