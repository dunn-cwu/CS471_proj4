/**
 * @file population.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the Population class.
 * Stores a population and fitness values.
 * @version 0.2
 * @date 2019-04-04
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "population.h"
#include "mem.h"
#include <new>
#include <fstream>

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
Population<T>::Population(size_t pSize, size_t dimensions) 
    : popMatrix(nullptr), popSize(pSize), popDim(dimensions), rdev(), rgen(rdev())
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

template <class T>
T* Population<T>::getBestPopulationPtr()
{
    return getPopulationPtr(getBestFitnessIndex());
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

    // Set up a uniform distribution for the random number generator with the correct function bounds
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

template <class T>
bool Population<T>::generateSingle(size_t popIndex, T minBound, T maxBound)
{
    if (popMatrix == nullptr || popIndex >= popSize) return false;

    // Set up a uniform distribution for the random number generator with the correct function bounds
    std::uniform_real_distribution<double> dist((double)minBound, (double)maxBound);

    for (size_t d = 0; d < popDim; d++)
    {
        T rand = (T)dist(rgen);
        popMatrix[popIndex][d] = rand;
    }

    popFitness[popIndex] = 0;

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

/**
 * @brief Uses the given function pointer to update the fitness value for the
 * population vector at the given index.
 * 
 * @tparam T Data type of the population.
 * @param popIndex Index of the population vector you wish to set the fitness for.
 * @param funcPtr Function pointer to the math function that will be used to calculate 
 * the fitness value.
 * @return Returns true on success, otherwise false.
 */
template<class T>
bool Population<T>::calcFitness(size_t popIndex, mfunc::mfuncPtr<T> funcPtr)
{
    if (popFitness == nullptr || popIndex >= popSize) return false;

    popFitness[popIndex] = funcPtr(popMatrix[popIndex], popDim);

    return true;
}

template<class T>
bool Population<T>::calcAllFitness(mfunc::mfuncPtr<T> funcPtr)
{
    for (size_t i = 0; i < popSize; i++)
    {
        if (!calcFitness(i, funcPtr))
            return false;
    }

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

/**
 * @brief Returns a pointer to the current best fitness value
 * 
 * @tparam T Data type of the population.
 * @return T* Pointer to the best fitness value
 */
template<class T>
T* Population<T>::getBestFitnessPtr()
{
    return &popFitness[getBestFitnessIndex()];
}

/**
 * @brief Returns the index of the current best fitness value
 * 
 * @tparam T Data type of the population.
 * @return size_t Index of the best fitness value
 */
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

template<class T>
T Population<T>::getBestFitness()
{
    return getFitness(getBestFitnessIndex());
}

template<class T>
size_t Population<T>::getWorstFitnessIndex()
{
    size_t worstIndex = 0;

    for (size_t i = 1; i < popSize; i++)
    {
        if (popFitness[i] > popFitness[worstIndex])
            worstIndex = i;
    }

    return worstIndex;
}

template<class T>
T Population<T>::getWorstFitness()
{
    return getFitness(getWorstFitnessIndex());
}

template<class T>
void Population<T>::sortFitnessAscend()
{
    qs_fit_ascend(0, popSize - 1);
}

template<class T>
void Population<T>::sortFitnessDescend()
{
    qs_fit_descend(0, popSize - 1);
}

template<class T>
bool Population<T>::copyFrom(Population<T>* srcPtr, size_t srcIndex, size_t destIndex)
{
    if (srcPtr == nullptr) return false;

    const size_t srcDim = srcPtr->getDimensionsSize();
    if (srcDim != popDim) return false;

    T* srcVector = srcPtr->getPopulationPtr(srcIndex);
    T* destVector = getPopulationPtr(destIndex);

    if (srcVector == nullptr || destVector == nullptr) return false;

    copyArray<T>(srcVector, destVector, popDim);
    setFitness(destIndex, srcPtr->getFitness(srcIndex));

    return true;
}

template<class T>
bool Population<T>::copyAllFrom(Population<T>* srcPtr)
{
    if (srcPtr == nullptr) return false;

    const size_t srcSize = srcPtr->getPopulationSize();
    const size_t srcDim = srcPtr->getDimensionsSize();

    if (srcSize != popSize || srcDim != popDim)
        return false;

    for (size_t i = 0; i < popSize; i++)
    {
        if (!copyFrom(srcPtr, i, i))
            return false;
    }

    return true;
}

template<class T>
bool Population<T>::copyPopulation(T* src, size_t destIndex)
{
    T* destVect = getPopulationPtr(destIndex);
    if (destVect == nullptr) 
        return false;

    for (size_t i = 0; i < popDim; i++)
    {
        destVect[i] = src[i];
    }

    return true;
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

template<class T>
bool Population<T>::outputPopulationCsv(std::string filePath)
{
    static const char* delim = ",";
    static const char* newline = "\n";

    std::ofstream file;
    file.open(filePath, std::ios::out | std::ios::trunc);
    if (!file.good()) return false;

    outputPopulation(file, delim, newline);
    file.close();

    return true;
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

// ==========================================
// Quicksort Implementation modified from:
// https://www.geeksforgeeks.org/quick-sort/
// ==========================================

template <class T>
void Population<T>::qs_swapval(T& a, T& b) 
{ 
    T t = a;
    a = b;
    b = t; 
}

template <class T>
void Population<T>::qs_swapptr(T*& a, T*& b) 
{ 
    T* t = a;
    a = b; 
    b = t; 
}

template <class T>
long Population<T>::part_fit_ascend(long low, long high) 
{ 
    T pivot = popFitness[high]; // pivot 
    long i = (low - 1);  // Index of smaller element 
  
    for (long j = low; j <= high- 1; j++) 
    { 
        if (popFitness[j] <= pivot) 
        { 
            i++;    // increment index of smaller element 
            qs_swapval(popFitness[i], popFitness[j]);
            qs_swapptr(popMatrix[i], popMatrix[j]);
        } 
    } 
    qs_swapval(popFitness[i + 1], popFitness[high]); 
    qs_swapptr(popMatrix[i + 1], popMatrix[high]);

    return (i + 1); 
} 
  
template <class T>
void Population<T>::qs_fit_ascend(long low, long high) 
{ 
    if (low < high) 
    {
        long pi = part_fit_ascend(low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        qs_fit_ascend(low, pi - 1); 
        qs_fit_ascend(pi + 1, high);
    } 
}

template <class T>
long Population<T>::part_fit_descend(long low, long high) 
{ 
    T pivot = popFitness[high]; // pivot 
    long i = (low - 1);  // Index of smaller element 
  
    for (long j = low; j <= high- 1; j++) 
    { 
        if (popFitness[j] > pivot) 
        { 
            i++;    // increment index of smaller element 
            qs_swapval(popFitness[i], popFitness[j]);
            qs_swapptr(popMatrix[i], popMatrix[j]);
        } 
    } 
    qs_swapval(popFitness[i + 1], popFitness[high]); 
    qs_swapptr(popMatrix[i + 1], popMatrix[high]);

    return (i + 1); 
} 
  
template <class T>
void Population<T>::qs_fit_descend(long low, long high) 
{ 
    if (low < high) 
    {
        long pi = part_fit_descend(low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        qs_fit_descend(low, pi - 1); 
        qs_fit_descend(pi + 1, high);
    } 
} 

// Explicit template specializations due to separate implementations in this CPP file
template class mdata::Population<float>;
template class mdata::Population<double>;
template class mdata::Population<long double>;

// =========================
// End of population.cpp
// =========================