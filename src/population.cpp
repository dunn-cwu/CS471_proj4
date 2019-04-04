#include "population.h"
#include "mem.h"
#include <new>

using namespace mdata;
using namespace util;

template <class T>
Population<T>::Population(size_t popSize, size_t dimensions) : popMatrix(nullptr), pop(popSize), dim(dimensions)
{
    if (!allocPopMatrix())
        throw std::bad_alloc();
}

template <class T>
Population<T>::~Population()
{
    releasePopMatrix();
}

template <class T>
bool Population<T>::generate(T minBound, T maxBound)
{

}

template <class T>
bool Population<T>::allocPopMatrix()
{
    if (pop == 0 || dim == 0) return false;

    popMatrix = allocMatrix<T>(pop, dim);
    return popMatrix != nullptr;
}

template <class T>
void Population<T>::releasePopMatrix()
{
    releaseMatrix<T>(popMatrix, pop);
}