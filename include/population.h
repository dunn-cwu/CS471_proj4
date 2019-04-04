#ifndef __POPULATION_H
#define __POPULATION_H

#include <cstddef> // size_t definition
#include <random>

namespace mdata
{
    template<class T>
    class Population
    {
    public:
        Population(size_t popSize, size_t dimensions);
        ~Population();

        bool generate(T minBound, T maxBound);
    protected:
        const size_t pop;
        const size_t dim;
        T** popMatrix;

        bool allocPopMatrix();
        void releasePopMatrix();
    };
}

#endif