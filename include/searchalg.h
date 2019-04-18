#ifndef __SEARCHALG_H
#define __SEARCHALG_H

#include <chrono>
#include "population.h"
#include "testresult.h"

using namespace std::chrono;

namespace mdata
{
    template<class T>
    class SearchAlgorithm
    {
    public:
        SearchAlgorithm() : timeDiff(0.0) {}
        virtual TestResult run(const unsigned int funcId, const T fMin, const T fMax, Population<T>* const pop, const size_t iterations, const T alpha) = 0;
    protected:
        double timeDiff;
        high_resolution_clock::time_point timer;

        inline void startTimer()
        {
            timer = high_resolution_clock::now();
        }

        inline double stopTimer()
        {
            high_resolution_clock::time_point t_end = high_resolution_clock::now();
            return (double)duration_cast<nanoseconds>(t_end - timer).count() / 1000000.0;
        }
    };
}

#endif