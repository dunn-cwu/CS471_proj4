#ifndef __SEARCHALG_H
#define __SEARCHALG_H

#include <chrono>
#include "population.h"
#include "testresult.h"
#include "mfuncPtr.h"

using namespace std::chrono;

namespace enums
{
    enum Algorithm
    {
        BlindSearch = 0,
        LocalSearch = 1,
        Count = 2
    };

    struct AlgorithmNames
    {
        static constexpr const char* BLIND_SEARCH = "Blind Search";
        static constexpr const char* LOCAL_SEARCH = "Local Search";

        static const char* get(Algorithm alg)
        {
            switch (alg)
            {
                case Algorithm::BlindSearch:
                    return BLIND_SEARCH;
                case Algorithm::LocalSearch:
                    return LOCAL_SEARCH;
                default:
                    return "";
                    break;
            }
        }
    };
}

namespace mdata
{
    template<class T>
    class SearchAlgorithm
    {
    public:
        SearchAlgorithm() : timeDiff(0.0) {}
        virtual TestResult<T> run(mfunc::mfuncPtr<T> funcPtr, const T fMin, const T fMax, Population<T>* const pop, const T alpha) = 0;
    protected:
        double timeDiff;
        high_resolution_clock::time_point timer;

        void startTimer()
        {
            timer = high_resolution_clock::now();
        }

        double stopTimer()
        {
            high_resolution_clock::time_point t_end = high_resolution_clock::now();
            return (double)duration_cast<nanoseconds>(t_end - timer).count() / 1000000.0;
        }
    };
}

#endif