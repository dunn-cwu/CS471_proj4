/**
 * @file searchalg.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Defines the SearchAlgorithm class, Algorithm enum,
 * and  AlgorithmNames struct. The SearchAlgorithm class serves
 * as a base class for implemented search algorithms.
 * @version 0.1
 * @date 2019-04-19
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef __SEARCHALG_H
#define __SEARCHALG_H

#include <chrono>
#include "population.h"
#include "testresult.h"
#include "mfuncptr.h"

using namespace std::chrono;

namespace enums
{
    /**
     * @brief Enum of different available search algorithms
     */
    enum Algorithm
    {
        BlindSearch = 0,
        LocalSearch = 1,
        Count = 2
    };

    /**
     * @brief Struct that contains constant string names for
     * the different search algorithms
     */
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
    /**
     * @brief The SearchAlgorithm class is used as a base class for
     * other implemented search algorithms. Provides a common interface
     * to run each algorithm.
     * 
     * @tparam T The data type used by the algorithm
     */
    template<class T>
    class SearchAlgorithm
    {
    public:
        SearchAlgorithm() : timeDiff(0.0) {}
        virtual ~SearchAlgorithm() = 0;
        virtual TestResult<T> run(mfunc::mfuncPtr<T> funcPtr, const T fMin, const T fMax, Population<T>* const pop, const T alpha) = 0;
    protected:
        double timeDiff;
        high_resolution_clock::time_point timer;

        /**
         * @brief Starts the execution time timer
         */
        void startTimer()
        {
            timer = high_resolution_clock::now();
        }

        /**
         * @brief Returns the amount of time that has passed since startTimer() was called in miliseconds.
         */
        double stopTimer()
        {
            high_resolution_clock::time_point t_end = high_resolution_clock::now();
            return static_cast<double>(duration_cast<nanoseconds>(t_end - timer).count()) / 1000000.0;
        }
    };
}

// Trivial implementation of pure-virtual destructor
// as required by the C++ language
template<class T>
mdata::SearchAlgorithm<T>::~SearchAlgorithm() { }

#endif

// =========================
// End of searchalg.h
// =========================
