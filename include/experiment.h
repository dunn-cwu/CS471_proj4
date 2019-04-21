/**
 * @file cs471.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the mfuncExperiment class.
 * Contains the basic logic and functions to run the cs471 project experiment.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __CS471_H
#define __CS471_H

#include <string>
#include <random>
#include <chrono>
#include <vector>
#include "mfunctions.h"
#include "inireader.h"
#include "population.h"
#include "threadpool.h"
#include "searchalg.h"
#include "testresult.h"
#include "testparam.h"

namespace enums
{
    enum class ThreadingMode
    {
        ByFunction = 0,
        ByFunctionIteration = 1,
        Count = 2
    };
}

namespace mfunc
{
    /**
     * @brief Simple struct for storing the minimum
     * and maximum input vector bounds for a function
     */
    template<class T>
    struct RandomBounds
    {
        T min = 0.0;
        T max = 0.0;
    };

    /**
     * @brief Contains classes for running the CS471 project experiment.
     * 
     * The mfuncExperiment class opens a given parameter .ini file
     * and executes the CS471 project 1 experiment with the specified
     * parameters. runAllFunc() runs all 18 functions defined in mfunc.cpp
     * a given number of times with vectors of random values that have a
     * given number of dimensions and collects all results/data. This
     * data is then entered into a DataTable and exported as a *.csv file.
     */
    template<class T>
    class Experiment
    {
    public:
        Experiment();
        ~Experiment();
        bool init(const char* paramFile);
        int testAllFunc();
        int testFunc(mdata::TestParameters<T>* tParams);
        int testFuncThreaded(mdata::TestParameters<T>* tParams);
    private:
        std::mutex popPoolMutex;
        util::IniReader iniParams; /** IniReader class instance for importing experiment parameters */
        std::vector<mdata::Population<T>*> populationsPool;
        std::string resultsFile;   /** The file path for the results output *.csv file */
        std::string execTimesFile;   /** The file path for the exec times output *.csv file */
        RandomBounds<T>* vBounds; /** An array of RandomBounds structs that holds the function bounds read from iniParams */
        ThreadPool* tPool;
        size_t iterations;
        T alpha;
        enums::Algorithm testAlg;
        enums::ThreadingMode threadMode;

        mdata::TestResult<T> asyncAlgIteration(mdata::TestParameters<T>* tParams, mdata::SearchAlgorithm<T>* alg);
        mdata::Population<T>* popPoolRemove();
        void popPoolAdd(mdata::Population<T>* popPtr);

        bool parseFuncBounds();

        bool allocatePopulationPool(size_t count, size_t popSize, size_t dimensions);
        void releasePopulationPool();

        bool allocateVBounds();
        void releaseVBounds();

        bool allocateThreadPool(size_t numThreads);
        void releaseThreadPool();
    };
} // proj1

#endif

// =========================
// End of experiment.h
// =========================
