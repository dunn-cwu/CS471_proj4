/**
 * @file experiment.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the Experiment class.
 * Contains the basic logic and functions to run the cs471 project experiment.
 * @version 0.2
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __EXPERIMENT_H
#define __EXPERIMENT_H

#include <string>
#include <random>
#include <chrono>
#include <vector>
#include "mfunctions.h"
#include "inireader.h"
#include "population.h"
#include "threadpool.h"
#include "testresult.h"
#include "testparam.h"
#include "geneticalg.h"

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

    enum class Algorithm
    {
        GeneticAlgorithm = 0,
        DifferentialEvolution = 1,
        Count = 2
    };

    /**
     * @brief Contains classes for running the CS471 project experiment.
     * 
     * The Experiment class opens a given parameter .ini file
     * and executes the CS471 project 2 experiment with the specified
     * parameters. runAllFunc() runs all 18 functions defined in mfunctions.h
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
        int testAllFunc_GA();
        int runGAThreaded(GAParams<T> gaParams);
    private:
        std::mutex popPoolMutex;
        util::IniReader iniParams; /** IniReader class instance for importing experiment parameters */
        std::vector<mdata::Population<T>*> populationsPool; /** Pool of population objects used by the worker threads. 1 per thread. */
        std::string resultsFile;   /** The file path for the results output *.csv file */
        std::string execTimesFile;   /** The file path for the exec times output *.csv file */
        RandomBounds<T>* vBounds; /** An array of RandomBounds structs that holds the function bounds read from iniParams */
        ThreadPool* tPool; /** Pool of worker threads which are used to run multiple tests in parallel */
        size_t iterations; /** Number of iterations for the selected test algorithm */
        Algorithm testAlg;

        bool loadGAParams(GAParams<T>& refParams);

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
} // mfunc

#endif

// =========================
// End of experiment.h
// =========================
