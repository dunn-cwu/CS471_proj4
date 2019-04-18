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
#include "mfunc.h"
#include "inireader.h"
#include "population.h"
#include "threadpool.h"
#include "testresult.h"

#define DEFAULT_NUM_THREADS 1

namespace cs471
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
    class mfuncExperiment
    {
    public:
        mfuncExperiment();
        ~mfuncExperiment();
        bool init(const char* paramFile);
        int runAllFunc();
        mdata::TestResult runFunc(unsigned int funcId);
    private:
        util::IniReader iniParams; /** IniReader class instance for importing experiment parameters */
        std::string resultsFile;   /** The file path for the results output *.csv file */
        mdata::Population<double>** populations; /** Array of population objects that contain matrices and fitness arrays for each function */
        RandomBounds<double>* vBounds; /** An array of RandomBounds structs that holds the function bounds read from iniParams */
        ThreadPool* tPool;
        size_t iterations;
        bool outputPop; /** If set to true, all population data will be exported to files */
        bool outputFitness; /** If set to true, all fitness data will be exported to files */

        bool genFuncVectors(unsigned int funcId);

        bool parseFuncBounds();

        void exportPop(unsigned int func);

        bool allocatePopulations(size_t popSize, size_t dimensions);
        void releasePopulations();

        bool allocateVBounds();
        void releaseVBounds();

        bool allocateThreadPool(size_t numThreads);
        void releaseThreadPool();
    };
} // proj1

#endif

// =========================
// End of proj1.h
// =========================
