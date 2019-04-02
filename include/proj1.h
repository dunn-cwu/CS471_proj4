/**
 * @file proj1.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the basic logic and functions to run
 * the cs471 project 1 experiment.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __PROJ1_H
#define __PROJ1_H

#include <string>
#include <random>
#include <chrono>
#include <vector>
#include "mfunc.h"
#include "inireader.h"

namespace proj1
{
    /**
     * @brief Simple struct for storing the minimum
     * and maximum input vector bounds for a function
     */
    struct RandomBounds
    {
        double min = 0.0;
        double max = 0.0;
    };

    /**
     * @brief Contains classes for running the CS471 project 1 experiment.
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
        int runFunc(unsigned int funcId, std::vector<double>& resultArrOut, double& timeOut);
    private:
        util::IniReader iniParams; /** IniReader class instance for importing experiment parameters */
        std::string resultsFile;   /** The file path for the output *.csv file */
        size_t nbrDim; /** The number of dimensions for the function vectors which is read from iniParams */
        size_t nbrSol; /** The number of solutions to generate for each function, which is read from iniParams */
        double** vMatrix; /** A two dimensional array of size nbrSol * nbrDim which stores function vectors */
        RandomBounds* vBounds; /** An array of RandomBounds structs that holds the function bounds read from iniParams */

        std::random_device rdev; /** Random seed for random number generator */
        std::mt19937 rgen; /** Mersenne twister random number generator engine */

        bool genFuncVectors(unsigned int funcId);

        bool parseFuncBounds();

        bool allocateVMatrix();
        void releaseVMatrix();

        bool allocateVBounds();
        void releaseVBounds();
    };
} // proj1

#endif

// =========================
// End of proj1.h
// =========================
