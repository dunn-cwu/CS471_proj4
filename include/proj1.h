#ifndef __PROJ1_H
#define __PROJ1_H

#include <random>
#include <chrono>
#include <vector>
#include "mfunc.h"

namespace proj1
{
    struct RandomBounds
    {
        double min = 0.0;
        double max = 0.0;
    };

    class mfuncExperiment
    {
    public:
        mfuncExperiment(size_t dimensions, size_t solutions);
        ~mfuncExperiment();
        int runAllFunc(const char* dataResultsFile);
        int runFunc(unsigned int funcId, std::vector<double>& resultArrOut, double& timeOut);
    private:
        const size_t nbrDim;
        const size_t nbrSol;
        double** vMatrix;
        RandomBounds* vBounds;

        std::random_device rdev;
        std::mt19937 rgen;

        bool genFuncVectors(unsigned int funcId);

        bool allocateVMatrix();
        void releaseVMatrix();
    };
} // proj1


#endif