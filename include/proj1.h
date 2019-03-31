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
    struct RandomBounds
    {
        double min = 0.0;
        double max = 0.0;
    };

    class mfuncExperiment
    {
    public:
        mfuncExperiment();
        ~mfuncExperiment();
        bool init(const char* paramFile);
        int runAllFunc();
        int runFunc(unsigned int funcId, std::vector<double>& resultArrOut, double& timeOut);
    private:
        util::IniReader iniParams;
        std::string resultsFile;
        size_t nbrDim;
        size_t nbrSol;
        double** vMatrix;
        RandomBounds* vBounds;

        std::random_device rdev;
        std::mt19937 rgen;

        bool genFuncVectors(unsigned int funcId);

        bool parseFuncBounds();

        bool allocateVMatrix();
        void releaseVMatrix();

        bool allocateVBounds();
        void releaseVBounds();
    };
} // proj1


#endif