#include "proj1.h"
#include "datatable.h"
#include "datastats.h"
#include <iostream>

using namespace std;
using namespace std::chrono;
using namespace proj1;

mfuncExperiment::mfuncExperiment(size_t dimensions, size_t solutions) : vMatrix(nullptr), vBounds(nullptr), nbrDim(dimensions), nbrSol(solutions)
{
    if (!allocateVMatrix())
        throw "Error allocating vector matrix";
}

mfuncExperiment::~mfuncExperiment()
{
    releaseVMatrix();
}

int mfuncExperiment::runAllFunc(const char* dataResultsFile)
{
    // function desc. | average | standard dev. | range | median | time
    mdata::DataTable resultsTable(6);
    resultsTable.setColLabel(0, "Function");
    resultsTable.setColLabel(1, "Average");
    resultsTable.setColLabel(2, "Standard Deviation");
    resultsTable.setColLabel(3, "Range");
    resultsTable.setColLabel(4, "Median");
    resultsTable.setColLabel(5, "Total Time (ms)");

    std::vector<double> fResults;
    double fTime = 0.0;

    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        int err = runFunc(f, fResults, fTime);
        if (err)
            return err;
        else
        {
            unsigned int rowIndex = resultsTable.addRow();
            resultsTable.setEntry(rowIndex, 0, mfunc::fDesc(f));
            resultsTable.setEntry(rowIndex, 1, mdata::average(fResults));
            resultsTable.setEntry(rowIndex, 2, mdata::standardDeviation(fResults));
            resultsTable.setEntry(rowIndex, 3, mdata::range(fResults));
            resultsTable.setEntry(rowIndex, 4, mdata::median(fResults));
            resultsTable.setEntry(rowIndex, 5, fTime);
        }
    }

    resultsTable.exportCSV(dataResultsFile);

    return 0;
}

int mfuncExperiment::runFunc(unsigned int funcId, std::vector<double>& resultArrOut, double& timeOut)
{
    if (!genFuncVectors(funcId)) return 1;

    resultArrOut.clear();
    resultArrOut.reserve(nbrSol);

    double fResult = 0;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    for (int i = 0; i < nbrSol; i++)
    {
        if (!mfunc::fExec(funcId, vMatrix[i], nbrDim, fResult))
            return 2;

        resultArrOut.push_back(fResult);
    }
    
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    timeOut = duration_cast<microseconds>(t_end - t_start).count() / 1000.0;

    return 0;
}

bool mfuncExperiment::genFuncVectors(unsigned int funcId)
{
    if (vMatrix == nullptr) return false;

    // Generate new seed for mersenne twister engine
    rgen = std::mt19937(rdev());

    std::uniform_real_distribution<> dist(0, M_PI);

    for (size_t s = 0; s < nbrSol; s++)
    {
        for (size_t d = 0; d < nbrDim; d++)
        {
            vMatrix[s][d] = dist(rgen);
        }
    }

    return true;
}

bool mfuncExperiment::allocateVMatrix()
{
    try
    {
        releaseVMatrix();
        vMatrix = new double*[nbrSol];
        
        for (size_t i = 0; i < nbrSol; i++)
        {
            vMatrix[i] = new double[nbrDim];
        }

        return true;
    }
    catch(const std::exception& e)
    {
        return false;
    }
}

void mfuncExperiment::releaseVMatrix()
{
    if (vMatrix == nullptr) return;

    for (size_t i = 0; i < nbrSol; i++)
    {
        if (vMatrix[i] != NULL)
        {
            delete vMatrix[i];
            vMatrix[i] = nullptr;
        }
    }

    delete vMatrix;
    vMatrix = nullptr;
}