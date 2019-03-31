#include "proj1.h"
#include "datatable.h"
#include "datastats.h"
#include "stringutils.h"
#include <iostream>

using namespace std;
using namespace std::chrono;
using namespace proj1;

mfuncExperiment::mfuncExperiment() : vMatrix(nullptr), vBounds(nullptr), nbrDim(0), nbrSol(0)
{
}

mfuncExperiment::~mfuncExperiment()
{
    releaseVMatrix();
    releaseVBounds();
}

bool mfuncExperiment::init(const char* paramFile)
{
    if (!iniParams.openFile(paramFile))
    {
        cout << "Experiment init failed: Unable to open param file: " << paramFile << endl;
        return false;
    }

    long numberSol;
    long numberDim;

    try 
    {   
        std::string entry;
        
        entry = iniParams.getEntry("test", "number");
        if (entry.empty())
        {
            cout << "Experiment init failed: Param file missing [test]->number entry: " << paramFile << endl;
            return false;
        }

        numberSol = std::atol(entry.c_str());

        entry = iniParams.getEntry("test", "dimensions");
        if (entry.empty())
        {
            cout << "Experiment init failed: Param file missing [test]->dimensions entry: " << paramFile << endl;
            return false;
        }

        numberDim = std::atol(entry.c_str());

        if (numberSol <= 0)
        {
            cout << "Experiment init failed: Param file [test]->number entry out of bounds: " << paramFile << endl;
            return false;
        }

        if (numberDim <= 0)
        {
            cout << "Experiment init failed: Param file [test]->dimensions entry out of bounds: " << paramFile << endl;
            return false;
        }
    }
    catch (const std::exception& ex)
    {
        cout << "Experiment init failed: Exception while parsing param file: " << paramFile << endl;
        return false;
    }

    nbrSol = (size_t)numberSol;
    nbrDim = (size_t)numberDim;
    resultsFile = iniParams.getEntry("test", "output_file");

    if (!allocateVMatrix()) 
    {
        cout << "Experiment init failed: Unable to allocate vector matrix." << endl;
        return false;
    }

    if (!allocateVBounds()) 
    {
        cout << "Experiment init failed: Unable to allocate vector bounds array." << endl;
        return false;
    }

    if (!parseFuncBounds())
    {
        cout << "Experiment init failed: Unable to parse vector bounds array." << endl;
        return false;
    }

    return true;
}

int mfuncExperiment::runAllFunc()
{
    if (vMatrix == nullptr || nbrDim == 0 || nbrSol == 0) return 1;

    // function desc. | average | standard dev. | range | median | time
    mdata::DataTable resultsTable(8);
    resultsTable.setColLabel(0, "Function");
    resultsTable.setColLabel(1, "Vector Min");
    resultsTable.setColLabel(2, "Vector Max");
    resultsTable.setColLabel(3, "Average");
    resultsTable.setColLabel(4, "Standard Deviation");
    resultsTable.setColLabel(5, "Range");
    resultsTable.setColLabel(6, "Median");
    resultsTable.setColLabel(7, "Total Time (ms)");

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
            resultsTable.setEntry(rowIndex, 1, to_string(vBounds[f-1].min));
            resultsTable.setEntry(rowIndex, 2, to_string(vBounds[f-1].max));
            resultsTable.setEntry(rowIndex, 3, mdata::average(fResults));
            resultsTable.setEntry(rowIndex, 4, mdata::standardDeviation(fResults));
            resultsTable.setEntry(rowIndex, 5, mdata::range(fResults));
            resultsTable.setEntry(rowIndex, 6, mdata::median(fResults));
            resultsTable.setEntry(rowIndex, 7, fTime);
        }
    }

    if (!resultsFile.empty())
    {
        cout << "Exporting results to: " << resultsFile << endl;
        resultsTable.exportCSV(resultsFile.c_str());
    }

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
    if (vMatrix == nullptr || vBounds == nullptr || funcId == 0 || funcId > mfunc::NUM_FUNCTIONS) return false;

    // Generate new seed for mersenne twister engine
    rgen = std::mt19937(rdev());

    std::uniform_real_distribution<> dist(vBounds[funcId - 1].min, vBounds[funcId - 1].max);

    for (size_t s = 0; s < nbrSol; s++)
    {
        for (size_t d = 0; d < nbrDim; d++)
        {
            vMatrix[s][d] = dist(rgen);
        }
    }

    return true;
}

bool mfuncExperiment::parseFuncBounds()
{
    if (vBounds == nullptr) return false;

    const string delim = ",";
    const string section = "function_range";
    string s_min;
    string s_max;

    for (unsigned int i = 1; i <= mfunc::NUM_FUNCTIONS; i++)
    {
        string entry = iniParams.getEntry(section, to_string(i));
        if (entry.empty())
        {
            cout << "Error parsing bounds for function: " << i << endl;
            return false;
        }

        auto delimPos = entry.find(delim);
        if (delimPos == string::npos || delimPos >= entry.length() - 1)
        {
            cout << "Error parsing bounds for function: " << i << endl;
            return false;
        }

        s_min = entry.substr((size_t)0, delimPos);
        s_max = entry.substr(delimPos + 1, entry.length());
        util::s_trim(s_min);
        util::s_trim(s_max);

        try
        {
            RandomBounds& b = vBounds[i - 1];
            b.min = atof(s_min.c_str());
            b.max = atof(s_max.c_str());
        }
        catch(const std::exception& e)
        {
            cout << "Error parsing bounds for function: " << i << endl;
            std::cerr << e.what() << '\n';
            return false;
        }
    }

    return true;
}

bool mfuncExperiment::allocateVMatrix()
{
    if (nbrSol == 0) return false;

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

bool mfuncExperiment::allocateVBounds()
{
    if (nbrSol == 0) return false;

    try
    {
        releaseVBounds();

        vBounds = new RandomBounds[nbrSol];
        
        return true;
    }
    catch(const std::exception& e)
    {
        return false;
    }
}

void mfuncExperiment::releaseVBounds()
{
    if (vBounds == nullptr) return;

    delete vBounds;
    vBounds = nullptr;
}