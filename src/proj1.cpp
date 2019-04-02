/**
 * @file proj1.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the basic logic and functions to run
 * the cs471 project 1 experiment.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "proj1.h"
#include "datatable.h"
#include "datastats.h"
#include "stringutils.h"
#include <iostream>

using namespace std;
using namespace std::chrono;
using namespace proj1;

/**
 * @brief Construct a new mfuncExperiment object
 */
mfuncExperiment::mfuncExperiment() : vMatrix(nullptr), vBounds(nullptr), nbrDim(0), nbrSol(0)
{
}

/**
 * @brief Destroys the mfuncExperiment object
 * 
 */
mfuncExperiment::~mfuncExperiment()
{
    releaseVMatrix();
    releaseVBounds();
}

/**
 * @brief Initializes the CS471 project 1 experiment. Opens the given parameter file
 * and extracts test parameters. Allocates memory for function vectors and function
 * bounds. Extracts all function bounds.
 * 
 * @param paramFile File path to the parameter ini file
 * @return Returns true if initialization was successful. Otherwise false.
 */
bool mfuncExperiment::init(const char* paramFile)
{
    // Open and parse parameters file
    if (!iniParams.openFile(paramFile))
    {
        cout << "Experiment init failed: Unable to open param file: " << paramFile << endl;
        return false;
    }

    long numberSol;
    long numberDim;

    // Attempt to parse number of solutions and vector dimensions size
    // from iniParams
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

    // Get csv output file path
    resultsFile = iniParams.getEntry("test", "output_file");

    // Allocate memory for vector * solutions matrix
    if (!allocateVMatrix()) 
    {
        cout << "Experiment init failed: Unable to allocate vector matrix." << endl;
        return false;
    }

    // Allocate memory for function bounds
    if (!allocateVBounds()) 
    {
        cout << "Experiment init failed: Unable to allocate vector bounds array." << endl;
        return false;
    }

    // Fill function bounds array with data parsed from iniParams
    if (!parseFuncBounds())
    {
        cout << "Experiment init failed: Unable to parse vector bounds array." << endl;
        return false;
    }

    return true;
}

/**
 * @brief Executes all functions as specified in the CS471 project 1
 * document, records results, computes statistics, and outputs the
 * data as a *.csv file.
 * 
 * @return Returns 0 on success. Returns a non-zero error code on failure.
 */
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

    // Create a vector which is used to store all function results
    std::vector<double> fResults;
    double fTime = 0.0;

    // Execute all functions
    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        int err = runFunc(f, fResults, fTime);
        if (err)
            return err;
        else
        {
            // Insert function result and statistics into results table as a new row
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
        // Export results table to a *.csv file
        cout << "Exporting results to: " << resultsFile << endl;
        resultsTable.exportCSV(resultsFile.c_str());
    }

    return 0;
}

/**
 * @brief Runs the specifed function given by it's function id a certain number of times, 
 * records the execution time, and appends all results to the resultArrOut reference vector.
 * 
 * @param funcId The id of the function to run
 * @param resultArrOut Out reference variable that function results are appended to
 * @param timeOut Out reference variable that the execution time in ms is set to.
 * @return Returns 0 on success. Returns a non-zero error code on failure.
 */
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
    timeOut = duration_cast<nanoseconds>(t_end - t_start).count() / 1000000.0;

    return 0;
}

/**
 * @brief Helper function called by mfuncExperiment::runFunc.
 * Generates 'nbrSol' random vectors of dimensions 'nbrDim' that
 * are within the correct bounds for the specified function. All
 * vectors are stored in the vMatrix two-dimensional array.
 * 
 * @param funcId The id of the function to retrieve the vector bounds for
 * @return Returns true if all vectors were succesfully generated. Otherwise false.
 */
bool mfuncExperiment::genFuncVectors(unsigned int funcId)
{
    if (vMatrix == nullptr || vBounds == nullptr || funcId == 0 || funcId > mfunc::NUM_FUNCTIONS) return false;

    // Generate a new seed for the mersenne twister engine
    rgen = std::mt19937(rdev());

    // Set up a uniform distribution for the random number generator with the correct function bounds
    std::uniform_real_distribution<> dist(vBounds[funcId - 1].min, vBounds[funcId - 1].max);

    // Generate values for all vectors in vMatrix
    for (size_t s = 0; s < nbrSol; s++)
    {
        for (size_t d = 0; d < nbrDim; d++)
        {
            vMatrix[s][d] = dist(rgen);
        }
    }

    return true;
}

/**
 * @brief Helper function that extracts and parses all function bounds
 * from the ini parameters. Function bounds are stored in vBounds.
 * 
 * @return Returns true if all function bounds were succesfully parsed. Otherwise false.
 */
bool mfuncExperiment::parseFuncBounds()
{
    if (vBounds == nullptr) return false;

    const string delim = ",";
    const string section = "function_range";
    string s_min;
    string s_max;

    // Extract the bounds for each function
    for (unsigned int i = 1; i <= mfunc::NUM_FUNCTIONS; i++)
    {
        // Get bounds entry from ini file for current function
        string entry = iniParams.getEntry(section, to_string(i));
        if (entry.empty())
        {
            cout << "Error parsing bounds for function: " << i << endl;
            return false;
        }

        // Find index of ',' delimeter in entry string
        auto delimPos = entry.find(delim);
        if (delimPos == string::npos || delimPos >= entry.length() - 1)
        {
            cout << "Error parsing bounds for function: " << i << endl;
            return false;
        }

        // Split string and extract min/max strings
        s_min = entry.substr((size_t)0, delimPos);
        s_max = entry.substr(delimPos + 1, entry.length());
        util::s_trim(s_min);
        util::s_trim(s_max);

        // Attempt to parse min and max strings into double values
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

/**
 * @brief Allocates the vector * solutions matrix, which is used to store
 * all input vectors for the current function being tested.
 * 
 * @return Returns true if the memory was successfully allocated. Otherwise false.
 */
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

/**
 * @brief Releases memory used by the vector matrix
 */
void mfuncExperiment::releaseVMatrix()
{
    if (vMatrix == nullptr) return;

    for (size_t i = 0; i < nbrSol; i++)
    {
        if (vMatrix[i] != NULL)
        {
            delete[] vMatrix[i];
            vMatrix[i] = nullptr;
        }
    }

    delete[] vMatrix;
    vMatrix = nullptr;
}

/**
 * @brief Allocates function bounds array, which
 * stores the min and max input values for all function
 * vectors.
 * 
 * @return Returns true if the memory was successfully allocated. Otherwise false.
 */
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

/**
 * @brief Releases memory used by the function bounds array
 */
void mfuncExperiment::releaseVBounds()
{
    if (vBounds == nullptr) return;

    delete[] vBounds;
    vBounds = nullptr;
}

// =========================
// End of proj1.cpp
// =========================
