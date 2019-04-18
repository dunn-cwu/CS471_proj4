/**
 * @file cs471.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the mfuncExperiment class.
 * Contains the basic logic and functions to run the cs471 project experiment.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "cs471.h"
#include "datatable.h"
#include "datastats.h"
#include "blindsearch.h"
#include "stringutils.h"
#include "mem.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace std::chrono;
using namespace cs471;

/**
 * @brief Construct a new mfuncExperiment object
 */
mfuncExperiment::mfuncExperiment() : populations(nullptr), vBounds(nullptr), outputPop(false), outputFitness(false), tPool(nullptr)
{
}

/**
 * @brief Destroys the mfuncExperiment object
 * 
 */
mfuncExperiment::~mfuncExperiment()
{
    releaseThreadPool();
    releasePopulations();
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
        cerr << "Experiment init failed: Unable to open param file: " << paramFile << endl;
        return false;
    }

    long numberSol;
    long numberDim;
    long numberThreads;

    // Attempt to parse number of solutions, vector dimensions size, and number of threads
    // from iniParams
    try 
    {   
        std::string entry;
        
        entry = iniParams.getEntry("test", "population");
        if (entry.empty())
        {
            cerr << "Experiment init failed: Param file missing [test]->population entry: " << paramFile << endl;
            return false;
        }

        numberSol = std::atol(entry.c_str());

        entry = iniParams.getEntry("test", "dimensions");
        if (entry.empty())
        {
            cerr << "Experiment init failed: Param file missing [test]->dimensions entry: " << paramFile << endl;
            return false;
        }

        numberDim = std::atol(entry.c_str());

        entry = iniParams.getEntry("test", "num_threads");
        if (entry.empty())
        {
            cerr << "Experiment init failed: Param file missing [test]->num_threads entry: " << paramFile << endl;
            return false;
        }

        numberThreads = std::atol(entry.c_str());

        if (numberSol <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->population entry out of bounds: " << paramFile << endl;
            return false;
        }
        else if (numberDim <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->dimensions entry out of bounds: " << paramFile << endl;
            return false;
        }
        else if (numberThreads <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->num_threads entry out of bounds: " << paramFile << endl;
            return false;
        }
    }
    catch (const std::exception& ex)
    {
        cerr << "Experiment init failed: Exception while parsing param file: " << paramFile << endl;
        return false;
    }

    cout << "Population size: " << numberSol << endl;
    cout << "Dimensions: " << numberDim << endl;

    // Get csv output file path
    resultsFile = iniParams.getEntry("test", "results_file");
    outputPop = iniParams.getEntry("test", "output_population") == "true";
    outputFitness = iniParams.getEntry("test", "output_fitness") == "true";

    // Allocate memory for vector * solutions matrix
    if (!allocatePopulations((size_t)numberSol, (size_t)numberDim))
    {
        cerr << "Experiment init failed: Unable to allocate populations." << endl;
        return false;
    }

    // Allocate memory for function bounds
    if (!allocateVBounds()) 
    {
        cerr << "Experiment init failed: Unable to allocate vector bounds array." << endl;
        return false;
    }

    // Allocate thread pool
    if (!allocateThreadPool((size_t)numberThreads))
    {
        cerr << "Experiment init failed: Unable to allocate thread pool." << endl;
        return false;
    }

    cout << "Started " << numberThreads << " worker threads." << endl;

    // Fill function bounds array with data parsed from iniParams
    if (!parseFuncBounds())
    {
        cerr << "Experiment init failed: Unable to parse vector bounds array." << endl;
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
    if (populations == nullptr) return 1;

    std::vector<std::future<mdata::TestResult>> testFutures;

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

    ofstream fitnessFile;
    if (outputFitness)
    {
        std::string fitFile = "fitness-dim_";
        fitFile += to_string(populations[0]->getDimensionsSize());
        fitFile += ".csv";
        fitnessFile.open(fitFile, ios::out | ios::trunc);
        if (!fitnessFile.good()) outputFitness = false;
    }

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    // Queue all functions to be ran by our thread pool
    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        testFutures.emplace_back(
            tPool->enqueue(&mfuncExperiment::runFunc, this, f)
        );
    }

    // Join and get all function return values from thread pool workers
    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        mdata::TestResult tRes = testFutures[f-1].get();
        if (tRes.err)
        {
            if (outputFitness) fitnessFile.close();
            return tRes.err;
        }
        else
        {
            mdata::Population<double>* curPopObj = populations[f - 1];

            // Export all population data if flag is set
            if (outputPop)
                exportPop(f);

            // Export all fitness data if flag is set
            if (outputFitness)
            {
                fitnessFile << mfunc::fDesc(f) << ",";
                curPopObj->outputFitness(fitnessFile, ",", "\n");
            }

            // Insert function results and statistics into the results table as a new row
            unsigned int rowIndex = resultsTable.addRow();
            resultsTable.setEntry(rowIndex, 0, mfunc::fDesc(f));
            resultsTable.setEntry(rowIndex, 1, to_string(vBounds[f-1].min));
            resultsTable.setEntry(rowIndex, 2, to_string(vBounds[f-1].max));
            resultsTable.setEntry(rowIndex, 3, curPopObj->getFitnessAverage());
            resultsTable.setEntry(rowIndex, 4, curPopObj->getFitnessStandardDev());
            resultsTable.setEntry(rowIndex, 5, curPopObj->getFitnessRange());
            resultsTable.setEntry(rowIndex, 6, curPopObj->getFitnessMedian());
            resultsTable.setEntry(rowIndex, 7, tRes.execTime);
        }
    }

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    double totalExecTime = (double)duration_cast<nanoseconds>(t_end - t_start).count() / 1000000.0;

    cout << "Total test time: " << totalExecTime << " miliseconds." << endl;

    if (!resultsFile.empty())
    {
        // Export results table to a *.csv file
        cout << "Exporting results to: " << resultsFile << endl;
        resultsTable.exportCSV(resultsFile.c_str());
    }

    if (outputFitness) fitnessFile.close();
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
mdata::TestResult mfuncExperiment::runFunc(unsigned int funcId)
{
    if (!genFuncVectors(funcId)) return mdata::TestResult(1, 0.0);

    mfunc::mfuncPtr fPtr = mfunc::fGet(funcId);
    if (fPtr == nullptr) return mdata::TestResult(2, 0.0);

    auto curPopObj = populations[funcId - 1];
    size_t nbrSol = curPopObj->getPopulationSize();

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    for (int i = 0; i < nbrSol; i++)
    {
        if (!curPopObj->setFitness(i, fPtr))
            return mdata::TestResult(4, 0.0);
    }
    
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    double execTime = (double)duration_cast<nanoseconds>(t_end - t_start).count() / 1000000.0;

    return mdata::TestResult(0, execTime);
}

/**
 * @brief Helper function called by mfuncExperiment::runFunc.
 * Generates a random population to be used as input for the functions.
 * 
 * @param funcId The id of the function to retrieve the vector bounds for
 * @return Returns true if the population was succesfully generated. Otherwise false.
 */
bool mfuncExperiment::genFuncVectors(unsigned int funcId)
{
    if (populations == nullptr || vBounds == nullptr || funcId == 0 || funcId > mfunc::NUM_FUNCTIONS) return false;

    return populations[funcId - 1]->generate(vBounds[funcId - 1].min, vBounds[funcId - 1].max);
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
            cerr << "Error parsing bounds for function: " << i << endl;
            return false;
        }

        // Find index of ',' delimeter in entry string
        auto delimPos = entry.find(delim);
        if (delimPos == string::npos || delimPos >= entry.length() - 1)
        {
            cerr << "Error parsing bounds for function: " << i << endl;
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
            RandomBounds<double>& b = vBounds[i - 1];
            b.min = atof(s_min.c_str());
            b.max = atof(s_max.c_str());
        }
        catch(const std::exception& e)
        {
            cerr << "Error parsing bounds for function: " << i << endl;
            std::cerr << e.what() << '\n';
            return false;
        }
    }

    return true;
}

/**
 * @brief Helper function that exports all current population data
 * to a file.
 * 
 * @param func The id of the function data being exported
 */
void mfuncExperiment::exportPop(unsigned int func)
{
    ofstream popFile;

    std::string fName = "pop-func_";
    fName += std::to_string(func);
    fName += "-dim_";
    fName += std::to_string(populations[func - 1]->getDimensionsSize());
    fName += ".csv";

    popFile.open(fName.c_str(), ios::out | ios::trunc);
    if (!popFile.good())
    {
        cerr << "Unable to open pop output file." << endl;
        return;
    }

    populations[func - 1]->outputPopulation(popFile, ",", "\n");
    popFile.close();
}

/**
 * @brief Allocates the vector * solutions matrix, which is used to store
 * all input vectors for the current function being tested.
 * 
 * @return Returns true if the memory was successfully allocated. Otherwise false.
 */
bool mfuncExperiment::allocatePopulations(size_t popSize, size_t dimensions)
{
    releasePopulations();

    try
    {
        populations = new(std::nothrow) mdata::Population<double>*[mfunc::NUM_FUNCTIONS];
        if (populations == nullptr) return false;

        for (int i = 0; i < mfunc::NUM_FUNCTIONS; i++)
        {
            populations[i] = new(std::nothrow) mdata::Population<double>(popSize, dimensions);
            if (populations[i] == nullptr)
            {
                releasePopulations();
                return false;
            }
        }

        return true;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return false;
    }
}

/**
 * @brief Releases memory used by the vector matrix
 */
void mfuncExperiment::releasePopulations()
{
    if (populations == nullptr) return;

    for (int i = 0; i < mfunc::NUM_FUNCTIONS; i++)
    {
        if (populations[i] != nullptr)
        {
            delete populations[i];
            populations[i] = nullptr;
        }
    }

    delete[] populations;
    populations = nullptr;
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
    if (populations == nullptr) return false;

    vBounds = util::allocArray<RandomBounds<double>>(mfunc::NUM_FUNCTIONS);
    return vBounds != nullptr;
}

/**
 * @brief Releases memory used by the function bounds array
 */
void mfuncExperiment::releaseVBounds()
{
    if (vBounds == nullptr) return;

    util::releaseArray<RandomBounds<double>>(vBounds);
}

bool mfuncExperiment::allocateThreadPool(size_t numThreads)
{
    releaseThreadPool();

    tPool = new(std::nothrow) ThreadPool(numThreads);
    return tPool != nullptr;
}

void mfuncExperiment::releaseThreadPool()
{
    if (tPool == nullptr) return;

    delete tPool;
    tPool = nullptr;
}

// =========================
// End of proj1.cpp
// =========================
