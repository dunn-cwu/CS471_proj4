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

#include "cs471.h"
#include "datatable.h"
#include "datastats.h"
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
mfuncExperiment::mfuncExperiment() : population(nullptr), vBounds(nullptr), outputPop(false), outputFitness(false)
{
}

/**
 * @brief Destroys the mfuncExperiment object
 * 
 */
mfuncExperiment::~mfuncExperiment()
{
    releasePopulation();
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

    // Attempt to parse number of solutions and vector dimensions size
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

        if (numberSol <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->population entry out of bounds: " << paramFile << endl;
            return false;
        }

        if (numberDim <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->dimensions entry out of bounds: " << paramFile << endl;
            return false;
        }
    }
    catch (const std::exception& ex)
    {
        cerr << "Experiment init failed: Exception while parsing param file: " << paramFile << endl;
        return false;
    }

    // Get csv output file path
    resultsFile = iniParams.getEntry("test", "results_file");
    outputPop = iniParams.getEntry("test", "output_population") == "true";
    outputFitness = iniParams.getEntry("test", "output_fitness") == "true";

    // Allocate memory for vector * solutions matrix
    if (!allocatePopulation((size_t)numberSol, (size_t)numberDim))
    {
        cerr << "Experiment init failed: Unable to allocate population matrix." << endl;
        return false;
    }

    // Allocate memory for function bounds
    if (!allocateVBounds()) 
    {
        cerr << "Experiment init failed: Unable to allocate vector bounds array." << endl;
        return false;
    }

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
    if (population == nullptr || !population->isReady()) return 1;

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

    double fTime = 0.0;
    ofstream fitnessFile;
    if (outputFitness)
    {
        std::string fitFile = "fitness-dim_";
        fitFile += to_string(population->getDimensionsSize());
        fitFile += ".csv";
        fitnessFile.open(fitFile, ios::out | ios::trunc);
        if (!fitnessFile.good()) outputFitness = false;
    }

    // Execute all functions
    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        int err = runFunc(f, fTime);
        if (err)
        {
            if (outputFitness) fitnessFile.close();
            return err;
        }
        else
        {
            // Export all population data if flag is set
            if (outputPop)
                exportPop(f);

            // Export all fitness data if flag is set
            if (outputFitness)
            {
                fitnessFile << mfunc::fDesc(f) << ",";
                population->outputFitness(fitnessFile, ",", "\n");
            }

            // Insert function results and statistics into the results table as a new row
            unsigned int rowIndex = resultsTable.addRow();
            resultsTable.setEntry(rowIndex, 0, mfunc::fDesc(f));
            resultsTable.setEntry(rowIndex, 1, to_string(vBounds[f-1].min));
            resultsTable.setEntry(rowIndex, 2, to_string(vBounds[f-1].max));
            resultsTable.setEntry(rowIndex, 3, population->getFitnessAverage());
            resultsTable.setEntry(rowIndex, 4, population->getFitnessStandardDev());
            resultsTable.setEntry(rowIndex, 5, population->getFitnessRange());
            resultsTable.setEntry(rowIndex, 6, population->getFitnessMedian());
            resultsTable.setEntry(rowIndex, 7, fTime);
        }
    }

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
int mfuncExperiment::runFunc(unsigned int funcId, double& timeOut)
{
    if (!genFuncVectors(funcId)) return 1;

    double fResult = 0;
    size_t nbrSol = population->getPopulationSize();
    size_t nbrDim = population->getDimensionsSize();
    double* curPop = nullptr;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    for (int i = 0; i < nbrSol; i++)
    {
        curPop = population->getPopulation(i);
        if (curPop == nullptr || !mfunc::fExec(funcId, curPop, nbrDim, fResult))
            return 2;

        if (!population->setFitness(i, fResult))
            return 3;
    }
    
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    timeOut = (double)duration_cast<nanoseconds>(t_end - t_start).count() / 1000000.0;

    return 0;
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
    if (population == nullptr || vBounds == nullptr || funcId == 0 || funcId > mfunc::NUM_FUNCTIONS) return false;

    return population->generate(vBounds[funcId - 1].min, vBounds[funcId - 1].max);
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
    fName += std::to_string(population->getDimensionsSize());
    fName += ".csv";

    popFile.open(fName.c_str(), ios::out | ios::trunc);
    if (!popFile.good())
    {
        cerr << "Unable to open pop output file." << endl;
        return;
    }

    population->outputPopulation(popFile, ",", "\n");
    popFile.close();
}

/**
 * @brief Allocates the vector * solutions matrix, which is used to store
 * all input vectors for the current function being tested.
 * 
 * @return Returns true if the memory was successfully allocated. Otherwise false.
 */
bool mfuncExperiment::allocatePopulation(size_t popSize, size_t dimensions)
{
    releasePopulation();

    try
    {
        population = new(std::nothrow) mdata::Population<double>(popSize, dimensions);
        return population != nullptr;
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
void mfuncExperiment::releasePopulation()
{
    if (population == nullptr) return;

    delete population;
    population = nullptr;
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
    if (population == nullptr) return false;

    vBounds = util::allocArray<RandomBounds<double>>(population->getPopulationSize());
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

// =========================
// End of proj1.cpp
// =========================
