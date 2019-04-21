/**
 * @file cs471.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the Experiment class.
 * Contains the basic logic and functions to run the cs471 project experiment.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "experiment.h"
#include "datatable.h"
#include "blindsearch.h"
#include "localsearch.h"
#include "stringutils.h"
#include "mem.h"
#include <iostream>
#include <fstream>

#define INI_TEST_SECTION "test"
#define INI_FUNC_RANGE_SECTION "function_range"
#define INI_TEST_POPULATION "population"
#define INI_TEST_DIMENSIONS "dimensions"
#define INI_TEST_ITERATIONS "iterations"
#define INI_TEST_NUMTHREADS "num_threads"
#define INI_TEST_ALPHA "alpha"
#define INI_TEST_ALGORITHM "algorithm"
#define INI_TEST_RESULTSFILE "results_file"
#define INI_TEST_EXECTIMESFILE "exec_times_file"

using namespace std;
using namespace std::chrono;
using namespace mfunc;

/**
 * @brief Construct a new Experiment object
 */
template<class T>
Experiment<T>::Experiment() 
    : vBounds(nullptr), tPool(nullptr), resultsFile(""), execTimesFile(""), iterations(0)
{
}

/**
 * @brief Destroys the Experiment object
 * 
 */
template<class T>
Experiment<T>::~Experiment()
{
    releaseThreadPool();
    releasePopulationPool();
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
template<class T>
bool Experiment<T>::init(const char* paramFile)
{
    try 
    {
        // Open and parse parameters file
        if (!iniParams.openFile(paramFile))
        {
            cerr << "Experiment init failed: Unable to open param file: " << paramFile << endl;
            return false;
        }

        long numberSol = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_POPULATION);
        long numberDim = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_DIMENSIONS);
        long numberIter = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_ITERATIONS);
        long numberThreads = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_NUMTHREADS);
        alpha = iniParams.getEntryAs<T>(INI_TEST_SECTION, INI_TEST_ALPHA);
        unsigned int selectedAlg = iniParams.getEntryAs<unsigned int>(INI_TEST_SECTION, INI_TEST_ALGORITHM);
        resultsFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_RESULTSFILE);
        execTimesFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_EXECTIMESFILE);

        if (numberSol <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_POPULATION << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }
        else if (numberDim <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_DIMENSIONS << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }
        else if (numberIter <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_ITERATIONS << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }
        else if (numberThreads <= 0)
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_NUMTHREADS << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }
        else if (alpha == 0)
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_ALPHA << " is missing or is equal to zero: " << paramFile << endl;
            return false;
        }
        else if (selectedAlg >= static_cast<unsigned int>(enums::Algorithm::Count))
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_ALGORITHM << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }

        iterations = (size_t)numberIter;
        testAlg = static_cast<enums::Algorithm>(selectedAlg);

        cout << "Population size: " << numberSol << endl;
        cout << "Dimensions: " << numberDim << endl;
        cout << "Iterations: " << iterations << endl;
        cout << "Alpha value: " << alpha << endl;
        cout << "Algorithm: " << enums::AlgorithmNames::get(testAlg) << endl;

        // Allocate memory for vector * solutions matrix
        if (!allocatePopulationPool((size_t)numberThreads, (size_t)numberSol, (size_t)numberDim))
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
    catch (const std::exception& ex)
    {
        cerr << "Exception occurred while initializing experiment: " << ex.what() << endl;
        return false;
    }
    catch (...)
    {
        cerr << "Unknown Exception occurred while initializing experiment." << endl;
        return false;
    }
}

/**
 * @brief Executes all functions as specified in the CS471 project 1
 * document, records results, computes statistics, and outputs the
 * data as a *.csv file.
 * 
 * @return Returns 0 on success. Returns a non-zero error code on failure.
 */
template<class T>
int Experiment<T>::testAllFunc()
{
    if (populationsPool.size() == 0) return 1;

    mdata::DataTable<T> resultsTable(iterations, (size_t)NUM_FUNCTIONS);
    mdata::DataTable<T> execTimesTable(iterations, (size_t)NUM_FUNCTIONS);

    mdata::TestParameters<T>* tParams = 
        new mdata::TestParameters<T>[NUM_FUNCTIONS];

    for (unsigned int i = 0; i < NUM_FUNCTIONS; i++)
    {
        mdata::TestParameters<T>& curParam = tParams[i];
        curParam.funcId = i + 1;
        curParam.alg = testAlg;
        curParam.iterations = iterations;
        curParam.resultsTable = &resultsTable;
        curParam.execTimesTable = &execTimesTable;
        curParam.resultsCol = i;
        curParam.execTimesCol = i;

        resultsTable.setColLabel((size_t)i, FunctionDesc::get(curParam.funcId));
        execTimesTable.setColLabel((size_t)i, FunctionDesc::get(curParam.funcId));
    }

    std::vector<std::future<int>> testFutures;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    // Queue all functions to be ran by our thread pool
    for (unsigned int f = 1; f <= NUM_FUNCTIONS; f++)
    {
        testFutures.emplace_back(
            tPool->enqueue(&Experiment::testFunc, this, &tParams[f - 1])
        );
    }

    // Join and get all function return values from thread pool workers
    for (unsigned int f = 1; f <= NUM_FUNCTIONS; f++)
    {
        int errCode = testFutures[f-1].get();
        if (errCode)
        {
            cerr << "Error occurred while testing function " << f << endl;
            return errCode;
        }
        
    }

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    double totalExecTime = (double)duration_cast<nanoseconds>(t_end - t_start).count() / 1000000000.0;

    cout << "Test finished. Total time: " << totalExecTime << " seconds." << endl;

    if (!resultsFile.empty())
    {
        // Export results table to a *.csv file
        cout << "Exporting results to: " << resultsFile << endl;
        resultsTable.exportCSV(resultsFile.c_str());
    }

    if (!execTimesFile.empty())
    {
        // Export exec times table to a *.csv file
        cout << "Exporting execution times to: " << execTimesFile << endl;
        execTimesTable.exportCSV(execTimesFile.c_str());
    }

    return 0;
}
 
template<class T>
int Experiment<T>::testFunc(mdata::TestParameters<T>* tParams)
{
    mdata::Population<T>* pop = popPoolRemove();
    mdata::SearchAlgorithm<T>* alg;

    switch (tParams->alg)
    {
        case enums::Algorithm::BlindSearch:
            alg = new mdata::BlindSearch<T>();
            break;
        case enums::Algorithm::LocalSearch:
            alg = new mdata::LocalSearch<T>();
            break;
        default:
            cerr << "Invalid algorithm selected." << endl;
            return 1;
    }

    RandomBounds<T>& funcBounds = vBounds[tParams->funcId - 1];
    int returnVal = 0;

    for (size_t i = 0; i < tParams->iterations; i++)
    {
        auto result = alg->run(Functions<T>::get(tParams->funcId), funcBounds.min, funcBounds.max, pop, alpha);
        if (result.err) 
        {
            returnVal = result.err;
            break;
        }

        tParams->resultsTable->setEntry(i, tParams->resultsCol, result.fitness);
        tParams->execTimesTable->setEntry(i, tParams->execTimesCol, result.execTime);
    }

    cout << "F" << tParams->funcId << " done." << endl << flush;

    delete alg;
    popPoolAdd(pop);

    return returnVal;
}

template<class T>
mdata::Population<T>* Experiment<T>::popPoolRemove()
{
    mdata::Population<T>* retPop = nullptr;
    std::chrono::microseconds waitTime(10);

    while (true)
    {
        {
            std::lock_guard<std::mutex> lk(popPoolMutex);
            if (populationsPool.size() > 0)
            {
                retPop = populationsPool.back();
                populationsPool.pop_back();
            }
        }

        if (retPop != nullptr)
            return retPop;
        else
            std::this_thread::sleep_for(waitTime);
    }
}

template<class T>
void Experiment<T>::popPoolAdd(mdata::Population<T>* popPtr)
{
    if (popPtr == nullptr) return;

    std::lock_guard<std::mutex> lk(popPoolMutex);

    populationsPool.push_back(popPtr);
}

/**
 * @brief Helper function that extracts and parses all function bounds
 * from the ini parameters. Function bounds are stored in vBounds.
 * 
 * @return Returns true if all function bounds were succesfully parsed. Otherwise false.
 */
template<class T>
bool Experiment<T>::parseFuncBounds()
{
    if (vBounds == nullptr) return false;

    const string delim = ",";
    const string section = "function_range";
    string s_min;
    string s_max;

    // Extract the bounds for each function
    for (unsigned int i = 1; i <= NUM_FUNCTIONS; i++)
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
            RandomBounds<T>& b = vBounds[i - 1];
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
 * @brief Allocates the vector * solutions matrix, which is used to store
 * all input vectors for the current function being tested.
 * 
 * @return Returns true if the memory was successfully allocated. Otherwise false.
 */
template<class T>
bool Experiment<T>::allocatePopulationPool(size_t count, size_t popSize, size_t dimensions)
{
    releasePopulationPool();

    std::lock_guard<std::mutex> lk(popPoolMutex);

    try
    {
        for (int i = 0; i < count; i++)
        {
            auto newPop = new(std::nothrow) mdata::Population<T>(popSize, dimensions);
            if (newPop == nullptr)
            {
                std::cerr << "Error allocating populations." << '\n';
                return false;
            }

            populationsPool.push_back(newPop);
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
template<class T>
void Experiment<T>::releasePopulationPool()
{
    std::lock_guard<std::mutex> lk(popPoolMutex);

    if (populationsPool.size() == 0) return;

    for (int i = 0; i < populationsPool.size(); i++)
    {
        if (populationsPool[i] != nullptr)
        {
            delete populationsPool[i];
            populationsPool[i] = nullptr;
        }
    }

    populationsPool.clear();
}

/**
 * @brief Allocates function bounds array, which
 * stores the min and max input values for all function
 * vectors.
 * 
 * @return Returns true if the memory was successfully allocated. Otherwise false.
 */
template<class T>
bool Experiment<T>::allocateVBounds()
{
    vBounds = util::allocArray<RandomBounds<T>>(NUM_FUNCTIONS);
    return vBounds != nullptr;
}

/**
 * @brief Releases memory used by the function bounds array
 */
template<class T>
void Experiment<T>::releaseVBounds()
{
    if (vBounds == nullptr) return;

    util::releaseArray<RandomBounds<T>>(vBounds);
}

template<class T>
bool Experiment<T>::allocateThreadPool(size_t numThreads)
{
    releaseThreadPool();

    tPool = new(std::nothrow) ThreadPool(numThreads);
    return tPool != nullptr;
}

template<class T>
void Experiment<T>::releaseThreadPool()
{
    if (tPool == nullptr) return;

    delete tPool;
    tPool = nullptr;
}

template class mfunc::Experiment<float>;
template class mfunc::Experiment<double>;
template class mfunc::Experiment<long double>;

// =========================
// End of proj1.cpp
// =========================
