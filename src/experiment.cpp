/**
 * @file experiment.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the Experiment class.
 * Contains the basic logic and functions to run the cs471 project experiment.
 * @version 0.2
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include "experiment.h"
#include "datatable.h"
#include "stringutils.h"
#include "mem.h"
#include "geneticalg.h"

// Ini file string sections and keys
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
 * @brief Initializes the CS471 project 2 experiment. Opens the given parameter file
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

        // Extract test parameters from ini file
        long numberSol = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_POPULATION);
        long numberDim = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_DIMENSIONS);
        long numberIter = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_ITERATIONS);
        long numberThreads = iniParams.getEntryAs<long>(INI_TEST_SECTION, INI_TEST_NUMTHREADS);
        alpha = iniParams.getEntryAs<T>(INI_TEST_SECTION, INI_TEST_ALPHA);
        unsigned int selectedAlg = iniParams.getEntryAs<unsigned int>(INI_TEST_SECTION, INI_TEST_ALGORITHM);
        resultsFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_RESULTSFILE);
        execTimesFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_EXECTIMESFILE);

        // Verify test parameters
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
        // else if (selectedAlg >= static_cast<unsigned int>(enums::Algorithm::Count))
        // {
        //     cerr << "Experiment init failed: Param file [test]->" 
        //         << INI_TEST_ALGORITHM << " entry missing or out of bounds: " << paramFile << endl;
        //     return false;
        // }

        // Cast iterations and test algorithm to correct types
        iterations = (size_t)numberIter;
        // testAlg = static_cast<enums::Algorithm>(selectedAlg);

        // Print test parameters to console
        cout << "Population size: " << numberSol << endl;
        cout << "Dimensions: " << numberDim << endl;
        cout << "Iterations: " << iterations << endl;
        cout << "Alpha value: " << alpha << endl;
        // cout << "Algorithm: " << enums::AlgorithmNames::get(testAlg) << endl;

        // Allocate memory for all population objects. We need one for each thread to prevent conflicts.
        if (!allocatePopulationPool((size_t)numberThreads, (size_t)numberSol, (size_t)numberDim))
        {
            cerr << "Experiment init failed: Unable to allocate populations." << endl;
            return false;
        }

        // Allocate memory for function vector bounds
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

        // Allocate thread pool
        if (!allocateThreadPool((size_t)numberThreads))
        {
            cerr << "Experiment init failed: Unable to allocate thread pool." << endl;
            return false;
        }

        cout << "Started " << numberThreads << " worker threads ..." << endl;

        // Ready to run an experiment
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
 * @brief Executes all functions as specified in the CS471 project 2
 * document, records results, and outputs the data as a *.csv file.
 * 
 * @return Returns 0 on success. Returns a non-zero error code on failure.
 */
template<class T>
int Experiment<T>::testAllFunc()
{
    if (populationsPool.size() == 0) return 1;

    mdata::DataTable<T> testTable(100, 50);
    mdata::Population<T> pop1(200, 30);
    mdata::Population<T> pop2(200, 30);

    for (size_t exp = 0; exp < 50; exp++)
    {
        GAParams<T> gParams;
        gParams.fitnessTable = &testTable;
        gParams.fitTableCol = exp;
        gParams.mainPop = &pop1;
        gParams.auxPop = &pop2;
        gParams.fPtr = Functions<T>::get(1);
        gParams.fMinBound = vBounds[0].min;
        gParams.fMaxBound = vBounds[0].max;
        gParams.generations = 100;
        gParams.crProb = 0.8;
        gParams.mutProb = 0.005;
        gParams.mutRange = 0.1;
        gParams.mutPrec = 1;
        gParams.elitismRate = 0.2;

        mfunc::GeneticAlgorithm<T>::run(gParams);
    }

    testTable.exportCSV("genalg_test.csv");


    /* testP.setFitnessNormalization(true);
    testP.generate(vBounds[0].min, vBounds[0].max);
    testP.calcAllFitness(Functions<T>::get(1));

    testP.debugOutputAll();

    cout << endl << "=======================" << endl;

    testP.sortDescendByFitness();
    testP.debugOutputAll();
 */
/*     // Construct results and execution times tables
    mdata::DataTable<T> resultsTable(iterations, (size_t)NUM_FUNCTIONS);
    mdata::DataTable<T> execTimesTable(iterations, (size_t)NUM_FUNCTIONS);

    // Prepare thread futures vector, used to ensure all async tasks complete
    // succesfully.
    std::vector<std::future<int>> testFutures;

    // Start recording total execution time
    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    // For each of the NUM_FUNCTIONS functions, prepare a TestParameters
    // struct and queue an asynchronous test that will be picked up and
    // executed by one of the threads in the thread pool.
    for (unsigned int i = 0; i < NUM_FUNCTIONS; i++)
    {
        // Update column labels for results and exec times tables
        resultsTable.setColLabel((size_t)i, FunctionDesc::get(i + 1));
        execTimesTable.setColLabel((size_t)i, FunctionDesc::get(i + 1));

        // Queue up a new function test for each iteration
        for (size_t iter = 0; iter < iterations; iter++)
        {
            mdata::TestParameters<T> curParam;
            curParam.funcId = i + 1;
            curParam.alpha = alpha;
            curParam.alg = testAlg;
            curParam.resultsTable = &resultsTable;
            curParam.execTimesTable = &execTimesTable;
            curParam.resultsCol = i;
            curParam.execTimesCol = i;
            curParam.resultsRow = iter;
            curParam.execTimesRow = iter;

            // Add function test to async queue
            testFutures.emplace_back(
                tPool->enqueue(&Experiment<T>::testFuncThreaded, this, curParam)
            );
        }
    }

    // Get the total number of async tasks queued
    const double totalFutures = static_cast<double>(testFutures.size());
    int tensPercentile = -1;
    std::chrono::microseconds waitTime(100);

    // Loop until all async tasks are completed and the thread futures
    // array is empty
    while (testFutures.size() > 0)
    {
        // Sleep a little bit since the async thread tasks are higher priority
        std::this_thread::sleep_for(waitTime);

        // Get iterator to first thread future
        auto it = testFutures.begin();

        // Loop through all thread futures
        while (it != testFutures.end())
        {
            if (!it->valid())
            {
                // An error occured with one of the threads
                cerr << "Error: Thread future invalid.";
                tPool->stopAndJoinAll();
                return 1;
            }

            // Get the status of the current thread future (async task)
            std::future_status status = it->wait_for(waitTime);
            if (status == std::future_status::ready)
            {
                // Task has completed, get return value
                int errCode = it->get();
                if (errCode)
                {
                    // An error occurred while running the task.
                    // Bail out of function
                    tPool->stopAndJoinAll();
                    return errCode;
                }

                // Remove processed task future from vector
                it = testFutures.erase(it);

                // Calculate the percent completed of all tasks, rounded to the nearest 10%
                int curPercentile = static_cast<int>(((totalFutures - testFutures.size()) / totalFutures) * 10);
                if (curPercentile > tensPercentile)
                {
                    // Print latest percent value to the console
                    tensPercentile = curPercentile;
                    cout << "~" << (tensPercentile * 10) << "% " << flush;
                }
            }
            else
            {
                // Async task has not yet completed, advance to the next one
                it++;
            }
        }
    }
    
    // Record total execution time and print it to the console
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    long double totalExecTime = static_cast<long double>(duration_cast<nanoseconds>(t_end - t_start).count()) / 1000000000.0L;

    cout << endl << "Test finished. Total time: " << std::setprecision(7) << totalExecTime << " seconds." << endl;

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

    cout << flush; */

    return 0;
}

/**
 * @brief Executes a single iteration of a test with the given parameters
 * 
 * @tparam T The data type used by the test
 * @param tParams The parameters used to set up the test
 * @return int An error code if any
 */
template<class T>
int Experiment<T>::testFuncThreaded(mdata::TestParameters<T> tParams)
{
    /* mdata::SearchAlgorithm<T>* alg;

    // Construct a search algorithm object for the selected alg
    switch (tParams.alg)
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

    // Retrieve the function bounds
    const RandomBounds<T>& funcBounds = vBounds[tParams.funcId - 1];

    // Retrieve the next available population object from the population pool
    mdata::Population<T>* pop = popPoolRemove();

    // Run the search algorithm one and record the results
    auto tResult = alg->run(Functions<T>::get(tParams.funcId), funcBounds.min, funcBounds.max, pop, tParams.alpha);

    // Place the population object back into the pool to be reused by anther thread
    popPoolAdd(pop);

    if (tResult.err)
    {
        cerr << "Error while testing function " << tParams.funcId << endl;
        return tResult.err;
    }

    // Update results table and execution times table with algorithm results
    tParams.resultsTable->setEntry(tParams.resultsRow, tParams.resultsCol, tResult.fitness);
    tParams.execTimesTable->setEntry(tParams.execTimesRow, tParams.execTimesCol, tResult.execTime);

    delete alg;
    return 0; */

    return 0;
}

/**
 * @brief Removes a single Population object from the pool,
 * or blocks until one is available. This function is thread-safe.
 * 
 * @tparam T Data type used by the Experiment class
 * @return mdata::Population<T>* Pointer to the removed Population object
 */
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

/**
 * @brief Adds a previously removed Population object back into
 * the population pool to be reused by another thread. This function
 * is thread-safe.
 * 
 * @tparam T Data type used by the Experiment class
 * @param popPtr Pointer to the Population object
 */
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
 * @brief Allocates all Population objects, which are used to store
 * population vectors and their associated fitness values. One population
 * object is created for each thread in the thread pool.
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
 * @brief Releases memory used by Population object pool.
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

/**
 * @brief Allocated the thread pool, which hosts all worker
 * threads and distributes all async tasks
 * 
 * @param numThreads Number of threads to create in the pool
 * @return Returns true if the thread pool was successfully created.
 * Otherwise returns false.
 */
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

// Explicit template specializations due to separate implementations in this CPP file
template class mfunc::Experiment<float>;
template class mfunc::Experiment<double>;
template class mfunc::Experiment<long double>;

// =========================
// End of experiment.cpp
// =========================
