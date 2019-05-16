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
#include <regex>
#include "experiment.h"
#include "datatable.h"
#include "stringutils.h"
#include "mem.h"

// Ini file string sections and keys
#define INI_TEST_SECTION "test"
#define INI_PSO_SECTION "particle_swarm"
#define INI_FA_SECTION "firefly"
#define INI_HS_SECTION "harmony_search"
#define INI_FUNC_RANGE_SECTION "function_range"

#define INI_TEST_POPULATION "population"
#define INI_TEST_DIMENSIONS "dimensions"
#define INI_TEST_ITERATIONS "iterations"
#define INI_TEST_NUMTHREADS "num_threads"
#define INI_TEST_ALGORITHM "algorithm"
#define INI_TEST_RESULTSFILE "results_file"
#define INI_TEST_WORSTFITNESSFILE "worst_fit_file"
#define INI_TEST_EXECTIMESFILE "exec_times_file"
#define INI_TEST_FUNCCALLSFILE "func_calls_file"
#define INI_TEST_POPULATIONFILE "population_file"

#define INI_PSO_C1 "c1"
#define INI_PSO_C2 "c2"
#define INI_PSO_K "k"

#define INI_FA_ALPHA "alpha"
#define INI_FA_BETAMIN "betamin"
#define INI_FA_GAMMA "gamma"

#define INI_HS_HMCR "hmcr"
#define INI_HS_PAR "par"
#define INI_HS_BW "bw"

#define PARAM_DEFAULT_PSO_C1 0.8
#define PARAM_DEFAULT_PSO_C2 1.2
#define PARAM_DEFAULT_PSO_K 1.0

#define PARAM_DEFAULT_FA_ALPHA 0.5
#define PARAM_DEFAULT_FA_BETAMIN 0.2
#define PARAM_DEFAULT_FA_GAMMA 0.1

#define PARAM_DEFAULT_HS_HMCR 0.9
#define PARAM_DEFAULT_HS_PAR 0.4
#define PARAM_DEFAULT_HS_BW 0.2

#define RESULTSFILE_ALG_PATTERN "%ALG%"

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
        unsigned int selectedAlg = iniParams.getEntryAs<unsigned int>(INI_TEST_SECTION, INI_TEST_ALGORITHM);
        resultsFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_RESULTSFILE);
        worstFitnessFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_WORSTFITNESSFILE);
        execTimesFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_EXECTIMESFILE);
        funcCallsFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_FUNCCALLSFILE);
        populationsFile = iniParams.getEntry(INI_TEST_SECTION, INI_TEST_POPULATIONFILE);

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
        else if (selectedAlg >= static_cast<unsigned int>(Algorithm::Count))
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_ALGORITHM << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }

        // Cast iterations and test algorithm to correct types
        iterations = (size_t)numberIter;
        selAlg = static_cast<Algorithm>(selectedAlg);

        // Print test parameters to console
        cout << "Population size: " << numberSol << endl;
        cout << "Dimensions: " << numberDim << endl;
        cout << "Iterations: " << iterations << endl;

        // Allocate memory for all population objects. We need one for each thread to prevent conflicts.
        if (!allocatePopulationPool((size_t)numberThreads * 2, (size_t)numberSol, (size_t)numberDim))
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
 * @brief Executes all functions as specified in the CS471 project 3
 * document, records results, and outputs the data as a *.csv file.
 * 
 * @return Returns 0 on success. Returns a non-zero error code on failure.
 */
template<class T>
int Experiment<T>::testAllFunc()
{
    switch (selAlg)
    {
    case Algorithm::ParticleSwarm:
        return testPS();
        break;
    case Algorithm::Firefly:
        return testFA();
        break;
    case Algorithm::HarmonySearch:
        return testHS();
        break;
    default:
        cout << "Error: Invalid algorithm selected." << endl;
        break;
    }

    return 1;
}

template<class T>
int Experiment<T>::testPS()
{
    const PSParams<T> paramTemplate = createPSParamsTemplate();
    mdata::DataTable<T> resultsTable(iterations, 18);
    mdata::DataTable<T> worstTable(iterations, 18);
    mdata::DataTable<T> execTimesTable(1, 18);
    mdata::DataTable<T> funcCallsTable(1, 18);
    std::vector<std::future<int>> testFutures;

    mfunc::Functions<T>::resetCallCounters();

    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        auto desc = mfunc::FunctionDesc::get(f);
        resultsTable.setColLabel(f - 1, desc);
        worstTable.setColLabel(f - 1, desc);
        execTimesTable.setColLabel(f - 1, desc);
        funcCallsTable.setColLabel(f - 1, desc);

        PSParams<T> params(paramTemplate);
        params.popFile = util::s_replace(populationsFile, "%FUNC%", std::to_string(f));
        params.bestFitnessTable = &resultsTable;
        params.worstFitnessTable = &worstTable;
        params.fitTableCol = f - 1;
        params.mainPop = nullptr;
        params.pbPop = nullptr;
        params.fPtr = mfunc::Functions<T>::get(f);
        params.fMinBound = vBounds[f-1].min;
        params.fMaxBound = vBounds[f-1].max;
        params.iterations = iterations;

        testFutures.emplace_back(
                tPool->enqueue(&Experiment<T>::runPSThreaded, this, params, &execTimesTable, 0, f - 1)
        );
    }

    cout << "Executing particle swarm ..." << endl << flush;

    // Wait for threads to finish running all functions
    waitThreadFutures(testFutures);
    testFutures.clear();

    cout << endl;

    if (!funcCallsFile.empty())
    {
        for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
            funcCallsTable.setEntry(0, f - 1, mfunc::Functions<T>::getCallCounter(f));

        std::string outFile = util::s_replace(funcCallsFile, RESULTSFILE_ALG_PATTERN, "PSO");
        if (funcCallsTable.exportCSV(outFile.c_str()))
            cout << "Function call counts written to: " << outFile << endl;
        else
            cout << "Unable to function call counts file: " << outFile << endl;
    }

    if (!resultsFile.empty())
    {
        std::string outFile = util::s_replace(resultsFile, RESULTSFILE_ALG_PATTERN, "PSO");
        if (resultsTable.exportCSV(outFile.c_str()))
            cout << "Best fitness results written to: " << outFile << endl;
        else
            cout << "Unable to open results file: " << outFile << endl;
    }

    if (!worstFitnessFile.empty())
    {
        std::string outFile = util::s_replace(worstFitnessFile, RESULTSFILE_ALG_PATTERN, "PSO");
        if (worstTable.exportCSV(outFile.c_str()))
            cout << "Worst fitness results written to: " << outFile << endl;
        else
            cout << "Unable to open worst fitness file: " << outFile << endl;
    }

    if (!execTimesFile.empty())
    {
        std::string outFile = util::s_replace(execTimesFile, RESULTSFILE_ALG_PATTERN, "PSO");
        if (execTimesTable.exportCSV(outFile.c_str()))
            cout << "Execution times written to: " << outFile << endl;
        else
            cout << "Unable to open execution times file: " << outFile << endl;
    }

    return 0;
}

template<class T>
int Experiment<T>::runPSThreaded(PSParams<T> params, mdata::DataTable<T>* timesTable, size_t tRow, size_t tCol)
{
    auto mainPop = popPoolRemove();
    auto pbPop = popPoolRemove();
    params.mainPop = mainPop;
    params.pbPop = pbPop;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    ParticleSwarm<T> pswarm;
    int ret = pswarm.run(params);

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    double execTimeMs = static_cast<double>(duration_cast<nanoseconds>(t_end - t_start).count()) / 1000000.0;

    // Record execution time
    if (timesTable != nullptr)
        timesTable->setEntry(tRow, tCol, execTimeMs);

    popPoolAdd(mainPop);
    popPoolAdd(pbPop);
    return ret;
}

template<class T>
int Experiment<T>::testFA()
{
    const FAParams<T> paramTemplate = createFAParamsTemplate();
    mdata::DataTable<T> resultsTable(iterations, 18);
    mdata::DataTable<T> worstTable(iterations, 18);
    mdata::DataTable<T> execTimesTable(1, 18);
    mdata::DataTable<T> funcCallsTable(1, 18);
    std::vector<std::future<int>> testFutures;

    mfunc::Functions<T>::resetCallCounters();

    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        auto desc = mfunc::FunctionDesc::get(f);
        resultsTable.setColLabel(f - 1, desc);
        worstTable.setColLabel(f - 1, desc);
        execTimesTable.setColLabel(f - 1, desc);
        funcCallsTable.setColLabel(f - 1, desc);

        FAParams<T> params(paramTemplate);
        params.popFile = util::s_replace(populationsFile, "%FUNC%", std::to_string(f));
        params.bestFitnessTable = &resultsTable;
        params.worstFitnessTable = &worstTable;
        params.fitTableCol = f - 1;
        params.mainPop = nullptr;
        params.fPtr = mfunc::Functions<T>::get(f);
        params.fMinBound = vBounds[f-1].min;
        params.fMaxBound = vBounds[f-1].max;
        params.iterations = iterations;

        testFutures.emplace_back(
                tPool->enqueue(&Experiment<T>::runFAThreaded, this, params, &execTimesTable, 0, f - 1)
        );
    }

    cout << "Executing firefly ..." << endl << flush;

    // Wait for all threads to finish
    waitThreadFutures(testFutures);
    testFutures.clear();

    cout << endl;

    if (!funcCallsFile.empty())
    {
        for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
            funcCallsTable.setEntry(0, f - 1, mfunc::Functions<T>::getCallCounter(f));

        std::string outFile = util::s_replace(funcCallsFile, RESULTSFILE_ALG_PATTERN, "FA");
        if (funcCallsTable.exportCSV(outFile.c_str()))
            cout << "Function call counts written to: " << outFile << endl;
        else
            cout << "Unable to function call counts file: " << outFile << endl;
    }

    if (!resultsFile.empty())
    {
        std::string outFile = util::s_replace(resultsFile, RESULTSFILE_ALG_PATTERN, "FA");
        if (resultsTable.exportCSV(outFile.c_str()))
            cout << "Best fitness results written to: " << outFile << endl;
        else
            cout << "Unable to open results file: " << outFile << endl;
    }

    if (!worstFitnessFile.empty())
    {
        std::string outFile = util::s_replace(worstFitnessFile, RESULTSFILE_ALG_PATTERN, "FA");
        if (worstTable.exportCSV(outFile.c_str()))
            cout << "Worst fitness results written to: " << outFile << endl;
        else
            cout << "Unable to open worst fitness file: " << outFile << endl;
    }

    if (!execTimesFile.empty())
    {
        std::string outFile = util::s_replace(execTimesFile, RESULTSFILE_ALG_PATTERN, "FA");
        if (execTimesTable.exportCSV(outFile.c_str()))
            cout << "Execution times written to: " << outFile << endl;
        else
            cout << "Unable to open execution times file: " << outFile << endl;
    }

    return 0;
}

template<class T>
int Experiment<T>::runFAThreaded(FAParams<T> params, mdata::DataTable<T>* timesTable, size_t tRow, size_t tCol)
{
    auto mainPop = popPoolRemove();
    auto nextPop = popPoolRemove();
    params.mainPop = mainPop;
    params.nextPop = nextPop;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    Firefly<T> ffly;
    int ret = ffly.run(params);

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    double execTimeMs = static_cast<double>(duration_cast<nanoseconds>(t_end - t_start).count()) / 1000000.0;

    // Record execution time
    if (timesTable != nullptr)
        timesTable->setEntry(tRow, tCol, execTimeMs);

    popPoolAdd(mainPop);
    popPoolAdd(nextPop);
    return ret;
}

template<class T>
int Experiment<T>::testHS()
{
    const HSParams<T> paramTemplate = createHSParamsTemplate();
    mdata::DataTable<T> resultsTable(iterations, 18);
    mdata::DataTable<T> worstTable(iterations, 18);
    mdata::DataTable<T> execTimesTable(1, 18);
    mdata::DataTable<T> funcCallsTable(1, 18);
    std::vector<std::future<int>> testFutures;

    mfunc::Functions<T>::resetCallCounters();

    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        auto desc = mfunc::FunctionDesc::get(f);
        resultsTable.setColLabel(f - 1, desc);
        worstTable.setColLabel(f - 1, desc);
        execTimesTable.setColLabel(f - 1, desc);
        funcCallsTable.setColLabel(f - 1, desc);

        HSParams<T> params(paramTemplate);
        params.popFile = util::s_replace(populationsFile, "%FUNC%", std::to_string(f));
        params.bestFitnessTable = &resultsTable;
        params.worstFitnessTable = &worstTable;
        params.fitTableCol = f - 1;
        params.mainPop = nullptr;
        params.fPtr = mfunc::Functions<T>::get(f);
        params.fMinBound = vBounds[f-1].min;
        params.fMaxBound = vBounds[f-1].max;
        params.iterations = iterations;

        testFutures.emplace_back(
                tPool->enqueue(&Experiment<T>::runHSThreaded, this, params, &execTimesTable, 0, f - 1)
        );
    }

    cout << "Executing harmony search ..." << endl << flush;

    waitThreadFutures(testFutures);

    // Clear thread futures
    testFutures.clear();

    cout << endl;

    if (!funcCallsFile.empty())
    {
        for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
            funcCallsTable.setEntry(0, f - 1, mfunc::Functions<T>::getCallCounter(f));

        std::string outFile = util::s_replace(funcCallsFile, RESULTSFILE_ALG_PATTERN, "HS");
        if (funcCallsTable.exportCSV(outFile.c_str()))
            cout << "Function call counts written to: " << outFile << endl;
        else
            cout << "Unable to function call counts file: " << outFile << endl;
    }

    if (!resultsFile.empty())
    {
        std::string outFile = util::s_replace(resultsFile, RESULTSFILE_ALG_PATTERN, "HS");
        if (resultsTable.exportCSV(outFile.c_str()))
            cout << "Best fitness results written to: " << outFile << endl;
        else
            cout << "Unable to open results file: " << outFile << endl;
    }

    if (!worstFitnessFile.empty())
    {
        std::string outFile = util::s_replace(worstFitnessFile, RESULTSFILE_ALG_PATTERN, "HS");
        if (worstTable.exportCSV(outFile.c_str()))
            cout << "Worst fitness results written to: " << outFile << endl;
        else
            cout << "Unable to open worst fitness file: " << outFile << endl;
    }

    if (!execTimesFile.empty())
    {
        std::string outFile = util::s_replace(execTimesFile, RESULTSFILE_ALG_PATTERN, "HS");
        if (execTimesTable.exportCSV(outFile.c_str()))
            cout << "Execution times written to: " << outFile << endl;
        else
            cout << "Unable to open execution times file: " << outFile << endl;
    }

    return 0;
}

template<class T>
int Experiment<T>::runHSThreaded(HSParams<T> params, mdata::DataTable<T>* timesTable, size_t tRow, size_t tCol)
{
    auto mainPop = popPoolRemove();
    params.mainPop = mainPop;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    HarmonySearch<T> hsearch;
    int ret = hsearch.run(params);

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    double execTimeMs = static_cast<double>(duration_cast<nanoseconds>(t_end - t_start).count()) / 1000000.0;

    // Record execution time
    if (timesTable != nullptr)
        timesTable->setEntry(tRow, tCol, execTimeMs);

    popPoolAdd(mainPop);
    return ret;
}

template<class T>
int Experiment<T>::waitThreadFutures(std::vector<std::future<int>>& testFutures)
{
    cout << "Waiting for threads to finish ..." << endl << flush;

    const size_t totalFutures = testFutures.size();

    // Join all thread futures and get result
    for (size_t futIndex = 0; futIndex < testFutures.size(); futIndex++)
    {
        auto& curFut = testFutures[futIndex];

        if (!curFut.valid())
        {
            // An error occured with one of the threads
            cerr << "Error: Thread future invalid.";
            tPool->stopAndJoinAll();
            return 1;
        }

        int errCode = curFut.get();
        if (errCode)
        {
            // An error occurred while running the task.
            // Bail out of function
            cerr << "Error: Threaded function returned error code: " << errCode << endl;
            tPool->stopAndJoinAll();
            return errCode;
        }

        cout << futIndex + 1 << ".." << flush;
    }

    return 0;
}

template<class T>
const PSParams<T> Experiment<T>::createPSParamsTemplate()
{
    PSParams<T> retParams;

    retParams.c1 = iniParams.getEntryAs<double>(INI_PSO_SECTION, INI_PSO_C1, PARAM_DEFAULT_PSO_C1);
    retParams.c2 = iniParams.getEntryAs<double>(INI_PSO_SECTION, INI_PSO_C2, PARAM_DEFAULT_PSO_C2);
    retParams.k = iniParams.getEntryAs<double>(INI_PSO_SECTION, INI_PSO_K, PARAM_DEFAULT_PSO_K);

    return retParams;
}

template<class T>
const FAParams<T> Experiment<T>::createFAParamsTemplate()
{
    FAParams<T> retParams;

    retParams.alpha = iniParams.getEntryAs<double>(INI_FA_SECTION, INI_FA_ALPHA, PARAM_DEFAULT_FA_ALPHA);
    retParams.betamin = iniParams.getEntryAs<double>(INI_FA_SECTION, INI_FA_BETAMIN, PARAM_DEFAULT_FA_BETAMIN);
    retParams.gamma = iniParams.getEntryAs<double>(INI_FA_SECTION, INI_FA_GAMMA, PARAM_DEFAULT_FA_GAMMA);

    return retParams;
}

template<class T>
const HSParams<T> Experiment<T>::createHSParamsTemplate()
{
    HSParams<T> retParams;

    retParams.hmcr = iniParams.getEntryAs<double>(INI_HS_SECTION, INI_HS_HMCR, PARAM_DEFAULT_HS_HMCR);
    retParams.par = iniParams.getEntryAs<double>(INI_HS_SECTION, INI_HS_PAR, PARAM_DEFAULT_HS_PAR);
    retParams.bw = iniParams.getEntryAs<double>(INI_HS_SECTION, INI_HS_BW, PARAM_DEFAULT_HS_BW);

    return retParams;
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
