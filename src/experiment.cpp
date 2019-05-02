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
#define INI_FUNC_RANGE_SECTION "function_range"
#define INI_GENALG_SECTION "genetic_alg"

#define INI_TEST_POPULATION "population"
#define INI_TEST_DIMENSIONS "dimensions"
#define INI_TEST_ITERATIONS "iterations"
#define INI_TEST_NUMTHREADS "num_threads"
#define INI_TEST_ALGORITHM "algorithm"
#define INI_TEST_RESULTSFILE "results_file"
#define INI_TEST_EXECTIMESFILE "exec_times_file"

#define INI_GENALG_GENERATIONS "generations"
#define INI_GENALG_CRPROB "crossover_prob"
#define INI_GENALG_MUTPROB "mutation_prob"
#define INI_GENALG_MUTRANGE "mutation_range"
#define INI_GENALG_MUTPREC "mutation_precision"
#define INI_GENALG_ELITISMRATE "elitism_rate"

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
        else if (selectedAlg >= static_cast<unsigned int>(Algorithm::Count))
        {
            cerr << "Experiment init failed: Param file [test]->" 
                << INI_TEST_ALGORITHM << " entry missing or out of bounds: " << paramFile << endl;
            return false;
        }

        // Cast iterations and test algorithm to correct types
        iterations = (size_t)numberIter;
        testAlg = static_cast<Algorithm>(selectedAlg);

        // Print test parameters to console
        cout << "Population size: " << numberSol << endl;
        cout << "Dimensions: " << numberDim << endl;
        cout << "Iterations: " << iterations << endl;
        // cout << "Algorithm: " << enums::AlgorithmNames::get(testAlg) << endl;

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
 * @brief Executes all functions as specified in the CS471 project 2
 * document, records results, and outputs the data as a *.csv file.
 * 
 * @return Returns 0 on success. Returns a non-zero error code on failure.
 */
template<class T>
int Experiment<T>::testAllFunc()
{
    switch (testAlg)
    {
        case Algorithm::GeneticAlgorithm:
            return testAllFunc_GA();
        default:
            return 1;
    }
}

template<class T>
int Experiment<T>::testAllFunc_GA()
{
    if (populationsPool.size() == 0) return 1;

    GAParams<T> _p;
    if (!loadGAParams(_p)) return 2;

    const GAParams<T>& paramTemplate = _p;

    mdata::DataTable<T> resultsTable(paramTemplate.generations, iterations);

    for (unsigned int c = 0; c < iterations; c++)
        resultsTable.setColLabel(c, std::string("Exp_") + std::to_string(c + 1));

    std::vector<std::future<int>> testFutures;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    for (unsigned int f = 1; f <= mfunc::NUM_FUNCTIONS; f++)
    {
        resultsTable.clearData();

        for (size_t exp = 0; exp < iterations; exp++)
        {
            GAParams<T> gaParams;
            gaParams.fitnessTable = &resultsTable;
            gaParams.fitTableCol = exp;
            gaParams.mainPop = nullptr;
            gaParams.auxPop = nullptr;
            gaParams.fPtr = Functions<T>::get(f);
            gaParams.fMinBound = vBounds[f-1].min;
            gaParams.fMaxBound = vBounds[f-1].max;
            gaParams.generations = paramTemplate.generations;
            gaParams.crProb = paramTemplate.crProb;
            gaParams.mutProb = paramTemplate.mutProb;
            gaParams.mutRange = paramTemplate.mutRange;
            gaParams.mutPrec = paramTemplate.mutPrec;
            gaParams.elitismRate = paramTemplate.elitismRate;

            testFutures.emplace_back(
                    tPool->enqueue(&Experiment<T>::runGAThreaded, this, gaParams)
                );
        }

        auto it = testFutures.begin();

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
                tPool->stopAndJoinAll();
                return errCode;
            }
        }

        testFutures.clear();

        std::string outFile = resultsFile;
        outFile = std::regex_replace(outFile, std::regex("\\%ALG%"), "GA");
        outFile = std::regex_replace(outFile, std::regex("\\%FUNC%"), std::to_string(f));

        if (!outFile.empty())
        {
            resultsTable.exportCSV(outFile.c_str());
            cout << "Exported function results to: " << outFile << endl << flush;
        }
    }

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    long double totalExecTime = static_cast<long double>(duration_cast<nanoseconds>(t_end - t_start).count()) / 1000000000.0L;

    cout << endl << "Test finished. Total time: " << std::setprecision(7) << totalExecTime << " seconds." << endl;

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
int Experiment<T>::runGAThreaded(GAParams<T> gaParams)
{
    // Retrieve the next two population objects from the population pool
    mdata::Population<T>* popMain = popPoolRemove();
    mdata::Population<T>* popAux = popPoolRemove();
    gaParams.mainPop = popMain;
    gaParams.auxPop = popAux;

    GeneticAlgorithm<T>* gaAlg = new GeneticAlgorithm<T>;
    int retVal = gaAlg->run(gaParams);
    delete gaAlg;

    popPoolAdd(popAux);
    popPoolAdd(popMain);

    return retVal;
}

template<class T>
bool Experiment<T>::loadGAParams(GAParams<T>& refParams)
{
    // Extract test parameters from ini file
    long generations = iniParams.getEntryAs<long>(INI_GENALG_SECTION, INI_GENALG_GENERATIONS);
    double crossover = iniParams.getEntryAs<double>(INI_GENALG_SECTION, INI_GENALG_CRPROB);
    double mutprob = iniParams.getEntryAs<double>(INI_GENALG_SECTION, INI_GENALG_MUTPROB);
    double mutrange = iniParams.getEntryAs<double>(INI_GENALG_SECTION, INI_GENALG_MUTRANGE);
    double mutprec = iniParams.getEntryAs<double>(INI_GENALG_SECTION, INI_GENALG_MUTPREC);
    double elitism = iniParams.getEntryAs<double>(INI_GENALG_SECTION, INI_GENALG_ELITISMRATE);

    // Verify test parameters
    if (generations <= 0)
    {
        cerr << "Experiment init failed: Param file [" << INI_GENALG_SECTION << "]->" 
            << INI_GENALG_GENERATIONS << " entry missing or out of bounds." << endl;
        return false;
    }
    else if (crossover <= 0)
    {
        cerr << "Experiment init failed: Param file [" << INI_GENALG_SECTION << "]->" 
            << INI_GENALG_CRPROB << " entry missing or out of bounds." << endl;
        return false;
    }
    else if (mutprob <= 0)
    {
        cerr << "Experiment init failed: Param file [" << INI_GENALG_SECTION << "]->" 
            << INI_GENALG_MUTPROB << " entry missing or out of bounds." << endl;
        return false;
    }
    else if (mutrange <= 0)
    {
        cerr << "Experiment init failed: Param file [" << INI_GENALG_SECTION << "]->" 
            << INI_GENALG_MUTRANGE << " entry missing or out of bounds." << endl;
        return false;
    }
    else if (mutprec <= 0)
    {
        cerr << "Experiment init failed: Param file [" << INI_GENALG_SECTION << "]->" 
            << INI_GENALG_MUTPREC << " entry missing or out of bounds." << endl;
        return false;
    }
    else if (elitism <= 0)
    {
        cerr << "Experiment init failed: Param file [" << INI_GENALG_SECTION << "]->" 
            << INI_GENALG_ELITISMRATE << " entry missing or out of bounds." << endl;
        return false;
    }

    refParams.generations = static_cast<unsigned int>(generations);
    refParams.crProb = crossover;
    refParams.mutProb = mutprob;
    refParams.mutRange = mutrange;
    refParams.mutPrec = mutprec;
    refParams.elitismRate = elitism;

    return true;
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
