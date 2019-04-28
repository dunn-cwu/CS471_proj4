/**
 * @file testparam.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the definition of the TestParameters struct,
 * which is a data type used to transfer test parameters between
 * functions.
 * @version 0.1
 * @date 2019-04-20
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __TESTPARAM_H
#define __TESTPARAM_H

#include <cstddef> // size_t definition
#include "datatable.h"

namespace mdata
{
    /**
     * @brief Packs together various test experiment parameters
     * 
     * @tparam T The data type that the test algorithm and functions are using
     */
    template <class T>
    struct TestParameters
    {
        unsigned int funcId; /** Id for the tested math function */
        T alpha; /** Alpha value if needed by the search alg */
        unsigned int resultsCol; /** DataTable column index to store the fitness result */
        unsigned int execTimesCol; /** DataTable column index to store the execution time */
        size_t resultsRow; /** DataTable row index to store the fitness result */
        size_t execTimesRow; /** DataTable row index to store the execution time */
        DataTable<T>* resultsTable; /** Pointer to the DataTable used to store the fitness results */
        DataTable<T>* execTimesTable; /** Pointer to the DataTable used to store the execution times */
        // enums::Algorithm alg; /** Selected search algorithm for this test */
        
        TestParameters()
        {
            funcId = 1;
            alpha = 0;
            resultsTable = nullptr;
            execTimesTable = nullptr;
            resultsCol = 0;
            execTimesCol = 0;
            resultsRow = 0;
            execTimesRow = 0;
        }
    };
}

#endif

// =========================
// End of testparam.h
// =========================
