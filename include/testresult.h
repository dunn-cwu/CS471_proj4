/**
 * @file testresult.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Simple structure that packs together various return values
 * for the search algorithms.
 * functions.
 * @version 0.1
 * @date 2019-04-19
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __TESTRESULT_H
#define __TESTRESULT_H

namespace mdata
{
    template<class T>
    struct TestResult
    {
        const int err; // Error code. 0 = no error.
        const T fitness; // Fitness result
        const double execTime; // Algorithm execution time in miliseconds

        TestResult(int _err, T _fitness, double _execTime) : err(_err), fitness(_fitness), execTime(_execTime)
        {
        }
    };
} // mdata

#endif

// =========================
// End of testresult.h
// =========================
