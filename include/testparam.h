#ifndef __TESTPARAM_H
#define __TESTPARAM_H

#include <cstddef> // size_t definition
#include "datatable.h"
#include "searchalg.h"

namespace mdata
{
    template <class T>
    struct TestParameters
    {
        unsigned int funcId;
        T alpha;
        unsigned int resultsCol;
        unsigned int execTimesCol;
        size_t resultsRow;
        size_t execTimesRow;
        DataTable<T>* resultsTable;
        DataTable<T>* execTimesTable;
        enums::Algorithm alg;
        
        TestParameters()
        {
            funcId = 1;
            alpha = 0;
            alg = enums::Algorithm::BlindSearch;
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
