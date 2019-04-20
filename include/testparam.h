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
        unsigned int resultsCol;
        unsigned int execTimesCol;
        DataTable<T>* resultsTable;
        DataTable<T>* execTimesTable;
        size_t iterations;
        enums::Algorithm alg;
        
        TestParameters()
        {
            funcId = 1;
            alg = enums::Algorithm::BlindSearch;
            iterations = 0;
            resultsTable = nullptr;
            execTimesTable = nullptr;
            resultsCol = 0;
            execTimesCol = 0;
        }
    };
}

#endif