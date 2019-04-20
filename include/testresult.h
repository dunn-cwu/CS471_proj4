#ifndef __TESTRESULT_H
#define __TESTRESULT_H

namespace mdata
{
    template<class T>
    struct TestResult
    {
        const int err;
        const T fitness;
        const double execTime;

        TestResult(int _err, T _fitness, double _execTime) : err(_err), fitness(_fitness), execTime(_execTime)
        {
        }
    };
} // mdata

#endif
