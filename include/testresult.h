#ifndef __TESTRESULT_H
#define __TESTRESULT_H

namespace mdata
{
    struct TestResult
    {
        const int err;
        const double execTime;

        TestResult(int _err, double _execTime) : err(_err), execTime(_execTime)
        {
        }
    };
} // mdata

#endif
