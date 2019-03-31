#ifndef __DATASTATS_H
#define __DATASTATS_H

#include <vector>
#include <cstddef> // size_t definition

namespace mdata
{
    double average(const std::vector<double>& v);

    double standardDeviation(const std::vector<double>& v);

    double range(const std::vector<double>& v);

    double median(const std::vector<double>& v);
}

#endif