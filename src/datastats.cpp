#include "datastats.h"
#include <cmath>
#include <bits/stdc++.h>

double mdata::average(const std::vector<double>& v)
{
    size_t vSize = v.size();
    double sum = 0.0;

    for (size_t i = 0; i < vSize; i++)
        sum += v[i];

    return sum / vSize;
}

double mdata::standardDeviation(const std::vector<double>& v)
{
    size_t vSize = v.size();
    double mean = average(v);
    double sum = 0;

    for (size_t i = 0; i < vSize; i++)
    {
        double subMean = v[i] - mean;
        sum += subMean * subMean;
    }

    return sqrt(sum / vSize);
}

double mdata::range(const std::vector<double>& v)
{
    size_t vSize = v.size();
    double min = v[0];
    double max = v[0];

    for (size_t i = 1; i < vSize; i++)
    {
        double cur = v[i];

        if (cur < min) min = cur;

        if (cur > max) max = cur;
    }

    return max - min;
}

double mdata::median(const std::vector<double>& v)
{
    size_t vSize = v.size();
    double* vSorted = new double[vSize];
    double retVal = 0;

    for (int i = 0; i < vSize; i++)
        vSorted[i] = v[i];

    std::sort(vSorted, vSorted + vSize);

    if (vSize % 2 != 0)
    {
        // Odd number of values
        retVal = vSorted[vSize / 2];
    }
    else
    {
        // Even number of values
        double low = vSorted[(vSize / 2) - 1];
        double high = vSorted[vSize / 2];
        retVal = (high + low) / 2;
    }

    delete vSorted;
    return retVal;
}