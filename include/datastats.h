/**
 * @file datastats.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for various data statistic functions
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __DATASTATS_H
#define __DATASTATS_H

#include <cmath>     // sqrt()
#include <algorithm> // std::sort()
#include <cstddef>   // size_t definition


namespace mdata
{
    /**
     * @brief Calculates the average for an array of values
     * 
     * @param v Array of values
     * @return The average value of the array
     */
    template<class T = double>
    T average(T* v, size_t vSize)
    {
        T sum = 0;

        for (size_t i = 0; i < vSize; i++)
            sum += v[i];

        return sum / vSize;
    }

    /**
     * @brief Calculates the standard deviation for an array of values
     * 
     * @param v Array of values
     * @return The standard deviation value of the array
     */
    template<class T = double>
    T standardDeviation(T* v, size_t vSize)
    {
        T mean = average<T>(v, vSize);
        T sum = 0;

        for (size_t i = 0; i < vSize; i++)
        {
            T subMean = v[i] - mean;
            sum += subMean * subMean;
        }

        return (T)sqrt((double)(sum / vSize));
    }

    /**
     * @brief Calculates the range for an array of values
     * 
     * @param v Array of values
     * @return The range of the array
     */
    template<class T = double>
    T range(T* v, size_t vSize)
    {
        T min = v[0];
        T max = v[0];

        for (size_t i = 1; i < vSize; i++)
        {
            T cur = v[i];

            if (cur < min) min = cur;

            if (cur > max) max = cur;
        }

        return max - min;
    }

    /**
     * @brief Calculates the median for an array of values
     * 
     * @param v Array of values
     * @return The median value of the array
     */
    template<class T = double>
    T median(T* v, size_t vSize)
    {
        T* vSorted = new T[vSize];
        T retVal = 0;

        for (size_t i = 0; i < vSize; i++)
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
            T low = vSorted[(vSize / 2) - 1];
            T high = vSorted[vSize / 2];
            retVal = (high + low) / 2;
        }

        delete[] vSorted;
        return retVal;
    }
}

#endif

// =========================
// End of datastats.h
// =========================
