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

#include <vector>

namespace mdata
{
    double average(const std::vector<double>& v);
    double standardDeviation(const std::vector<double>& v);
    double range(const std::vector<double>& v);
    double median(const std::vector<double>& v);
}

#endif

// =========================
// End of datastats.h
// =========================
