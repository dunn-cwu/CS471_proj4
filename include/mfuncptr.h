/**
 * @file mfuncptr.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Contains the type definition for mfuncPtr, 
 * a templated function pointer to one of the math functions
 * in mfunctions.h
 * @version 0.1
 * @date 2019-04-19
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __MFUNCPTR_H
#define __MFUNCPTR_H

#include <cstddef> // size_t definition

namespace mfunc
{
    /**
     * @brief Function pointer that takes two arguments T* and size_t,
     * and returns a T value.
     * 
     * @tparam T Data type for vector and return value
     */
    template <class T>
    using mfuncPtr = T (*)(T*, size_t);
}

#endif

// =========================
// End of mfuncptr.h
// =========================
