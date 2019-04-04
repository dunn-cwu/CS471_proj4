/**
 * @file mem.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for various memory utility functions
 * @version 0.1
 * @date 2019-04-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <new> // std::nothrow
#include <cstddef> // size_t definition

namespace util
{
    template <class T>
    inline T* allocArray(size_t size)
    {
        return new(std::nothrow) T[size];
    }

    template <class T>
    inline T** allocMatrix(size_t rows, size_t cols)
    {
        T** m = (T**)allocArray<T*>(rows);
        if (m == nullptr) return nullptr;
    }
}

// =========================
// End of mem.h
// =========================
