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
    template <class T = double>
    void releaseArray(T*& a)
    {
        if (a == nullptr) return;

        delete[] a;
        a = nullptr;
    }

    template <class T = double>
    void releaseMatrix(T**& m, size_t rows)
    {
        if (m == nullptr) return;

        for (size_t i = 0; i < rows; i++)
        {
            if (m[i] != nullptr)
            {
                releaseArray<T>(m[i]);
                m[i] = nullptr;
            }
        }

        delete[] m;
        m = nullptr;
    }

    template <class T = double>
    inline T* allocArray(size_t size)
    {
        return new(std::nothrow) T[size];
    }

    template <class T = double>
    inline T** allocMatrix(size_t rows, size_t cols)
    {
        T** m = (T**)allocArray<T*>(rows);
        if (m == nullptr) return nullptr;

        for (size_t i = 0; i < rows; i++)
        {
            m[i] = allocArray<T>(cols);
            if (m[i] == nullptr)
            {
                releaseMatrix<T>(m, rows);
                return nullptr;
            }
        }

        return m;
    }
}

// =========================
// End of mem.h
// =========================
