/**
 * @file mem.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for various memory utility functions
 * @version 0.2
 * @date 2019-04-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __MEM_H
#define __MEM_H

#include <new> // std::nothrow
#include <cstddef> // size_t definition

namespace util
{
    /**
     * @brief Initializes an array with some set value
     * 
     * @tparam Data type of array
     * @param a Pointer to array
     * @param size Size of the array
     * @param val Value to initialize the array to
     */
    template <class T = double>
    inline void initArray(T* a, size_t size, T val)
    {
        if (a == nullptr) return;

        for (size_t i = 0; i < size; i++)
        {
            a[i] = val;
        }
    }

    /**
     * @brief Initializes a matrix with a set value for each entry
     * 
     * @tparam Data type of matrix entries
     * @param m Pointer to a matrix
     * @param rows Number of rows in matrix
     * @param cols Number of columns in matrix
     * @param val Value to initialize the matrix to
     */
    template <class T = double>
    inline void initMatrix(T** m, size_t rows, size_t cols, T val)
    {
        if (m == nullptr) return;

        for (size_t i = 0; i < rows; i++)
        {
            initArray(m[i], cols, val);
        }
    }

    /**
     * @brief Releases an allocated array's memory and sets the pointer to nullptr
     * 
     * @tparam Data type of array
     * @param a Pointer to array
     */
    template <class T = double>
    bool releaseArray(T*& a)
    {
        if (a == nullptr) return true;

        try
        {
            delete[] a;
            a = nullptr;
            return true;
        }
        catch(...)
        {
            return false;
        }
    }

    /**
     * @brief Releases an allocated matrix's memory and sets the pointer to nullptr
     * 
     * @tparam Data type of the matrix
     * @param m Pointer th the matrix
     * @param rows The number of rows in the matrix
     */
    template <class T = double>
    void releaseMatrix(T**& m, size_t rows)
    {
        if (m == nullptr) return;

        for (size_t i = 0; i < rows; i++)
        {
            if (m[i] != nullptr)
            {
                // Release each row
                releaseArray<T>(m[i]);
            }
        }

        // Release columns
        delete[] m;
        m = nullptr;
    }

    /**
     * @brief Allocates a new array of the given data type
     * 
     * @tparam Data type of the array
     * @param size Number of elements in the array
     * @return Returns a pointer to the new array, or nullptr allocation fails
     */
    template <class T = double>
    inline T* allocArray(size_t size)
    {
        return new(std::nothrow) T[size];
    }

    /**
     * @brief Allocates a new matrix of the given data type
     * 
     * @tparam Data type of the matrix entries
     * @param rows The number of rows
     * @param cols The number of columns
     * @return Returns a pointer to the new matrix, or nullptr if allocation fails
     */
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

    /**
     * @brief Copies the elements from one equal-sized array to another
     * 
     * @tparam Data type of the array
     * @param src Source array from where the elements will be copied from
     * @param dest Destination array from where the elements will be copied to
     * @param size Number of elements in the array
     */
    template <class T = double>
    inline void copyArray(T* src, T* dest, size_t size)
    {
        for (size_t i = 0; i < size; i++)
            dest[i] = src[i];
    }
}

#endif

// =========================
// End of mem.h
// =========================
