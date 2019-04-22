/**
 * @file datatable.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the DataTable class, which represents a
 * spreadsheet/table of values that can easily be exported to a *.csv file
 * @version 0.2
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __DATATABLE_H
#define __DATATABLE_H

#include <vector>
#include <string>
#include <stdexcept>
#include "mem.h"

namespace mdata
{
    /**
     * @brief The DataTable class is a simple table of values with labeled columns.
     * 
     * --
     * Initialize a DataTable object with a specified number of rows and columns:
     * DataTable table(rows, columns);
     * 
     * --
     * Set a column's label:
     * 
     * table.setColLabel(0, "Column 1");
     * 
     * --
     * Set an entry in the table:
     * 
     * table.setEntry(n, m, value);
     * 
     * Where 'n' is the row, 'm' is the column, and 'value' is the value
     * of the entry
     * 
     * --
     * Export the table to a *.csv file:
     * 
     * bool success = table.exportCSV("my_file.csv");
     */
    template <class T>
    class DataTable
    {
    public:
        DataTable(size_t _rows, size_t _cols) : rows(_rows), cols(_cols), dataMatrix(nullptr)
        {
            if (rows == 0)
                throw std::length_error("Table rows must be greater than 0.");
            else if (cols == 0)
                throw std::length_error("Table columns must be greater than 0.");

            dataMatrix = util::allocMatrix<T>(rows, cols);
            if (dataMatrix == nullptr)
                throw std::bad_alloc();

            colLabels.resize(_cols, std::string());
        }

        ~DataTable()
        {
            util::releaseMatrix(dataMatrix, rows);
        }
        
        std::string getColLabel(size_t colIndex)
        {
            if (colIndex >= colLabels.size())
                throw std::out_of_range("Column index out of range");

            return colLabels[colIndex];
        }

        void setColLabel(size_t colIndex, std::string newLabel)
        {
             if (colIndex >= colLabels.size())
                throw std::out_of_range("Column index out of range");

            colLabels[colIndex] = newLabel;
        }

        T getEntry(size_t row, size_t col)
        {
            if (dataMatrix == nullptr)
                throw std::runtime_error("Data matrix not allocated");
            if (row >= rows)
                throw std::out_of_range("Table row out of range");
            else if (col >= cols)
                throw std::out_of_range("Table column out of range");

            return dataMatrix[row][col];
        }

        void setEntry(size_t row, size_t col, T val)
        {
            if (dataMatrix == nullptr)
                throw std::runtime_error("Data matrix not allocated");
            if (row >= rows)
                throw std::out_of_range("Table row out of range");
            else if (col >= cols)
                throw std::out_of_range("Table column out of range");

            dataMatrix[row][col] = val;
        }

        bool exportCSV(const char* filePath)
        {
            if (dataMatrix == nullptr) return false;

            using namespace std;
            ofstream outFile;
            outFile.open(filePath, ofstream::out | ofstream::trunc);
            if (!outFile.good()) return false;

            // Print column labels
            for (unsigned int c = 0; c < cols; c++)
            {
                outFile << colLabels[c];
                if (c < cols - 1) outFile << ",";
            }

            outFile << endl;

            // Print data rows
            for (unsigned int r = 0; r < rows; r++)
            {
                for (unsigned int c = 0; c < cols; c++)
                {
                    outFile << dataMatrix[r][c];
                    if (c < cols - 1) outFile << ",";
                }
                outFile << endl;
            }

            outFile.close();
            return true;
        }
    private:
        size_t rows; /** Number of rows in the table. */
        size_t cols; /** Number of columns in the table. */
        std::vector<std::string> colLabels; /** Vector of column labels. Index n = Col n. */
        T** dataMatrix;
        
    };
} // mdata

#endif

// =========================
// End of datatable.h
// =========================
