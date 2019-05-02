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
#include <iomanip>
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
        /**
         * @brief Construct a new Data Table object
         * Throws std::length_error and std::bad_alloc
         * 
         * @param _rows Number of rows in table
         * @param _cols Number of columns in table
         */
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

        /**
         * @brief Destroy the Data Table object
         */
        ~DataTable()
        {
            util::releaseMatrix(dataMatrix, rows);
        }
        
        void clearData()
        {
            util::initMatrix<T>(dataMatrix, rows, cols, 0);
        }

        /**
         * @brief Gets the string label for the column with the given index
         * 
         * @param colIndex Index of the column
         * @return std::string String value of the column label
         */
        std::string getColLabel(size_t colIndex)
        {
            if (colIndex >= colLabels.size())
                throw std::out_of_range("Column index out of range");

            return colLabels[colIndex];
        }

        /**
         * @brief Sets the string label for the column with the given index
         * 
         * @param colIndex Index of the column
         * @param newLabel New string label for the column
         */
        void setColLabel(size_t colIndex, std::string newLabel)
        {
             if (colIndex >= colLabels.size())
                throw std::out_of_range("Column index out of range");

            colLabels[colIndex] = newLabel;
        }

        /**
         * @brief Returns the value in the table at the given row and column
         * 
         * @param row Row index of the table
         * @param col Column index of the table
         * @return T Value of the entry at the given row and column
         */
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

        /**
         * @brief Set the value for the table entry at the given row and column
         * 
         * @param row Row index of the table
         * @param col Column index of the table
         * @param val New value for the entry
         */
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

        /**
         * @brief Exports the contents of this DataTable to a .csv file
         * 
         * @param filePath Path to the file that will be filled with this table's values
         * @return true If the file was successfully written to
         * @return false If there was an error opening the file
         */
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
                    outFile << std::setprecision(8) << dataMatrix[r][c];
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
        T** dataMatrix; /** Matrix of table data values */
        
    };
} // mdata

#endif

// =========================
// End of datatable.h
// =========================
