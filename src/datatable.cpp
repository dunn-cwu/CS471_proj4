/**
 * @file datatable.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the DataTable class, which represents a
 * spreadsheet/table of values that can be exported to a *.csv file
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "datatable.h"
#include <iostream>
#include <fstream>

using namespace mdata;

/**
 * @brief Constructs a new DataTable object with a specified number of columns.
 * 
 * @param columns The number of columns to be created for the table
 */
DataTable::DataTable(unsigned int columns) : rows(0), cols(columns), colLabels(columns)
{
    for (int i = 0; i < cols; i++)
    {
        colLabels.push_back("(No label)");
    }
}

/**
 * @brief Destroys the DataTable object
 */
DataTable::~DataTable()
{
    colLabels.clear();
    tableData.clear();
}

/**
 * @brief Returns the label for the column at the specified index.
 * The first column = index 0.
 * 
 * @param colIndex Column index
 * @return A std::string containing the column label
 */
std::string DataTable::getColLabel(unsigned int colIndex)
{
    if (colIndex >= cols) throw "Invalid Column Index";

    return colLabels[colIndex];
}

/**
 * @brief Sets the label for the column at the specified index.
 * The first column = index 0.
 * 
 * @param colIndex Column index
 * @param newLabel std::string containing the new column label
 * @return true If the column label was succesfully changed.
 * @return false If the column index was invalid.
 */
bool DataTable::setColLabel(unsigned int colIndex, std::string newLabel)
{
    if (colIndex >= cols) return false;

    colLabels[colIndex] = newLabel;
    return true;
}

/**
 * @brief Adds a new row to the end of the table
 * 
 * @return The index of the newly added row
 */
unsigned int DataTable::addRow()
{
    unsigned int newRowIndex = rows;
    rows++;

    auto& tableRow = tableData[newRowIndex];
    tableRow.clear();

    for (int i = 0; i < cols; i++)
    {
        tableRow.push_back("");
    }

    return newRowIndex;
}

/**
 * @brief Adds a new row to the end of the table and fills the row
 * with the data given in the vector of strings rowData.
 * 
 * @param rowData Vector of strings to be entered into the table. rowData[n] = Column[n]
 * @return The index of the newly added row
 */
unsigned int DataTable::addRow(const std::vector<std::string>& rowData)
{
    unsigned int newRowIndex = rows;
    rows++;
    setRow(newRowIndex, rowData);

    return newRowIndex;
}

/**
 * @brief Returns a reference to the string vector that contains the
 * entries for the given row index.
 * 
 * @param row Index of the row that you wish to access.
 * @return std::vector<std::string>& 
 */
std::vector<std::string>& DataTable::getRow(unsigned int row)
{   
    if (row >= rows) throw "Invalid row index";

    return tableData[row];
}

/**
 * @brief Sets the data entries for the row at the given index.
 * 
 * @param row Index of the row that you wish to update.
 * @param rowData Vector of strings that contain the new row data entries.
 */
void DataTable::setRow(unsigned int row, const std::vector<std::string>& rowData)
{
    if (row >= rows) throw "Invalid row index";

    auto& tableRow = tableData[row];
    tableRow.clear();

    for (unsigned int i = 0; i < cols; i++)
    {
        if (i < rowData.size())
            tableRow.push_back(rowData[i]);
        else
            tableRow.push_back("(No data)");
    }
}

/**
 * @brief Returns the string value of the entry at the given row
 * and column indices.
 * 
 * @param row Index of the row you wish to access
 * @param col Index of the column you wish to access
 * @return The value of the given row and column. Throws a string
 * exception of the row or column is out of bounds.
 */
std::string DataTable::getEntry(unsigned int row, unsigned int col)
{
    if (row >= rows || col >= cols) throw "Invalid row/column";

    return tableData[row][col];
}

/**
 * @brief Sets the value of the entry at the given row and column indices.
 * 
 * @param row Index of the row you wish to access
 * @param col Index of the column you wish to access
 * @param val The new value for the entry
 */
void DataTable::setEntry(unsigned int row, unsigned int col, std::string val)
{
    if (row >= rows || col >= cols) throw "Invalid row/column";

    tableData[row][col] = val;
}

/**
 * @brief Sets the value of the entry at the given row and column indices.
 * 
 * @param row Index of the row you wish to access
 * @param col Index of the column you wish to access
 * @param val The new value for the entry
 */
void DataTable::setEntry(unsigned int row, unsigned int col, int val)
{
    setEntry(row, col, std::to_string(val));
}

/**
 * @brief Sets the value of the entry at the given row and column indices.
 * 
 * @param row Index of the row you wish to access
 * @param col Index of the column you wish to access
 * @param val The new value for the entry
 */
void DataTable::setEntry(unsigned int row, unsigned int col, long val)
{
    setEntry(row, col, std::to_string(val));
}

/**
 * @brief Sets the value of the entry at the given row and column indices.
 * 
 * @param row Index of the row you wish to access
 * @param col Index of the column you wish to access
 * @param val The new value for the entry
 */
void DataTable::setEntry(unsigned int row, unsigned int col, float val)
{
    setEntry(row, col, std::to_string(val));
}

/**
 * @brief Sets the value of the entry at the given row and column indices.
 * 
 * @param row Index of the row you wish to access
 * @param col Index of the column you wish to access
 * @param val The new value for the entry
 */
void DataTable::setEntry(unsigned int row, unsigned int col, double val)
{
    setEntry(row, col, std::to_string(val));
}

/**
 * @brief Exports the current data table to the given file path
 * in the *.csv format. If the file already exists, it is replaced.
 * 
 * @param filePath File path to be exported to
 * @return Returns true if the file was succesfully exported. Otherwise false.
 */
bool DataTable::exportCSV(const char* filePath)
{
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
            outFile << tableData[r][c];
            if (c < cols - 1) outFile << ",";
        }
        outFile << endl;
    }

    outFile.close();
    return true;
}

// =========================
// End of datatable.cpp
// =========================
