/**
 * @file datatable.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the DataTable class, which represents a
 * spreadsheet/table of values that can be exported to a *.csv file
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __DATATABLE_H
#define __DATATABLE_H

#include <map>
#include <vector>
#include <string>

namespace mdata
{
    /**
     * @brief The DataTable class is a simple table of values with labeled columns.
     * 
     * --
     * Initialize a DataTable object with a specified number of columns:
     * DataTable table(n);
     * 
     * --
     * Set a column's label:
     * 
     * table.setColLabel(0, "Column 1");
     * 
     * --
     * Add a row to the table:
     * int rowIndex = table.addRow();
     * 
     * or
     * 
     * int rowIndex = table.addRow((std::vector<std::string>)dataVector);
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
    class DataTable
    {
    public:
        DataTable(unsigned int cols);
        ~DataTable();
        
        std::string getColLabel(unsigned int colIndex);
        bool setColLabel(unsigned int colIndex, std::string newLabel);

        unsigned int addRow();
        unsigned int addRow(const std::vector<std::string>& rowData);
        std::vector<std::string>& getRow(unsigned int row);
        void setRow(unsigned int row, const std::vector<std::string>& rowData);

        std::string getEntry(unsigned int row, unsigned int col);
        void setEntry(unsigned int row, unsigned int col, std::string val);
        void setEntry(unsigned int row, unsigned int col, int val);
        void setEntry(unsigned int row, unsigned int col, long val);
        void setEntry(unsigned int row, unsigned int col, float val);
        void setEntry(unsigned int row, unsigned int col, double val);

        bool exportCSV(const char* filePath);
    private:
        unsigned int cols; /** Number of columns in the table. */
        unsigned int rows; /** Number of rows in the table. */
        std::vector<std::string> colLabels; /** Vector of column labels. Index n = Col n. */
        std::map<unsigned int, std::vector<std::string>> tableData; /** Table data composed of a map of vectors. map[n] = row vector. map[n][m] = entry. */
    };
} // mdata

#endif

// =========================
// End of datatable.h
// =========================
