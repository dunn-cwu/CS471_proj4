#include "datatable.h"

using namespace mdata;

DataTable::DataTable(unsigned int columns) : rows(0), cols(columns), colLabels(columns)
{
    for (int i = 0; i < cols; i++)
    {
        colLabels.push_back("(No label)");
    }
}

std::string DataTable::getColLabel(unsigned int colIndex)
{
    if (colIndex >= cols) throw "Invalid Column Index";

    return colLabels[colIndex];
}

bool DataTable::setColLabel(unsigned int colIndex, std::string newLabel)
{
    if (colIndex >= cols) return false;

    colLabels[colIndex] = newLabel;
    return true;
}

unsigned int DataTable::addRow()
{
    unsigned int newRowIndex = rows;
    rows++;

    auto& tableRow = tableData[newRowIndex];
    tableRow.clear();

    for (int i = 0; i < cols; i++)
    {
        tableRow.push_back("(No data)");
    }

    return newRowIndex;
}

unsigned int DataTable::addRow(const std::vector<std::string>& rowData)
{
    unsigned int newRowIndex = rows;
    rows++;
    setRow(newRowIndex, rowData);

    return newRowIndex;
}

std::vector<std::string>& DataTable::getRow(unsigned int row)
{   
    if (row >= rows) throw "Invalid row index";

    return tableData[row];
}

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

std::string DataTable::getEntry(unsigned int row, unsigned int col)
{
    if (row >= rows || col >= cols) throw "Invalid row/column";

    return tableData[row][col];
}

void DataTable::setEntry(unsigned int row, unsigned int col, std::string val)
{
    if (row >= rows || col >= cols) throw "Invalid row/column";

    tableData[row][col] = val;
}

