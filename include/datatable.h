#ifndef __DATATABLE_H
#define __DATATABLE_H

#include <map>
#include <vector>
#include <string>

namespace mdata
{
    class DataTable
    {
    public:
        DataTable(unsigned int cols);

        std::string getColLabel(unsigned int colIndex);
        bool setColLabel(unsigned int colIndex, std::string newLabel);

        unsigned int addRow();
        unsigned int addRow(const std::vector<std::string>& rowData);
        std::vector<std::string>& getRow(unsigned int row);
        void setRow(unsigned int row, const std::vector<std::string>& rowData);

        std::string getEntry(unsigned int row, unsigned int col);
        void setEntry(unsigned int row, unsigned int col, std::string val);
    protected:
        unsigned int cols;
        unsigned int rows;
        std::vector<std::string> colLabels;
        std::map<unsigned int, std::vector<std::string>> tableData;
    };
} // mdata

#endif