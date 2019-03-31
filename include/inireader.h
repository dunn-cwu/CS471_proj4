#ifndef __INIREADER_H
#define __INIREADER_H

#include <string>
#include <map>
#include <iostream>
#include <fstream>

namespace util
{
    class IniReader
    {
    public:
        IniReader();
        bool openFile(std::string filePath);
        bool sectionExists(std::string section);
        bool entryExists(std::string section, std::string entry);
        std::string getEntry(std::string section, std::string entry);
    protected:
        std::string file;
        std::map<std::string, std::map<std::string, std::string>> iniMap;

        bool parseFile();
        void parseEntry(const std::string& sectionName, const std::string& entry);
    };
}

#endif