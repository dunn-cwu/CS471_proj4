/**
 * @file inireader.h
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Header file for the IniReader class, which can 
 * open and parse simple *.ini files.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef __INIREADER_H
#define __INIREADER_H

#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <fstream>

namespace util
{
    /**
     * @brief The IniReader class is a simple *.ini file reader and parser
     * 
     * --
     * Initialize an IniReader object:
     * 
     * IniReader ini;
     * 
     * --
     * Open and parse an *.ini file:
     * 
     * ini.openFile("my_ini_file.ini");
     * 
     * Note that the file is immediately closed after parsing, and the
     * file data is retained in memory.
     * 
     * --
     * Retrieve an entry from the ini file:
     * 
     * std::string value = ini.getEntry("My Section", "entryKey");
     * 
     */
    class IniReader
    {
    public:
        IniReader();
        ~IniReader();
        bool openFile(std::string filePath);
        bool sectionExists(std::string section);
        bool entryExists(std::string section, std::string entry);
        std::string getEntry(std::string section, std::string entry, std::string defVal = "");

        template <class T>
        T getEntryAs(std::string section, std::string entry, T defVal = {})
        {
            std::stringstream ss(getEntry(section, entry, std::to_string(defVal)));
            T retVal;
            ss >> retVal;
            return retVal;
        }
    private:
        std::string file; /** The file path for the current ini file data. */
        std::map<std::string, std::map<std::string, std::string>> iniMap; /** The parsed ini file data. */

        bool parseFile();
        void parseEntry(const std::string& sectionName, const std::string& entry);
    };
}

#endif

// =========================
// End of inireader.h
// =========================
