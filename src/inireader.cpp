/**
 * @file inireader.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Implementation file for the IniReader class, which can 
 * open and parse simple *.ini files.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "inireader.h"
#include "stringutils.h"

using namespace util;

/**
 * @brief Construct a new IniReader object
 */
IniReader::IniReader() : file(""), iniMap()
{
}

/**
 * @brief Destroys the IniReader object
 */
IniReader::~IniReader()
{
    iniMap.clear();
}

/**
 * @brief Opens the given ini file and parses all sections/entries.
 * The all file data is stored in memory and the file is closed.
 * 
 * @param filePath Path to the ini file you wish to open
 * @return Returns true if the file was succesfully opened and parsed. Otherwise false.
 */
bool IniReader::openFile(std::string filePath)
{
    file = filePath;
    if (!parseFile()) 
        return false;

    return true;
}

/**
 * @brief Returns true if the given section exists in the current ini file.
 * 
 * @param section std::string containing the section name
 * @return Returns true if the section exists in the ini file, otherwise false.
 */
bool IniReader::sectionExists(std::string section)
{
    return iniMap.find(section) != iniMap.end();
}

/**
 * @brief Returns true if the given section and entry key exists in the current ini file.
 * 
 * @param section std::string containing the section name
 * @param entry std::string containing the entry key name
 * @return Returns true if the section and entry key exist in the ini file, otherwise false.
 */
bool IniReader::entryExists(std::string section, std::string entry)
{
    auto it = iniMap.find(section);
    if (it == iniMap.end()) return false;

    return it->second.find(entry) != it->second.end();
}

/**
 * @brief Returns the value for the entry that has the given entry key
 * within the given section
 * 
 * @param section std::string containing the section name
 * @param entry std::string containing the entry key name
 * @return The value of the entry with the given entry key and section. 
 * Returns an empty string if the entry does not exist.
 */
std::string IniReader::getEntry(std::string section, std::string entry, std::string defVal)
{
    if (!entryExists(section, entry)) return defVal;

    return iniMap[section][entry];
}

/**
 * @brief Protected helper function that is called by IniReader::openFile().
 * Parses the complete ini file and stores all sections and entries in memory.
 * 
 * @return Returns true if the file was succesfully opened and parsed.
 */
bool IniReader::parseFile()
{
    iniMap.clear();

    using namespace std;

    ifstream inputF(file, ifstream::in);
    if (!inputF.good()) return false;

    string curSection;
    string line;

    while (getline(inputF, line))
    {
        // Trim whitespace on both ends of the line
        s_trim(line);

        // Ignore empty lines and comments
        if (line.empty() || line.front() == '#')
        {
            continue;
        }
        else if (line.front() == '[' && line.back() == ']')
        {
            // Line is a section definition
            // Erase brackets and trim to get section name
            line.erase(0, 1);
            line.erase(line.length() - 1, 1);
            s_trim(line);
            curSection = line;
        }
        else if (!curSection.empty())
        {
            // Line is an entry, parse the key and value
            parseEntry(curSection, line);
        }
    }

    // Close input file
    inputF.close();
    return true;
}

/**
 * @brief Protected helper function that is called by IniReader::parseFile().
 * Parses a single entry by extracting the key and value.
 */
void IniReader::parseEntry(const std::string& sectionName, const std::string& entry)
{
    using namespace std;

    // Split string around equals sign character
    const string delim = "=";
    string entryName;
    string entryValue;

    // Find index of '='
    auto delimPos = entry.find(delim);

    if (delimPos == string::npos || delimPos >= entry.length() - 1) 
        return; // '=' is missing, or is last char in string

    // Extract entry name/key and value
    entryName = entry.substr((size_t)0, delimPos);
    entryValue = entry.substr(delimPos + 1, entry.length());

    // Remove leading and trailing whitespace
    s_trim(entryName);
    s_trim(entryValue);

    // We cannot have entries with empty keys
    if (entryName.empty()) return;

    // Add entry to cache
    iniMap[sectionName][entryName] = entryValue;
}

// =========================
// End of inireader.cpp
// =========================
