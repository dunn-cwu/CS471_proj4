#include "inireader.h"
#include "stringutils.h"

using namespace util;

IniReader::IniReader() : file(""), iniMap()
{
}

bool IniReader::openFile(std::string filePath)
{
    file = filePath;
    if (!parseFile()) 
        return false;

    return true;
}

bool IniReader::sectionExists(std::string section)
{
    return iniMap.find(section) != iniMap.end();
}

bool IniReader::entryExists(std::string section, std::string entry)
{
    auto it = iniMap.find(section);
    if (it == iniMap.end()) return false;

    return it->second.find(entry) != it->second.end();
}

std::string IniReader::getEntry(std::string section, std::string entry)
{
    if (!entryExists(section, entry)) return std::string();

    return iniMap[section][entry];
}

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
            // Line is an entry
            parseEntry(curSection, line);
        }
    }

    inputF.close();
    return true;
}

void IniReader::parseEntry(const std::string& sectionName, const std::string& entry)
{
    using namespace std;

    const string delim = "=";
    string entryName;
    string entryValue;

    auto delimPos = entry.find(delim);
    if (delimPos == string::npos || delimPos >= entry.length() - 1) return;

    entryName = entry.substr((size_t)0, delimPos);
    entryValue = entry.substr(delimPos + 1, entry.length());

    s_trim(entryName);
    s_trim(entryValue);

    if (entryName.empty()) return;

    iniMap[sectionName][entryName] = entryValue;
}