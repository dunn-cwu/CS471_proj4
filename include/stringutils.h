#ifndef __STRINGUTILS_H
#define __STRINGUTILS_H

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

namespace util
{
    // =======================================================
    // The string functions below were written by Evan Teran
    // from Stack Overflow:
    // https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
    // =======================================================

    // trim from start (in place)
    static inline void s_ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                std::not1(std::ptr_fun<int, int>(std::isspace))));
    }

    // trim from end (in place)
    static inline void s_rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    }

    // trim from both ends (in place)
    static inline void s_trim(std::string &s) {
        s_ltrim(s);
        s_rtrim(s);
    }

    // trim from start (copying)
    static inline std::string s_ltrim_copy(std::string s) {
        s_ltrim(s);
        return s;
    }

    // trim from end (copying)
    static inline std::string s_rtrim_copy(std::string s) {
        s_rtrim(s);
        return s;
    }

    // trim from both ends (copying)
    static inline std::string s_trim_copy(std::string s) {
        s_trim(s);
        return s;
    }
}
#endif