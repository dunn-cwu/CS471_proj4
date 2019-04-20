#ifndef __MFUNCPTR_H
#define __MFUNCPTR_H

#include <cstddef> // size_t definition

namespace mfunc
{
    template <class T>
    using mfuncPtr = T (*)(T*, size_t);
}

#endif