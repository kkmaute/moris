/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Cell.cpp
 *
 */

#include "cl_Cell.hpp"

// C++ header files.
#include <vector>
#include <algorithm> // for unique
#include <iostream>

// MORIS library header files.
#include "moris_typedefs.hpp" // COR/src
#include "assert.hpp"

namespace moris
{
    //------------------------------------------------------------------

    inline
    moris::Cell<char>
    string_to_char(moris::Cell<std::string>& strings)
    {
        moris::Cell<char> cstrings;
        cstrings.reserve(strings.size());
        for(std::string s: strings)
        {
            for(size_t i = 0; i < strlen(s.c_str()); ++i)
            {
                cstrings.push_back(s.c_str()[i]);
            }
        }

        return cstrings;
    }

    //------------------------------------------------------------------

} // namespace moris

