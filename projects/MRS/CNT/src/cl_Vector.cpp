/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Vector.cpp
 *
 */

#include "cl_Vector.hpp"

// C++ header files.
#include <vector>
#include <algorithm>    // for unique
#include <iostream>

// MORIS library header files.
#include "moris_typedefs.hpp"    // COR/src
#include "assert.hpp"

namespace moris
{
    //------------------------------------------------------------------

    inline Vector< char >
    string_to_char( Vector< std::string >& strings )
    {
        Vector< char > cstrings;
        cstrings.reserve( strings.size() );
        for ( const std::string& s : strings )
        {
            for ( size_t i = 0; i < strlen( s.c_str() ); ++i )
            {
                cstrings.push_back( s.c_str()[ i ] );
            }
        }

        return cstrings;
    }

    //------------------------------------------------------------------

}    // namespace moris
