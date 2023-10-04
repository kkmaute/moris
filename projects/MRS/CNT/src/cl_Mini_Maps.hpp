/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Mini_Maps.hpp
 *
 */

#pragma once

// C++ header files.
#include <exception>
#include <map>
#include <iostream>

// MORIS library header files.
#include "typedefs.hpp"    // COR/src
#include "assert.hpp"

/**
 * @brief This is here to define how (!) SMALL (!) maps are defined 
 * An example application for this would be the global to local node index relation on a single element
 *
 */

namespace moris
{
    // using std::unordered_map
    template< typename A, typename B > using Mini_Map = std::unordered_map< A, B >;

    // using std::map
    // template< typename A, typename B > using Mini_Map = std::map< A, B >;

    // define the IndexMap for convenience as it is used very often
    typedef Mini_Map< moris::moris_index, moris::moris_index > IndexMap;

}    // namespace moris
