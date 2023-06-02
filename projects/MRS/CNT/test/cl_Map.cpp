/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map.cpp
 *
 */

#include <catch.hpp>

// MORIS project header files.
#include "core.hpp"
#include "algorithms.hpp"
#include "cl_Map.hpp" // CON/src

// ----------------------------------------------------------------------------

TEST_CASE( "moris::map" )
{
    #include "containers/cl_Map/cl_Map.inc"
    #include "containers/cl_Map/cl_Map_size.inc"

    SECTION("moris::map size")
    {
        REQUIRE( map_size == 3 );
    }

    SECTION("moris::map clear")
    {
        #include "containers/cl_Map/cl_Map_clear.inc"
        REQUIRE( myMap.size() == 2 );
    }

    SECTION("moris::map find")
    {
        #include "containers/cl_Map/cl_Map_clear.inc"
        #include "containers/cl_Map/cl_Map_find.inc"
        REQUIRE( aValue == 12 );
    }

    SECTION("moris::map empty")
    {
        #include "containers/cl_Map/cl_Map_empty.inc"
        REQUIRE( aEmpty == 0 );
    }
}

