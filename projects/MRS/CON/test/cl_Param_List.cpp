/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Param_List.cpp
 *
 */

#include <catch.hpp>

// MORIS project header files.
#include "cl_Param_List.hpp" // CON/src

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::Param_List" )
{
    SECTION( "moris::Param_List.insert( aKey, aVal )" )
    {
        // test the ability to create maps of different types
        #include "containers/cl_param_list/cl_param_list_insert.inc"

        REQUIRE( list.get< moris::sint >( "max_its" )   == 1000 );
        REQUIRE( list.get< moris::real >( "step_size" ) == 0.01 );
        REQUIRE( list.get< bool        >( "is_true" )   == true );

        // test the ability to change values of different types
        #include "../../../snippets/containers/cl_param_list/cl_param_list_access.inc"

        REQUIRE( list.get< moris::sint >( "max_its" )   == 100   );
        REQUIRE( list.get< moris::real >( "step_size" ) == 0.05  );
        REQUIRE( list.get< bool        >( "is_true" )   == false );
    }

    SECTION( "moris::Param_List.which( aKey )" )
    {
        #include "containers/cl_param_list/cl_param_list_insert.inc"

        // test the ability to get the underlying entry type
        REQUIRE( list.which( "max_its" )   == 1   );
        REQUIRE( list.which( "step_size" ) == 2   );
        REQUIRE( list.which( "is_true" )   == 0   );
    }
}

