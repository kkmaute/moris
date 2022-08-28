/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Array.cpp
 *
 */

#include <catch.hpp>

// MORIS project header files.
#include "cl_Array.hpp" // CON/src

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::Array" )
{
    SECTION( "moris::Array.begin()",
             "moris::Array.end()" )
    {
        moris::Array< moris::uint, 3 > myarray1 = { {10, 2, 3} };

        moris::Array< moris::uint, 3 > myarray2{ {10, 2, 3} };

        auto first1 = myarray1.begin();
        auto first2 = myarray2.begin();
        auto end1   = myarray1.end()-1;

        REQUIRE( *first1 == 10 );
        REQUIRE( *first2 == 10 );
        REQUIRE( *end1   == 3 );
    }

    SECTION( "moris::Array.size()",
             "moris::Array.empty()" )
    {
        moris::Array< moris::uint, 3 > myarray1 = { {10, 2, 3} };

        moris::Array< moris::uint, 0 > myarray2 = { {} };

        moris::size_t array_size = myarray1.size();

        bool isempty = myarray2.empty();

        REQUIRE( array_size == 3 );
        REQUIRE( isempty == true );
    }

    SECTION( "moris::Array.front()",
             "moris::Array.back()" )
    {
        moris::Array< moris::uint, 3 > myarray = { {10, 2, 3} };

        moris::uint front = myarray.front();
        moris::uint back  = myarray.back();

        REQUIRE( front == 10 );
        REQUIRE( back  == 3 );
    }

    SECTION( "moris::Array.fill( constant )" )
    {
        moris::Array< moris::uint, 5 > myarray1 = { {} };
        myarray1.fill(5);

        REQUIRE( myarray1[0] == 5 );
        REQUIRE( myarray1[4] == 5 );
    }
}

