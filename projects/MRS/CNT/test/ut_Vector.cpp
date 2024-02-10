/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Cell.cpp
 *
 */

#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "core.hpp"
#include "cl_Vector.hpp"    // CON/src

using namespace moris;

// ----------------------------------------------------------------------------

TEST_CASE( "moris::Vector" )
{

    SECTION( "Vector operator=" )
    {
        Vector< moris::uint > myval   = { 4 };
        Vector< moris::uint > myCell1 = { 10, 2, 3 };

        // Test with single value
        Vector< moris::uint > val1;
        val1 = myval;

        REQUIRE( val1( 0 ) == 4 );

        // Test with row Cell
        Vector< moris::uint > Cell1;
        Cell1 = myCell1;

        REQUIRE( Cell1( 0 ) == 10 );
        REQUIRE( Cell1( 1 ) == 2 );
        REQUIRE( Cell1( 2 ) == 3 );
    }

    SECTION( "Vector size constructor" )
    {
        Vector< moris::uint > myCell1( 10, 4 );
        Vector< moris::real > myCell2( 10, 3.0 );

        REQUIRE( myCell1( 5 ) == 4 );
        REQUIRE( myCell2( 5 ) == 3.0 );
    }

    SECTION( "Vector.size()",
            "Vector.empty()" )
    {
        Vector< moris::uint > myCell1 = { 10, 2, 3 };

        Vector< moris::uint > myCell2;

        moris::size_t Cell_size = myCell1.size();

        bool isempty = myCell2.empty();

        REQUIRE( Cell_size == 3 );
        REQUIRE( isempty == true );
    }

    SECTION( "Vector.append()" )
    {
        Vector< moris::uint > myCell1 = { 4, 5, 9, 8 };
        Vector< moris::uint > myCell2 = { 10, 2, 3 };

        myCell1.append( myCell2 );

        REQUIRE( myCell1( 0 ) == 4 );
        REQUIRE( myCell1( 1 ) == 5 );
        REQUIRE( myCell1( 2 ) == 9 );
        REQUIRE( myCell1( 3 ) == 8 );
        REQUIRE( myCell1( 4 ) == 10 );
        REQUIRE( myCell1( 5 ) == 2 );
        REQUIRE( myCell1( 6 ) == 3 );
    }


    SECTION( "Vector print_as_row_vector" )
    {
        Vector< moris::uint >                     myCell1 = { 4, 5, 9, 8 };
        Vector< Vector< moris::uint > >           myCell2 = { { 10, 2, 3 }, { 4, 5, 9, 8 } };
        Vector< Vector< Vector< moris::uint > > > myCell3 = { { { 1, 2, 3 }, { 4, 5, 6 } }, { { 7, 8, 9 } } };
        ;

        print_as_row_vector( myCell1 );
        print_as_row_vector( myCell2 );
        print_as_row_vector( myCell3 );
    }
}
