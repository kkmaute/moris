/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Vector.cpp
 *
 */

#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "core.hpp"
#include "cl_Vector.hpp" // CON/src

// ----------------------------------------------------------------------------

TEST_CASE(
        "Vector" )
{

    SECTION( "Vector operator=")
    {
        moris::Vector< moris::uint > myval = { 4 };
        moris::Vector< moris::uint > myCell1 = { 10, 2, 3 };

        // Test with single value
        moris::Vector< moris::uint > val1;
        val1 = myval;

        REQUIRE( val1(0) == 4 );

        // Test with row Cell
        moris::Vector< moris::uint > Cell1;
        Cell1 = myCell1;

        REQUIRE( Cell1(0) == 10 );
        REQUIRE( Cell1(1) == 2 );
        REQUIRE( Cell1(2) == 3 );
    }

    SECTION( "Vector size constructor" )
    {
        moris::Vector< moris::uint > myCell1( 10, 4  );
        moris::Vector< moris::real > myCell2( 10, 3.0 );

        REQUIRE( myCell1( 5 ) == 4   );
        REQUIRE( myCell2( 5 ) == 3.0 );
    }

    SECTION( "Vector.size()",
             "Vector.empty()" )
    {
        moris::Vector< moris::uint > myCell1 = {10, 2, 3};

        moris::Vector< moris::uint > myCell2;

        moris::size_t Cell_size = myCell1.size();

        bool isempty = myCell2.empty();

        REQUIRE( Cell_size == 3 );
        REQUIRE( isempty == true );
    }

    SECTION( "Vector.append()")
    {
        moris::Vector< moris::uint > myCell1 = { 4, 5, 9, 8 };
        moris::Vector< moris::uint > myCell2 = { 10, 2, 3 };

        myCell1.append(myCell2);

        REQUIRE( myCell1(0) == 4 );  REQUIRE( myCell1(1) == 5 );
        REQUIRE( myCell1(2) == 9 );  REQUIRE( myCell1(3) == 8 );
        REQUIRE( myCell1(4) == 10 ); REQUIRE( myCell1(5) == 2 );
        REQUIRE( myCell1(6) == 3 );
    }

    
    SECTION( "Vector print_as_row_vector")
    {
        moris::Vector< moris::uint > myCell1 = { 4, 5, 9, 8 };
        moris::Vector< moris::Vector<moris::uint> > myCell2 = { {10, 2, 3 }, { 4, 5, 9, 8 } };
        moris::Vector< moris::Vector<moris::Vector<moris::uint> > > myCell3 = {{{1,2,3}, {4,5,6}}, {{7,8,9}}};;

        print_as_row_vector(myCell1);
        print_as_row_vector(myCell2);
        print_as_row_vector(myCell3);

    }

}

