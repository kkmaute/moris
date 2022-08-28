/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_unique.cpp
 *
 */

#include <catch.hpp>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_unique.hpp"

namespace moris
{
TEST_CASE("unique",
          "[linalgebra],[unique]")
{
    //Demonstrates the functionality of "unique". It provides a vector, in which every value from a given vector is only once in there.
    // Unique first sorts the vector and takes then the unique entries.
    SECTION( "unique of uint row vector" )
    {
        //Uniqueness of row vector

        Matrix< DDUMat > A( 7, 1 );

        A( 0 ) = 1;
        A( 1 ) = 2;
        A( 2 ) = 3;
        A( 3 ) = 2;
        A( 4 ) = 2;
        A( 5 ) = 4;
        A( 6 ) = 1;

        Matrix< DDUMat > C;
        unique( A, C );
        REQUIRE( C( 0 ) == 1 );
        REQUIRE( C( 1 ) == 2 );
        REQUIRE( C( 2 ) == 3 );
    }

    SECTION( "unique of real row vector" )
    {

        //Uniqueness of row vector
        Matrix< DDRMat >A( 7, 1 );

        A( 0 ) = 1.1;
        A( 1 ) = 2.1;
        A( 2 ) = 2.1;
        A( 3 ) = 3.1;
        A( 4 ) = 2.1;
        A( 5 ) = 4.1;
        A( 6 ) = 1.1;

        Matrix< DDRMat > C;
        unique( A, C );

        REQUIRE( C( 0 ) == 1.1 );
        REQUIRE( C( 1 ) == 2.1 );
        REQUIRE( C( 2 ) == 3.1 );
    }

    SECTION( "unique of uint col vector" )
    {
        //Uniqueness of col vector

        Matrix< DDUMat > A( 1,7 );

        A( 0 ) = 1;
        A( 1 ) = 2;
        A( 2 ) = 2;
        A( 3 ) = 3;
        A( 4 ) = 1;
        A( 5 ) = 2;
        A( 6 ) = 4;

        Matrix< DDUMat > C;
        unique( A, C );

        REQUIRE( C( 0 ) == 1 );
        REQUIRE( C( 1 ) == 2 );
        REQUIRE( C( 2 ) == 3 );
        REQUIRE( C( 3 ) == 4 );
    }

    SECTION( "unique of real col vector" )
    {
        //Uniqueness of col vector

        Matrix< DDRMat > A( 1,8 );

        A( 0 ) = 1.1;
        A( 1 ) = 2.1;
        A( 2 ) = 2.1;
        A( 3 ) = 3.1;
        A( 4 ) = 1.1;
        A( 5 ) = 4.1;
        A( 6 ) = 2.3;
        A( 7 ) = 5.3;

        Matrix< DDRMat > C;
        unique( A, C );

        REQUIRE( C( 0 ) == 1.1 );
        REQUIRE( C( 1 ) == 2.1 );
        REQUIRE( C( 2 ) == 2.3 );
        REQUIRE( C( 3 ) == 3.1 );
        REQUIRE( C( 4 ) == 4.1 );
    }
}
}

