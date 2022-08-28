/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find_unique.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_find_unique.hpp"

namespace moris
{
TEST_CASE( "moris::find_unique", "[linalgebra],[find_unique]" )
    {
    SECTION( "find unique values of uint row vector" )
    {
        Matrix< DDUMat > A( 6, 1 );

        A( 0, 0 ) = 2;
        A( 1, 0 ) = 2;
        A( 2, 0 ) = 4;
        A( 3, 0 ) = 4;
        A( 4, 0 ) = 6;
        A( 5, 0 ) = 6;

        Matrix< DDNIMat > C = find_unique(A);

        REQUIRE( C(0,0) == 0 );
        REQUIRE( C(1,0) == 2 );
        REQUIRE( C(2,0) == 4 );
    }
    SECTION( "find unique values of uint col vector" )
    {
        Matrix< DDUMat > A( 1, 6 );

        A( 0, 0 ) = 2;
        A( 0, 1 ) = 2;
        A( 0, 2 ) = 4;
        A( 0, 3 ) = 4;
        A( 0, 4 ) = 6;
        A( 0, 5 ) = 6;

        Matrix< DDNIMat > C = find_unique(A);

        REQUIRE( C(0,0) == 0 );
        REQUIRE( C(1,0) == 2 );
        REQUIRE( C(2,0) == 4 );
    }

    SECTION( "find unique values of real row vector" )
    {
        Matrix< DDRMat > A( 6, 1 );

        A( 0, 0 ) = 2.2;
        A( 1, 0 ) = 2.2;
        A( 2, 0 ) = 4.5;
        A( 3, 0 ) = 4.5;
        A( 4, 0 ) = 6.3;
        A( 5, 0 ) = 6.3;

        Matrix< DDNIMat > C = find_unique(A);

        REQUIRE( C(0,0) == 0 );
        REQUIRE( C(1,0) == 2 );
        REQUIRE( C(2,0) == 4 );
    }

    SECTION( "find unique values of real col vector" )
    {
        Matrix< DDRMat > A( 1, 6 );

        A( 0, 0 ) = 2.2;
        A( 0, 1 ) = 2.2;
        A( 0, 2 ) = 4.5;
        A( 0, 3 ) = 4.5;
        A( 0, 4 ) = 6.3;
        A( 0, 5 ) = 6.3;

        Matrix< DDNIMat > C = find_unique(A);

        REQUIRE( C(0,0) == 0 );
        REQUIRE( C(1,0) == 2 );
        REQUIRE( C(2,0) == 4 );
    }
    }
}

