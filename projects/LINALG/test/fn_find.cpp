/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_find.hpp"

namespace moris
{
TEST_CASE( "moris::find", "[linalgebra],[find]" )
    {
    SECTION( "find of uint row vector" )
    {
        Matrix< DDUMat > A( 7, 1 );

        A( 0, 0 ) = 0;
        A( 1, 0 ) = 1;
        A( 2, 0 ) = 0;
        A( 3, 0 ) = 2;
        A( 4, 0 ) = 2;
        A( 5, 0 ) = 0;
        A( 6, 0 ) = 0;

        Matrix< DDNIMat > C = find(A);

        REQUIRE( C(0,0) == 1 );
        REQUIRE( C(1,0) == 3 );
        REQUIRE( C(2,0) == 4 );
    }
    SECTION( "find of real row vector" )
    {
        Matrix< DDRMat > A( 7, 1 );

        A( 0, 0 ) = 0;
        A( 1, 0 ) = 1.1;
        A( 2, 0 ) = 0;
        A( 3, 0 ) = 2.1;
        A( 4, 0 ) = 2.2;
        A( 5, 0 ) = 0;
        A( 6, 0 ) = 0;

        Matrix< DDNIMat > C = find(A);

        REQUIRE( C(0,0) == 1 );
        REQUIRE( C(1,0) == 3 );
        REQUIRE( C(2,0) == 4 );
    }
    SECTION( "find of uint col vector" )
    {
        Matrix< DDUMat > A( 1, 7 );

        A( 0, 0 ) = 0;
        A( 0, 1 ) = 1;
        A( 0, 2 ) = 0;
        A( 0, 3 ) = 2;
        A( 0, 4 ) = 2;
        A( 0, 5 ) = 0;
        A( 0, 6 ) = 0;

        Matrix< DDNIMat > C = find(A);

        REQUIRE( C(0,0) == 1 );
        REQUIRE( C(1,0) == 3 );
        REQUIRE( C(2,0) == 4 );
    }
    SECTION( "find first 2 nonzeros of uint row vector" )
    {
        Matrix< DDUMat > A( 7, 1 );

        A( 0, 0 ) = 0;
        A( 1, 0 ) = 1;
        A( 2, 0 ) = 0;
        A( 3, 0 ) = 2;
        A( 4, 0 ) = 2;
        A( 5, 0 ) = 0;
        A( 6, 0 ) = 1;

        Matrix< DDNIMat > C = find( A, 2 );

        REQUIRE( C(0,0) == 1 );
        REQUIRE( C(1,0) == 3 );
        REQUIRE( C.length() == 2 );
    }
    }
}

