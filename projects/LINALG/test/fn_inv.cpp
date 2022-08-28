/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_inv.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_inv.hpp"

TEST_CASE("moris::inv",
             "[linalgebra],[inv]" )
{
    // test of 3x3 matrix
    moris::Matrix< moris::DDRMat > A( 3, 3 );

    A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
    A( 1, 0 ) = 0.0; A( 1, 1 ) = 1.0; A( 1, 2 ) = 4.0;
    A( 2, 0 ) = 5.0; A( 2, 1 ) = 6.0; A( 2, 2 ) = 0.0;

    moris::Matrix< moris::DDRMat > B = moris::inv( A );

    REQUIRE( moris::equal_to( B( 0, 0 ), -24.0 ) );
    REQUIRE( moris::equal_to( B( 0, 1 ),  18.0 ) );
    REQUIRE( moris::equal_to( B( 0, 2 ),  5.0 ) );

    REQUIRE( moris::equal_to( B( 1, 0 ),  20.0 ) );
    REQUIRE( moris::equal_to( B( 1, 1 ), -15.0 ) );
    REQUIRE( moris::equal_to( B( 1, 2 ), -4.0 ) );

    REQUIRE( moris::equal_to( B( 2, 0 ),  -5.0 ) );
    REQUIRE( moris::equal_to( B( 2, 1 ),   4.0 ) );
    REQUIRE( moris::equal_to( B( 2, 2 ),  1.0 ) );

}

