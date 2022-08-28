/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_trans.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"

TEST_CASE("moris::trans",
          "[linalgebra],[trans]")
{
    moris::Matrix< moris::DDRMat > A( 3, 3 );

    A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
    A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
    A( 2, 0 ) = 7.0; A( 2, 1 ) = 8.0; A( 2, 2 ) = 9.0;

    moris::Matrix< moris::DDRMat > B( 3, 3 );
    B = moris::trans( A );

    REQUIRE( B( 0,0 ) == 1.0 ); REQUIRE( B( 0,1 ) == 4.0 ); REQUIRE( B( 0,2 ) == 7.0 );
    REQUIRE( B( 1,0 ) == 2.0 ); REQUIRE( B( 1,1 ) == 5.0 ); REQUIRE( B( 1,2 ) == 8.0 );
    REQUIRE( B( 2,0 ) == 3.0 ); REQUIRE( B( 2,1 ) == 6.0 ); REQUIRE( B( 2,2 ) == 9.0 );
}

