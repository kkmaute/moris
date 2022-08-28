/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sum.cpp
 *
 */

#include <catch.hpp>

// MORIS project header files.

#include "cl_Matrix.hpp"
#include "fn_equal_to.hpp" // ALG/src
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_sum.hpp"

TEST_CASE(
        "moris::sum",
        "[linalgebra],[sum]" )
{
    //Demonstrates the functionality of "sum". You can sum a col or row vector or a matrix.
    SECTION( "Sum of row vector" )
    {
    //Sum of row vector

    moris::Matrix< moris::DDRMat > A( 3, 1 );

    A( 0, 0 ) = 1.0;
    A( 1, 0 ) = 2.0;
    A( 2, 0 ) = 3.0;

    moris::real C = sum(A);

    REQUIRE( moris::equal_to( C, 6 ) );
    }

    SECTION( "Sum of col vector" )
    {
    //Sum of col vector
    moris::Matrix< moris::DDRMat > B( 1, 3 );

    B( 0, 0 ) = 1.0;
    B( 0, 1 ) = 2.0;
    B( 0, 2 ) = 3.0;

    moris::real C = sum(B);

    REQUIRE( moris::equal_to( C, 6 ) );
    }

    SECTION( "Sum of a matrix" )
    {
    //Sum of col vector
    moris::Matrix< moris::DDRMat > B( 3, 3 );

    B( 0, 0 ) = 1.0;    B( 0, 1 ) = 1.0;     B( 0, 2 ) = 1.0;
    B( 1, 0 ) = 2.0;    B( 1, 1 ) = 2.0;     B( 1, 2 ) = 2.0;
    B( 2, 0 ) = 3.0;    B( 2, 1 ) = 3.0;     B( 2, 2 ) = 3.0;

    moris::real C = sum(B);

    REQUIRE( moris::equal_to( C, 18 ) );
    }

}

