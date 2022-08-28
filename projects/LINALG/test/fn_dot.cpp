/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_dot.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_dot.hpp"

//calculates the dot product of two vectors

TEST_CASE("moris::dot",
             "[linalgebra],[dot]" )
{
    moris::Matrix< moris::DDRMat > A( 3, 1 );

    A( 0, 0 ) = 1.0;
    A( 1, 0 ) = 3.0;
    A( 2, 0 ) = -5.0;

    moris::Matrix< moris::DDRMat > B( 3, 1 );

    B(0, 0) = 4;
    B(1, 0) = -2;
    B(2, 0) = -1;

    moris::real C = dot(A, B);

    REQUIRE( moris::equal_to( C, 3 ) );

}

