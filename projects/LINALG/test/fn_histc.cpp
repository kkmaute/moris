/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_histc.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_histc.hpp"
#include "fn_unique.hpp"

namespace moris
{
TEST_CASE( "moris::histc", "[linalgebra],[histc]" )
    {
    //Demonstrates the functionality of "histc". It provides a vector, that contains the counts of the number of values that fall between the elements in a given range.
    // histc first sorts the vector and takes then the counts of the number of values.
    // Only a column vector as input is allowed
    SECTION( "histc of uint col vector" )
    {
    Matrix< DDUMat > A( 7, 1 );
    Matrix< DDUMat > B;

    A( 0, 0 ) = 1;
    A( 1, 0 ) = 2;
    A( 2, 0 ) = 3;
    A( 3, 0 ) = 2;
    A( 4, 0 ) = 2;
    A( 5, 0 ) = 4;
    A( 6, 0 ) = 1;

    unique( A, B );

    Matrix< DDNIMat > C = histc( A, B );
    REQUIRE( C(0,0) == 2 );
    REQUIRE( C(1,0) == 3 );
    REQUIRE( C(2,0) == 1 );
    REQUIRE( C(3,0) == 1 );
    }
}
}

