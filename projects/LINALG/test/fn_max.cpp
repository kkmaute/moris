/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_max.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_max.hpp"

namespace moris
{
TEST_CASE(
        "moris::max",
        "[linalgebra],[max]" )
{
    Matrix< DDRMat > tA{ {1.0,2.0,3.0}, {4.0,.0,6.0},{7.0,8.0,9.0}} ; 
    Matrix< DDRMat > tB{ {2.0,6.0,1.0}, {8.0,5.0,-2.0},{12.0,6.0,9.0}} ; 

    Matrix< DDRMat > tMax= max( tA, tB) ; 

    REQUIRE( tMax( 0,1 ) == 6.0 );
    REQUIRE( tMax( 0,0 ) == 2.0 );
    REQUIRE( tMax( 0,2 ) == 3.0 );

    REQUIRE( tMax( 1,0 ) == 8.0 );
    REQUIRE( tMax( 1,1 ) == 5.0 );
    REQUIRE( tMax( 1,2 ) == 6.0 );

    REQUIRE( tMax( 2,0 ) == 12.0 );
    REQUIRE( tMax( 2,1 ) == 8.0 );
    REQUIRE( tMax( 2,2 ) == 9.0 );
}
}