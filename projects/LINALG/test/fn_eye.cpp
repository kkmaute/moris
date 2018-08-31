/*
 * fn_eye.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_eye.hpp"

namespace moris
{
TEST_CASE(
        "moris::eye",
        "[linalgebra],[eye]" )
{
    Matrix< real, DDRMat > a;
    eye( 3, 3, a );

    REQUIRE( a( 0,0 ) == 1.0 );
    REQUIRE( a( 0,1 ) == 0.0 );
    REQUIRE( a( 0,2 ) == 0.0 );

    REQUIRE( a( 1,0 ) == 0.0 );
    REQUIRE( a( 1,1 ) == 1.0 );
    REQUIRE( a( 1,2 ) == 0.0 );

    REQUIRE( a( 2,0 ) == 0.0 );
    REQUIRE( a( 2,1 ) == 0.0 );
    REQUIRE( a( 2,2 ) == 1.0 );
}
}
