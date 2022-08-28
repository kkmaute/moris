/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_elemwise_div.cpp
 *
 */

#include <catch.hpp>
#include <op_elemwise_mult.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_elemwise_div.hpp"
#include "op_plus.hpp"
#include "fn_print.hpp"

namespace moris
{
    TEST_CASE(
            "moris::elemwise_div",
            "[linalgebra],[elemwise_div]" )
    {
        moris::Matrix< moris::DDRMat > tA = { { 4 }, { 64 }, { 162 }, { 196 } };
        moris::Matrix< moris::DDRMat > tB = { { 2 }, { 8 }, { 18 }, { 14 } };

        moris::Matrix< moris::DDRMat > tC = tA / tB;

        REQUIRE( tC( 0 ) == 2 );
        REQUIRE( tC( 1 ) == 8 );
        REQUIRE( tC( 2 ) == 9 );
        REQUIRE( tC( 3 ) == 14 );

        moris::Matrix< moris::DDRMat > tD = tA / ( tB + tC );

        REQUIRE( tD( 0 ) == 1 );
        REQUIRE( tD( 1 ) == 4 );
        REQUIRE( tD( 2 ) == 6 );
        REQUIRE( tD( 3 ) == 7 );

        moris::Matrix< moris::DDRMat > tE = ( tA + tB ) / ( tC );

        REQUIRE( tE( 0 ) == 3 );
        REQUIRE( tE( 1 ) == 9 );
        REQUIRE( tE( 2 ) == 20 );
        REQUIRE( tE( 3 ) == 15 );

        moris::Matrix< moris::DDRMat > tF = ( tA + tB ) / ( tC + tD );

        REQUIRE( tF( 0 ) == 2 );
        REQUIRE( tF( 1 ) == 6 );
        REQUIRE( tF( 2 ) == 12 );
        REQUIRE( tF( 3 ) == 10 );
    }
}

