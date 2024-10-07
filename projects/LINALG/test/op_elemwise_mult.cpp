/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_elemwise_mult.cpp
 *
 */

#include <catch.hpp>
#include <op_elemwise_mult.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_plus.hpp"

namespace moris
{
    TEST_CASE(
            "moris::elemwise_mult",
            "[linalgebra],[elemwise_mult]" )
    {
        moris::Matrix< moris::DDRMat > tA = { { 1 }, { 2 }, { 3 }, { 5 } };
        moris::Matrix< moris::DDRMat > tB = { { 7 }, { 11 }, { 13 }, { 17 } };

        moris::Matrix< moris::DDRMat > tC = tA % tB;

        REQUIRE( tC( 0 ) == 7 );
        REQUIRE( tC( 1 ) == 22 );
        REQUIRE( tC( 2 ) == 39 );
        REQUIRE( tC( 3 ) == 85 );

        moris::Matrix< moris::DDRMat > tD = tA % ( tB + tC );

        REQUIRE( tD( 0 ) == 14 );
        REQUIRE( tD( 1 ) == 66 );
        REQUIRE( tD( 2 ) == 156 );
        REQUIRE( tD( 3 ) == 510 );

        moris::Matrix< moris::DDRMat > aE = ( tA + tB ) % ( tC );

        REQUIRE( aE( 0 ) == 56 );
        REQUIRE( aE( 1 ) == 286 );
        REQUIRE( aE( 2 ) == 624 );
        REQUIRE( aE( 3 ) == 1870 );

        moris::Matrix< moris::DDRMat > tF = ( tA + tB ) % ( tC + tD );

        REQUIRE( tF( 0 ) == 168 );
        REQUIRE( tF( 1 ) == 1144 );
        REQUIRE( tF( 2 ) == 3120 );
        REQUIRE( tF( 3 ) == 13090 );
    }
}

