/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_chol_u.cpp
 *
 */

#include <catch.hpp>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_chol_u.hpp"

TEST_CASE("moris::chol_u",
          "[linalgebra],[chol_u]")
{
    moris::Matrix< moris::DDRMat >  tA = { { 4, 12, -16 } ,  { 12, 37, -43 } , { -16, -43, 98 } };

    moris::Matrix< moris::DDRMat > tU = chol_u( tA );

    REQUIRE( tU( 0, 0 ) ==  2 );
    REQUIRE( tU( 1, 0 ) ==  0 );
    REQUIRE( tU( 2, 0 ) ==  0 );
    REQUIRE( tU( 0, 1 ) ==  6 );
    REQUIRE( tU( 1, 1 ) ==  1 );
    REQUIRE( tU( 2, 1 ) ==  0 );
    REQUIRE( tU( 0, 2 ) == -8 );
    REQUIRE( tU( 1, 2 ) ==  5 );
    REQUIRE( tU( 2, 2 ) ==  3 );

}

