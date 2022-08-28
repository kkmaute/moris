/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_div.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_div.hpp"

TEST_CASE(
         "moris::op_div",
         "[linalgebra],[op_div]" )
{
    SECTION( "scalar division" )
    {
        moris::Matrix< moris::DDRMat > tValues = { { -1 }, { 0 }, { 1 }, { 2 } };

        moris::Matrix< moris::DDRMat > tDivide = tValues / 2.0;

        REQUIRE( tDivide( 0, 0 ) == -0.5 );
        REQUIRE( tDivide( 1, 0 ) ==  0.0 );
        REQUIRE( tDivide( 2, 0 ) ==  0.5 );
        REQUIRE( tDivide( 3, 0 ) ==  1.0 );
    }
}

