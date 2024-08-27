/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isempty.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_isempty.hpp"

namespace moris
{
    TEST_CASE( "moris::is_empty", "[linalgebra],[is_empty]" )
    {
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b;
        Matrix< DDRMat > c( 0, 3 );

        bool tIsEmpty_1 = isempty( a );
        bool tIsEmpty_2 = isempty( b );
        bool tIsEmpty_3 = isempty( b );

        REQUIRE( !tIsEmpty_1 );
        REQUIRE( tIsEmpty_2 );
        REQUIRE( tIsEmpty_3 );
    }
}    // namespace moris
