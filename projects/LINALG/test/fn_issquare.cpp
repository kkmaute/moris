/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_issquare.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_issquare.hpp"

namespace moris
{
    TEST_CASE( "moris::issquare", "[linalgebra],[issquare]" )
    {
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b( 1, 3 );
        Matrix< DDRMat > c( 1, 1 );
        Matrix< DDRMat > d;

        bool tIsSquare_1 = issquare( a );
        bool tIsSquare_2 = issquare( b );
        bool tIsSquare_3 = issquare( c );
        bool tIsSquare_4 = issquare( d );

        REQUIRE( tIsSquare_1 );
        REQUIRE( !tIsSquare_2 );
        REQUIRE( tIsSquare_3 );
        REQUIRE( tIsSquare_4 );
    }
}    // namespace moris
