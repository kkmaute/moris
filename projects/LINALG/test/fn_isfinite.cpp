/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isfinite.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_isfinite.hpp"

namespace moris
{
    TEST_CASE( "moris::isfinite", "[linalgebra],[isfinite]" )
    {
        moris::real      tZero = 0.0;
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b( 3, 3 );

        a( 0, 0 ) = 1.0;
        a( 0, 1 ) = 2.0;
        a( 0, 2 ) = 3.0;
        a( 1, 0 ) = 0.0;
        a( 1, 1 ) = 1.0;
        a( 1, 2 ) = 4.0;
        a( 2, 0 ) = 5.0;
        a( 2, 1 ) = 6.0;
        a( 2, 2 ) = 0.0;

        b( 0, 0 ) = 1.0;
        b( 0, 1 ) = 2.0;
        b( 0, 2 ) = 3.0;
        b( 1, 0 ) = 0.0;
        b( 1, 1 ) = 1.0;
        b( 1, 2 ) = 4.0;
        b( 2, 0 ) = 5.0;
        b( 2, 1 ) = 6.0;
        b( 2, 2 ) = 1.0 / tZero;

        bool tIsFinite_1 = isfinite( a );
        bool tIsFinite_2 = isfinite( b );

        REQUIRE( tIsFinite_1 );
        REQUIRE( !tIsFinite_2 );
    }
}    // namespace moris
