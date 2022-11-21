/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_unique.cpp
 *
 */

#include <catch.hpp>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_intersect.hpp"

namespace moris
{
    TEST_CASE( "intersect",
            "[linalgebra],[intersect]" )
    {
        SECTION( "intersect of uint row vectors" )
        {
            // Uniqueness of row vector
            Matrix< DDUMat > A( 5, 1 );
            Matrix< DDUMat > B( 3, 1 );

            A( 0 ) = 1;
            A( 1 ) = 2;
            A( 2 ) = 3;
            A( 3 ) = 4;
            A( 4 ) = 5;

            B( 0 ) = 5;
            B( 1 ) = 6;
            B( 2 ) = 7;

            Matrix< DDUMat > C = intersect( A, B );

            REQUIRE( C.numel() == 1 );
            REQUIRE( C( 0 ) == 5 );

            Matrix< DDUMat > Aind;
            Matrix< DDUMat > Bind;

            intersect( A, B, C, Aind, Bind );

            REQUIRE( C.numel() == 1 );
            REQUIRE( C( 0 ) == 5 );

            REQUIRE( Aind.numel() == 1 );
            REQUIRE( Bind.numel() == 1 );

            REQUIRE( Aind( 0 ) == 4 );
            REQUIRE( Bind( 0 ) == 0 );
        }

        SECTION( "intersect of sint col vectors" )
        {
            // Uniqueness of row vector
            Matrix< DDSMat > A( 1, 5 );
            Matrix< DDSMat > B( 1, 3 );

            A( 0 ) = 1;
            A( 1 ) = -2;
            A( 2 ) = 3;
            A( 3 ) = -4;
            A( 4 ) = 5;

            B( 0 ) = 5;
            B( 1 ) = 6;
            B( 2 ) = 7;

            Matrix< DDSMat > C = intersect( A, B );

            REQUIRE( C.numel() == 1 );
            REQUIRE( C( 0 ) == 5 );
        }

        SECTION( "intersect of sint col vectors with no column elements" )
        {
            // Uniqueness of row vector
            Matrix< DDUMat > A( 1, 5 );
            Matrix< DDUMat > B( 1, 3 );

            A( 0 ) = 1;
            A( 1 ) = -2;
            A( 2 ) = 3;
            A( 3 ) = -4;
            A( 4 ) = 5;

            B( 0 ) = -5;
            B( 1 ) = 6;
            B( 2 ) = 7;

            Matrix< DDUMat > C = intersect( A, B );

            REQUIRE( C.numel() == 0 );
        }
    }
}    // namespace moris
