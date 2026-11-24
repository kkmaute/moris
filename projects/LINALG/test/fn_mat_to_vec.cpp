/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_vectorize.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_mat_to_vec.hpp"

namespace moris
{
    TEST_CASE(
            "mat_to_vec",
            "[linalgebra],[mat_to_vec]" )
    {
        // create rectangular matrix of size 4x2
        Matrix< DDRMat > A( 4, 2 );
        A( 0, 0 ) = 1.1;
        A( 1, 0 ) = 2.1;
        A( 2, 0 ) = 3.1;
        A( 3, 0 ) = 4.1;
        A( 0, 1 ) = 1.2;
        A( 1, 1 ) = 2.2;
        A( 2, 1 ) = 3.2;
        A( 3, 1 ) = 4.2;

        // flatten matrix into column vector
        moris::Vector< real > B = mat_to_vec( A, false );

        REQUIRE( equal_to( B( 0 ), 1.1 ) );
        REQUIRE( equal_to( B( 1 ), 2.1 ) );
        REQUIRE( equal_to( B( 2 ), 3.1 ) );
        REQUIRE( equal_to( B( 3 ), 4.1 ) );
        REQUIRE( equal_to( B( 4 ), 1.2 ) );
        REQUIRE( equal_to( B( 5 ), 2.2 ) );
        REQUIRE( equal_to( B( 6 ), 3.2 ) );
        REQUIRE( equal_to( B( 7 ), 4.2 ) );
    }
}    // namespace moris
