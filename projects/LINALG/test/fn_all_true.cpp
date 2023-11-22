/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_all_true.cpp
 *
 */

#include <catch.hpp>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_print.hpp"

namespace moris
{
    TEST_CASE( "moris::all_true",
            "[linalgebra],[all_true]" )
    {
        Matrix< DDRMat > A( 3, 1 );
        Matrix< DDRMat > B( 3, 1 );

        A( 0, 0 ) = 1.0;
        A( 1, 0 ) = 2.0;
        A( 2, 0 ) = 3.0;

        B( 0, 0 ) = 1.0;
        B( 1, 0 ) = 0.0;
        B( 2, 0 ) = 3.0;

        REQUIRE( !all_true( A == B ) );

        B( 1, 0 ) = 2.0;

        REQUIRE( all_true( A == B ) );
    }
}    // namespace moris
