/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isrow.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_isrow.hpp"

namespace moris
{
    TEST_CASE( "moris::isrow", "[linalgebra],[isrow]" )
    {
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b( 1, 3 );
        Matrix< DDRMat > c( 1, 1 );
        Matrix< DDRMat > d;

        bool tIsRow_1 = isrow( a );
        bool tIsRow_2 = isrow( b );
        bool tIsRow_3 = isrow( b );
        bool tIsRow_4 = isrow( d );

        REQUIRE( !tIsRow_1 );
        REQUIRE( tIsRow_2 );
        REQUIRE( tIsRow_3 );
        REQUIRE( !tIsRow_4 );
    }
}    // namespace moris
