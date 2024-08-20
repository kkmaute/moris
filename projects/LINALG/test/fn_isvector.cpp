/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isvector.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_isvector.hpp"

namespace moris
{
    TEST_CASE( "moris::isvector", "[linalgebra],[isvector]" )
    {
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b( 1, 3 );
        Matrix< DDRMat > c( 1, 1 );
        Matrix< DDRMat > d;

        bool tIsVector_1 = isvector( a );
        bool tIsVector_2 = isvector( b );
        bool tIsVector_3 = isvector( b );
        bool tIsVector_4 = isvector( d );

        REQUIRE( !tIsVector_1 );
        REQUIRE( tIsVector_2 );
        REQUIRE( tIsVector_3 );
        REQUIRE( !tIsVector_4 );
    }
}    // namespace moris
