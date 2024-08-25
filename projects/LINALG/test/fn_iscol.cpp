/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_iscol.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_iscol.hpp"

namespace moris
{
    TEST_CASE( "moris::iscol", "[linalgebra],[iscol]" )
    {
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b( 3, 1 );
        Matrix< DDRMat > c( 1, 1 );
        Matrix< DDRMat > d;

        bool tIsCol_1 = iscol( a );
        bool tIsCol_2 = iscol( b );
        bool tIsCol_3 = iscol( c );
        bool tIsCol_4 = iscol( d );

        REQUIRE( !tIsCol_1 );
        REQUIRE( tIsCol_2 );
        REQUIRE( tIsCol_3 );
        REQUIRE( !tIsCol_4 );
    }
}    // namespace moris
