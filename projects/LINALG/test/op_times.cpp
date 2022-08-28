/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_times.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_times.hpp"

TEST_CASE(
         "moris::op_times",
         "[linalgebra],[op_times]" )
{

     moris::Matrix< moris::DDRMat > Am( 3, 3 );
     moris::Matrix< moris::DDRMat > Bm( 3, 3 );
     moris::Matrix< moris::DDRMat > Cm( 3, 3 );

    Am( 0, 0 ) = 1.0; Am( 0, 1 ) = 1.0; Am( 0, 2 ) = 1.0;
    Am( 1, 0 ) = 1.0; Am( 1, 1 ) = 1.0; Am( 1, 2 ) = 1.0;
    Am( 2, 0 ) = 1.0; Am( 2, 1 ) = 1.0; Am( 2, 2 ) = 1.0;

    Bm( 0, 0 ) = 2.0; Bm( 0, 1 ) = 2.0; Bm( 0, 2 ) = 2.0;
    Bm( 1, 0 ) = 2.0; Bm( 1, 1 ) = 2.0; Bm( 1, 2 ) = 2.0;
    Bm( 2, 0 ) = 2.0; Bm( 2, 1 ) = 2.0; Bm( 2, 2 ) = 2.0;

    Cm = Am * Bm;

    SECTION( "moris::Mat * moris::Mat" )
    {
        REQUIRE( Cm( 0, 0 ) == 6.0 );
        REQUIRE( Cm( 0, 1 ) == 6.0 );
        REQUIRE( Cm( 0, 2 ) == 6.0 );

        REQUIRE( Cm( 1, 0 ) == 6.0 );
        REQUIRE( Cm( 1, 1 ) == 6.0 );
        REQUIRE( Cm( 1, 2 ) == 6.0 );

        REQUIRE( Cm( 2, 0 ) == 6.0 );
        REQUIRE( Cm( 2, 1 ) == 6.0 );
        REQUIRE( Cm( 2, 2 ) == 6.0 );
    }
}

