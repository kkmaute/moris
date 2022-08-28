/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_minus.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"

namespace moris
{
TEST_CASE(
        "moris::op_minus",
        "[linalgebra],[op_minus]" )
{
    Matrix< DDRMat > Am( 3, 3 );
    Matrix< DDRMat > Bm( 3, 3 );
    Matrix< DDRMat > Cm( 3, 3 );
    Matrix< DDRMat > Dm( 3, 3 );

    Am( 0, 0 ) = 1.0; Am( 0, 1 ) = 1.0; Am( 0, 2 ) = 1.0;
    Am( 1, 0 ) = 1.0; Am( 1, 1 ) = 1.0; Am( 1, 2 ) = 1.0;
    Am( 2, 0 ) = 1.0; Am( 2, 1 ) = 1.0; Am( 2, 2 ) = 1.0;

    Bm( 0, 0 ) = 2.0; Bm( 0, 1 ) = 2.0; Bm( 0, 2 ) = 2.0;
    Bm( 1, 0 ) = 2.0; Bm( 1, 1 ) = 2.0; Bm( 1, 2 ) = 2.0;
    Bm( 2, 0 ) = 2.0; Bm( 2, 1 ) = 2.0; Bm( 2, 2 ) = 2.0;

    Cm( 0, 0 ) = 3.0; Cm( 0, 1 ) = 3.0; Cm( 0, 2 ) = 3.0;
    Cm( 1, 0 ) = 3.0; Cm( 1, 1 ) = 3.0; Cm( 1, 2 ) = 3.0;
    Cm( 2, 0 ) = 3.0; Cm( 2, 1 ) = 3.0; Cm( 2, 2 ) = 3.0;

    Dm = Am - Bm - Cm;
    REQUIRE( Dm( 0, 0 ) == -4.0 );
    REQUIRE( Dm( 0, 1 ) == -4.0 );
    REQUIRE( Dm( 0, 2 ) == -4.0 );

    REQUIRE( Dm( 1, 0 ) == -4.0 );
    REQUIRE( Dm( 1, 1 ) == -4.0 );
    REQUIRE( Dm( 1, 2 ) == -4.0 );

    REQUIRE( Dm( 2, 0 ) == -4.0 );
    REQUIRE( Dm( 2, 1 ) == -4.0 );
    REQUIRE( Dm( 2, 2 ) == -4.0 );
}
}

