/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linspace.cpp
 *
 */

#include <catch.hpp>
#include <cmath>

// MORIS project header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_equal_to.hpp"
#include "fn_linspace.hpp"

#include "fn_print.hpp"
// ----------------------------------------------------------------------------
namespace moris
{
TEST_CASE("moris::linspace", "[linalgebra],[linspace]" )
{

    SECTION( "linspace" )
    {
        moris::Matrix< moris::DDRMat > aVecNr  = moris::linspace( 0.0, 1.0, 6 );
        Matrix< DDSMat > aVecNsi = moris::linspace( 0, 10, 6 );
        Matrix< DDLMat > aVecNli = moris::linspace( 0, 10, 6 );

        SECTION( "moris::linspace real" )
        {
            REQUIRE( moris::equal_to( aVecNr( 0, 0 ), 0.00 ) );
            REQUIRE( moris::equal_to( aVecNr( 1, 0 ), 0.20 ) );
            REQUIRE( moris::equal_to( aVecNr( 2, 0 ), 0.40 ) );
            REQUIRE( moris::equal_to( aVecNr( 3, 0 ), 0.60 ) );
            REQUIRE( moris::equal_to( aVecNr( 4, 0 ), 0.80 ) );
            REQUIRE( moris::equal_to( aVecNr( 5, 0 ), 1.00 ) );
        }

        SECTION( "moris::linspace sint" )
        {
            REQUIRE( moris::equal_to( aVecNsi( 0, 0 ), 0  ) );
            REQUIRE( moris::equal_to( aVecNsi( 1, 0 ), 2  ) );
            REQUIRE( moris::equal_to( aVecNsi( 2, 0 ), 4  ) );
            REQUIRE( moris::equal_to( aVecNsi( 3, 0 ), 6  ) );
            REQUIRE( moris::equal_to( aVecNsi( 4, 0 ), 8  ) );
            REQUIRE( moris::equal_to( aVecNsi( 5, 0 ), 10 ) );
        }

        SECTION( "moris::linspace lint" )
        {
            REQUIRE( moris::equal_to( aVecNli( 0, 0 ), 0  ) );
            REQUIRE( moris::equal_to( aVecNli( 1, 0 ), 2  ) );
            REQUIRE( moris::equal_to( aVecNli( 2, 0 ), 4  ) );
            REQUIRE( moris::equal_to( aVecNli( 3, 0 ), 6  ) );
            REQUIRE( moris::equal_to( aVecNli( 4, 0 ), 8  ) );
            REQUIRE( moris::equal_to( aVecNli( 5, 0 ), 10 ) );
        }

    }
}
}

