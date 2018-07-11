// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::eye",
        "[linalgebra],[eye]" )
{
#include "linalg/fn_eye.inc"

    REQUIRE( a( 0,0 ) == 1.0 );
    REQUIRE( a( 0,1 ) == 0.0 );
    REQUIRE( a( 0,2 ) == 0.0 );

    REQUIRE( a( 1,0 ) == 0.0 );
    REQUIRE( a( 1,1 ) == 1.0 );
    REQUIRE( a( 1,2 ) == 0.0 );

    REQUIRE( a( 2,0 ) == 0.0 );
    REQUIRE( a( 2,1 ) == 0.0 );
    REQUIRE( a( 2,2 ) == 1.0 );
}
