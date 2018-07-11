// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::choleskydecomposition_lower",
        "[linalgebra],[choleskydecomposition_lower]" )
{
    SECTION( "choll( moris::Mat, moris::Mat )" )
    {
        #include "linalg/fn_chol_l.inc"

        REQUIRE( moris::equal_to( L( 0,0 ), +1.414213562373095e+00 ) );
        REQUIRE( moris::equal_to( L( 0,1 ), +0.000000000000000e+00 ) );
        REQUIRE( moris::equal_to( L( 0,2 ), +0.000000000000000e+00 ) );

        REQUIRE( moris::equal_to( L( 1,0 ), -7.071067811865475e-01 ) );
        REQUIRE( moris::equal_to( L( 1,1 ), +1.224744871391589e+00 ) );
        REQUIRE( moris::equal_to( L( 1,2 ), +0.000000000000000e+00 ) );

        REQUIRE( moris::equal_to( L( 2,0 ), +0.000000000000000e+00 ) );
        REQUIRE( moris::equal_to( L( 2,1 ), -8.164965809277261e-01 ) );
        REQUIRE( moris::equal_to( L( 2,2 ), +1.154700538379251e+00 ) );
    }
}
