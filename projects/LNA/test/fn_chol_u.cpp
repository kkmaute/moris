// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::choleskydecompositionupper",
        "[linalgebra],[choleskydecompositionupper]" )
{
    SECTION( "cholu( moris::Mat, moris::Mat )" )
    {
        #include "linalg/fn_chol_u.inc"
        REQUIRE( U( 0,0 ) - 1.414213562373095e+00 < 1.0e-12 );
        REQUIRE( U( 0,1 ) + 7.071067811865475e-01 < 1.0e-12 );
        REQUIRE( U( 0,2 ) - 0.000000000000000e+00 < 1.0e-12 );

        REQUIRE( U( 1,0 ) - 0.000000000000000e+00 < 1.0e-12 );
        REQUIRE( U( 1,1 ) - 1.224744871391589e+00 < 1.0e-12 );
        REQUIRE( U( 1,2 ) + 8.164965809277261e-01 < 1.0e-12 );

        REQUIRE( U( 2,0 ) - 0.000000000000000e+00 < 1.0e-12 );
        REQUIRE( U( 2,1 ) - 0.000000000000000e+00 < 1.0e-12 );
        REQUIRE( U( 2,2 ) - 1.154700538379251e+00 < 1.0e-12 );
    }
}
