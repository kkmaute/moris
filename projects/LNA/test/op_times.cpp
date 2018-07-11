// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
         "moris::op_times",
         "[linalgebra],[op_times]" )
{

#include "linalg/op_times.inc"

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

    SECTION( "moris::Sp_Mat * moris::Sp_Mat" )
    {
        REQUIRE( Cs( 0, 0 ) == 6.0 );
        REQUIRE( Cs( 0, 1 ) == 6.0 );
        REQUIRE( Cs( 0, 2 ) == 6.0 );

        REQUIRE( Cs( 1, 0 ) == 6.0 );
        REQUIRE( Cs( 1, 1 ) == 6.0 );
        REQUIRE( Cs( 1, 2 ) == 6.0 );

        REQUIRE( Cs( 2, 0 ) == 6.0 );
        REQUIRE( Cs( 2, 1 ) == 6.0 );
        REQUIRE( Cs( 2, 2 ) == 6.0 );
    }
}
