// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::op_minus",
        "[linalgebra],[op_minus]" )
{

#include "linalg/op_minus.inc"

    SECTION( "moris::Mat - moris::Mat - moris::Mat" )
    {
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

    SECTION( "moris::Sp_Mat - moris::Sp_Mat - moris::Sp_Mat" )
    {
       REQUIRE( Ds( 0, 0 ) == -4.0 );
       REQUIRE( Ds( 0, 1 ) == -4.0 );
       REQUIRE( Ds( 0, 2 ) == -4.0 );

       REQUIRE( Ds( 1, 0 ) == -4.0 );
       REQUIRE( Ds( 1, 1 ) == -4.0 );
       REQUIRE( Ds( 1, 2 ) == -4.0 );

       REQUIRE( Ds( 2, 0 ) == -4.0 );
       REQUIRE( Ds( 2, 1 ) == -4.0 );
       REQUIRE( Ds( 2, 2 ) == -4.0 );
    }

    SECTION( "moris::Col - moris::Col - moris::Col" )
    {
        REQUIRE( Dc( 0, 0 ) == -4.0 );
        REQUIRE( Dc( 1, 0 ) == -4.0 );
        REQUIRE( Dc( 2, 0 ) == -4.0 );
    }
}
