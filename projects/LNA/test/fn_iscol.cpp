// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::iscol",
        "[linalgebra],[iscol]" )
{
    SECTION( "iscol( moris::Mat )" )
    {
        #include "linalg/fn_iscol.inc"
        REQUIRE( e == false );
        REQUIRE( f == true  );
        REQUIRE( g == false );
        REQUIRE( h == true  );

    }

    SECTION( "iscol( moris::Mat cplx )" )
    {
        moris::Mat< moris::cplx > A( 3, 3 );
        moris::Mat< moris::cplx > B( 3, 1 );
        moris::Mat< moris::cplx > C( 1, 3 );
        moris::Mat< moris::cplx > D( 1, 1 );

        A( 0, 0 ) = { 1.0, 3.0 };
        A( 0, 1 ) = { 2.0, 2.0 };
        A( 0, 2 ) = { 3.0, 2.5 };

        A( 1, 0 ) = { 4.0, 5.0 };
        A( 1, 1 ) = { 5.0, 3.0 };
        A( 1, 2 ) = { 6.0, 4.5 };

        A( 2, 0 ) = { 7.0, 3.0 };
        A( 2, 1 ) = { 8.0, 1.5 };
        A( 2, 2 ) = { 9.0, 0.0 };

        B( 0, 0 ) = { 1.0, 2.5 };
        B( 1, 0 ) = { 2.0, 3.5 };
        B( 2, 0 ) = { 3.0, 3.4 };

        C( 0, 0 ) = { 1.0, 3.6 };
        C( 0, 1 ) = { 2.0, 1.7 };
        C( 0, 2 ) = { 3.0, 1.3 };

        D( 0, 0 ) = {1.0, 0.0 };

        bool e = moris::iscol( A );
        bool f = moris::iscol( B );
        bool g = moris::iscol( C );
        bool h = moris::iscol( D );

        REQUIRE( e == false );
        REQUIRE( f == true  );
        REQUIRE( g == false );
        REQUIRE( h == true  );

    }
}
