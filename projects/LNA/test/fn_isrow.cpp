// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::isrow",
        "[linalgebra],[isrow]" )
{
    SECTION( "isrow( moris::Mat )" )
    {
        #include "linalg/fn_isrow.inc"
        REQUIRE( e == false );
        REQUIRE( f == false );
        REQUIRE( g == true  );
        REQUIRE( h == true  );

    }

    SECTION( "isrow( moris::Mat cplx )" )
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

        bool e = moris::isrow( A );
        bool f = moris::isrow( B );
        bool g = moris::isrow( C );
        bool h = moris::isrow( D );

        REQUIRE( e == false );
        REQUIRE( f == false );
        REQUIRE( g == true  );
        REQUIRE( h == true  );

    }
}
