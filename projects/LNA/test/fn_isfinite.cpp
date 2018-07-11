// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::isfinite",
        "[linalgebra],[isfinite]" )
{
// Testing the functionality of isfinite(A) to see if matrix is finite or not.

    SECTION( "moris::matrix" )
    {
        moris::Mat< moris::real > A( 3, 3 );
        A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
        A( 2, 0 ) = 7.0; A( 2, 1 ) = 8.0; A( 2, 2 ) = 9.0;

        bool e = moris::isfinite ( A );
        REQUIRE( e == true );

        moris::Mat< moris::real > B( 3, 3 );
        double inf = 1.0/0.0;
        B( 0, 0 ) = inf; B( 0, 1 ) = 2.0; B( 0, 2 ) = 3.0;
        B( 1, 0 ) = 4.0; B( 1, 1 ) = 5.0; B( 1, 2 ) = 6.0;
        B( 2, 0 ) = 7.0; B( 2, 1 ) = 8.0; B( 2, 2 ) = 9.0;

        bool f = moris::isfinite( B );
        REQUIRE_FALSE( f );
    }

    SECTION( "moris::Mat cplx" )
    {
        moris::Mat< moris::cplx > C( 3, 3 );
        C( 0, 0 ) = {0.0,1.0}; C( 0, 1 ) = {2.0,3.0}; C( 0, 2 ) = {4.0,5.0};
        C( 1, 0 ) = {6.0,7.0}; C( 1, 1 ) = {8.0,9.0}; C( 1, 2 ) = {0.0,1.0};
        C( 2, 0 ) = {2.0,3.0}; C( 2, 1 ) = {4.0,5.0}; C( 2, 2 ) = {6.0,7.0};

        bool g = moris::isfinite( C );
        REQUIRE( g == true );

        moris::Mat< moris::cplx > D( 3, 3 );
        double inf = 1.0/0.0;
        D( 0, 0 ) = {0.0,1.0}; D( 0, 1 ) = {2.0,3.0}; D( 0, 2 ) = {4.0,5.0};
        D( 1, 0 ) = {6.0,7.0}; D( 1, 1 ) = {inf,9.0}; D( 1, 2 ) = {0.0,1.0};
        D( 2, 0 ) = {2.0,3.0}; D( 2, 1 ) = {4.0,5.0}; D( 2, 2 ) = {6.0,7.0};

        bool h = moris::isfinite( D );
        REQUIRE_FALSE( h );
    }

    //Testing the functionality of isfinite(a) to see if vector is finite or not
    SECTION( "moris::matrix" )
    {
        moris::Mat< moris::real > v( 3, 1 );
        v( 0, 0 ) = 1.0; v( 1, 0 ) = 2.0; v( 2, 0 ) = 3.0;

        bool g = moris::isfinite( v );
        REQUIRE( g == true );

        moris::Mat< moris::real > w( 3, 1 );
        double nan = 0.0/0.0;
        w( 0, 0 ) = 1.0; w( 1, 0 ) = 2.0; w( 2, 0 ) = nan;

        bool h = moris::isfinite( w );
        REQUIRE_FALSE( h );
    }

    SECTION( "moris::Cell cplx" )
    {
        moris::Mat< moris::cplx > x( 3, 1 );
        x( 0, 0 ) = {0.0,1.0}; x( 1, 0 ) = {2.0,3.0}; x( 2, 0 ) = {4.0,5.0};

        bool k = moris::isfinite( x );
        REQUIRE( k == true );

        moris::Mat< moris::cplx > w( 3, 1 );
        double nan = 0.0/0.0;
        w( 0, 0 ) = {0.0,1.0}; w( 1, 0 ) = {nan,3.0}; w( 2, 0 ) = {4.0,5.0};

        bool l = moris::isfinite( w );
        REQUIRE_FALSE( l );
    }
}
