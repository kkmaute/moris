// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::isempty",
        "[linalgebra],[isempty]" )
{
    //Testing the functionality of isempty(A) to see if a matrix is empty or not.

    SECTION( "moris::matrix" )
    {
        #include "linalg/fn_isempty.inc"
        REQUIRE_FALSE( e );

        moris::Mat< moris::real > B;
        bool f = moris::isempty( B );
        REQUIRE( f == true );
    }


    SECTION( "moris::Mat cplx" )
    {
        moris::Mat< moris::cplx > C( 3, 3 );

        C( 0, 0 ) = {0.0,1.0}; C( 0, 1 ) = {2.0,3.0}; C( 0, 2 ) = {4.0,5.0};
        C( 1, 0 ) = {6.0,7.0}; C( 1, 1 ) = {8.0,9.0}; C( 1, 2 ) = {0.0,1.0};
        C( 2, 0 ) = {2.0,3.0}; C( 2, 1 ) = {4.0,5.0}; C( 2, 2 ) = {6.0,7.0};

        bool g = moris::isempty( C );
        REQUIRE( g == false );

        moris::Mat< moris::cplx > D( 3, 3 );
        bool h = moris::isempty( D );
        REQUIRE_FALSE( h );
    }

    //Testing the functionality of isempty(a) to see if a vector is empty or not
    SECTION( "moris::Cell" )
    {
        moris::Mat< moris::real > v( 3, 1 );
        v( 0, 0 ) = 1.0; v( 1, 0 ) = 2.0; v( 2, 0 ) = 3.0;

        bool i = moris::isempty( v );
        REQUIRE_FALSE( i );

        moris::Mat< moris::real > w;
        bool j = moris::isempty( w );
        REQUIRE( j == true );
    }

    SECTION( "moris::Cell cplx" )
    {
        moris::Mat< moris::cplx > x( 3, 1 );
        x( 0, 0 ) = {0.0, 1.0}; x( 1, 0 ) = {2.0, 3.0}; x( 2, 0 ) = {4.0, 5.0};

        bool k = moris::isempty( x );
        REQUIRE_FALSE( k );

        moris::Mat< moris::cplx > y;
        bool l = moris::isempty( y );
        REQUIRE( l == true );
    }
}
