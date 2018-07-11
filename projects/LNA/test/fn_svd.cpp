// Third-party header files.
#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::svd",
        "[linalgebra],[svd]" )
{
#include "linalg/fn_svd/fn_svd_short.inc"
#include "linalg/fn_svd/fn_svd_long.inc"

    /*
     Due to relative numerical deviation in the svd in EIGEN, the unit tests
     might fail if the range of values in the matrix (e.g if greater than 1e6)
     is getting to big. The equal_to tolerance therefore might have to be
     adjusted to not cause the unit tests to fail.
     (e.g. 1.0e+7 * machine precision)
    */

    SECTION( "moris::svd1" )
    {
        REQUIRE( moris::equal_to( aS2( 0, 0 ), aS( 0, 0 ) ) );
        REQUIRE( moris::equal_to( aS2( 1, 0 ), aS( 1, 0 ) ) );
        REQUIRE( moris::equal_to( aS2( 2, 0 ), aS( 2, 0 ) ) );
    }

    SECTION( "moris::svd1 cplx" )
    {
        REQUIRE( moris::equal_to( aSc2( 0, 0 ), aSc( 0, 0 ) ) );
        REQUIRE( moris::equal_to( aSc2( 1, 0 ), aSc( 1, 0 ) ) );
        REQUIRE( moris::equal_to( aSc2( 2, 0 ), aSc( 2, 0 ) ) );
    }


    SECTION( "moris::svd2" )
    {
        moris::Mat< moris::real > aAtest( 3, 3 );
        aAtest = aU * moris::diag( aS2 ) * moris::trans( aV );
        REQUIRE( moris::equal_to( aA( 0, 0 ), aAtest( 0, 0 ) ) );
        REQUIRE( moris::equal_to( aA( 1, 0 ), aAtest( 1, 0 ) ) );
        REQUIRE( moris::equal_to( aA( 2, 0 ), aAtest( 2, 0 ) ) );
        REQUIRE( moris::equal_to( aA( 0, 1 ), aAtest( 0, 1 ) ) );
        REQUIRE( moris::equal_to( aA( 1, 1 ), aAtest( 1, 1 ) ) );
        REQUIRE( moris::equal_to( aA( 2, 1 ), aAtest( 2, 1 ) ) );
        REQUIRE( moris::equal_to( aA( 0, 2 ), aAtest( 0, 2 ) ) );
        REQUIRE( moris::equal_to( aA( 1, 2 ), aAtest( 1, 2 ) ) );
        REQUIRE( moris::equal_to( aA( 2, 2 ), aAtest( 2, 2 ) ) );
    }

    SECTION( "moris::svd2 cplx" )
    {
        moris::Mat< moris::cplx > aActest( 3, 3 );
        aActest = aUc * moris::diag( aSc2 ) * moris::ctrans( aVc );
        REQUIRE( aActest( 0, 0 ).real() - aAc2( 0, 0 ).real() < 1.0e-10 );
        REQUIRE( aActest( 1, 0 ).real() - aAc2( 1, 0 ).real() < 1.0e-10 );
        REQUIRE( aActest( 2, 0 ).real() - aAc2( 2, 0 ).real() < 1.0e-10 );
        REQUIRE( aActest( 0, 1 ).real() - aAc2( 0, 1 ).real() < 1.0e-10 );
        REQUIRE( aActest( 1, 1 ).real() - aAc2( 1, 1 ).real() < 1.0e-10 );
        REQUIRE( aActest( 2, 1 ).real() - aAc2( 2, 1 ).real() < 1.0e-10 );
        REQUIRE( aActest( 0, 2 ).real() - aAc2( 0, 2 ).real() < 1.0e-10 );
        REQUIRE( aActest( 1, 2 ).real() - aAc2( 1, 2 ).real() < 1.0e-10 );
        REQUIRE( aActest( 2, 2 ).real() - aAc2( 2, 2 ).real() < 1.0e-10 );

        REQUIRE( aActest( 0, 0 ).imag() - aAc2( 0, 0 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 1, 0 ).imag() - aAc2( 1, 0 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 2, 0 ).imag() - aAc2( 2, 0 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 0, 1 ).imag() - aAc2( 0, 1 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 1, 1 ).imag() - aAc2( 1, 1 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 2, 1 ).imag() - aAc2( 2, 1 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 0, 2 ).imag() - aAc2( 0, 2 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 1, 2 ).imag() - aAc2( 1, 2 ).imag() < 1.0e-10 );
        REQUIRE( aActest( 2, 2 ).imag() - aAc2( 2, 2 ).imag() < 1.0e-10 );

    }

}
