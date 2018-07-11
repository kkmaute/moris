// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "fn_sum.hpp" // LNA/src

TEST_CASE(
        "moris::sum",
        "[linalgebra],[sum]" )
{
    //Demonstrates the functionality of "sum". You can sum a col or row vector or a matrix.
    SECTION( "Sum of row vector" )
    {
    //Sum of row vector

    moris::Mat< moris::real > A( 3, 1 );

    A( 0, 0 ) = 1.0;
    A( 1, 0 ) = 2.0;
    A( 2, 0 ) = 3.0;

    moris::uint C = sum(A);

    REQUIRE( moris::equal_to( C, 6 ) );
    }

    SECTION( "Sum of col vector" )
    {
    //Sum of col vector
    moris::Mat< moris::real > B( 1, 3 );

    B( 0, 0 ) = 1.0;
    B( 0, 1 ) = 2.0;
    B( 0, 2 ) = 3.0;

    moris::uint C = sum(B);

    REQUIRE( moris::equal_to( C, 6 ) );
    }

    SECTION( "Sum of a matrix" )
    {
    //Sum of col vector
    moris::Mat< moris::real > B( 3, 3 );

    B( 0, 0 ) = 1.0;    B( 0, 1 ) = 1.0;     B( 0, 2 ) = 1.0;
    B( 1, 0 ) = 2.0;    B( 1, 1 ) = 2.0;     B( 1, 2 ) = 2.0;
    B( 2, 0 ) = 3.0;    B( 2, 1 ) = 3.0;     B( 2, 2 ) = 3.0;

    moris::uint C = sum(B);

    REQUIRE( moris::equal_to( C, 18 ) );
    }

}



