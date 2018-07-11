// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "fn_find.hpp" // LNA/src

TEST_CASE(
        "moris::find",
        "[linalgebra],[find]" )
{
    //Demonstrates the functionality of "find". It provides a vector, in which every value from a condition is found.
    // Find requires a vector which contains zeros and ones. They come from a condition, for example: A==B, A>=B, ...
    // The output is always a col vector
    SECTION( "find of uint row vector" )
    {
    //Uniqueness of row vector

    moris::Mat< moris::uint > A( 7, 1 );

    A( 0, 0 ) = 1;
    A( 1, 0 ) = 2;
    A( 2, 0 ) = 3;
    A( 3, 0 ) = 2;
    A( 4, 0 ) = 2;
    A( 5, 0 ) = 4;
    A( 6, 0 ) = 1;

    moris::Mat< moris::uint > B = (A==2);
    moris::Mat< moris::uint > C = find(B);
    REQUIRE( moris::equal_to( C(0,0), 1 ) );
    REQUIRE( moris::equal_to( C(1,0), 3 ) );
    REQUIRE( moris::equal_to( C(2,0), 4 ) );
    }

    SECTION( "find of real row vector" )
    {
    //Uniqueness of row vector

    moris::Mat< moris::real > A( 7, 1 );

    A( 0, 0 ) = 1.1;
    A( 1, 0 ) = 2.1;
    A( 2, 0 ) = 2.1;
    A( 3, 0 ) = 3.1;
    A( 4, 0 ) = 2.1;
    A( 5, 0 ) = 4.1;
    A( 6, 0 ) = 1.1;

    moris::Mat< moris::uint > B = (A==3.1);
    moris::Mat< moris::uint > C = find(B);
    REQUIRE( moris::equal_to( C(0,0), 3 ) );
    }

    SECTION( "find of uint col vector" )
    {
    //Uniqueness of col vector

    moris::Mat< moris::uint > A( 1,7 );

    A( 0, 0 ) = 1;
    A( 0, 1 ) = 2;
    A( 0, 2 ) = 2;
    A( 0, 3 ) = 3;
    A( 0, 4 ) = 1;
    A( 0, 5 ) = 2;
    A( 0, 6 ) = 4;

    moris::Mat< moris::uint > B = (A==2);
    moris::Mat< moris::uint > C = find(B);
    REQUIRE( moris::equal_to( C(0,0), 1 ) );
    REQUIRE( moris::equal_to( C(1,0), 2 ) );
    REQUIRE( moris::equal_to( C(2,0), 5 ) );
    }

    SECTION( "find of real col vector" )
    {
    //Uniqueness of col vector

    moris::Mat< moris::real > A( 1,8 );

    A( 0, 0 ) = 1.1;
    A( 0, 1 ) = 2.1;
    A( 0, 2 ) = 2.1;
    A( 0, 3 ) = 3.1;
    A( 0, 4 ) = 1.1;
    A( 0, 5 ) = 4.1;
    A( 0, 6 ) = 2.3;
    A( 0, 7 ) = 5.3;
    moris::Mat< moris::uint > B = (A==2.3);
    moris::Mat< moris::uint > C = find(B);
    REQUIRE( moris::equal_to( C(0,0), 6 ) );
    }

}



