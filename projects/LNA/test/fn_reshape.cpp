// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "fn_reshape.hpp" // LNA/src

TEST_CASE(
        "moris::reshape",
        "[linalgebra],[reshape]" )
{
    //Demonstrates the functionality of "reshape". You can reshape a vector or matrix into different configurations and it preserves the data
    //The elements in the generated object are placed column-wise (ie. the first column is filled up before filling the second column)
    SECTION( "reshape of uint row vector" )
    {
    moris::Mat< moris::uint > A( 6, 1 );
    A( 0, 0 ) = 1;
    A( 1, 0 ) = 2;
    A( 2, 0 ) = 3;
    A( 3, 0 ) = 2;
    A( 4, 0 ) = 2;
    A( 5, 0 ) = 4;

    moris::Mat< moris::uint > C = reshape(A,2,3);
    REQUIRE( moris::equal_to( C(0,0), 1 ) );
    REQUIRE( moris::equal_to( C(1,0), 2 ) );
    REQUIRE( moris::equal_to( C(0,1), 3 ) );
    REQUIRE( moris::equal_to( C(1,1), 2 ) );
    REQUIRE( moris::equal_to( C(0,2), 2 ) );
    REQUIRE( moris::equal_to( C(1,2), 4 ) );
    }

    SECTION( "reshape of real row vector" )
    {
    moris::Mat< moris::real > A( 6, 1 );
    A( 0, 0 ) = 1.1;
    A( 1, 0 ) = 2.1;
    A( 2, 0 ) = 2.1;
    A( 3, 0 ) = 3.1;
    A( 4, 0 ) = 2.1;
    A( 5, 0 ) = 4.1;

    moris::Mat< moris::real > C = reshape(A,3,2);
    REQUIRE( moris::equal_to( C(0,0), 1.1 ) );
    REQUIRE( moris::equal_to( C(1,0), 2.1 ) );
    REQUIRE( moris::equal_to( C(2,0), 2.1 ) );
    REQUIRE( moris::equal_to( C(0,1), 3.1 ) );
    REQUIRE( moris::equal_to( C(1,1), 2.1 ) );
    REQUIRE( moris::equal_to( C(2,1), 4.1 ) );
    }

    SECTION( "reshape of uint col vector" )
    {
    moris::Mat< moris::uint > A( 1,6 );
    A( 0, 0 ) = 1;
    A( 0, 1 ) = 2;
    A( 0, 2 ) = 2;
    A( 0, 3 ) = 3;
    A( 0, 4 ) = 1;
    A( 0, 5 ) = 2;

    moris::Mat< moris::uint > C = reshape(A,2,3);
    REQUIRE( moris::equal_to( C(0,0), 1 ) );
    REQUIRE( moris::equal_to( C(1,0), 2 ) );
    REQUIRE( moris::equal_to( C(0,1), 2 ) );
    REQUIRE( moris::equal_to( C(1,1), 3 ) );
    REQUIRE( moris::equal_to( C(0,2), 1 ) );
    REQUIRE( moris::equal_to( C(1,2), 2 ) );
    }

    SECTION( "reshape of real col vector" )
    {
    moris::Mat< moris::real > A( 1,8 );
    A( 0, 0 ) = 1.1;
    A( 0, 1 ) = 2.1;
    A( 0, 2 ) = 2.1;
    A( 0, 3 ) = 3.1;
    A( 0, 4 ) = 1.1;
    A( 0, 5 ) = 4.1;
    A( 0, 6 ) = 2.3;
    A( 0, 7 ) = 5.3;

    moris::uint col = 2;
    moris::uint row = 4;
    moris::Mat< moris::real > C = reshape(A,row,col);
    REQUIRE( moris::equal_to( C(0,0), 1.1 ) );
    REQUIRE( moris::equal_to( C(1,0), 2.1 ) );
    REQUIRE( moris::equal_to( C(2,0), 2.1 ) );
    REQUIRE( moris::equal_to( C(3,0), 3.1 ) );
    REQUIRE( moris::equal_to( C(0,1), 1.1 ) );
    REQUIRE( moris::equal_to( C(1,1), 4.1 ) );
    REQUIRE( moris::equal_to( C(2,1), 2.3 ) );
    REQUIRE( moris::equal_to( C(3,1), 5.3 ) );
    }

    SECTION( "reshape of a mat" )
     {
     moris::Mat< moris::real > A( 2,2 );

     A( 0, 0 ) = 1.1;
     A( 0, 1 ) = 2.1;
     A( 1, 0 ) = 2.1;
     A( 1, 1 ) = 3.1;

     moris::Mat< moris::real > C = reshape(A,1,4);
     REQUIRE( moris::equal_to( C(0,0), 1.1 ) );
     REQUIRE( moris::equal_to( C(0,1), 2.1 ) );
     REQUIRE( moris::equal_to( C(0,2), 2.1 ) );
     REQUIRE( moris::equal_to( C(0,3), 3.1 ) );
     }

}



