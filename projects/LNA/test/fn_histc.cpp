// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "fn_histc.hpp" // LNA/src
#include "fn_unique.hpp" // LNA/src

TEST_CASE(
        "moris::histc",
        "[linalgebra],[histc]" )
{
    //Demonstrates the functionality of "histc". It provides a vector, that contains the counts of the number of values that fall between the elements in a given range.
    // histc first sorts the vector and takes then the counts of the number of values.
    // Only a coloumn vector as input is allowed
    SECTION( "histc of uint col vector" )
    {
    moris::Mat< moris::uint > A( 7, 1 );

    A( 0, 0 ) = 1;
    A( 1, 0 ) = 2;
    A( 2, 0 ) = 3;
    A( 3, 0 ) = 2;
    A( 4, 0 ) = 2;
    A( 5, 0 ) = 4;
    A( 6, 0 ) = 1;

    moris::Mat< moris::uint > C = histc(A,unique(A));
    REQUIRE( moris::equal_to( C(0,0), 2 ) );
    REQUIRE( moris::equal_to( C(1,0), 3 ) );
    REQUIRE( moris::equal_to( C(2,0), 1 ) );
    REQUIRE( moris::equal_to( C(3,0), 1 ) );
    }


}



