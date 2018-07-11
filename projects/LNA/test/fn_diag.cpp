// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::diag",
        "[linalgebra],[diag]" )
{
    // Testing the functionality of diag(A) to extract
    // the diagonal entries of a matrix A.
    SECTION("moris::diagvec" )
    {
#include "linalg/fn_diag/diagmat.inc"
        REQUIRE( v( 0, 0 ) == 1.0 );
        REQUIRE( v( 1, 0 ) == 5.0 );
        REQUIRE( v( 2, 0 ) == 9.0 );
    }


    // Testing the functionality of diag(v) to interpret
    // the vector v as a diagonal matrix.
    SECTION("moris::diagmat" )

    {
#include "linalg/fn_diag/diagvec.inc"

        REQUIRE( C( 0,0 ) == 1.0 );
        REQUIRE( C( 0,1 ) == 0.0 );
        REQUIRE( C( 0,2 ) == 0.0 );

        REQUIRE( C( 1,0 ) == 0.0 );
        REQUIRE( C( 1,1 ) == 2.0 );
        REQUIRE( C( 1,2 ) == 0.0 );

        REQUIRE( C( 2,0 ) == 0.0 );
        REQUIRE( C( 2,1 ) == 0.0 );
        REQUIRE( C( 2,2 ) == 3.0 );
    }

}
