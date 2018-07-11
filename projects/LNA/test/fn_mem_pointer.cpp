// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// C++ header files
#include <stdio.h>

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::topointer",
        "[linalgebra],[mem_pointer]" )
{
    SECTION( "mem_pointer( moris::Mat )" )
    {
        moris::Mat< moris::real > MatA( 3, 3 );

        MatA( 0, 0 ) = 1.1; MatA( 0, 1 ) = 2.2; MatA( 0, 2 ) = 3.3;
        MatA( 1, 0 ) = 4.4; MatA( 1, 1 ) = 5.5; MatA( 1, 2 ) = 6.6;
        MatA( 2, 0 ) = 7.7; MatA( 2, 1 ) = 8.8; MatA( 2, 2 ) = 9.9;

        moris::real* PointerA = moris::mem_pointer( MatA );

        REQUIRE( std::abs( PointerA[0] - MatA( 0, 0 ) ) < 1.0e-12 );
        REQUIRE( std::abs( PointerA[3] - MatA( 0, 1 ) ) < 1.0e-12 );
        REQUIRE( std::abs( PointerA[8] - MatA( 2, 2 ) ) < 1.0e-12 );
    }
}
