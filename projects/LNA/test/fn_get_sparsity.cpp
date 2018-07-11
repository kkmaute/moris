// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// C++ header files
#include <stdio.h>

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::get_sparsity",
        "[linalgebra],[moris::Sp_Mat]" )
{
    #include "linalg/fn_get_sparsity.inc"

    REQUIRE( Sparsity( 0, 0 ) == 4 ); REQUIRE( Sparsity( 0, 1 ) == 0 );
    REQUIRE( Sparsity( 1, 0 ) == 1 ); REQUIRE( Sparsity( 1, 1 ) == 2 );
    REQUIRE( Sparsity( 2, 0 ) == 2 ); REQUIRE( Sparsity( 2, 1 ) == 3 );
    REQUIRE( Sparsity( 3, 0 ) == 3 ); REQUIRE( Sparsity( 3, 1 ) == 4 );
}
