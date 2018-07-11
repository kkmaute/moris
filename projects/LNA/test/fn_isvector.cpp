// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::isvector",
        "[linalgebra],[isvector]" )
{
    SECTION( "isvector( moris::Mat )" )
    {
        #include "linalg/fn_isvector.inc"
        REQUIRE( e == false );
        REQUIRE( f == true  );
        REQUIRE( g == true  );

    }
}
