// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::issquare",
        "[linalgebra],[issquare]" )
{
    SECTION( "issquare( moris::Mat )" )
    {
        #include "linalg/fn_issquare.inc"

        REQUIRE( e == true  );
    }
}
