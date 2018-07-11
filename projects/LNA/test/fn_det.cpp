// Third-party header files.
#include <catch.hpp>
#include <iostream>
#include <stdlib.h>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
                "moris::det",
                "[linalgebra],[det]" )
{

    #include "linalg/fn_det.inc"

    SECTION( "moris::det real" )
    {
        REQUIRE( moris::equal_to( detr, 1.2649337e+01 ) );
    }

}
