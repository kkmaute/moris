// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::cond",
        "[linalgebra],[cond]" )
{
    #include "linalg/fn_cond.inc"

    REQUIRE( moris::equal_to( aCond, 4.984802690419728e+01 ) );
    REQUIRE( moris::equal_to( bCond, 9.312148994138860e+01 ) );

    /*
     Due to relative numerical differences in the svd (used for computation of
     the condition number of a matrix in EIGEN), the equal_to tolerance is set
     to 1.0e+7 * machine precision. Depending on the range of values in the
     matrix (e.g if greater than 1e6) this tolerance might not be big enough
     and cause the unit test to fail.
    */

    REQUIRE( moris::equal_to( cCond, 3.634808466391376e+07, 1.0e+07 ) );

}
