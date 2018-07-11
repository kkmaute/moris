// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::lu",
        "[linalgebra],[lu]" )
        {
    SECTION( "lu( moris::Mat, moris::Mat, moris::Mat, moris::Mat )" )
            {
        #include "linalg/fn_lu.inc"

        REQUIRE( moris::equal_to(L( 0,0 ) , 1.000000000000000e+00 ));
        REQUIRE( moris::equal_to(L( 0,1 ) , 0.000000000000000e+00 ));
        REQUIRE( moris::equal_to(L( 0,2 ) , 0.000000000000000e+00 ));

        REQUIRE( moris::equal_to(L( 1,0 ) , 6.666666666666666e-01 ));
        REQUIRE( moris::equal_to(L( 1,1 ) , 1.000000000000000e+00 ));
        REQUIRE( moris::equal_to(L( 1,2 ) , 0.000000000000000e+00 ));

        REQUIRE( moris::equal_to(L( 2,0 ) , 3.333333333333333e-01 ));
        REQUIRE( moris::equal_to(L( 2,1 ) , -9.999999999999993e-01 ));
        REQUIRE( moris::equal_to(L( 2,2 ) , 1.000000000000000e+00 ));

        REQUIRE( moris::equal_to(U( 0,0 ) , 3.000000000000000e+00 ));
        REQUIRE( moris::equal_to(U( 0,1 ) , 8.000000000000000e+00 ));
        REQUIRE( moris::equal_to(U( 0,2 ) , 1.400000000000000e+01 ));

        REQUIRE( moris::equal_to(U( 1,0 ) , 0.000000000000000e+00 ));
        REQUIRE( moris::equal_to(U( 1,1 ) , 6.666666666666670e-01 ));
        REQUIRE( moris::equal_to(U( 1,2 ) , 3.666666666666668e+00 ));

        REQUIRE( moris::equal_to(U( 2,0 ) , 0.000000000000000e+00 ));
        REQUIRE( moris::equal_to(U( 2,1 ) , 0.000000000000000e+00 ));
        REQUIRE( moris::equal_to(U( 2,2 ) , 2.999999999999999e+00 ));

        REQUIRE( moris::equal_to(P( 0,0 ) , 0.0000e+00 ));
        REQUIRE( moris::equal_to(P( 0,1 ) , 1.0000e+00 ));
        REQUIRE( moris::equal_to(P( 0,2 ) , 0.0000e+00 ));

        REQUIRE( moris::equal_to(P( 1,0 ) , 0.0000e+00 ));
        REQUIRE( moris::equal_to(P( 1,1 ) , 0.0000e+00 ));
        REQUIRE( moris::equal_to(P( 1,2 ) , 1.0000e+00 ));

        REQUIRE( moris::equal_to(P( 2,0 ) , 1.0000e+00 ));
        REQUIRE( moris::equal_to(P( 2,1 ) , 0.0000e+00 ));
        REQUIRE( moris::equal_to(P( 2,2 ) , 0.0000e+00 ));
            }
        }
