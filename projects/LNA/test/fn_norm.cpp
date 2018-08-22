/*
 * fn_norm.cpp
 *
 *  Created on: Aug 22, 2018
 *      Author: doble
 */

// Third-party header files.
#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src

#include "cl_Mat.hpp"
#include "fn_norm.hpp"

TEST_CASE(
        "moris::norm",
        "[linalgebra],[norm]" )
{
    #include "linalg/cl_Mat/Mat.inc"

    SECTION( "real moris::Mat::norm" )
    {
        moris::Mat< moris::real >aVec = a.col(0);

        moris::real VecNorm = norm(aVec);

        REQUIRE( moris::equal_to( VecNorm,  9.899494, 1.0e+11 ));
    }
}


