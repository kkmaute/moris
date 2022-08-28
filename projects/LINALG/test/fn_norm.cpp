/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_norm.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_norm.hpp"

namespace moris
{
TEST_CASE(
        "moris::norm",
        "[linalgebra],[norm]" )
    {
    Matrix< DDRMat > a( 3, 3 );

    a( 0, 0 ) = 1.0; a( 0, 1 ) = 2.0; a( 0, 2 ) = 3.0;
    a( 1, 0 ) = 4.0; a( 1, 1 ) = 5.0; a( 1, 2 ) = 6.0;
    a( 2, 0 ) = 9.0; a( 2, 1 ) = 8.0; a( 2, 2 ) = 9.0;

    Matrix< DDRMat > rMat( 3, 3, 123.5 );
    Matrix< DDRMat >aVec = a.get_column(0);

    moris::real VecNorm = norm(aVec);

    REQUIRE( moris::equal_to( VecNorm,  9.899494, 1.0e+11 ));
    }
}

