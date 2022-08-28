/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cond.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_cond.hpp"

namespace moris
{
TEST_CASE(
        "cond",
        "[linalgebra],[cond]" )
{
    Matrix< DDRMat > A( 3, 3 );
    Matrix< DDRMat > B( 3, 3 );
    Matrix< DDRMat > C( 3, 3 );

    A( 0, 0 ) = 1.0; A( 0, 1 ) = 1.0; A( 0, 2 ) = 2.0;
    A( 1, 0 ) = 3.0; A( 1, 1 ) = 4.0; A( 1, 2 ) = 5.0;
    A( 2, 0 ) = 6.0; A( 2, 1 ) = 7.0; A( 2, 2 ) = 7.0;

    B( 0, 0 ) = 1.0e-12; B( 0, 1 ) = 1.0; B( 0, 2 ) = 2.0;
    B( 1, 0 ) = 3.0;     B( 1, 1 ) = 4.0; B( 1, 2 ) = 5.0;
    B( 2, 0 ) = 6.0;     B( 2, 1 ) = 7.0; B( 2, 2 ) = 7.0;

    C( 0, 0 ) = 1.0; C( 0, 1 ) = 1.0; C( 0, 2 ) = 2.0;
    C( 1, 0 ) = 3.0; C( 1, 1 ) = 4.0; C( 1, 2 ) = 5.0;
    C( 2, 0 ) = 6.0; C( 2, 1 ) = 7.0; C( 2, 2 ) = 7.0e6;

    real aCond = cond(A);
    real bCond = cond(B);
    real cCond = cond(C);

    REQUIRE( equal_to( aCond, 4.984802690419728e+01 ) );
    REQUIRE( equal_to( bCond, 9.312148994138860e+01 ) );

    /*
     Due to relative numerical differences in the svd (used for computation of
     the condition number of a matrix in EIGEN), the equal_to tolerance is set
     to 1.0e+7 * machine precision. Depending on the range of values in the
     matrix (e.g if greater than 1e6) this tolerance might not be big enough
     and cause the unit test to fail.
    */

    REQUIRE( equal_to( cCond, 3.634808466391376e+07, 1.0e+07 ) );

}
}

