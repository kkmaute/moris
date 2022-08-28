/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_diag.cpp
 *
 */

#include "catch.hpp"

#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_diag_mat.hpp"
#include "fn_diag_vec.hpp"

namespace moris
{

TEST_CASE(
        "moris::diag",
        "[linalgebra],[diag]" )
{
    // Testing the functionality of diag(A) to extract
    // the diagonal entries of a matrix A.
    SECTION("moris::diagvec" )
    {
        Matrix< DDRMat > A( 3, 3 );
        Matrix< DDRMat > B( 3, 3 );

        A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
        A( 2, 0 ) = 7.0; A( 2, 1 ) = 8.0; A( 2, 2 ) = 9.0;

        Matrix< DDRMat > v = diag_vec( A );

        REQUIRE( v( 0, 0 ) == 1.0 );
        REQUIRE( v( 1, 0 ) == 5.0 );
        REQUIRE( v( 2, 0 ) == 9.0 );
    }

    // Testing the functionality of diag(v) to interpret
    // the vector v as a diagonal matrix.
    SECTION("moris::diagmat" )

    {
        Matrix< DDRMat > w( 3,1 );
        w( 0,0 ) = 1.0; w( 1,0 ) = 2.0; w( 2,0 ) = 3.0;

        Matrix< DDRMat > C = diag_mat( w );

        REQUIRE( C( 0,0 ) == 1.0 );
        REQUIRE( C( 0,1 ) == 0.0 );
        REQUIRE( C( 0,2 ) == 0.0 );

        REQUIRE( C( 1,0 ) == 0.0 );
        REQUIRE( C( 1,1 ) == 2.0 );
        REQUIRE( C( 1,2 ) == 0.0 );

        REQUIRE( C( 2,0 ) == 0.0 );
        REQUIRE( C( 2,1 ) == 0.0 );
        REQUIRE( C( 2,2 ) == 3.0 );
    }

}
}

