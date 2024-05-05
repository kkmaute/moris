/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eye.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_diag_mat.hpp"

namespace moris
{
    TEST_CASE(
            "moris::diag_mat",
            "[linalgebra],[eye]" )
    {
        // generate diagonal matrix from vector

        // create vector
        Matrix< moris::DDRMat > A( 3, 1 );

        A( 0, 0 ) = -1.0;
        A( 1, 0 ) = 1.0;
        A( 2, 0 ) = 3.0;

        // generate matrix
        Matrix< moris::DDRMat > Adiag = diag_mat( A );

        // check result
        REQUIRE( Adiag( 0, 0 ) == -1.0 );
        REQUIRE( Adiag( 0, 1 ) == 0.0 );
        REQUIRE( Adiag( 0, 2 ) == 0.0 );

        REQUIRE( Adiag( 1, 0 ) == 0.0 );
        REQUIRE( Adiag( 1, 1 ) == 1.0 );
        REQUIRE( Adiag( 1, 2 ) == 0.0 );

        REQUIRE( Adiag( 2, 0 ) == 0.0 );
        REQUIRE( Adiag( 2, 1 ) == 0.0 );
        REQUIRE( Adiag( 2, 2 ) == 3.0 );

        // generate matrix using vector on 1st off-diagonal
        Adiag = diag_mat( A, 1 );

        // check result
        REQUIRE( Adiag( 0, 0 ) == 0.0 );
        REQUIRE( Adiag( 0, 1 ) == -1.0 );
        REQUIRE( Adiag( 0, 2 ) == 0.0 );
        REQUIRE( Adiag( 0, 3 ) == 0.0 );

        REQUIRE( Adiag( 1, 0 ) == 0.0 );
        REQUIRE( Adiag( 1, 1 ) == 0.0 );
        REQUIRE( Adiag( 1, 2 ) == 1.0 );
        REQUIRE( Adiag( 1, 3 ) == 0.0 );

        REQUIRE( Adiag( 2, 0 ) == 0.0 );
        REQUIRE( Adiag( 2, 1 ) == 0.0 );
        REQUIRE( Adiag( 2, 2 ) == 0.0 );
        REQUIRE( Adiag( 2, 3 ) == 3.0 );

        REQUIRE( Adiag( 3, 0 ) == 0.0 );
        REQUIRE( Adiag( 3, 1 ) == 0.0 );
        REQUIRE( Adiag( 3, 2 ) == 0.0 );
        REQUIRE( Adiag( 3, 3 ) == 0.0 );

        // generate diagonal matrix from matrix

        // create vector
        A.set_size( 3, 3 );

        A( 0, 0 ) = -1.0;
        A( 0, 1 ) = 1.0;
        A( 0, 2 ) = 3.0;

        A( 1, 0 ) = 1.0;
        A( 1, 1 ) = 2.0;
        A( 1, 2 ) = 0.0;

        A( 2, 0 ) = 3.0;
        A( 2, 1 ) = 0.0;
        A( 2, 2 ) = 2.0;

        // generate matrix
        Adiag = diag_mat( A );

        // check result
        REQUIRE( Adiag( 0, 0 ) == -1.0 );
        REQUIRE( Adiag( 0, 1 ) == 0.0 );
        REQUIRE( Adiag( 0, 2 ) == 0.0 );

        REQUIRE( Adiag( 1, 0 ) == 0.0 );
        REQUIRE( Adiag( 1, 1 ) == 2.0 );
        REQUIRE( Adiag( 1, 2 ) == 0.0 );

        REQUIRE( Adiag( 2, 0 ) == 0.0 );
        REQUIRE( Adiag( 2, 1 ) == 0.0 );
        REQUIRE( Adiag( 2, 2 ) == 2.0 );

        // generate matrix
        Adiag = diag_mat( A, 1 );

        // check result
        REQUIRE( Adiag( 0, 0 ) == 0.0 );
        REQUIRE( Adiag( 0, 1 ) == 1.0 );
        REQUIRE( Adiag( 0, 2 ) == 0.0 );

        REQUIRE( Adiag( 1, 0 ) == 0.0 );
        REQUIRE( Adiag( 1, 1 ) == 0.0 );
        REQUIRE( Adiag( 1, 2 ) == 0.0 );

        REQUIRE( Adiag( 2, 0 ) == 0.0 );
        REQUIRE( Adiag( 2, 1 ) == 0.0 );
        REQUIRE( Adiag( 2, 2 ) == 0.0 );
    }
}    // namespace moris
