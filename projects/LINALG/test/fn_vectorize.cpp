/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_vectorize.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_vectorize.hpp"
#include "fn_reshape.hpp"
#include "fn_trans.hpp"

namespace moris
{
    TEST_CASE(
            "vectorize",
            "[linalgebra],[vectorize]" )
    {
        // Demonstrates the functionality of "vectorize". You can vectorize a vector or matrix into either a row or column vector.

        SECTION( "vectorize of rectangular matrix" )
        {
            // create rectangular matrix of size 4x2
            Matrix< DDRMat > A( 4, 2 );
            A( 0, 0 ) = 1.1;
            A( 1, 0 ) = 2.1;
            A( 2, 0 ) = 3.1;
            A( 3, 0 ) = 4.1;
            A( 0, 1 ) = 1.2;
            A( 1, 1 ) = 2.2;
            A( 2, 1 ) = 3.2;
            A( 3, 1 ) = 4.2;

            // flatten matrix into column vector
            Matrix< DDRMat > B = vectorize( A, false );

            // check for proper dimensions
            REQUIRE( B.n_rows() == 8 );
            REQUIRE( B.n_cols() == 1 );

            // check for proper values
            REQUIRE( equal_to( B( 0 ), 1.1 ) );
            REQUIRE( equal_to( B( 1 ), 2.1 ) );
            REQUIRE( equal_to( B( 2 ), 3.1 ) );
            REQUIRE( equal_to( B( 3 ), 4.1 ) );
            REQUIRE( equal_to( B( 4 ), 1.2 ) );
            REQUIRE( equal_to( B( 5 ), 2.2 ) );
            REQUIRE( equal_to( B( 6 ), 3.2 ) );
            REQUIRE( equal_to( B( 7 ), 4.2 ) );

            // flatten matrix into row vector
            Matrix< DDRMat > C = vectorize( A, true );

            // check for proper dimensions
            REQUIRE( C.n_rows() == 1 );
            REQUIRE( C.n_cols() == 8 );

            // check for proper values
            REQUIRE( equal_to( C( 0 ), 1.1 ) );
            REQUIRE( equal_to( C( 1 ), 1.2 ) );
            REQUIRE( equal_to( C( 2 ), 2.1 ) );
            REQUIRE( equal_to( C( 3 ), 2.2 ) );
            REQUIRE( equal_to( C( 4 ), 3.1 ) );
            REQUIRE( equal_to( C( 5 ), 3.2 ) );
            REQUIRE( equal_to( C( 6 ), 4.1 ) );
            REQUIRE( equal_to( C( 7 ), 4.2 ) );
        }

        SECTION( "vectorize of a column vector into column and row vector" )
        {
            // create column vector of size 4
            Matrix< DDRMat > A( 4, 1 );

            A( 0, 0 ) = 1.1;
            A( 1, 0 ) = 2.1;
            A( 2, 0 ) = 3.1;
            A( 3, 0 ) = 4.1;

            // flatten matrix into column vector - nothing should happen
            Matrix< DDRMat > B = vectorize( A );

            // check for proper dimensions
            REQUIRE( B.n_rows() == 4 );
            REQUIRE( B.n_cols() == 1 );

            // check for proper values
            REQUIRE( equal_to( B( 0 ), 1.1 ) );
            REQUIRE( equal_to( B( 1 ), 2.1 ) );
            REQUIRE( equal_to( B( 2 ), 3.1 ) );
            REQUIRE( equal_to( B( 3 ), 4.1 ) );

            // flatten matrix into row vector
            Matrix< DDRMat > C = vectorize( A, true );

            // check for proper dimensions
            REQUIRE( C.n_rows() == 1 );
            REQUIRE( C.n_cols() == 4 );

            // check for proper values
            REQUIRE( equal_to( C( 0 ), 1.1 ) );
            REQUIRE( equal_to( C( 1 ), 2.1 ) );
            REQUIRE( equal_to( C( 2 ), 3.1 ) );
            REQUIRE( equal_to( C( 3 ), 4.1 ) );
        }

        SECTION( "vectorize of a product of two matrices" )
        {
            // create two matrices: 2x1 and 1x2
            Matrix< DDRMat > A( 2, 1 );
            A( 0, 0 ) = 1.0;
            A( 1, 0 ) = 3.0;

            Matrix< DDRMat > B( 1, 2 );
            B( 0, 0 ) = 1.0;
            B( 0, 1 ) = 2.0;

            // flatten the matrix product into column vector
            Matrix< DDRMat > C = vectorize( A * B );

            // check for proper dimensions
            REQUIRE( C.n_rows() == 4 );
            REQUIRE( C.n_cols() == 1 );

            // check for proper values
            REQUIRE( equal_to( C( 0 ), 1.0 ) );
            REQUIRE( equal_to( C( 1 ), 3.0 ) );
            REQUIRE( equal_to( C( 2 ), 2.0 ) );
            REQUIRE( equal_to( C( 3 ), 6.0 ) );
        }

        SECTION( "vectorize versus reshape of matrices" )
        {
            // create two matrices: 3x4 and 1x2
            Matrix< DDRMat > A( 3, 4 );
            A( 0, 0 ) = 1.1;
            A( 0, 1 ) = 1.2;
            A( 0, 2 ) = 1.3;
            A( 0, 3 ) = 1.4;
            A( 1, 0 ) = 2.1;
            A( 1, 1 ) = 2.2;
            A( 1, 2 ) = 2.3;
            A( 1, 3 ) = 2.4;
            A( 2, 0 ) = 3.1;
            A( 2, 1 ) = 3.2;
            A( 2, 2 ) = 3.3;
            A( 2, 3 ) = 3.4;

            // create a 2x6 matrix using reshape
            Matrix< DDRMat > B = reshape( A, 2, 6 );

            // check for proper dimensions
            REQUIRE( B.n_rows() == 2 );
            REQUIRE( B.n_cols() == 6 );

            // check for proper values
            REQUIRE( equal_to( B( 0, 0 ), 1.1 ) );
            REQUIRE( equal_to( B( 1, 0 ), 2.1 ) );
            REQUIRE( equal_to( B( 0, 1 ), 3.1 ) );
            REQUIRE( equal_to( B( 1, 1 ), 1.2 ) );
            REQUIRE( equal_to( B( 0, 2 ), 2.2 ) );
            REQUIRE( equal_to( B( 1, 2 ), 3.2 ) );
            REQUIRE( equal_to( B( 0, 3 ), 1.3 ) );
            REQUIRE( equal_to( B( 1, 3 ), 2.3 ) );
            REQUIRE( equal_to( B( 0, 4 ), 3.3 ) );
            REQUIRE( equal_to( B( 1, 4 ), 1.4 ) );
            REQUIRE( equal_to( B( 0, 5 ), 2.4 ) );
            REQUIRE( equal_to( B( 1, 5 ), 3.4 ) );

            // create a 12x1 column vector of A using reshape
            Matrix< DDRMat > C = reshape( A, 12, 1 );

            // check for proper dimensions
            REQUIRE( C.n_rows() == 12 );
            REQUIRE( C.n_cols() == 1 );

            // check for proper values
            REQUIRE( equal_to( C( 0 ), 1.1 ) );
            REQUIRE( equal_to( C( 1 ), 2.1 ) );
            REQUIRE( equal_to( C( 2 ), 3.1 ) );
            REQUIRE( equal_to( C( 3 ), 1.2 ) );
            REQUIRE( equal_to( C( 4 ), 2.2 ) );
            REQUIRE( equal_to( C( 5 ), 3.2 ) );
            REQUIRE( equal_to( C( 6 ), 1.3 ) );
            REQUIRE( equal_to( C( 7 ), 2.3 ) );
            REQUIRE( equal_to( C( 8 ), 3.3 ) );
            REQUIRE( equal_to( C( 9 ), 1.4 ) );
            REQUIRE( equal_to( C( 10 ), 2.4 ) );
            REQUIRE( equal_to( C( 11 ), 3.4 ) );

            // create a 12x1 column vector of A using vectorize
            Matrix< DDRMat > D = vectorize( A );

            // check the result from reshape and vectorize are identical
            REQUIRE( norm( D - C ) < 1e-12 );

            // create a 1x12 row vector using reshape
            Matrix< DDRMat > E = reshape( A, 1, 12 );

            // check the results from reshape are identical
            REQUIRE( norm( D - trans( E ) ) < 1e-12 );

            // create a 1x12 row vector using vectorize
            Matrix< DDRMat > F = vectorize( A, true );

            // check the result from reshape and vectorize are different !
            REQUIRE( norm( E - F ) > 1e-1 );

            // check for proper values
            REQUIRE( equal_to( F( 0 ), 1.1 ) );
            REQUIRE( equal_to( F( 1 ), 1.2 ) );
            REQUIRE( equal_to( F( 2 ), 1.3 ) );
            REQUIRE( equal_to( F( 3 ), 1.4 ) );
            REQUIRE( equal_to( F( 4 ), 2.1 ) );
            REQUIRE( equal_to( F( 5 ), 2.2 ) );
            REQUIRE( equal_to( F( 6 ), 2.3 ) );
            REQUIRE( equal_to( F( 7 ), 2.4 ) );
            REQUIRE( equal_to( F( 8 ), 3.1 ) );
            REQUIRE( equal_to( F( 9 ), 3.2 ) );
            REQUIRE( equal_to( F( 10 ), 3.3 ) );
            REQUIRE( equal_to( F( 11 ), 3.4 ) );
        }
    }
}    // namespace moris
