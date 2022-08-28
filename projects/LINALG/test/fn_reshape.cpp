/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_reshape.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_reshape.hpp"

namespace moris
{
    TEST_CASE(
            "moris::reshape",
            "[linalgebra],[reshape]" )
    {
        //Demonstrates the functionality of "reshape". You can reshape a vector or matrix into different configurations and it preserves the data
        //The elements in the generated object are placed column-wise (ie. the first column is filled up before filling the second column)
        SECTION( "reshape of uint row vector" )
                    {
            Matrix< DDUMat > A( 6, 1 );
            A( 0, 0 ) = 1;
            A( 1, 0 ) = 2;
            A( 2, 0 ) = 3;
            A( 3, 0 ) = 2;
            A( 4, 0 ) = 2;
            A( 5, 0 ) = 4;

            Matrix< DDUMat > C = moris::reshape(A,2,3);
            REQUIRE( C(0,0) == 1 );
            REQUIRE( C(1,0) == 2 );
            REQUIRE( C(0,1) == 3 );
            REQUIRE( C(1,1) == 2 );
            REQUIRE( C(0,2) == 2 );
            REQUIRE( C(1,2) == 4 );

            A.reshape(2,3);
            REQUIRE( A(0,0) == 1 );
            REQUIRE( A(1,0) == 2 );
            REQUIRE( A(0,1) == 3 );
            REQUIRE( A(1,1) == 2 );
            REQUIRE( A(0,2) == 2 );
            REQUIRE( A(1,2) == 4 );
                    }

        SECTION( "reshape of real row vector" )
        {
            moris::Matrix< moris::DDRMat > A( 6, 1 );
            A( 0, 0 ) = 1.1;
            A( 1, 0 ) = 2.1;
            A( 2, 0 ) = 2.1;
            A( 3, 0 ) = 3.1;
            A( 4, 0 ) = 2.1;
            A( 5, 0 ) = 4.1;

            moris::Matrix< moris::DDRMat > C = reshape(A,3,2);
            REQUIRE( moris::equal_to( C(0,0), 1.1 ) );
            REQUIRE( moris::equal_to( C(1,0), 2.1 ) );
            REQUIRE( moris::equal_to( C(2,0), 2.1 ) );
            REQUIRE( moris::equal_to( C(0,1), 3.1 ) );
            REQUIRE( moris::equal_to( C(1,1), 2.1 ) );
            REQUIRE( moris::equal_to( C(2,1), 4.1 ) );

            A.reshape(3,2);
            REQUIRE( moris::equal_to( A(0,0), 1.1 ) );
            REQUIRE( moris::equal_to( A(1,0), 2.1 ) );
            REQUIRE( moris::equal_to( A(2,0), 2.1 ) );
            REQUIRE( moris::equal_to( A(0,1), 3.1 ) );
            REQUIRE( moris::equal_to( A(1,1), 2.1 ) );
            REQUIRE( moris::equal_to( A(2,1), 4.1 ) );
        }

        SECTION( "reshape of uint col vector" )
        {
            Matrix< DDUMat > A( 1,6 );
            A( 0, 0 ) = 1;
            A( 0, 1 ) = 2;
            A( 0, 2 ) = 2;
            A( 0, 3 ) = 3;
            A( 0, 4 ) = 1;
            A( 0, 5 ) = 2;

            Matrix< DDUMat > C = reshape(A,2,3);
            REQUIRE( C(0,0) == 1 );
            REQUIRE( C(1,0) == 2 );
            REQUIRE( C(0,1) == 2 );
            REQUIRE( C(1,1) == 3 );
            REQUIRE( C(0,2) == 1 );
            REQUIRE( C(1,2) == 2 );

            A.reshape(2,3);
            REQUIRE( A(0,0) == 1 );
            REQUIRE( A(1,0) == 2 );
            REQUIRE( A(0,1) == 2 );
            REQUIRE( A(1,1) == 3 );
            REQUIRE( A(0,2) == 1 );
            REQUIRE( A(1,2) == 2 );

        }

        SECTION( "reshape of real col vector" )
        {
            moris::Matrix< moris::DDRMat > A( 1,8 );
            A( 0, 0 ) = 1.1;
            A( 0, 1 ) = 2.1;
            A( 0, 2 ) = 2.1;
            A( 0, 3 ) = 3.1;
            A( 0, 4 ) = 1.1;
            A( 0, 5 ) = 4.1;
            A( 0, 6 ) = 2.3;
            A( 0, 7 ) = 5.3;

            moris::uint col = 2;
            moris::uint row = 4;

            moris::Matrix< moris::DDRMat > C = reshape(A,row,col);
            REQUIRE( moris::equal_to( C(0,0), 1.1 ) );
            REQUIRE( moris::equal_to( C(1,0), 2.1 ) );
            REQUIRE( moris::equal_to( C(2,0), 2.1 ) );
            REQUIRE( moris::equal_to( C(3,0), 3.1 ) );
            REQUIRE( moris::equal_to( C(0,1), 1.1 ) );
            REQUIRE( moris::equal_to( C(1,1), 4.1 ) );
            REQUIRE( moris::equal_to( C(2,1), 2.3 ) );
            REQUIRE( moris::equal_to( C(3,1), 5.3 ) );

            A.reshape(row,col);
            REQUIRE( moris::equal_to( A(0,0), 1.1 ) );
            REQUIRE( moris::equal_to( A(1,0), 2.1 ) );
            REQUIRE( moris::equal_to( A(2,0), 2.1 ) );
            REQUIRE( moris::equal_to( A(3,0), 3.1 ) );
            REQUIRE( moris::equal_to( A(0,1), 1.1 ) );
            REQUIRE( moris::equal_to( A(1,1), 4.1 ) );
            REQUIRE( moris::equal_to( A(2,1), 2.3 ) );
            REQUIRE( moris::equal_to( A(3,1), 5.3 ) );
        }

        SECTION( "reshape of a mat" )
        {
            moris::Matrix< moris::DDRMat > A( 2,2 );

            A( 0, 0 ) = 1.1;
            A( 0, 1 ) = 2.1;
            A( 1, 0 ) = 2.1;
            A( 1, 1 ) = 3.1;

            moris::Matrix< moris::DDRMat > C = reshape(A,1,4);
            REQUIRE( moris::equal_to( C(0,0), 1.1 ) );
            REQUIRE( moris::equal_to( C(0,1), 2.1 ) );
            REQUIRE( moris::equal_to( C(0,2), 2.1 ) );
            REQUIRE( moris::equal_to( C(0,3), 3.1 ) );

            A.reshape(1,4);
            REQUIRE( moris::equal_to( A(0,0), 1.1 ) );
            REQUIRE( moris::equal_to( A(0,1), 2.1 ) );
            REQUIRE( moris::equal_to( A(0,2), 2.1 ) );
            REQUIRE( moris::equal_to( A(0,3), 3.1 ) );
        }

        SECTION( "reshape of a product of two matrices" )
        {
            // create two matrices: 2x1 and 1x2
            Matrix< DDRMat > A( 2,1 );
            A( 0, 0 ) = 1.0;
            A( 1, 0 ) = 3.0;

            Matrix< DDRMat > B( 1,2 );
            B( 0, 0 ) = 1.0;
            B( 0, 1 ) = 2.0;

            // flatten the matrix product into column vector
            Matrix< DDRMat > C = reshape(A*B, 4, 1);

            // check for proper dimensions
            REQUIRE( C.n_rows() == 4 );
            REQUIRE( C.n_cols() == 1 );

            // check for proper values
            REQUIRE( equal_to( C(0), 1.0 ) );
            REQUIRE( equal_to( C(1), 3.0 ) );
            REQUIRE( equal_to( C(2), 2.0 ) );
            REQUIRE( equal_to( C(3), 6.0 ) );
        }
    }
}

