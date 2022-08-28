/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_less_equal.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_less_equal.hpp"

namespace moris
{
TEST_CASE(
        "moris::op_less_equal",
        "[linalgebra],[Mat],[Col],[op_less_equal]" )
{
    SECTION( "moris::Col <= moris::Col" )
    {
        Matrix< DDRMat > A( 3, 1 );
        Matrix< DDRMat > B( 3, 1 );

        A( 0, 0 ) = 1.0;
        A( 1, 0 ) = 0.0;
        A( 2, 0 ) = 3.0;

        B( 0, 0 ) = 1.0;
        B( 1, 0 ) = 2.0;
        B( 2, 0 ) = 1.0;

        Matrix< DDBMat > C = ( A <= B );

        REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
    }

    SECTION( "moris::Mat <= moris::Mat" )
    {
        Matrix< DDRMat > A( 3, 3 );
        Matrix< DDRMat > B( 3, 3 );

        A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
        A( 2, 0 ) = 7.0; A( 2, 1 ) = 8.0; A( 2, 2 ) = 9.0;

        B( 0, 0 ) = 3.0; B( 0, 1 ) = 2.0; B( 0, 2 ) = 1.0;
        B( 1, 0 ) = 6.0; B( 1, 1 ) = 5.0; B( 1, 2 ) = 4.0;
        B( 2, 0 ) = 9.0; B( 2, 1 ) = 8.0; B( 2, 2 ) = 7.0;

        Matrix< DDBMat > C = ( A <= B );        REQUIRE( C( 0, 0 ) == 1 );

        REQUIRE( C( 0, 1 ) == 1 );
        REQUIRE( C( 0, 2 ) == 0 );

        REQUIRE( C( 1, 0 ) == 1 );
        REQUIRE( C( 1, 1 ) == 1 );
        REQUIRE( C( 1, 2 ) == 0 );

        REQUIRE( C( 2, 0 ) == 1 );
        REQUIRE( C( 2, 1 ) == 1 );
        REQUIRE( C( 2, 2 ) == 0 );
    }

    SECTION( "moris::Col <= scalar" )
     {
         Matrix< DDRMat > A( 3, 1 );
         moris::real B;

         A( 0, 0 ) = 1.0;
         A( 1, 0 ) = 2.0;
         A( 2, 0 ) = 3.0;

         B = 2.0;

         Matrix< DDBMat > C = ( A <= B );

         REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
         REQUIRE( moris::equal_to( C.n_cols(), 1 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
     }

     SECTION( "moris::Mat <= scalar" )
     {
         Matrix< DDRMat > A( 3, 2 );
         moris::real B;

         A( 0, 0 ) = 1.0;    A( 0, 1 ) = 4.0;
         A( 1, 0 ) = 2.0;    A( 1, 1 ) = 0.0;
         A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

         B = 2.0;

         Matrix< DDBMat > C = ( A <= B );

         REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
         REQUIRE( moris::equal_to( C.n_cols(), 2 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 0, 1 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 1 ), 1 ) );
     }

     SECTION( "scalar <= moris::Col" )
     {
         Matrix< DDRMat > A( 3, 1 );
         moris::real B;

         A( 0, 0 ) = 1.0;
         A( 1, 0 ) = 2.0;
         A( 2, 0 ) = 3.0;

         B = 2.0;

         Matrix< DDBMat > C = ( B <= A );

         REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
         REQUIRE( moris::equal_to( C.n_cols(), 1 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
     }

     SECTION( "scalar <= moris::Mat " )
     {
         Matrix< DDRMat > A( 3, 2 );
         moris::real B;

         A( 0, 0 ) = 1.0;    A( 0, 1 ) = 1.0;
         A( 1, 0 ) = 2.0;    A( 1, 1 ) = 3.0;
         A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

         B = 2.0;

         Matrix< DDBMat > C = ( B <= A );

         REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
         REQUIRE( moris::equal_to( C.n_cols(), 2 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 0, 1 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 1 ), 1 ) );
     }
}
}

