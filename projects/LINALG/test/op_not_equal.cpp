/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_not_equal.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_not_equal.hpp"

namespace moris
{
TEST_CASE(
        "moris::op_not_equal",
        "[linalgebra],[op_not_equal]" )
{

    Matrix< DDRMat > Am( 3, 3 );
    Matrix< DDRMat > Bm( 3, 3 );

    Am( 0, 0 ) = 1.0; Am( 0, 1 ) = 1.0; Am( 0, 2 ) = 1.0;
    Am( 1, 0 ) = 1.0; Am( 1, 1 ) = 1.0; Am( 1, 2 ) = 1.0;
    Am( 2, 0 ) = 1.0; Am( 2, 1 ) = 1.0; Am( 2, 2 ) = 1.0;

    Bm( 0, 0 ) = 1.0; Bm( 0, 1 ) = 2.0; Bm( 0, 2 ) = 2.0;
    Bm( 1, 0 ) = 2.0; Bm( 1, 1 ) = 1.0; Bm( 1, 2 ) = 2.0;
    Bm( 2, 0 ) = 2.0; Bm( 2, 1 ) = 2.0; Bm( 2, 2 ) = 1.0;

    Matrix< DDBMat > Cm = ( Am != Bm );

    Matrix< DDRMat > Ac( 3, 1 );
    Matrix< DDRMat > Bc( 3, 1 );

    Ac( 0, 0 ) = 1.0;
    Ac( 1, 0 ) = 2.0;
    Ac( 2, 0 ) = 3.0;

    Bc( 0, 0 ) = 1.0;
    Bc( 1, 0 ) = 0.0;
    Bc( 2, 0 ) = 3.0;

    Matrix< DDBMat > Cc = ( Ac != Bc );

    SECTION( "moris::Mat != moris::Mat" )
    {
       REQUIRE( Cm( 0, 0 ) == 0 );
       REQUIRE( Cm( 0, 1 ) == 1 );
       REQUIRE( Cm( 0, 2 ) == 1 );

       REQUIRE( Cm( 1, 0 ) == 1 );
       REQUIRE( Cm( 1, 1 ) == 0 );
       REQUIRE( Cm( 1, 2 ) == 1 );

       REQUIRE( Cm( 2, 0 ) == 1 );
       REQUIRE( Cm( 2, 1 ) == 1 );
       REQUIRE( Cm( 2, 2 ) == 0 );
    }

    SECTION( "moris::Col != moris::Col" )
    {
        REQUIRE( Cc( 0, 0 ) == 0 );
        REQUIRE( Cc( 1, 0 ) == 1 );
        REQUIRE( Cc( 2, 0 ) == 0 );
    }

    SECTION( "moris::Col != scalar" )
    {
        Matrix< DDRMat > A( 3, 1 );
        moris::real B;

        A( 0, 0 ) = 1.0;
        A( 1, 0 ) = 2.0;
        A( 2, 0 ) = 3.0;

        B = 2.0;

        Matrix< DDBMat > C = ( A != B );

        REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
        REQUIRE( moris::equal_to( C.n_cols(), 1 ) );

        REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
        REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
    }

    SECTION( "moris::Mat != scalar" )
    {
        Matrix< DDRMat > A( 3, 2 );
        moris::real B;

        A( 0, 0 ) = 1.0;    A( 0, 1 ) = 1.0;
        A( 1, 0 ) = 2.0;    A( 1, 1 ) = 0.0;
        A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

        B = 2.0;

        Matrix< DDBMat > C = ( A != B );

        REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
        REQUIRE( moris::equal_to( C.n_cols(), 2 ) );

        REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
        REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 0, 1 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
        REQUIRE( moris::equal_to( C( 2, 1 ), 0 ) );
    }

    SECTION( "scalar != moris::Col" )
    {
        Matrix< DDRMat > A( 3, 1 );
        moris::real B;

        A( 0, 0 ) = 1.0;
        A( 1, 0 ) = 2.0;
        A( 2, 0 ) = 3.0;

        B = 2.0;

        Matrix< DDBMat > C = ( B != A );

        REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
        REQUIRE( moris::equal_to( C.n_cols(), 1 ) );

        REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
        REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
    }

    SECTION( "scalar != moris::Mat " )
    {
        Matrix< DDRMat > A( 3, 2 );
        moris::real B;

        A( 0, 0 ) = 1.0;    A( 0, 1 ) = 1.0;
        A( 1, 0 ) = 2.0;    A( 1, 1 ) = 0.0;
        A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

        B = 2.0;

        Matrix< DDBMat > C = ( B != A );

        REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
        REQUIRE( moris::equal_to( C.n_cols(), 2 ) );

        REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
        REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 0, 1 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
        REQUIRE( moris::equal_to( C( 2, 1 ), 0 ) );
    }
}
}

