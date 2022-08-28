/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_greater.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_greater.hpp"

namespace moris
{
TEST_CASE(
        "moris::op_greater",
        "[linalgebra], [op_greater]" )
{
    SECTION(
                "moris::Col > moris::Col" )
        {
            Matrix< DDRMat > a( 3, 1 );
            Matrix< DDRMat > b( 3, 1 );

            a( 0, 0 ) = 1.0;
            a( 1, 0 ) = 0.0;
            a( 2, 0 ) = 3.0;

            b( 0, 0 ) = 1.0;
            b( 1, 0 ) = 2.0;
            b( 2, 0 ) = 1.0;

            Matrix< DDBMat > c = ( a > b );

            REQUIRE( c( 0, 0 ) == 0 );
            REQUIRE( c( 1, 0 ) == 0 );
            REQUIRE( c( 2, 0 ) == 1 );
        }

        SECTION(
                "moris::Mat > moris::Mat" )
        {
            Matrix< DDRMat > a( 3, 3 );
            Matrix< DDRMat > b( 3, 3 );

            a( 0, 0 ) = 1.0; a( 0, 1 ) = 2.0; a( 0, 2 ) = 3.0;
            a( 1, 0 ) = 4.0; a( 1, 1 ) = 5.0; a( 1, 2 ) = 6.0;
            a( 2, 0 ) = 7.0; a( 2, 1 ) = 9.0; a( 2, 2 ) = 9.0;

            b( 0, 0 ) = 3.0; b( 0, 1 ) = 2.0; b( 0, 2 ) = 1.0;
            b( 1, 0 ) = 6.0; b( 1, 1 ) = 5.0; b( 1, 2 ) = 4.0;
            b( 2, 0 ) = 9.0; b( 2, 1 ) = 4.0; b( 2, 2 ) = 7.0;

            Matrix< DDBMat > c = ( a > b );

            REQUIRE( c( 0, 0 ) == 0 );
            REQUIRE( c( 0, 1 ) == 0 );
            REQUIRE( c( 0, 2 ) == 1 );

            REQUIRE( c( 1, 0 ) == 0 );
            REQUIRE( c( 1, 1 ) == 0 );
            REQUIRE( c( 1, 2 ) == 1 );

            REQUIRE( c( 2, 0 ) == 0 );
            REQUIRE( c( 2, 1 ) == 1 );
            REQUIRE( c( 2, 2 ) == 1 );
        }

    //    greater_than operator is currently not yet working for cplx numbers
    //    MATLAB only compares the real parts of a cplx number
    //
    //    SECTION(
    //            "moris::Mat cplx > moris::Mat cplx" )
    //    {
    //        moris::Mat< moris::cplx > d( 3, 3 );
    //        moris::Mat< moris::cplx > e( 3, 3 );
    //
    //        d( 0, 0 ) = {0.0,1.0}; d( 0, 1 ) = {2.0,3.0}; d( 0, 2 ) = {4.0,5.0};
    //        d( 1, 0 ) = {6.0,7.0}; d( 1, 1 ) = {8.0,9.0}; d( 1, 2 ) = {0.0,1.0};
    //        d( 2, 0 ) = {2.0,3.0}; d( 2, 1 ) = {4.0,5.0}; d( 2, 2 ) = {6.0,7.0};
    //
    //        e( 0, 0 ) = {9.0,8.0}; e( 0, 1 ) = {7.0,6.0}; e( 0, 2 ) = {5.0,4.0};
    //        e( 1, 0 ) = {3.0,2.0}; e( 1, 1 ) = {1.0,0.0}; e( 1, 2 ) = {9.0,8.0};
    //        e( 2, 0 ) = {7.0,6.0}; e( 2, 1 ) = {5.0,4.0}; e( 2, 2 ) = {3.0,2.0};
    //
    //        Matrix< DDBMat > f = ( d > e );
    //
    //        REQUIRE( f( 0, 0 ) == 0 ); REQUIRE( f( 0, 1 ) == 0 ); REQUIRE( f( 0, 2 ) == 1 );
    //        REQUIRE( f( 1, 0 ) == 0 ); REQUIRE( f( 1, 1 ) == 0 ); REQUIRE( f( 1, 2 ) == 1 );
    //        REQUIRE( f( 2, 0 ) == 0 ); REQUIRE( f( 2, 1 ) == 1 ); REQUIRE( f( 2, 2 ) == 1 );
    //    }

        SECTION( "moris::Col > scalar" )
         {
             Matrix< DDRMat > A( 3, 1 );
             moris::real B;

             A( 0, 0 ) = 1.0;
             A( 1, 0 ) = 2.0;
             A( 2, 0 ) = 3.0;

             B = 2.0;

             Matrix< DDBMat > C = ( A > B );

             REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
             REQUIRE( moris::equal_to( C.n_cols(), 1 ) );

             REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
         }

         SECTION( "moris::Mat > scalar" )
         {
             Matrix< DDRMat > A( 3, 2 );
             moris::real B;

             A( 0, 0 ) = 1.0;    A( 0, 1 ) = 4.0;
             A( 1, 0 ) = 2.0;    A( 1, 1 ) = 0.0;
             A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

             B = 2.0;

             Matrix< DDBMat > C = ( A > B );

             REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
             REQUIRE( moris::equal_to( C.n_cols(), 2 ) );

             REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
             REQUIRE( moris::equal_to( C( 0, 1 ), 1 ) );
             REQUIRE( moris::equal_to( C( 1, 1 ), 0 ) );
             REQUIRE( moris::equal_to( C( 2, 1 ), 0 ) );
         }

         SECTION( "scalar > moris::Col" )
         {
             Matrix< DDRMat > A( 3, 1 );
             moris::real B;

             A( 0, 0 ) = 1.0;
             A( 1, 0 ) = 2.0;
             A( 2, 0 ) = 3.0;

             B = 2.0;

             Matrix< DDBMat > C = ( B > A );

             REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
             REQUIRE( moris::equal_to( C.n_cols(), 1 ) );

             REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
             REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
         }

         SECTION( "scalar > moris::Mat " )
         {
             Matrix< DDRMat > A( 3, 2 );
             moris::real B;

             A( 0, 0 ) = 1.0;    A( 0, 1 ) = 1.0;
             A( 1, 0 ) = 2.0;    A( 1, 1 ) = 3.0;
             A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

             B = 2.0;

             Matrix< DDBMat > C = ( B > A );

             REQUIRE( moris::equal_to( C.n_rows(), 3 ) );
             REQUIRE( moris::equal_to( C.n_cols(), 2 ) );

             REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
             REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
             REQUIRE( moris::equal_to( C( 0, 1 ), 1 ) );
             REQUIRE( moris::equal_to( C( 1, 1 ), 0 ) );
             REQUIRE( moris::equal_to( C( 2, 1 ), 0 ) );
         }
}
}

