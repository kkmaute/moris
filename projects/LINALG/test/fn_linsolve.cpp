/*
 * fn_linsolve.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

// Third-party header files.
#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_linsolve.hpp"

namespace moris
{
TEST_CASE( "moris::linsolve", "[linalgebra],[linsolve]" )
    {
    SECTION("linsolve_real" )
    {
        Matrix< DDRMat > tA( 3, 3 );
        Matrix< DDRMat > tB( 3, 1 );

        tA( 0, 0 ) = 3.0; tA( 0, 1 ) = 0.0; tA( 0, 2 ) = 3.0;
        tA( 1, 0 ) = 4.0; tA( 1, 1 ) = 5.0; tA( 1, 2 ) = 6.0;
        tA( 2, 0 ) = 9.0; tA( 2, 1 ) = 1.0; tA( 2, 2 ) = 9.0;

        tB( 0, 0 ) = 0.0;
        tB( 1, 0 ) = 1.0;
        tB( 2, 0 ) = 2.0;

        Matrix< DDRMat > tC = moris::solve( tA, tB );

        REQUIRE( moris::equal_to ( tC( 0,0 ),  4.5 ) );
        REQUIRE( moris::equal_to ( tC( 1,0 ),  2.0 ) );
        REQUIRE( moris::equal_to ( tC( 2,0 ), -4.5 ) );
    }

//    SECTION( "linsolve_complexe" )
//    {
//        moris::Matrix< moris::DDCMat > A( 2, 2 );
//        moris::Matrix< moris::DDCMat > B( 2, 1 );
//
//        A( 0, 0 ) = { 1.0, 1.0 }; A( 0, 1 ) = { 2.0, -1.0 };
//        A( 1, 0 ) = { 7.0, 0.0 }; A( 1, 1 ) = { 8.0, -2.0 };
//
//        B( 0, 0 ) = { 2.0,  7.0 };
//        B( 1, 0 ) = { 4.0, -9.0 };
//
//        moris::Matrix< moris::DDCMat > C = moris::solve( A,B );
//
//        REQUIRE( moris::equal_to( C( 0,0 ).real(),  4.529729729729730 ) );
//        REQUIRE( moris::equal_to( C( 0,0 ).imag(), -3.778378378378378 ) );
//        REQUIRE( moris::equal_to( C( 1,0 ).real(), -3.772972972972973 ) );
//        REQUIRE( moris::equal_to( C( 1,0 ).imag(),  1.237837837837838 ) );
//    }
    }
}
