/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_det.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_det.hpp"

TEST_CASE("moris::det",
          "[linalgebra],[det]")
{
    moris::Matrix< moris::DDRMat > aMatr ( 3, 3 );

    aMatr( 0, 0 ) = 2.33; aMatr( 0, 1 ) = 0.00; aMatr( 0, 2 ) = 0.00;
    aMatr( 1, 0 ) = 0.00; aMatr( 1, 1 ) = 2.33; aMatr( 1, 2 ) = 0.00;
    aMatr( 2, 0 ) = 0.00; aMatr( 2, 1 ) = 0.00; aMatr( 2, 2 ) = 2.33;

    moris::real detr = moris::det( aMatr );
    SECTION( "moris::det real" )
    {
        REQUIRE( moris::equal_to( detr, 1.2649337e+01 ) );
    }

}

