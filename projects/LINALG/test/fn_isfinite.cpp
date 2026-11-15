/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isfinite.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_isfinite.hpp"
#include <limits>

namespace moris
{
    TEST_CASE( "moris::isfinite", "[linalgebra],[isfinite]" )
    {
        moris::real      tZero = 0.0;
        Matrix< DDRMat > a( 3, 3 );
        Matrix< DDRMat > b( 3, 3 );

        a( 0, 0 ) = 1.0;
        a( 0, 1 ) = 2.0;
        a( 0, 2 ) = 3.0;
        a( 1, 0 ) = 0.0;
        a( 1, 1 ) = 1.0;
        a( 1, 2 ) = 4.0;
        a( 2, 0 ) = 5.0;
        a( 2, 1 ) = 6.0;
        a( 2, 2 ) = 0.0;

        b( 0, 0 ) = 1.0;
        b( 0, 1 ) = 2.0;
        b( 0, 2 ) = 3.0;
        b( 1, 0 ) = 0.0;
        b( 1, 1 ) = 1.0;
        b( 1, 2 ) = 4.0;
        b( 2, 0 ) = 5.0;
        b( 2, 1 ) = 6.0;
        b( 2, 2 ) = 1.0 / tZero;

        bool tIsFinite_1 = isfinite( a );
        bool tIsFinite_2 = isfinite( b );

        REQUIRE( tIsFinite_1 );
        REQUIRE( !tIsFinite_2 );
    }
}    // namespace moris

TEST_CASE( "LINALG finite checks handle extremes, views and expressions", "[linalgebra],[isfinite]" )
{
    using namespace moris;
    const real       tInfinity = std::numeric_limits< real >::infinity();
    const real       tNaN      = std::numeric_limits< real >::quiet_NaN();
    Matrix< DDRMat > tValues( 3, 3, 1.0 );
    CHECK( moris::isfinite( tValues ) );
    CHECK( moris::isfinite( Matrix< DDRMat >( 0, 3 ) ) );
    tValues( 2, 2 ) = tInfinity;
    CHECK_FALSE( moris::isfinite( tValues ) );
    CHECK( moris::isfinite( tValues( { 0, 1 }, { 0, 1 } ) ) );
    tValues( 0, 1 ) = -tInfinity;
    CHECK_FALSE( moris::isfinite( tValues( { 0, 1 }, { 0, 1 } ) ) );
    tValues.fill( tNaN );
    CHECK_FALSE( moris::isfinite( tValues ) );
    tValues.fill( std::numeric_limits< real >::max() );
    CHECK( moris::isfinite( tValues ) );
    CHECK_FALSE( moris::isfinite( tValues + tValues ) );
    CHECK( moris::isfinite( tValues - tValues ) );

    Matrix< DDCMat > tComplex( 2, 1, cplx( 1.0, 2.0 ) );
    CHECK( moris::isfinite( tComplex ) );
    tComplex( 0 ) = cplx( 1.0, tNaN );
    CHECK_FALSE( moris::isfinite( tComplex ) );
    tComplex( 0 ) = cplx( tInfinity, 0.0 );
    CHECK_FALSE( moris::isfinite( tComplex ) );
}
