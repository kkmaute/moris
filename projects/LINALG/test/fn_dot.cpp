/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_dot.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_dot.hpp"
#include "fn_trans.hpp"

// calculates the dot product of two vectors

TEST_CASE( "moris::dot",
        "[linalgebra],[dot]" )
{
    moris::Matrix< moris::DDRMat > A( 3, 1 );

    A( 0, 0 ) = 1.0;
    A( 1, 0 ) = 3.0;
    A( 2, 0 ) = -5.0;

    moris::Matrix< moris::DDRMat > B( 3, 1 );

    B( 0, 0 ) = 4;
    B( 1, 0 ) = -2;
    B( 2, 0 ) = -1;

    moris::real C = dot( A, B );

    REQUIRE( moris::equal_to( C, 3 ) );
}


TEST_CASE( "moris::dot ignores vector orientation and matrix shape", "[linalgebra],[dot]" )
{
    using namespace moris;
    const Matrix< DDRMat > tColumn = { { 1.0 }, { 3.0 }, { -5.0 } };
    const Matrix< DDRMat > tRow    = { { 4.0, -2.0, -1.0 } };

    CHECK( moris::dot( tColumn, tRow ) == Approx( 3.0 ) );
    CHECK( moris::dot( tRow, tColumn ) == Approx( 3.0 ) );
    CHECK( moris::dot( tColumn, 2.0 * tRow ) == Approx( 6.0 ) );
    CHECK( moris::dot( 2.0 * tRow, tColumn ) == Approx( 6.0 ) );
    CHECK( moris::dot( 2.0 * tColumn, tRow + tRow ) == Approx( 12.0 ) );
    CHECK( moris::dot( trans( tColumn ), trans( tRow ) ) == Approx( 3.0 ) );

    // Flatten in column-major order even when the shapes differ.
    const Matrix< DDRMat > tMatrix = { { 1.0, 3.0, 5.0 }, { 2.0, 4.0, 6.0 } };
    const Matrix< DDRMat > tOther  = { { 1.0, 4.0 }, { 2.0, 5.0 }, { 3.0, 6.0 } };
    CHECK( moris::dot( tMatrix, tOther ) == Approx( 91.0 ) );
    CHECK( moris::dot( trans( tMatrix ), tOther ) == Approx( 86.0 ) );

    // dot is bilinear; it does not conjugate complex entries.
    Matrix< DDCMat > tComplexColumn( 2, 1 );
    Matrix< DDCMat > tComplexRow( 1, 2 );
    tComplexColumn( 0 ) = cplx( 1.0, 2.0 );
    tComplexColumn( 1 ) = cplx( 3.0, -4.0 );
    tComplexRow( 0 )    = cplx( 2.0, 1.0 );
    tComplexRow( 1 )    = cplx( 4.0, 3.0 );
    const cplx tDot     = moris::dot( tComplexColumn, tComplexRow );
    CHECK( tDot.real() == Approx( 24.0 ) );
    CHECK( tDot.imag() == Approx( -2.0 ) );
}
