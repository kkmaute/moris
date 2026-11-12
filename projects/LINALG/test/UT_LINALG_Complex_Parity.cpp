/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include <catch.hpp>
#include "cl_Matrix.hpp"
#include "fn_chol_l.hpp"
#include "fn_cross.hpp"
#include "fn_ctrans.hpp"
#include "fn_eig_sym.hpp"
#include "fn_linsolve.hpp"
#include "fn_max.hpp"
#include "fn_min.hpp"
#include "fn_norm.hpp"

namespace moris
{
    TEST_CASE( "LINALG complex factorizations preserve imaginary components", "[linalgebra],[backend_parity]" )
    {
        // This Hermitian positive definite matrix requires complex factors.
        const Matrix< DDCMat > tA = { { cplx( 3.0 ), cplx( 0.0, 1.0 ) }, { cplx( 0.0, -1.0 ), cplx( 3.0 ) } };
        const Matrix< DDCMat > tL = chol_l( tA );
        CHECK( norm( tL * ctrans( tL ) - tA ) < 1e-12 );
        // The shared solve interface reserves systems up to 10x10 for inv().
        Matrix< DDCMat > tSystem( 12, 12, cplx( 0.0 ) );
        for ( uint i = 0; i < 12; ++i )
        {
            tSystem( i, i ) = 3.0;
        }
        tSystem( 0, 1 ) = cplx( 0.0, 1.0 );
        tSystem( 1, 0 ) = cplx( 0.0, -1.0 );
        const Matrix< DDCMat > tB( 12, 1, cplx( 1.0, 2.0 ) );
        const Matrix< DDCMat > tX = solve( tSystem, tB );
        CHECK( norm( tSystem * tX - tB ) < 1e-12 );

        Matrix< DDRMat > tValues;
        Matrix< DDCMat > tVectors;
        eig_sym( tValues, tVectors, tA );
        REQUIRE( tValues.numel() == 2 );
        CHECK( tValues( 0 ) == Approx( 2.0 ) );
        CHECK( tValues( 1 ) == Approx( 4.0 ) );
        for ( uint j = 0; j < 2; ++j )
        {
            for ( uint i = 0; i < 2; ++i )
            {
                cplx tResidual = -tValues( j ) * tVectors( i, j );
                for ( uint k = 0; k < 2; ++k )
                {
                    tResidual += tA( i, k ) * tVectors( k, j );
                }
                CHECK( std::abs( tResidual ) < 1e-12 );
            }
        }
    }

    TEST_CASE( "LINALG complex extrema use magnitude and cross preserves orientation", "[linalgebra],[backend_parity]" )
    {
        const Matrix< DDCMat > tA       = { { cplx( 0.0, 2.0 ), cplx( 3.0 ) }, { cplx( 1.0 ), cplx( 0.0, -4.0 ) } };
        Matrix< DDCMat >       tMinimum = min( tA );
        Matrix< DDCMat >       tMaximum = max( tA );
        CHECK( tMinimum( 0 ) == cplx( 1.0 ) );
        CHECK( tMinimum( 1 ) == cplx( 3.0 ) );
        CHECK( tMaximum( 0 ) == cplx( 0.0, 2.0 ) );
        CHECK( tMaximum( 1 ) == cplx( 0.0, -4.0 ) );
        const Matrix< DDCMat > tB    = tA + tA;
        Matrix< DDCMat >       tPair = min( tA, tB );
        CHECK( norm( tPair - tA ) < 1e-12 );

        const Matrix< DDCMat > tX     = { { cplx( 1.0, 1.0 ), cplx( 0.0 ), cplx( 0.0 ) } };
        const Matrix< DDCMat > tY     = { { cplx( 0.0 ), cplx( 1.0 ), cplx( 0.0 ) } };
        Matrix< DDCMat >       tCross = cross( tX, tY );
        REQUIRE( tCross.n_rows() == 1 );
        REQUIRE( tCross.n_cols() == 3 );
        CHECK( tCross( 0 ) == cplx( 0.0 ) );
        CHECK( tCross( 1 ) == cplx( 0.0 ) );
        CHECK( tCross( 2 ) == cplx( 1.0, 1.0 ) );
    }
}    // namespace moris
