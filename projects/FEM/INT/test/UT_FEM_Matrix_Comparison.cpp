/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include "fn_FEM_Matrix_Comparison.hpp"

namespace moris::fem::test
{
    TEST_CASE( "FEM matrix comparison bounds preserve spectral norm tolerances", "[fem],[Matrix_Comparison]" )
    {
        const Matrix< DDRMat > tReference = { { 1.0, 2.0 }, { -3.0, 4.0 }, { 5.0, -6.0 } };
        const Matrix< DDRMat > tError     = { { 0.25, -0.5 }, { 0.0, 0.75 }, { -1.0, 0.125 } };
        for ( real tScale : { 0.0, 0.4, 2.5 } )
        {
            const Matrix< DDRMat > tActual        = tScale * tReference + tError;
            const real             tSpectralError = norm( Matrix< DDRMat >( tActual - tScale * tReference ) ) / ( 1.0 + norm( tReference ) );
            CHECK( matrix_difference_bound( tActual, tReference, tScale ) >= tSpectralError );
            // Materialization rounds the product before a potentially fused
            // multiply/subtract in the comparison, leaving roundoff-sized error.
            CHECK( matrix_difference_bound( Matrix< DDRMat >( tScale * tReference ), tReference, tScale ) < 1e-15 );
        }

        const Matrix< DDRMat > tVector       = { { 1.0 }, { -2.0 }, { 3.0 } };
        const Matrix< DDRMat > tActualVector = { { 1.25 }, { -2.0 }, { 3.0 } };
        CHECK( matrix_difference_bound( tActualVector, tVector )
                == Approx( norm( Matrix< DDRMat >( tActualVector - tVector ) ) / ( 1.0 + norm( tVector ) ) ) );

        const Matrix< DDRMat > tSymmetric = { { 2.0, 1.0 }, { 1.0, 3.0 } };
        CHECK( matrix_symmetry_bound( tSymmetric ) == 0.0 );
        Matrix< DDRMat > tAsymmetric = tSymmetric;
        tAsymmetric( 0, 1 ) += 1e-6;
        const real tSpectralAsymmetry = norm( Matrix< DDRMat >( tAsymmetric - trans( tAsymmetric ) ) ) / norm( tAsymmetric );
        CHECK( matrix_symmetry_bound( tAsymmetric ) >= tSpectralAsymmetry );
        CHECK( matrix_symmetry_bound( tAsymmetric ) > 1e-12 );

        // Bounds must remain meaningful when unscaled sums of squares would
        // underflow or overflow, especially for the relative symmetry check.
        for ( real tScale : { 1e-200, 1e200 } )
        {
            const Matrix< DDRMat > tScaled              = tScale * tAsymmetric;
            const real             tScaledSpectralError = norm( Matrix< DDRMat >( tScaled - trans( tScaled ) ) ) / norm( tScaled );
            CHECK( matrix_symmetry_bound( tScaled ) >= tScaledSpectralError );
            CHECK( matrix_symmetry_bound( tScaled ) > 1e-12 );
        }

        const Matrix< DDRMat > tZero( 2, 2, 0.0 );
        CHECK( matrix_difference_bound( tZero, tZero ) == 0.0 );
        CHECK( matrix_symmetry_bound( tZero ) == 0.0 );
        CHECK( matrix_difference_bound( tSymmetric, tZero ) > 1.0 );
    }
}    // namespace moris::fem::test
