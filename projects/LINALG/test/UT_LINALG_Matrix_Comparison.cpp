/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include "catch.hpp"
#include "fn_frobenius_norm.hpp"
#include "fn_matrix_difference_bound.hpp"
#include "fn_matrix_symmetry_bound.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include <cmath>
#include <limits>

namespace moris
{
    // Validate the comparison helpers themselves, independently of FEM assembly.
    // A bound must never understate the original spectral-norm error: otherwise
    // replacing an expensive SVD-based check could hide a numerical regression.
    TEST_CASE( "LINALG matrix comparison bounds preserve spectral norm tolerances", "[linalgebra],[Matrix_Comparison]" )
    {
        const Matrix< DDRMat > tReference = { { 1.0, 2.0 }, { -3.0, 4.0 }, { 5.0, -6.0 } };
        const Matrix< DDRMat > tError     = { { 0.25, -0.5 }, { 0.0, 0.75 }, { -1.0, 0.125 } };
        // Exercise absent, fractional, and amplified reference contributions.
        for ( real tScale : { 0.0, 0.4, 2.5 } )
        {
            const Matrix< DDRMat > tActual        = tScale * tReference + tError;
            const real             tSpectralError = norm( Matrix< DDRMat >( tActual - tScale * tReference ) ) / ( 1.0 + norm( tReference ) );
            CHECK( matrix_difference_bound( tActual, tReference, tScale ) >= tSpectralError );
            // Materialization rounds the product before a potentially fused
            // multiply/subtract in the comparison, leaving roundoff-sized error.
            CHECK( matrix_difference_bound( Matrix< DDRMat >( tScale * tReference ), tReference, tScale ) < 1e-15 );
        }

        // For a vector the Frobenius and spectral norms coincide, so the
        // conservative matrix bound must reduce to the original comparison.
        const Matrix< DDRMat > tVector       = { { 1.0 }, { -2.0 }, { 3.0 } };
        const Matrix< DDRMat > tActualVector = { { 1.25 }, { -2.0 }, { 3.0 } };
        CHECK( matrix_difference_bound( tActualVector, tVector )
                == Approx( norm( Matrix< DDRMat >( tActualVector - tVector ) ) / ( 1.0 + norm( tVector ) ) ) );

        // Accept exact symmetry and detect a small, deliberately introduced
        // off-diagonal mismatch without understating its relative size.
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

        // Zero references must not cause division by zero or accept a nonzero
        // discrepancy; the relative symmetry check must accept the zero matrix.
        const Matrix< DDRMat > tZero( 2, 2, 0.0 );
        CHECK( matrix_difference_bound( tZero, tZero ) == 0.0 );
        CHECK( matrix_symmetry_bound( tZero ) == 0.0 );
        CHECK( matrix_difference_bound( tSymmetric, tZero ) > 1.0 );
    }

    TEST_CASE( "LINALG Frobenius norm handles extreme and nonfinite entries", "[linalgebra],[Matrix_Comparison]" )
    {
        // Rectangular data distinguish the Frobenius norm from the spectral norm.
        const Matrix< DDRMat > tMatrix = { { 3.0, 0.0, 0.0 }, { 0.0, 4.0, 0.0 } };
        CHECK( frobenius_norm( tMatrix ) == Approx( 5.0 ) );
        for ( real tScale : { 1e-200, 1e200 } )
        {
            CHECK( frobenius_norm( Matrix< DDRMat >( tScale * tMatrix ) ) / tScale == Approx( 5.0 ) );
        }

        // A reciprocal of a subnormal scale would overflow, even though this
        // norm is representable. Exact equality prevents an absolute tolerance from accepting zero.
        const real tTiny = std::numeric_limits< real >::denorm_min();
        CHECK( frobenius_norm( Matrix< DDRMat >{ { tTiny } } ) == tTiny );

        const Matrix< DDRMat > tEmpty;
        CHECK( frobenius_norm( tEmpty ) == 0.0 );
        CHECK( matrix_difference_bound( tEmpty, tEmpty ) == 0.0 );
        CHECK( matrix_symmetry_bound( tEmpty ) == 0.0 );

        // Nonfinite values must never make an error or symmetry check pass.
        const Matrix< DDRMat > tZero( 1, 1, 0.0 );
        for ( real tValue : { std::numeric_limits< real >::infinity(), std::numeric_limits< real >::quiet_NaN() } )
        {
            const Matrix< DDRMat > tNonfinite = { { tValue } };
            CHECK( std::isinf( frobenius_norm( tNonfinite ) ) );
            CHECK( std::isinf( matrix_difference_bound( tNonfinite, tZero ) ) );
            CHECK( std::isinf( matrix_difference_bound( tZero, tNonfinite ) ) );
            CHECK( std::isinf( matrix_symmetry_bound( tNonfinite ) ) );
        }
    }
}    // namespace moris
