/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "fn_frobenius_norm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace moris
{
    // Conservative, SVD-free comparison for real matrices of matching dimensions.
    // Passing a tolerance guarantees the spectral-norm comparison also passes;
    // the bound can reject a case that the spectral-norm comparison accepts.
    // Nonfinite errors or normalization scales return infinity. A zero error
    // with a zero denominator returns zero; a nonzero error returns infinity.
    // Bound ||actual - referenceScale * reference||_2 /
    //       (absoluteFloor + ||reference||_2).
    // referenceScale models an expected physical multiplier (e.g. thickness);
    // the denominator deliberately uses the unscaled reference. The default
    // floor keeps comparisons meaningful near zero; zero requests a relative check.
    inline real matrix_difference_bound(
            const Matrix< DDRMat >& aActual,
            const Matrix< DDRMat >& aReference,
            real                    aReferenceScale = 1.0,
            real                    aAbsoluteFloor  = 1.0 )
    {
        MORIS_ASSERT( aActual.n_rows() == aReference.n_rows() && aActual.n_cols() == aReference.n_cols(),
                "matrix_difference_bound() requires matrices with matching dimensions." );
        MORIS_ASSERT( std::isfinite( aAbsoluteFloor ) && aAbsoluteFloor >= 0.0,
                "matrix_difference_bound() requires a finite, nonnegative absolute floor." );

        const Matrix< DDRMat > tDifference = aActual - aReferenceScale * aReference;
        const real             tError      = frobenius_norm( tDifference );
        const size_t           tRankBound  = std::min( aReference.n_rows(), aReference.n_cols() );

        // ||E||_2 <= ||E||_F and ||A||_F/sqrt(min(m,n)) <= ||A||_2.
        // This ratio therefore bounds the original spectral-norm comparison
        // from above, preserving its tolerance without weakening acceptance.
        const real tReferenceBound = tRankBound == 0 ? 0.0
                                                     : frobenius_norm( aReference ) / std::sqrt( static_cast< real >( tRankBound ) );
        const real tScale          = aAbsoluteFloor + tReferenceBound;
        if ( !std::isfinite( tError ) || !std::isfinite( tScale ) )
        {
            return std::numeric_limits< real >::infinity();
        }
        return tScale > 0.0 ? tError / tScale
                            : ( tError == 0.0 ? 0.0 : std::numeric_limits< real >::infinity() );
    }

}    // namespace moris
