/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace moris::fem::test
{
    inline real frobenius_norm( const Matrix< DDRMat >& aMatrix )
    {
        // Vector 2-norm equals matrix Frobenius norm and avoids an SVD.
        Matrix< DDRMat > tVector = vectorize( aMatrix );
        real             tScale  = 0.0;
        for ( uint i = 0; i < tVector.numel(); ++i )
        {
            const real tMagnitude = std::abs( tVector( i ) );
            if ( !std::isfinite( tMagnitude ) )
            {
                return std::numeric_limits< real >::infinity();
            }
            tScale = std::max( tScale, tMagnitude );
        }
        if ( tScale == 0.0 )
        {
            return 0.0;
        }

        // Normalize before squaring. Explicit scalar division also handles
        // subnormal scales whose reciprocal would overflow in a vector kernel.
        for ( uint i = 0; i < tVector.numel(); ++i )
        {
            tVector( i ) = tVector( i ) / tScale;
        }
        return tScale * norm( tVector );
    }

    inline real matrix_difference_bound(
            const Matrix< DDRMat >& aActual,
            const Matrix< DDRMat >& aReference,
            real                    aReferenceScale = 1.0,
            real                    aAbsoluteFloor  = 1.0 )
    {
        REQUIRE( aActual.n_rows() == aReference.n_rows() );
        REQUIRE( aActual.n_cols() == aReference.n_cols() );

        const Matrix< DDRMat > tDifference = aActual - aReferenceScale * aReference;
        const real             tError      = frobenius_norm( tDifference );
        const uint             tRankBound  = std::min( aReference.n_rows(), aReference.n_cols() );

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

    inline real matrix_symmetry_bound( const Matrix< DDRMat >& aMatrix )
    {
        const Matrix< DDRMat > tTranspose = trans( aMatrix );
        return matrix_difference_bound( aMatrix, tTranspose, 1.0, 0.0 );
    }
}    // namespace moris::fem::test
