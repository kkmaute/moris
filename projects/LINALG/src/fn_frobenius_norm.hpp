/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_vectorize.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace moris
{
    /**
     * Compute the Frobenius norm of a real matrix without an SVD.
     * Scaling before squaring avoids intermediate overflow and underflow.
     * Empty and zero matrices return zero; nonfinite entries return infinity.
     */
    inline real frobenius_norm( const Matrix< DDRMat >& aMatrix )
    {
        // Vector 2-norm equals matrix Frobenius norm and avoids an SVD.
        Matrix< DDRMat > tVector = vectorize( aMatrix );
        real             tScale  = 0.0;
        for ( size_t i = 0; i < tVector.numel(); ++i )
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
        for ( size_t i = 0; i < tVector.numel(); ++i )
        {
            tVector( i ) = tVector( i ) / tScale;
        }
        return tScale * norm( tVector );
    }

}    // namespace moris
