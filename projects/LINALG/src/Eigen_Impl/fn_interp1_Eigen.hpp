/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"

namespace moris
{
    template< class X, class Y, class Q, class M >
    void interp1( const Eigen::MatrixBase< X >& aX, const Eigen::MatrixBase< Y >& aY, const Eigen::MatrixBase< Q >& aQueries, Matrix< M >& aResult )
    {
        MORIS_ERROR( aX.size() >= 2 && aX.size() == aY.size(), "interp1 requires at least two paired values." );
        auto tIndices = linalg_internal::eigen_sorted_indices( aX.derived() );
        for ( size_t i = 1; i < tIndices.size(); ++i )
        {
            MORIS_ERROR( aX.derived().coeff( tIndices[ i ] ) > aX.derived().coeff( tIndices[ i - 1 ] ), "interp1 requires distinct abscissae." );
        }
        using Scalar = typename Y::Scalar;
        Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > tValues( aQueries.rows(), aQueries.cols() );
        for ( Eigen::Index j = 0; j < aQueries.cols(); ++j )
        {
            for ( Eigen::Index i = 0; i < aQueries.rows(); ++i )
            {
                const auto tQuery = aQueries( i, j );
                if ( !( tQuery >= aX.derived().coeff( tIndices.front() ) && tQuery <= aX.derived().coeff( tIndices.back() ) ) )
                {
                    tValues( i, j ) = Scalar( std::numeric_limits< typename Eigen::NumTraits< Scalar >::Real >::quiet_NaN() );
                    continue;
                }
                size_t tLow = 0, tHigh = tIndices.size() - 1;
                while ( tLow + 1 < tHigh )
                {
                    const auto tMiddle = ( tLow + tHigh ) / 2;
                    if ( tQuery < aX.derived().coeff( tIndices[ tMiddle ] ) )
                    {
                        tHigh = tMiddle;
                    }
                    else
                    {
                        tLow = tMiddle;
                    }
                }
                const auto tFraction = ( tQuery - aX.derived().coeff( tIndices[ tLow ] ) ) / ( aX.derived().coeff( tIndices[ tHigh ] ) - aX.derived().coeff( tIndices[ tLow ] ) );
                tValues( i, j )      = ( 1 - tFraction ) * aY.derived().coeff( tIndices[ tLow ] ) + tFraction * aY.derived().coeff( tIndices[ tHigh ] );
            }
        }
        aResult = tValues;
    }
}    // namespace moris
