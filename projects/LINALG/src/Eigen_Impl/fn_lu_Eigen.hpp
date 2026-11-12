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
    template< class L, class U, class P, class E >
    void lu( Matrix< L >& l, Matrix< U >& u, Matrix< P >& p, const Eigen::MatrixBase< E >& a )
    {
        using Plain = linalg_internal::eigen_plain< E >;
        if ( a.rows() == a.cols() && a.rows() )
        {
            Eigen::PartialPivLU< Plain > tFactor( a );
            Plain                        tLower = tFactor.matrixLU().template triangularView< Eigen::StrictlyLower >();
            tLower.diagonal().setOnes();
            l = tLower;
            u = tFactor.matrixLU().template triangularView< Eigen::Upper >();
            p = tFactor.permutationP().toDenseMatrix().template cast< typename E::Scalar >();
            return;
        }

        // Explicit partial pivoting supports rectangular matrices, unlike PartialPivLU.
        Plain      work        = a;
        Plain      lower       = Plain::Identity( a.rows(), std::min( a.rows(), a.cols() ) );
        Plain      permutation = Plain::Identity( a.rows(), a.rows() );
        const auto kmax        = std::min( a.rows(), a.cols() );
        for ( Eigen::Index k = 0; k < kmax; ++k )
        {
            Eigen::Index pivot = k;
            for ( Eigen::Index i = k + 1; i < a.rows(); ++i )
            {
                if ( std::abs( work( i, k ) ) > std::abs( work( pivot, k ) ) )
                {
                    pivot = i;
                }
            }
            if ( pivot != k )
            {
                work.row( k ).swap( work.row( pivot ) );
                permutation.row( k ).swap( permutation.row( pivot ) );
                if ( k )
                {
                    lower.row( k ).head( k ).swap( lower.row( pivot ).head( k ) );
                }
            }
            if ( std::abs( work( k, k ) ) == 0 )
            {
                continue;
            }
            for ( Eigen::Index i = k + 1; i < a.rows(); ++i )
            {
                lower( i, k ) = work( i, k ) / work( k, k );
                work.row( i ).segment( k, a.cols() - k ) -= lower( i, k ) * work.row( k ).segment( k, a.cols() - k );
            }
        }
        l = lower;
        u = work.topRows( kmax ).template triangularView< Eigen::Upper >();
        p = permutation;
    }
}    // namespace moris
