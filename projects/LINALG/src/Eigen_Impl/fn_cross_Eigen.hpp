/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"

namespace moris
{
    template< class A, class B >
    auto cross( const Eigen::MatrixBase< A >& a, const Eigen::MatrixBase< B >& b )
    {
        MORIS_ERROR( a.size() == 3 && b.size() == 3, "cross requires three-component vectors." );
        typename Eigen::internal::nested_eval< A, 1 >::type tA = a.derived();
        typename Eigen::internal::nested_eval< B, 1 >::type tB = b.derived();

        // Armadillo preserves the first vector's orientation and uses the
        // algebraic cross product for complex values, without conjugation.
        Eigen::Matrix< typename A::Scalar, Eigen::Dynamic, Eigen::Dynamic > result( a.rows(), a.cols() );
        for ( Eigen::Index i = 0; i < 3; ++i )
        {
            const auto j = ( i + 1 ) % 3;
            const auto k = ( i + 2 ) % 3;
            result( i )  = tA.coeff( j ) * tB.coeff( k )
                        - tA.coeff( k ) * tB.coeff( j );
        }
        return result;
    }

    template< class M, class E >
    auto cross( const Matrix< M >& a, const Eigen::MatrixBase< E >& b )
    {
        return cross( a.matrix_data(), b );
    }

    template< class M, class E >
    auto cross( const Eigen::MatrixBase< E >& a, const Matrix< M >& b )
    {
        return cross( a, b.matrix_data() );
    }
}    // namespace moris
