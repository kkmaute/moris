/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "fn_Eigen_Utilities.hpp"

namespace moris::linalg_internal
{
    template< typename E >
    struct eigen_reshape
    {
        typename Eigen::internal::nested_eval< E, 1 >::type values;
        Eigen::Index                                        rows;
        bool                                                rowOrder;

        typename E::Scalar operator()( Eigen::Index i, Eigen::Index j ) const
        {
            const auto tIndex = i + j * rows;
            if ( tIndex >= values.size() )
            {
                return typename E::Scalar( 0 );
            }
            return rowOrder ? values.coeff( tIndex / values.cols(), tIndex % values.cols() )
                            : values.coeff( tIndex % values.rows(), tIndex / values.rows() );
        }
    };
}    // namespace moris::linalg_internal

namespace moris
{
    template< typename E >
    auto reshape( const Eigen::MatrixBase< E >& aA, size_t aRows, size_t aCols )
    {
        // Retain the operand expression; padding/truncation follows Armadillo.
        return linalg_internal::eigen_plain< E >::NullaryExpr( aRows, aCols, linalg_internal::eigen_reshape< E >{ aA.derived(), static_cast< Eigen::Index >( aRows ), false } );
    }
}    // namespace moris
