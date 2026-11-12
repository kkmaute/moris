/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"

namespace moris::linalg_internal
{
    template< class E >
    struct eigen_diagonal
    {
        typename Eigen::internal::nested_eval< E, 1 >::type values;
        Eigen::Index                                        offset;
        bool                                                vector;
        typename E::Scalar                                  operator()( Eigen::Index i, Eigen::Index j ) const
        {
            if ( j - i != offset )
            {
                return typename E::Scalar( 0 );
            }
            return vector ? values.coeff( offset >= 0 ? i : j ) : values.coeff( i, j );
        }
    };
}    // namespace moris::linalg_internal
namespace moris
{
    template< class E >
    auto diag_mat( const Eigen::MatrixBase< E >& a, size_t k = 0 )
    {
        const auto offset = static_cast< Eigen::Index >( k );
        const bool vector = a.rows() == 1 || a.cols() == 1;
        const auto n      = a.size() + std::abs( offset );
        return linalg_internal::eigen_plain< E >::NullaryExpr( vector ? n : a.rows(), vector ? n : a.cols(), linalg_internal::eigen_diagonal< E >{ a.derived(), offset, vector } );
    }
}    // namespace moris
