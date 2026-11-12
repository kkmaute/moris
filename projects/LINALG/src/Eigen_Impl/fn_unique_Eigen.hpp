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
    template< class E, class M >
    void unique( const Eigen::MatrixBase< E >& a, Matrix< M >& result )
    {
        std::vector< typename E::Scalar > values;
        values.reserve( a.size() );
        for ( Eigen::Index j = 0; j < a.cols(); ++j )
        {
            for ( Eigen::Index i = 0; i < a.rows(); ++i )
            {
                values.push_back( a( i, j ) );
            }
        }
        std::sort( values.begin(), values.end(), linalg_internal::eigen_less< typename E::Scalar > );
        values.erase( std::unique( values.begin(), values.end() ), values.end() );
        result.set_size( a.rows() == 1 ? 1 : values.size(), a.rows() == 1 ? values.size() : 1 );
        std::copy( values.begin(), values.end(), result.data() );
    }
}    // namespace moris
