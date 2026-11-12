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
    template< class E >
    size_t rank( const Eigen::MatrixBase< E >& a )
    {
        using Real = typename Eigen::NumTraits< typename E::Scalar >::Real;
        if ( !a.size() )
        {
            return 0;
        }
        Eigen::JacobiSVD< linalg_internal::eigen_plain< E > > svd( a );
        if ( svd.info() != Eigen::Success )
        {
            throw std::runtime_error( "rank: SVD failed" );
        }
        const Real tolerance = std::max( a.rows(), a.cols() ) * svd.singularValues()( 0 ) * std::numeric_limits< Real >::epsilon();
        return ( svd.singularValues().array() > tolerance ).count();
    }
}    // namespace moris
