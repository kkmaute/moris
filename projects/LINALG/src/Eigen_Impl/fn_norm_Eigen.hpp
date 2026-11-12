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
    typename Eigen::NumTraits< typename E::Scalar >::Real norm( const Eigen::MatrixBase< E >& a )
    {
        if ( a.size() == 0 )
        {
            return 0;
        }
        if ( a.rows() == 1 || a.cols() == 1 )
        {
            return a.norm();
        }
        Eigen::JacobiSVD< linalg_internal::eigen_plain< E > > svd( a );
        if ( svd.info() != Eigen::Success )
        {
            throw std::runtime_error( "norm: SVD failed" );
        }
        return svd.singularValues()( 0 );
    }
}    // namespace moris
