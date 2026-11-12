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
    auto cond( const Eigen::MatrixBase< E >& a )
    {
        using Real = typename Eigen::NumTraits< typename E::Scalar >::Real;
        Eigen::JacobiSVD< linalg_internal::eigen_plain< E > > solver( a );
        if ( solver.info() != Eigen::Success )
        {
            throw std::runtime_error( "cond: SVD failed" );
        }
        const auto& s = solver.singularValues();
        return s.size() && s( s.size() - 1 ) > 0 ? Real( s( 0 ) / s( s.size() - 1 ) ) : std::numeric_limits< Real >::infinity();
    }
}    // namespace moris
