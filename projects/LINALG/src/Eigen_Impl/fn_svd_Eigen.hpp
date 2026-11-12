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
    template< class U, class S, class V, class E >
    void svd( Matrix< U >& u, Matrix< S >& s, Matrix< V >& v, const Eigen::MatrixBase< E >& a )
    {
        // The Armadillo wrapper retains full U but trims V to the singular-value count.
        Eigen::JacobiSVD< linalg_internal::eigen_plain< E > > solver( a, Eigen::ComputeFullU | Eigen::ComputeThinV );
        if ( solver.info() != Eigen::Success )
        {
            throw std::runtime_error( "svd failed" );
        }
        u = solver.matrixU();
        s = solver.singularValues();
        v = solver.matrixV();
    }

    template< class S, class E >
    void svd( Matrix< S >& s, const Eigen::MatrixBase< E >& a )
    {
        Eigen::JacobiSVD< linalg_internal::eigen_plain< E > > solver( a );
        if ( solver.info() != Eigen::Success )
        {
            throw std::runtime_error( "svd failed" );
        }
        s = solver.singularValues();
    }
}    // namespace moris
