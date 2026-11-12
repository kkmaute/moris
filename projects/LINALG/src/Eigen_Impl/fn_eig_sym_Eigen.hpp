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
    template< class V, class W, class E >
    void eig_sym( Matrix< V >& values, Matrix< W >& vectors, const Eigen::MatrixBase< E >& a )
    {
        MORIS_ERROR( a.rows() == a.cols(), "eig_sym requires a square matrix." );
        Eigen::SelfAdjointEigenSolver< linalg_internal::eigen_plain< E > > solver( a );
        if ( solver.info() != Eigen::Success )
        {
            throw std::runtime_error( "eig_sym failed" );
        }
        values  = solver.eigenvalues();
        vectors = solver.eigenvectors();
    }
}    // namespace moris
