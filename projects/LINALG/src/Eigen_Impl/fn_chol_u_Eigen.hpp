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
    linalg_internal::eigen_plain< E > chol_u( const Eigen::MatrixBase< E >& a )
    {
        MORIS_ERROR( a.rows() == a.cols(), "Cholesky requires a square matrix." );
        Eigen::LLT< linalg_internal::eigen_plain< E > > solver( a );
        if ( solver.info() != Eigen::Success )
        {
            throw std::runtime_error( "Cholesky factorization failed" );
        }
        return solver.matrixU();
    }
}    // namespace moris
