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
    void eig_gen( Matrix< V >& values, Matrix< W >& vectors, const Eigen::MatrixBase< E >& a )
    {
        using C     = std::complex< typename Eigen::NumTraits< typename E::Scalar >::Real >;
        using Plain = Eigen::Matrix< C, Eigen::Dynamic, Eigen::Dynamic >;
        MORIS_ERROR( a.rows() == a.cols(), "eig_gen requires a square matrix." );
        Eigen::ComplexEigenSolver< Plain > solver( a.template cast< C >() );
        if ( solver.info() != Eigen::Success )
        {
            throw std::runtime_error( "eig_gen failed" );
        }
        values  = solver.eigenvalues();
        vectors = solver.eigenvectors();
    }
}    // namespace moris
