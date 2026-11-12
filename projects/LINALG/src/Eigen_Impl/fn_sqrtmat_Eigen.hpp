/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"
#include <unsupported/Eigen/MatrixFunctions>

namespace moris
{
    template< class E >
    linalg_internal::eigen_plain< E > sqrtmat( const Eigen::MatrixBase< E >& a )
    {
        MORIS_ERROR( a.rows() == a.cols(), "sqrtmat requires a square matrix." );
        using Real                                                      = typename Eigen::NumTraits< typename E::Scalar >::Real;
        using Complex                                                   = std::complex< Real >;
        Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic > result = a.template cast< Complex >().eval().sqrt();
        if constexpr ( Eigen::NumTraits< typename E::Scalar >::IsComplex )
        {
            return result;
        }
        else
        {
            MORIS_ERROR( result.imag().norm() <= 1e-12 * std::max< Real >( 1, a.norm() ), "sqrtmat: square root is not real." );
            return result.real();
        }
    }

    template< class E, class M >
    void sqrtmat( const Eigen::MatrixBase< E >& a, Matrix< M >& x )
    {
        x = sqrtmat( a );
    }
}    // namespace moris
