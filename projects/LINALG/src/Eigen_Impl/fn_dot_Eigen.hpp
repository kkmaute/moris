/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_dot_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DOT_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DOT_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include <Eigen/Dense>
#include "fn_vectorize_Eigen.hpp"

namespace moris
{
    // Match Armadillo: dot pairs column-major entries, independent of shape.
    // Reduce the lazy vectorizations immediately while both operands are alive.
    template< typename T2, typename ET >
    auto
    dot( Matrix< T2 > const               &aA,
            const Eigen::MatrixBase< ET > &aB )
            -> decltype( moris::vectorize( aA.matrix_data() ).cwiseProduct( moris::vectorize( aB ) ).sum() )
    {
        return moris::vectorize( aA.matrix_data() ).cwiseProduct( moris::vectorize( aB ) ).sum();
    }

    template< typename T2, typename ET >
    auto
    dot( const Eigen::MatrixBase< ET > &aA,
            Matrix< T2 > const         &aB )
            -> decltype( moris::vectorize( aA ).cwiseProduct( moris::vectorize( aB.matrix_data() ) ).sum() )
    {
        return moris::vectorize( aA ).cwiseProduct( moris::vectorize( aB.matrix_data() ) ).sum();
    }

    template< typename ETA, typename ETB >
    auto
    dot( const Eigen::MatrixBase< ETA >    &aA,
            const Eigen::MatrixBase< ETB > &aB )
            -> decltype( moris::vectorize( aA ).cwiseProduct( moris::vectorize( aB ) ).sum() )
    {
        return moris::vectorize( aA ).cwiseProduct( moris::vectorize( aB ) ).sum();
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DOT_EIGEN_HPP_ */
