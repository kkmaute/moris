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

namespace moris
{
template< typename T1, typename T2, typename ET >
auto
dot( Matrix<T2> const & aA,
     const Eigen::MatrixBase<ET> &  aB)
->decltype(aA.matrix_data().cwiseProduct(aB).sum())
{
    return aA.matrix_data().cwiseProduct(aB).sum();
}

template< typename T1, typename T2, typename ET >
auto
dot( const Eigen::MatrixBase<ET> &  aA,
     Matrix<T2> const & aB)
->decltype(aA.cwiseProduct(aB.matrix_data()).sum())
{
    return aA.cwiseProduct(aB.matrix_data()).sum();
}

template< typename ET >
auto
dot( const Eigen::MatrixBase<ET> &  aA,
     const Eigen::MatrixBase<ET> & aB)
->decltype(aA.cwiseProduct(aB).sum())
{
    return aA.cwiseProduct(aB).sum();
}

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DOT_EIGEN_HPP_ */

