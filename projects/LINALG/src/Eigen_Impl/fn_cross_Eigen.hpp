/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cross_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CROSS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CROSS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template<typename M1, typename ET1 >
auto
cross( const Matrix<M1> & aA,
       const Eigen::MatrixBase<ET1> & aB)
->decltype( aA.matrix_data().cross(aB) )
{
    return aA.matrix_data().cross(aB);
}

template<typename M1, typename ET1 >
auto
cross( const Eigen::MatrixBase<ET1> & aA,
        Matrix<M1> & aB)
->decltype( aA.cross(aB.matrix_data()) )
{
    return aA.cross(aB.matrix_data());
}

template<typename ET1, typename ET2 >
auto
cross( const Eigen::MatrixBase<ET1> & aA,
       const Eigen::MatrixBase<ET2> & aB)
->decltype( aA.cross(aB) )
{
    return aA.cross(aB);
}
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CROSS_EIGEN_HPP_ */

