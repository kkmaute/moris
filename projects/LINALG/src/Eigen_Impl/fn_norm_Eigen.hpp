/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_norm_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_NORM_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_NORM_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template<typename ET >
auto
norm( const Eigen::MatrixBase<ET> &  A)
->decltype( A.norm() )
{
    return  A.norm();
}

template<typename ET >
auto
norm( Eigen::MatrixBase<ET> &  A)
->decltype( A.norm() )
{
    return  A.norm();
}

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_NORM_EIGEN_HPP_ */

