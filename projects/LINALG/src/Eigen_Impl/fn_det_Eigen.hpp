/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_det_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DET_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DET_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{

template<typename ET >
auto
det( const Eigen::MatrixBase<ET> &  A)
->decltype( A.determinant() )
{
    return A.determinant();
}

template<typename ET >
auto
det( Eigen::MatrixBase<ET> &  A)
->decltype( A.determinant() )
{
    return A.determinant();
}

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DET_EIGEN_HPP_ */

