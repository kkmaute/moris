/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_trans_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_TRANS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_TRANS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
namespace linalg_internal
{
template<typename ET >
auto
trans( const Eigen::MatrixBase<ET> &  A)
->decltype( A.transpose() )
{
    return A.transpose();
}

template<typename ET >
auto
trans( Eigen::MatrixBase<ET> &  A)
->decltype( A.transpose() )
{
    return A.transpose();
}

}}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_TRANS_EIGEN_HPP_ */

