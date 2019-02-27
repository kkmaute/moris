/*
 * op_trans_Eigen.hpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
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
