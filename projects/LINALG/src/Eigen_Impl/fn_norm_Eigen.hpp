/*
 * fn_norm_Eigen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
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
