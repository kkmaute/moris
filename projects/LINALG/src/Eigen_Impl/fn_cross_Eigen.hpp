/*
 * fn_cross_Eigen.hpp
 *
 *  Created on: Jan 16, 2019
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CROSS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CROSS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

//FIXME: Verify this builds once eigen cmake issue is fixed.
namespace moris
{
template<typename M1, typename ET1 >
auto
cross( const Matrix<M1> & aA,
       const Eigen::MatrixBase<ET1> & aB)
->decltype( aA.matrix_data().cross(aB) )
{
    return aA.cross(aB);
}

template<typename M1, typename ET1 >
auto
cross( const Eigen::MatrixBase<ET1> & aA,
       const Matrix<M1> & aB)
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
