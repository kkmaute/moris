/*
 * op_greater_equal_Eigen.hpp
 *
 *  Created on: Aug 31, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_GREATER_EQUAL_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_GREATER_EQUAL_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template< typename T1, typename T2, typename ET >
auto
operator>=( const Eigen::MatrixBase<ET> &  aA,
           Matrix< T1, T2 > & aB )
->decltype( aA.array() >= aB.matrix_data().array() )
{
    return  aA.array() >= aB.matrix_data().array();
}

template< typename T1, typename T2, typename ET >
auto
operator>=( Matrix< T1, T2 > & aA,
           const Eigen::MatrixBase<ET> &  aB)
->decltype( aA.matrix_data().array() >= aB.array() )
{
    return  aA.matrix_data().array() >= aB.array();
}

template< typename ET1, typename ET2 >
auto
operator>=( const Eigen::MatrixBase<ET1> &  aA,
           const Eigen::MatrixBase<ET2> &  aB)
->decltype( aA.array() >= aB.array() )
{
    return  aA.array() >= aB.array();
}

// scalar to matrix comparisons
template< typename ET1, typename Scalar >
auto
operator>=(Eigen::MatrixBase<ET1> const & aA,
            Scalar const & aB )
->decltype( aA.array() >= aB  )
{
    return aA.array() >=  aB ;
}

template< typename Scalar, typename ET1 >
auto
operator>=(
        Scalar & aA,
        Eigen::MatrixBase<ET1> const & aB)
->decltype( aA >= aB.array() )
{
    return aA >= aB.array();
}
}
#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_GREATER_EQUAL_EIGEN_HPP_ */
