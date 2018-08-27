/*
 * op_times.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_TIMES_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_TIMES_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"


namespace moris
{
template< typename T1, typename T2, typename ET >
auto
operator*( const Eigen::MatrixBase<ET> &  aA,
           Mat_New< T1, T2 > & aB )
->decltype( aA * aB.matrix_data() )
{
    return  aA * aB.matrix_data();
}

template< typename T1, typename T2, typename ET >
auto
operator*( Mat_New< T1, T2 > & aA,
           const Eigen::MatrixBase<ET> &  aB)
->decltype( aA.matrix_data() * aB )
{
    return  aA.matrix_data() * aB;
}

}



#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_TIMES_EIGEN_HPP_ */
