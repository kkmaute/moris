/*
 * op_less_Eigen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_LESS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_LESS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
//template< typename T1, typename T2, typename ET >
//auto
//operator<( const Eigen::MatrixBase<ET> &  aA,
//           Matrix< T1, T2 > & aB )
//->decltype( aA < aB.matrix_data() )
//{
//    return  aA < aB.matrix_data();
//}
//
//template< typename T1, typename T2, typename ET >
//auto
//operator<( Matrix< T1, T2 > & aA,
//           const Eigen::MatrixBase<ET> &  aB)
//->decltype( aA.matrix_data() < aB )
//{
//    return  aA.matrix_data() < aB;
//}

template< typename ET1, typename ET2 >
auto
operator<( const Eigen::MatrixBase<ET1> &  aA,
           const Eigen::MatrixBase<ET2> &  aB)
->decltype( aA.array() < aB.array() )
{
    return  aA.array() < aB.array();
}

}


#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_LESS_EIGEN_HPP_ */
