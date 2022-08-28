/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_times_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_TIMES_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_TIMES_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template< typename Matrix_Type, typename ET >
auto
operator*( const Eigen::MatrixBase<ET> &  aA,
           Matrix< Matrix_Type > & aB )
->decltype( aA * aB.matrix_data() )
{
    return  aA * aB.matrix_data();
}

template< typename Matrix_Type, typename ET >
auto
operator*( Matrix< Matrix_Type > & aA,
           const Eigen::MatrixBase<ET> &  aB)
->decltype( aA.matrix_data() * aB )
{
    return  aA.matrix_data() * aB;
}

template< typename Scalar, typename ET >
auto
scalar_times( Scalar & aA,
           const Eigen::MatrixBase<ET> &  aB)
->decltype( aA * aB )
{
    return  aA * aB;
}

}

/*template< typename Matrix_Type >
auto
operator*( const Matrix< Matrix_Type > & aA, const typename Matrix< Matrix_Type >::Data_Type & aB )
->decltype( aA.matrix_data() * aB )
{
    return  aA.matrix_data() * aB;
}

template< typename Matrix_Type  >
auto
operator*( const typename Matrix< Matrix_Type >::Data_Type & aA,  Matrix< Matrix_Type > & aB )
->decltype( aA * aB.matrix_data() )
{
    return  aA * aB.matrix_data();
}*/

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_TIMES_EIGEN_HPP_ */

